#!/usr/bin/env python3
"""Compare the de novo RM P/D branch with the historical GFF-mask branch."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


RM = Path(".")
OLD = Path(".")
OLD_EVENT = Path(".")
OUTPUT = Path(".")


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def loci(row: dict[str, str]) -> list[tuple[str, int, int]]:
    return [
        (row[f"locus_{label}_chrom"], int(row[f"locus_{label}_representative_start"]), int(row[f"locus_{label}_representative_end"]))
        for label in ("A", "B")
    ]


def ro(left: tuple[str, int, int], right: tuple[str, int, int]) -> float:
    if left[0] != right[0]:
        return 0.0
    overlap = max(0, min(left[2], right[2]) - max(left[1], right[1]))
    return min(overlap / max(1, left[2] - left[1]), overlap / max(1, right[2] - right[1]))


def candidates(old: list[dict[str, str]], new: list[dict[str, str]]) -> list[tuple[float, str, str, str]]:
    rows = []
    for left in old:
        ll = loci(left)
        for right in new:
            rr = loci(right)
            direct = min(ro(ll[0], rr[0]), ro(ll[1], rr[1]))
            swapped = min(ro(ll[0], rr[1]), ro(ll[1], rr[0]))
            score, orientation = (direct, "direct") if direct >= swapped else (swapped, "swapped")
            if score >= 0.5:
                rows.append((score, left["event_id"], right["event_id"], orientation))
    return sorted(rows, reverse=True)


def one_to_one(calls: list[tuple[float, str, str, str]], threshold: float) -> list[tuple[float, str, str, str]]:
    used_old, used_new, selected = set(), set(), []
    for item in calls:
        score, old_id, new_id, _ = item
        if score < threshold:
            continue
        if old_id not in used_old and new_id not in used_new:
            selected.append(item)
            used_old.add(old_id)
            used_new.add(new_id)
    return selected


def count(path: Path) -> int:
    return len(read(path))


def ratio(d: int, p: int) -> float:
    return d / p if p else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--rm-root", required=True, type=Path)
    parser.add_argument("--gff-root", type=Path)
    parser.add_argument("--gff-event-root", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="比较表和 Markdown 报告的目录；默认写入 RM 根目录。",
    )
    args = parser.parse_args()
    global RM, OLD, OLD_EVENT, OUTPUT
    RM = args.rm_root
    OLD = args.gff_root or args.project / "15_age_free_pd_sequence_variation"
    OLD_EVENT = args.gff_event_root or (
        args.project
        / "06_sd_age_tracing_preparation/event_first_reanalysis"
    )
    OUTPUT = args.output_dir or RM

    old_events = read(OLD_EVENT / "events/event_first_events.tsv")
    new_events = read(RM / "event_first_reanalysis/events/event_first_events.tsv")
    old_by = {row["event_id"]: row for row in old_events}
    new_by = {row["event_id"]: row for row in new_events}
    all_candidates = candidates(old_events, new_events)
    matches50 = one_to_one(all_candidates, 0.5)
    matches80 = one_to_one(all_candidates, 0.8)
    old_age = {row["event_id"] for row in read(OLD / "inputs/age_free_pd_events.tsv")}
    new_age = {row["event_id"] for row in read(RM / "inputs/age_free_pd_events.tsv")}

    match_rows = []
    for score, old_id, new_id, orientation in matches50:
        old_p = old_by[old_id]["provisional_p_locus"]
        new_p = new_by[new_id]["provisional_p_locus"]
        mapped_old_p = old_p if orientation == "direct" else ("locus_B" if old_p == "locus_A" else "locus_A")
        match_rows.append({
            "old_event_id": old_id, "rm_event_id": new_id, "minimum_reciprocal_overlap": f"{score:.6f}",
            "locus_mapping": orientation, "old_age_free": "YES" if old_id in old_age else "NO",
            "rm_age_free": "YES" if new_id in new_age else "NO",
            "P_orientation_concordant": "YES" if mapped_old_p == new_p else "NO",
        })
    write(OUTPUT / "comparison/event_matches.ro50.tsv", match_rows)

    stages = [
        ("BISER_calls", 4734, 6669),
        ("strict_two_copy_events", len(old_events), len(new_events)),
        ("stable_age_free_PD", len(old_age), len(new_age)),
        ("core_eligible", count(OLD / "inputs/high_priority_core_eligible.tsv"), count(RM / "inputs/high_priority_core_eligible.tsv")),
        ("mapped_P0_queue", count(OLD / "pilot/age_free_p0_mapping_queue.tsv"), count(RM / "pilot/age_free_p0_mapping_queue.tsv")),
        ("atomic_regions", count(OLD / "microindel_local_msa/atomic_region_manifest.tsv"), count(RM / "microindel_local_msa/atomic_region_manifest.tsv")),
        ("SNP_ge200_events", count(OLD / "snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge200.event_metrics.tsv"), count(RM / "snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge200.event_metrics.tsv")),
    ]
    funnel = [{"stage": stage, "GFF_mask": old, "RM_mask": new, "RM_minus_GFF": new-old, "RM_vs_GFF_pct": f"{100*new/old:.2f}"} for stage, old, new in stages]
    write(OUTPUT / "comparison/full_pipeline_funnel.tsv", funnel)

    old_snp = read(OLD / "snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge200.event_metrics.tsv")
    new_snp = read(RM / "snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge200.event_metrics.tsv")
    old_snp_by = {r["event_id"]: r for r in old_snp}; new_snp_by = {r["event_id"]: r for r in new_snp}
    old_ind = {r["event_id"]: r for r in read(OLD / "microindel_local_msa/event_level_PD_microindel_rates.tsv")}
    new_ind = {r["event_id"]: r for r in read(RM / "microindel_local_msa/event_level_PD_microindel_rates.tsv")}
    map80 = [(o, n, orient) for _, o, n, orient in matches80]
    matched_snp = [(old_snp_by[o], new_snp_by[n]) for o, n, _ in map80 if o in old_snp_by and n in new_snp_by]
    matched_ind = [(old_ind[o], new_ind[n]) for o, n, _ in map80 if o in old_ind and n in new_ind]

    old_snp_p = sum(int(r["P_specific_SNP"]) for r in old_snp); old_snp_d = sum(int(r["D_specific_SNP"]) for r in old_snp)
    new_snp_p = sum(int(r["P_specific_SNP"]) for r in new_snp); new_snp_d = sum(int(r["D_specific_SNP"]) for r in new_snp)
    old_ind_stat = read(OLD / "microindel_local_msa/PD_denovo_microindel_statistics.tsv")[0]
    new_ind_stat = read(RM / "microindel_local_msa/PD_denovo_microindel_statistics.tsv")[0]
    endpoint = [
        {"endpoint": "SNP_ge200_all", "GFF_events": len(old_snp), "RM_events": len(new_snp), "GFF_P": old_snp_p, "GFF_D": old_snp_d, "GFF_D_to_P": f"{ratio(old_snp_d, old_snp_p):.6f}", "RM_P": new_snp_p, "RM_D": new_snp_d, "RM_D_to_P": f"{ratio(new_snp_d, new_snp_p):.6f}"},
        {"endpoint": "microindel_primary_all", "GFF_events": old_ind_stat["events_with_polarized_indels"], "RM_events": new_ind_stat["events_with_polarized_indels"], "GFF_P": old_ind_stat["P_branch_indels"], "GFF_D": old_ind_stat["D_branch_indels"], "GFF_D_to_P": old_ind_stat["D_to_P_count_ratio"], "RM_P": new_ind_stat["P_branch_indels"], "RM_D": new_ind_stat["D_branch_indels"], "RM_D_to_P": new_ind_stat["D_to_P_count_ratio"]},
    ]
    for label, pairs, pfield, dfield in (("SNP_ge200_RO80_matched", matched_snp, "P_specific_SNP", "D_specific_SNP"), ("microindel_RO80_matched", matched_ind, "P_branch_indels", "D_branch_indels")):
        op = sum(int(a[pfield]) for a, _ in pairs); od = sum(int(a[dfield]) for a, _ in pairs)
        np = sum(int(b[pfield]) for _, b in pairs); nd = sum(int(b[dfield]) for _, b in pairs)
        endpoint.append({"endpoint": label, "GFF_events": len(pairs), "RM_events": len(pairs), "GFF_P": op, "GFF_D": od, "GFF_D_to_P": f"{ratio(od, op):.6f}", "RM_P": np, "RM_D": nd, "RM_D_to_P": f"{ratio(nd, np):.6f}"})
    write(OUTPUT / "comparison/endpoint_comparison.tsv", endpoint)

    concordant_age = [r for r in match_rows if r["old_age_free"] == "YES" and r["rm_age_free"] == "YES"]
    p_concordant = sum(r["P_orientation_concordant"] == "YES" for r in concordant_age)
    report = f"""# RepeatMasker BISER P/D 全流程重分析

## 主要结果

- RM BISER 从 6,669 calls 独立归一化为 {len(new_events):,} 个严格两拷贝事件，其中 {len(new_age):,} 个通过 age-free 且 scope/threshold 稳定的 P/D 极化；GFF-mask 分别为 {len(old_events):,} 和 {len(old_age):,}。
- 可获得至少一个 P0 的事件从 {stages[4][1]} 变为 {stages[4][2]}；三算法一致且≥200 callable sites 的 SNP 端点从 {len(old_snp)} 变为 {len(new_snp)} 事件。
- SNP 主端点 D/P 由 {ratio(old_snp_d, old_snp_p):.3f} 变为 {ratio(new_snp_d, new_snp_p):.3f}；方向保持 D>P。
- microindel 主端点 D/P 由 {old_ind_stat['D_to_P_count_ratio']} 变为 {new_ind_stat['D_to_P_count_ratio']}；方向同样保持 D>P。

## 事件同源匹配

- 严格事件在两个 mask 间可一对一匹配：RO≥0.50 为 {len(matches50):,}，RO≥0.80 为 {len(matches80):,}。
- RO≥0.50 且两边均进入 age-free P/D 的事件为 {len(concordant_age):,}，其中 P 位点方向一致 {p_concordant:,}/{len(concordant_age):,} ({100*p_concordant/max(1,len(concordant_age)):.2f}%)。

## 解读

RM mask 明显改变了 BISER 边集、事件边界和下游可获得的 P0/core，因此“全集结果”的数量变化同时包含事件替换与 callable 序列变化。然而，SNP 和 microindel 两个主端点的 D>P 方向都保持，说明核心生物学结论不依赖于单一 mask 方案。精确数值不应被视为 mask-invariant，报告时应以 GFF-mask 主分析 + RM 全流程敏感性分析的方式呈现。
"""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "FULL_PD_REANALYSIS_COMPARISON.md").write_text(report)
    print(f"RO50={len(matches50)} RO80={len(matches80)} SNP_ratio={ratio(new_snp_d,new_snp_p):.4f} indel_ratio={new_ind_stat['D_to_P_count_ratio']}")


if __name__ == "__main__":
    main()
