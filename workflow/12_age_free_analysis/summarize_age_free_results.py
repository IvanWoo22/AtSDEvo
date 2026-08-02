#!/usr/bin/env python3
"""Freeze the age-free P/D funnel, old/new comparison, and paired statistics."""

from __future__ import annotations

import csv
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


PROJECT = Path(
    os.environ.get("SD_PROJECT_ROOT", Path(__file__).resolve().parents[2])
)
ROOT = Path(
    os.environ.get(
        "SD_AGE_FREE_ROOT", PROJECT / "15_age_free_pd_sequence_variation"
    )
)
OLD = PROJECT / "12_inclusive_pd_sequence_variation"
GENOME_BP = 142_481_245
RNG = np.random.default_rng(20260729)


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def interval_metrics(rows: list[dict[str, str]]) -> tuple[int, int]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    summed = 0
    for row in rows:
        for locus in ("A", "B"):
            chrom = row[f"locus_{locus}_chrom"]
            start = int(row[f"locus_{locus}_representative_start"])
            end = int(row[f"locus_{locus}_representative_end"])
            by_chrom[chrom].append((start, end))
            summed += end - start
    union = 0
    for intervals in by_chrom.values():
        intervals.sort()
        left, right = intervals[0]
        for start, end in intervals[1:]:
            if start <= right:
                right = max(right, end)
            else:
                union += right - left
                left, right = start, end
        union += right - left
    return summed, union


def event_ids(path: Path, predicate=lambda row: True) -> set[str]:
    return {row["event_id"] for row in read(path) if predicate(row)}


def bootstrap_paired(
    rows: list[dict[str, str]],
    p_field: str,
    d_field: str,
    denominator_field: str,
    replicates: int = 20_000,
) -> dict[str, object]:
    p = np.array([int(row[p_field]) for row in rows], dtype=float)
    d = np.array([int(row[d_field]) for row in rows], dtype=float)
    den = np.array([int(row[denominator_field]) for row in rows], dtype=float)
    p_rate, d_rate = p / den, d / den
    difference = d_rate - p_rate
    nonzero = difference != 0
    wilcoxon = (
        stats.wilcoxon(d_rate, p_rate, zero_method="wilcox").pvalue
        if np.any(nonzero)
        else math.nan
    )
    boot_diff = np.empty(replicates)
    boot_ratio = np.empty(replicates)
    for index in range(replicates):
        sample = RNG.integers(0, len(rows), len(rows))
        boot_diff[index] = np.mean(d_rate[sample] - p_rate[sample])
        p_sum, d_sum = p[sample].sum(), d[sample].sum()
        boot_ratio[index] = d_sum / p_sum if p_sum else np.nan
    return {
        "events": len(rows),
        "P_count": int(p.sum()),
        "D_count": int(d.sum()),
        "D_to_P_count_ratio": f"{d.sum() / p.sum():.8f}" if p.sum() else "NA",
        "mean_event_D_minus_P_rate": f"{difference.mean():.10f}",
        "mean_event_D_minus_P_rate_bootstrap95_low": (
            f"{np.quantile(boot_diff, 0.025):.10f}"
        ),
        "mean_event_D_minus_P_rate_bootstrap95_high": (
            f"{np.quantile(boot_diff, 0.975):.10f}"
        ),
        "D_to_P_count_ratio_bootstrap95_low": (
            f"{np.nanquantile(boot_ratio, 0.025):.8f}"
        ),
        "D_to_P_count_ratio_bootstrap95_high": (
            f"{np.nanquantile(boot_ratio, 0.975):.8f}"
        ),
        "paired_wilcoxon_p": f"{wilcoxon:.8g}",
        "bootstrap_replicates": replicates,
        "bootstrap_seed": 20260729,
    }


all_events = read(
    PROJECT
    / "06_sd_age_tracing_preparation/event_first_reanalysis/events/"
    "event_first_events.tsv"
)
by_id = {row["event_id"]: row for row in all_events}
age_free = read(ROOT / "inputs/high_priority_core_eligible.tsv")
queue = read(ROOT / "pilot/age_free_p0_mapping_queue.tsv")
atom_ids = event_ids(ROOT / "microindel_local_msa/atomic_region_manifest.tsv")
snp20_ids = event_ids(
    ROOT / "snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge20.event_metrics.tsv"
)
snp200 = read(
    ROOT / "snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge200.event_metrics.tsv"
)
snp200_ids = {row["event_id"] for row in snp200}
indel_ids = event_ids(
    ROOT / "microindel_local_msa/denovo_microindel_inference.tsv",
    lambda row: row["evidence_tier"] == "PRIMARY",
)

funnel_rows = []
stages = (
    ("normalized_two_copy_events", {row["event_id"] for row in all_events}),
    ("scope_threshold_stable_age_free_PD", {row["event_id"] for row in age_free}),
    ("at_least_one_event_matched_P0", {row["event_id"] for row in queue}),
    ("continuous_PD_P0_atoms", atom_ids),
    ("three_aligner_local_MSA_ge20", snp20_ids),
    ("three_aligner_local_MSA_ge200", snp200_ids),
    ("at_least_one_primary_microindel", indel_ids),
)
previous = len(all_events)
for order, (stage, ids) in enumerate(stages, 1):
    rows = [by_id[event_id] for event_id in ids]
    summed, union = interval_metrics(rows)
    funnel_rows.append(
        {
            "stage_order": order,
            "stage": stage,
            "events": len(ids),
            "retained_from_previous_pct": f"{100 * len(ids) / previous:.4f}",
            "two_copy_interval_sum_bp": summed,
            "genome_union_covered_bp": union,
            "genome_union_coverage_pct": f"{100 * union / GENOME_BP:.5f}",
        }
    )
    previous = len(ids)
write(ROOT / "statistics/age_free_filter_funnel_with_coverage.tsv", funnel_rows)

old_queue = read(OLD / "pilot/inclusive_p0_mapping_queue.tsv")
old_snp = next(
    row
    for row in read(OLD / "snp_local_msa/local_MSA_endpoint_summary.tsv")
    if row["endpoint"] == "MAFFT_MUSCLE_PRANK_local_MSA"
    and row["minimum_local_MSA_callable_sites"] == "200"
)
new_snp = next(
    row
    for row in read(ROOT / "snp_local_msa/local_MSA_endpoint_summary.tsv")
    if row["endpoint"] == "MAFFT_MUSCLE_PRANK_local_MSA"
    and row["minimum_local_MSA_callable_sites"] == "200"
)
old_indel = read(OLD / "microindel_local_msa/PD_denovo_microindel_statistics.tsv")[0]
new_indel = read(ROOT / "microindel_local_msa/PD_denovo_microindel_statistics.tsv")[0]
old_atom_ids = event_ids(OLD / "microindel_local_msa/atomic_region_manifest.tsv")

comparison = [
    {
        "endpoint": "events_with_mapped_P0",
        "former_time_gated": len(old_queue),
        "age_free_PD": len(queue),
        "fold_change": f"{len(queue) / len(old_queue):.4f}",
    },
    {
        "endpoint": "events_with_continuous_atoms",
        "former_time_gated": len(old_atom_ids),
        "age_free_PD": len(atom_ids),
        "fold_change": f"{len(atom_ids) / len(old_atom_ids):.4f}",
    },
]
for label, field in (
    ("SNP_events_ge200", "events"),
    ("SNP_callable_sites_ge200", "callable_sites"),
    ("P_specific_SNP", "P_specific_SNP"),
    ("D_specific_SNP", "D_specific_SNP"),
):
    old_value, new_value = int(old_snp[field]), int(new_snp[field])
    comparison.append(
        {
            "endpoint": label,
            "former_time_gated": old_value,
            "age_free_PD": new_value,
            "fold_change": f"{new_value / old_value:.4f}",
        }
    )
for label, field in (
    ("primary_microindels", "total_indels"),
    ("P_branch_microindels", "P_branch_indels"),
    ("D_branch_microindels", "D_branch_indels"),
    ("events_with_primary_microindels", "events_with_polarized_indels"),
):
    old_value, new_value = int(old_indel[field]), int(new_indel[field])
    comparison.append(
        {
            "endpoint": label,
            "former_time_gated": old_value,
            "age_free_PD": new_value,
            "fold_change": f"{new_value / old_value:.4f}",
        }
    )
write(ROOT / "statistics/former_vs_age_free_endpoint_comparison.tsv", comparison)

old_ids, new_ids = (
    {row["event_id"] for row in old_queue},
    {row["event_id"] for row in queue},
)
write(
    ROOT / "statistics/P0_queue_overlap_audit.tsv",
    [
        {"category": "former_queue", "events": len(old_ids)},
        {"category": "age_free_queue", "events": len(new_ids)},
        {"category": "shared", "events": len(old_ids & new_ids)},
        {"category": "former_only", "events": len(old_ids - new_ids)},
        {"category": "age_free_only", "events": len(new_ids - old_ids)},
    ],
)

former_age = Counter(row["former_age_bin"] for row in queue)
write(
    ROOT / "statistics/age_free_P0_queue_former_age_labels.tsv",
    [
        {"former_age_label_audit_only": label, "events": count}
        for label, count in sorted(former_age.items())
    ],
)

snp_stats = bootstrap_paired(
    snp200,
    "P_specific_SNP",
    "D_specific_SNP",
    "local_MSA_callable_sites",
)
snp_stats["analysis"] = "three_aligner_local_MSA_SNP_ge200"
indel_event_rows = read(
    ROOT / "microindel_local_msa/event_level_PD_microindel_rates.tsv"
)
indel_stats = bootstrap_paired(
    indel_event_rows,
    "P_branch_indels",
    "D_branch_indels",
    "ancestral_callable_bp",
)
indel_stats["analysis"] = "primary_denovo_microindel_1_10bp"
write(ROOT / "statistics/paired_event_bootstrap_statistics.tsv", [snp_stats, indel_stats])

# Compact advisor-facing overview.
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
labels = ["events", "P0", "atoms", "SNP≥200"]
values = [1655, len(age_free), len(queue), len(atom_ids), len(snp200_ids)]
axes[0].plot(range(len(values)), values, marker="o", color="#365f91", linewidth=2)
axes[0].set_xticks(range(len(values)), ["events", "P/D", "P0", "atoms", "SNP≥200"], rotation=25)
axes[0].set_ylabel("Nonredundant SD events")
axes[0].set_title("A  Age-free event funnel", loc="left")
for x, value in enumerate(values):
    axes[0].text(x, value, str(value), ha="center", va="bottom", fontsize=9)

axes[1].bar(
    np.arange(2) - 0.18,
    [int(old_snp["P_specific_SNP"]), int(new_snp["P_specific_SNP"])],
    0.36,
    label="P",
    color="#3973ac",
)
axes[1].bar(
    np.arange(2) + 0.18,
    [int(old_snp["D_specific_SNP"]), int(new_snp["D_specific_SNP"])],
    0.36,
    label="D",
    color="#c45145",
)
axes[1].set_xticks(range(2), ["time-gated\n70 events", "age-free\n195 events"])
axes[1].set_ylabel("Polarized SNPs")
axes[1].set_title("B  Three-aligner SNP endpoint", loc="left")
axes[1].legend(frameon=False)

axes[2].bar(
    np.arange(2) - 0.18,
    [int(old_indel["P_branch_indels"]), int(new_indel["P_branch_indels"])],
    0.36,
    label="P",
    color="#3973ac",
)
axes[2].bar(
    np.arange(2) + 0.18,
    [int(old_indel["D_branch_indels"]), int(new_indel["D_branch_indels"])],
    0.36,
    label="D",
    color="#c45145",
)
axes[2].set_xticks(range(2), ["time-gated\n110 events", "age-free\n239 events"])
axes[2].set_ylabel("Primary 1–10 bp micro-indels")
axes[2].set_title("C  Independent local-MSA indels", loc="left")
axes[2].legend(frameon=False)

figure_dir = ROOT / "figures"
figure_dir.mkdir(exist_ok=True)
for suffix in ("png", "pdf"):
    fig.savefig(figure_dir / f"age_free_PD_result_overview.{suffix}", dpi=300)
plt.close(fig)

print(
    f"age_free={len(age_free)} P0={len(queue)} atoms={len(atom_ids)} "
    f"SNP_ge200={len(snp200_ids)} microindel_events={len(indel_ids)}"
)
