#!/usr/bin/env python3
"""Map all callable P/D positions to genomic contexts for proper rate denominators."""

from __future__ import annotations

import argparse
import bisect
import csv
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats


CONTEXTS = ("CDS", "exon_nonCDS", "intron_or_gene_body", "promoter_2kb", "intergenic")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        if not rows:
            return
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    expansion = args.project / "08_event_inclusion_sensitivity"
    pilot = expansion / "core_500"
    controlled = read_tsv(expansion / "controlled_expansion_endpoint_events.tsv")
    tiers = {row["event_id"]: row["admission_tier"] for row in controlled}
    sources = {
        "core_500_strict": pilot / "sequence_variation/polarized_sites.tsv",
        "partial_postduplication": (
            pilot / "sequence_variation_partial_postdup/polarized_sites.tsv"
        ),
        "deeper_P0_fallback": (
            pilot / "sequence_variation_deeper_P0/polarized_sites.tsv"
        ),
    }
    rows_by_event: dict[str, list[dict[str, str]]] = defaultdict(list)
    for tier, path in sources.items():
        wanted = {event for event, value in tiers.items() if value == tier}
        for row in read_tsv(path):
            if row["event_id"] in wanted:
                rows_by_event[row["event_id"]].append(row)
    blocks: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(pilot / "core/homologous_core_blocks.tsv"):
        if row["event_id"] in tiers:
            blocks[row["event_id"]].append(row)
    for rows in blocks.values():
        rows.sort(key=lambda row: int(row["core_block_index"]))

    positions = []
    for event in controlled:
        event_id = event["event_id"]
        event_blocks = blocks[event_id]
        for source in rows_by_event[event_id]:
            core_position = int(source["core_position_1based"]) - 1
            cumulative = 0
            selected = None
            offset = -1
            for block in event_blocks:
                length = int(block["core_block_bp"])
                if cumulative <= core_position < cumulative + length:
                    selected = block
                    offset = core_position - cumulative
                    break
                cumulative += length
            if selected is None:
                raise AssertionError(f"Missing block for {event_id}/{core_position}")
            for role in ("P", "D"):
                copy = "copy1" if selected["copy1_role"] == role else "copy2"
                strand = selected[f"{copy}_strand"]
                start, end = int(selected[f"{copy}_start"]), int(selected[f"{copy}_end"])
                position = start + offset if strand == "+" else end - 1 - offset
                positions.append(
                    {
                        "event_id": event_id,
                        "age_bin": event["age_bin"],
                        "copy_role": role,
                        "chrom": selected[f"{copy}_chrom"],
                        "position0": position,
                        "is_role_specific_change": (
                            source["primary_class"] == f"{role}_specific"
                        ),
                        "flags": set(),
                    }
                )
    by_chrom: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in positions:
        by_chrom[str(row["chrom"])].append(row)
    sorted_positions = {}
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda row: int(row["position0"]))
        sorted_positions[chrom] = [int(row["position0"]) for row in rows]
    gff = args.project / "01_reference/prepared_data/TAIR12.Col-CC.annotation.gff3"
    with gff.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, feature = fields[0], fields[2]
            if chrom not in by_chrom or feature not in {"gene", "exon", "CDS"}:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            values = sorted_positions[chrom]
            left, right = bisect.bisect_left(values, start), bisect.bisect_left(values, end)
            for row in by_chrom[chrom][left:right]:
                row["flags"].add(feature)  # type: ignore[union-attr]
            if feature == "gene":
                pstart, pend = (
                    (max(0, start - 2000), start)
                    if fields[6] == "+"
                    else (end, end + 2000)
                )
                left = bisect.bisect_left(values, pstart)
                right = bisect.bisect_left(values, pend)
                for row in by_chrom[chrom][left:right]:
                    row["flags"].add("promoter_2kb")  # type: ignore[union-attr]
    repeat_path = (
        args.project / "01_reference/prepared_data/TAIR12.annotated_repeats.merged.bed"
    )
    with repeat_path.open() as handle:
        for line in handle:
            chrom, start_text, end_text, *_ = line.rstrip().split("\t")
            if chrom not in by_chrom:
                continue
            values = sorted_positions[chrom]
            left = bisect.bisect_left(values, int(start_text))
            right = bisect.bisect_left(values, int(end_text))
            for row in by_chrom[chrom][left:right]:
                row["flags"].add("repeat")  # type: ignore[union-attr]
    for row in positions:
        flags: set[str] = row["flags"]  # type: ignore[assignment]
        row["context"] = (
            "CDS"
            if "CDS" in flags
            else "exon_nonCDS"
            if "exon" in flags
            else "intron_or_gene_body"
            if "gene" in flags
            else "promoter_2kb"
            if "promoter_2kb" in flags
            else "intergenic"
        )
        row["repeat_overlap"] = "PASS" if "repeat" in flags else "FAIL"

    aggregate = Counter()
    event_counts = Counter()
    for row in positions:
        key = (str(row["copy_role"]), str(row["context"]))
        aggregate[(*key, "callable")] += 1
        if row["is_role_specific_change"]:
            aggregate[(*key, "change")] += 1
        event_key = (str(row["event_id"]), *key)
        event_counts[(*event_key, "callable")] += 1
        if row["is_role_specific_change"]:
            event_counts[(*event_key, "change")] += 1
    aggregate_rows = []
    for role in ("P", "D"):
        for context in CONTEXTS:
            callable_count = aggregate[(role, context, "callable")]
            change_count = aggregate[(role, context, "change")]
            aggregate_rows.append(
                {
                    "copy_role": role,
                    "gene_context": context,
                    "callable_sites": callable_count,
                    "role_specific_changes": change_count,
                    "role_specific_change_rate": (
                        f"{change_count / callable_count:.8f}"
                        if callable_count
                        else "NA"
                    ),
                }
            )
    write_tsv(args.output / "snp_context_callable_rate_summary.tsv", aggregate_rows)

    statistical = []
    rng = np.random.default_rng(20260724)
    for context in CONTEXTS:
        p_rates, d_rates = [], []
        usable = 0
        for event in controlled:
            event_id = event["event_id"]
            p_den = event_counts[(event_id, "P", context, "callable")]
            d_den = event_counts[(event_id, "D", context, "callable")]
            if not p_den or not d_den:
                continue
            usable += 1
            p_rates.append(event_counts[(event_id, "P", context, "change")] / p_den)
            d_rates.append(event_counts[(event_id, "D", context, "change")] / d_den)
        p = np.array(p_rates)
        d = np.array(d_rates)
        differences = d - p
        if usable:
            samples = rng.integers(0, usable, size=(20_000, usable))
            boot = np.median(differences[samples], axis=1)
        else:
            boot = np.array([np.nan])
        nonzero = differences[differences != 0]
        statistical.append(
            {
                "gene_context": context,
                "events_with_both_copy_denominators": usable,
                "median_P_rate": f"{np.median(p):.8f}" if usable else "NA",
                "median_D_rate": f"{np.median(d):.8f}" if usable else "NA",
                "median_D_minus_P_rate": (
                    f"{np.median(differences):.8f}" if usable else "NA"
                ),
                "median_difference_CI95_low": (
                    f"{np.quantile(boot, 0.025):.8f}" if usable else "NA"
                ),
                "median_difference_CI95_high": (
                    f"{np.quantile(boot, 0.975):.8f}" if usable else "NA"
                ),
                "events_D_gt_P": int(np.sum(differences > 0)),
                "events_P_gt_D": int(np.sum(differences < 0)),
                "event_sign_test_p": (
                    f"{stats.binomtest(int(np.sum(nonzero > 0)), len(nonzero), 0.5).pvalue:.12g}"
                    if len(nonzero)
                    else "1"
                ),
                "paired_wilcoxon_p": (
                    f"{stats.wilcoxon(d, p, zero_method='wilcox').pvalue:.12g}"
                    if np.any(differences)
                    else "1"
                ),
            }
        )
    write_tsv(args.output / "snp_context_event_statistical_summary.tsv", statistical)
    repeat_summary = Counter()
    for row in positions:
        repeat_summary[
            (str(row["copy_role"]), str(row["repeat_overlap"]), "callable")
        ] += 1
        if row["is_role_specific_change"]:
            repeat_summary[
                (str(row["copy_role"]), str(row["repeat_overlap"]), "change")
            ] += 1
    write_tsv(
        args.output / "snp_repeat_callable_rate_summary.tsv",
        [
            {
                "copy_role": role,
                "repeat_overlap": status,
                "callable_sites": repeat_summary[(role, status, "callable")],
                "role_specific_changes": repeat_summary[(role, status, "change")],
                "role_specific_change_rate": (
                    f"{repeat_summary[(role, status, 'change')] / repeat_summary[(role, status, 'callable')]:.8f}"
                    if repeat_summary[(role, status, "callable")]
                    else "NA"
                ),
            }
            for role in ("P", "D")
            for status in ("PASS", "FAIL")
        ],
    )
    print(f"Annotated {len(positions)} role-specific callable genomic positions")


if __name__ == "__main__":
    main()
