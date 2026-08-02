#!/usr/bin/env python3
"""Summarize SNP types and genomic-context rates for the local-MSA endpoint."""

from __future__ import annotations

import argparse
import bisect
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


COMPLEMENT = str.maketrans("ACGT", "TGCA")
SIX_CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
CONTEXTS = (
    "CDS",
    "exon_nonCDS",
    "intron_or_gene_body",
    "promoter_2kb",
    "intergenic",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
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


def normalize_change(ancestral: str, derived: str) -> str:
    if ancestral in "AG":
        ancestral = ancestral.translate(COMPLEMENT)
        derived = derived.translate(COMPLEMENT)
    return f"{ancestral}>{derived}"


def exact_binomial(successes: int, trials: int) -> float:
    if not trials:
        return math.nan
    observed = math.comb(trials, successes)
    return min(
        1.0,
        sum(
            math.comb(trials, value)
            for value in range(trials + 1)
            if math.comb(trials, value) <= observed
        )
        / (2**trials),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--sites", required=True, type=Path)
    parser.add_argument("--events", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    events = {row["event_id"]: row for row in read_tsv(args.events)}
    sites = [
        row for row in read_tsv(args.sites) if row["event_id"] in events
    ]

    spectrum = Counter()
    variant_rows = []
    for row in sites:
        category = row["polarized_class"]
        if category not in {"P_specific", "D_specific"}:
            continue
        role = category[0]
        derived = row[f"{role}_base"]
        ancestral = row["ancestral_base"]
        change = normalize_change(ancestral, derived)
        spectrum[(role, change)] += 1
        spectrum[(role, "transition" if change in {"C>T", "T>C"} else "transversion")] += 1
        variant_rows.append(
            {
                "event_id": row["event_id"],
                "age_bin": row["age_bin"],
                "copy_role": role,
                "P_chrom": row["P_chrom"],
                "P_position_0based": row["P_position_0based"],
                "D_chrom": row["D_chrom"],
                "D_position_0based": row["D_position_0based"],
                "ancestral_base": ancestral,
                "derived_base": derived,
                "six_class_change": change,
                "transition_transversion": (
                    "transition" if change in {"C>T", "T>C"} else "transversion"
                ),
            }
        )
    write_tsv(args.output / "polarized_local_MSA_SNP_sites.tsv", variant_rows)
    total_callable = len(sites)
    spectrum_rows = []
    for role in ("P", "D"):
        for change in SIX_CLASSES:
            count = spectrum[(role, change)]
            spectrum_rows.append(
                {
                    "copy_role": role,
                    "six_class_change": change,
                    "SNP_count": count,
                    "callable_sites": total_callable,
                    "SNP_per_callable_site": f"{count / total_callable:.8f}",
                }
            )
    write_tsv(args.output / "six_class_SNP_spectrum.tsv", spectrum_rows)

    positions = []
    for row in sites:
        for role in ("P", "D"):
            positions.append(
                {
                    "event_id": row["event_id"],
                    "age_bin": row["age_bin"],
                    "copy_role": role,
                    "chrom": row[f"{role}_chrom"],
                    "position0": int(row[f"{role}_position_0based"]),
                    "is_role_specific_change": (
                        row["polarized_class"] == f"{role}_specific"
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
                left, right = bisect.bisect_left(values, pstart), bisect.bisect_left(values, pend)
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
        row["repeat_overlap"] = "YES" if "repeat" in flags else "NO"
        row.pop("flags")
    write_tsv(args.output / "local_MSA_callable_genomic_positions.tsv", positions)

    aggregate = Counter()
    event_counts = Counter()
    for row in positions:
        role, context = str(row["copy_role"]), str(row["context"])
        aggregate[(role, context, "callable")] += 1
        event_counts[(row["event_id"], role, context, "callable")] += 1
        if row["is_role_specific_change"]:
            aggregate[(role, context, "change")] += 1
            event_counts[(row["event_id"], role, context, "change")] += 1
    context_rows = []
    for role in ("P", "D"):
        for context in CONTEXTS:
            denominator = aggregate[(role, context, "callable")]
            numerator = aggregate[(role, context, "change")]
            context_rows.append(
                {
                    "copy_role": role,
                    "genomic_context": context,
                    "callable_sites": denominator,
                    "role_specific_SNP": numerator,
                    "SNP_per_callable_site": (
                        f"{numerator / denominator:.8f}" if denominator else "NA"
                    ),
                }
            )
    write_tsv(args.output / "genomic_context_SNP_rates.tsv", context_rows)

    statistical_rows = []
    for context in CONTEXTS:
        p_rates, d_rates = [], []
        d_greater = p_greater = ties = 0
        for event_id in events:
            p_den = event_counts[(event_id, "P", context, "callable")]
            d_den = event_counts[(event_id, "D", context, "callable")]
            if not p_den or not d_den:
                continue
            p_rate = event_counts[(event_id, "P", context, "change")] / p_den
            d_rate = event_counts[(event_id, "D", context, "change")] / d_den
            p_rates.append(p_rate)
            d_rates.append(d_rate)
            d_greater += d_rate > p_rate
            p_greater += p_rate > d_rate
            ties += p_rate == d_rate
        non_ties = d_greater + p_greater
        statistical_rows.append(
            {
                "genomic_context": context,
                "paired_events": len(p_rates),
                "median_P_rate": f"{np.median(p_rates):.8f}" if p_rates else "NA",
                "median_D_rate": f"{np.median(d_rates):.8f}" if d_rates else "NA",
                "events_D_greater": d_greater,
                "events_P_greater": p_greater,
                "events_tied": ties,
                "event_sign_test_p": (
                    f"{exact_binomial(d_greater, non_ties):.8g}"
                    if non_ties else "NA"
                ),
                "paired_wilcoxon_p": (
                    f"{stats.wilcoxon(d_rates, p_rates, zero_method='wilcox').pvalue:.8g}"
                    if non_ties else "NA"
                ),
            }
        )
    write_tsv(args.output / "genomic_context_event_statistics.tsv", statistical_rows)

    event_rows = list(events.values())
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    age_labels = [
        age
        for age in ("time1", "time2", "time3", "age_free")
        if any(row["age_bin"] == age for row in event_rows)
    ]
    x = np.arange(len(age_labels))
    width = 0.36
    age_values = {}
    for age in age_labels:
        subset = [row for row in event_rows if row["age_bin"] == age]
        denominator = sum(
            int(row["local_MSA_callable_sites"]) for row in subset
        )
        age_values[(age, "P")] = (
            sum(int(row["P_specific_SNP"]) for row in subset) / denominator
            if denominator
            else 0
        )
        age_values[(age, "D")] = (
            sum(int(row["D_specific_SNP"]) for row in subset) / denominator
            if denominator
            else 0
        )
    axes[0].bar(
        x - width / 2,
        [age_values[(age, "P")] for age in age_labels],
        width,
        label="P",
    )
    axes[0].bar(
        x + width / 2,
        [age_values[(age, "D")] for age in age_labels],
        width,
        label="D",
    )
    axes[0].set_xticks(x, age_labels)
    axes[0].set_ylabel("SNP / local-MSA callable site")
    axes[0].set_title("A  Age-stratified endpoint", loc="left")
    axes[0].legend(frameon=False)

    sx = np.arange(len(SIX_CLASSES))
    axes[1].bar(sx - width / 2, [spectrum[("P", c)] / total_callable for c in SIX_CLASSES], width, label="P")
    axes[1].bar(sx + width / 2, [spectrum[("D", c)] / total_callable for c in SIX_CLASSES], width, label="D")
    axes[1].set_xticks(sx, SIX_CLASSES, rotation=45)
    axes[1].set_ylabel("SNP / callable site")
    axes[1].set_title("B  Six-class spectrum", loc="left")

    cx = np.arange(len(CONTEXTS))
    def context_rate(role: str, context: str) -> float:
        den = aggregate[(role, context, "callable")]
        return aggregate[(role, context, "change")] / den if den else 0
    axes[2].bar(cx - width / 2, [context_rate("P", c) for c in CONTEXTS], width, label="P")
    axes[2].bar(cx + width / 2, [context_rate("D", c) for c in CONTEXTS], width, label="D")
    axes[2].set_xticks(cx, ("CDS", "exon", "intron/gene", "promoter", "intergenic"), rotation=35, ha="right")
    axes[2].set_ylabel("Context-specific SNP rate")
    axes[2].set_title("C  Genomic context", loc="left")
    for suffix in ("png", "pdf"):
        fig.savefig(args.output / f"local_MSA_SNP_types_and_context.{suffix}", dpi=300)
    plt.close(fig)
    print(f"events={len(events)} callable={len(sites)} variants={len(variant_rows)}")


if __name__ == "__main__":
    main()
