#!/usr/bin/env python3
"""Measure masked, N, and unmasked-callable sequence in event representative arms."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from statistics import median


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            else:
                sequences[name].append(line.strip())
    return {name: "".join(parts) for name, parts in sequences.items()}


def arm_qc(sequence: str) -> dict[str, int]:
    counts = Counter(sequence)
    return {
        "bp": len(sequence),
        "unmasked_acgt_bp": sum(counts[base] for base in "ACGT"),
        "softmasked_acgt_bp": sum(counts[base] for base in "acgt"),
        "n_bp": counts["N"] + counts["n"],
        "other_bp": len(sequence)
        - sum(counts[base] for base in "ACGTacgtNn"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--events", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    genome = read_fasta(args.fasta)
    with args.events.open() as handle:
        events = list(csv.DictReader(handle, delimiter="\t"))

    output_rows = []
    for event in events:
        row: dict[str, object] = {
            "event_id": event["event_id"],
            "strict_pd_status": event["strict_pd_status"],
            "strict_age_bin": event["strict_age_bin"],
        }
        for copy in (1, 2):
            chrom = event[f"representative_copy{copy}_chrom"]
            start = int(event[f"representative_copy{copy}_start"])
            end = int(event[f"representative_copy{copy}_end"])
            values = arm_qc(genome[chrom][start:end])
            for key, value in values.items():
                row[f"copy{copy}_{key}"] = value
            row[f"copy{copy}_unmasked_callable_fraction"] = (
                f"{values['unmasked_acgt_bp'] / values['bp']:.6f}"
            )
            row[f"copy{copy}_non_n_fraction"] = (
                f"{(values['bp'] - values['n_bp']) / values['bp']:.6f}"
            )
        minimum_callable = min(
            float(row["copy1_unmasked_callable_fraction"]),
            float(row["copy2_unmasked_callable_fraction"]),
        )
        minimum_non_n = min(
            float(row["copy1_non_n_fraction"]),
            float(row["copy2_non_n_fraction"]),
        )
        row["minimum_unmasked_callable_fraction"] = f"{minimum_callable:.6f}"
        row["minimum_non_n_fraction"] = f"{minimum_non_n:.6f}"
        output_rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(output_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    passed = [row for row in output_rows if row["strict_pd_status"] == "PASS"]
    fractions = [float(row["minimum_unmasked_callable_fraction"]) for row in passed]
    summary_rows = [
        {"metric": "strict_pd_candidates", "value": len(passed)},
        {
            "metric": "median_minimum_unmasked_callable_fraction",
            "value": f"{median(fractions):.6f}",
        },
    ]
    for threshold in (0.5, 0.7, 0.8, 0.9):
        summary_rows.append(
            {
                "metric": f"candidates_with_both_arms_unmasked_callable_fraction_ge_{threshold}",
                "value": sum(value >= threshold for value in fractions),
            }
        )
    summary_path = args.output.with_name("event_sequence_qc_summary.tsv")
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summary_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Annotated representative-arm sequence QC for {len(events)} events")


if __name__ == "__main__":
    main()
