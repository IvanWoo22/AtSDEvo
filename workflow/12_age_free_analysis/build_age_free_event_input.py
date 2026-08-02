#!/usr/bin/env python3
"""Materialize the scope/threshold-stable, age-free P/D event input."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


ROOT = Path(os.environ.get("SD_PROJECT_ROOT", Path(__file__).resolve().parents[2]))
OUT = Path(
    os.environ.get(
        "SD_AGE_FREE_ROOT", ROOT / "15_age_free_pd_sequence_variation"
    )
)


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", type=Path)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    events_path = args.events or (
        ROOT
        / "06_sd_age_tracing_preparation/event_first_reanalysis/events/"
        "event_first_events.tsv"
    )
    audit_path = args.audit or (
        ROOT
        / "14_pd_polarization_reassessment/events/"
        "age_free_pd_polarization_audit.tsv"
    )
    output_root = args.output_root or OUT

    events = {row["event_id"]: row for row in read(events_path)}
    audit = read(audit_path)

    rows: list[dict[str, str]] = []
    for decision in audit:
        if decision["scope_and_threshold_stable"] != "PASS":
            continue
        row = dict(events[decision["event_id"]])
        p_locus = decision["age_free_P_locus"]
        row["former_strict_age_bin"] = row["strict_age_bin"]
        row["former_strict_pd_status"] = row["strict_pd_status"]
        row["former_provisional_p_locus"] = row["provisional_p_locus"]
        row["strict_pd_status"] = "PASS"
        row["strict_age_bin"] = "age_free"
        row["provisional_p_locus"] = p_locus
        row["provisional_d_locus"] = (
            "locus_B" if p_locus == "locus_A" else "locus_A"
        )
        row["classification_stable_at_overlap_0.25_0.5_0.75"] = "PASS"
        row["age_free_PD_evidence_class"] = decision["age_free_PD_status"]
        row["age_free_single_copy_support_nodes"] = decision[
            "single_copy_support_nodes"
        ]
        row["age_free_scope_threshold_rule"] = (
            "same_P_locus_in_representative_common_core_envelope_at_"
            "overlap_0.25_0.50_0.75"
        )
        rows.append(row)

    rows.sort(key=lambda row: row["event_id"])
    write(output_root / "inputs/age_free_pd_events.tsv", rows)
    print(f"age_free_pd_events={len(rows)}")


if __name__ == "__main__":
    main()
