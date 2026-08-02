#!/usr/bin/env python3
"""Build an age-free queue requiring any mapped singleton-state P0."""

from __future__ import annotations

import argparse
import csv
import os
from collections import Counter
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
ROOT = Path(
    os.environ.get(
        "SD_AGE_FREE_ROOT",
        Path(os.environ.get("SD_PROJECT_ROOT", Path.cwd()))
        / "15_age_free_pd_sequence_variation",
    )
)


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path)
    args = parser.parse_args()
    root = args.root or ROOT
    events = read(root / "inputs/high_priority_core_eligible.tsv")
    qc = {row["event_id"]: row for row in read(root / "core/core_event_qc.tsv")}
    mapping = {
        (row["event_id"], row["species"]): row
        for row in read(root / "outgroup_mapping/event_species_blastn_summary.tsv")
    }

    queue: list[dict[str, object]] = []
    audit: list[dict[str, object]] = []
    for event in events:
        event_id = event["event_id"]
        candidates = [
            species for species in SPECIES if event[f"{species}_state"] in {"1", "2"}
        ]
        passed = [
            species
            for species in candidates
            if mapping.get((event_id, species), {}).get("mapping_status") == "PASS"
        ]
        postdup_pass = [
            species
            for species in SPECIES
            if event[f"{species}_state"] == "3"
            and mapping.get((event_id, species), {}).get("mapping_status") == "PASS"
        ]
        audit.append(
            {
                "event_id": event_id,
                "primary_pattern": event["primary_pattern"],
                "candidate_P0_species": ",".join(candidates),
                "mapped_P0_species": ",".join(passed),
                "mapped_P0_species_count": len(passed),
                "queue_status": "PASS" if passed else "FAIL",
            }
        )
        if not passed:
            continue
        queue.append(
            {
                "event_id": event_id,
                "age_bin": "age_free",
                "former_age_bin": event["former_strict_age_bin"],
                "analysis_tier": "MULTI_P0" if len(passed) >= 2 else "SINGLE_P0",
                "boundary_P0_species": passed[0],
                "mapped_P0_species": ",".join(passed),
                "mapped_P0_species_count": len(passed),
                "mapped_postduplication_species": ",".join(postdup_pass),
                "mapped_postduplication_species_count": len(postdup_pass),
                "paired_M_core_bp": qc[event_id]["paired_M_core_bp"],
                "jointly_callable_uppercase_acgt_bp": qc[event_id]["jointly_callable_uppercase_acgt_bp"],
                "present_day_PD_mismatch_pct_callable": qc[event_id]["present_day_PD_mismatch_pct_callable"],
                "age_filter_policy": "NOT_GATING",
                "legacy_1kb_filters_policy": "DIAGNOSTIC_NOT_GATING",
            }
        )

    write(root / "pilot/age_free_p0_mapping_queue.tsv", queue)
    write(root / "statistics/age_free_p0_event_audit.tsv", audit)
    counts = Counter(row["analysis_tier"] for row in queue)
    write(
        root / "statistics/age_free_p0_queue_summary.tsv",
        [
            {"metric": "age_free_core_events", "count": len(events)},
            {"metric": "at_least_one_mapped_P0", "count": len(queue)},
            {"metric": "single_mapped_P0", "count": counts["SINGLE_P0"]},
            {"metric": "multiple_mapped_P0", "count": counts["MULTI_P0"]},
        ],
    )
    print(f"age_free_core_events={len(events)} mapped_P0_queue={len(queue)}")


if __name__ == "__main__":
    main()
