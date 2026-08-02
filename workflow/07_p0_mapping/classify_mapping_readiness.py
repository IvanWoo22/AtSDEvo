#!/usr/bin/env python3
"""Classify event readiness from boundary and post-duplication locus mapping."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
BOUNDARY = {"time1": "Alyrata", "time2": "Bstricta", "time3": "Dstrictus"}
POSTDUP = {
    "time1": (),
    "time2": ("Alyrata",),
    "time3": ("Alyrata", "Bstricta"),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", required=True, type=Path)
    args = parser.parse_args()

    events = read_tsv(args.pilot / "inputs/high_priority_core_eligible.tsv")
    mapping = {
        (row["event_id"], row["species"]): row
        for row in read_tsv(
            args.pilot / "outgroup_mapping/event_species_blastn_summary.tsv"
        )
    }
    output = []
    for event in events:
        event_id = event["event_id"]
        age = event["strict_age_bin"]
        boundary_species = BOUNDARY[age]
        boundary_pass = (
            mapping.get((event_id, boundary_species), {}).get("mapping_status")
            == "PASS"
        )
        expected_postdup = POSTDUP[age]
        postdup_pass = [
            species
            for species in expected_postdup
            if mapping.get((event_id, species), {}).get("mapping_status") == "PASS"
        ]
        single_expected = [
            species
            for species in SPECIES
            if event[f"{species}_state"] in {"1", "2"}
        ]
        single_pass = [
            species
            for species in single_expected
            if mapping.get((event_id, species), {}).get("mapping_status") == "PASS"
        ]
        strict = boundary_pass and len(postdup_pass) == len(expected_postdup)
        relaxed = boundary_pass and (
            not expected_postdup or len(postdup_pass) >= 1
        )
        output.append(
            {
                "event_id": event_id,
                "age_bin": age,
                "boundary_species": boundary_species,
                "boundary_single_copy_mapping": (
                    "PASS" if boundary_pass else "FAIL"
                ),
                "expected_postduplication_species": ",".join(expected_postdup),
                "mapped_postduplication_species": ",".join(postdup_pass),
                "postduplication_species_expected_count": len(expected_postdup),
                "postduplication_species_mapped_count": len(postdup_pass),
                "single_copy_species_expected": ",".join(single_expected),
                "single_copy_species_mapped": ",".join(single_pass),
                "single_copy_species_expected_count": len(single_expected),
                "single_copy_species_mapped_count": len(single_pass),
                "strict_mapping_ready": "PASS" if strict else "FAIL",
                "relaxed_mapping_ready": "PASS" if relaxed else "FAIL",
            }
        )
    write_tsv(
        args.pilot / "outgroup_mapping/event_mapping_readiness.tsv", output
    )

    summary = []
    for age in ("time1", "time2", "time3"):
        subset = [row for row in output if row["age_bin"] == age]
        counts = Counter(
            (
                row["boundary_single_copy_mapping"],
                row["strict_mapping_ready"],
                row["relaxed_mapping_ready"],
            )
            for row in subset
        )
        summary.append(
            {
                "age_bin": age,
                "core_eligible_events": len(subset),
                "boundary_mapping_pass": sum(
                    row["boundary_single_copy_mapping"] == "PASS"
                    for row in subset
                ),
                "strict_mapping_ready": sum(
                    row["strict_mapping_ready"] == "PASS" for row in subset
                ),
                "relaxed_mapping_ready": sum(
                    row["relaxed_mapping_ready"] == "PASS" for row in subset
                ),
            }
        )
    write_tsv(
        args.pilot / "outgroup_mapping/event_mapping_readiness_summary.tsv",
        summary,
    )
    print(
        f"Classified mapping readiness for {len(output)} core-eligible events"
    )


if __name__ == "__main__":
    main()
