#!/usr/bin/env python3
"""Build the inclusive, length-ungated queue with at least one mapped P0."""

import csv
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
BOUNDARY_INDEX = {"time1": 0, "time2": 1, "time3": 2}


def read(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


root = Path(__file__).resolve().parent
events = read(root / "inputs/high_priority_core_eligible.tsv")
qc = {row["event_id"]: row for row in read(root / "core/core_event_qc.tsv")}
mapping = {
    (row["event_id"], row["species"]): row
    for row in read(root / "outgroup_mapping/event_species_blastn_summary.tsv")
}
rows = []
for event in events:
    event_id, age = event["event_id"], event["strict_age_bin"]
    mapped_p0 = [
        species
        for species in SPECIES[BOUNDARY_INDEX[age] :]
        if event[f"{species}_state"] in {"1", "2"}
        and mapping.get((event_id, species), {}).get("mapping_status") == "PASS"
    ]
    if not mapped_p0:
        continue
    boundary = mapped_p0[0]
    mapped_postdup = [
        species
        for species in SPECIES[: BOUNDARY_INDEX[age]]
        if event[f"{species}_state"] == "3"
        and mapping.get((event_id, species), {}).get("mapping_status") == "PASS"
    ]
    rows.append(
        {
            "event_id": event_id,
            "age_bin": age,
            "analysis_tier": (
                "EXPECTED_BOUNDARY_P0"
                if SPECIES[BOUNDARY_INDEX[age]] == boundary
                else "DEEPER_P0_FALLBACK"
            ),
            "boundary_P0_species": boundary,
            "mapped_P0_species": ",".join(mapped_p0),
            "mapped_P0_species_count": len(mapped_p0),
            "mapped_postduplication_species": ",".join(mapped_postdup),
            "mapped_postduplication_species_count": len(mapped_postdup),
            "paired_M_core_bp": qc[event_id]["paired_M_core_bp"],
            "jointly_callable_uppercase_acgt_bp": qc[event_id][
                "jointly_callable_uppercase_acgt_bp"
            ],
            "present_day_PD_mismatch_pct_callable": qc[event_id][
                "present_day_PD_mismatch_pct_callable"
            ],
            "legacy_1kb_filters_policy": "DIAGNOSTIC_NOT_GATING",
        }
    )
write(root / "pilot/inclusive_p0_mapping_queue.tsv", rows)
print(f"inclusive_events={len(events)} mapped_P0_queue={len(rows)}")
