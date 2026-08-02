#!/usr/bin/env python3
"""Summarize event-matched minimap2 mappings to synteny candidate regions."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")


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


def covered_bp(intervals: list[tuple[int, int]]) -> int:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", required=True, type=Path)
    parser.add_argument("--min-query-coverage", type=float, default=0.5)
    parser.add_argument("--min-identity", type=float, default=0.6)
    args = parser.parse_args()

    events = {
        row["event_id"]: row
        for row in read_tsv(args.pilot / "inputs/high_priority_core_eligible.tsv")
    }
    candidates = {
        row["candidate_id"]: row
        for row in read_tsv(
            args.pilot / "outgroup_mapping/candidate_region_manifest.tsv"
        )
    }
    groups: dict[tuple[str, str], list[list[str]]] = defaultdict(list)
    for code in SPECIES:
        paf = args.pilot / f"outgroup_mapping/minimap2/{code}.paf"
        with paf.open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                query_event = fields[0].split("|", 1)[0]
                target_event = fields[5].split("|", 1)[0]
                if query_event == target_event and fields[5] in candidates:
                    groups[(fields[0], fields[5])].append(fields)

    mapping_rows = []
    for (query, target), alignments in groups.items():
        query_length = int(alignments[0][1])
        query_intervals = [(int(row[2]), int(row[3])) for row in alignments]
        query_covered = covered_bp(query_intervals)
        matches = sum(int(row[9]) for row in alignments)
        aligned_columns = sum(int(row[10]) for row in alignments)
        candidate = candidates[target]
        query_role = query.split("|")[1]
        mapping_rows.append(
            {
                "event_id": candidate["event_id"],
                "age_bin": candidate["age_bin"],
                "species": candidate["species"],
                "species_state": candidate["species_state"],
                "query_role": query_role,
                "target_TAIR12_copy": candidate["TAIR12_copy"],
                "target_copy_role": candidate["TAIR12_copy_role"],
                "candidate_id": target,
                "outgroup_scaffold": candidate["outgroup_scaffold"],
                "candidate_start": candidate["candidate_start"],
                "candidate_end": candidate["candidate_end"],
                "query_length": query_length,
                "query_covered_bp": query_covered,
                "query_coverage": f"{query_covered / query_length:.6f}",
                "aggregate_matches": matches,
                "aggregate_aligned_columns": aligned_columns,
                "aggregate_identity": (
                    f"{matches / aligned_columns:.6f}"
                    if aligned_columns
                    else "0"
                ),
                "max_mapq": max(int(row[11]) for row in alignments),
                "alignment_segments": len(alignments),
            }
        )
    mapping_rows.sort(
        key=lambda row: (
            row["event_id"],
            row["species"],
            row["query_role"],
            -float(row["query_coverage"]),
            -float(row["aggregate_identity"]),
        )
    )
    write_tsv(
        args.pilot / "outgroup_mapping/event_matched_minimap2_mappings.tsv",
        mapping_rows,
    )

    by_key: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in mapping_rows:
        by_key[
            (
                str(row["event_id"]),
                str(row["species"]),
                str(row["query_role"]),
                str(row["target_copy_role"]),
            )
        ].append(row)

    summary_rows = []
    for event_id, event in events.items():
        for code in SPECIES:
            state = int(event[f"{code}_state"])
            if state == 0:
                continue
            requirements = [("P", "P")]
            if state == 3:
                requirements.append(("D", "D"))
            results = []
            detail = {}
            for query_role, target_role in requirements:
                candidates_for_role = by_key.get(
                    (event_id, code, query_role, target_role), []
                )
                best = max(
                    candidates_for_role,
                    key=lambda row: (
                        float(row["query_coverage"]),
                        float(row["aggregate_identity"]),
                    ),
                    default=None,
                )
                passed = (
                    best is not None
                    and float(best["query_coverage"]) >= args.min_query_coverage
                    and float(best["aggregate_identity"]) >= args.min_identity
                )
                results.append(passed)
                label = f"{query_role}_to_{target_role}"
                detail[f"{label}_status"] = "PASS" if passed else "FAIL"
                detail[f"{label}_query_coverage"] = (
                    best["query_coverage"] if best else "NA"
                )
                detail[f"{label}_identity"] = (
                    best["aggregate_identity"] if best else "NA"
                )
                detail[f"{label}_candidate_id"] = (
                    best["candidate_id"] if best else "NA"
                )
            summary_rows.append(
                {
                    "event_id": event_id,
                    "age_bin": event["strict_age_bin"],
                    "species": code,
                    "species_state": state,
                    "expected_locus_class": (
                        "postduplication_two_copy" if state == 3 else "preduplication_single_copy"
                    ),
                    "mapping_status": "PASS" if all(results) else "FAIL",
                    "min_query_coverage_rule": args.min_query_coverage,
                    "min_identity_rule": args.min_identity,
                    **detail,
                }
            )
    write_tsv(
        args.pilot / "outgroup_mapping/event_species_mapping_summary.tsv",
        summary_rows,
    )

    aggregate = []
    for code in SPECIES:
        for expected in ("preduplication_single_copy", "postduplication_two_copy"):
            subset = [
                row
                for row in summary_rows
                if row["species"] == code and row["expected_locus_class"] == expected
            ]
            if subset:
                pass_count = sum(
                    row["mapping_status"] == "PASS" for row in subset
                )
                aggregate.append(
                    {
                        "species": code,
                        "expected_locus_class": expected,
                        "event_species_tests": len(subset),
                        "mapping_pass": pass_count,
                        "mapping_pass_pct": (
                            f"{100 * pass_count / len(subset):.4f}"
                        ),
                    }
                )
    write_tsv(
        args.pilot / "outgroup_mapping/mapping_aggregate_summary.tsv",
        aggregate,
    )
    print(
        f"Summarized {len(mapping_rows)} event-matched mappings across "
        f"{len(summary_rows)} event-species tests"
    )


if __name__ == "__main__":
    main()
