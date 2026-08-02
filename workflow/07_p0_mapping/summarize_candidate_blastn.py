#!/usr/bin/env python3
"""Summarize event-matched sensitive BLASTN hits in synteny candidate regions."""

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


def newly_covered(start: int, end: int, accepted: list[tuple[int, int]]) -> int:
    covered = 0
    cursor = start
    for other_start, other_end in sorted(accepted):
        if other_end <= cursor:
            continue
        if other_start >= end:
            break
        if other_start > cursor:
            covered += other_start - cursor
        cursor = max(cursor, other_end)
        if cursor >= end:
            break
    if cursor < end:
        covered += end - cursor
    return covered


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
        hits = args.pilot / f"outgroup_mapping/blastn/{code}.hits.tsv"
        with hits.open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                query_event = fields[0].split("|", 1)[0]
                target_event = fields[1].split("|", 1)[0]
                if query_event == target_event and fields[1] in candidates:
                    groups[(fields[0], fields[1])].append(fields)

    mapping_rows = []
    for (query, target), hsps in groups.items():
        query_length = int(hsps[0][2])
        candidate = candidates[target]
        candidate_start = int(candidate["candidate_start"])
        accepted: list[tuple[int, int]] = []
        covered = 0
        estimated_matches = 0.0
        accepted_hsps = 0
        target_starts = []
        target_ends = []
        for hsp in sorted(hsps, key=lambda row: float(row[12]), reverse=True):
            start, end = sorted((int(hsp[7]) - 1, int(hsp[8])))
            new_bp = newly_covered(start, end, accepted)
            if not new_bp:
                continue
            accepted.append((start, end))
            covered += new_bp
            estimated_matches += new_bp * float(hsp[4]) / 100
            accepted_hsps += 1
            subject_start, subject_end = sorted((int(hsp[9]), int(hsp[10])))
            target_starts.append(candidate_start + subject_start - 1)
            target_ends.append(candidate_start + subject_end)
        mapping_rows.append(
            {
                "event_id": candidate["event_id"],
                "age_bin": candidate["age_bin"],
                "species": candidate["species"],
                "species_state": candidate["species_state"],
                "query_role": query.split("|")[1],
                "target_TAIR12_copy": candidate["TAIR12_copy"],
                "target_copy_role": candidate["TAIR12_copy_role"],
                "two_sided_ordered_anchor_status": candidate[
                    "two_sided_ordered_anchor_status"
                ],
                "candidate_id": target,
                "outgroup_scaffold": candidate["outgroup_scaffold"],
                "candidate_start": candidate["candidate_start"],
                "candidate_end": candidate["candidate_end"],
                "mapped_genomic_start": min(target_starts),
                "mapped_genomic_end": max(target_ends),
                "query_length": query_length,
                "query_covered_bp": covered,
                "query_coverage": f"{covered / query_length:.6f}",
                "nonoverlap_weighted_identity": (
                    f"{estimated_matches / covered:.6f}" if covered else "0"
                ),
                "accepted_nonoverlap_hsps": accepted_hsps,
                "best_evalue": min(float(row[11]) for row in hsps),
                "best_bitscore": max(float(row[12]) for row in hsps),
            }
        )
    mapping_rows.sort(
        key=lambda row: (
            row["event_id"],
            row["species"],
            row["query_role"],
            -float(row["query_coverage"]),
            -float(row["nonoverlap_weighted_identity"]),
        )
    )
    write_tsv(
        args.pilot / "outgroup_mapping/event_matched_blastn_mappings.tsv",
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
            detail: dict[str, object] = {}
            passed_requirements = []
            best_by_label: dict[str, dict[str, object] | None] = {}
            for query_role, target_role in requirements:
                choices = by_key.get((event_id, code, query_role, target_role), [])
                best = max(
                    choices,
                    key=lambda row: (
                        row["two_sided_ordered_anchor_status"] == "PASS",
                        float(row["query_coverage"]),
                        float(row["nonoverlap_weighted_identity"]),
                    ),
                    default=None,
                )
                passed = (
                    best is not None
                    and best["two_sided_ordered_anchor_status"] == "PASS"
                    and float(best["query_coverage"]) >= args.min_query_coverage
                    and float(best["nonoverlap_weighted_identity"])
                    >= args.min_identity
                )
                passed_requirements.append(passed)
                label = f"{query_role}_to_{target_role}"
                best_by_label[label] = best
                detail[f"{label}_status"] = "PASS" if passed else "FAIL"
                detail[f"{label}_query_coverage"] = (
                    best["query_coverage"] if best else "NA"
                )
                detail[f"{label}_identity"] = (
                    best["nonoverlap_weighted_identity"] if best else "NA"
                )
                detail[f"{label}_candidate_id"] = (
                    best["candidate_id"] if best else "NA"
                )
                detail[f"{label}_two_sided_ordered_anchor_status"] = (
                    best["two_sided_ordered_anchor_status"] if best else "NA"
                )
            if state == 3:
                best_p = best_by_label["P_to_P"]
                best_d = best_by_label["D_to_D"]
                distinct = False
                if best_p is not None and best_d is not None:
                    distinct = (
                        best_p["outgroup_scaffold"] != best_d["outgroup_scaffold"]
                        or int(best_p["mapped_genomic_end"])
                        <= int(best_d["mapped_genomic_start"])
                        or int(best_d["mapped_genomic_end"])
                        <= int(best_p["mapped_genomic_start"])
                    )
                detail["state3_target_loci_distinct"] = (
                    "PASS" if distinct else "FAIL"
                )
                passed_requirements.append(distinct)
            else:
                detail["state3_target_loci_distinct"] = "NA"
            summary_rows.append(
                {
                    "event_id": event_id,
                    "age_bin": event["strict_age_bin"],
                    "species": code,
                    "species_state": state,
                    "expected_locus_class": (
                        "postduplication_two_copy"
                        if state == 3
                        else "preduplication_single_copy"
                    ),
                    "mapping_status": (
                        "PASS" if all(passed_requirements) else "FAIL"
                    ),
                    "min_query_coverage_rule": args.min_query_coverage,
                    "min_identity_rule": args.min_identity,
                    **detail,
                }
            )
    write_tsv(
        args.pilot / "outgroup_mapping/event_species_blastn_summary.tsv",
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
        args.pilot / "outgroup_mapping/blastn_mapping_aggregate_summary.tsv",
        aggregate,
    )
    print(
        f"Summarized {len(mapping_rows)} event-matched BLASTN mappings across "
        f"{len(summary_rows)} event-species tests"
    )


if __name__ == "__main__":
    main()
