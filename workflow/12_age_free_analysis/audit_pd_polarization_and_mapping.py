#!/usr/bin/env python3
"""Audit age-free P/D polarization, attrition, coverage, and P0 mapping loss."""

from __future__ import annotations

import csv
import os
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(os.environ.get("SD_PROJECT_ROOT", Path(__file__).resolve().parents[2]))
OUT = Path(
    os.environ.get(
        "SD_REASSESSMENT_ROOT", ROOT / "14_pd_polarization_reassessment"
    )
)
SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
BOUNDARY_INDEX = {"time1": 0, "time2": 1, "time3": 2}
GENOME_BP = 142_481_245


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


def orient(pattern: str) -> tuple[str, str, int]:
    singletons = [state for state in pattern if state in "12"]
    distinct = set(singletons)
    if not singletons:
        return "UNINFORMATIVE_NO_SINGLE_COPY", "NA", 0
    if len(distinct) > 1:
        return "CONFLICT_BOTH_SINGLE_COPY_LOCI", "NA", len(singletons)
    state = singletons[0]
    return (
        "POLARIZABLE_STRONG" if len(singletons) >= 2 else "POLARIZABLE_SINGLE_NODE",
        "locus_A" if state == "1" else "locus_B",
        len(singletons),
    )


def interval_metrics(events: list[dict[str, str]]) -> tuple[int, int]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    summed = 0
    for row in events:
        for locus in ("A", "B"):
            chrom = row[f"locus_{locus}_chrom"]
            start = int(row[f"locus_{locus}_representative_start"])
            end = int(row[f"locus_{locus}_representative_end"])
            intervals[chrom].append((start, end))
            summed += end - start
    union = 0
    for values in intervals.values():
        values.sort()
        current_start, current_end = values[0]
        for start, end in values[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                union += current_end - current_start
                current_start, current_end = start, end
        union += current_end - current_start
    return summed, union


def biser_interval_metrics(path: Path, wanted: set[int] | None = None) -> tuple[int, int]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    summed = 0
    with path.open() as handle:
        for call_id, line in enumerate(handle, 1):
            if wanted is not None and call_id not in wanted:
                continue
            fields = line.rstrip().split("\t")
            for chrom, start, end in (
                (fields[0], int(fields[1]), int(fields[2])),
                (fields[3], int(fields[4]), int(fields[5])),
            ):
                intervals[chrom].append((start, end))
                summed += end - start
    union = 0
    for values in intervals.values():
        values.sort()
        left, right = values[0]
        for start, end in values[1:]:
            if start <= right:
                right = max(right, end)
            else:
                union += right - left
                left, right = start, end
        union += right - left
    return summed, union


def mapping_failure(row: dict[str, str]) -> str:
    if row.get("P_to_P_candidate_id", "NA") in {"", "NA"}:
        return "NO_SEQUENCE_CANDIDATE"
    if row.get("P_to_P_two_sided_ordered_anchor_status") != "PASS":
        return "TWO_SIDED_ORDERED_ANCHOR_FAIL"
    qcov = float(row.get("P_to_P_query_coverage", "nan"))
    identity = float(row.get("P_to_P_identity", "nan"))
    if qcov < float(row["min_query_coverage_rule"]):
        return "QUERY_COVERAGE_FAIL"
    if identity < float(row["min_identity_rule"]):
        return "IDENTITY_FAIL"
    return "OTHER_OR_COMPOSITE_FAIL"


def main() -> None:
    event_root = ROOT / "06_sd_age_tracing_preparation/event_first_reanalysis"
    all_events = read(event_root / "events/event_first_events.tsv")
    by_id = {row["event_id"]: row for row in all_events}
    matrix = read(event_root / "statistics/event_classification_by_scope_threshold.tsv")
    patterns: dict[str, dict[tuple[str, str], str]] = defaultdict(dict)
    for row in matrix:
        patterns[row["event_id"]][
            (row["interval_scope"], row["overlap_threshold"])
        ] = row["primary_pattern"]

    audit_rows = []
    for event in all_events:
        event_id = event["event_id"]
        primary_status, primary_p, single_nodes = orient(event["primary_pattern"])
        rep_orientations = [
            orient(patterns[event_id][("representative", threshold)])
            for threshold in ("0.25", "0.5", "0.75")
        ]
        all_orientations = [
            orient(pattern)
            for pattern in patterns[event_id].values()
        ]
        rep_stable = (
            all(value[1] == primary_p and value[1] != "NA" for value in rep_orientations)
        )
        scope_stable = (
            len(all_orientations) == 9
            and all(value[1] == primary_p and value[1] != "NA" for value in all_orientations)
        )
        audit_rows.append(
            {
                "event_id": event_id,
                "primary_pattern": event["primary_pattern"],
                "age_free_PD_status": primary_status,
                "age_free_P_locus": primary_p,
                "single_copy_support_nodes": single_nodes,
                "representative_threshold_stable": "PASS" if rep_stable else "FAIL",
                "scope_and_threshold_stable": "PASS" if scope_stable else "FAIL",
                "former_age_PD_status": event["strict_pd_status"],
                "former_age_bin": event["strict_age_bin"],
                "former_P_locus": event["provisional_p_locus"],
                "locus_A_chrom": event["locus_A_chrom"],
                "locus_A_start": event["locus_A_representative_start"],
                "locus_A_end": event["locus_A_representative_end"],
                "locus_B_chrom": event["locus_B_chrom"],
                "locus_B_start": event["locus_B_representative_start"],
                "locus_B_end": event["locus_B_representative_end"],
            }
        )
    write(OUT / "events/age_free_pd_polarization_audit.tsv", audit_rows)

    summary = []
    for label, selected in (
        ("all_normalized_two_copy_events", audit_rows),
        (
            "primary_pattern_age_free_PD",
            [row for row in audit_rows if str(row["age_free_PD_status"]).startswith("POLARIZABLE")],
        ),
        (
            "primary_pattern_strong_ge2_single_copy_nodes",
            [row for row in audit_rows if row["age_free_PD_status"] == "POLARIZABLE_STRONG"],
        ),
        (
            "representative_threshold_stable_age_free_PD",
            [row for row in audit_rows if row["representative_threshold_stable"] == "PASS"],
        ),
        (
            "scope_threshold_stable_age_free_PD",
            [row for row in audit_rows if row["scope_and_threshold_stable"] == "PASS"],
        ),
        (
            "former_age_PD_pass",
            [row for row in audit_rows if row["former_age_PD_status"] == "PASS"],
        ),
    ):
        original = [by_id[row["event_id"]] for row in selected]
        summed, union = interval_metrics(original)
        summary.append(
            {
                "candidate_definition": label,
                "events": len(selected),
                "two_copy_interval_sum_bp": summed,
                "genome_union_covered_bp": union,
                "genome_union_coverage_pct": f"{100 * union / GENOME_BP:.4f}",
            }
        )
    write(OUT / "statistics/age_free_pd_summary.tsv", summary)
    rescued = Counter()
    for row in audit_rows:
        if (
            row["scope_and_threshold_stable"] == "PASS"
            and row["former_age_PD_status"] != "PASS"
        ):
            rescued[str(row["former_age_bin"])] += 1
    write(
        OUT / "statistics/age_free_pd_rescued_from_former_exclusions.tsv",
        [
            {"former_exclusion_reason": reason, "rescued_events": count}
            for reason, count in rescued.most_common()
        ],
    )

    inclusive_root = ROOT / "12_inclusive_pd_sequence_variation"
    inclusive = read(inclusive_root / "inputs/high_priority_events.tsv")
    queue = read(inclusive_root / "pilot/inclusive_p0_mapping_queue.tsv")
    atom_events = {
        row["event_id"]
        for row in read(inclusive_root / "microindel_local_msa/atomic_region_manifest.tsv")
    }
    snp20 = read(
        inclusive_root
        / "snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge20.event_metrics.tsv"
    )
    snp200 = read(
        inclusive_root
        / "snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge200.event_metrics.tsv"
    )
    indel_events = {
        row["event_id"]
        for row in read(inclusive_root / "microindel_local_msa/denovo_microindel_inference.tsv")
        if row["evidence_tier"] == "PRIMARY"
    }
    stable_calls = {
        int(row["call_id"])
        for row in read(event_root / "network_age_blind/call_network_tiers.tsv")
        if row["network_tier"] == "strict_stable_two_copy"
    }
    biser_path = (
        ROOT
        / "03_biser_segmental_duplication/runs/annotation_extended_softmask/biser_out"
    )
    all_call_sum, all_call_union = biser_interval_metrics(biser_path)
    stable_call_sum, stable_call_union = biser_interval_metrics(biser_path, stable_calls)

    funnel = [
        {
            "stage_order": 1,
            "branch": "common",
            "parent_stage": "START",
            "stage": "annotation_extended_BISER_calls",
            "units": "calls",
            "count": 4734,
            "two_copy_interval_sum_bp": all_call_sum,
            "genome_union_covered_bp": all_call_union,
            "genome_union_coverage_pct": f"{100 * all_call_union / GENOME_BP:.4f}",
            "retained_from_previous_pct": "100.000",
        },
        {
            "stage_order": 2,
            "branch": "common",
            "parent_stage": "annotation_extended_BISER_calls",
            "stage": "network_stable_source_calls",
            "units": "calls",
            "count": len(stable_calls),
            "two_copy_interval_sum_bp": stable_call_sum,
            "genome_union_covered_bp": stable_call_union,
            "genome_union_coverage_pct": f"{100 * stable_call_union / GENOME_BP:.4f}",
            "retained_from_previous_pct": f"{100 * len(stable_calls) / 4734:.3f}",
        },
    ]
    event_stages = [
        ("common", "network_stable_source_calls", "normalized_two_copy_events", all_events, len(stable_calls)),
        (
            "proposed_age_free",
            "normalized_two_copy_events",
            "age_free_PD_scope_threshold_stable",
            [
                by_id[row["event_id"]]
                for row in audit_rows
                if row["scope_and_threshold_stable"] == "PASS"
            ],
            len(all_events),
        ),
        ("current_branch", "normalized_two_copy_events", "current_time1_time3_input", inclusive, len(all_events)),
        ("current_branch", "current_time1_time3_input", "at_least_one_mapped_P0", [by_id[row["event_id"]] for row in queue], len(inclusive)),
        ("current_branch", "at_least_one_mapped_P0", "continuous_PD_P0_atoms", [by_id[event_id] for event_id in atom_events], len(queue)),
        ("current_branch", "continuous_PD_P0_atoms", "three_aligner_local_MSA_ge20", [by_id[row["event_id"]] for row in snp20], len(atom_events)),
        ("current_branch", "three_aligner_local_MSA_ge20", "three_aligner_local_MSA_ge200", [by_id[row["event_id"]] for row in snp200], len(snp20)),
        ("current_branch", "three_aligner_local_MSA_ge200", "at_least_one_primary_microindel", [by_id[event_id] for event_id in indel_events], len(snp200)),
    ]
    for order, (branch, parent, label, selected, previous) in enumerate(event_stages, 3):
        summed, union = interval_metrics(selected)
        funnel.append(
            {
                "stage_order": order,
                "branch": branch,
                "parent_stage": parent,
                "stage": label,
                "units": "events",
                "count": len(selected),
                "two_copy_interval_sum_bp": summed,
                "genome_union_covered_bp": union,
                "genome_union_coverage_pct": f"{100 * union / GENOME_BP:.4f}",
                "retained_from_previous_pct": f"{100 * len(selected) / previous:.3f}",
            }
        )
    write(OUT / "statistics/filter_funnel_with_genome_coverage.tsv", funnel)

    mapping_rows = read(inclusive_root / "outgroup_mapping/event_species_blastn_summary.tsv")
    mapping = {(row["event_id"], row["species"]): row for row in mapping_rows}
    candidate_manifest = read(
        inclusive_root / "outgroup_mapping/candidate_region_manifest.tsv"
    )
    candidate_event_species = {
        (row["event_id"], row["species"]) for row in candidate_manifest
    }
    queue_ids = {row["event_id"] for row in queue}
    failure_species = Counter()
    failure_event_profile = Counter()
    pattern_summary: dict[str, Counter[str]] = defaultdict(Counter)
    details = []
    for event in inclusive:
        event_id, age = event["event_id"], event["strict_age_bin"]
        candidate_species = [
            species
            for species in SPECIES[BOUNDARY_INDEX[age] :]
            if event[f"{species}_state"] in {"1", "2"}
        ]
        reasons = []
        passed = []
        for species in candidate_species:
            row = mapping.get((event_id, species))
            if row and row["mapping_status"] == "PASS":
                passed.append(species)
                continue
            reason = mapping_failure(row) if row else "NO_SPECIES_MAPPING_ROW"
            reasons.append(reason)
            failure_species[(species, reason)] += 1
        status = "PASS_ANY_P0" if event_id in queue_ids else "FAIL_ALL_P0"
        pattern_summary[event["primary_pattern"]][status] += 1
        if status == "FAIL_ALL_P0":
            profile = reasons[0] if len(set(reasons)) == 1 else "MIXED_FAILURES"
            failure_event_profile[profile] += 1
        details.append(
            {
                "event_id": event_id,
                "primary_pattern": event["primary_pattern"],
                "age_bin": age,
                "P_locus": event["provisional_p_locus"],
                "candidate_P0_species": ",".join(candidate_species),
                "passed_P0_species": ",".join(passed),
                "failed_P0_reasons": ",".join(reasons),
                "P0_queue_status": status,
            }
        )
    write(OUT / "mapping/current_362_P0_mapping_audit.tsv", details)
    failed_details = [row for row in details if row["P0_queue_status"] == "FAIL_ALL_P0"]
    failed_with_candidate = sum(
        any(
            (row["event_id"], species) in candidate_event_species
            for species in str(row["candidate_P0_species"]).split(",")
        )
        for row in failed_details
    )
    write(
        OUT / "statistics/P0_candidate_region_availability.tsv",
        [
            {
                "P0_failed_events": len(failed_details),
                "with_event_specific_candidate_region_for_any_P0_species": failed_with_candidate,
                "without_event_specific_candidate_region_for_all_P0_species": (
                    len(failed_details) - failed_with_candidate
                ),
            }
        ],
    )
    write(
        OUT / "statistics/P0_mapping_failure_event_profiles.tsv",
        [
            {"event_failure_profile": key, "events": value}
            for key, value in failure_event_profile.most_common()
        ],
    )
    write(
        OUT / "statistics/P0_mapping_failure_species_reasons.tsv",
        [
            {"species": key[0], "failure_reason": key[1], "event_species_pairs": value}
            for key, value in sorted(failure_species.items())
        ],
    )
    write(
        OUT / "statistics/P0_mapping_by_primary_pattern.tsv",
        [
            {
                "primary_pattern": pattern,
                "input_events": counts["PASS_ANY_P0"] + counts["FAIL_ALL_P0"],
                "events_with_any_P0": counts["PASS_ANY_P0"],
                "events_without_P0": counts["FAIL_ALL_P0"],
                "P0_success_pct": (
                    f"{100 * counts['PASS_ANY_P0'] / (counts['PASS_ANY_P0'] + counts['FAIL_ALL_P0']):.3f}"
                ),
            }
            for pattern, counts in sorted(
                pattern_summary.items(),
                key=lambda item: -(item[1]["PASS_ANY_P0"] + item[1]["FAIL_ALL_P0"]),
            )
        ],
    )


if __name__ == "__main__":
    main()
