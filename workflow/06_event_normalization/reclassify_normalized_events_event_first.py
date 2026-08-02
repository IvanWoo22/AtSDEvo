#!/usr/bin/env python3
"""Classify normalized two-locus SD events against primary synteny nodes.

This workflow deliberately separates event construction from age/P-D
classification.  It ignores all call-level age labels while reconstructing
canonical physical loci A/B, then scores synteny on an event representative
core.  Call-level classifications are retained only as a QC comparison after
the event result has been assigned.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
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


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def interval_overlap(
    start: int, end: int, intervals: list[tuple[int, int]]
) -> int:
    total = 0
    for other_start, other_end in intervals:
        if other_start >= end:
            break
        if other_end > start:
            total += max(0, min(end, other_end) - max(start, other_start))
    return total


def overlapping_blocks(
    chrom: str,
    start: int,
    end: int,
    blocks: dict[str, list[tuple[int, int, int]]],
) -> list[int]:
    return sorted(
        {
            block_id
            for block_start, block_end, block_id in blocks.get(chrom, [])
            if block_start < end and block_end > start
        }
    )


def strict_classify(states: tuple[int, ...]) -> tuple[str, str, str]:
    if not any(states):
        return "EXCLUDE", "all_zero", "NA"
    prefix_both = 0
    while prefix_both < len(states) and states[prefix_both] == 3:
        prefix_both += 1
    if any(state == 3 for state in states[prefix_both:]):
        return "EXCLUDE", "nonmonotonic_both_reappears", "NA"
    if prefix_both == len(states):
        return "EXCLUDE", "older_than_N4_unpolarized", "NA"
    boundary = states[prefix_both]
    if boundary == 0:
        return "EXCLUDE", "boundary_node_uninformative", "NA"
    if boundary not in (1, 2):
        return "EXCLUDE", "invalid_boundary_state", "NA"
    if any(state not in (0, boundary) for state in states[prefix_both + 1 :]):
        return "EXCLUDE", "older_nodes_copy_conflict", "NA"
    return (
        "PASS",
        f"time{prefix_both + 1}",
        "locus_A" if boundary == 1 else "locus_B",
    )


def representative_call(calls: list[dict[str, str]]) -> dict[str, str]:
    return max(
        calls,
        key=lambda row: (
            min(
                int(row["biser_max_mate_length_bp"]),
                int(row["biser_alignment_span_bp"]),
            ),
            -float(row["biser_error"]),
            -int(row["call_id"]),
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--primary-node-analysis",
        type=Path,
        help="所选 BISER 分支对应的 call-level 共线性证据目录。",
    )
    parser.add_argument("--minimum-length", type=int, default=1000)
    parser.add_argument("--thresholds", default="0.25,0.5,0.75")
    args = parser.parse_args()

    old = args.primary_node_analysis or (
        args.project / "06_sd_age_tracing_preparation/primary_node_analysis"
    )
    network = args.output / "network_age_blind"
    membership = read_tsv(network / "event_call_membership.tsv")
    call_rows = read_tsv(old / "call_evidence.threshold_0.5.tsv")
    calls = {row["call_id"]: row for row in call_rows}
    network_events = {
        row["event_id"]: row
        for row in read_tsv(network / "strict_two_copy_events.tsv")
    }

    event_members: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in membership:
        event_members[row["event_id"]].append(row)

    projected: dict[str, dict[str, list[tuple[int, int]]]] = {}
    block_index: dict[str, dict[str, list[tuple[int, int, int]]]] = {}
    for species in SPECIES:
        projected[species] = defaultdict(list)
        bed = old / f"projected_intervals/{species}.TAIR12_collinear.bed"
        with bed.open() as handle:
            for line in handle:
                chrom, start, end = line.rstrip().split("\t")
                projected[species][chrom].append((int(start), int(end)))
        for chrom in projected[species]:
            projected[species][chrom] = merge_intervals(
                projected[species][chrom]
            )
        block_index[species] = defaultdict(list)
        for row in read_tsv(old / f"blocks/{species}.projected_blocks.tsv"):
            block_index[species][row["chromosome"]].append(
                (int(row["start"]), int(row["end"]), int(row["block_id"]))
            )
        for chrom in block_index[species]:
            block_index[species][chrom].sort()

    catalog = []
    arm_map_rows = []
    event_intervals: dict[
        str, dict[str, dict[str, tuple[str, int, int]]]
    ] = {}
    representative_by_event: dict[str, dict[str, str]] = {}
    qc = Counter()

    for event_id in sorted(event_members):
        members = sorted(
            event_members[event_id], key=lambda row: int(row["call_id"])
        )
        source_calls = [calls[row["call_id"]] for row in members]
        arms_by_locus: dict[str, list[tuple[str, int, int, str, int]]] = (
            defaultdict(list)
        )
        member_by_call = {row["call_id"]: row for row in members}
        for member in members:
            call = calls[member["call_id"]]
            for copy in (1, 2):
                locus_id = member[f"locus80_copy{copy}"]
                arms_by_locus[locus_id].append(
                    (
                        call[f"copy{copy}_chrom"],
                        int(call[f"copy{copy}_start"]),
                        int(call[f"copy{copy}_end"]),
                        call["call_id"],
                        copy,
                    )
                )
        if len(arms_by_locus) != 2:
            raise SystemExit(
                f"{event_id}: expected two normalized loci, "
                f"observed {len(arms_by_locus)}"
            )

        physical = []
        for locus_id, arms in arms_by_locus.items():
            chroms = {arm[0] for arm in arms}
            if len(chroms) != 1:
                raise SystemExit(f"{event_id}/{locus_id}: multi-chrom locus")
            chrom = next(iter(chroms))
            envelope_start = min(arm[1] for arm in arms)
            envelope_end = max(arm[2] for arm in arms)
            core_start = max(arm[1] for arm in arms)
            core_end = min(arm[2] for arm in arms)
            if core_start >= core_end:
                raise SystemExit(f"{event_id}/{locus_id}: empty common core")
            physical.append(
                (
                    (chrom, envelope_start, envelope_end, int(locus_id)),
                    locus_id,
                    chrom,
                    envelope_start,
                    envelope_end,
                    core_start,
                    core_end,
                )
            )
        physical.sort()
        locus_for_label = {
            "locus_A": physical[0],
            "locus_B": physical[1],
        }
        label_for_locus = {
            value[1]: label for label, value in locus_for_label.items()
        }

        representative = representative_call(source_calls)
        representative_by_event[event_id] = representative
        representative_member = member_by_call[representative["call_id"]]
        representative_intervals = {}
        for copy in (1, 2):
            locus_id = representative_member[f"locus80_copy{copy}"]
            label = label_for_locus[locus_id]
            representative_intervals[label] = (
                representative[f"copy{copy}_chrom"],
                int(representative[f"copy{copy}_start"]),
                int(representative[f"copy{copy}_end"]),
            )

        intervals = {}
        for label, value in locus_for_label.items():
            _, locus_id, chrom, env_start, env_end, core_start, core_end = value
            intervals[label] = {
                "representative": representative_intervals[label],
                "common_core": (chrom, core_start, core_end),
                "envelope": (chrom, env_start, env_end),
            }
        event_intervals[event_id] = intervals

        for member in members:
            call = calls[member["call_id"]]
            for copy in (1, 2):
                locus_id = member[f"locus80_copy{copy}"]
                arm_map_rows.append(
                    {
                        "event_id": event_id,
                        "call_id": call["call_id"],
                        "source_copy_label": f"copy{copy}",
                        "locus80_id": locus_id,
                        "canonical_event_locus": label_for_locus[locus_id],
                        "chromosome": call[f"copy{copy}_chrom"],
                        "start": call[f"copy{copy}_start"],
                        "end": call[f"copy{copy}_end"],
                    }
                )

        row: dict[str, object] = {
            "event_id": event_id,
            "source_call_count": len(source_calls),
            "source_call_ids": ",".join(call["call_id"] for call in source_calls),
            "representative_call_id": representative["call_id"],
            "relative_orientation": network_events[event_id][
                "relative_orientation"
            ],
            "network_stable_at_reciprocal_overlap_0.5_and_0.8": "PASS",
        }
        for label in ("locus_A", "locus_B"):
            value = locus_for_label[label]
            _, locus_id, chrom, env_start, env_end, core_start, core_end = value
            rep_chrom, rep_start, rep_end = intervals[label]["representative"]
            prefix = label
            row.update(
                {
                    f"{prefix}_network_id": locus_id,
                    f"{prefix}_chrom": chrom,
                    f"{prefix}_representative_start": rep_start,
                    f"{prefix}_representative_end": rep_end,
                    f"{prefix}_representative_bp": rep_end - rep_start,
                    f"{prefix}_common_core_start": core_start,
                    f"{prefix}_common_core_end": core_end,
                    f"{prefix}_common_core_bp": core_end - core_start,
                    f"{prefix}_envelope_start": env_start,
                    f"{prefix}_envelope_end": env_end,
                    f"{prefix}_envelope_bp": env_end - env_start,
                }
            )
            if rep_chrom != chrom:
                raise SystemExit(f"{event_id}/{label}: representative mismatch")
        row.update(
            {
                "representative_max_mate_length_bp": representative[
                    "biser_max_mate_length_bp"
                ],
                "representative_alignment_span_bp": representative[
                    "biser_alignment_span_bp"
                ],
            }
        )
        catalog.append(row)
        qc["normalized_events"] += 1
        qc[f"events_with_{len(source_calls)}_source_calls"] += 1

    thresholds = [float(value) for value in args.thresholds.split(",")]
    scope_classifications: dict[
        tuple[str, str, str], tuple[str, str, str, str]
    ] = {}
    detail_by_main: dict[str, dict[str, object]] = {}
    classification_rows = []

    for scope in ("representative", "common_core", "envelope"):
        for threshold in thresholds:
            for event in catalog:
                event_id = str(event["event_id"])
                states = []
                detail: dict[str, object] = {}
                for species in SPECIES:
                    present = []
                    for index, label in enumerate(("locus_A", "locus_B"), 1):
                        chrom, start, end = event_intervals[event_id][label][scope]
                        overlap_bp = interval_overlap(
                            start, end, projected[species].get(chrom, [])
                        )
                        fraction = overlap_bp / max(1, end - start)
                        detail[f"{species}_{label}_overlap_bp"] = overlap_bp
                        detail[f"{species}_{label}_overlap_fraction"] = (
                            f"{fraction:.6f}"
                        )
                        detail[f"{species}_{label}_block_ids"] = ",".join(
                            map(
                                str,
                                overlapping_blocks(
                                    chrom,
                                    start,
                                    end,
                                    block_index[species],
                                ),
                            )
                        )
                        present.append(fraction >= threshold)
                    state = (
                        3
                        if all(present)
                        else 1
                        if present[0]
                        else 2
                        if present[1]
                        else 0
                    )
                    states.append(state)
                    detail[f"{species}_state"] = state
                status, age_bin, p_locus = strict_classify(tuple(states))
                pattern = "".join(map(str, states))
                scope_classifications[(event_id, scope, str(threshold))] = (
                    status,
                    age_bin,
                    p_locus,
                    pattern,
                )
                classification_rows.append(
                    {
                        "event_id": event_id,
                        "interval_scope": scope,
                        "overlap_threshold": threshold,
                        "primary_pattern": pattern,
                        "strict_status_before_length": status,
                        "strict_age_bin_before_length": age_bin,
                        "provisional_p_locus_before_length": p_locus,
                    }
                )
                if scope == "representative" and threshold == 0.5:
                    detail_by_main[event_id] = detail

    event_rows = []
    for event in catalog:
        event_id = str(event["event_id"])
        status, age_bin, p_locus, pattern = scope_classifications[
            (event_id, "representative", "0.5")
        ]
        length_pass = (
            int(event["locus_A_representative_bp"]) >= args.minimum_length
            and int(event["locus_B_representative_bp"]) >= args.minimum_length
            and int(event["representative_max_mate_length_bp"])
            >= args.minimum_length
            and int(event["representative_alignment_span_bp"])
            >= args.minimum_length
        )
        # Length and callable-sequence metrics are descriptive QC only.
        # Do not overwrite a successful synteny-derived age/P-D assignment:
        # short events are retained for downstream local-comparability analysis.
        d_locus = (
            "locus_B"
            if p_locus == "locus_A"
            else "locus_A"
            if p_locus == "locus_B"
            else "NA"
        )
        stable = True
        main_without_length = scope_classifications[
            (event_id, "representative", "0.5")
        ][:3]
        for threshold in thresholds:
            if (
                scope_classifications[
                    (event_id, "representative", str(threshold))
                ][:3]
                != main_without_length
            ):
                stable = False
        scope_agreement = all(
            scope_classifications[(event_id, scope, "0.5")][:3]
            == main_without_length
            for scope in ("common_core", "envelope")
        )
        row = {
            **event,
            **detail_by_main[event_id],
            "primary_pattern": pattern,
            "strict_minimum_1kb_core": "PASS" if length_pass else "FAIL",
            "minimum_1kb_filter_policy": "DIAGNOSTIC_NOT_GATING",
            "strict_pd_status": status,
            "strict_age_bin": age_bin,
            "provisional_p_locus": p_locus,
            "provisional_d_locus": d_locus,
            "classification_stable_at_overlap_0.25_0.5_0.75": (
                "PASS" if stable else "FAIL"
            ),
            "classification_agrees_representative_core_envelope": (
                "PASS" if scope_agreement else "FAIL"
            ),
        }
        event_rows.append(row)

    source_qc_rows = []
    false_label_swaps = 0
    for event in event_rows:
        event_id = str(event["event_id"])
        members = event_members[event_id]
        physical_results = []
        call_label_results = []
        for member in members:
            call = calls[member["call_id"]]
            call_status = call["strict_status"]
            call_age = call["strict_age_bin"]
            call_p = call["provisional_p_copy"]
            if call_p in ("copy1", "copy2"):
                copy = 1 if call_p == "copy1" else 2
                network_id = member[f"locus80_copy{copy}"]
                event_p = next(
                    label
                    for label in ("locus_A", "locus_B")
                    if str(event[f"{label}_network_id"]) == network_id
                )
            else:
                event_p = "NA"
            physical_results.append((call_status, call_age, event_p))
            call_label_results.append((call_status, call_age, call_p))
            source_qc_rows.append(
                {
                    "event_id": event_id,
                    "call_id": call["call_id"],
                    "call_level_status": call_status,
                    "call_level_age_bin": call_age,
                    "call_level_p_copy_label": call_p,
                    "call_level_p_physical_event_locus": event_p,
                    "event_first_status": event["strict_pd_status"],
                    "event_first_age_bin": event["strict_age_bin"],
                    "event_first_p_locus": event["provisional_p_locus"],
                }
            )
        if len(set(call_label_results)) > 1 and len(set(physical_results)) == 1:
            false_label_swaps += 1

    disposition = Counter(
        (row["strict_pd_status"], row["strict_age_bin"]) for row in event_rows
    )
    pass_events = [
        row for row in event_rows if row["strict_pd_status"] == "PASS"
    ]
    summary = [
        {"metric": "input_biser_calls", "value": len(call_rows)},
        {"metric": "network_stable_source_calls", "value": len(membership)},
        {"metric": "normalized_two_locus_events", "value": len(catalog)},
        {
            "metric": "events_passing_event_first_age_PD_without_length_gating",
            "value": len(pass_events),
        },
        {
            "metric": "pass_events_stable_at_overlap_0.25_0.5_0.75",
            "value": sum(
                row["strict_pd_status"] == "PASS"
                and row[
                    "classification_stable_at_overlap_0.25_0.5_0.75"
                ]
                == "PASS"
                for row in event_rows
            ),
        },
        {
            "metric": "pass_events_agree_across_interval_scopes",
            "value": sum(
                row["strict_pd_status"] == "PASS"
                and row[
                    "classification_agrees_representative_core_envelope"
                ]
                == "PASS"
                for row in event_rows
            ),
        },
        {
            "metric": "events_with_false_source_copy_label_disagreement",
            "value": false_label_swaps,
        },
    ]
    for (status, age_bin), count in sorted(disposition.items()):
        summary.append(
            {
                "metric": f"event_disposition_{status}_{age_bin}",
                "value": count,
            }
        )

    write_tsv(args.output / "events/event_first_catalog.tsv", catalog)
    write_tsv(args.output / "events/event_first_events.tsv", event_rows)
    write_tsv(args.output / "events/source_arm_to_event_locus.tsv", arm_map_rows)
    write_tsv(
        args.output / "statistics/event_classification_by_scope_threshold.tsv",
        classification_rows,
    )
    write_tsv(
        args.output / "statistics/source_call_event_concordance.tsv",
        source_qc_rows,
    )
    write_tsv(args.output / "statistics/event_first_summary.tsv", summary)
    print(
        f"Classified {len(catalog)} normalized events; "
        f"{len(pass_events)} pass age/P-D; {args.minimum_length}-bp length is "
        "reported as diagnostic QC only"
    )


if __name__ == "__main__":
    main()
