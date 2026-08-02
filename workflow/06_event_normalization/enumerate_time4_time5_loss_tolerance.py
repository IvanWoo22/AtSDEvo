#!/usr/bin/env python3
"""Enumerate loss-tolerant time4/time5 patterns on normalized SD events."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


THRESHOLDS = ("0.25", "0.5", "0.75")
SCOPES = ("representative", "common_core", "envelope")


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def classify(pattern: str) -> tuple[str, str, str, str]:
    """Return age, P locus, loss class, and loss-node description."""
    if len(pattern) != 4:
        return "other", "NA", "other", "NA"
    a, b, d, c = pattern

    # time4: duplication predates D. strictus but postdates C. violacea.
    # B. stricta may retain both copies or show an isolated lineage loss.
    if a == "3" and b in "123" and d == "3" and c in "12":
        return (
            "time4_loss_tolerant",
            "locus_A" if c == "1" else "locus_B",
            "strict_no_loss" if b == "3" else "isolated_postdup_loss",
            "none" if b == "3" else "Bstricta",
        )

    # time5: duplication predates all four primary nodes. The two internal
    # lineages may independently retain A, B, or both; their retained sides
    # need not agree. P/D cannot be polarized with these four nodes.
    if a == "3" and b in "123" and d in "123" and c == "3":
        internal_losses = [
            node
            for node, state in (("Bstricta", b), ("Dstrictus", d))
            if state in "12"
        ]
        return (
            "time5_older_than_N4",
            "NA",
            (
                "strict_all_both"
                if not internal_losses
                else "one_internal_postdup_loss"
                if len(internal_losses) == 1
                else "two_internal_postdup_losses"
            ),
            ",".join(internal_losses) if internal_losses else "none",
        )
    return "other", "NA", "other", "NA"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    args = parser.parse_args()

    root = (
        args.project
        / "06_sd_age_tracing_preparation/event_first_reanalysis"
    )
    events = {
        row["event_id"]: row
        for row in read(root / "events/event_first_events.tsv")
    }
    effective_qc_path = (
        root / "events/event_effective_unmasked_length_qc.tsv"
    )
    effective_qc = (
        {
            row["event_id"]: row
            for row in read(effective_qc_path)
        }
        if effective_qc_path.exists()
        else {}
    )
    classification_rows = read(
        root / "statistics/event_classification_by_scope_threshold.tsv"
    )
    patterns = {
        (
            row["event_id"],
            row["interval_scope"],
            row["overlap_threshold"],
        ): row["primary_pattern"]
        for row in classification_rows
    }

    pattern_universe = []
    for b in "123":
        for c in "12":
            pattern_universe.append(
                {
                    "extended_age_bin": "time4_loss_tolerant",
                    "pattern": f"3{b}3{c}",
                    "design": (
                        "strict"
                        if b == "3"
                        else "Bstricta_isolated_postdup_loss"
                    ),
                }
            )
    for b in "123":
        for d in "123":
            losses = sum(state in "12" for state in (b, d))
            pattern_universe.append(
                {
                    "extended_age_bin": "time5_older_than_N4",
                    "pattern": f"3{b}{d}3",
                    "design": (
                        "strict_all_both"
                        if losses == 0
                        else "one_internal_postdup_loss"
                        if losses == 1
                        else "two_internal_postdup_losses"
                    ),
                }
            )

    pattern_counts = Counter(
        row["primary_pattern"] for row in events.values()
    )
    length_counts = Counter(
        row["primary_pattern"]
        for row in events.values()
        if row["strict_minimum_1kb_core"] == "PASS"
    )

    selected_rows = []
    extended_stability: dict[str, str] = {}
    scope_agreement: dict[str, str] = {}
    for event_id, event in events.items():
        main_pattern = patterns[(event_id, "representative", "0.5")]
        age, p_locus, loss_class, loss_nodes = classify(main_pattern)
        if age == "other":
            continue
        threshold_results = [
            classify(patterns[(event_id, "representative", threshold)])[:2]
            for threshold in THRESHOLDS
        ]
        stable = len(set(threshold_results)) == 1
        scope_results = [
            classify(patterns[(event_id, scope, "0.5")])[:2]
            for scope in SCOPES
        ]
        scopes_agree = len(set(scope_results)) == 1
        extended_stability[event_id] = "PASS" if stable else "FAIL"
        scope_agreement[event_id] = "PASS" if scopes_agree else "FAIL"
        d_locus = (
            "locus_B"
            if p_locus == "locus_A"
            else "locus_A"
            if p_locus == "locus_B"
            else "NA"
        )
        selected_rows.append(
            {
                "event_id": event_id,
                "primary_pattern": main_pattern,
                "extended_age_bin": age,
                "loss_tolerance_class": loss_class,
                "postdup_loss_nodes": loss_nodes,
                "provisional_p_locus": p_locus,
                "provisional_d_locus": d_locus,
                "source_call_ids": event["source_call_ids"],
                "representative_call_id": event["representative_call_id"],
                "strict_minimum_1kb_core": event[
                    "strict_minimum_1kb_core"
                ],
                "extended_class_stable_at_0.25_0.5_0.75": (
                    "PASS" if stable else "FAIL"
                ),
                "extended_class_agrees_across_interval_scopes": (
                    "PASS" if scopes_agree else "FAIL"
                ),
                "former_strict_status": event["strict_pd_status"],
                "former_strict_age_bin": event["strict_age_bin"],
            }
        )

    selected_by_pattern = Counter(row["primary_pattern"] for row in selected_rows)
    stable_by_pattern = Counter(
        row["primary_pattern"]
        for row in selected_rows
        if row["strict_minimum_1kb_core"] == "PASS"
        and row["extended_class_stable_at_0.25_0.5_0.75"] == "PASS"
    )
    scope_stable_by_pattern = Counter(
        row["primary_pattern"]
        for row in selected_rows
        if row["strict_minimum_1kb_core"] == "PASS"
        and row["extended_class_stable_at_0.25_0.5_0.75"] == "PASS"
        and row["extended_class_agrees_across_interval_scopes"] == "PASS"
    )
    effective_by_pattern = Counter(
        row["primary_pattern"]
        for row in selected_rows
        if effective_qc
        and effective_qc[row["event_id"]][
            "joint_callable_uppercase_ACGT_aligned_ge1kb"
        ]
        == "PASS"
    )
    effective_stable_by_pattern = Counter(
        row["primary_pattern"]
        for row in selected_rows
        if effective_qc
        and effective_qc[row["event_id"]][
            "joint_callable_uppercase_ACGT_aligned_ge1kb"
        ]
        == "PASS"
        and row["extended_class_stable_at_0.25_0.5_0.75"] == "PASS"
        and row["extended_class_agrees_across_interval_scopes"] == "PASS"
    )
    enumeration_rows = []
    for design in pattern_universe:
        pattern = design["pattern"]
        enumeration_rows.append(
            {
                **design,
                "normalized_events": pattern_counts[pattern],
                "events_ge1kb": length_counts[pattern],
                "events_ge1kb_extended_threshold_stable": (
                    stable_by_pattern[pattern]
                ),
                "events_ge1kb_threshold_and_scope_stable": (
                    scope_stable_by_pattern[pattern]
                ),
                "events_effective_joint_uppercase_aligned_ge1kb": (
                    effective_by_pattern[pattern]
                ),
                "events_effective_ge1kb_threshold_and_scope_stable": (
                    effective_stable_by_pattern[pattern]
                ),
            }
        )

    summary_rows = []
    for age in ("time4_loss_tolerant", "time5_older_than_N4"):
        age_rows = [
            row for row in selected_rows if row["extended_age_bin"] == age
        ]
        for loss_class in sorted(
            {row["loss_tolerance_class"] for row in age_rows}
        ) + ["TOTAL"]:
            subset = (
                age_rows
                if loss_class == "TOTAL"
                else [
                    row
                    for row in age_rows
                    if row["loss_tolerance_class"] == loss_class
                ]
            )
            summary_rows.append(
                {
                    "extended_age_bin": age,
                    "loss_tolerance_class": loss_class,
                    "normalized_events": len(subset),
                    "events_ge1kb": sum(
                        row["strict_minimum_1kb_core"] == "PASS"
                        for row in subset
                    ),
                    "events_ge1kb_extended_threshold_stable": sum(
                        row["strict_minimum_1kb_core"] == "PASS"
                        and row[
                            "extended_class_stable_at_0.25_0.5_0.75"
                        ]
                        == "PASS"
                        for row in subset
                    ),
                    "events_ge1kb_threshold_and_scope_stable": sum(
                        row["strict_minimum_1kb_core"] == "PASS"
                        and row[
                            "extended_class_stable_at_0.25_0.5_0.75"
                        ]
                        == "PASS"
                        and row[
                            "extended_class_agrees_across_interval_scopes"
                        ]
                        == "PASS"
                        for row in subset
                    ),
                    "events_effective_joint_uppercase_aligned_ge1kb": sum(
                        bool(effective_qc)
                        and effective_qc[row["event_id"]][
                            "joint_callable_uppercase_ACGT_aligned_ge1kb"
                        ]
                        == "PASS"
                        for row in subset
                    ),
                    "events_effective_ge1kb_threshold_and_scope_stable": sum(
                        bool(effective_qc)
                        and effective_qc[row["event_id"]][
                            "joint_callable_uppercase_ACGT_aligned_ge1kb"
                        ]
                        == "PASS"
                        and row[
                            "extended_class_stable_at_0.25_0.5_0.75"
                        ]
                        == "PASS"
                        and row[
                            "extended_class_agrees_across_interval_scopes"
                        ]
                        == "PASS"
                        for row in subset
                    ),
                }
            )

    write(
        root / "events/time4_time5_loss_tolerant_events.tsv",
        selected_rows,
    )
    write(
        root / "statistics/time4_time5_pattern_enumeration.tsv",
        enumeration_rows,
    )
    write(
        root / "statistics/time4_time5_loss_tolerant_summary.tsv",
        summary_rows,
    )
    print(
        "Enumerated "
        f"{sum(row['extended_age_bin'] == 'time4_loss_tolerant' for row in selected_rows)} "
        "time4 and "
        f"{sum(row['extended_age_bin'] == 'time5_older_than_N4' for row in selected_rows)} "
        "time5 normalized events"
    )


if __name__ == "__main__":
    main()
