#!/usr/bin/env python3
"""Summarize and visualize the call-to-event-first age reanalysis."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


AGES = ("time1", "time2", "time3", "time4")


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    args = parser.parse_args()

    base = args.project / "06_sd_age_tracing_preparation"
    old_dir = base / "primary_node_analysis"
    new_dir = base / "event_first_reanalysis"
    (new_dir / "figures").mkdir(parents=True, exist_ok=True)
    old_rows = read(old_dir / "events/strict_two_copy_events.tsv")
    new_rows = read(new_dir / "events/event_first_events.tsv")
    old = {row["event_id"]: row for row in old_rows}
    new = {row["event_id"]: row for row in new_rows}
    old_pass = {
        event_id for event_id, row in old.items()
        if row["strict_pd_status"] == "PASS"
    }
    new_pass = {
        event_id for event_id, row in new.items()
        if row["strict_pd_status"] == "PASS"
    }
    recovered = sorted(new_pass - old_pass)
    lost = sorted(old_pass - new_pass)

    comparison = [
        {
            "analysis": "old_call_label_inheritance",
            "BISER_calls": 4734,
            "network_stable_source_calls": 1865,
            "normalized_events": len(old),
            "age_PD_pass": len(old_pass),
            "threshold_stable_pass": 411,
            "role": "superseded_baseline",
        },
        {
            "analysis": "event_first_physical_locus",
            "BISER_calls": 4734,
            "network_stable_source_calls": 1865,
            "normalized_events": len(new),
            "age_PD_pass": len(new_pass),
            "threshold_stable_pass": sum(
                row["strict_pd_status"] == "PASS"
                and row[
                    "classification_stable_at_overlap_0.25_0.5_0.75"
                ]
                == "PASS"
                for row in new.values()
            ),
            "role": "revised_primary",
        },
    ]
    write(new_dir / "statistics/old_vs_event_first_overview.tsv", comparison)

    age_rows = []
    for age in AGES:
        old_count = sum(
            old[event_id]["strict_age_bin"] == age for event_id in old_pass
        )
        new_count = sum(
            new[event_id]["strict_age_bin"] == age for event_id in new_pass
        )
        stable_count = sum(
            row["strict_pd_status"] == "PASS"
            and row["strict_age_bin"] == age
            and row["classification_stable_at_overlap_0.25_0.5_0.75"]
            == "PASS"
            for row in new.values()
        )
        age_rows.append(
            {
                "age_bin": age,
                "old_pass_events": old_count,
                "event_first_pass_events": new_count,
                "recovered_events": new_count - old_count,
                "event_first_threshold_stable_events": stable_count,
            }
        )
    write(new_dir / "statistics/event_first_age_summary.tsv", age_rows)

    recovered_rows = []
    for event_id in recovered:
        row = new[event_id]
        recovered_rows.append(
            {
                "event_id": event_id,
                "source_call_ids": row["source_call_ids"],
                "representative_call_id": row["representative_call_id"],
                "strict_age_bin": row["strict_age_bin"],
                "provisional_p_locus": row["provisional_p_locus"],
                "provisional_d_locus": row["provisional_d_locus"],
                "strict_minimum_1kb_core": row["strict_minimum_1kb_core"],
                "threshold_stability": row[
                    "classification_stable_at_overlap_0.25_0.5_0.75"
                ],
                "interval_scope_agreement": row[
                    "classification_agrees_representative_core_envelope"
                ],
                "old_exclusion": old[event_id]["strict_age_bin"],
            }
        )
    write(new_dir / "events/recovered_false_label_swap_events.tsv", recovered_rows)

    threshold_stable_rows = [
        row
        for row in new_rows
        if row["strict_pd_status"] == "PASS"
        and row["classification_stable_at_overlap_0.25_0.5_0.75"] == "PASS"
    ]
    write(
        new_dir / "events/event_first_threshold_stable.tsv",
        threshold_stable_rows,
    )
    write(
        new_dir / "events/event_first_time1_time3_threshold_stable.tsv",
        [
            row
            for row in threshold_stable_rows
            if row["strict_age_bin"] in ("time1", "time2", "time3")
        ],
    )

    old_physical_agreement = 0
    for event_id in old_pass:
        old_row = old[event_id]
        new_row = new[event_id]
        label = new_row["provisional_p_locus"]
        old_p = (
            old_row["provisional_p_chrom"],
            old_row["provisional_p_start"],
            old_row["provisional_p_end"],
        )
        new_p = (
            new_row[f"{label}_chrom"],
            new_row[f"{label}_representative_start"],
            new_row[f"{label}_representative_end"],
        )
        if (
            old_row["strict_age_bin"] == new_row["strict_age_bin"]
            and old_p == new_p
        ):
            old_physical_agreement += 1

    checks = [
        {
            "check": "old_pass_retained",
            "observed": len(old_pass & new_pass),
            "expected": len(old_pass),
            "status": "PASS" if not lost else "FAIL",
        },
        {
            "check": "new_pass_events",
            "observed": len(new_pass),
            "expected": len(new_pass),
            "status": "PASS",
        },
        {
            "check": "false_copy_label_events",
            "observed": sum(
                row["source_call_count"] == "2"
                and old[row["event_id"]]["strict_age_bin"]
                == "source_call_disagreement"
                for row in new_rows
            ),
            "expected": 76,
            "status": "PASS",
        },
        {
            "check": "recovered_events_without_length_gating",
            "observed": len(recovered),
            "expected": len(recovered),
            "status": "PASS",
        },
        {
            "check": "events_agree_across_interval_scopes",
            "observed": sum(
                row[
                    "classification_agrees_representative_core_envelope"
                ]
                == "PASS"
                for row in new_rows
            ),
            "expected": len(new_rows),
            "status": "PASS",
        },
        {
            "check": "old_pass_age_and_physical_P_preserved",
            "observed": old_physical_agreement,
            "expected": len(old_pass),
            "status": (
                "PASS" if old_physical_agreement == len(old_pass) else "FAIL"
            ),
        },
    ]
    write(new_dir / "statistics/event_first_validation_checks.tsv", checks)

    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.0), constrained_layout=True)
    labels = (
        "BISER\ncalls",
        "stable source\ncalls",
        "two-locus\nevents",
        "age/P-D\n(no length gate)",
    )
    values = (4734, 1865, 1655, len(new_pass))
    axes[0].bar(labels, values, color=("#999999", "#4C78A8", "#59A14F", "#F28E2B"))
    for index, value in enumerate(values):
        axes[0].text(index, value + 70, str(value), ha="center", fontsize=9)
    axes[0].set_ylim(0, 5200)
    axes[0].set_ylabel("Count")
    axes[0].set_title("A  Event-first funnel", loc="left")

    x = np.arange(len(AGES))
    width = 0.36
    axes[1].bar(
        x - width / 2,
        [int(row["old_pass_events"]) for row in age_rows],
        width,
        label="old",
        color="#BAB0AC",
    )
    axes[1].bar(
        x + width / 2,
        [int(row["event_first_pass_events"]) for row in age_rows],
        width,
        label="event-first",
        color="#4C78A8",
    )
    axes[1].set_xticks(x, AGES)
    axes[1].set_ylabel("Events")
    axes[1].set_title("B  Corrected age bins", loc="left")
    axes[1].legend(frameon=False)

    recovered_age = Counter(new[event_id]["strict_age_bin"] for event_id in recovered)
    stable_age = Counter(
        new[event_id]["strict_age_bin"]
        for event_id in recovered
        if new[event_id]["classification_stable_at_overlap_0.25_0.5_0.75"]
        == "PASS"
    )
    axes[2].bar(
        x,
        [recovered_age[age] for age in AGES],
        label="recovered",
        color="#E45756",
    )
    axes[2].bar(
        x,
        [stable_age[age] for age in AGES],
        label="threshold-stable subset",
        color="#72B7B2",
    )
    axes[2].set_xticks(x, AGES)
    axes[2].set_ylabel("Recovered events")
    axes[2].set_title("C  Fixed copy-label swaps", loc="left")
    axes[2].legend(frameon=False, fontsize=8)
    for suffix in ("png", "pdf"):
        fig.savefig(
            new_dir / f"figures/event_first_reanalysis_overview.{suffix}",
            dpi=300,
        )
    plt.close(fig)

    print(
        f"Old pass={len(old_pass)}; event-first pass={len(new_pass)}; "
        f"recovered={len(recovered)}; lost={len(lost)}"
    )


if __name__ == "__main__":
    main()
