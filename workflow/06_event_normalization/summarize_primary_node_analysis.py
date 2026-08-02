#!/usr/bin/env python3
"""Summarize and plot the frozen four-primary-node SD analysis."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
AGE_BINS = ("time1", "time2", "time3", "time4")


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
    parser.add_argument("--project", required=True, type=Path)
    args = parser.parse_args()

    synteny_dir = args.project / "05_mcscanx_synteny"
    analysis_dir = (
        args.project / "06_sd_age_tracing_preparation/primary_node_analysis"
    )
    event_dir = analysis_dir / "events"
    statistics_dir = analysis_dir / "statistics"
    figure_dir = analysis_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    synteny = read_tsv(synteny_dir / "statistics/synteny_summary.tsv")
    calls = read_tsv(analysis_dir / "call_evidence.threshold_0.5.tsv")
    events = read_tsv(event_dir / "strict_two_copy_events.tsv")
    threshold_rows = read_tsv(
        statistics_dir / "call_classification_by_threshold.tsv"
    )
    threshold_classification = {
        (row["call_id"], row["threshold"]): (
            row["strict_status"],
            row["strict_age_bin"],
            row["provisional_p_copy"],
        )
        for row in threshold_rows
    }
    normalization = {
        row["metric"]: int(row["count"])
        for row in read_tsv(event_dir / "normalization_summary.tsv")
    }

    event_status = Counter(
        (row["strict_pd_status"], row["strict_age_bin"]) for row in events
    )
    pass_events = [row for row in events if row["strict_pd_status"] == "PASS"]
    sensitivity_rows = []
    for event in events:
        call_ids = event["source_call_ids"].split(",")
        main = (
            event["strict_pd_status"],
            event["strict_age_bin"],
            event["provisional_p_copy"],
        )
        stable = all(
            threshold_classification[(call_id, threshold)] == main
            for call_id in call_ids
            for threshold in ("0.25", "0.5", "0.75")
        )
        sensitivity_rows.append(
            {
                "event_id": event["event_id"],
                "main_strict_pd_status": event["strict_pd_status"],
                "main_strict_age_bin": event["strict_age_bin"],
                "main_provisional_p_copy": event["provisional_p_copy"],
                "classification_stable_at_0.25_0.5_0.75": (
                    "PASS" if stable else "FAIL"
                ),
            }
        )
    write_tsv(
        statistics_dir / "event_node_overlap_threshold_sensitivity.tsv",
        sensitivity_rows,
    )
    stable_pass_events = sum(
        row["main_strict_pd_status"] == "PASS"
        and row["classification_stable_at_0.25_0.5_0.75"] == "PASS"
        for row in sensitivity_rows
    )
    age_rows = []
    for age in AGE_BINS:
        selected = [row for row in pass_events if row["strict_age_bin"] == age]
        p_copy = Counter(row["provisional_p_copy"] for row in selected)
        orientation = Counter(row["relative_orientation"] for row in selected)
        relation = Counter(
            "intrachromosomal"
            if row["representative_copy1_chrom"]
            == row["representative_copy2_chrom"]
            else "interchromosomal"
            for row in selected
        )
        lengths = [int(row["representative_alignment_span_bp"]) for row in selected]
        age_rows.append(
            {
                "age_bin": age,
                "strict_events": len(selected),
                "p_is_copy1": p_copy["copy1"],
                "p_is_copy2": p_copy["copy2"],
                "same_orientation": orientation["same"],
                "opposite_orientation": orientation["opposite"],
                "intrachromosomal": relation["intrachromosomal"],
                "interchromosomal": relation["interchromosomal"],
                "median_alignment_span_bp": f"{median(lengths):.1f}",
                "min_alignment_span_bp": min(lengths),
                "max_alignment_span_bp": max(lengths),
            }
        )
    write_tsv(statistics_dir / "strict_event_age_summary.tsv", age_rows)

    exclusion_order = (
        "older_than_N4_unpolarized",
        "boundary_node_uninformative",
        "nonmonotonic_both_reappears",
        "all_zero",
        "source_call_disagreement",
        "below_1kb_strict_sd_core",
        "older_nodes_copy_conflict",
    )
    disposition_rows = [
        {
            "status": "PASS",
            "category": age,
            "events": event_status[("PASS", age)],
            "fraction_of_1655_events": (
                f"{100 * event_status[('PASS', age)] / len(events):.4f}"
            ),
        }
        for age in AGE_BINS
    ]
    disposition_rows.extend(
        {
            "status": "EXCLUDE",
            "category": reason,
            "events": event_status[("EXCLUDE", reason)],
            "fraction_of_1655_events": (
                f"{100 * event_status[('EXCLUDE', reason)] / len(events):.4f}"
            ),
        }
        for reason in exclusion_order
    )
    write_tsv(statistics_dir / "strict_event_disposition.tsv", disposition_rows)

    call_status = Counter((row["strict_status"], row["strict_age_bin"]) for row in calls)
    overview = [
        {"metric": "annotation_extended_biser_calls", "value": len(calls)},
        {
            "metric": "call_level_strict_primary_node_pass",
            "value": sum(status == "PASS" for status, _ in call_status.elements()),
        },
        {
            "metric": "network_stable_source_calls",
            "value": normalization["strict_calls_stable_at_both"],
        },
        {
            "metric": "network_stable_two_copy_events",
            "value": len(events),
        },
        {"metric": "strict_pd_events_after_1kb_filter", "value": len(pass_events)},
        {
            "metric": "strict_pd_events_stable_at_node_overlap_0.25_0.5_0.75",
            "value": stable_pass_events,
        },
        {
            "metric": "strict_pd_event_fraction_of_biser_calls_pct",
            "value": f"{100 * len(pass_events) / len(calls):.4f}",
        },
        {
            "metric": "strict_pd_event_fraction_of_normalized_events_pct",
            "value": f"{100 * len(pass_events) / len(events):.4f}",
        },
    ]
    write_tsv(statistics_dir / "primary_analysis_overview.tsv", overview)

    plt.rcParams.update({"font.size": 9, "figure.dpi": 150})
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))

    coverage = {row["species"]: float(row["atha_collinear_coverage_pct"]) for row in synteny}
    axes[0].bar(SPECIES, [coverage[x] for x in SPECIES], color="#4472C4")
    axes[0].set_ylabel("TAIR12 genes in synteny (%)")
    axes[0].set_ylim(0, 90)
    axes[0].set_title("Primary-node collinearity")
    axes[0].tick_params(axis="x", rotation=35)

    attrition_labels = ("BISER\ncalls", "Stable\ncalls", "2-copy\nevents", "Strict P/D\n≥1 kb")
    attrition_values = (
        len(calls),
        normalization["strict_calls_stable_at_both"],
        len(events),
        len(pass_events),
    )
    axes[1].bar(attrition_labels, attrition_values, color=("#888888", "#70AD47", "#5B9BD5", "#ED7D31"))
    axes[1].set_ylabel("Count")
    axes[1].set_title("Conservative event attrition")
    for index, value in enumerate(attrition_values):
        axes[1].text(index, value + 70, str(value), ha="center", fontsize=8)
    axes[1].set_ylim(0, 5200)

    age_counts = [event_status[("PASS", age)] for age in AGE_BINS]
    axes[2].bar(AGE_BINS, age_counts, color=("#4E79A7", "#59A14F", "#F28E2B", "#E15759"))
    axes[2].set_ylabel("Strict events")
    axes[2].set_title("Primary-node age/P-D bins")
    for index, value in enumerate(age_counts):
        axes[2].text(index, value + 3, str(value), ha="center", fontsize=8)
    axes[2].set_ylim(0, max(age_counts) * 1.18)

    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(figure_dir / f"primary_node_strict_sd_overview.{suffix}", bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote summaries and figures for {len(pass_events)} strict P/D events")


if __name__ == "__main__":
    main()
