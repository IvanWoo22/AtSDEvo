#!/usr/bin/env python3
"""Plot the revised attrition funnel with event count and genomic coverage."""

import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROJECT = Path(
    os.environ.get("SD_PROJECT_ROOT", Path(__file__).resolve().parents[2])
)
OUT = Path(
    os.environ.get(
        "SD_REASSESSMENT_ROOT", PROJECT / "14_pd_polarization_reassessment"
    )
)


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    rows = read(OUT / "statistics/filter_funnel_with_genome_coverage.tsv")
    label_map = {
        "annotation_extended_BISER_calls": "BISER calls",
        "network_stable_source_calls": "network-stable calls",
        "normalized_two_copy_events": "two-copy events",
        "age_free_PD_scope_threshold_stable": "age-free stable P/D",
        "current_time1_time3_input": "current time1–3",
        "at_least_one_mapped_P0": "≥1 mapped P0",
        "continuous_PD_P0_atoms": "continuous atoms",
        "three_aligner_local_MSA_ge20": "local MSA ≥20",
        "three_aligner_local_MSA_ge200": "local MSA ≥200",
        "at_least_one_primary_microindel": "primary micro-indel",
    }
    current = [
        row for row in rows if row["branch"] in {"common", "current_branch"}
    ]
    proposed = [
        row for row in rows if row["branch"] in {"common", "proposed_age_free"}
    ]
    labels = [label_map[row["stage"]] for row in current]
    counts = np.array([int(row["count"]) for row in current])
    coverage = np.array(
        [int(row["genome_union_covered_bp"]) / 1e6 for row in current]
    )
    retained = [float(row["retained_from_previous_pct"]) for row in current]
    x = np.arange(len(current))

    fig, axes = plt.subplots(2, 1, figsize=(12.5, 7.2), constrained_layout=True)
    axes[0].plot(x, counts, marker="o", linewidth=2, color="#3973ac")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Calls/events (log scale)")
    axes[0].set_xticks(x, labels, rotation=32, ha="right")
    axes[0].set_title("A  Current branch and proposed age-free P/D branch", loc="left")
    axes[0].grid(axis="y", alpha=0.2)
    for index, value in enumerate(counts):
        axes[0].annotate(
            f"{value:,}\n({retained[index]:.1f}%)",
            (index, value),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            fontsize=8,
        )
    proposed_x = np.array([0, 1, 2, 3])
    proposed_counts = np.array([int(row["count"]) for row in proposed])
    axes[0].plot(
        proposed_x,
        proposed_counts,
        marker="s",
        linewidth=2,
        linestyle="--",
        color="#c45145",
        label="proposed age-free branch",
    )
    axes[0].annotate(
        f"{proposed_counts[-1]:,}\n({float(proposed[-1]['retained_from_previous_pct']):.1f}%)",
        (proposed_x[-1], proposed_counts[-1]),
        xytext=(0, -34),
        textcoords="offset points",
        ha="center",
        fontsize=8,
        color="#c45145",
    )
    axes[0].legend(frameon=False, loc="upper right")

    axes[1].bar(x, coverage, color="#3b8c6e", label="current branch")
    proposed_coverage = int(proposed[-1]["genome_union_covered_bp"]) / 1e6
    axes[1].bar(
        3,
        proposed_coverage,
        color="#c45145",
        alpha=0.78,
        width=0.55,
        label="proposed age-free P/D",
    )
    axes[1].set_ylabel("Union genomic coverage (Mb)")
    axes[1].set_xticks(x, labels, rotation=32, ha="right")
    axes[1].set_title(
        "B  Nonredundant union of both-copy genomic intervals", loc="left"
    )
    for index, value in enumerate(coverage):
        axes[1].text(index, value + 0.2, f"{value:.2f}", ha="center", fontsize=8)
    axes[1].set_ylim(0, max(coverage) * 1.15)
    axes[1].text(
        3,
        proposed_coverage + 0.2,
        f"{proposed_coverage:.2f}",
        ha="center",
        fontsize=8,
        color="#8f3028",
    )
    axes[1].legend(frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"figures/revised_filter_funnel_with_coverage.{suffix}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
