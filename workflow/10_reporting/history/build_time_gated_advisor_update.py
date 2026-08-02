#!/usr/bin/env python3
"""Build the current advisor-report summary tables and figures."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUT = Path(__file__).resolve().parents[1]
FIGURES = OUT / "figures"
TABLES = OUT / "tables"

BLUE = "#3973ac"
RED = "#c45145"
GREY = "#6b7280"
GREEN = "#3b8c6e"


def save(fig: plt.Figure, name: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)

    # Current analysis funnel and the two direct primary endpoints.
    stages = [
        ("BISER\ncalls", 4734),
        ("stable\ncalls", 1865),
        ("two-copy\nevents", 1655),
        ("age/P-D\npass", 585),
        ("threshold\nstable", 542),
        ("time1–3\ninput", 362),
        ("P0\nmapped", 110),
        ("local MSA\n≥20 bp", 77),
        ("main SNP\n≥200 bp", 70),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.4), constrained_layout=True)
    x = np.arange(len(stages))
    counts = [value for _, value in stages]
    axes[0].plot(x, counts, marker="o", color=BLUE, linewidth=2)
    axes[0].fill_between(x, counts, 1, color=BLUE, alpha=0.08)
    axes[0].set_yscale("log")
    axes[0].set_xticks(x, [label for label, _ in stages], fontsize=7)
    axes[0].set_ylabel("Number (log scale)")
    axes[0].set_title("A  Event-first analysis funnel", loc="left")
    axes[0].grid(axis="y", alpha=0.2)

    ages = ["time1", "time2", "time3"]
    p_snp = np.array([789, 188, 1072])
    d_snp = np.array([1691, 490, 2153])
    width = 0.36
    x = np.arange(3)
    axes[1].bar(x - width / 2, p_snp, width, label="P", color=BLUE)
    axes[1].bar(x + width / 2, d_snp, width, label="D", color=RED)
    axes[1].set_xticks(x, ages)
    axes[1].set_ylabel("Polarized SNPs")
    axes[1].set_title("B  Direct local-MSA SNP endpoint", loc="left")
    axes[1].legend(frameon=False)
    for i, ratio in enumerate((2.143, 2.606, 2.008)):
        axes[1].text(i, max(p_snp[i], d_snp[i]) * 1.04, f"D/P={ratio:.2f}",
                     ha="center", fontsize=8)
    axes[1].set_ylim(0, 2450)

    p_indel = np.array([10, 0, 3])
    d_indel = np.array([24, 0, 16])
    axes[2].bar(x - width / 2, p_indel, width, label="P", color=BLUE)
    axes[2].bar(x + width / 2, d_indel, width, label="D", color=RED)
    axes[2].set_xticks(x, ages)
    axes[2].set_ylabel("Primary 1–10 bp micro-indels")
    axes[2].set_title("C  Independent local-MSA micro-indels", loc="left")
    axes[2].legend(frameon=False)
    axes[2].text(
        1, 2.0, "no event met\nfull evidence rule",
        ha="center", va="center", fontsize=8, color=GREY,
    )
    save(fig, "current_evidence_overview")

    # Direct vs historical vs topology-aware SNP endpoints.
    labels = [
        "Historical exact\nprojection (55)",
        "Old local MSA\n≥200 (71)",
        "Current direct\n≥200 (70)",
        "Direct + concordant\ndeeper P0 (35)",
        "Fixed-topology ASR\n≥0.90 (70)",
    ]
    p_counts = np.array([2184, 1873, 2049, 603, 2016])
    d_counts = np.array([5235, 4452, 4334, 1756, 5085])
    ratios = d_counts / p_counts
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.5), constrained_layout=True)
    x = np.arange(len(labels))
    axes[0].bar(x - width / 2, p_counts, width, label="P", color=BLUE)
    axes[0].bar(x + width / 2, d_counts, width, label="D", color=RED)
    axes[0].set_xticks(x, labels, fontsize=8)
    axes[0].set_ylabel("Polarized SNPs")
    axes[0].set_title("A  Endpoints use different callable denominators", loc="left")
    axes[0].legend(frameon=False)
    for i, ratio in enumerate(ratios):
        axes[0].text(i, max(p_counts[i], d_counts[i]) + 130,
                     f"{ratio:.2f}", ha="center", fontsize=8)
    axes[0].set_ylim(0, 5900)

    thresholds = np.array([0.80, 0.90, 0.95, 0.99])
    asr_ratio = np.array([2.383, 2.522, 2.777, 3.599])
    callable_bp = np.array([79991, 79469, 78624, 75746])
    axes[1].plot(thresholds, asr_ratio, marker="o", color=GREEN, linewidth=2)
    axes[1].set_xlabel("Minimum JC and K2P DUP posterior")
    axes[1].set_ylabel("D/P SNP count ratio", color=GREEN)
    axes[1].tick_params(axis="y", labelcolor=GREEN)
    axes[1].set_xticks(thresholds)
    axes[1].grid(alpha=0.2)
    ax2 = axes[1].twinx()
    ax2.plot(thresholds, callable_bp, marker="s", color=GREY, linestyle="--")
    ax2.set_ylabel("Callable sites", color=GREY)
    ax2.tick_params(axis="y", labelcolor=GREY)
    axes[1].set_title("B  ASR posterior sensitivity", loc="left")
    save(fig, "endpoint_and_asr_sensitivity")

    # Window-parameter sensitivity for the independent micro-indel analysis.
    labels = ["Historical 55", "Current main", "Window sensitivity"]
    p = np.array([8, 13, 11])
    d = np.array([30, 40, 39])
    fig, ax = plt.subplots(figsize=(7.6, 4.2), constrained_layout=True)
    x = np.arange(3)
    ax.bar(x - width / 2, p, width, label="P", color=BLUE)
    ax.bar(x + width / 2, d, width, label="D", color=RED)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Primary 1–10 bp micro-indels")
    ax.set_title("Micro-indel result expands without reversing direction", loc="left")
    ax.legend(frameon=False)
    for i, ratio in enumerate(d / p):
        ax.text(i, d[i] + 1.3, f"D/P={ratio:.2f}", ha="center", fontsize=9)
    ax.set_ylim(0, 46)
    save(fig, "microindel_endpoint_comparison")


if __name__ == "__main__":
    main()
