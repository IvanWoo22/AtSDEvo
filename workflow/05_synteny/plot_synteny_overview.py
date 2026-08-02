#!/usr/bin/env python3
"""Plot a compact overview of TAIR12-outgroup synteny results."""

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-tair12-synteny")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = Path(args.root)
    with open(root / "statistics" / "synteny_summary.tsv") as handle:
        data = {row["species"]: row for row in csv.DictReader(handle, delimiter="\t")}

    order = ["Alyrata", "Bstricta", "Dstrictus", "Cviolacea"]
    labels = ["A. lyrata", "B. stricta", "D. strictus", "C. violacea"]
    groups = ["congener", "close", "core", "sister"]
    palette = {
        "congener": "#2166ac", "hybrid": "#762a83", "close": "#4393c3",
        "core": "#92c5de", "sister": "#f4a582", "brassicales": "#d6604d",
        "malvales": "#b2182b",
    }
    colors = [palette[group] for group in groups]
    coverage = [float(data[code]["atha_collinear_coverage_pct"]) for code in order]
    blocks = [int(data[code]["syntenic_blocks"]) for code in order]
    median_size = [float(data[code]["median_gene_pairs_per_block"]) for code in order]
    y = np.arange(len(order))

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.4), sharey=True, gridspec_kw={"width_ratios": [1.5, 1, 1]})
    axes[0].barh(y, coverage, color=colors, edgecolor="white")
    axes[0].set_xlabel("TAIR12 genes in syntenic blocks (%)")
    axes[0].set_yticks(y, labels)
    axes[0].set_xlim(0, max(coverage) * 1.12)
    axes[1].barh(y, blocks, color=colors, edgecolor="white")
    axes[1].set_xlabel("MCScanX blocks")
    axes[2].barh(y, median_size, color=colors, edgecolor="white")
    axes[2].set_xlabel("Median gene pairs per block")
    for ax in axes:
        ax.grid(axis="x", alpha=0.2)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.invert_yaxis()
    axes[0].set_title("Conserved TAIR12 gene coverage")
    axes[1].set_title("Block fragmentation / multiplicity")
    axes[2].set_title("Typical block size")
    fig.suptitle("TAIR12–primary outgroup synteny (BLASTP 1e-10, top 5; MCScanX -s 3 -m 2 -w 0)", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    figures = root / "figures"
    figures.mkdir(exist_ok=True)
    fig.savefig(figures / "pairwise_synteny_overview.png", dpi=220)
    fig.savefig(figures / "pairwise_synteny_overview.pdf")

    with open(root / "statistics" / "atha_chromosome_coverage.tsv") as handle:
        chrom_data = {(row["species"], row["chromosome"]): float(row["coverage_pct"]) for row in csv.DictReader(handle, delimiter="\t")}
    chromosomes = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
    matrix = np.array([[chrom_data[(code, chrom)] for chrom in chromosomes] for code in order])
    heat_fig, heat_ax = plt.subplots(figsize=(7.2, 4.2))
    image = heat_ax.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=90, aspect="auto")
    heat_ax.set_xticks(np.arange(len(chromosomes)), chromosomes)
    heat_ax.set_yticks(np.arange(len(labels)), labels)
    heat_ax.set_xlabel("TAIR12 chromosome")
    heat_ax.set_title("TAIR12 chromosome-level syntenic gene coverage (%)")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            heat_ax.text(j, i, f"{matrix[i, j]:.1f}", ha="center", va="center", fontsize=8,
                         color="white" if matrix[i, j] > 55 else "black")
    heat_fig.colorbar(image, ax=heat_ax, label="Coverage (%)", shrink=0.8)
    heat_fig.tight_layout()
    heat_fig.savefig(figures / "atha_chromosome_synteny_heatmap.png", dpi=220)
    heat_fig.savefig(figures / "atha_chromosome_synteny_heatmap.pdf")
    print("Wrote overview and chromosome heatmap PNG/PDF figures")


if __name__ == "__main__":
    main()
