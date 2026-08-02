#!/usr/bin/env python3
"""Build compact, reproducible summary tables/figure for thesis part I."""

from __future__ import annotations

import csv
import math
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parents[1]
TABLES = OUT / "tables"
FIGURES = OUT / "figures"


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


def exact_binomial(k: int, n: int) -> float:
    observed = math.comb(n, k)
    return sum(
        math.comb(n, value)
        for value in range(n + 1)
        if math.comb(n, value) <= observed
    ) / (2**n)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    overview = {
        row["metric"]: int(float(row["value"]))
        for row in read(
            ROOT
            / "06_sd_age_tracing_preparation/primary_node_analysis/statistics"
            / "primary_analysis_overview.tsv"
        )
    }
    pipeline = [
        {
            "stage_order": 1,
            "stage": "annotation_extended_BISER_calls",
            "units": "calls",
            "count": overview["annotation_extended_biser_calls"],
            "role": "SD discovery only",
        },
        {
            "stage_order": 2,
            "stage": "network_stable_source_calls",
            "units": "calls",
            "count": overview["network_stable_source_calls"],
            "role": "stable at reciprocal overlap 0.5 and 0.8",
        },
        {
            "stage_order": 3,
            "stage": "normalized_two_copy_events",
            "units": "events",
            "count": overview["network_stable_two_copy_events"],
            "role": "nonredundant two-locus events",
        },
        {
            "stage_order": 4,
            "stage": "strict_age_PD_events_ge1kb",
            "units": "events",
            "count": overview["strict_pd_events_after_1kb_filter"],
            "role": "time1-time4 plus strict P/D",
        },
        {
            "stage_order": 5,
            "stage": "node_threshold_stable_strict_events",
            "units": "events",
            "count": overview[
                "strict_pd_events_stable_at_node_overlap_0.25_0.5_0.75"
            ],
            "role": "age and P stable across node-overlap thresholds",
        },
        {
            "stage_order": 6,
            "stage": "time1_time3_threshold_stable_candidates",
            "units": "events",
            "count": 264,
            "role": "time4 excluded by study scope",
        },
        {
            "stage_order": 7,
            "stage": "frozen_primary_endpoint",
            "units": "events",
            "count": 44,
            "role": "1 kb strict sequence endpoint",
        },
        {
            "stage_order": 8,
            "stage": "controlled_expanded_endpoint",
            "units": "events",
            "count": 55,
            "role": "500 bp plus prespecified controlled sensitivity tiers",
        },
    ]
    write(TABLES / "analysis_pipeline_attrition.tsv", pipeline)

    controlled = read(
        ROOT / "08_event_inclusion_sensitivity/controlled_expansion_endpoint_events.tsv"
    )
    snp_age = []
    for age in ("time1", "time2", "time3"):
        rows = [row for row in controlled if row["age_bin"] == age]
        p_count = sum(int(row["P_specific_changes"]) for row in rows)
        d_count = sum(int(row["D_specific_changes"]) for row in rows)
        d_greater = sum(
            int(row["D_specific_changes"]) > int(row["P_specific_changes"])
            for row in rows
        )
        p_greater = sum(
            int(row["P_specific_changes"]) > int(row["D_specific_changes"])
            for row in rows
        )
        ties = len(rows) - d_greater - p_greater
        callable_bp = sum(int(row["primary_bidir_callable_sites"]) for row in rows)
        snp_age.append(
            {
                "age_bin": age,
                "events": len(rows),
                "callable_bp": callable_bp,
                "P_specific_SNP": p_count,
                "D_specific_SNP": d_count,
                "D_to_P_ratio": f"{d_count / p_count:.6f}",
                "P_rate": f"{p_count / callable_bp:.8f}",
                "D_rate": f"{d_count / callable_bp:.8f}",
                "events_D_greater": d_greater,
                "events_P_greater": p_greater,
                "events_tied": ties,
                "event_sign_test_p": (
                    f"{exact_binomial(d_greater, d_greater + p_greater):.8g}"
                ),
            }
        )
    write(TABLES / "SNP_age_stratified_summary.tsv", snp_age)

    event_info = {
        row["event_id"]: row
        for row in read(
            ROOT
            / "08_event_inclusion_sensitivity/core_500/inputs"
            / "high_priority_core_eligible.tsv"
        )
    }
    copy_order_rows = []
    for p_copy in ("copy1", "copy2"):
        rows = [
            row
            for row in controlled
            if event_info[row["event_id"]]["provisional_p_copy"] == p_copy
        ]
        p_count = sum(int(row["P_specific_changes"]) for row in rows)
        d_count = sum(int(row["D_specific_changes"]) for row in rows)
        d_greater = sum(
            int(row["D_specific_changes"]) > int(row["P_specific_changes"])
            for row in rows
        )
        p_greater = sum(
            int(row["P_specific_changes"]) > int(row["D_specific_changes"])
            for row in rows
        )
        copy_order_rows.append(
            {
                "P_copy_index": p_copy,
                "events": len(rows),
                "P_specific_SNP": p_count,
                "D_specific_SNP": d_count,
                "D_to_P_ratio": f"{d_count / p_count:.6f}",
                "events_D_greater": d_greater,
                "events_P_greater": p_greater,
                "event_sign_test_p": f"{exact_binomial(d_greater, d_greater + p_greater):.8g}",
            }
        )
    write(TABLES / "SNP_copy_order_bias_check.tsv", copy_order_rows)

    indel_root = (
        ROOT
        / "09_variant_type_analysis_55/microindel_denovo_msa"
    )
    indels = [
        row
        for row in read(indel_root / "denovo_microindel_inference.tsv")
        if row["evidence_tier"] == "PRIMARY"
    ]
    direction_rows = []
    for age in ("time1", "time2", "time3", "all"):
        rows = indels if age == "all" else [
            row for row in indels if row["age_bin"] == age
        ]
        counts = Counter(row["parsimonious_direction"] for row in rows)
        direction_rows.append(
            {
                "age_bin": age,
                "P_insertion": counts["P_insertion"],
                "P_deletion": counts["P_deletion"],
                "D_insertion": counts["D_insertion"],
                "D_deletion": counts["D_deletion"],
                "primary_microindels": len(rows),
            }
        )
    write(TABLES / "microindel_direction_type_summary.tsv", direction_rows)

    indel_order_rows = []
    for p_copy in ("copy1", "copy2"):
        rows = [
            row
            for row in indels
            if event_info[row["event_id"]]["provisional_p_copy"] == p_copy
        ]
        indel_order_rows.append(
            {
                "P_copy_index": p_copy,
                "primary_microindels": len(rows),
                "P_branch_microindels": sum(
                    row["parsimonious_direction"].startswith("P_") for row in rows
                ),
                "D_branch_microindels": sum(
                    row["parsimonious_direction"].startswith("D_") for row in rows
                ),
                "independent_events": len({row["event_id"] for row in rows}),
            }
        )
    write(TABLES / "microindel_copy_order_bias_check.tsv", indel_order_rows)

    synteny = read(ROOT / "05_mcscanx_synteny/statistics/synteny_summary.tsv")
    write(
        TABLES / "outgroup_synteny_ladder.tsv",
        [
            {
                "species": row["species"],
                "syntenic_blocks": row["syntenic_blocks"],
                "collinear_gene_pairs": row["collinear_gene_pairs"],
                "unique_TAIR12_collinear_genes": row[
                    "unique_atha_collinear_genes"
                ],
                "TAIR12_collinear_coverage_pct": row[
                    "atha_collinear_coverage_pct"
                ],
            }
            for row in synteny
        ],
    )

    evidence = [
        {
            "question": "Can stable nonredundant SD events be defined?",
            "result": "1655 normalized stable two-copy events; 446 strict age/P-D passes",
            "status": "SUPPORTED",
            "main_limit": "complex/nonmonotonic events remain excluded",
        },
        {
            "question": "Can an ordered outgroup synteny ladder be recovered?",
            "result": "TAIR12 gene coverage declines 81.23% to 43.31% across four nodes",
            "status": "SUPPORTED",
            "main_limit": "block count itself is not a molecular clock",
        },
        {
            "question": "Can strict P/D-polarized sequence endpoints be formed?",
            "result": "44 frozen primary and 55 controlled events; 81,220 callable bp",
            "status": "SUPPORTED",
            "main_limit": "time2 remains the smallest stratum",
        },
        {
            "question": "Do P and D differ in derived SNP accumulation?",
            "result": "P=2184, D=5235, D/P=2.397; all six classes D>P",
            "status": "SUPPORTED_ROBUST",
            "main_limit": "does not yet separate mutation, selection, and conversion",
        },
        {
            "question": "Do P and D differ in de novo micro-indel accumulation?",
            "result": "main P=8, D=30; 28/38 exact calls robust to window parameters",
            "status": "SUPPORTED_PRELIMINARY",
            "main_limit": "time2 has no primary micro-indel; coordinate sensitivity remains",
        },
        {
            "question": "Is the observed asymmetry a copy1/copy2 ordering artifact?",
            "result": "D excess occurs in both P=copy1 and P=copy2 strata",
            "status": "NOT_EXPLAINED_BY_COPY_ORDER",
            "main_limit": "other ascertainment effects still require monitoring",
        },
    ]
    write(TABLES / "first_part_evidence_matrix.tsv", evidence)

    # Compact presentation figure: event funnel, SNP by age, micro-indel by age.
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    funnel = [4734, 1865, 1655, 446, 411, 264, 55]
    labels = ["BISER", "stable calls", "events", "age/P-D", "stable P-D", "t1–3", "endpoint"]
    axes[0].plot(range(len(funnel)), funnel, marker="o", color="#4c78a8")
    axes[0].set_yscale("log")
    axes[0].set_xticks(range(len(funnel)), labels, rotation=45, ha="right")
    axes[0].set_ylabel("Count (log scale)")
    axes[0].set_title("A  SD-to-endpoint funnel", loc="left")

    x = np.arange(3)
    width = 0.36
    axes[1].bar(
        x - width / 2,
        [float(row["P_rate"]) for row in snp_age],
        width,
        label="P",
        color="#3973ac",
    )
    axes[1].bar(
        x + width / 2,
        [float(row["D_rate"]) for row in snp_age],
        width,
        label="D",
        color="#c45145",
    )
    axes[1].set_xticks(x, ("time1", "time2", "time3"))
    axes[1].set_ylabel("Polarized SNP / callable bp")
    axes[1].set_title("B  SNP asymmetry", loc="left")
    axes[1].legend(frameon=False)

    indel_age = read(indel_root / "age_stratified_PD_microindel_summary.tsv")
    axes[2].bar(
        x - width / 2,
        [int(row["P_branch_microindels"]) for row in indel_age],
        width,
        label="P",
        color="#3973ac",
    )
    axes[2].bar(
        x + width / 2,
        [int(row["D_branch_microindels"]) for row in indel_age],
        width,
        label="D",
        color="#c45145",
    )
    axes[2].set_xticks(x, ("time1", "time2", "time3"))
    axes[2].set_ylabel("Primary de novo micro-indels")
    axes[2].set_title("C  Micro-indel asymmetry", loc="left")
    axes[2].legend(frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"first_part_evidence_overview.{suffix}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
