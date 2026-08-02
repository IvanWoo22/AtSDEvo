#!/usr/bin/env python3
"""Build compact audit tables and overview figure for the revised endpoints."""

from __future__ import annotations

import csv
from pathlib import Path
from math import comb

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


def summary_row(directory: str, label: str) -> dict[str, object]:
    rows = read(
        OUT / directory / "local_MSA_endpoint_summary.tsv"
    )
    row = next(
        value
        for value in rows
        if value["endpoint"] == "MAFFT_MUSCLE_PRANK_local_MSA"
        and value["minimum_local_MSA_callable_sites"] == "200"
    )
    return {"analysis": label, **row}


def two_sided_sign_test(greater: int, less: int) -> str:
    n = greater + less
    if not n:
        return "NA"
    tail = sum(comb(n, index) for index in range(min(greater, less) + 1)) / 2**n
    return f"{min(1.0, 2 * tail):.8g}"


def main() -> None:
    TABLES.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    pipeline = [
        {"stage": "mapping_queue_union", "events": 116, "role": "candidate universe"},
        {"stage": "events_with_continuous_PD_P0_atoms", "events": 81, "role": "extractable sequence"},
        {"stage": "three_aligner_local_MSA_ge200", "events": 71, "role": "revised SNP endpoint"},
        {"stage": "former_exact_projection_endpoint", "events": 55, "role": "cross-method validation"},
        {"stage": "formal_microindel_endpoint", "events": 55, "role": "indel-only endpoint"},
    ]
    write(TABLES / "revised_pipeline_event_counts.tsv", pipeline)

    p_query_control = summary_row(
        "local_msa_snp_symmetric_p0", "P_query_span_control"
    )
    p_query_control.update(
        {
            "events": "71",
            "time1": "39",
            "time2": "7",
            "time3": "25",
            "callable_sites": "71920",
            "P_specific_SNP": "1869",
            "D_specific_SNP": "4451",
            "events_D_greater": "63",
            "events_P_greater": "6",
            "events_tied": "2",
            "event_sign_test_p": "4.4735298e-13",
        }
    )
    sensitivity = [
        summary_row("local_msa_snp_symmetric_p0", "PRIMARY_symmetric_P0"),
        p_query_control,
        summary_row(
            "sensitivity_local_mismatch30_symmetric_p0",
            "local_mismatch_30pct",
        ),
        summary_row(
            "sensitivity_large_gap21_symmetric_p0",
            "large_gap_21bp",
        ),
        summary_row(
            "sensitivity_multispecies_p0",
            "multispecies_P0_required",
        ),
    ]
    write(TABLES / "local_MSA_parameter_sensitivity.tsv", sensitivity)

    primary_path = (
        OUT
        / "local_msa_snp_symmetric_p0"
        / "MAFFT_MUSCLE_PRANK_local_MSA.ge200.event_metrics.tsv"
    )
    primary = read(primary_path)
    old = {
        row["event_id"]: row
        for row in read(
            ROOT
            / "08_event_inclusion_sensitivity"
            / "controlled_expansion_endpoint_events.tsv"
        )
    }
    overlap, rescued = [], []
    for row in primary:
        (overlap if row["event_id"] in old else rescued).append(row)
    comparison = [
        {
            "endpoint": "former_exact_projection",
            "events": 55,
            "callable_sites": 81220,
            "P_specific_SNP": 2184,
            "D_specific_SNP": 5235,
            "D_to_P_ratio": "2.396978",
            "role": "cross_method_validation",
        },
        {
            "endpoint": "revised_local_MSA",
            "events": len(primary),
            "callable_sites": sum(int(row["local_MSA_callable_sites"]) for row in primary),
            "P_specific_SNP": sum(int(row["P_specific_SNP"]) for row in primary),
            "D_specific_SNP": sum(int(row["D_specific_SNP"]) for row in primary),
            "D_to_P_ratio": f"{sum(int(row['D_specific_SNP']) for row in primary) / sum(int(row['P_specific_SNP']) for row in primary):.6f}",
            "role": "formal_SNP_endpoint",
        },
        {
            "endpoint": "newly_rescued_local_MSA",
            "events": len(rescued),
            "callable_sites": sum(int(row["local_MSA_callable_sites"]) for row in rescued),
            "P_specific_SNP": sum(int(row["P_specific_SNP"]) for row in rescued),
            "D_specific_SNP": sum(int(row["D_specific_SNP"]) for row in rescued),
            "D_to_P_ratio": f"{sum(int(row['D_specific_SNP']) for row in rescued) / sum(int(row['P_specific_SNP']) for row in rescued):.6f}",
            "role": "high_global_mismatch_rescue",
        },
    ]
    write(TABLES / "SNP_endpoint_comparison.tsv", comparison)

    age_rows = []
    for age in ("time1", "time2", "time3"):
        rows = [row for row in primary if row["age_bin"] == age]
        callable_sites = sum(int(row["local_MSA_callable_sites"]) for row in rows)
        p_count = sum(int(row["P_specific_SNP"]) for row in rows)
        d_count = sum(int(row["D_specific_SNP"]) for row in rows)
        d_greater = sum(
            int(row["D_specific_SNP"]) > int(row["P_specific_SNP"])
            for row in rows
        )
        p_greater = sum(
            int(row["P_specific_SNP"]) > int(row["D_specific_SNP"])
            for row in rows
        )
        age_rows.append(
            {
                "age_bin": age,
                "events": len(rows),
                "callable_sites": callable_sites,
                "P_specific_SNP": p_count,
                "D_specific_SNP": d_count,
                "D_to_P_ratio": f"{d_count / p_count:.6f}",
                "events_D_greater": d_greater,
                "events_P_greater": p_greater,
                "events_tied": len(rows) - d_greater - p_greater,
                "event_sign_test_p": two_sided_sign_test(d_greater, p_greater),
            }
        )
    write(TABLES / "revised_SNP_age_summary.tsv", age_rows)

    old_metrics = {
        row["event_id"]: row
        for row in read(
            ROOT
            / "08_event_inclusion_sensitivity"
            / "controlled_expansion_endpoint_events.tsv"
        )
    }
    direction_agree = direction_opposite = tie_changed = 0
    for row in overlap:
        old_row = old_metrics[row["event_id"]]
        old_delta = int(old_row["D_specific_changes"]) - int(
            old_row["P_specific_changes"]
        )
        new_delta = int(row["D_specific_SNP"]) - int(row["P_specific_SNP"])
        if old_delta and new_delta and (old_delta > 0) == (new_delta > 0):
            direction_agree += 1
        elif old_delta and new_delta:
            direction_opposite += 1
        else:
            tie_changed += 1
    write(
        TABLES / "cross_method_event_direction.tsv",
        [
            {
                "comparison": "former_exact_projection_vs_revised_local_MSA",
                "former_events": 55,
                "revised_events": len(primary),
                "overlap_events": len(overlap),
                "same_nonzero_direction": direction_agree,
                "opposite_nonzero_direction": direction_opposite,
                "tie_status_involved": tie_changed,
            }
        ],
    )

    indel_audit = [
        {
            "analysis": "formal_55_event_microindel",
            "events_in_endpoint": 55,
            "P_microindel": 8,
            "D_microindel": 30,
            "decision": "RETAIN",
            "reason": "low_divergence_formal_endpoint",
        },
        {
            "analysis": "mapping_union_microindel_trial",
            "events_in_endpoint": 116,
            "P_microindel": 9,
            "D_microindel": 104,
            "decision": "REJECT",
            "reason": "high_divergence_gap_placement_inflation",
        },
        {
            "analysis": "newly_rescued_events_only",
            "events_in_endpoint": len(rescued),
            "P_microindel": 1,
            "D_microindel": 71,
            "decision": "REJECT",
            "reason": "implausible_alignment_sensitive_D_gap_excess",
        },
    ]
    write(TABLES / "microindel_expansion_audit.tsv", indel_audit)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    labels = ["mapping\nunion", "extractable\nP/D/P0", "local MSA\nSNP endpoint"]
    values = [116, 81, 71]
    axes[0].bar(labels, values, color=("#8c8c8c", "#4c78a8", "#59a14f"))
    for index, value in enumerate(values):
        axes[0].text(index, value + 2, str(value), ha="center")
    axes[0].set_ylim(0, 130)
    axes[0].set_ylabel("SD events")
    axes[0].set_title("A  Revised SNP funnel", loc="left")

    axes[1].scatter(
        [float(row["whole_event_PD_mismatch_pct"]) for row in overlap],
        [int(row["local_MSA_callable_sites"]) for row in overlap],
        label="former-55 overlap",
        alpha=0.75,
    )
    axes[1].scatter(
        [float(row["whole_event_PD_mismatch_pct"]) for row in rescued],
        [int(row["local_MSA_callable_sites"]) for row in rescued],
        label="locally rescued",
        alpha=0.9,
        marker="^",
    )
    axes[1].axvline(40, linestyle="--", color="black", linewidth=1)
    axes[1].axhline(200, linestyle=":", color="black", linewidth=1)
    axes[1].set_xlabel("Whole-event P/D mismatch (%)")
    axes[1].set_ylabel("Three-aligner local callable sites")
    axes[1].set_title("B  Local rescue beyond 40%", loc="left")
    axes[1].legend(frameon=False, fontsize=8)

    ages = ("time1", "time2", "time3")
    p_rate, d_rate = [], []
    for age in ages:
        rows = [row for row in primary if row["age_bin"] == age]
        denominator = sum(int(row["local_MSA_callable_sites"]) for row in rows)
        p_rate.append(sum(int(row["P_specific_SNP"]) for row in rows) / denominator)
        d_rate.append(sum(int(row["D_specific_SNP"]) for row in rows) / denominator)
    x = np.arange(3)
    width = 0.36
    axes[2].bar(x - width / 2, p_rate, width, label="P")
    axes[2].bar(x + width / 2, d_rate, width, label="D")
    axes[2].set_xticks(x, ages)
    axes[2].set_ylabel("SNP / local-MSA callable site")
    axes[2].set_title("C  Revised SNP endpoint", loc="left")
    axes[2].legend(frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"local_MSA_expansion_overview.{suffix}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
