#!/usr/bin/env python3
"""Statistical analysis and figures for polarized P/D substitution rates."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


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


def bootstrap(
    rows: list[dict[str, str]], prefix: str, rng: np.random.Generator
) -> tuple[float, float, float, float]:
    ratios = []
    median_differences = []
    for _ in range(20_000):
        indices = rng.integers(0, len(rows), len(rows))
        sampled = [rows[index] for index in indices]
        p_count = sum(int(row[f"{prefix}_P_specific"]) for row in sampled)
        d_count = sum(int(row[f"{prefix}_D_specific"]) for row in sampled)
        if p_count:
            ratios.append(d_count / p_count)
        median_differences.append(
            np.median(
                [
                    float(row[f"{prefix}_D_minus_P_rate"])
                    for row in sampled
                ]
            )
        )
    ratio_low, ratio_high = np.quantile(ratios, [0.025, 0.975])
    diff_low, diff_high = np.quantile(
        median_differences, [0.025, 0.975]
    )
    return ratio_low, ratio_high, diff_low, diff_high


def summarize(
    name: str,
    rows: list[dict[str, str]],
    prefix: str,
    rng: np.random.Generator,
) -> dict[str, object]:
    p_counts = np.array(
        [int(row[f"{prefix}_P_specific"]) for row in rows], dtype=int
    )
    d_counts = np.array(
        [int(row[f"{prefix}_D_specific"]) for row in rows], dtype=int
    )
    callable_sites = np.array(
        [int(row[f"{prefix}_callable_sites"]) for row in rows], dtype=int
    )
    p_rates = p_counts / callable_sites
    d_rates = d_counts / callable_sites
    differences = d_rates - p_rates
    p_total = int(p_counts.sum())
    d_total = int(d_counts.sum())
    non_ties = differences[differences != 0]
    binomial = stats.binomtest(
        d_total, p_total + d_total, p=0.5, alternative="two-sided"
    )
    sign = stats.binomtest(
        int((non_ties > 0).sum()),
        len(non_ties),
        p=0.5,
        alternative="two-sided",
    )
    wilcoxon = stats.wilcoxon(
        d_rates,
        p_rates,
        zero_method="wilcox",
        alternative="two-sided",
        method="auto",
    )
    ratio_low, ratio_high, diff_low, diff_high = bootstrap(
        rows, prefix, rng
    )
    return {
        "analysis": name,
        "metric_prefix": prefix,
        "events": len(rows),
        "callable_sites": int(callable_sites.sum()),
        "P_specific_changes": p_total,
        "D_specific_changes": d_total,
        "D_to_P_count_ratio": f"{d_total / p_total:.6f}",
        "D_to_P_ratio_bootstrap_CI95_low": f"{ratio_low:.6f}",
        "D_to_P_ratio_bootstrap_CI95_high": f"{ratio_high:.6f}",
        "median_P_specific_rate": f"{np.median(p_rates):.8f}",
        "median_D_specific_rate": f"{np.median(d_rates):.8f}",
        "median_D_minus_P_rate": f"{np.median(differences):.8f}",
        "median_difference_bootstrap_CI95_low": f"{diff_low:.8f}",
        "median_difference_bootstrap_CI95_high": f"{diff_high:.8f}",
        "events_D_gt_P": int((differences > 0).sum()),
        "events_P_gt_D": int((differences < 0).sum()),
        "events_tied": int((differences == 0).sum()),
        "aggregate_binomial_p": f"{binomial.pvalue:.12g}",
        "event_sign_test_p": f"{sign.pvalue:.12g}",
        "paired_wilcoxon_statistic": f"{wilcoxon.statistic:.6f}",
        "paired_wilcoxon_p": f"{wilcoxon.pvalue:.12g}",
    }


def summarize_terminal(
    name: str,
    rows: list[dict[str, str]],
    rng: np.random.Generator,
) -> dict[str, object]:
    p_counts = np.array(
        [int(row["P_terminal_mismatches"]) for row in rows], dtype=int
    )
    d_counts = np.array(
        [int(row["D_terminal_mismatches"]) for row in rows], dtype=int
    )
    callable_sites = np.array(
        [int(row["joint_callable_sites"]) for row in rows], dtype=int
    )
    p_rates = p_counts / callable_sites
    d_rates = d_counts / callable_sites
    differences = d_rates - p_rates
    non_ties = differences[differences != 0]
    bootstrap_ratios = []
    bootstrap_differences = []
    for _ in range(20_000):
        indices = rng.integers(0, len(rows), len(rows))
        p_sum = int(p_counts[indices].sum())
        d_sum = int(d_counts[indices].sum())
        if p_sum:
            bootstrap_ratios.append(d_sum / p_sum)
        bootstrap_differences.append(np.median(differences[indices]))
    ratio_ci = np.quantile(bootstrap_ratios, [0.025, 0.975])
    diff_ci = np.quantile(bootstrap_differences, [0.025, 0.975])
    binomial = stats.binomtest(
        int(d_counts.sum()),
        int(p_counts.sum() + d_counts.sum()),
        0.5,
        alternative="two-sided",
    )
    sign = stats.binomtest(
        int((non_ties > 0).sum()),
        len(non_ties),
        0.5,
        alternative="two-sided",
    )
    wilcoxon = stats.wilcoxon(
        d_rates, p_rates, zero_method="wilcox", method="auto"
    )
    return {
        "analysis": name,
        "events": len(rows),
        "callable_sites": int(callable_sites.sum()),
        "P_terminal_mismatches": int(p_counts.sum()),
        "D_terminal_mismatches": int(d_counts.sum()),
        "D_to_P_count_ratio": f"{d_counts.sum() / p_counts.sum():.6f}",
        "D_to_P_ratio_bootstrap_CI95_low": f"{ratio_ci[0]:.6f}",
        "D_to_P_ratio_bootstrap_CI95_high": f"{ratio_ci[1]:.6f}",
        "median_D_minus_P_terminal_rate": f"{np.median(differences):.8f}",
        "median_difference_bootstrap_CI95_low": f"{diff_ci[0]:.8f}",
        "median_difference_bootstrap_CI95_high": f"{diff_ci[1]:.8f}",
        "events_D_gt_P": int((differences > 0).sum()),
        "events_P_gt_D": int((differences < 0).sum()),
        "events_tied": int((differences == 0).sum()),
        "aggregate_binomial_p": f"{binomial.pvalue:.12g}",
        "event_sign_test_p": f"{sign.pvalue:.12g}",
        "paired_wilcoxon_p": f"{wilcoxon.pvalue:.12g}",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", required=True, type=Path)
    args = parser.parse_args()
    source = read_tsv(
        args.pilot / "sequence_variation/event_pd_substitution_metrics.tsv"
    )
    rng = np.random.default_rng(20260724)

    no_high_divergence = lambda row: (
        row["endpoint_sensitivity_flags"] == "NONE"
    )
    definitions = [
        (
            "PRIMARY_bidir_ge500_no_high_divergence",
            "primary_bidir",
            lambda row: int(row["primary_bidir_callable_sites"]) >= 500
            and no_high_divergence(row),
        ),
        (
            "all_strict_bidir_ge500",
            "primary_bidir",
            lambda row: int(row["primary_bidir_callable_sites"]) >= 500,
        ),
        (
            "stringent_bidir_ge1000_no_high_divergence",
            "primary_bidir",
            lambda row: int(row["primary_bidir_callable_sites"]) >= 1000
            and no_high_divergence(row),
        ),
        (
            "corroborator_agree_bidir_ge500_no_high_divergence",
            "primary_bidir",
            lambda row: int(row["primary_bidir_callable_sites"]) >= 500
            and no_high_divergence(row)
            and row["boundary_corroborator_exact_state_agreement"] == "PASS",
        ),
        (
            "multi_outgroup_ge500_no_high_divergence",
            "multi_outgroup",
            lambda row: int(row["multi_outgroup_callable_sites"]) >= 500
            and no_high_divergence(row),
        ),
        (
            "Pquery_only_bias_control_ge500_no_high_divergence",
            "boundary_Pquery",
            lambda row: int(row["boundary_Pquery_callable_sites"]) >= 500
            and no_high_divergence(row),
        ),
    ]
    primary = [
        row
        for row in source
        if int(row["primary_bidir_callable_sites"]) >= 500
        and no_high_divergence(row)
    ]
    for age in ("time1", "time2", "time3"):
        definitions.append(
            (
                f"PRIMARY_{age}",
                "primary_bidir",
                lambda row, age=age: row in primary and row["age_bin"] == age,
            )
        )

    summaries = []
    for name, prefix, predicate in definitions:
        subset = [row for row in source if predicate(row)]
        if subset:
            summaries.append(summarize(name, subset, prefix, rng))
    write_tsv(
        args.pilot / "sequence_variation/pd_rate_statistical_summary.tsv",
        summaries,
    )

    terminal_source = read_tsv(
        args.pilot
        / "sequence_variation/postduplication_terminal_branch_metrics.tsv"
    )
    terminal_summaries = []
    for species in ("Alyrata", "Bstricta"):
        subset = [
            row
            for row in terminal_source
            if row["postduplication_species"] == species
            and int(row["joint_callable_sites"]) >= 500
            and row["endpoint_sensitivity_flags"] == "NONE"
        ]
        if subset:
            terminal_summaries.append(
                summarize_terminal(
                    f"{species}_terminal_ge500_no_high_divergence",
                    subset,
                    rng,
                )
            )
    write_tsv(
        args.pilot
        / "sequence_variation/postduplication_terminal_statistical_summary.tsv",
        terminal_summaries,
    )

    eligibility = []
    primary_ids = {row["event_id"] for row in primary}
    for row in source:
        reasons = []
        if int(row["primary_bidir_callable_sites"]) < 500:
            reasons.append("bidirectional_callable_lt_500")
        if not no_high_divergence(row):
            reasons.append(row["endpoint_sensitivity_flags"])
        eligibility.append(
            {
                "event_id": row["event_id"],
                "age_bin": row["age_bin"],
                "primary_endpoint_eligible": (
                    "PASS" if row["event_id"] in primary_ids else "FAIL"
                ),
                "exclusion_or_sensitivity_reason": (
                    ",".join(reasons) if reasons else "NONE"
                ),
                "primary_bidir_callable_sites": row[
                    "primary_bidir_callable_sites"
                ],
                "P_specific_rate": row["primary_bidir_P_specific_rate"],
                "D_specific_rate": row["primary_bidir_D_specific_rate"],
                "D_minus_P_rate": row["primary_bidir_D_minus_P_rate"],
            }
        )
    write_tsv(
        args.pilot / "sequence_variation/primary_endpoint_eligibility.tsv",
        eligibility,
    )

    core_qc = {
        row["event_id"]: row
        for row in read_tsv(args.pilot / "core/core_event_qc.tsv")
    }
    gap_rows = []
    for row in source:
        qc = core_qc[row["event_id"]]
        if qc["p_copy"] == "copy1":
            p_unique = int(qc["mate1_insertion_I_bp"])
            d_unique = int(qc["mate2_insertion_D_bp"])
            p_masked = int(qc["mate1_softmasked_S_bp"])
            d_masked = int(qc["mate2_softmasked_N_bp"])
        else:
            p_unique = int(qc["mate2_insertion_D_bp"])
            d_unique = int(qc["mate1_insertion_I_bp"])
            p_masked = int(qc["mate2_softmasked_N_bp"])
            d_masked = int(qc["mate1_softmasked_S_bp"])
        denominator = (
            int(qc["paired_M_core_bp"]) + p_unique + d_unique
        )
        gap_rows.append(
            {
                "event_id": row["event_id"],
                "age_bin": row["age_bin"],
                "primary_endpoint_eligible": (
                    "PASS" if row["event_id"] in primary_ids else "FAIL"
                ),
                "P_unique_unmasked_gap_bp": p_unique,
                "D_unique_unmasked_gap_bp": d_unique,
                "P_unique_gap_rate": f"{p_unique / denominator:.8f}",
                "D_unique_gap_rate": f"{d_unique / denominator:.8f}",
                "D_minus_P_unique_gap_rate": (
                    f"{(d_unique - p_unique) / denominator:.8f}"
                ),
                "P_softmasked_omission_bp": p_masked,
                "D_softmasked_omission_bp": d_masked,
                "relative_gap_denominator_bp": denominator,
                "interpretation_limit": (
                    "relative_copy_unique_sequence_not_polarized_insertion_or_deletion"
                ),
            }
        )
    write_tsv(
        args.pilot / "sequence_variation/pd_relative_gap_metrics.tsv",
        gap_rows,
    )

    primary_gaps = [
        row for row in gap_rows
        if row["primary_endpoint_eligible"] == "PASS"
    ]
    p_gap = np.array(
        [int(row["P_unique_unmasked_gap_bp"]) for row in primary_gaps]
    )
    d_gap = np.array(
        [int(row["D_unique_unmasked_gap_bp"]) for row in primary_gaps]
    )
    denominators = np.array(
        [int(row["relative_gap_denominator_bp"]) for row in primary_gaps]
    )
    gap_differences = d_gap / denominators - p_gap / denominators
    gap_non_ties = gap_differences[gap_differences != 0]
    gap_sign = stats.binomtest(
        int((gap_non_ties > 0).sum()), len(gap_non_ties), 0.5
    )
    gap_wilcoxon = stats.wilcoxon(
        d_gap / denominators,
        p_gap / denominators,
        zero_method="wilcox",
        method="auto",
    )
    gap_bootstrap_ratio = []
    gap_bootstrap_difference = []
    for _ in range(20_000):
        indices = rng.integers(0, len(primary_gaps), len(primary_gaps))
        p_sum = int(p_gap[indices].sum())
        d_sum = int(d_gap[indices].sum())
        if p_sum:
            gap_bootstrap_ratio.append(d_sum / p_sum)
        gap_bootstrap_difference.append(
            np.median(gap_differences[indices])
        )
    gap_ratio_ci = np.quantile(gap_bootstrap_ratio, [0.025, 0.975])
    gap_diff_ci = np.quantile(gap_bootstrap_difference, [0.025, 0.975])
    write_tsv(
        args.pilot / "sequence_variation/pd_relative_gap_summary.tsv",
        [
            {
                "analysis": "PRIMARY_relative_unmasked_gap_asymmetry",
                "events": len(primary_gaps),
                "P_unique_unmasked_gap_bp": int(p_gap.sum()),
                "D_unique_unmasked_gap_bp": int(d_gap.sum()),
                "D_to_P_gap_bp_ratio": f"{d_gap.sum() / p_gap.sum():.6f}",
                "ratio_bootstrap_CI95_low": f"{gap_ratio_ci[0]:.6f}",
                "ratio_bootstrap_CI95_high": f"{gap_ratio_ci[1]:.6f}",
                "median_D_minus_P_gap_rate": (
                    f"{np.median(gap_differences):.8f}"
                ),
                "median_difference_CI95_low": f"{gap_diff_ci[0]:.8f}",
                "median_difference_CI95_high": f"{gap_diff_ci[1]:.8f}",
                "events_D_gt_P": int((gap_differences > 0).sum()),
                "events_P_gt_D": int((gap_differences < 0).sum()),
                "events_tied": int((gap_differences == 0).sum()),
                "event_sign_test_p": f"{gap_sign.pvalue:.12g}",
                "paired_wilcoxon_p": f"{gap_wilcoxon.pvalue:.12g}",
                "interpretation_limit": (
                    "BISER_I_D_are_relative_copy_unique_bases;"
                    "ancestral_insertion_vs_deletion_direction_not_resolved"
                ),
            }
        ],
    )

    colors = {"time1": "#2b8cbe", "time2": "#f28e2b", "time3": "#7b3294"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for age in colors:
        subset = [row for row in primary if row["age_bin"] == age]
        axes[0].scatter(
            [float(row["primary_bidir_P_specific_rate"]) for row in subset],
            [float(row["primary_bidir_D_specific_rate"]) for row in subset],
            label=f"{age} (n={len(subset)})",
            color=colors[age],
            alpha=0.8,
            s=38,
        )
    limit = max(
        max(
            float(row["primary_bidir_P_specific_rate"]),
            float(row["primary_bidir_D_specific_rate"]),
        )
        for row in primary
    )
    axes[0].plot([0, limit], [0, limit], "--", color="black", linewidth=1)
    axes[0].set_xlabel("P-specific substitutions / callable site")
    axes[0].set_ylabel("D-specific substitutions / callable site")
    axes[0].legend(frameon=False)
    axes[0].set_title("Event-level polarized rates")

    age_data = [
        [
            float(row["primary_bidir_D_minus_P_rate"])
            for row in primary
            if row["age_bin"] == age
        ]
        for age in ("time1", "time2", "time3")
    ]
    box = axes[1].boxplot(
        age_data,
        tick_labels=["time1", "time2", "time3"],
        patch_artist=True,
        showfliers=True,
    )
    for patch, age in zip(box["boxes"], ("time1", "time2", "time3")):
        patch.set_facecolor(colors[age])
        patch.set_alpha(0.65)
    axes[1].axhline(0, linestyle="--", color="black", linewidth=1)
    axes[1].set_ylabel("D-specific rate − P-specific rate")
    axes[1].set_title("Within-event rate difference")
    fig.tight_layout()
    figure_dir = args.pilot / "sequence_variation/figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / "PD_substitution_rate_asymmetry.png", dpi=240)
    fig.savefig(figure_dir / "PD_substitution_rate_asymmetry.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    terminal_colors = {"Alyrata": "#4daf4a", "Bstricta": "#e41a1c"}
    terminal_limit = 0.0
    for species, color in terminal_colors.items():
        subset = [
            row
            for row in terminal_source
            if row["postduplication_species"] == species
            and int(row["joint_callable_sites"]) >= 500
            and row["endpoint_sensitivity_flags"] == "NONE"
        ]
        p_values = [
            float(row["P_terminal_mismatch_rate"]) for row in subset
        ]
        d_values = [
            float(row["D_terminal_mismatch_rate"]) for row in subset
        ]
        if p_values:
            terminal_limit = max(
                terminal_limit, max(p_values), max(d_values)
            )
        axes[0].scatter(
            p_values,
            d_values,
            label=f"{species} (n={len(subset)})",
            color=color,
            alpha=0.8,
            s=40,
        )
    axes[0].plot(
        [0, terminal_limit],
        [0, terminal_limit],
        "--",
        color="black",
        linewidth=1,
    )
    axes[0].set_xlabel("P terminal mismatch rate")
    axes[0].set_ylabel("D terminal mismatch rate")
    axes[0].set_title("Independent post-duplication comparison")
    axes[0].legend(frameon=False)

    p_gap_rates = p_gap / denominators
    d_gap_rates = d_gap / denominators
    gap_limit = max(float(p_gap_rates.max()), float(d_gap_rates.max()))
    axes[1].scatter(
        p_gap_rates,
        d_gap_rates,
        color="#666666",
        alpha=0.75,
        s=38,
    )
    axes[1].plot(
        [0, gap_limit], [0, gap_limit], "--", color="black", linewidth=1
    )
    axes[1].set_xlabel("P relative unique-gap rate")
    axes[1].set_ylabel("D relative unique-gap rate")
    axes[1].set_title("BISER relative gap asymmetry")
    fig.tight_layout()
    fig.savefig(
        figure_dir / "PD_terminal_and_gap_sensitivity.png", dpi=240
    )
    fig.savefig(figure_dir / "PD_terminal_and_gap_sensitivity.pdf")
    plt.close(fig)

    print(
        f"Primary endpoint: {len(primary)} events; "
        f"wrote {len(summaries)} statistical comparisons"
    )


if __name__ == "__main__":
    main()
