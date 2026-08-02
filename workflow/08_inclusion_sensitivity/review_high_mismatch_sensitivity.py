#!/usr/bin/env python3
"""Independent review of the >40% present-day P/D mismatch sensitivity tier.

This script consumes the already projected, bidirectionally exact P0-coordinate
metrics.  It does not realign sequences, relax the comparable-site definition,
or pool high-mismatch events with either the 44-event primary endpoint or the
55-event controlled sensitivity endpoint.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import platform
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/drmegd_high_mismatch_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats


SEED = 20260725
BOOTSTRAP_REPLICATES = 20_000


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mismatch_stratum(value: float) -> str:
    if not value > 40:
        raise ValueError(f"High-mismatch event has mismatch={value}")
    if value <= 50:
        return "(40,50]"
    if value <= 60:
        return "(50,60]"
    if value <= 70:
        return "(60,70]"
    return ">70"


def finite(value: float) -> str:
    return "NA" if not np.isfinite(value) else f"{value:.8f}"


def pvalue(value: float) -> str:
    return "NA" if not np.isfinite(value) else f"{value:.12g}"


def safe_ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else np.nan


def summarize(
    label: str,
    analysis_role: str,
    rows: list[dict[str, object]],
    events_total: int,
    rng: np.random.Generator,
) -> dict[str, object]:
    """Summarize events with >=1 exact bidirectional comparable site."""
    if not rows:
        return {
            "analysis_set": label,
            "analysis_role": analysis_role,
            "events_total": events_total,
            "events_analyzable": 0,
            "events_zero_comparable": events_total,
            "callable_sites": 0,
            "P_specific_changes": 0,
            "D_specific_changes": 0,
            "D_to_P_count_ratio": "NA",
            "D_to_P_ratio_bootstrap_CI95_low": "NA",
            "D_to_P_ratio_bootstrap_CI95_high": "NA",
            "median_event_D_minus_P_rate": "NA",
            "median_D_minus_P_rate_bootstrap_CI95_low": "NA",
            "median_D_minus_P_rate_bootstrap_CI95_high": "NA",
            "events_D_gt_P": 0,
            "events_P_gt_D": 0,
            "events_tied": 0,
            "event_sign_test_p": "NA",
            "paired_wilcoxon_statistic": "NA",
            "paired_wilcoxon_p": "NA",
        }

    p = np.array([int(row["P_specific_changes"]) for row in rows], dtype=float)
    d = np.array([int(row["D_specific_changes"]) for row in rows], dtype=float)
    callable_sites = np.array(
        [int(row["comparable_sites"]) for row in rows], dtype=float
    )
    differences = (d - p) / callable_sites
    samples = rng.integers(
        0, len(rows), size=(BOOTSTRAP_REPLICATES, len(rows))
    )
    boot_p = p[samples].sum(axis=1)
    boot_d = d[samples].sum(axis=1)
    boot_ratio = np.divide(
        boot_d,
        boot_p,
        out=np.full(BOOTSTRAP_REPLICATES, np.nan),
        where=boot_p != 0,
    )
    boot_median_difference = np.median(differences[samples], axis=1)
    non_ties = differences != 0
    greater = int((differences > 0).sum())
    sign_p = (
        stats.binomtest(greater, int(non_ties.sum()), 0.5).pvalue
        if non_ties.any()
        else np.nan
    )
    try:
        wilcoxon = stats.wilcoxon(
            differences,
            zero_method="wilcox",
            alternative="two-sided",
            method="auto",
        )
        wilcoxon_statistic = float(wilcoxon.statistic)
        wilcoxon_p = float(wilcoxon.pvalue)
    except ValueError:
        wilcoxon_statistic = np.nan
        wilcoxon_p = np.nan
    ratio_ci = np.nanquantile(boot_ratio, [0.025, 0.975])
    diff_ci = np.quantile(boot_median_difference, [0.025, 0.975])
    return {
        "analysis_set": label,
        "analysis_role": analysis_role,
        "events_total": events_total,
        "events_analyzable": len(rows),
        "events_zero_comparable": events_total - len(rows),
        "callable_sites": int(callable_sites.sum()),
        "P_specific_changes": int(p.sum()),
        "D_specific_changes": int(d.sum()),
        "D_to_P_count_ratio": finite(safe_ratio(d.sum(), p.sum())),
        "D_to_P_ratio_bootstrap_CI95_low": finite(float(ratio_ci[0])),
        "D_to_P_ratio_bootstrap_CI95_high": finite(float(ratio_ci[1])),
        "median_event_D_minus_P_rate": finite(float(np.median(differences))),
        "median_D_minus_P_rate_bootstrap_CI95_low": finite(float(diff_ci[0])),
        "median_D_minus_P_rate_bootstrap_CI95_high": finite(float(diff_ci[1])),
        "events_D_gt_P": greater,
        "events_P_gt_D": int((differences < 0).sum()),
        "events_tied": int((differences == 0).sum()),
        "event_sign_test_p": pvalue(float(sign_p)),
        "paired_wilcoxon_statistic": finite(wilcoxon_statistic),
        "paired_wilcoxon_p": pvalue(wilcoxon_p),
    }


def metric_row(
    event_id: str,
    age_bin: str,
    callable_sites: int,
    p_count: int,
    d_count: int,
) -> dict[str, object]:
    p_rate = safe_ratio(p_count, callable_sites)
    d_rate = safe_ratio(d_count, callable_sites)
    return {
        "event_id": event_id,
        "age_bin": age_bin,
        "comparable_sites": callable_sites,
        "P_specific_changes": p_count,
        "D_specific_changes": d_count,
        "P_specific_rate": p_rate,
        "D_specific_rate": d_rate,
        "D_minus_P_rate": d_rate - p_rate,
    }


def make_figure(
    detail: list[dict[str, object]],
    summary: list[dict[str, object]],
    output_dir: Path,
) -> None:
    strata = ["(40,50]", "(50,60]", "(60,70]", ">70"]
    colors = {
        "(40,50]": "#4c78a8",
        "(50,60]": "#72b7b2",
        "(60,70]": "#f2cf5b",
        ">70": "#e45756",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.1))

    ax = axes[0]
    rng = np.random.default_rng(SEED)
    for index, stratum in enumerate(strata):
        values = [
            float(row["D_minus_P_rate"])
            for row in detail
            if row["mismatch_stratum"] == stratum
            and row["included_in_event_level_statistics"] == "YES"
        ]
        if values:
            jitter = rng.uniform(-0.14, 0.14, len(values))
            ax.scatter(
                np.full(len(values), index) + jitter,
                values,
                s=36,
                alpha=0.82,
                color=colors[stratum],
                edgecolor="white",
                linewidth=0.45,
            )
            ax.plot(
                [index - 0.22, index + 0.22],
                [np.median(values), np.median(values)],
                color="black",
                linewidth=2,
            )
    ax.axhline(0, color="#555555", linestyle="--", linewidth=1)
    ax.set_xticks(range(len(strata)), strata)
    ax.set_xlabel("Present-day P/D mismatch stratum (%)")
    ax.set_ylabel("Event-level D − P substitution rate")
    ax.set_title("A  High-mismatch events only\n(exact bidirectional P0 coordinates)")

    ax = axes[1]
    selected_labels = [
        "primary_44",
        "controlled_55",
        "high_mismatch_all",
        "high_mismatch_(40,50]",
        "high_mismatch_(50,60]",
        "high_mismatch_(60,70]",
        "high_mismatch_>70",
    ]
    selected = {
        str(row["analysis_set"]): row
        for row in summary
        if row["analysis_set"] in selected_labels
    }
    ylabels = [
        "Primary (44)",
        "Controlled (55)",
        "High mismatch (32; 25 evaluable)",
        "(40,50]",
        "(50,60]",
        "(60,70]",
        ">70",
    ]
    estimates = np.array(
        [float(selected[label]["D_to_P_count_ratio"]) for label in selected_labels]
    )
    lows = np.array(
        [
            float(selected[label]["D_to_P_ratio_bootstrap_CI95_low"])
            for label in selected_labels
        ]
    )
    highs = np.array(
        [
            float(selected[label]["D_to_P_ratio_bootstrap_CI95_high"])
            for label in selected_labels
        ]
    )
    y = np.arange(len(selected_labels))
    point_colors = ["#1f4e79", "#4c78a8"] + [
        "#9c3d3d",
        colors["(40,50]"],
        colors["(50,60]"],
        colors["(60,70]"],
        colors[">70"],
    ]
    for index, color in enumerate(point_colors):
        ax.errorbar(
            estimates[index],
            y[index],
            xerr=np.array(
                [
                    [estimates[index] - lows[index]],
                    [highs[index] - estimates[index]],
                ]
            ),
            fmt="none",
            ecolor=color,
            elinewidth=1.6,
            capsize=3,
        )
    ax.scatter(estimates, y, color=point_colors, s=44, zorder=3)
    ax.axvline(1, color="#555555", linestyle="--", linewidth=1)
    ax.axhline(1.5, color="#999999", linestyle=":", linewidth=1)
    ax.set_yticks(y, ylabels)
    ax.invert_yaxis()
    ax.set_xlabel("Aggregate D/P substitution count ratio\n(95% event-bootstrap CI)")
    ax.set_title("B  Parallel reporting; no pooling")
    ax.grid(axis="x", color="#eeeeee")
    fig.tight_layout()
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures / "high_mismatch_independent_review.png", dpi=240)
    fig.savefig(figures / "high_mismatch_independent_review.pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
    )
    args = parser.parse_args()
    project = args.project_root.resolve()
    output_dir = Path(__file__).resolve().parent
    pilot = project / "07_pd_sequence_variation_pilot"
    inclusion = project / "08_event_inclusion_sensitivity"
    queue_path = pilot / "pilot/strict_primary_event_queue.tsv"
    metrics_path = pilot / "sequence_variation/event_pd_substitution_metrics.tsv"
    eligibility_path = pilot / "sequence_variation/primary_endpoint_eligibility.tsv"
    controlled_path = inclusion / "controlled_expansion_endpoint_events.tsv"

    queue = {row["event_id"]: row for row in read_tsv(queue_path)}
    metrics = {row["event_id"]: row for row in read_tsv(metrics_path)}
    if set(queue) != set(metrics):
        raise ValueError("Strict queue and projected metric event IDs differ")

    high_ids = {
        event_id
        for event_id, row in queue.items()
        if float(row["present_day_PD_mismatch_pct_callable"]) > 40
    }
    if len(high_ids) != 32:
        raise ValueError(f"Expected 32 high-mismatch events, found {len(high_ids)}")

    detail: list[dict[str, object]] = []
    high_stats: list[dict[str, object]] = []
    for event_id in sorted(
        high_ids,
        key=lambda value: (
            float(queue[value]["present_day_PD_mismatch_pct_callable"]),
            value,
        ),
    ):
        queue_row = queue[event_id]
        metric = metrics[event_id]
        callable_sites = int(metric["primary_bidir_callable_sites"])
        p_count = int(metric["primary_bidir_P_specific"])
        d_count = int(metric["primary_bidir_D_specific"])
        base = metric_row(
            event_id, metric["age_bin"], callable_sites, p_count, d_count
        )
        mismatch = float(queue_row["present_day_PD_mismatch_pct_callable"])
        stratum = mismatch_stratum(mismatch)
        analyzable = callable_sites > 0
        detail.append(
            {
                "event_id": event_id,
                "age_bin": metric["age_bin"],
                "present_day_PD_mismatch_bp": queue_row[
                    "present_day_PD_mismatch_bp"
                ],
                "present_day_PD_comparable_bp": queue_row[
                    "jointly_callable_uppercase_acgt_bp"
                ],
                "present_day_PD_mismatch_pct": f"{mismatch:.6f}",
                "mismatch_stratum": stratum,
                "P0_projection_rule": "bidirectional_same_coordinate",
                "comparable_site_rule": "P0,P,D_all_ACGT_at_same_projected_coordinate",
                "primary_bidir_comparable_sites": callable_sites,
                "P_specific_changes": p_count,
                "D_specific_changes": d_count,
                "P_specific_rate": (
                    finite(float(base["P_specific_rate"])) if analyzable else "NA"
                ),
                "D_specific_rate": (
                    finite(float(base["D_specific_rate"])) if analyzable else "NA"
                ),
                "D_minus_P_rate": (
                    finite(float(base["D_minus_P_rate"])) if analyzable else "NA"
                ),
                "included_in_event_level_statistics": "YES" if analyzable else "NO",
                "nonanalysis_reason": (
                    "NONE" if analyzable else "zero_exact_bidirectional_comparable_sites"
                ),
            }
        )
        if analyzable:
            base["mismatch_stratum"] = stratum
            high_stats.append(base)

    primary_ids = {
        row["event_id"]
        for row in read_tsv(eligibility_path)
        if row["primary_endpoint_eligible"] == "PASS"
    }
    if len(primary_ids) != 44:
        raise ValueError(f"Expected 44 primary events, found {len(primary_ids)}")
    controlled_source = read_tsv(controlled_path)
    controlled_ids = {row["event_id"] for row in controlled_source}
    if len(controlled_ids) != 55:
        raise ValueError(f"Expected 55 controlled events, found {len(controlled_ids)}")
    if high_ids & primary_ids or high_ids & controlled_ids:
        raise ValueError("High-mismatch tier overlaps a reference endpoint set")

    primary_rows = [
        metric_row(
            event_id,
            metrics[event_id]["age_bin"],
            int(metrics[event_id]["primary_bidir_callable_sites"]),
            int(metrics[event_id]["primary_bidir_P_specific"]),
            int(metrics[event_id]["primary_bidir_D_specific"]),
        )
        for event_id in sorted(primary_ids)
    ]
    controlled_rows = [
        metric_row(
            row["event_id"],
            row["age_bin"],
            int(row["primary_bidir_callable_sites"]),
            int(row["P_specific_changes"]),
            int(row["D_specific_changes"]),
        )
        for row in controlled_source
    ]

    rng = np.random.default_rng(SEED)
    summary = [
        summarize("primary_44", "REFERENCE_PRIMARY_NOT_POOLED", primary_rows, 44, rng),
        summarize(
            "controlled_55",
            "REFERENCE_CONTROLLED_NOT_POOLED",
            controlled_rows,
            55,
            rng,
        ),
        summarize(
            "high_mismatch_all",
            "INDEPENDENT_HIGH_MISMATCH_SENSITIVITY",
            high_stats,
            32,
            rng,
        ),
    ]
    for stratum in ("(40,50]", "(50,60]", "(60,70]", ">70"):
        stratum_rows = [
            row for row in high_stats if row["mismatch_stratum"] == stratum
        ]
        stratum_total = sum(
            row["mismatch_stratum"] == stratum for row in detail
        )
        summary.append(
            summarize(
                f"high_mismatch_{stratum}",
                "INDEPENDENT_HIGH_MISMATCH_STRATUM",
                stratum_rows,
                stratum_total,
                rng,
            )
        )

    write_tsv(output_dir / "high_mismatch_event_detail.tsv", detail)
    write_tsv(output_dir / "parallel_statistical_summary.tsv", summary)
    make_figure(detail, summary, output_dir)

    manifest = [
        {
            "input": str(path.relative_to(project)),
            "sha256": sha256(path),
        }
        for path in (queue_path, metrics_path, eligibility_path, controlled_path)
    ]
    write_tsv(output_dir / "input_manifest.tsv", manifest)
    log = (
        "# Reproduction log\n\n"
        "Run from this directory with:\n\n"
        "```bash\n"
        "python3 review_high_mismatch_sensitivity.py\n"
        "```\n\n"
        f"- Random seed: `{SEED}`\n"
        f"- Event-bootstrap replicates: `{BOOTSTRAP_REPLICATES}`\n"
        f"- Python: `{platform.python_version()}`\n"
        f"- NumPy: `{np.__version__}`\n"
        f"- SciPy: `{scipy.__version__}`\n"
        f"- Matplotlib: `{matplotlib.__version__}`\n"
        "- Statistical unit: event\n"
        "- Bootstrap: event resampling with replacement within each reported set\n"
        "- Sign test: exact two-sided binomial test after removing D−P rate ties\n"
        "- Wilcoxon: paired two-sided signed-rank test; zero differences omitted\n"
        "- High-mismatch events are never unioned with the 44- or 55-event sets\n"
    )
    (output_dir / "RUN_LOG.md").write_text(log)


if __name__ == "__main__":
    main()
