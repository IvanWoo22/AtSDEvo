#!/usr/bin/env python3
"""Summarize nested, controlled expansions of the P/D substitution endpoint."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, wilcoxon


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    root = Path(__file__).resolve().parent
    project = root.parent
    primary = project / "07_pd_sequence_variation_pilot"
    pilot = root / "core_500"

    def metrics(directory: Path) -> dict[str, dict[str, str]]:
        return {
            row["event_id"]: row
            for row in read_tsv(directory / "event_pd_substitution_metrics.tsv")
        }

    old_metrics = metrics(primary / "sequence_variation")
    old_ids = {
        row["event_id"]
        for row in read_tsv(
            primary / "sequence_variation/primary_endpoint_eligibility.tsv"
        )
        if row["primary_endpoint_eligible"] == "PASS"
    }
    strict_queue = {
        row["event_id"]: row
        for row in read_tsv(pilot / "pilot/strict_primary_event_queue.tsv")
    }
    strict_metrics = metrics(pilot / "sequence_variation")
    partial_queue = {
        row["event_id"]: row
        for row in read_tsv(pilot / "pilot/partial_postdup_event_queue.tsv")
    }
    partial_metrics = metrics(pilot / "sequence_variation_partial_postdup")
    fallback_queue = {
        row["event_id"]: row
        for row in read_tsv(root / "deeper_P0_sensitivity_queue.tsv")
    }
    fallback_metrics = metrics(pilot / "sequence_variation_deeper_P0")

    def eligible(queue_row: dict[str, str], metric: dict[str, str]) -> bool:
        return (
            int(metric["primary_bidir_callable_sites"]) >= 500
            and float(queue_row["present_day_PD_mismatch_pct_callable"]) <= 40
        )

    strict_ids = {
        event_id
        for event_id, row in strict_queue.items()
        if eligible(row, strict_metrics[event_id])
    }
    partial_added = {
        event_id
        for event_id, row in partial_queue.items()
        if event_id not in strict_queue
        and eligible(row, partial_metrics[event_id])
    }
    fallback_added = {
        event_id
        for event_id, row in fallback_queue.items()
        if event_id not in strict_queue
        and eligible(row, fallback_metrics[event_id])
    }

    event_source: dict[str, tuple[str, dict[str, str]]] = {
        event_id: ("core_500_strict", strict_metrics[event_id])
        for event_id in strict_ids
    }
    event_source.update(
        {
            event_id: ("partial_postduplication", partial_metrics[event_id])
            for event_id in partial_added
        }
    )
    event_source.update(
        {
            event_id: ("deeper_P0_fallback", fallback_metrics[event_id])
            for event_id in fallback_added
        }
    )
    event_rows = []
    for event_id, (tier, row) in sorted(event_source.items()):
        event_rows.append(
            {
                "event_id": event_id,
                "age_bin": row["age_bin"],
                "admission_tier": tier,
                "primary_bidir_callable_sites": row[
                    "primary_bidir_callable_sites"
                ],
                "P_specific_changes": row["primary_bidir_P_specific"],
                "D_specific_changes": row["primary_bidir_D_specific"],
                "P_specific_rate": row["primary_bidir_P_specific_rate"],
                "D_specific_rate": row["primary_bidir_D_specific_rate"],
                "D_minus_P_rate": row["primary_bidir_D_minus_P_rate"],
            }
        )
    write_tsv(root / "controlled_expansion_endpoint_events.tsv", event_rows)

    nested = [
        ("original_1kb_primary", old_ids, old_metrics),
        ("core_500_strict_mapping", strict_ids, strict_metrics),
        (
            "core_500_plus_partial_postdup",
            strict_ids | partial_added,
            {**strict_metrics, **partial_metrics},
        ),
        (
            "core_500_plus_deeper_P0",
            strict_ids | fallback_added,
            {**strict_metrics, **fallback_metrics},
        ),
        (
            "controlled_union",
            set(event_source),
            {event_id: row for event_id, (_, row) in event_source.items()},
        ),
    ]
    rng = np.random.default_rng(20260724)
    summary = []
    for label, event_ids, source in nested:
        rows = [source[event_id] for event_id in sorted(event_ids)]
        p = np.array(
            [int(row["primary_bidir_P_specific"]) for row in rows],
            dtype=float,
        )
        d = np.array(
            [int(row["primary_bidir_D_specific"]) for row in rows],
            dtype=float,
        )
        rates = np.array(
            [float(row["primary_bidir_D_minus_P_rate"]) for row in rows]
        )
        samples = rng.integers(0, len(rows), size=(20_000, len(rows)))
        boot_p = p[samples].sum(axis=1)
        boot_d = d[samples].sum(axis=1)
        boot_ratio = np.divide(
            boot_d,
            boot_p,
            out=np.full_like(boot_d, np.nan),
            where=boot_p != 0,
        )
        non_ties = int(np.sum(d != p))
        greater = int(np.sum(d > p))
        wilcox = wilcoxon(rates, zero_method="wilcox", alternative="two-sided")
        summary.append(
            {
                "analysis_set": label,
                "events": len(rows),
                "time1": sum(row["age_bin"] == "time1" for row in rows),
                "time2": sum(row["age_bin"] == "time2" for row in rows),
                "time3": sum(row["age_bin"] == "time3" for row in rows),
                "callable_sites": sum(
                    int(row["primary_bidir_callable_sites"]) for row in rows
                ),
                "P_specific_changes": int(p.sum()),
                "D_specific_changes": int(d.sum()),
                "D_to_P_count_ratio": f"{d.sum() / p.sum():.6f}",
                "bootstrap_ratio_CI95_low": (
                    f"{np.nanquantile(boot_ratio, 0.025):.6f}"
                ),
                "bootstrap_ratio_CI95_high": (
                    f"{np.nanquantile(boot_ratio, 0.975):.6f}"
                ),
                "events_D_gt_P": greater,
                "events_P_gt_D": int(np.sum(d < p)),
                "events_tied": int(np.sum(d == p)),
                "event_sign_test_p": (
                    f"{binomtest(greater, non_ties, 0.5).pvalue:.12g}"
                ),
                "paired_wilcoxon_p": f"{wilcox.pvalue:.12g}",
            }
        )
    write_tsv(root / "controlled_expansion_statistical_summary.tsv", summary)


if __name__ == "__main__":
    main()
