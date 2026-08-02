#!/usr/bin/env python3
"""Compare exact calls and event-level P/D direction between two de novo runs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row["event_id"],
        row["unique_role"],
        row["chrom"],
        row["start_0based"],
        row["end_0based"],
        row["parsimonious_direction"],
    )


def event_counts(rows: list[dict[str, str]]) -> dict[str, tuple[int, int]]:
    result: dict[str, list[int]] = {}
    for row in rows:
        values = result.setdefault(row["event_id"], [0, 0])
        values[0 if row["parsimonious_direction"].startswith("P_") else 1] += 1
    return {event: (values[0], values[1]) for event, values in result.items()}


def direction(values: tuple[int, int]) -> str:
    return "D" if values[1] > values[0] else "P" if values[0] > values[1] else "TIE"


def exact_binomial(k: int, n: int) -> float:
    observed = math.comb(n, k)
    return sum(
        math.comb(n, value)
        for value in range(n + 1)
        if math.comb(n, value) <= observed
    ) / (2**n)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument("--sensitivity", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    main_rows = [row for row in read(args.main) if row["evidence_tier"] == "PRIMARY"]
    sensitivity_rows = [
        row for row in read(args.sensitivity) if row["evidence_tier"] == "PRIMARY"
    ]
    main_calls = {key(row): row for row in main_rows}
    sensitivity_calls = {key(row): row for row in sensitivity_rows}
    event_age = {
        row["event_id"]: row["age_bin"] for row in main_rows + sensitivity_rows
    }
    call_rows = []
    for call in sorted(set(main_calls) | set(sensitivity_calls)):
        source = main_calls.get(call) or sensitivity_calls[call]
        call_rows.append(
            {
                "event_id": call[0],
                "unique_role": call[1],
                "chrom": call[2],
                "start_0based": call[3],
                "end_0based": call[4],
                "direction": call[5],
                "fragment_bp": source["fragment_bp"],
                "main_present": "YES" if call in main_calls else "NO",
                "sensitivity_present": "YES" if call in sensitivity_calls else "NO",
                "exact_call_robust": (
                    "YES" if call in main_calls and call in sensitivity_calls else "NO"
                ),
            }
        )
    write(args.output / "window_parameter_call_robustness.tsv", call_rows)
    main_events = event_counts(main_rows)
    sensitivity_events = event_counts(sensitivity_rows)
    event_rows = []
    for event in sorted(set(main_events) | set(sensitivity_events)):
        m = main_events.get(event, (0, 0))
        s = sensitivity_events.get(event, (0, 0))
        event_rows.append(
            {
                "event_id": event,
                "age_bin": event_age[event],
                "main_P": m[0],
                "main_D": m[1],
                "main_direction": direction(m),
                "sensitivity_P": s[0],
                "sensitivity_D": s[1],
                "sensitivity_direction": direction(s),
                "D_direction_robust": (
                    "YES" if direction(m) == direction(s) == "D" else "NO"
                ),
                "P_direction_robust": (
                    "YES" if direction(m) == direction(s) == "P" else "NO"
                ),
            }
        )
    write(args.output / "window_parameter_event_robustness.tsv", event_rows)
    exact = sum(row["exact_call_robust"] == "YES" for row in call_rows)
    robust_d = sum(row["D_direction_robust"] == "YES" for row in event_rows)
    robust_p = sum(row["P_direction_robust"] == "YES" for row in event_rows)
    write(
        args.output / "window_parameter_robustness_summary.tsv",
        [
            {
                "main_primary_calls": len(main_rows),
                "sensitivity_primary_calls": len(sensitivity_rows),
                "exact_shared_calls": exact,
                "events_D_direction_in_both": robust_d,
                "events_P_direction_in_both": robust_p,
                "robust_event_sign_test_two_sided_p": (
                    f"{exact_binomial(robust_d, robust_d + robust_p):.8g}"
                    if robust_d + robust_p
                    else "NA"
                ),
            }
        ],
    )
    age_rows = []
    for age in ("time1", "time2", "time3"):
        subset = [row for row in event_rows if row["age_bin"] == age]
        d_count = sum(row["D_direction_robust"] == "YES" for row in subset)
        p_count = sum(row["P_direction_robust"] == "YES" for row in subset)
        age_rows.append(
            {
                "age_bin": age,
                "events_D_direction_in_both": d_count,
                "events_P_direction_in_both": p_count,
                "robust_non_tied_events": d_count + p_count,
                "exact_two_sided_sign_test_p": (
                    f"{exact_binomial(d_count, d_count + p_count):.8g}"
                    if d_count + p_count
                    else "NA"
                ),
            }
        )
    write(args.output / "window_parameter_age_robustness.tsv", age_rows)


if __name__ == "__main__":
    main()
