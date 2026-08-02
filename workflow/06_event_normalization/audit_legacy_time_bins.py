#!/usr/bin/env python3
"""Audit the hard-coded five-node SD age classification used in the legacy workflow.

This script is deliberately read-only with respect to the legacy results.  It
reproduces colink2time_point.sh, records patterns that the shell workflow drops,
and flags patterns that are not monotonic under a simple species-tree model.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def classify_legacy(states: tuple[int, ...]) -> tuple[str | None, int | None]:
    """Return the legacy time bin and orientation (0/1/2), or (None, None).

    State meanings inherited from classify_links.pl:
      0 neither SD mate covered; 1 mate 1; 2 mate 2; 3 both mates.
    """
    if len(states) != 5:
        raise ValueError(f"expected five node states, got {states!r}")

    # time1: no node has both copies; one-sided evidence may only support one mate.
    if 3 not in states:
        sides = {state for state in states if state in (1, 2)}
        if len(sides) <= 1:
            return "time1", next(iter(sides), 0)
        return None, None

    # time2..time5: the current node has both copies and all more distant nodes
    # contain only zero or evidence for the same single mate.
    for current_index in range(4):
        if states[current_index] != 3:
            continue
        older = states[current_index + 1 :]
        if 3 in older:
            continue
        sides = {state for state in older if state in (1, 2)}
        if len(sides) <= 1:
            return f"time{current_index + 2}", next(iter(sides), 0)

    # time6: the most distant node has both copies.  The old workflow does not
    # attempt P/D orientation because it has no still older outgroup.
    if states[-1] == 3:
        return "time6", 0
    return None, None


def monotonicity_flag(states: tuple[int, ...]) -> str:
    """Flag hard-state patterns that contradict a simple gain-only ladder."""
    both = [i for i, state in enumerate(states) if state == 3]
    if not both:
        return "no_both_copy_node"
    farthest_both = max(both)
    if any(states[i] != 3 for i in range(farthest_both)):
        return "nearer_node_not_both"
    older_sides = {state for state in states[farthest_both + 1 :] if state in (1, 2)}
    if len(older_sides) > 1:
        return "older_nodes_switch_supported_mate"
    return "simple_monotonic"


def read_patterns(path: Path) -> Counter[tuple[int, ...]]:
    patterns: Counter[tuple[int, ...]] = Counter()
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_number}: fewer than eight tab fields")
            try:
                states = tuple(int(value) for value in fields[3:8])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid node states") from exc
            patterns[states] += 1
    return patterns


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--patterns", required=True, type=Path)
    args = parser.parse_args()

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.patterns.parent.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []
    pattern_rows: list[dict[str, object]] = []

    for path in args.inputs:
        patterns = read_patterns(path)
        bin_counts: Counter[str] = Counter()
        orientation_counts: Counter[int] = Counter()
        flag_counts: Counter[str] = Counter()
        classified = 0

        for states, count in sorted(patterns.items()):
            time_bin, orientation = classify_legacy(states)
            flag = monotonicity_flag(states)
            flag_counts[flag] += count
            if time_bin is not None and orientation is not None:
                classified += count
                bin_counts[time_bin] += count
                orientation_counts[orientation] += count
            pattern_rows.append(
                {
                    "dataset": path.parent.name,
                    "source_file": str(path),
                    "pattern": "".join(map(str, states)),
                    "count": count,
                    "legacy_time_bin": time_bin or "UNCLASSIFIED",
                    "legacy_orientation": (
                        {0: "unpolarized", 1: "mate1_is_P", 2: "mate2_is_P"}.get(
                            orientation, "UNCLASSIFIED"
                        )
                    ),
                    "monotonicity": flag,
                }
            )

        total = sum(patterns.values())
        summary_rows.append(
            {
                "dataset": path.parent.name,
                "source_file": str(path),
                "total_pairs": total,
                "unique_patterns": len(patterns),
                "legacy_classified": classified,
                "legacy_unclassified": total - classified,
                "legacy_classified_pct": f"{100 * classified / total:.3f}" if total else "NA",
                "legacy_polarized": orientation_counts[1] + orientation_counts[2],
                "legacy_unpolarized": orientation_counts[0],
                **{f"time{i}_pairs": bin_counts[f"time{i}"] for i in range(1, 7)},
                "simple_monotonic": flag_counts["simple_monotonic"],
                "nearer_node_not_both": flag_counts["nearer_node_not_both"],
                "older_nodes_switch_supported_mate": flag_counts[
                    "older_nodes_switch_supported_mate"
                ],
                "no_both_copy_node": flag_counts["no_both_copy_node"],
            }
        )

    summary_fields = list(summary_rows[0]) if summary_rows else []
    with args.summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    pattern_fields = list(pattern_rows[0]) if pattern_rows else []
    with args.patterns.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pattern_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(pattern_rows)


if __name__ == "__main__":
    main()
