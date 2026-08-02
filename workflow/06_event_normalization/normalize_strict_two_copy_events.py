#!/usr/bin/env python3
"""Build a conservative threshold-stable two-copy event catalogue from BISER calls."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


class DSU:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.size = [1] * size

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        left, right = self.find(left), self.find(right)
        if left == right:
            return
        if self.size[left] < self.size[right]:
            left, right = right, left
        self.parent[right] = left
        self.size[left] += self.size[right]


def reciprocal_overlap(
    left: tuple[str, int, int], right: tuple[str, int, int]
) -> float:
    if left[0] != right[0]:
        return 0.0
    overlap = max(0, min(left[2], right[2]) - max(left[1], right[1]))
    if not overlap:
        return 0.0
    return min(overlap / (left[2] - left[1]), overlap / (right[2] - right[1]))


def cluster_arms(
    arms: list[tuple[str, int, int]], threshold: float
) -> tuple[list[int], set[int]]:
    dsu = DSU(len(arms))
    by_chrom: dict[str, list[int]] = defaultdict(list)
    for index, arm in enumerate(arms):
        by_chrom[arm[0]].append(index)
    for indices in by_chrom.values():
        indices.sort(key=lambda index: arms[index][1])
        active: list[int] = []
        for index in indices:
            current = arms[index]
            active = [
                other for other in active if arms[other][2] > current[1]
            ]
            for other in active:
                if reciprocal_overlap(current, arms[other]) >= threshold:
                    dsu.union(index, other)
            active.append(index)
    roots = [dsu.find(index) for index in range(len(arms))]
    remap = {root: number for number, root in enumerate(sorted(set(roots)), 1)}
    locus_ids = [remap[root] for root in roots]
    members: dict[int, list[int]] = defaultdict(list)
    for index, locus_id in enumerate(locus_ids):
        members[locus_id].append(index)
    complete_link_loci = {
        locus_id
        for locus_id, indices in members.items()
        if all(
            reciprocal_overlap(arms[left], arms[right]) >= threshold
            for offset, left in enumerate(indices)
            for right in indices[offset + 1 :]
        )
    }
    return locus_ids, complete_link_loci


def strict_network_calls(
    locus_ids: list[int], relative_orientations: list[str]
) -> tuple[set[int], dict[int, tuple[int, int]], dict[int, int]]:
    edges: dict[tuple[int, int], list[int]] = defaultdict(list)
    adjacency: dict[int, set[int]] = defaultdict(set)
    for call_index, orientation in enumerate(relative_orientations):
        left, right = locus_ids[2 * call_index], locus_ids[2 * call_index + 1]
        edge = tuple(sorted((left, right)))
        edges[edge].append(call_index)
        adjacency[left].add(right)
        adjacency[right].add(left)

    all_loci = set(locus_ids)
    seen = set()
    strict_calls = set()
    call_edge = {}
    call_component_size = {}
    for locus in all_loci:
        if locus in seen:
            continue
        stack = [locus]
        seen.add(locus)
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        component_set = set(component)
        component_edges = [
            edge for edge in edges if edge[0] in component_set and edge[1] in component_set
        ]
        component_call_indices = [
            index for edge in component_edges for index in edges[edge]
        ]
        for index in component_call_indices:
            call_component_size[index] = len(component)
        if (
            len(component) == 2
            and len(component_edges) == 1
            and component_edges[0][0] != component_edges[0][1]
        ):
            orientations = {
                relative_orientations[index] for index in component_call_indices
            }
            if len(orientations) == 1:
                strict_calls.update(component_call_indices)
                for index in component_call_indices:
                    call_edge[index] = component_edges[0]
    return strict_calls, call_edge, call_component_size


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--biser", required=True, type=Path)
    parser.add_argument(
        "--evidence",
        type=Path,
        help=(
            "Optional call-level age/P-D evidence. Omit for an age-blind "
            "call-to-event catalogue."
        ),
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--network-policy",
        choices=("stable", "0.5", "0.8"),
        default="stable",
        help="Default retains calls strict at both 0.5 and 0.8.",
    )
    parser.add_argument("--minimum-length", type=int, default=1000)
    args = parser.parse_args()

    biser_rows = []
    arms = []
    orientations = []
    with args.biser.open() as handle:
        for call_id, line in enumerate(handle, 1):
            fields = line.rstrip().split("\t")
            left = (fields[0], *sorted((int(fields[1]), int(fields[2]))))
            right = (fields[3], *sorted((int(fields[4]), int(fields[5]))))
            arms.extend((left, right))
            orientations.append("same" if fields[8] == fields[9] else "opposite")
            biser_rows.append(
                {
                    "call_id": call_id,
                    "copy1_chrom": left[0],
                    "copy1_start": left[1],
                    "copy1_end": left[2],
                    "copy2_chrom": right[0],
                    "copy2_start": right[1],
                    "copy2_end": right[2],
                    "biser_error": float(fields[7]),
                    "max_mate_length_bp": int(fields[10]),
                    "alignment_span_bp": int(fields[11]),
                }
            )

    evidence = {}
    if args.evidence:
        with args.evidence.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                evidence[int(row["call_id"])] = row
        if set(evidence) != set(range(1, len(biser_rows) + 1)):
            raise SystemExit("Call IDs in evidence and BISER output do not match")

    threshold_results = {}
    for threshold in (0.5, 0.8):
        loci, complete_link_loci = cluster_arms(arms, threshold)
        strict, edges, component_sizes = strict_network_calls(loci, orientations)
        complete_link_calls = {
            index
            for index in range(len(biser_rows))
            if loci[2 * index] in complete_link_loci
            and loci[2 * index + 1] in complete_link_loci
        }
        strict &= complete_link_calls
        threshold_results[threshold] = (
            loci,
            strict,
            edges,
            component_sizes,
            complete_link_calls,
        )

    loci80, strict80, edges80, sizes80, complete80 = threshold_results[0.8]
    loci50, strict50, edges50, sizes50, complete50 = threshold_results[0.5]
    stable_calls = strict80 & strict50
    if args.network_policy == "stable":
        selected_calls = stable_calls
        selected_loci = loci80
        selected_edges = edges80
    elif args.network_policy == "0.5":
        selected_calls = strict50
        selected_loci = loci50
        selected_edges = edges50
    else:
        selected_calls = strict80
        selected_loci = loci80
        selected_edges = edges80
    grouped: dict[tuple[int, int], list[int]] = defaultdict(list)
    for call_index in selected_calls:
        grouped[selected_edges[call_index]].append(call_index)

    event_rows = []
    membership_rows = []
    reason_counts = Counter()
    for event_number, (edge, indices) in enumerate(
        sorted(grouped.items(), key=lambda item: min(item[1])), 1
    ):
        indices.sort()
        calls = [biser_rows[index] for index in indices]
        classifications = (
            {
                (
                    evidence[call["call_id"]]["strict_status"],
                    evidence[call["call_id"]]["strict_age_bin"],
                    evidence[call["call_id"]]["provisional_p_copy"],
                )
                for call in calls
            }
            if evidence
            else set()
        )
        if not evidence:
            status, age_bin, p_copy = (
                "NOT_EVALUATED",
                "event_first_pending",
                "NA",
            )
        elif len(classifications) == 1:
            status, age_bin, p_copy = next(iter(classifications))
        else:
            status, age_bin, p_copy = "EXCLUDE", "source_call_disagreement", "NA"
        representative = max(
            calls,
            key=lambda row: (
                min(
                    row["max_mate_length_bp"],
                    row["alignment_span_bp"],
                ),
                -row["biser_error"],
                -row["call_id"],
            ),
        )
        representative_index = int(representative["call_id"]) - 1
        representative_copy1_bp = (
            representative["copy1_end"] - representative["copy1_start"]
        )
        representative_copy2_bp = (
            representative["copy2_end"] - representative["copy2_start"]
        )
        if status == "PASS" and (
            representative_copy1_bp < args.minimum_length
            or representative_copy2_bp < args.minimum_length
            or representative["max_mate_length_bp"] < args.minimum_length
            or representative["alignment_span_bp"] < args.minimum_length
        ):
            status = "EXCLUDE"
            age_bin = (
                "below_1kb_strict_sd_core"
                if args.minimum_length == 1000
                else "below_minimum_length_strict_sd_core"
            )
            p_copy = "NA"
        if status != "PASS":
            reason_counts[age_bin] += 1

        if status == "PASS":
            p_arm = 0 if p_copy == "copy1" else 1
            d_arm = 1 - p_arm
            p_locus = selected_loci[2 * representative_index + p_arm]
            d_locus = selected_loci[2 * representative_index + d_arm]
            p_chrom = representative[f"copy{p_arm + 1}_chrom"]
            p_start = representative[f"copy{p_arm + 1}_start"]
            p_end = representative[f"copy{p_arm + 1}_end"]
            d_chrom = representative[f"copy{d_arm + 1}_chrom"]
            d_start = representative[f"copy{d_arm + 1}_start"]
            d_end = representative[f"copy{d_arm + 1}_end"]
        else:
            p_locus = d_locus = "NA"
            p_chrom = p_start = p_end = "NA"
            d_chrom = d_start = d_end = "NA"
        event_id = f"SDJGI{event_number:05d}"
        event_rows.append(
            {
                "event_id": event_id,
                "source_call_count": len(calls),
                "source_call_ids": ",".join(str(row["call_id"]) for row in calls),
                "representative_call_id": representative["call_id"],
                "copy_locus_1": edge[0],
                "copy_locus_2": edge[1],
                "representative_copy1_chrom": representative["copy1_chrom"],
                "representative_copy1_start": representative["copy1_start"],
                "representative_copy1_end": representative["copy1_end"],
                "representative_copy2_chrom": representative["copy2_chrom"],
                "representative_copy2_start": representative["copy2_start"],
                "representative_copy2_end": representative["copy2_end"],
                "representative_copy1_bp": representative_copy1_bp,
                "representative_copy2_bp": representative_copy2_bp,
                "representative_max_mate_length_bp": representative[
                    "max_mate_length_bp"
                ],
                "representative_alignment_span_bp": representative[
                    "alignment_span_bp"
                ],
                "strict_minimum_length_core": (
                    "PASS"
                    if representative_copy1_bp >= args.minimum_length
                    and representative_copy2_bp >= args.minimum_length
                    and representative["max_mate_length_bp"] >= args.minimum_length
                    and representative["alignment_span_bp"] >= args.minimum_length
                    else "FAIL"
                ),
                "strict_minimum_1kb_core": (
                    "PASS"
                    if representative_copy1_bp >= 1000
                    and representative_copy2_bp >= 1000
                    and representative["max_mate_length_bp"] >= 1000
                    and representative["alignment_span_bp"] >= 1000
                    else "FAIL"
                ),
                "copy1_envelope_start": min(row["copy1_start"] for row in calls),
                "copy1_envelope_end": max(row["copy1_end"] for row in calls),
                "copy2_envelope_start": min(row["copy2_start"] for row in calls),
                "copy2_envelope_end": max(row["copy2_end"] for row in calls),
                "relative_orientation": orientations[indices[0]],
                "network_stable_at_reciprocal_overlap_0.5_and_0.8": (
                    "PASS"
                    if all(index in stable_calls for index in indices)
                    else "FAIL"
                ),
                "network_selection_policy": args.network_policy,
                "minimum_length_bp": args.minimum_length,
                "strict_pd_status": status,
                "strict_age_bin": age_bin,
                "provisional_p_copy": p_copy,
                "provisional_p_locus": p_locus,
                "provisional_d_locus": d_locus,
                "provisional_p_chrom": p_chrom,
                "provisional_p_start": p_start,
                "provisional_p_end": p_end,
                "provisional_d_chrom": d_chrom,
                "provisional_d_start": d_start,
                "provisional_d_end": d_end,
            }
        )
        for index in indices:
            membership_rows.append(
                {
                    "event_id": event_id,
                    "call_id": biser_rows[index]["call_id"],
                    "locus80_copy1": loci80[2 * index],
                    "locus80_copy2": loci80[2 * index + 1],
                    "component_loci_at_0.8": sizes80[index],
                    "component_loci_at_0.5": sizes50[index],
                    "strict_network_at_0.8": "PASS" if index in strict80 else "FAIL",
                    "strict_network_at_0.5": "PASS" if index in strict50 else "FAIL",
                    "complete_link_arms_at_0.8": (
                        "PASS" if index in complete80 else "FAIL"
                    ),
                    "complete_link_arms_at_0.5": (
                        "PASS" if index in complete50 else "FAIL"
                    ),
                }
            )

    call_tiers = []
    for index, row in enumerate(biser_rows):
        tier = (
            "strict_stable_two_copy"
            if index in stable_calls
            else "strict_at_0.8_only"
            if index in strict80
            else "complex_or_self_network"
        )
        call_tiers.append(
            {
                "call_id": row["call_id"],
                "network_tier": tier,
                "component_loci_at_0.8": sizes80.get(index, 0),
                "component_loci_at_0.5": sizes50.get(index, 0),
                "complete_link_arms_at_0.8": (
                    "PASS" if index in complete80 else "FAIL"
                ),
                "complete_link_arms_at_0.5": (
                    "PASS" if index in complete50 else "FAIL"
                ),
            }
        )

    args.output.mkdir(parents=True, exist_ok=True)

    def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
        with path.open("w", newline="") as handle:
            if not rows:
                return
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    write_tsv(args.output / "strict_two_copy_events.tsv", event_rows)
    write_tsv(args.output / "event_call_membership.tsv", membership_rows)
    write_tsv(args.output / "call_network_tiers.tsv", call_tiers)
    summary = [
        {"metric": "input_biser_calls", "count": len(biser_rows)},
        {
            "metric": "calls_with_complete_link_arms_at_overlap_0.8",
            "count": len(complete80),
        },
        {
            "metric": "calls_with_complete_link_arms_at_overlap_0.5",
            "count": len(complete50),
        },
        {"metric": "strict_calls_at_overlap_0.8", "count": len(strict80)},
        {"metric": "strict_calls_at_overlap_0.5", "count": len(strict50)},
        {"metric": "strict_calls_stable_at_both", "count": len(stable_calls)},
        {
            "metric": "selected_two_copy_events",
            "count": len(event_rows),
        },
        {
            "metric": "strict_stable_two_copy_events",
            "count": len(event_rows),
        },
        {
            "metric": "selected_network_policy",
            "count": args.network_policy,
        },
        {
            "metric": "selected_minimum_length_bp",
            "count": args.minimum_length,
        },
        {
            "metric": "events_passing_strict_primary_node_pd",
            "count": sum(row["strict_pd_status"] == "PASS" for row in event_rows),
        },
    ]
    write_tsv(args.output / "normalization_summary.tsv", summary)
    print(
        f"Normalized {len(biser_rows)} calls into {len(event_rows)} "
        f"two-copy events under network policy {args.network_policy}"
    )


if __name__ == "__main__":
    main()
