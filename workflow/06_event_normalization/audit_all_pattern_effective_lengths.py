#!/usr/bin/env python3
"""Audit uppercase effective lengths and enumerate all 256 node patterns."""

from __future__ import annotations

import argparse
import csv
import itertools
import re
from collections import Counter
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNS])")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


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


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            else:
                sequences[name].append(line.strip())
    return {name: "".join(parts) for name, parts in sequences.items()}


def advance(cursor: int, length: int, strand: str) -> tuple[int, int, int]:
    if strand == "+":
        return cursor, cursor + length, cursor + length
    return cursor - length, cursor, cursor - length


def oriented_sequence(
    genome: dict[str, str],
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> str:
    sequence = genome[chrom][start:end]
    return sequence if strand == "+" else sequence.translate(COMPLEMENT)[::-1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    args = parser.parse_args()

    root = (
        args.project
        / "06_sd_age_tracing_preparation/event_first_reanalysis"
    )
    events = read(root / "events/event_first_events.tsv")
    event_by_id = {row["event_id"]: row for row in events}
    arm_map = {
        (row["call_id"], row["source_copy_label"]): row[
            "canonical_event_locus"
        ]
        for row in read(root / "events/source_arm_to_event_locus.tsv")
    }
    biser_path = (
        args.project
        / "03_biser_segmental_duplication/runs/"
        "annotation_extended_softmask/biser_out"
    )
    biser = {
        str(call_id): line.rstrip().split("\t")
        for call_id, line in enumerate(biser_path.open(), 1)
    }
    genome = read_fasta(
        args.project
        / "01_reference/prepared_data/"
        "TAIR12.Col-CC.annotation_softmasked.fa"
    )

    qc_rows = []
    for event in events:
        event_id = event["event_id"]
        call_id = event["representative_call_id"]
        fields = biser[call_id]
        chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
        chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
        strand1, strand2 = fields[8], fields[9]
        cigar = fields[12]
        operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
        if "".join(f"{length}{op}" for length, op in operations) != cigar:
            raise ValueError(f"Unsupported CIGAR for call {call_id}: {cigar}")

        interval_sequences = {
            "copy1": genome[chrom1][start1:end1],
            "copy2": genome[chrom2][start2:end2],
        }
        uppercase_by_locus = {}
        interval_bp_by_locus = {}
        for copy in ("copy1", "copy2"):
            locus = arm_map[(call_id, copy)]
            sequence = interval_sequences[copy]
            uppercase_by_locus[locus] = sum(base in "ACGT" for base in sequence)
            interval_bp_by_locus[locus] = len(sequence)

        cursor1 = start1 if strand1 == "+" else end1
        cursor2 = start2 if strand2 == "+" else end2
        operation_bp = Counter()
        joint_callable = mismatch_bp = 0
        for length, op in operations:
            operation_bp[op] += length
            interval1 = interval2 = None
            if op in "MIS":
                left, right, cursor1 = advance(cursor1, length, strand1)
                interval1 = (left, right)
            if op in "MDN":
                left, right, cursor2 = advance(cursor2, length, strand2)
                interval2 = (left, right)
            if op != "M":
                continue
            if interval1 is None or interval2 is None:
                raise AssertionError("M must consume both copies")
            seq1 = oriented_sequence(
                genome, chrom1, interval1[0], interval1[1], strand1
            )
            seq2 = oriented_sequence(
                genome, chrom2, interval2[0], interval2[1], strand2
            )
            if len(seq1) != length or len(seq2) != length:
                raise AssertionError(f"Sequence length mismatch: {event_id}")
            for base1, base2 in zip(seq1, seq2):
                if base1 in "ACGT" and base2 in "ACGT":
                    joint_callable += 1
                    mismatch_bp += base1 != base2

        expected1 = end1 if strand1 == "+" else start1
        expected2 = end2 if strand2 == "+" else start2
        if cursor1 != expected1 or cursor2 != expected2:
            raise AssertionError(f"CIGAR coordinate mismatch: {event_id}")

        geometric_pass = event["strict_minimum_1kb_core"] == "PASS"
        both_uppercase_pass = (
            uppercase_by_locus["locus_A"] >= 1000
            and uppercase_by_locus["locus_B"] >= 1000
        )
        paired_m_pass = operation_bp["M"] >= 1000
        joint_callable_pass = joint_callable >= 1000
        qc_rows.append(
            {
                "event_id": event_id,
                "primary_pattern": event["primary_pattern"],
                "representative_call_id": call_id,
                "locus_A_interval_bp": interval_bp_by_locus["locus_A"],
                "locus_A_uppercase_ACGT_bp": uppercase_by_locus["locus_A"],
                "locus_A_uppercase_fraction": (
                    f"{uppercase_by_locus['locus_A'] / interval_bp_by_locus['locus_A']:.6f}"
                ),
                "locus_B_interval_bp": interval_bp_by_locus["locus_B"],
                "locus_B_uppercase_ACGT_bp": uppercase_by_locus["locus_B"],
                "locus_B_uppercase_fraction": (
                    f"{uppercase_by_locus['locus_B'] / interval_bp_by_locus['locus_B']:.6f}"
                ),
                "paired_CIGAR_M_bp": operation_bp["M"],
                "joint_callable_uppercase_ACGT_aligned_bp": joint_callable,
                "joint_callable_fraction_of_M": (
                    f"{joint_callable / operation_bp['M']:.6f}"
                    if operation_bp["M"]
                    else "0"
                ),
                "present_day_PD_mismatch_bp": mismatch_bp,
                "present_day_PD_mismatch_pct_callable": (
                    f"{100 * mismatch_bp / joint_callable:.6f}"
                    if joint_callable
                    else "NA"
                ),
                "former_geometric_BISER_core_ge1kb": (
                    "PASS" if geometric_pass else "FAIL"
                ),
                "both_locus_interval_uppercase_ACGT_ge1kb": (
                    "PASS" if both_uppercase_pass else "FAIL"
                ),
                "paired_CIGAR_M_ge1kb": (
                    "PASS" if paired_m_pass else "FAIL"
                ),
                "joint_callable_uppercase_ACGT_aligned_ge1kb": (
                    "PASS" if joint_callable_pass else "FAIL"
                ),
                "event_first_strict_pd_status": event["strict_pd_status"],
                "event_first_strict_age_bin": event["strict_age_bin"],
                "event_first_threshold_stability": event[
                    "classification_stable_at_overlap_0.25_0.5_0.75"
                ],
            }
        )

    qc = {row["event_id"]: row for row in qc_rows}
    enumeration_rows = []
    for states in itertools.product("0123", repeat=4):
        pattern = "".join(states)
        subset = [
            event for event in events if event["primary_pattern"] == pattern
        ]
        enumeration_rows.append(
            {
                "primary_pattern": pattern,
                "Alyrata_state": states[0],
                "Bstricta_state": states[1],
                "Dstrictus_state": states[2],
                "Cviolacea_state": states[3],
                "normalized_events": len(subset),
                "former_geometric_BISER_core_ge1kb": sum(
                    qc[event["event_id"]][
                        "former_geometric_BISER_core_ge1kb"
                    ]
                    == "PASS"
                    for event in subset
                ),
                "both_locus_interval_uppercase_ACGT_ge1kb": sum(
                    qc[event["event_id"]][
                        "both_locus_interval_uppercase_ACGT_ge1kb"
                    ]
                    == "PASS"
                    for event in subset
                ),
                "paired_CIGAR_M_ge1kb": sum(
                    qc[event["event_id"]]["paired_CIGAR_M_ge1kb"] == "PASS"
                    for event in subset
                ),
                "joint_callable_uppercase_ACGT_aligned_ge1kb": sum(
                    qc[event["event_id"]][
                        "joint_callable_uppercase_ACGT_aligned_ge1kb"
                    ]
                    == "PASS"
                    for event in subset
                ),
            }
        )

    summary = [
        {
            "criterion": "all_normalized_events",
            "events": len(qc_rows),
            "definition": "network-stable two-locus events",
        },
        {
            "criterion": "former_geometric_BISER_core_ge1kb",
            "events": sum(
                row["former_geometric_BISER_core_ge1kb"] == "PASS"
                for row in qc_rows
            ),
            "definition": (
                "both representative intervals, max-mate length, and "
                "alignment span >=1000 bp"
            ),
        },
        {
            "criterion": "both_locus_interval_uppercase_ACGT_ge1kb",
            "events": sum(
                row["both_locus_interval_uppercase_ACGT_ge1kb"] == "PASS"
                for row in qc_rows
            ),
            "definition": (
                "each representative physical locus contains >=1000 "
                "uppercase A/C/G/T bases"
            ),
        },
        {
            "criterion": "paired_CIGAR_M_ge1kb",
            "events": sum(
                row["paired_CIGAR_M_ge1kb"] == "PASS" for row in qc_rows
            ),
            "definition": "BISER paired M operations total >=1000 bp",
        },
        {
            "criterion": "joint_callable_uppercase_ACGT_aligned_ge1kb",
            "events": sum(
                row["joint_callable_uppercase_ACGT_aligned_ge1kb"] == "PASS"
                for row in qc_rows
            ),
            "definition": (
                ">=1000 aligned M columns where both copies are uppercase "
                "A/C/G/T; recommended effective unmasked criterion"
            ),
        },
    ]
    age_rows = []
    for age in ("time1", "time2", "time3", "time4"):
        selected = [
            event
            for event in events
            if event["strict_pd_status"] == "PASS"
            and event["strict_age_bin"] == age
        ]
        age_rows.append(
            {
                "age_bin": age,
                "event_first_geometric_ge1kb_events": len(selected),
                "joint_callable_uppercase_aligned_ge1kb_events": sum(
                    qc[event["event_id"]][
                        "joint_callable_uppercase_ACGT_aligned_ge1kb"
                    ]
                    == "PASS"
                    for event in selected
                ),
                "effective_ge1kb_and_node_threshold_stable_events": sum(
                    qc[event["event_id"]][
                        "joint_callable_uppercase_ACGT_aligned_ge1kb"
                    ]
                    == "PASS"
                    and event[
                        "classification_stable_at_overlap_0.25_0.5_0.75"
                    ]
                    == "PASS"
                    for event in selected
                ),
            }
        )

    write(root / "events/event_effective_unmasked_length_qc.tsv", qc_rows)
    write(
        root / "statistics/all_256_pattern_counts.tsv",
        enumeration_rows,
    )
    write(
        root / "statistics/observed_pattern_counts.tsv",
        [row for row in enumeration_rows if row["normalized_events"]],
    )
    write(
        root / "statistics/effective_length_criterion_summary.tsv",
        summary,
    )
    write(
        root / "statistics/effective_length_age_summary.tsv",
        age_rows,
    )
    print(
        f"Audited {len(qc_rows)} events; "
        f"{summary[-1]['events']} have >=1 kb jointly callable uppercase "
        "aligned sequence"
    )


if __name__ == "__main__":
    main()
