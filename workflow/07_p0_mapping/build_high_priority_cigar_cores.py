#!/usr/bin/env python3
"""Select P/D-polarized events and parse paired BISER CIGAR cores.

The default retains the historical time1-time3 selection.  ``--age-free``
removes the age-bin admission rule; the supplied event table must already
encode the accepted P/D orientation and its stability.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNS])")
SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        if not rows:
            return
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
    genome: dict[str, str], chrom: str, start: int, end: int, strand: str
) -> str:
    sequence = genome[chrom][start:end]
    return sequence if strand == "+" else sequence.translate(COMPLEMENT)[::-1]


def fasta_record(name: str, sequence: str, width: int = 80) -> str:
    return f">{name}\n" + "\n".join(
        sequence[index : index + width]
        for index in range(0, len(sequence), width)
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--events", type=Path)
    parser.add_argument("--sensitivity", type=Path)
    parser.add_argument(
        "--age-free",
        action="store_true",
        help="Admit stable P/D events regardless of strict_age_bin.",
    )
    parser.add_argument(
        "--minimum-core-bp",
        type=int,
        default=1,
        help="Technical non-empty minimum; no 1-kb biological gate is applied.",
    )
    parser.add_argument(
        "--minimum-callable-bp",
        type=int,
        default=1,
        help="Technical non-empty minimum; local comparable sequence is retained.",
    )
    args = parser.parse_args()

    event_first = (
        args.project
        / "06_sd_age_tracing_preparation/event_first_reanalysis"
    )
    legacy_analysis = (
        args.project / "06_sd_age_tracing_preparation/primary_node_analysis"
    )
    event_path = args.events or event_first / "events/event_first_events.tsv"
    events = read_tsv(event_path)
    is_event_first = bool(events and "locus_A_representative_start" in events[0])
    sensitivity = {}
    if not is_event_first:
        sensitivity_path = (
            args.sensitivity
            or legacy_analysis
            / "statistics/event_node_overlap_threshold_sensitivity.tsv"
        )
        sensitivity = {
            row["event_id"]: row for row in read_tsv(sensitivity_path)
        }
    arm_map = {}
    effective_qc = {}
    if is_event_first:
        arm_map = {
            (row["event_id"], row["call_id"], row["canonical_event_locus"]):
            row["source_copy_label"]
            for row in read_tsv(
                event_first / "events/source_arm_to_event_locus.tsv"
            )
        }
        effective_qc = {
            row["event_id"]: row
            for row in read_tsv(
                event_first / "events/event_effective_unmasked_length_qc.tsv"
            )
        }
    biser_path = (
        args.project
        / "03_biser_segmental_duplication/runs/annotation_extended_softmask/biser_out"
    )
    biser = {
        str(call_id): line.rstrip().split("\t")
        for call_id, line in enumerate(biser_path.open(), 1)
    }
    genome = read_fasta(
        args.project
        / "01_reference/prepared_data/TAIR12.Col-CC.annotation_softmasked.fa"
    )

    selected = []
    for event in events:
        stable = (
            event["classification_stable_at_overlap_0.25_0.5_0.75"]
            if is_event_first
            else sensitivity[event["event_id"]][
                "classification_stable_at_0.25_0.5_0.75"
            ]
        )
        if event["strict_pd_status"] != "PASS" or stable != "PASS":
            continue
        if not args.age_free and event["strict_age_bin"] not in {
            "time1",
            "time2",
            "time3",
        }:
            continue
        selected.append(event)

    event_rows: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    fasta_parts: list[str] = []

    for event in selected:
        call_id = event["representative_call_id"]
        fields = biser[call_id]
        chrom1, start1, end1 = fields[0], int(fields[1]), int(fields[2])
        chrom2, start2, end2 = fields[3], int(fields[4]), int(fields[5])
        strand1, strand2 = fields[8], fields[9]
        cigar = fields[12]
        operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
        if "".join(f"{length}{op}" for length, op in operations) != cigar:
            raise ValueError(f"Unsupported CIGAR for call {call_id}: {cigar}")

        p_copy = (
            arm_map[
                (
                    event["event_id"],
                    call_id,
                    event["provisional_p_locus"],
                )
            ]
            if is_event_first
            else event["provisional_p_copy"]
        )
        event_row: dict[str, object] = dict(event)
        event_row["provisional_p_copy"] = p_copy
        event_row["inclusive_sequence_analysis_status"] = "INCLUDE"
        event_row["legacy_length_filters_policy"] = "DIAGNOSTIC_NOT_GATING"
        if is_event_first:
            event_row.update(
                {
                    "representative_copy1_chrom": chrom1,
                    "representative_copy1_start": start1,
                    "representative_copy1_end": end1,
                    "representative_copy2_chrom": chrom2,
                    "representative_copy2_start": start2,
                    "representative_copy2_end": end2,
                }
            )
            for field in (
                "former_geometric_BISER_core_ge1kb",
                "both_locus_interval_uppercase_ACGT_ge1kb",
                "paired_CIGAR_M_ge1kb",
                "joint_callable_uppercase_ACGT_aligned_ge1kb",
            ):
                event_row[f"legacy_{field}_diagnostic"] = effective_qc[
                    event["event_id"]
                ][field]
            for code in SPECIES:
                for copy in (1, 2):
                    locus = next(
                        locus
                        for locus in ("locus_A", "locus_B")
                        if arm_map[(event["event_id"], call_id, locus)]
                        == f"copy{copy}"
                    )
                    event_row[f"{code}_copy{copy}_block_ids"] = event[
                        f"{code}_{locus}_block_ids"
                    ]
        event_rows.append(event_row)

        cursor1 = start1 if strand1 == "+" else end1
        cursor2 = start2 if strand2 == "+" else end2
        operation_bp = Counter()
        copy_sequences: dict[str, list[str]] = {"copy1": [], "copy2": []}
        callable_bp = mismatch_bp = 0
        core_index = 0

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
                raise AssertionError("M operation must consume both mates")
            core_index += 1
            seq1 = oriented_sequence(
                genome, chrom1, interval1[0], interval1[1], strand1
            )
            seq2 = oriented_sequence(
                genome, chrom2, interval2[0], interval2[1], strand2
            )
            if len(seq1) != length or len(seq2) != length:
                raise AssertionError(f"Sequence length mismatch for call {call_id}")
            copy_sequences["copy1"].append(seq1)
            copy_sequences["copy2"].append(seq2)
            for base1, base2 in zip(seq1, seq2):
                if base1 in "ACGT" and base2 in "ACGT":
                    callable_bp += 1
                    mismatch_bp += base1 != base2
            role1 = "P" if p_copy == "copy1" else "D"
            role2 = "P" if p_copy == "copy2" else "D"
            block_rows.append(
                {
                    "event_id": event["event_id"],
                    "call_id": call_id,
                    "core_block_index": core_index,
                    "core_block_bp": length,
                    "copy1_role": role1,
                    "copy1_chrom": chrom1,
                    "copy1_start": interval1[0],
                    "copy1_end": interval1[1],
                    "copy1_strand": strand1,
                    "copy2_role": role2,
                    "copy2_chrom": chrom2,
                    "copy2_start": interval2[0],
                    "copy2_end": interval2[1],
                    "copy2_strand": strand2,
                }
            )

        expected1 = end1 if strand1 == "+" else start1
        expected2 = end2 if strand2 == "+" else start2
        if cursor1 != expected1 or cursor2 != expected2:
            raise AssertionError(
                f"CIGAR coordinate mismatch for call {call_id}: "
                f"{cursor1}/{expected1}, {cursor2}/{expected2}"
            )

        core_bp = operation_bp["M"]
        eligible = (
            core_bp >= args.minimum_core_bp
            and callable_bp >= args.minimum_callable_bp
        )
        role_sequences = {
            ("P" if p_copy == copy else "D"): "".join(parts)
            for copy, parts in copy_sequences.items()
        }
        if eligible:
            fasta_parts.append(
                fasta_record(
                    f"{event['event_id']}|P|{event['strict_age_bin']}",
                    role_sequences["P"],
                )
            )
            fasta_parts.append(
                fasta_record(
                    f"{event['event_id']}|D|{event['strict_age_bin']}",
                    role_sequences["D"],
                )
            )
        qc_rows.append(
            {
                "event_id": event["event_id"],
                "age_bin": event["strict_age_bin"],
                "representative_call_id": call_id,
                "primary_pattern": event["primary_pattern"],
                "p_copy": p_copy,
                "cigar_core_blocks": core_index,
                "paired_M_core_bp": core_bp,
                "jointly_callable_uppercase_acgt_bp": callable_bp,
                "jointly_callable_fraction_of_M": (
                    f"{callable_bp / core_bp:.6f}" if core_bp else "0"
                ),
                "present_day_PD_mismatch_bp": mismatch_bp,
                "present_day_PD_mismatch_pct_callable": (
                    f"{100 * mismatch_bp / callable_bp:.6f}"
                    if callable_bp
                    else "NA"
                ),
                "mate1_insertion_I_bp": operation_bp["I"],
                "mate2_insertion_D_bp": operation_bp["D"],
                "mate1_softmasked_S_bp": operation_bp["S"],
                "mate2_softmasked_N_bp": operation_bp["N"],
                "minimum_core_bp_rule": args.minimum_core_bp,
                "minimum_callable_bp_rule": args.minimum_callable_bp,
                "legacy_1kb_filters_policy": "DIAGNOSTIC_NOT_GATING",
                "core_minimum_rule": (
                    "PASS" if core_bp >= args.minimum_core_bp else "FAIL"
                ),
                "joint_callable_minimum_rule": (
                    "PASS"
                    if callable_bp >= args.minimum_callable_bp
                    else "FAIL"
                ),
                "core_minimum_1kb": "PASS" if core_bp >= 1000 else "FAIL",
                "joint_callable_minimum_1kb": (
                    "PASS" if callable_bp >= 1000 else "FAIL"
                ),
                "pilot_core_eligible": "PASS" if eligible else "FAIL",
            }
        )

    args.output.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output / "inputs/high_priority_events.tsv", event_rows)
    write_tsv(args.output / "core/homologous_core_blocks.tsv", block_rows)
    write_tsv(args.output / "core/core_event_qc.tsv", qc_rows)
    eligible_ids = {
        row["event_id"] for row in qc_rows if row["pilot_core_eligible"] == "PASS"
    }
    write_tsv(
        args.output / "inputs/high_priority_core_eligible.tsv",
        [row for row in event_rows if row["event_id"] in eligible_ids],
    )
    fasta_path = args.output / "core/TAIR12_PD_homologous_cores.fa"
    fasta_path.write_text("".join(fasta_parts))

    summary = Counter(row["age_bin"] for row in qc_rows)
    eligible_summary = Counter(
        row["age_bin"] for row in qc_rows if row["pilot_core_eligible"] == "PASS"
    )
    age_labels = sorted(summary) if args.age_free else ["time1", "time2", "time3"]
    write_tsv(
        args.output / "core/core_selection_summary.tsv",
        [
            {
                "age_bin": age,
                "threshold_stable_events": summary[age],
                "core_eligible_events": eligible_summary[age],
                "excluded_after_cigar_core_qc": summary[age]
                - eligible_summary[age],
            }
            for age in age_labels
        ],
    )
    print(
        f"Selected {len(selected)} stable "
        f"{'age-free P/D' if args.age_free else 'time1-time3'} events; "
        f"{len(eligible_ids)} contain non-empty paired and jointly callable "
        "sequence; legacy 1-kb filters are diagnostic only"
    )


if __name__ == "__main__":
    main()
