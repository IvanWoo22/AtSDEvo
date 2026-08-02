#!/usr/bin/env python3
"""Score annotation-extended TAIR12 BISER calls against four primary nodes."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
PAIR_RE = re.compile(r"^\s*\d+-\s*\d+:\s+(\S+)\s+(\S+)\s+\S+")
BLOCK_RE = re.compile(r"^## Alignment (\d+):")


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


def read_atha_genes(path: Path) -> dict[str, tuple[str, int, int]]:
    genes = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            pos1, pos2 = int(fields[2]), int(fields[3])
            genes[fields[1]] = (
                fields[0].removeprefix("Atha_"),
                min(pos1, pos2) - 1,
                max(pos1, pos2),
            )
    return genes


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def project_blocks(
    path: Path,
    atha_genes: dict[str, tuple[str, int, int]],
    flank: int,
) -> tuple[
    dict[str, list[tuple[int, int]]],
    dict[str, list[tuple[int, int, int]]],
    list[dict[str, object]],
]:
    raw: dict[str, list[tuple[int, int]]] = {}
    indexed: dict[str, list[tuple[int, int, int]]] = {}
    block_rows: list[dict[str, object]] = []
    current_id = None
    current_genes: list[str] = []

    def finish() -> None:
        if current_id is None or not current_genes:
            return
        by_chrom: dict[str, list[tuple[int, int]]] = {}
        for gene in dict.fromkeys(current_genes):
            chrom, start, end = atha_genes[gene]
            by_chrom.setdefault(chrom, []).append((start, end))
        for chrom, spans in by_chrom.items():
            start = max(0, min(span[0] for span in spans) - flank)
            end = max(span[1] for span in spans) + flank
            raw.setdefault(chrom, []).append((start, end))
            indexed.setdefault(chrom, []).append((start, end, int(current_id)))
            block_rows.append(
                {
                    "block_id": current_id,
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "atha_anchor_genes": len(spans),
                    "flank_bp": flank,
                }
            )

    with path.open() as handle:
        for line in handle:
            block = BLOCK_RE.match(line)
            if block:
                finish()
                current_id = int(block.group(1))
                current_genes = []
                continue
            pair = PAIR_RE.match(line)
            if pair:
                gene1, gene2 = pair.groups()
                if gene1 in atha_genes:
                    current_genes.append(gene1)
                elif gene2 in atha_genes:
                    current_genes.append(gene2)
        finish()

    merged = {chrom: merge_intervals(values) for chrom, values in raw.items()}
    for chrom in indexed:
        indexed[chrom].sort()
    return merged, indexed, block_rows


def interval_overlap(
    start: int, end: int, intervals: list[tuple[int, int]]
) -> int:
    total = 0
    for other_start, other_end in intervals:
        if other_start >= end:
            break
        if other_end > start:
            total += max(0, min(end, other_end) - max(start, other_start))
    return total


def overlapping_blocks(
    start: int, end: int, blocks: list[tuple[int, int, int]]
) -> list[int]:
    return sorted(
        {
            block_id
            for block_start, block_end, block_id in blocks
            if block_start < end and block_end > start
        }
    )


def read_biser(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open() as handle:
        for call_id, line in enumerate(handle, 1):
            fields = line.rstrip().split("\t")
            start1, end1 = sorted((int(fields[1]), int(fields[2])))
            start2, end2 = sorted((int(fields[4]), int(fields[5])))
            rows.append(
                {
                    "call_id": call_id,
                    "copy1_chrom": fields[0],
                    "copy1_start": start1,
                    "copy1_end": end1,
                    "copy2_chrom": fields[3],
                    "copy2_start": start2,
                    "copy2_end": end2,
                    "relative_orientation": (
                        "same" if fields[8] == fields[9] else "opposite"
                    ),
                    "biser_error": fields[7],
                    "biser_max_mate_length_bp": fields[10],
                    "biser_alignment_span_bp": fields[11],
                    "biser_cigar": fields[12],
                }
            )
    return rows


def strict_classify(states: tuple[int, ...]) -> tuple[str, str, str]:
    if not any(states):
        return "EXCLUDE", "all_zero", "NA"

    prefix_both = 0
    while prefix_both < len(states) and states[prefix_both] == 3:
        prefix_both += 1
    if any(state == 3 for state in states[prefix_both:]):
        return "EXCLUDE", "nonmonotonic_both_reappears", "NA"
    if prefix_both == len(states):
        return "EXCLUDE", "older_than_N4_unpolarized", "NA"

    boundary = states[prefix_both]
    if boundary == 0:
        return "EXCLUDE", "boundary_node_uninformative", "NA"
    if boundary not in (1, 2):
        return "EXCLUDE", "invalid_boundary_state", "NA"
    if any(state not in (0, boundary) for state in states[prefix_both + 1 :]):
        return "EXCLUDE", "older_nodes_copy_conflict", "NA"

    return "PASS", f"time{prefix_both + 1}", f"copy{boundary}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--biser",
        type=Path,
        help="覆盖默认的注释扩展 soft-mask BISER 结果路径。",
    )
    parser.add_argument("--flank", type=int, default=2000)
    parser.add_argument("--thresholds", default="0.25,0.5,0.75")
    args = parser.parse_args()

    mcscan = args.project / "05_mcscanx_synteny"
    biser_path = args.biser or (
        args.project
        / "03_biser_segmental_duplication/runs/annotation_extended_softmask/biser_out"
    )
    atha_genes = read_atha_genes(mcscan / "inputs/Atha.gff")
    calls = read_biser(biser_path)
    thresholds = [float(value) for value in args.thresholds.split(",")]

    projected = {}
    block_index = {}
    projection_summary = []
    for code in SPECIES:
        collinearity = (
            mcscan / f"results/Atha_{code}/Atha_{code}.collinearity"
        )
        merged, indexed, blocks = project_blocks(
            collinearity, atha_genes, args.flank
        )
        projected[code] = merged
        block_index[code] = indexed
        for block in blocks:
            block["species"] = code
        write_tsv(args.output / f"blocks/{code}.projected_blocks.tsv", blocks)
        bed = args.output / f"projected_intervals/{code}.TAIR12_collinear.bed"
        bed.parent.mkdir(parents=True, exist_ok=True)
        with bed.open("w") as handle:
            for chrom in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
                for start, end in merged.get(chrom, []):
                    handle.write(f"{chrom}\t{start}\t{end}\n")
                projection_summary.append(
                    {
                        "species": code,
                        "chromosome": chrom,
                        "merged_intervals": len(merged.get(chrom, [])),
                        "projected_bp": sum(
                            end - start for start, end in merged.get(chrom, [])
                        ),
                        "source_blocks": sum(
                            block["chromosome"] == chrom for block in blocks
                        ),
                        "flank_bp": args.flank,
                    }
                )
    write_tsv(args.output / "statistics/projected_interval_summary.tsv", projection_summary)

    threshold_summary = []
    threshold_classifications = []
    main_rows: list[dict[str, object]] = []
    for threshold in thresholds:
        class_counts = Counter()
        pattern_counts = Counter()
        species_counts = {code: Counter() for code in SPECIES}
        scored_rows = []
        for call in calls:
            row = dict(call)
            states = []
            for code in SPECIES:
                pass_copy = []
                for copy in (1, 2):
                    chrom = str(call[f"copy{copy}_chrom"])
                    start = int(call[f"copy{copy}_start"])
                    end = int(call[f"copy{copy}_end"])
                    overlap_bp = interval_overlap(
                        start, end, projected[code].get(chrom, [])
                    )
                    fraction = overlap_bp / max(1, end - start)
                    blocks = overlapping_blocks(
                        start, end, block_index[code].get(chrom, [])
                    )
                    row[f"{code}_copy{copy}_overlap_bp"] = overlap_bp
                    row[f"{code}_copy{copy}_overlap_fraction"] = f"{fraction:.6f}"
                    row[f"{code}_copy{copy}_block_ids"] = ",".join(map(str, blocks))
                    pass_copy.append(fraction >= threshold)
                state = (
                    3 if all(pass_copy)
                    else 1 if pass_copy[0]
                    else 2 if pass_copy[1]
                    else 0
                )
                row[f"{code}_state"] = state
                states.append(state)
                species_counts[code][state] += 1
            status, age_bin, p_copy = strict_classify(tuple(states))
            row["primary_pattern"] = "".join(map(str, states))
            row["strict_status"] = status
            row["strict_age_bin"] = age_bin
            row["provisional_p_copy"] = p_copy
            threshold_classifications.append(
                {
                    "call_id": call["call_id"],
                    "threshold": threshold,
                    "primary_pattern": row["primary_pattern"],
                    "strict_status": status,
                    "strict_age_bin": age_bin,
                    "provisional_p_copy": p_copy,
                }
            )
            class_counts[status] += 1
            pattern_counts[(status, age_bin)] += 1
            scored_rows.append(row)
        for code in SPECIES:
            threshold_summary.append(
                {
                    "threshold": threshold,
                    "summary_type": "species_state",
                    "category": code,
                    "count0": species_counts[code][0],
                    "count1": species_counts[code][1],
                    "count2": species_counts[code][2],
                    "count3": species_counts[code][3],
                    "total": len(calls),
                }
            )
        for (status, age_bin), count in sorted(pattern_counts.items()):
            threshold_summary.append(
                {
                    "threshold": threshold,
                    "summary_type": status,
                    "category": age_bin,
                    "count0": "NA",
                    "count1": "NA",
                    "count2": "NA",
                    "count3": "NA",
                    "total": count,
                }
            )
        if threshold == 0.5:
            main_rows = scored_rows

    write_tsv(args.output / "call_evidence.threshold_0.5.tsv", main_rows)
    write_tsv(
        args.output / "statistics/call_classification_by_threshold.tsv",
        threshold_classifications,
    )
    write_tsv(args.output / "statistics/threshold_and_strict_summary.tsv", threshold_summary)
    print(
        f"Scored {len(calls)} annotation-extended BISER calls against "
        f"{len(SPECIES)} primary nodes"
    )


if __name__ == "__main__":
    main()
