#!/usr/bin/env python3
"""Compare TAIR12 source, GFF-extended, and de-novo RepeatMasker masks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import statistics
from collections import Counter, defaultdict
from pathlib import Path


ACCESSION_TO_CHR = {
    "OZ408683.1": "Chr1",
    "OZ408684.1": "Chr2",
    "OZ408685.1": "Chr3",
    "OZ408686.1": "Chr4",
    "OZ408687.1": "Chr5",
}


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                raw = line[1:].split()[0]
                name = ACCESSION_TO_CHR.get(raw, raw)
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {name: "".join(parts) for name, parts in seqs.items()}


def runs_from_predicate(seq: str, predicate) -> list[tuple[int, int]]:
    runs = []
    start = None
    for i, base in enumerate(seq):
        hit = predicate(base)
        if hit and start is None:
            start = i
        elif not hit and start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, len(seq)))
    return runs


def overlap_bp(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    i = j = total = 0
    while i < len(a) and j < len(b):
        start = max(a[i][0], b[j][0])
        end = min(a[i][1], b[j][1])
        if end > start:
            total += end - start
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return total


def interval_overlap(intervals: list[tuple[int, int]], start: int, end: int) -> int:
    """Return overlap of one half-open interval with sorted disjoint intervals."""
    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        if intervals[mid][1] <= start:
            lo = mid + 1
        else:
            hi = mid
    i = lo
    total = 0
    while i < len(intervals) and intervals[i][0] < end:
        total += max(0, min(end, intervals[i][1]) - max(start, intervals[i][0]))
        i += 1
    return total


def interval_stats(intervals: list[tuple[int, int]]) -> tuple[int, int, float, int]:
    lengths = [end - start for start, end in intervals]
    return len(lengths), sum(lengths), statistics.median(lengths), max(lengths)


def parse_repeatmasker_out(path: Path):
    by_class: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            chrom = ACCESSION_TO_CHR.get(fields[4], fields[4])
            start, end = int(fields[5]) - 1, int(fields[6])
            repeat_class = fields[10].split("/")[0]
            by_class[repeat_class][chrom].append((start, end))
    return by_class


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return [(start, end) for start, end in merged]


def write_tsv(path: Path, header: list[str], rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    ref = args.reference_dir
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    paths = {
        "source": ref / "prepared_data/TAIR12.Col-CC.source_softmasked.fa",
        "gff_extended": ref / "prepared_data/TAIR12.Col-CC.annotation_softmasked.fa",
        "repeatmasker": ref
        / "repeatmasker_arabidopsis/GCA_978657495.1_TAIR12_genomic.uppercase.fna.masked",
    }
    seqs = {label: read_fasta(path) for label, path in paths.items()}
    chroms = [f"Chr{i}" for i in range(1, 6)]
    if any(set(genome) != set(chroms) for genome in seqs.values()):
        raise ValueError("Unexpected chromosome names")

    masks: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(dict)
    case_rows = []
    identity_rows = []
    for chrom in chroms:
        canonical = seqs["source"][chrom].upper()
        for label in paths:
            seq = seqs[label][chrom]
            if len(seq) != len(canonical) or seq.upper() != canonical:
                raise ValueError(f"Case-insensitive sequence mismatch: {label} {chrom}")
            intervals = runs_from_predicate(seq, str.islower)
            masks[label][chrom] = intervals
            n_runs, bp, median_len, max_len = interval_stats(intervals)
            case_rows.append(
                [label, chrom, len(seq), bp, f"{100 * bp / len(seq):.6f}",
                 n_runs, median_len, max_len]
            )
        identity_rows.append(
            [chrom, len(canonical), hashlib.sha256(canonical.encode()).hexdigest()]
        )

    write_tsv(
        out / "mask_by_chromosome.tsv",
        ["mask", "chromosome", "length_bp", "masked_bp", "masked_pct",
         "mask_runs", "median_run_bp", "max_run_bp"],
        case_rows,
    )
    write_tsv(out / "sequence_identity.tsv", ["chromosome", "length_bp", "uppercase_sha256"], identity_rows)

    pattern_rows = []
    pattern_totals = Counter()
    labels = ["source", "gff_extended", "repeatmasker"]
    for chrom in chroms:
        strings = [seqs[label][chrom] for label in labels]
        counts = Counter(
            (int(a.islower()) << 2) | (int(b.islower()) << 1) | int(c.islower())
            for a, b, c in zip(*strings)
        )
        for code in range(8):
            bp = counts[code]
            pattern_totals[code] += bp
            pattern_rows.append(
                [chrom, (code >> 2) & 1, (code >> 1) & 1, code & 1, bp]
            )
    for code in range(8):
        pattern_rows.append(
            ["TOTAL", (code >> 2) & 1, (code >> 1) & 1, code & 1,
             pattern_totals[code]]
        )
    write_tsv(
        out / "mask_venn_patterns.tsv",
        ["chromosome", "source_masked", "gff_extended_masked",
         "repeatmasker_masked", "bp"],
        pattern_rows,
    )

    pair_rows = []
    for left, right in [("source", "gff_extended"), ("source", "repeatmasker"),
                        ("gff_extended", "repeatmasker")]:
        left_bp = right_bp = shared = 0
        for chrom in chroms:
            l_bp = sum(e - s for s, e in masks[left][chrom])
            r_bp = sum(e - s for s, e in masks[right][chrom])
            ov = overlap_bp(masks[left][chrom], masks[right][chrom])
            left_bp += l_bp
            right_bp += r_bp
            shared += ov
        union = left_bp + right_bp - shared
        pair_rows.append(
            [left, right, left_bp, right_bp, shared, left_bp - shared,
             right_bp - shared, union, f"{shared / union:.6f}",
             f"{shared / left_bp:.6f}", f"{shared / right_bp:.6f}"]
        )
    write_tsv(
        out / "pairwise_mask_overlap.tsv",
        ["mask_a", "mask_b", "a_bp", "b_bp", "intersection_bp", "a_only_bp",
         "b_only_bp", "union_bp", "jaccard", "fraction_a_recovered",
         "fraction_b_recovered"],
        pair_rows,
    )

    annotated: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with (ref / "prepared_data/TAIR12.annotated_repeats.merged.bed").open() as handle:
        for line in handle:
            chrom, start, end = line.split()[:3]
            annotated[chrom].append((int(start), int(end)))
    annotated_bp = sum(e - s for chrom in chroms for s, e in annotated[chrom])
    rm_annotated_overlap = sum(
        overlap_bp(masks["repeatmasker"][chrom], annotated[chrom]) for chrom in chroms
    )
    rm_bp = sum(
        e - s for chrom in chroms for s, e in masks["repeatmasker"][chrom]
    )
    union_bp = rm_bp + annotated_bp - rm_annotated_overlap
    write_tsv(
        out / "candidate_mask_combinations.tsv",
        ["mask_definition", "masked_bp", "masked_pct"],
        [
            ["repeatmasker_only", rm_bp, f"{100 * rm_bp / 142481245:.6f}"],
            ["gff_selected_repeats_only", annotated_bp,
             f"{100 * annotated_bp / 142481245:.6f}"],
            ["repeatmasker_union_gff_selected_repeats", union_bp,
             f"{100 * union_bp / 142481245:.6f}"],
        ],
    )

    rm_classes = parse_repeatmasker_out(
        ref / "repeatmasker_arabidopsis/GCA_978657495.1_TAIR12_genomic.uppercase.fna.out"
    )
    class_rows = []
    for repeat_class, by_chrom in rm_classes.items():
        raw_count = sum(len(v) for v in by_chrom.values())
        merged_by_chrom = {chrom: merge(by_chrom.get(chrom, [])) for chrom in chroms}
        union_bp = sum(sum(e - s for s, e in v) for v in merged_by_chrom.values())
        source_overlap = sum(
            overlap_bp(merged_by_chrom[chrom], masks["source"][chrom]) for chrom in chroms
        )
        gff_overlap = sum(
            overlap_bp(merged_by_chrom[chrom], masks["gff_extended"][chrom]) for chrom in chroms
        )
        class_rows.append(
            [repeat_class, raw_count, union_bp, source_overlap,
             f"{source_overlap / union_bp:.6f}", gff_overlap,
             f"{gff_overlap / union_bp:.6f}"]
        )
    class_rows.sort(key=lambda row: row[2], reverse=True)
    write_tsv(
        out / "repeatmasker_class_overlap.tsv",
        ["repeat_class", "raw_hits", "class_union_bp", "overlap_source_bp",
         "fraction_in_source", "overlap_gff_extended_bp", "fraction_in_gff_extended"],
        class_rows,
    )

    contrast_masks: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(dict)
    for chrom in chroms:
        gff_seq = seqs["gff_extended"][chrom]
        rm_seq = seqs["repeatmasker"][chrom]
        contrast_masks["rm_newly_masked"][chrom] = runs_from_predicate(
            "".join(
                "x" if rm.islower() and not gff.islower() else "X"
                for gff, rm in zip(gff_seq, rm_seq)
            ),
            str.islower,
        )
        contrast_masks["rm_reopened"][chrom] = runs_from_predicate(
            "".join(
                "x" if gff.islower() and not rm.islower() else "X"
                for gff, rm in zip(gff_seq, rm_seq)
            ),
            str.islower,
        )

    biser_root = ref.parent / "03_biser_segmental_duplication/runs"
    call_rows = []
    for condition, call_path in [
        ("source_softmask", biser_root / "source_softmask/biser_out"),
        ("gff_extended_softmask", biser_root / "annotation_extended_softmask/biser_out"),
    ]:
        calls = touch_new = touch_new_both = ge10 = ge50 = 0
        current_both_ge1k = rm_both_ge1k = current_ge1k_rm_lt1k = 0
        arm_span_bp = rm_new_bp = rm_reopened_bp = 0
        with call_path.open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                arms = [
                    (fields[0], int(fields[1]), int(fields[2])),
                    (fields[3], int(fields[4]), int(fields[5])),
                ]
                calls += 1
                new_per_arm = []
                current_visible = []
                rm_visible = []
                call_span = 0
                for chrom, start, end in arms:
                    span = end - start
                    call_span += span
                    new_bp = interval_overlap(
                        contrast_masks["rm_newly_masked"][chrom], start, end
                    )
                    reopened_bp = interval_overlap(
                        contrast_masks["rm_reopened"][chrom], start, end
                    )
                    current_masked = interval_overlap(
                        masks["gff_extended"][chrom], start, end
                    )
                    rm_masked = interval_overlap(
                        masks["repeatmasker"][chrom], start, end
                    )
                    new_per_arm.append(new_bp)
                    current_visible.append(span - current_masked)
                    rm_visible.append(span - rm_masked)
                    rm_new_bp += new_bp
                    rm_reopened_bp += reopened_bp
                arm_span_bp += call_span
                new_total = sum(new_per_arm)
                touch_new += new_total > 0
                touch_new_both += all(value > 0 for value in new_per_arm)
                ge10 += new_total / call_span >= 0.10
                ge50 += new_total / call_span >= 0.50
                current_ok = all(value >= 1000 for value in current_visible)
                rm_ok = all(value >= 1000 for value in rm_visible)
                current_both_ge1k += current_ok
                rm_both_ge1k += rm_ok
                current_ge1k_rm_lt1k += current_ok and not rm_ok
        call_rows.append(
            [condition, calls, arm_span_bp, touch_new, f"{100*touch_new/calls:.4f}",
             touch_new_both, ge10, ge50, rm_new_bp, rm_reopened_bp,
             current_both_ge1k, rm_both_ge1k, current_ge1k_rm_lt1k]
        )
    write_tsv(
        out / "existing_biser_call_mask_sensitivity.tsv",
        ["condition", "calls", "total_two_arm_span_bp", "calls_touching_rm_new_mask",
         "pct_calls_touching_rm_new_mask", "calls_both_arms_touching_rm_new_mask",
         "calls_ge10pct_span_newly_masked", "calls_ge50pct_span_newly_masked",
         "rm_newly_masked_bp_in_arms", "rm_reopened_bp_in_arms",
         "calls_both_arms_current_visible_ge1kb", "calls_both_arms_rm_visible_ge1kb",
         "current_ge1kb_but_rm_lt1kb_calls"],
        call_rows,
    )


if __name__ == "__main__":
    main()
