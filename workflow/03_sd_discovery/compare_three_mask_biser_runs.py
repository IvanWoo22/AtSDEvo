#!/usr/bin/env python3
"""Compare BISER calls from source, GFF-extended, and RepeatMasker masks."""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from scipy.stats import binomtest, wilcoxon


CONDITIONS = {
    "source_softmask": "source_softmask",
    "gff_extended_softmask": "annotation_extended_softmask",
    "repeatmasker_arabidopsis_softmask": "repeatmasker_arabidopsis_softmask",
}
ACCESSION_TO_CHR = {
    "OZ408683.1": "Chr1", "OZ408684.1": "Chr2", "OZ408685.1": "Chr3",
    "OZ408686.1": "Chr4", "OZ408687.1": "Chr5",
}


@dataclass(frozen=True)
class Call:
    call_id: int
    chr1: str
    start1: int
    end1: int
    chr2: str
    start2: int
    end2: int
    strand1: str
    strand2: str
    error: float
    max_len: int
    aln_len: int
    raw_exact_key: tuple

    @property
    def arm1(self) -> tuple[str, int, int]:
        return self.chr1, self.start1, self.end1

    @property
    def arm2(self) -> tuple[str, int, int]:
        return self.chr2, self.start2, self.end2

    @property
    def span(self) -> int:
        return self.end1 - self.start1 + self.end2 - self.start2

    @property
    def relative_orientation(self) -> str:
        return "same" if self.strand1 == self.strand2 else "opposite"

    @property
    def group_key(self) -> tuple[str, str, str]:
        return self.chr1, self.chr2, self.relative_orientation

    @property
    def exact_key(self) -> tuple:
        return (
            self.chr1, self.start1, self.end1, self.chr2, self.start2, self.end2,
            self.relative_orientation,
        )


def read_calls(path: Path) -> list[Call]:
    calls = []
    with path.open() as handle:
        for call_id, line in enumerate(handle, 1):
            f = line.rstrip().split("\t")
            arm1 = (f[0], int(f[1]), int(f[2]), f[8])
            arm2 = (f[3], int(f[4]), int(f[5]), f[9])
            if arm2[:3] < arm1[:3]:
                arm1, arm2 = arm2, arm1
            calls.append(Call(
                call_id, arm1[0], arm1[1], arm1[2], arm2[0], arm2[1], arm2[2],
                arm1[3], arm2[3], float(f[7]), int(f[10]), int(f[11]),
                tuple(f[:6] + f[8:10]),
            ))
    return calls


def overlap(a: tuple[str, int, int], b: tuple[str, int, int]) -> int:
    if a[0] != b[0]:
        return 0
    return max(0, min(a[2], b[2]) - max(a[1], b[1]))


def reciprocal_score(a: Call, b: Call) -> float:
    if a.group_key != b.group_key:
        return 0.0
    ov1 = overlap(a.arm1, b.arm1)
    ov2 = overlap(a.arm2, b.arm2)
    if ov1 == 0 or ov2 == 0:
        return 0.0
    return min(
        ov1 / (a.end1 - a.start1), ov1 / (b.end1 - b.start1),
        ov2 / (a.end2 - a.start2), ov2 / (b.end2 - b.start2),
    )


def pairwise_matches(a_calls: list[Call], b_calls: list[Call]):
    b_groups: dict[tuple[str, str, str], list[Call]] = defaultdict(list)
    for call in b_calls:
        b_groups[call.group_key].append(call)
    a_to_b: dict[int, list[tuple[int, float]]] = defaultdict(list)
    b_to_a: dict[int, list[tuple[int, float]]] = defaultdict(list)
    for a in a_calls:
        for b in b_groups.get(a.group_key, []):
            score = reciprocal_score(a, b)
            if score >= 0.25:
                a_to_b[a.call_id].append((b.call_id, score))
                b_to_a[b.call_id].append((a.call_id, score))
    return a_to_b, b_to_a


def merge_intervals(items: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(items):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return [(start, end) for start, end in merged]


def call_union(calls: list[Call]) -> dict[str, list[tuple[int, int]]]:
    by_chr: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for c in calls:
        by_chr[c.chr1].append((c.start1, c.end1))
        by_chr[c.chr2].append((c.start2, c.end2))
    return {chrom: merge_intervals(v) for chrom, v in by_chr.items()}


def interval_lists_overlap(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    i = j = total = 0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return total


def union_bp(intervals: dict[str, list[tuple[int, int]]]) -> int:
    return sum(e - s for values in intervals.values() for s, e in values)


def read_mask_runs(path: Path) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = defaultdict(list)
    chrom = None
    pos = 0
    run_start = None
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if chrom is not None and run_start is not None:
                    result[chrom].append((run_start, pos))
                raw = line[1:].split()[0]
                chrom = ACCESSION_TO_CHR.get(raw, raw)
                pos = 0
                run_start = None
                continue
            for base in line.strip():
                if base.islower() and run_start is None:
                    run_start = pos
                elif not base.islower() and run_start is not None:
                    result[chrom].append((run_start, pos))
                    run_start = None
                pos += 1
    if chrom is not None and run_start is not None:
        result[chrom].append((run_start, pos))
    return result


def subtract_intervals(a: list[tuple[int, int]], b: list[tuple[int, int]]):
    result = []
    j = 0
    for start, end in a:
        pos = start
        while j < len(b) and b[j][1] <= pos:
            j += 1
        k = j
        while k < len(b) and b[k][0] < end:
            if b[k][0] > pos:
                result.append((pos, min(end, b[k][0])))
            pos = max(pos, b[k][1])
            if pos >= end:
                break
            k += 1
        if pos < end:
            result.append((pos, end))
    return result


def interval_overlap(intervals: list[tuple[int, int]], start: int, end: int) -> int:
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


def call_overlap(call: Call, intervals: dict[str, list[tuple[int, int]]]) -> int:
    return (
        interval_overlap(intervals.get(call.chr1, []), call.start1, call.end1)
        + interval_overlap(intervals.get(call.chr2, []), call.start2, call.end2)
    )


def parse_repeat_classes(path: Path) -> dict[str, dict[str, list[tuple[int, int]]]]:
    classes: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    with path.open() as handle:
        for line in handle:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            chrom = ACCESSION_TO_CHR.get(f[4], f[4])
            classes[f[10].split("/")[0]][chrom].append((int(f[5]) - 1, int(f[6])))
    return {
        cls: {chrom: merge_intervals(items) for chrom, items in by_chr.items()}
        for cls, by_chr in classes.items()
    }


def parse_runtime(path: Path) -> tuple[float, int]:
    text = path.read_text(errors="replace")
    match = re.search(r"Elapsed \(wall clock\) time.*?:\s*([0-9:.]+)", text)
    raw = match.group(1) if match else "0"
    parts = [float(x) for x in raw.split(":")]
    seconds = parts[-1] + (parts[-2] * 60 if len(parts) >= 2 else 0)
    if len(parts) == 3:
        seconds += parts[0] * 3600
    rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    return seconds, int(rss.group(1)) if rss else 0


def putative_alignments(path: Path) -> int:
    match = re.search(r"Total alignments:\s*([0-9,]+)", path.read_text())
    return int(match.group(1).replace(",", ""))


def write_tsv(path: Path, header: list[str], rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def has_match(mapping, call_id: int, threshold: float) -> bool:
    return any(score >= threshold for _, score in mapping.get(call_id, []))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_dir
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    analysis = project / "03_biser_segmental_duplication"
    runs = analysis / "runs"

    calls = {
        label: read_calls(runs / dirname / "biser_out")
        for label, dirname in CONDITIONS.items()
    }
    unions = {label: call_union(values) for label, values in calls.items()}

    summary_rows = []
    chromosome_rows = []
    for label, dirname in CONDITIONS.items():
        values = calls[label]
        lengths = [x for c in values for x in (c.end1 - c.start1, c.end2 - c.start2)]
        errors = [c.error for c in values]
        chrom_counts = Counter((c.chr1, c.chr2) for c in values)
        for (chrom1, chrom2), count in sorted(chrom_counts.items()):
            chromosome_rows.append([label, chrom1, chrom2, count])
        runtime, rss = parse_runtime(runs / dirname / "resource_usage.log")
        summary_rows.append([
            label, len(values), sum(1 for c in values if c.chr1 == c.chr2),
            sum(1 for c in values if c.chr1 != c.chr2),
            len((runs / dirname / "biser_out.elem.txt").read_text().splitlines()),
            putative_alignments(runs / dirname / "stdout.log"),
            f"{statistics.mean(lengths):.2f}", f"{statistics.median(lengths):.1f}",
            f"{statistics.quantiles(lengths, n=4)[0]:.1f}",
            f"{statistics.quantiles(lengths, n=4)[2]:.1f}",
            f"{statistics.mean(errors):.4f}", f"{statistics.median(errors):.4f}",
            union_bp(unions[label]), f"{runtime:.2f}", rss,
        ])
    write_tsv(out / "three_run_summary.tsv", [
        "condition", "reported_sd_pairs", "intrachromosomal", "interchromosomal",
        "decomposition_elements", "putative_alignments", "mean_copy_length_bp",
        "median_copy_length_bp", "copy_length_q1_bp", "copy_length_q3_bp",
        "mean_total_error_pct", "median_total_error_pct", "union_copy_coverage_bp",
        "wall_seconds", "max_rss_kb",
    ], summary_rows)
    write_tsv(out / "three_run_chromosome_pair_counts.tsv",
              ["condition", "chromosome_1", "chromosome_2", "events"], chromosome_rows)

    labels = list(CONDITIONS)
    exact_sets = {label: {c.exact_key for c in calls[label]} for label in labels}
    raw_exact_sets = {label: {c.raw_exact_key for c in calls[label]} for label in labels}
    exact_rows = []
    for code in range(1, 8):
        members = [labels[i] for i in range(3) if code & (1 << i)]
        excluded = [labels[i] for i in range(3) if not code & (1 << i)]
        keys = set.intersection(*(exact_sets[x] for x in members))
        if excluded:
            keys -= set.union(*(exact_sets[x] for x in excluded))
        exact_rows.append([int(bool(code & 1)), int(bool(code & 2)), int(bool(code & 4)), len(keys)])
    write_tsv(out / "exact_call_venn.tsv",
              ["source", "gff_extended", "repeatmasker", "exact_calls"], exact_rows)
    raw_exact_rows = []
    for code in range(1, 8):
        members = [labels[i] for i in range(3) if code & (1 << i)]
        excluded = [labels[i] for i in range(3) if not code & (1 << i)]
        keys = set.intersection(*(raw_exact_sets[x] for x in members))
        if excluded:
            keys -= set.union(*(raw_exact_sets[x] for x in excluded))
        raw_exact_rows.append([
            int(bool(code & 1)), int(bool(code & 2)), int(bool(code & 4)), len(keys)
        ])
    write_tsv(out / "raw_exact_call_venn.tsv",
              ["source", "gff_extended", "repeatmasker", "raw_exact_calls"], raw_exact_rows)

    pair_rows = []
    all_maps = {}
    for i, a_label in enumerate(labels):
        for b_label in labels[i + 1:]:
            a_to_b, b_to_a = pairwise_matches(calls[a_label], calls[b_label])
            all_maps[(a_label, b_label)] = (a_to_b, b_to_a)
            a_bp, b_bp = union_bp(unions[a_label]), union_bp(unions[b_label])
            shared_bp = sum(
                interval_lists_overlap(unions[a_label].get(chrom, []), unions[b_label].get(chrom, []))
                for chrom in set(unions[a_label]) | set(unions[b_label])
            )
            exact_raw = len(raw_exact_sets[a_label] & raw_exact_sets[b_label])
            exact_canonical = len(exact_sets[a_label] & exact_sets[b_label])
            row = [a_label, b_label, len(calls[a_label]), len(calls[b_label]),
                   exact_raw, exact_canonical]
            for threshold in (0.5, 0.8):
                a_match = sum(has_match(a_to_b, c.call_id, threshold) for c in calls[a_label])
                b_match = sum(has_match(b_to_a, c.call_id, threshold) for c in calls[b_label])
                row.extend([a_match, f"{100*a_match/len(calls[a_label]):.4f}",
                            b_match, f"{100*b_match/len(calls[b_label]):.4f}"])
            row.extend([a_bp, b_bp, shared_bp, f"{shared_bp/(a_bp+b_bp-shared_bp):.6f}"])
            pair_rows.append(row)
    write_tsv(out / "pairwise_call_stability.tsv", [
        "condition_a", "condition_b", "calls_a", "calls_b", "exact_shared_raw",
        "exact_shared_canonical",
        "a_matched_ro50", "a_matched_ro50_pct", "b_matched_ro50", "b_matched_ro50_pct",
        "a_matched_ro80", "a_matched_ro80_pct", "b_matched_ro80", "b_matched_ro80_pct",
        "a_union_bp", "b_union_bp", "union_coverage_intersection_bp", "coverage_jaccard",
    ], pair_rows)

    gff_label = "gff_extended_softmask"
    rm_label = "repeatmasker_arabidopsis_softmask"
    source_label = "source_softmask"
    gff_to_rm, rm_to_gff = all_maps[(gff_label, rm_label)]
    source_to_gff, gff_to_source = all_maps[(source_label, gff_label)]
    mapping_rows = []
    for c in calls[gff_label]:
        matches = sorted(gff_to_rm.get(c.call_id, []), key=lambda x: x[1], reverse=True)
        best_id, best_score = matches[0] if matches else ("", 0.0)
        mapping_rows.append([
            c.call_id, best_id, f"{best_score:.6f}", int(c.exact_key in exact_sets[rm_label]),
            int(best_score >= 0.5), int(best_score >= 0.8),
            len([1 for _, score in matches if score >= 0.5]),
            len([1 for _, score in matches if score >= 0.8]),
        ])
    write_tsv(out / "gff_to_repeatmasker_call_mapping.tsv", [
        "gff_call_id", "best_repeatmasker_call_id", "best_reciprocal_score", "exact_match",
        "matched_ro50", "matched_ro80", "repeatmasker_matches_ro50", "repeatmasker_matches_ro80",
    ], mapping_rows)

    source_fa = project / "01_reference/prepared_data/TAIR12.Col-CC.annotation_softmasked.fa"
    rm_fa = project / "01_reference/repeatmasker_arabidopsis/GCA_978657495.1_TAIR12_genomic.uppercase.fna.masked"
    gff_mask = read_mask_runs(source_fa)
    rm_mask = read_mask_runs(rm_fa)
    chroms = [f"Chr{i}" for i in range(1, 6)]
    rm_new = {chrom: subtract_intervals(rm_mask[chrom], gff_mask[chrom]) for chrom in chroms}
    rm_reopened = {chrom: subtract_intervals(gff_mask[chrom], rm_mask[chrom]) for chrom in chroms}
    repeat_classes = parse_repeat_classes(
        project / "01_reference/repeatmasker_arabidopsis/GCA_978657495.1_TAIR12_genomic.uppercase.fna.out"
    )

    context_rows = []
    class_rows = []
    for label, values, mapping in [
        (gff_label, calls[gff_label], gff_to_rm),
        (rm_label, calls[rm_label], rm_to_gff),
    ]:
        groups: dict[str, list[Call]] = defaultdict(list)
        for c in values:
            best = max((score for _, score in mapping.get(c.call_id, [])), default=0.0)
            if c.exact_key in exact_sets[gff_label] & exact_sets[rm_label]:
                group = "exact"
            elif best >= 0.8:
                group = "nonexact_ro80"
            elif best >= 0.5:
                group = "ro50_only"
            elif best >= 0.25:
                group = "ro25_only"
            else:
                group = "unmatched_ro25"
            groups[group].append(c)
        for group, subset in sorted(groups.items()):
            span = sum(c.span for c in subset)
            new_bp = sum(call_overlap(c, rm_new) for c in subset)
            reopened_bp = sum(call_overlap(c, rm_reopened) for c in subset)
            context_rows.append([
                label, group, len(subset), span, new_bp, f"{100*new_bp/span:.4f}",
                reopened_bp, f"{100*reopened_bp/span:.4f}",
                f"{statistics.median(c.error for c in subset):.4f}",
                sum(c.chr1 == c.chr2 for c in subset), sum(c.chr1 != c.chr2 for c in subset),
            ])
            for repeat_class, intervals in repeat_classes.items():
                bp = sum(call_overlap(c, intervals) for c in subset)
                if bp:
                    class_rows.append([
                        label, group, repeat_class, bp, f"{100*bp/span:.4f}"
                    ])
    write_tsv(out / "call_stability_mask_context.tsv", [
        "condition", "stability_group", "calls", "two_arm_span_bp", "rm_new_mask_bp",
        "rm_new_mask_pct_span", "rm_reopened_bp", "rm_reopened_pct_span",
        "median_total_error_pct", "intrachromosomal", "interchromosomal",
    ], context_rows)
    write_tsv(out / "call_stability_repeat_class_overlap.tsv", [
        "condition", "stability_group", "repeat_class", "overlap_bp",
        "overlap_pct_two_arm_span",
    ], class_rows)

    gff_source_rm50 = sum(
        has_match(gff_to_source, c.call_id, 0.5) and has_match(gff_to_rm, c.call_id, 0.5)
        for c in calls[gff_label]
    )
    gff_source_rm80 = sum(
        has_match(gff_to_source, c.call_id, 0.8) and has_match(gff_to_rm, c.call_id, 0.8)
        for c in calls[gff_label]
    )
    write_tsv(out / "three_way_gff_anchored_stability.tsv",
              ["criterion", "gff_calls", "pct_of_gff_calls"], [
                  ["exact_in_all_three", len(exact_sets[source_label] & exact_sets[gff_label] & exact_sets[rm_label]),
                   f"{100*len(exact_sets[source_label] & exact_sets[gff_label] & exact_sets[rm_label])/len(calls[gff_label]):.4f}"],
                  ["matched_to_source_and_repeatmasker_ro50", gff_source_rm50,
                   f"{100*gff_source_rm50/len(calls[gff_label]):.4f}"],
                  ["matched_to_source_and_repeatmasker_ro80", gff_source_rm80,
                   f"{100*gff_source_rm80/len(calls[gff_label]):.4f}"],
              ])

    event_files = [
        ("strict_two_copy_events", project / "06_sd_age_tracing_preparation/event_first_reanalysis/network_age_blind/strict_two_copy_events.tsv"),
        ("event_first_threshold_stable", project / "06_sd_age_tracing_preparation/event_first_reanalysis/events/event_first_threshold_stable.tsv"),
        ("event_first_time1_time3_stable", project / "06_sd_age_tracing_preparation/event_first_reanalysis/events/event_first_time1_time3_threshold_stable.tsv"),
        ("age_free_pd_events", project / "15_age_free_pd_sequence_variation/inputs/age_free_pd_events.tsv"),
    ]
    event_rows = []
    strata_rows = []
    age_free_context: dict[str, list[Call]] = defaultdict(list)
    for dataset, path in event_files:
        with path.open() as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        metrics = Counter()
        for rec in records:
            ids = [int(x) for x in rec["source_call_ids"].split(",")]
            rep = int(rec["representative_call_id"])
            metrics["representative_exact"] += calls[gff_label][rep - 1].exact_key in exact_sets[rm_label]
            for threshold, suffix in [(0.5, "ro50"), (0.8, "ro80")]:
                flags = [has_match(gff_to_rm, call_id, threshold) for call_id in ids]
                metrics[f"representative_{suffix}"] += has_match(gff_to_rm, rep, threshold)
                metrics[f"any_source_call_{suffix}"] += any(flags)
                metrics[f"all_source_calls_{suffix}"] += all(flags)
            if dataset == "age_free_pd_events":
                if has_match(gff_to_rm, rep, 0.8):
                    stability = "representative_ro80"
                elif has_match(gff_to_rm, rep, 0.5):
                    stability = "representative_ro50_only"
                else:
                    stability = "representative_unmatched_ro50"
                age_free_context[stability].append(calls[gff_label][rep - 1])
        row = [dataset, len(records)]
        for key in ["representative_exact", "representative_ro50", "representative_ro80",
                    "any_source_call_ro50", "all_source_calls_ro50",
                    "any_source_call_ro80", "all_source_calls_ro80"]:
            row.extend([metrics[key], f"{100*metrics[key]/len(records):.4f}"])
        event_rows.append(row)
        for stratifier in ["strict_age_bin", "former_strict_age_bin",
                           "age_free_PD_evidence_class", "age_free_single_copy_support_nodes"]:
            if not records or stratifier not in records[0]:
                continue
            levels: dict[str, list[dict[str, str]]] = defaultdict(list)
            for rec in records:
                levels[rec[stratifier]].append(rec)
            for level, subset in sorted(levels.items()):
                rep50 = sum(
                    has_match(gff_to_rm, int(rec["representative_call_id"]), 0.5)
                    for rec in subset
                )
                rep80 = sum(
                    has_match(gff_to_rm, int(rec["representative_call_id"]), 0.8)
                    for rec in subset
                )
                strata_rows.append([
                    dataset, stratifier, level, len(subset), rep50,
                    f"{100*rep50/len(subset):.4f}", rep80,
                    f"{100*rep80/len(subset):.4f}",
                ])
    write_tsv(out / "downstream_event_retention.tsv", [
        "dataset", "events", "representative_exact", "representative_exact_pct",
        "representative_ro50", "representative_ro50_pct", "representative_ro80",
        "representative_ro80_pct", "any_source_call_ro50", "any_source_call_ro50_pct",
        "all_source_calls_ro50", "all_source_calls_ro50_pct", "any_source_call_ro80",
        "any_source_call_ro80_pct", "all_source_calls_ro80", "all_source_calls_ro80_pct",
    ], event_rows)
    write_tsv(out / "downstream_event_retention_stratified.tsv", [
        "dataset", "stratifier", "level", "events", "representative_ro50",
        "representative_ro50_pct", "representative_ro80", "representative_ro80_pct",
    ], strata_rows)

    age_context_rows = []
    age_class_rows = []
    for stability, subset in sorted(age_free_context.items()):
        span = sum(c.span for c in subset)
        new_bp = sum(call_overlap(c, rm_new) for c in subset)
        reopened_bp = sum(call_overlap(c, rm_reopened) for c in subset)
        age_context_rows.append([
            stability, len(subset), span, new_bp, f"{100*new_bp/span:.4f}",
            reopened_bp, f"{100*reopened_bp/span:.4f}",
            f"{statistics.median(c.error for c in subset):.4f}",
        ])
        for repeat_class, intervals in repeat_classes.items():
            bp = sum(call_overlap(c, intervals) for c in subset)
            if bp:
                age_class_rows.append([
                    stability, repeat_class, bp, f"{100*bp/span:.4f}"
                ])
    write_tsv(out / "age_free_event_mask_context.tsv", [
        "stability_group", "events", "representative_two_arm_span_bp", "rm_new_mask_bp",
        "rm_new_mask_pct_span", "rm_reopened_bp", "rm_reopened_pct_span",
        "median_total_error_pct",
    ], age_context_rows)
    write_tsv(out / "age_free_event_repeat_class_overlap.tsv", [
        "stability_group", "repeat_class", "overlap_bp", "overlap_pct_span",
    ], age_class_rows)

    age_free_path = project / "15_age_free_pd_sequence_variation/inputs/age_free_pd_events.tsv"
    with age_free_path.open() as handle:
        age_free_records = list(csv.DictReader(handle, delimiter="\t"))
    event_mapping_rows = []
    for rec in age_free_records:
        rep = int(rec["representative_call_id"])
        call = calls[gff_label][rep - 1]
        matches = sorted(gff_to_rm.get(rep, []), key=lambda x: x[1], reverse=True)
        best_id, best_score = matches[0] if matches else ("", 0.0)
        new_bp = call_overlap(call, rm_new)
        reopened_bp = call_overlap(call, rm_reopened)
        class_values = {
            cls: call_overlap(call, repeat_classes.get(cls, {}))
            for cls in ["rRNA", "Satellite", "LTR", "DNA", "RC", "LINE",
                        "Simple_repeat", "Low_complexity"]
        }
        event_mapping_rows.append([
            rec["event_id"], rep, best_id, f"{best_score:.6f}",
            int(best_score >= 0.5), int(best_score >= 0.8),
            rec.get("former_strict_age_bin", ""), rec.get("age_free_PD_evidence_class", ""),
            rec.get("age_free_single_copy_support_nodes", ""), rec.get("primary_pattern", ""),
            call.span, new_bp, f"{100*new_bp/call.span:.4f}", reopened_bp,
            f"{100*reopened_bp/call.span:.4f}", f"{call.error:.4f}",
            *[class_values[x] for x in ["rRNA", "Satellite", "LTR", "DNA", "RC", "LINE",
                                           "Simple_repeat", "Low_complexity"]],
        ])
    event_mapping_header = [
        "event_id", "representative_gff_call_id", "best_repeatmasker_call_id",
        "best_reciprocal_score", "mask_stable_ro50", "mask_stable_ro80",
        "former_strict_age_bin", "age_free_PD_evidence_class",
        "age_free_single_copy_support_nodes", "primary_pattern",
        "representative_two_arm_span_bp", "rm_new_mask_bp", "rm_new_mask_pct_span",
        "rm_reopened_bp", "rm_reopened_pct_span", "gff_biser_total_error_pct",
        "rRNA_bp", "Satellite_bp", "LTR_bp", "DNA_bp", "RC_bp", "LINE_bp",
        "Simple_repeat_bp", "Low_complexity_bp",
    ]
    write_tsv(out / "age_free_event_repeatmasker_mapping.tsv",
              event_mapping_header, event_mapping_rows)
    write_tsv(out / "age_free_mask_stable_ro50.tsv", event_mapping_header,
              [row for row in event_mapping_rows if row[4] == 1])
    write_tsv(out / "age_free_mask_stable_ro80.tsv", event_mapping_header,
              [row for row in event_mapping_rows if row[5] == 1])
    write_tsv(out / "age_free_mask_sensitive_lt_ro50.tsv", event_mapping_header,
              [row for row in event_mapping_rows if row[4] == 0])

    stable_sets = {
        "all_current_events": {row[0] for row in event_mapping_rows},
        "mask_stable_ro50": {row[0] for row in event_mapping_rows if row[4] == 1},
        "mask_stable_ro80": {row[0] for row in event_mapping_rows if row[5] == 1},
        "mask_sensitive_lt_ro50": {row[0] for row in event_mapping_rows if row[4] == 0},
    }
    endpoint_inputs = [
        (
            "three_aligner_SNP_ge200",
            project / "15_age_free_pd_sequence_variation/snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge200.event_metrics.tsv",
            "local_MSA_callable_sites", "P_specific_SNP", "D_specific_SNP",
        ),
        (
            "independent_microindel_primary",
            project / "15_age_free_pd_sequence_variation/microindel_local_msa/event_level_PD_microindel_rates.tsv",
            "ancestral_callable_bp", "P_branch_indels", "D_branch_indels",
        ),
    ]
    endpoint_rows = []
    for endpoint, path, callable_col, p_col, d_col in endpoint_inputs:
        with path.open() as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        for subset_name, wanted in stable_sets.items():
            subset = [rec for rec in records if rec["event_id"] in wanted]
            p_values = [int(rec[p_col]) for rec in subset]
            d_values = [int(rec[d_col]) for rec in subset]
            p_total, d_total = sum(p_values), sum(d_values)
            d_greater = sum(d > p for p, d in zip(p_values, d_values))
            p_greater = sum(p > d for p, d in zip(p_values, d_values))
            tied = len(subset) - d_greater - p_greater
            sign_p = (
                binomtest(d_greater, d_greater + p_greater, 0.5).pvalue
                if d_greater + p_greater else 1.0
            )
            site_p = binomtest(d_total, d_total + p_total, 0.5).pvalue if d_total + p_total else 1.0
            differences = [d - p for p, d in zip(p_values, d_values)]
            wilcoxon_p = (
                wilcoxon(
                    differences,
                    zero_method="wilcox",
                    alternative="two-sided",
                    method="approx",
                ).pvalue
                if any(differences) else 1.0
            )
            endpoint_rows.append([
                endpoint, subset_name, len(subset), sum(int(rec[callable_col]) for rec in subset),
                p_total, d_total, f"{d_total/p_total:.6f}" if p_total else "Inf",
                d_greater, p_greater, tied, d_greater + p_greater,
                f"{sign_p:.10g}", f"{wilcoxon_p:.10g}", f"{site_p:.10g}",
            ])
    write_tsv(out / "mask_stable_downstream_endpoints.tsv", [
        "endpoint", "mask_subset", "events", "callable_bp_or_sites", "P_count",
        "D_count", "D_to_P_ratio", "events_D_greater", "events_P_greater",
        "events_tied", "non_tied_events", "event_sign_test_p", "paired_wilcoxon_p",
        "site_count_binomial_p_descriptive",
    ], endpoint_rows)


if __name__ == "__main__":
    main()
