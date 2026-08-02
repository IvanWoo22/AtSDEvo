#!/usr/bin/env python3
"""Validate MCScanX non-detections with nucleotide hits and flanking anchors.

The script deliberately uses three-valued missing-side calls:
  R: an independent second target locus is rescued by sequence evidence;
  A: a distinct, assembly-callable orthologous window exists but has no hit;
  U: absence cannot be tested (ambiguous/no anchors, assembly gap, or only a
     shared ancestral homolog).

An A call is a candidate assembly-level absence, not proof of a biological
deletion; validation in an independent assembly/read set remains desirable.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
GENOMES = {
    "Alyrata": "Phytozome/PhytozomeV12/Alyrata/assembly/Alyrata_384_v1.fa.gz",
    "Bstricta": "Phytozome/PhytozomeV10/Bstricta/assembly/Bstricta_278_v1.fa.gz",
    "Dstrictus": "Phytozome/PhytozomeV13/Dstrictus/v2.1/assembly/Dstrictus_582_v2.0.fa.gz",
    "Cviolacea": "Phytozome/PhytozomeV13/Cviolacea/v2.1/assembly/Cviolacea_585_v2.0.fa.gz",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def percentile(values: list[float], p: float) -> float:
    if not values:
        return 0.0
    values = sorted(values)
    position = (len(values) - 1) * p
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    return values[lower] * (upper - position) + values[upper] * (position - lower)


def merged_bp(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    left, right = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= right:
            right = max(right, end)
        else:
            total += right - left
            left, right = start, end
    return total + right - left


def parse_blast(path: Path) -> dict[str, list[dict[str, object]]]:
    """Cluster HSPs by query, subject and target proximity."""
    hsps: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 14:
                continue
            qid, qlen, sid, slen = f[:4]
            s1, s2 = int(f[10]), int(f[11])
            hsps[(qid, sid)].append(
                {
                    "q0": min(int(f[8]), int(f[9])) - 1,
                    "q1": max(int(f[8]), int(f[9])),
                    "s0": min(s1, s2) - 1,
                    "s1": max(s1, s2),
                    "pident": float(f[4]),
                    "alen": int(f[5]),
                    "bitscore": float(f[13]),
                    "strand": "+" if s1 <= s2 else "-",
                }
            )
    output: dict[str, list[dict[str, object]]] = defaultdict(list)
    for (qid, sid), records in hsps.items():
        records.sort(key=lambda x: (x["s0"], x["s1"]))
        clusters: list[list[dict[str, object]]] = []
        for record in records:
            if (
                not clusters
                or record["s0"] - max(x["s1"] for x in clusters[-1]) > 50000
            ):
                clusters.append([record])
            else:
                clusters[-1].append(record)
        for cluster in clusters:
            aligned = sum(x["alen"] for x in cluster)
            output[qid].append(
                {
                    "target": sid,
                    "start": min(x["s0"] for x in cluster),
                    "end": max(x["s1"] for x in cluster),
                    "query_union_bp": merged_bp(
                        [(x["q0"], x["q1"]) for x in cluster]
                    ),
                    "identity": sum(x["pident"] * x["alen"] for x in cluster)
                    / max(1, aligned),
                    "bitscore": sum(x["bitscore"] for x in cluster),
                    "hsp_count": len(cluster),
                    "method": "dc-megablast",
                }
            )
    return output


def parse_paf(path: Path) -> dict[str, list[dict[str, object]]]:
    output: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            output[f[0]].append(
                {
                    "target": f[5],
                    "start": int(f[7]),
                    "end": int(f[8]),
                    "identity": 100 * int(f[9]) / max(1, int(f[10])),
                    "query_union_bp": int(f[3]) - int(f[2]),
                    "mapq": int(f[11]),
                    "method": "minimap2-asm20",
                }
            )
    return output


def overlaps(a: dict[str, object], b: dict[str, object], pad=50000) -> bool:
    return (
        a["target"] == b["target"]
        and int(a["start"]) <= int(b["end"]) + pad
        and int(b["start"]) <= int(a["end"]) + pad
    )


def read_bed(path: Path) -> tuple[dict[str, tuple[str, int, int]], dict[str, list]]:
    genes = {}
    chromosomes: dict[str, list] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            chrom, start, end, gene = f[0], int(f[1]), int(f[2]), f[3]
            genes[gene] = (chrom, start, end)
            chromosomes[chrom].append((start, end, gene))
    for chrom in chromosomes:
        chromosomes[chrom].sort()
    return genes, chromosomes


def parse_collinear_anchors(
    path: Path, atha_genes: set[str], out_genes: set[str]
) -> dict[str, set[str]]:
    anchors: dict[str, set[str]] = defaultdict(set)
    pair = re.compile(r"^\s*\d+-\s*\d+:\s+(\S+)\s+(\S+)\s+")
    with path.open() as handle:
        for line in handle:
            match = pair.match(line)
            if not match:
                continue
            one, two = match.groups()
            if one in atha_genes and two in out_genes:
                anchors[one].add(two)
            elif two in atha_genes and one in out_genes:
                anchors[two].add(one)
    return anchors


def nearest_anchor_genes(
    chrom_genes: list[tuple[int, int, str]],
    locus_start: int,
    locus_end: int,
    anchors: dict[str, set[str]],
    max_distance=500000,
    each_side=3,
) -> tuple[list[str], list[str]]:
    starts = [x[0] for x in chrom_genes]
    index = bisect_left(starts, locus_start)
    left, right = [], []
    for item in reversed(chrom_genes[:index]):
        if locus_start - item[1] > max_distance:
            break
        if item[2] in anchors:
            left.append(item[2])
            if len(left) == each_side:
                break
    for item in chrom_genes[index:]:
        if item[0] < locus_end:
            continue
        if item[0] - locus_end > max_distance:
            break
        if item[2] in anchors:
            right.append(item[2])
            if len(right) == each_side:
                break
    return left, right


def infer_window(
    left: list[str],
    right: list[str],
    anchors: dict[str, set[str]],
    out_genes: dict[str, tuple[str, int, int]],
) -> dict[str, object]:
    candidates = []
    for li, left_gene in enumerate(left):
        for ri, right_gene in enumerate(right):
            for out_left in anchors[left_gene]:
                for out_right in anchors[right_gene]:
                    if out_left not in out_genes or out_right not in out_genes:
                        continue
                    lc, ls, le = out_genes[out_left]
                    rc, rs, re = out_genes[out_right]
                    if lc != rc:
                        continue
                    start, end = min(ls, rs), max(le, re)
                    span = end - start
                    if span > 5_000_000:
                        continue
                    candidates.append(
                        {
                            "target": lc,
                            "start": start,
                            "end": end,
                            "left_atha": left_gene,
                            "right_atha": right_gene,
                            "left_out": out_left,
                            "right_out": out_right,
                            "rank": li + ri,
                            "span": span,
                        }
                    )
    if not candidates:
        return {"window_status": "NO_PAIRED_FLANK_ANCHORS"}
    candidates.sort(key=lambda x: (x["rank"], x["span"]))
    best = candidates[0]
    competing = [
        x
        for x in candidates[1:]
        if x["rank"] == best["rank"]
        and not (
            x["target"] == best["target"]
            and x["start"] <= best["end"] + 100000
            and best["start"] <= x["end"] + 100000
        )
    ]
    if competing:
        best["window_status"] = "AMBIGUOUS_PAIRED_ANCHORS"
    else:
        best["window_status"] = "PAIRED_ANCHORS"
    return best


def read_fasta_gz(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = ""
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in seqs.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    analysis = args.project / "06_sd_age_tracing_preparation/event_first_reanalysis"
    work = analysis / "missing_state_validation"
    manifest = read_tsv(work / "query_manifest.tsv")
    candidates = read_tsv(work / "candidate_event_species.tsv")
    events = {
        row["event_id"]: row
        for row in read_tsv(analysis / "events/event_first_events.tsv")
    }
    manifest_by_species = defaultdict(dict)
    for row in manifest:
        manifest_by_species[row["species"]][row["query_id"]] = row

    atha_gene_coords, atha_chromosomes = read_bed(
        args.project / "05_mcscanx_synteny/inputs/Atha.bed"
    )
    all_hit_rows, all_locus_rows, all_event_rows, threshold_rows = [], [], [], []
    for species in SPECIES:
        species_manifest = manifest_by_species[species]
        blast = parse_blast(work / f"dc_megablast/{species}.tsv")
        paf = parse_paf(work / f"minimap2_asm20/{species}.paf")

        # Calibrate a high-confidence tier to the 10th percentile of state-3
        # positive controls, while retaining a deliberately sensitive tier.
        control_cov, control_ident = [], []
        for qid, meta in species_manifest.items():
            if meta["candidate_role"] != "state3_positive_control":
                continue
            denom = max(1, int(meta["query_uppercase_ACGT_bp"]))
            if blast.get(qid):
                best = max(blast[qid], key=lambda x: x["bitscore"])
                control_cov.append(min(1.0, best["query_union_bp"] / denom))
                control_ident.append(best["identity"])
            else:
                control_cov.append(0.0)
                control_ident.append(0.0)
        strong_cov = max(0.10, min(0.50, percentile(control_cov, 0.10)))
        strong_ident = max(65.0, min(85.0, percentile(control_ident, 0.10)))
        threshold_rows.append(
            {
                "species": species,
                "control_queries": len(control_cov),
                "sensitive_min_coverage": "0.100000",
                "sensitive_min_identity_pct": "65.000",
                "strong_min_coverage": f"{strong_cov:.6f}",
                "strong_min_identity_pct": f"{strong_ident:.3f}",
                "control_blast_detected": sum(x > 0 for x in control_cov),
                "control_strong_pass": sum(
                    c >= strong_cov and i >= strong_ident
                    for c, i in zip(control_cov, control_ident)
                ),
            }
        )

        accepted: dict[str, list[dict[str, object]]] = defaultdict(list)
        for qid, meta in species_manifest.items():
            denom = max(1, int(meta["query_uppercase_ACGT_bp"]))
            for hit_index, hit in enumerate(
                sorted(blast.get(qid, []), key=lambda x: x["bitscore"], reverse=True),
                1,
            ):
                coverage = min(1.0, hit["query_union_bp"] / denom)
                sensitive = (
                    coverage >= 0.10
                    and hit["identity"] >= 65.0
                    and hit["query_union_bp"] >= 100
                    and hit["bitscore"] >= 100
                )
                strong = (
                    coverage >= strong_cov
                    and hit["identity"] >= strong_ident
                    and hit["query_union_bp"] >= 100
                    and hit["bitscore"] >= 100
                )
                mm_support = any(
                    overlaps(hit, mm, pad=10000)
                    and mm["mapq"] >= 10
                    and mm["identity"] >= 65
                    for mm in paf.get(qid, [])
                )
                row = {
                    "species": species,
                    "query_id": qid,
                    "event_id": meta["event_id"],
                    "event_locus": meta["event_locus"],
                    "hit_rank": hit_index,
                    "target": hit["target"],
                    "target_start": hit["start"],
                    "target_end": hit["end"],
                    "query_coverage_of_unmasked_bp": f"{coverage:.6f}",
                    "weighted_identity_pct": f"{hit['identity']:.3f}",
                    "bitscore_sum": f"{hit['bitscore']:.1f}",
                    "hsp_count": hit["hsp_count"],
                    "sensitive_pass": "PASS" if sensitive else "FAIL",
                    "control_calibrated_strong_pass": "PASS" if strong else "FAIL",
                    "minimap2_strict_support": "PASS" if mm_support else "FAIL",
                }
                all_hit_rows.append(row)
                if sensitive:
                    kept = dict(hit)
                    kept["coverage"] = coverage
                    kept["strong"] = strong
                    kept["minimap"] = mm_support
                    accepted[qid].append(kept)

        out_gene_coords, _ = read_bed(
            args.project / f"05_mcscanx_synteny/inputs/{species}.bed"
        )
        anchors = parse_collinear_anchors(
            args.project
            / f"05_mcscanx_synteny/results/Atha_{species}/Atha_{species}.collinearity",
            set(atha_gene_coords),
            set(out_gene_coords),
        )
        genome = read_fasta_gz(
            args.project
            / "04_outgroup_resources/primary_nodes_jgi_20260724/raw"
            / GENOMES[species]
        )
        windows: dict[tuple[str, str], dict[str, object]] = {}
        for candidate in (x for x in candidates if x["species"] == species):
            event = events[candidate["event_id"]]
            for locus in ("locus_A", "locus_B"):
                chrom = event[f"{locus}_chrom"]
                start = int(event[f"{locus}_representative_start"])
                end = int(event[f"{locus}_representative_end"])
                left, right = nearest_anchor_genes(
                    atha_chromosomes[chrom], start, end, anchors
                )
                window = infer_window(left, right, anchors, out_gene_coords)
                if "target" in window and window["target"] in genome:
                    seq = genome[window["target"]][window["start"] : window["end"]]
                    non_n = sum(x in "ACGT" for x in seq) / max(1, len(seq))
                    window["non_N_fraction"] = non_n
                    window["assembly_callable"] = (
                        window["window_status"] == "PAIRED_ANCHORS"
                        and non_n >= 0.90
                        and window.get("rank") == 0
                        and len(seq) <= 1_000_000
                    )
                else:
                    window["non_N_fraction"] = 0.0
                    window["assembly_callable"] = False
                windows[(candidate["event_id"], locus)] = window

        for candidate in (x for x in candidates if x["species"] == species):
            eid, state = candidate["event_id"], candidate["mcscan_state"]
            locus_hits = {
                locus: accepted.get(f"{eid}__{locus}", [])
                for locus in ("locus_A", "locus_B")
            }
            detected = {
                "locus_A": state in ("1", "3"),
                "locus_B": state in ("2", "3"),
            }
            statuses = {}
            for locus, other in (("locus_A", "locus_B"), ("locus_B", "locus_A")):
                window = windows[(eid, locus)]
                hits = locus_hits[locus]
                hits_in_window = [
                    h
                    for h in hits
                    if "target" in window and overlaps(h, window, pad=0)
                ]
                independent_hits = [
                    h
                    for h in hits
                    if not any(overlaps(h, x) for x in locus_hits[other])
                ]
                high_independent_hits = [
                    h for h in independent_hits if h["strong"] or h["minimap"]
                ]
                anchored_independent_hits = [
                    h
                    for h in high_independent_hits
                    if "target" in window and overlaps(h, window, pad=0)
                ]
                shared_hits = [
                    h
                    for h in hits
                    if any(overlaps(h, x) for x in locus_hits[other])
                ]
                if detected[locus]:
                    status = "P_MCSCAN"
                elif anchored_independent_hits:
                    status = "R_ANCHORED_INDEPENDENT_SEQUENCE_LOCUS"
                elif window.get("assembly_callable") and not hits_in_window:
                    other_window = windows[(eid, other)]
                    if (
                        "target" in other_window
                        and overlaps(window, other_window, pad=100000)
                    ):
                        status = "U_SAME_ORTHOLOGOUS_WINDOW"
                    else:
                        status = "A_CANDIDATE_CALLABLE_WINDOW_NO_HIT"
                elif high_independent_hits:
                    status = "U_UNANCHORED_SECOND_LOCUS_CANDIDATE"
                elif shared_hits:
                    status = "U_SHARED_ANCESTRAL_HOMOLOG_ONLY"
                else:
                    status = "U_NO_CALLABLE_ORTHOLOGOUS_WINDOW"
                statuses[locus] = status
                all_locus_rows.append(
                    {
                        "event_id": eid,
                        "species": species,
                        "mcscan_state": state,
                        "event_locus": locus,
                        "mcscan_detected": detected[locus],
                        "validation_status": status,
                        "sensitive_hit_clusters": len(hits),
                        "strong_hit_clusters": sum(h["strong"] for h in hits),
                        "minimap_supported_clusters": sum(h["minimap"] for h in hits),
                        "independent_from_other_query_clusters": len(independent_hits),
                        "high_confidence_independent_clusters": len(
                            high_independent_hits
                        ),
                        "anchored_independent_clusters": len(
                            anchored_independent_hits
                        ),
                        "shared_with_other_query_clusters": len(shared_hits),
                        "window_status": window.get("window_status", ""),
                        "window_target": window.get("target", ""),
                        "window_start": window.get("start", ""),
                        "window_end": window.get("end", ""),
                        "window_non_N_fraction": (
                            f"{window.get('non_N_fraction', 0):.6f}"
                        ),
                        "window_anchor_rank_sum": window.get("rank", ""),
                        "window_span_bp": window.get("span", ""),
                        "left_atha_anchor": window.get("left_atha", ""),
                        "right_atha_anchor": window.get("right_atha", ""),
                        "left_outgroup_anchor": window.get("left_out", ""),
                        "right_outgroup_anchor": window.get("right_out", ""),
                    }
                )
            missing = [
                statuses[x]
                for x in ("locus_A", "locus_B")
                if not detected[x]
            ]
            all_event_rows.append(
                {
                    "event_id": eid,
                    "species": species,
                    "candidate_role": candidate["candidate_role"],
                    "original_mcscan_state": state,
                    "locus_A_validation": statuses["locus_A"],
                    "locus_B_validation": statuses["locus_B"],
                    "missing_loci_n": len(missing),
                    "missing_rescued_n": sum(x.startswith("R_") for x in missing),
                    "missing_candidate_absent_n": sum(x.startswith("A_") for x in missing),
                    "missing_uncallable_n": sum(x.startswith("U_") for x in missing),
                    "event_validation": (
                        "MCSCAN_STATE3_CONTROL"
                        if state == "3"
                        else "HAS_SEQUENCE_RESCUE"
                        if any(x.startswith("R_") for x in missing)
                        else "HAS_CANDIDATE_ABSENCE"
                        if any(x.startswith("A_") for x in missing)
                        else "UNCALLABLE_MISSING_STATE"
                    ),
                }
            )

    write_tsv(args.output / "species_control_calibrated_thresholds.tsv", threshold_rows)
    write_tsv(args.output / "sequence_hit_clusters.tsv", all_hit_rows)
    write_tsv(args.output / "locus_missing_validation.tsv", all_locus_rows)
    write_tsv(args.output / "event_species_missing_validation.tsv", all_event_rows)
    summary = []
    for species in SPECIES:
        rows = [
            x
            for x in all_event_rows
            if x["species"] == species and x["candidate_role"] == "missing_state_target"
        ]
        labels = sorted({x["event_validation"] for x in rows})
        for label in labels:
            summary.append(
                {
                    "species": species,
                    "event_validation": label,
                    "event_species_n": sum(x["event_validation"] == label for x in rows),
                }
            )
    write_tsv(args.output / "validation_summary.tsv", summary)

    validation_lookup = {
        (row["event_id"], row["species"]): row for row in all_event_rows
    }
    matrix_rows = []
    for event_id, event in events.items():
        output = {
            "event_id": event_id,
            "original_primary_pattern": event["primary_pattern"],
        }
        pattern = []
        contains_u = False
        for species in SPECIES:
            original = event[f"{species}_state"]
            if original == "3":
                pair = "P/P"
            else:
                validation = validation_lookup[(event_id, species)]
                symbols = []
                for locus in ("locus_A", "locus_B"):
                    status = validation[f"{locus}_validation"]
                    symbol = (
                        "P"
                        if status.startswith(("P_", "R_"))
                        else "A"
                        if status.startswith("A_")
                        else "U"
                    )
                    symbols.append(symbol)
                pair = "/".join(symbols)
            output[f"{species}_validated_pair_state"] = pair
            pattern.append(pair)
            contains_u |= "U" in pair
        output["validated_pattern"] = "|".join(pattern)
        output["validation_completeness"] = (
            "COMPLETE_NO_U" if not contains_u else "INCOMPLETE_HAS_U"
        )
        matrix_rows.append(output)
    write_tsv(args.output / "validated_event_species_matrix.tsv", matrix_rows)
    pattern_counts: dict[tuple[str, str], int] = defaultdict(int)
    for row in matrix_rows:
        pattern_counts[
            (row["validated_pattern"], row["validation_completeness"])
        ] += 1
    write_tsv(
        args.output / "validated_pattern_counts.tsv",
        [
            {
                "validated_pattern": key[0],
                "validation_completeness": key[1],
                "event_n": count,
            }
            for key, count in sorted(
                pattern_counts.items(), key=lambda x: (-x[1], x[0])
            )
        ],
    )


if __name__ == "__main__":
    main()
