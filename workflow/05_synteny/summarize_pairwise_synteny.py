#!/usr/bin/env python3
"""Summarize pairwise TAIR12-outgroup MCScanX results."""

import argparse
import csv
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


BLOCK_RE = re.compile(
    r"^## Alignment (\d+): score=(\S+) e_value=(\S+) N=(\d+) (\S+)&(\S+) (plus|minus)"
)
PAIR_RE = re.compile(r"^\s*\d+-\s*\d+:\s+(\S+)\s+(\S+)\s+(\S+)")


def gff_ids(path):
    ids = []
    chrom = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip().split("\t")
            if len(parts) >= 2:
                ids.append(parts[1])
                chrom[parts[1]] = parts[0].removeprefix("Atha_")
    return ids, chrom


def time_value(path, label):
    if not path.exists():
        return "NA"
    value = "NA"
    with open(path) as handle:
        for line in handle:
            if label in line and ": " in line:
                value = line.rsplit(": ", 1)[1].strip()
    return value


def elapsed_seconds(value):
    if value == "NA":
        return "NA"
    fields = value.split(":")
    try:
        if len(fields) == 2:
            return f"{int(fields[0]) * 60 + float(fields[1]):.2f}"
        if len(fields) == 3:
            return f"{int(fields[0]) * 3600 + int(fields[1]) * 60 + float(fields[2]):.2f}"
    except ValueError:
        pass
    return "NA"


def command_value(path):
    if not path.exists():
        return ""
    with open(path) as handle:
        for line in handle:
            if "Command being timed:" in line:
                return line.split('"', 1)[-1].rsplit('"', 1)[0]
    return ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = Path(args.root)
    inputs = root / "inputs"
    results = root / "results"
    stats = root / "statistics"
    logs = root / "logs"
    stats.mkdir(parents=True, exist_ok=True)

    codes = ["Alyrata", "Bstricta", "Dstrictus", "Cviolacea"]
    atha_ids, atha_chrom = gff_ids(inputs / "Atha.gff")
    atha_set = set(atha_ids)
    atha_total = len(atha_set)

    summary_rows = []
    block_rows = []
    chrom_rows = []
    runtime_rows = []
    for code in codes:
        outgroup_ids, _ = gff_ids(inputs / f"{code}.gff")
        outgroup_total = len(set(outgroup_ids))
        prefix = results / f"Atha_{code}" / f"Atha_{code}"
        blast_hits = sum(1 for _ in open(prefix.with_suffix(".blast")))
        atha_collinear = set()
        outgroup_collinear = set()
        blocks = []
        current = None
        reported_collinear = None
        reported_all = None
        with open(prefix.with_suffix(".collinearity")) as handle:
            for line in handle:
                if line.startswith("# Number of collinear genes:"):
                    match = re.search(r"genes: (\d+)", line)
                    if match:
                        reported_collinear = int(match.group(1))
                elif line.startswith("# Number of all genes:"):
                    match = re.search(r"genes: (\d+)", line)
                    if match:
                        reported_all = int(match.group(1))
                block_match = BLOCK_RE.match(line)
                if block_match:
                    number, score, evalue, n_pairs, seq1, seq2, orientation = block_match.groups()
                    current = {
                        "pair": f"Atha_{code}", "block_id": int(number),
                        "score": float(score), "e_value": evalue,
                        "gene_pairs": int(n_pairs), "sequence_1": seq1,
                        "sequence_2": seq2, "orientation": orientation,
                    }
                    blocks.append(current)
                    continue
                pair_match = PAIR_RE.match(line)
                if pair_match:
                    gene1, gene2, _ = pair_match.groups()
                    if gene1 in atha_set:
                        atha_collinear.add(gene1)
                        outgroup_collinear.add(gene2)
                    elif gene2 in atha_set:
                        atha_collinear.add(gene2)
                        outgroup_collinear.add(gene1)

        sizes = [block["gene_pairs"] for block in blocks]
        orientations = Counter(block["orientation"] for block in blocks)
        summary_rows.append({
            "species": code,
            "blast_hits": blast_hits,
            "syntenic_blocks": len(blocks),
            "collinear_gene_pairs": sum(sizes),
            "unique_atha_collinear_genes": len(atha_collinear),
            "atha_total_genes": atha_total,
            "atha_collinear_coverage_pct": f"{100 * len(atha_collinear) / atha_total:.4f}",
            "unique_outgroup_collinear_genes": len(outgroup_collinear),
            "outgroup_total_genes": outgroup_total,
            "outgroup_collinear_coverage_pct": f"{100 * len(outgroup_collinear) / outgroup_total:.4f}",
            "mean_gene_pairs_per_block": f"{statistics.mean(sizes):.4f}" if sizes else "0",
            "median_gene_pairs_per_block": f"{statistics.median(sizes):.4f}" if sizes else "0",
            "max_gene_pairs_per_block": max(sizes, default=0),
            "plus_blocks": orientations["plus"],
            "minus_blocks": orientations["minus"],
            "mcscanx_reported_collinear_genes": reported_collinear if reported_collinear is not None else "NA",
            "mcscanx_reported_all_genes": reported_all if reported_all is not None else "NA",
        })
        block_rows.extend(blocks)
        per_chrom = Counter(atha_chrom[gene] for gene in atha_collinear)
        total_per_chrom = Counter(atha_chrom.values())
        for chrom in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
            chrom_rows.append({
                "species": code, "chromosome": chrom,
                "collinear_atha_genes": per_chrom[chrom],
                "total_atha_genes": total_per_chrom[chrom],
                "coverage_pct": f"{100 * per_chrom[chrom] / total_per_chrom[chrom]:.4f}",
            })
        blast_elapsed = time_value(logs / f"Atha_{code}.blast.time.log", "Elapsed (wall clock)")
        mcscan_elapsed = time_value(logs / f"Atha_{code}.mcscanx.time.log", "Elapsed (wall clock)")
        blast_command = command_value(logs / f"Atha_{code}.blast.time.log")
        thread_match = re.search(r"-num_threads (\d+)", blast_command)
        runtime_rows.append({
            "species": code,
            "blast_threads": thread_match.group(1) if thread_match else "NA",
            "blast_exit_status": time_value(logs / f"Atha_{code}.blast.time.log", "Exit status"),
            "blast_wall_seconds": elapsed_seconds(blast_elapsed),
            "blast_max_rss_kb": time_value(logs / f"Atha_{code}.blast.time.log", "Maximum resident set size"),
            "mcscanx_exit_status": time_value(logs / f"Atha_{code}.mcscanx.time.log", "Exit status"),
            "mcscanx_wall_seconds": elapsed_seconds(mcscan_elapsed),
            "mcscanx_max_rss_kb": time_value(logs / f"Atha_{code}.mcscanx.time.log", "Maximum resident set size"),
        })

    def write_table(path, rows):
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    write_table(stats / "synteny_summary.tsv", summary_rows)
    write_table(stats / "block_summary.tsv", block_rows)
    write_table(stats / "atha_chromosome_coverage.tsv", chrom_rows)
    write_table(stats / "runtime_summary.tsv", runtime_rows)
    print(f"Summarized {len(summary_rows)} species and {len(block_rows)} syntenic blocks")


if __name__ == "__main__":
    main()
