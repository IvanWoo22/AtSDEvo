#!/usr/bin/env python3
"""Shared, dependency-light utilities for the formal micro-indel workflow."""

from __future__ import annotations

import csv
import math
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f"duplicate FASTA identifier in {path}: {name}")
                records[name] = []
            elif not name:
                raise ValueError(f"sequence before FASTA header in {path}")
            else:
                records[name].append(line)
    return {key: "".join(value) for key, value in records.items()}


def fasta_text(records: dict[str, str], width: int = 80) -> str:
    chunks = []
    for name, sequence in records.items():
        chunks.append(f">{name}")
        chunks.extend(
            sequence[start : start + width]
            for start in range(0, len(sequence), width)
        )
    return "\n".join(chunks) + "\n"


def interval_overlap(
    intervals: dict[str, list[tuple[int, int]]],
    chrom: str,
    start: int,
    end: int,
) -> bool:
    return any(left < end and start < right for left, right in intervals.get(chrom, []))


def longest_homopolymer(sequence: str) -> int:
    longest = current = 0
    previous = ""
    for base in sequence.upper():
        if base not in "ACGT":
            previous, current = "", 0
        elif base == previous:
            current += 1
        else:
            previous, current = base, 1
        longest = max(longest, current)
    return longest


def sequence_entropy(sequence: str) -> float:
    counts = [sequence.upper().count(base) for base in "ACGT"]
    total = sum(counts)
    if not total:
        return 0.0
    return -sum(
        (count / total) * math.log2(count / total) for count in counts if count
    )


def exact_two_sided_binomial(successes: int, trials: int) -> float:
    """Exact two-sided binomial P value for p=0.5."""
    if trials == 0:
        return math.nan
    observed = math.comb(trials, successes) / (2**trials)
    return min(
        1.0,
        sum(
            math.comb(trials, value) / (2**trials)
            for value in range(trials + 1)
            if math.comb(trials, value) / (2**trials) <= observed + 1e-15
        ),
    )


def parse_hsp(fields: list[str]) -> dict[str, object]:
    """Parse the 16-column BLASTN format used by the P0 projection workflow."""
    if len(fields) < 16:
        raise ValueError(f"expected 16 BLASTN fields, observed {len(fields)}")
    return {
        "query": fields[0],
        "target": fields[1],
        "qlen": int(fields[2]),
        "length": int(fields[3]),
        "pident": float(fields[4]),
        "qstart": int(fields[7]),
        "qend": int(fields[8]),
        "sstart": int(fields[9]),
        "send": int(fields[10]),
        "evalue": float(fields[11]),
        "bitscore": float(fields[12]),
        # NCBI BLAST+ omits the unsupported qstrand token; query coordinates
        # are in forward orientation and field 14 is the subject strand.
        "qstrand": "plus",
        "sstrand": fields[13],
        "qseq": fields[14],
        "sseq": fields[15],
    }


def _interval_coverage(intervals: list[tuple[int, int]]) -> int:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def _effective_subject_strand(hsp: dict[str, object]) -> str:
    strand = str(hsp["sstrand"])
    if hsp["qstrand"] == "plus":
        return strand
    return "minus" if strand == "plus" else "plus"


def cluster_hsps(
    hsps: list[dict[str, object]], max_subject_gap: int = 5000
) -> list[list[dict[str, object]]]:
    """Cluster HSPs by effective strand and proximity on one target region."""
    clusters = []
    for strand in ("plus", "minus"):
        subset = [
            hsp for hsp in hsps if _effective_subject_strand(hsp) == strand
        ]
        subset.sort(key=lambda hsp: min(int(hsp["sstart"]), int(hsp["send"])))
        current: list[dict[str, object]] = []
        current_end = -1
        for hsp in subset:
            start = min(int(hsp["sstart"]), int(hsp["send"]))
            end = max(int(hsp["sstart"]), int(hsp["send"]))
            if current and start > current_end + max_subject_gap:
                clusters.append(current)
                current = []
            current.append(hsp)
            current_end = max(current_end, end)
        if current:
            clusters.append(current)
    return clusters


def cluster_score(cluster: list[dict[str, object]]) -> tuple[int, float]:
    coverage = _interval_coverage(
        [
            (
                min(int(hsp["qstart"]), int(hsp["qend"])) - 1,
                max(int(hsp["qstart"]), int(hsp["qend"])),
            )
            for hsp in cluster
        ]
    )
    return coverage, sum(float(hsp["bitscore"]) for hsp in cluster)


def project_cluster(
    cluster: list[dict[str, object]], query_sequence: str
) -> tuple[dict[int, tuple[str, int]], int]:
    """Project ungapped query positions to subject bases and 0-based coordinates."""
    projected: dict[int, tuple[float, str, int]] = {}
    query_validation_failures = 0
    for hsp in sorted(
        cluster, key=lambda row: float(row["bitscore"]), reverse=True
    ):
        if hsp["qstrand"] == "plus":
            position = int(hsp["qstart"]) - 1
            step = 1
            normalize = lambda base: base
        else:
            position = int(hsp["qstart"]) - 1
            step = -1
            normalize = lambda base: base.translate(COMPLEMENT)
        subject_position = int(hsp["sstart"]) - 1
        subject_step = 1 if int(hsp["send"]) >= int(hsp["sstart"]) else -1
        for query_base, subject_base in zip(str(hsp["qseq"]), str(hsp["sseq"])):
            if query_base == "-":
                if subject_base != "-":
                    subject_position += subject_step
                continue
            normalized_query = normalize(query_base).upper()
            if (
                0 <= position < len(query_sequence)
                and normalized_query in "ACGT"
                and query_sequence[position].upper() in "ACGT"
                and normalized_query != query_sequence[position].upper()
            ):
                query_validation_failures += 1
            if subject_base != "-":
                base = normalize(subject_base).upper()
                if base in "ACGT" and position not in projected:
                    projected[position] = (
                        float(hsp["bitscore"]),
                        base,
                        subject_position,
                    )
                subject_position += subject_step
            position += step
    return (
        {
            position: (value[1], value[2])
            for position, value in projected.items()
        },
        query_validation_failures,
    )
