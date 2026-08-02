#!/usr/bin/env python3
"""Normalize TAIR12 RepeatMasker accessions and export its repeat BED."""

from __future__ import annotations

import argparse
from pathlib import Path


ACCESSION = {f"OZ40868{i}.1": f"Chr{i - 2}" for i in range(3, 8)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--masked-fasta", required=True, type=Path)
    parser.add_argument("--repeat-out", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    fasta_out = args.output_root / "TAIR12.Col-CC.repeatmasker_softmasked.fa"
    with args.masked_fasta.open() as source, fasta_out.open("w") as target:
        for line in source:
            if line.startswith(">"):
                accession = line[1:].split()[0]
                if accession not in ACCESSION:
                    raise SystemExit(f"Unexpected sequence: {accession}")
                target.write(f">{ACCESSION[accession]}\n")
            else:
                target.write(line)

    intervals: dict[str, list[tuple[int, int]]] = {chrom: [] for chrom in ACCESSION.values()}
    with args.repeat_out.open() as source:
        for line in source:
            fields = line.split()
            if len(fields) < 15 or fields[4] not in ACCESSION:
                continue
            intervals[ACCESSION[fields[4]]].append((int(fields[5]) - 1, int(fields[6])))
    bed_out = args.output_root / "TAIR12.repeatmasker.merged.bed"
    with bed_out.open("w") as target:
        for chrom in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
            values = sorted(intervals[chrom])
            if not values:
                continue
            left, right = values[0]
            for start, end in values[1:]:
                if start <= right:
                    right = max(right, end)
                else:
                    target.write(f"{chrom}\t{left}\t{right}\n")
                    left, right = start, end
            target.write(f"{chrom}\t{left}\t{right}\n")
    print(f"reference={fasta_out}")
    print(f"repeat_bed={bed_out}")


if __name__ == "__main__":
    main()
