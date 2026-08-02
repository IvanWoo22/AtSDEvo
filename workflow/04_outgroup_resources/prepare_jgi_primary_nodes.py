#!/usr/bin/env python3
"""Validate and normalize the four version-matched JGI primary-node datasets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from collections import Counter
from pathlib import Path


PRIMARY_ORDER = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_attrs(text: str) -> dict[str, str]:
    result = {}
    for field in text.rstrip().split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
    return result


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name = None
    length = 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                if name in lengths:
                    raise ValueError(f"Duplicate FASTA ID in {path}: {name}")
                length = 0
            else:
                length += len(line.strip())
        if name is not None:
            lengths[name] = length
    return lengths


def read_primary_peptides(path: Path) -> dict[str, str]:
    peptides: dict[str, str] = {}
    transcript = None
    chunks: list[str] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if transcript is not None:
                    peptides[transcript] = "".join(chunks).replace("*", "")
                match = re.search(r"(?:^| )transcript=([^ ]+)", line.rstrip())
                if not match:
                    raise ValueError(f"Missing transcript= in peptide header: {line.rstrip()}")
                transcript = match.group(1)
                if transcript in peptides:
                    raise ValueError(f"Duplicate peptide transcript: {transcript}")
                chunks = []
            else:
                chunks.append(line.strip())
        if transcript is not None:
            peptides[transcript] = "".join(chunks).replace("*", "")
    return peptides


def read_primary_mrnas(path: Path) -> dict[str, tuple[str, int, int, str]]:
    mrnas: dict[str, tuple[str, int, int, str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "mRNA":
                continue
            attrs = parse_attrs(fields[8])
            if attrs.get("longest") != "1":
                continue
            transcript = attrs.get("Name")
            if not transcript:
                raise ValueError(f"Primary mRNA without Name in {path}: {line.rstrip()}")
            if transcript in mrnas:
                raise ValueError(f"Duplicate primary mRNA Name in {path}: {transcript}")
            mrnas[transcript] = (
                fields[0],
                int(fields[3]),
                int(fields[4]),
                fields[6],
            )
    return mrnas


def find_single(root: Path, pattern: str) -> Path:
    matches = sorted(root.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"Expected one match for {pattern}, found {len(matches)}: {matches}")
    return matches[0]


def write_wrapped(handle, sequence: str, width: int = 60) -> None:
    for index in range(0, len(sequence), width):
        handle.write(sequence[index:index + width] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    manifest_path = args.raw_root / "Download_2214609_File_Manifest.csv"
    with manifest_path.open(newline="") as handle:
        manifest = list(csv.DictReader(handle))
    by_name = {row["filename"]: row for row in manifest}
    if len(by_name) != 12:
        raise SystemExit(f"Expected 12 biological files in manifest, found {len(by_name)}")

    args.output.mkdir(parents=True, exist_ok=True)
    qc_rows = []
    resource_rows = []

    for code in PRIMARY_ORDER:
        species_root = find_single(
            args.raw_root / "Phytozome", f"**/{code}"
        )
        genome = find_single(species_root, "**/assembly/*.fa.gz")
        gff = find_single(species_root, "**/annotation/*.gene_exons.gff3.gz")
        peptide = find_single(
            species_root, "**/annotation/*.protein_primaryTranscriptOnly.fa.gz"
        )

        for source in (genome, gff, peptide):
            expected = by_name[source.name]["md5 checksum"]
            observed = md5sum(source)
            if observed != expected:
                raise SystemExit(
                    f"MD5 mismatch for {source}: expected {expected}, observed {observed}"
                )

        genome_lengths = read_fasta_lengths(genome)
        peptides = read_primary_peptides(peptide)
        mrnas = read_primary_mrnas(gff)
        peptide_ids = set(peptides)
        mrna_ids = set(mrnas)
        missing_peptide = sorted(mrna_ids - peptide_ids)
        missing_mrna = sorted(peptide_ids - mrna_ids)
        missing_seqids = sorted({value[0] for value in mrnas.values()} - set(genome_lengths))
        out_of_bounds = [
            transcript
            for transcript, (seqid, start, end, _) in mrnas.items()
            if seqid in genome_lengths and (start < 1 or end > genome_lengths[seqid])
        ]
        if missing_peptide or missing_mrna or missing_seqids or out_of_bounds:
            raise SystemExit(
                f"{code} validation failed: missing_peptide={missing_peptide[:5]}, "
                f"missing_mrna={missing_mrna[:5]}, missing_seqids={missing_seqids[:5]}, "
                f"out_of_bounds={out_of_bounds[:5]}"
            )

        records = sorted(
            mrnas.items(),
            key=lambda item: (item[1][0], item[1][1], item[0]),
        )
        with (
            (args.output / f"{code}.pep").open("w") as pep_out,
            (args.output / f"{code}.gff").open("w") as gff_out,
            (args.output / f"{code}.bed").open("w") as bed_out,
            (args.output / f"{code}.transcript_map.tsv").open("w") as map_out,
        ):
            map_out.write(
                "mcscan_id\tjgi_transcript_id\tassembly_sequence\tstart\tend\tstrand\n"
            )
            for transcript, (seqid, start, end, strand) in records:
                mcscan_id = f"{code}_{transcript}"
                pep_out.write(f">{mcscan_id}\n")
                write_wrapped(pep_out, peptides[transcript])
                left, right = (start, end) if strand == "+" else (end, start)
                gff_out.write(
                    f"{code}_{seqid}\t{mcscan_id}\t{left}\t{right}\n"
                )
                bed_out.write(
                    f"{seqid}\t{start - 1}\t{end}\t{mcscan_id}\t0\t{strand}\n"
                )
                map_out.write(
                    f"{mcscan_id}\t{transcript}\t{seqid}\t{start}\t{end}\t{strand}\n"
                )

        lengths = [len(sequence) for sequence in peptides.values()]
        qc_rows.append(
            {
                "code": code,
                "assembly_sequences": len(genome_lengths),
                "assembly_bp": sum(genome_lengths.values()),
                "primary_mrnas": len(mrnas),
                "primary_peptides": len(peptides),
                "median_peptide_aa": sorted(lengths)[len(lengths) // 2],
                "minus_strand_genes": Counter(value[3] for value in mrnas.values())["-"],
                "gff_seqids_missing_from_assembly": len(missing_seqids),
                "out_of_bounds_mrnas": len(out_of_bounds),
                "pep_gff_id_match": "PASS",
                "md5_match": "PASS",
            }
        )
        for kind, source in (("assembly", genome), ("gff3", gff), ("peptide", peptide)):
            resource_rows.append(
                {
                    "code": code,
                    "resource": kind,
                    "filename": source.name,
                    "relative_path": str(source.relative_to(args.raw_root)),
                    "md5": md5sum(source),
                    "jgi_grouping_id": by_name[source.name]["JGI Grouping ID"],
                    "dataset_name": by_name[source.name]["Genome/Metagenome Name"],
                }
            )

    def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)

    write_tsv(args.output / "input_qc.tsv", qc_rows)
    write_tsv(args.output / "resource_manifest.tsv", resource_rows)
    print(
        "Prepared and validated "
        f"{len(PRIMARY_ORDER)} JGI primary nodes; "
        f"{sum(row['primary_mrnas'] for row in qc_rows)} representative proteins"
    )


if __name__ == "__main__":
    main()
