#!/usr/bin/env python3
"""Create canonical TAIR12 peptide, MCScanX coordinate and BED files."""

import argparse
import re
import subprocess
from collections import defaultdict
from pathlib import Path


def attrs(text):
    result = {}
    for field in text.rstrip().split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
    return result


def read_fasta(path):
    records = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).replace("*", "")
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            records[name] = "".join(chunks).replace("*", "")
    return records


def isoform_key(name):
    match = re.search(r"\.(\d+)$", name)
    return (int(match.group(1)) if match else 10**9, name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--gffread", default="gffread")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    all_pep_path = outdir / ".Atha.all_transcripts.pep"
    subprocess.run(
        [args.gffread, args.gff, "-g", args.genome, "-y", str(all_pep_path)],
        check=True,
    )
    all_pep = read_fasta(all_pep_path)

    genes = {}
    transcripts = defaultdict(list)
    with open(args.gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) != 9:
                continue
            seqid, _, feature, start, end, _, strand, _, raw_attrs = parts
            data = attrs(raw_attrs)
            if feature == "gene" and data.get("gene_biotype") == "protein_coding":
                internal_id = data["ID"]
                locus = data.get("locus_tag", internal_id)
                gene_id = re.sub(r"^(?:gene-)?TAIR12_TAIR12_", "", locus)
                genes[internal_id] = (gene_id, seqid, int(start), int(end), strand)
            elif feature == "mRNA":
                parent = data.get("Parent")
                transcript_id = data.get("ID")
                standard_name = data.get("standard_name")
                if parent and transcript_id and standard_name:
                    transcripts[parent].append((standard_name, transcript_id))

    selected = []
    missing = []
    for internal_id, (gene_id, seqid, start, end, strand) in genes.items():
        candidates = sorted(transcripts.get(internal_id, []), key=lambda x: isoform_key(x[0]))
        candidates = [item for item in candidates if item[1] in all_pep]
        if not candidates:
            missing.append(gene_id)
            continue
        standard_name, transcript_id = candidates[0]
        selected.append((seqid, start, end, strand, gene_id, standard_name, transcript_id))

    selected.sort(key=lambda x: (int(x[0][3:]) if x[0].startswith("Chr") and x[0][3:].isdigit() else 999, x[1], x[4]))
    if missing:
        raise SystemExit(f"No translated mRNA for {len(missing)} protein-coding genes; first: {missing[:10]}")

    with open(outdir / "Atha.pep", "w") as pep_out, \
         open(outdir / "Atha.gff", "w") as gff_out, \
         open(outdir / "Atha.bed", "w") as bed_out, \
         open(outdir / "Atha.transcript_map.tsv", "w") as map_out:
        map_out.write("gene_id\tselected_transcript\tgff_transcript_id\n")
        for seqid, start, end, strand, gene_id, standard_name, transcript_id in selected:
            sequence = all_pep[transcript_id]
            pep_out.write(f">{gene_id}\n")
            for i in range(0, len(sequence), 60):
                pep_out.write(sequence[i:i + 60] + "\n")
            left, right = (start, end) if strand == "+" else (end, start)
            gff_out.write(f"Atha_{seqid}\t{gene_id}\t{left}\t{right}\n")
            bed_out.write(f"{seqid}\t{start - 1}\t{end}\t{gene_id}\t0\t{strand}\n")
            map_out.write(f"{gene_id}\t{standard_name}\t{transcript_id}\n")

    all_pep_path.unlink()
    print(f"TAIR12 canonical protein-coding genes: {len(selected)}")


if __name__ == "__main__":
    main()
