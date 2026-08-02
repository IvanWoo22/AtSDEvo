#!/usr/bin/env python3
"""Validate ID consistency and write compact MCScanX input QC."""

import argparse
from pathlib import Path


def fasta_ids(path):
    with open(path) as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">")]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = Path(args.root)
    inputs = root / "inputs"
    codes = ["Atha", "Alyrata", "Bstricta", "Dstrictus", "Cviolacea"]
    rows = []
    for code in codes:
        pep = fasta_ids(inputs / f"{code}.pep")
        with open(inputs / f"{code}.gff") as handle:
            gff = [line.rstrip().split("\t")[1] for line in handle if line.strip()]
        with open(inputs / f"{code}.bed") as handle:
            bed = [line.rstrip().split("\t")[3] for line in handle if line.strip()]
        problems = []
        if len(pep) != len(set(pep)): problems.append("duplicate_pep_id")
        if len(gff) != len(set(gff)): problems.append("duplicate_gff_id")
        if set(pep) != set(gff): problems.append("pep_gff_id_mismatch")
        if set(gff) != set(bed): problems.append("gff_bed_id_mismatch")
        rows.append((code, len(pep), len(gff), len(bed), "PASS" if not problems else ";".join(problems)))
    statistics = root / "statistics"
    statistics.mkdir(parents=True, exist_ok=True)
    with open(statistics / "input_qc.tsv", "w") as out:
        out.write("species\tpep_records\tgff_records\tbed_records\tstatus\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")
    failed = [row for row in rows if row[-1] != "PASS"]
    atha_count = rows[0][2]
    pair_rows = []
    for code, _, gff_count, _, _ in rows[1:]:
        pair_path = root / "pair_inputs" / f"Atha_{code}.gff"
        with open(pair_path) as handle:
            pair_ids = [line.rstrip().split("\t")[1] for line in handle if line.strip()]
        expected = atha_count + gff_count
        problems = []
        if len(pair_ids) != expected: problems.append("wrong_record_count")
        if len(pair_ids) != len(set(pair_ids)): problems.append("cross_species_duplicate_id")
        pair_rows.append((f"Atha_{code}", len(pair_ids), expected, "PASS" if not problems else ";".join(problems)))
    with open(statistics / "pair_input_qc.tsv", "w") as out:
        out.write("pair\trecords\texpected_records\tstatus\n")
        for row in pair_rows:
            out.write("\t".join(map(str, row)) + "\n")
    failed.extend(row for row in pair_rows if row[-1] != "PASS")
    if failed:
        raise SystemExit(f"Validation failed: {failed}")
    print(f"Validated {len(rows)} species and {len(pair_rows)} pairs; all IDs and counts match")


if __name__ == "__main__":
    main()
