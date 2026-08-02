#!/usr/bin/env python3
"""Prepare targeted non-state3 SD locus queries for outgroup genome mapping."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            else:
                sequences[name].append(line.strip())
    return {name: "".join(parts) for name, parts in sequences.items()}


def fasta_record(name: str, sequence: str, width: int = 80) -> str:
    return f">{name}\n" + "\n".join(
        sequence[index : index + width]
        for index in range(0, len(sequence), width)
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--state3-controls", type=int, default=50)
    parser.add_argument("--minimum-uppercase-bp", type=int, default=200)
    args = parser.parse_args()

    analysis = (
        args.project
        / "06_sd_age_tracing_preparation/event_first_reanalysis"
    )
    events = read(analysis / "events/event_first_events.tsv")
    genome = read_fasta(
        args.project
        / "01_reference/prepared_data/"
        "TAIR12.Col-CC.annotation_softmasked.fa"
    )

    candidate_rows = []
    query_rows = []
    fasta_parts: dict[str, list[str]] = defaultdict(list)
    selected_controls: dict[str, set[str]] = {}
    for species in SPECIES:
        controls = [
            row["event_id"]
            for row in events
            if row[f"{species}_state"] == "3"
            and row["strict_minimum_1kb_core"] == "PASS"
        ][: args.state3_controls]
        selected_controls[species] = set(controls)

    for event in events:
        event_id = event["event_id"]
        sequences = {}
        for locus in ("locus_A", "locus_B"):
            chrom = event[f"{locus}_chrom"]
            start = int(event[f"{locus}_representative_start"])
            end = int(event[f"{locus}_representative_end"])
            sequence = genome[chrom][start:end]
            uppercase_bp = sum(base in "ACGT" for base in sequence)
            # Minimap2 does not honor soft-mask case. Replace lowercase bases
            # with N so repeat-masked sequence cannot drive a rescue.
            uppercase_only = "".join(
                base if base in "ACGT" else "N" for base in sequence
            )
            sequences[locus] = (
                chrom,
                start,
                end,
                sequence,
                uppercase_only,
                uppercase_bp,
            )
        for species in SPECIES:
            state = event[f"{species}_state"]
            role = (
                "missing_state_target"
                if state != "3"
                else "state3_positive_control"
                if event_id in selected_controls[species]
                else "skip"
            )
            if role == "skip":
                continue
            candidate_rows.append(
                {
                    "event_id": event_id,
                    "species": species,
                    "mcscan_state": state,
                    "candidate_role": role,
                    "locus_A_mcscan_detected": state in ("1", "3"),
                    "locus_B_mcscan_detected": state in ("2", "3"),
                    "locus_A_block_ids": event[
                        f"{species}_locus_A_block_ids"
                    ],
                    "locus_B_block_ids": event[
                        f"{species}_locus_B_block_ids"
                    ],
                }
            )
            for locus in ("locus_A", "locus_B"):
                chrom, start, end, raw, masked, uppercase_bp = sequences[locus]
                query_id = f"{event_id}__{locus}"
                eligible = uppercase_bp >= args.minimum_uppercase_bp
                query_rows.append(
                    {
                        "query_id": query_id,
                        "event_id": event_id,
                        "species": species,
                        "event_locus": locus,
                        "candidate_role": role,
                        "mcscan_state": state,
                        "mcscan_detected": (
                            state in ("1", "3")
                            if locus == "locus_A"
                            else state in ("2", "3")
                        ),
                        "TAIR12_chrom": chrom,
                        "TAIR12_start": start,
                        "TAIR12_end": end,
                        "query_interval_bp": end - start,
                        "query_uppercase_ACGT_bp": uppercase_bp,
                        "query_uppercase_fraction": (
                            f"{uppercase_bp / max(1, end - start):.6f}"
                        ),
                        "mapping_query_eligible": "PASS" if eligible else "FAIL",
                    }
                )
                if eligible:
                    fasta_parts[species].append(
                        fasta_record(query_id, masked)
                    )

    args.output.mkdir(parents=True, exist_ok=True)
    write(args.output / "candidate_event_species.tsv", candidate_rows)
    write(args.output / "query_manifest.tsv", query_rows)
    for species in SPECIES:
        (args.output / f"queries/{species}.non3_plus_controls.fa").parent.mkdir(
            parents=True, exist_ok=True
        )
        (args.output / f"queries/{species}.non3_plus_controls.fa").write_text(
            "".join(fasta_parts[species])
        )
    summary = []
    for species in SPECIES:
        species_candidates = [
            row for row in candidate_rows if row["species"] == species
        ]
        species_queries = [
            row
            for row in query_rows
            if row["species"] == species
            and row["mapping_query_eligible"] == "PASS"
        ]
        summary.append(
            {
                "species": species,
                "non3_event_species": sum(
                    row["candidate_role"] == "missing_state_target"
                    for row in species_candidates
                ),
                "state3_positive_controls": sum(
                    row["candidate_role"] == "state3_positive_control"
                    for row in species_candidates
                ),
                "eligible_locus_queries": len(species_queries),
            }
        )
    write(args.output / "preparation_summary.tsv", summary)
    print(
        f"Prepared {len(candidate_rows)} event-species candidates and "
        f"{sum(row['mapping_query_eligible'] == 'PASS' for row in query_rows)} "
        "locus queries"
    )


if __name__ == "__main__":
    main()
