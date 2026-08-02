#!/usr/bin/env python3
"""Prepare contiguous P/D/outgroup blocks for conservative GENECONV screening."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def fasta_record(name: str, sequence: str) -> str:
    return f">{name}\n" + "\n".join(
        sequence[index : index + 80]
        for index in range(0, len(sequence), 80)
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", required=True, type=Path)
    parser.add_argument("--minimum-block-bp", type=int, default=100)
    parser.add_argument("--minimum-polymorphic-sites", type=int, default=5)
    args = parser.parse_args()

    eligible = {
        row["event_id"]
        for row in read_tsv(
            args.pilot
            / "sequence_variation/primary_endpoint_eligibility.tsv"
        )
        if row["primary_endpoint_eligible"] == "PASS"
    }
    sites: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(
        args.pilot / "sequence_variation/polarized_sites.tsv"
    ):
        if row["event_id"] in eligible:
            sites[row["event_id"]].append(row)

    input_dir = args.pilot / "sequence_variation/geneconv/input"
    input_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for event_id, event_sites in sorted(sites.items()):
        event_sites.sort(key=lambda row: int(row["core_position_1based"]))
        blocks: list[list[dict[str, str]]] = []
        current = []
        previous = None
        for row in event_sites:
            position = int(row["core_position_1based"])
            if current and position != previous + 1:
                blocks.append(current)
                current = []
            current.append(row)
            previous = position
        if current:
            blocks.append(current)
        retained_index = 0
        for block in blocks:
            p_sequence = "".join(row["P_base"] for row in block)
            d_sequence = "".join(row["D_base"] for row in block)
            out_sequence = "".join(
                row["boundary_bidirectional_ancestral_base"]
                for row in block
            )
            polymorphic = sum(
                len({p, d, out}) > 1
                for p, d, out in zip(p_sequence, d_sequence, out_sequence)
            )
            if (
                len(block) < args.minimum_block_bp
                or polymorphic < args.minimum_polymorphic_sites
            ):
                continue
            retained_index += 1
            filename = f"{event_id}.block{retained_index:03d}.fa"
            path = input_dir / filename
            path.write_text(
                fasta_record("P", p_sequence)
                + fasta_record("D", d_sequence)
                + fasta_record("OUT", out_sequence)
            )
            manifest.append(
                {
                    "event_id": event_id,
                    "block_id": f"{event_id}.block{retained_index:03d}",
                    "core_start_1based": block[0]["core_position_1based"],
                    "core_end_1based": block[-1]["core_position_1based"],
                    "block_bp": len(block),
                    "polymorphic_sites": polymorphic,
                    "input_fasta": str(path),
                    "screen_interpretation": (
                        "conversion_risk_flag_not_standalone_proof"
                    ),
                }
            )
    write_tsv(
        args.pilot / "sequence_variation/geneconv/block_manifest.tsv",
        manifest,
    )
    print(
        f"Prepared {len(manifest)} contiguous blocks from "
        f"{len({row['event_id'] for row in manifest})} events"
    )


if __name__ == "__main__":
    main()
