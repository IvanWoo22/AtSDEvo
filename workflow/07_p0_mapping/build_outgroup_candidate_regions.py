#!/usr/bin/env python3
"""Build anchor-centred outgroup candidate regions for high-priority events."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
BLOCK_RE = re.compile(r"^## Alignment (\d+):")
PAIR_RE = re.compile(r"^\s*\d+-\s*\d+:\s+(\S+)\s+(\S+)\s+\S+")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        if not rows:
            return
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_coordinates(path: Path) -> dict[str, tuple[str, int, int]]:
    coordinates = {}
    with path.open() as handle:
        for line in handle:
            seqid, gene, first, second = line.rstrip().split("\t")
            coordinates[gene] = (
                seqid,
                min(int(first), int(second)) - 1,
                max(int(first), int(second)),
            )
    return coordinates


def read_blocks(
    path: Path,
    atha: dict[str, tuple[str, int, int]],
    outgroup: dict[str, tuple[str, int, int]],
) -> dict[int, list[tuple[str, str]]]:
    blocks: dict[int, list[tuple[str, str]]] = defaultdict(list)
    block_id = None
    with path.open() as handle:
        for line in handle:
            match = BLOCK_RE.match(line)
            if match:
                block_id = int(match.group(1))
                continue
            pair = PAIR_RE.match(line)
            if not pair or block_id is None:
                continue
            left, right = pair.groups()
            if left in atha and right in outgroup:
                blocks[block_id].append((left, right))
            elif right in atha and left in outgroup:
                blocks[block_id].append((right, left))
    return blocks


def read_gzip_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    name = ""
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            else:
                sequences[name].append(line.strip())
    return {name: "".join(parts) for name, parts in sequences.items()}


def write_fasta_record(handle, name: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for index in range(0, len(sequence), width):
        handle.write(sequence[index : index + width] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--pilot", required=True, type=Path)
    parser.add_argument("--radius", type=int, default=100_000)
    args = parser.parse_args()

    mcscan = args.project / "05_mcscanx_synteny"
    resources = (
        args.project
        / "04_outgroup_resources/primary_nodes_jgi_20260724/prepared"
    )
    raw_root = (
        args.project / "04_outgroup_resources/primary_nodes_jgi_20260724/raw"
    )
    events = read_tsv(args.pilot / "inputs/high_priority_core_eligible.tsv")
    atha = read_coordinates(mcscan / "inputs/Atha.gff")
    resource_manifest = read_tsv(resources / "resource_manifest.tsv")
    genome_path = {
        row["code"]: raw_root / row["relative_path"]
        for row in resource_manifest
        if row["resource"] == "assembly"
    }

    manifest_rows: list[dict[str, object]] = []
    fasta_dir = args.pilot / "outgroup_mapping/candidate_regions"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    for code in SPECIES:
        outgroup = read_coordinates(mcscan / f"inputs/{code}.gff")
        transcript_map = {
            row["mcscan_id"]: row
            for row in read_tsv(resources / f"{code}.transcript_map.tsv")
        }
        blocks = read_blocks(
            mcscan / f"results/Atha_{code}/Atha_{code}.collinearity",
            atha,
            outgroup,
        )
        genome = read_gzip_fasta(genome_path[code])
        fasta_path = fasta_dir / f"{code}.candidate_regions.fa"
        with fasta_path.open("w") as fasta_out:
            for event in events:
                state = int(event[f"{code}_state"])
                supported_copies = (
                    (1, 2) if state == 3 else (1,) if state == 1 else (2,) if state == 2 else ()
                )
                for copy in supported_copies:
                    block_text = event[f"{code}_copy{copy}_block_ids"]
                    if not block_text:
                        continue
                    event_chrom = event[f"representative_copy{copy}_chrom"]
                    event_start = int(event[f"representative_copy{copy}_start"])
                    event_end = int(event[f"representative_copy{copy}_end"])
                    event_midpoint = (event_start + event_end) / 2
                    for block_id in sorted({int(value) for value in block_text.split(",")}):
                        anchors = []
                        for atha_gene, outgroup_gene in blocks.get(block_id, []):
                            atha_seqid, atha_start, atha_end = atha[atha_gene]
                            if atha_seqid.removeprefix("Atha_") != event_chrom:
                                continue
                            out_seqid_prefixed, out_start, out_end = outgroup[outgroup_gene]
                            source = transcript_map[outgroup_gene]
                            source_seqid = source["assembly_sequence"]
                            if out_seqid_prefixed != f"{code}_{source_seqid}":
                                raise AssertionError(
                                    f"Sequence-ID translation mismatch: {out_seqid_prefixed}"
                                )
                            distance = abs((atha_start + atha_end) / 2 - event_midpoint)
                            anchors.append(
                                (
                                    distance,
                                    atha_gene,
                                    outgroup_gene,
                                    source_seqid,
                                    atha_start,
                                    atha_end,
                                    out_start,
                                    out_end,
                                )
                            )
                        if not anchors:
                            continue
                        (
                            distance,
                            atha_gene,
                            outgroup_gene,
                            source_seqid,
                            _,
                            _,
                            out_start,
                            out_end,
                        ) = min(anchors)
                        sequence = genome[source_seqid]
                        left_anchors = [
                            anchor
                            for anchor in anchors
                            if anchor[5] <= event_start
                        ]
                        right_anchors = [
                            anchor
                            for anchor in anchors
                            if anchor[4] >= event_end
                        ]
                        left_anchor = (
                            min(
                                left_anchors,
                                key=lambda anchor: event_start - anchor[5],
                            )
                            if left_anchors
                            else None
                        )
                        right_anchor = (
                            min(
                                right_anchors,
                                key=lambda anchor: anchor[4] - event_end,
                            )
                            if right_anchors
                            else None
                        )
                        two_sided = (
                            left_anchor is not None
                            and right_anchor is not None
                            and left_anchor[3] == right_anchor[3]
                        )
                        if two_sided:
                            start = max(
                                0,
                                min(left_anchor[6], right_anchor[6]) - 20_000,
                            )
                            end = min(
                                len(sequence),
                                max(left_anchor[7], right_anchor[7]) + 20_000,
                            )
                            outgroup_orientation = (
                                "plus"
                                if (left_anchor[6] + left_anchor[7])
                                < (right_anchor[6] + right_anchor[7])
                                else "minus"
                            )
                        else:
                            center = (out_start + out_end) // 2
                            start = max(0, center - args.radius)
                            end = min(len(sequence), center + args.radius)
                            outgroup_orientation = "unresolved"
                        candidate_id = (
                            f"{event['event_id']}|{code}|copy{copy}|"
                            f"block{block_id}|{source_seqid}:{start}-{end}"
                        )
                        write_fasta_record(
                            fasta_out, candidate_id, sequence[start:end]
                        )
                        manifest_rows.append(
                            {
                                "candidate_id": candidate_id,
                                "event_id": event["event_id"],
                                "age_bin": event["strict_age_bin"],
                                "species": code,
                                "species_state": state,
                                "TAIR12_copy": f"copy{copy}",
                                "TAIR12_copy_role": (
                                    "P"
                                    if event["provisional_p_copy"] == f"copy{copy}"
                                    else "D"
                                ),
                                "block_id": block_id,
                                "nearest_atha_anchor_gene": atha_gene,
                                "nearest_outgroup_anchor_gene": outgroup_gene,
                                "atha_anchor_distance_to_event_midpoint_bp": (
                                    f"{distance:.1f}"
                                ),
                                "two_sided_ordered_anchor_status": (
                                    "PASS" if two_sided else "FAIL"
                                ),
                                "left_atha_anchor_gene": (
                                    left_anchor[1] if left_anchor else "NA"
                                ),
                                "right_atha_anchor_gene": (
                                    right_anchor[1] if right_anchor else "NA"
                                ),
                                "left_outgroup_anchor_gene": (
                                    left_anchor[2] if left_anchor else "NA"
                                ),
                                "right_outgroup_anchor_gene": (
                                    right_anchor[2] if right_anchor else "NA"
                                ),
                                "outgroup_anchor_orientation": outgroup_orientation,
                                "outgroup_scaffold": source_seqid,
                                "candidate_start": start,
                                "candidate_end": end,
                                "candidate_bp": end - start,
                                "candidate_radius_bp": args.radius,
                            }
                        )
        del genome

    write_tsv(
        args.pilot / "outgroup_mapping/candidate_region_manifest.tsv",
        manifest_rows,
    )
    print(
        f"Built {len(manifest_rows)} anchor-centred candidate regions for "
        f"{len(events)} core-eligible events"
    )


if __name__ == "__main__":
    main()
