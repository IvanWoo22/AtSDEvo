#!/usr/bin/env python3
"""De novo 1-10 bp indel discovery from continuous syntenic haplotypes.

BISER I/D operations are not read. BISER supplies only the already accepted SD
locus and paired homologous block coordinates. Continuous sequence between
core-coordinate anchors is realigned independently with MAFFT, MUSCLE5, and
fixed-species-tree PRANK; gap blocks are discovered anew from each MSA.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from microindel_utils import (
    cluster_hsps,
    cluster_score,
    exact_two_sided_binomial,
    fasta_text,
    interval_overlap,
    longest_homopolymer,
    parse_hsp,
    project_cluster,
    read_fasta,
    read_tsv,
    sequence_entropy,
    write_tsv,
)


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def block_ranges(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    result = []
    cursor = 0
    for row in sorted(rows, key=lambda value: int(value["core_block_index"])):
        length = int(row["core_block_bp"])
        normalized: dict[str, object] = {
            **row,
            "core_start": cursor,
            "core_end": cursor + length,
        }
        for copy in ("copy1", "copy2"):
            role = row[f"{copy}_role"]
            for field in ("chrom", "start", "end", "strand"):
                normalized[f"{role.lower()}_{field}"] = row[f"{copy}_{field}"]
        result.append(normalized)
        cursor += length
    return result


def core_cut(
    blocks: list[dict[str, object]], role: str, position: int, side: str
) -> tuple[str, int, str]:
    total = int(blocks[-1]["core_end"])
    selected = None
    if side == "start":
        if position == total:
            selected = blocks[-1]
        else:
            selected = next(
                row
                for row in blocks
                if int(row["core_start"]) <= position < int(row["core_end"])
            )
    else:
        if position == 0:
            selected = blocks[0]
        else:
            selected = next(
                row
                for row in blocks
                if int(row["core_start"]) < position <= int(row["core_end"])
            )
    offset = position - int(selected["core_start"])
    chrom = str(selected[f"{role.lower()}_chrom"])
    strand = str(selected[f"{role.lower()}_strand"])
    start = int(selected[f"{role.lower()}_start"])
    end = int(selected[f"{role.lower()}_end"])
    coordinate = start + offset if strand == "+" else end - offset
    return chrom, coordinate, strand


def extract_atha_atom(
    genome: dict[str, str],
    blocks: list[dict[str, object]],
    role: str,
    core_start: int,
    core_end: int,
) -> tuple[str, dict[str, object]] | None:
    chrom1, cut1, strand1 = core_cut(blocks, role, core_start, "start")
    chrom2, cut2, strand2 = core_cut(blocks, role, core_end, "end")
    if chrom1 != chrom2 or strand1 != strand2:
        return None
    start, end = sorted((cut1, cut2))
    if end <= start or end - start > 1200:
        return None
    sequence = genome[chrom1][start:end]
    if strand1 == "-":
        sequence = reverse_complement(sequence)
    return sequence, {
        "chrom": chrom1,
        "start": start,
        "end": end,
        "strand": strand1,
    }


def target_atom(
    event_id: str,
    species: str,
    candidate_id: str,
    query_roles: tuple[str, ...],
    core_start: int,
    core_end: int,
    core_sequences: dict[str, str],
    hit_index: dict[tuple[str, str, str, str], list[dict[str, object]]],
    candidate_sequences: dict[str, dict[str, str]],
    dual: bool,
    endpoint_anchor: int = 20,
    dual_span_mode: str = "first",
) -> tuple[str | None, str]:
    projections = []
    for role in query_roles:
        clusters = cluster_hsps(
            hit_index.get((event_id, species, role, candidate_id), [])
        )
        if not clusters:
            return None, f"no_{role}_cluster"
        projection, failures = project_cluster(
            max(clusters, key=cluster_score), core_sequences[role]
        )
        if failures:
            return None, f"{role}_query_validation_failure"
        left = range(core_start, min(core_end, core_start + endpoint_anchor))
        right = range(max(core_start, core_end - endpoint_anchor), core_end)
        positions = list(left) + list(right)
        covered = [position for position in positions if position in projection]
        if not positions or len(covered) / len(positions) < 0.80:
            return None, f"{role}_endpoint_coverage_below_80pct"
        identity = sum(
            projection[position][0] == core_sequences[role][position].upper()
            for position in covered
        ) / len(covered)
        if identity < 0.60:
            return None, f"{role}_endpoint_identity_below_60pct"
        left_subject = [
            projection[position][1] for position in left if position in projection
        ]
        right_subject = [
            projection[position][1] for position in right if position in projection
        ]
        strand = "plus" if np.median(right_subject) > np.median(left_subject) else "minus"
        all_subject = left_subject + right_subject
        projections.append(
            (min(all_subject), max(all_subject) + 1, strand, identity)
        )
    if dual:
        if len(projections) != 2 or projections[0][2] != projections[1][2]:
            return None, "P0_dual_strand_discordant"
        if (
            abs(projections[0][0] - projections[1][0]) > 20
            or abs(projections[0][1] - projections[1][1]) > 20
        ):
            return None, "P0_dual_endpoint_discordant"
        if dual_span_mode == "union":
            start = min(row[0] for row in projections)
            end = max(row[1] for row in projections)
        elif dual_span_mode == "intersection":
            start = max(row[0] for row in projections)
            end = min(row[1] for row in projections)
        else:
            start, end = projections[0][0], projections[0][1]
        strand = projections[0][2]
    else:
        start, end, strand, _ = projections[0]
    sequence = candidate_sequences[species][candidate_id]
    if start < 0 or end > len(sequence) or end <= start or end - start > 1200:
        return None, "target_span_invalid"
    value = sequence[start:end]
    if strand == "minus":
        value = reverse_complement(value)
    return value.upper(), "PASS"


def nested_clade(names: list[str]) -> str:
    if not names:
        raise ValueError("empty clade")
    current = names[0]
    for name in names[1:]:
        current = f"({current},{name})"
    return current


def fixed_tree(names: set[str]) -> str:
    def copy_clade(role: str) -> str:
        ordered = [
            f"Atha_{role}",
            f"Alyrata_{role}",
            f"Bstricta_{role}",
            f"Dstrictus_{role}",
        ]
        return nested_clade([name for name in ordered if name in names])

    duplication = f"({copy_clade('P')},{copy_clade('D')})"
    # Wrap from the closest available P0 to progressively deeper P0 species.
    for species in ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea"):
        name = f"{species}_P0"
        if name in names:
            duplication = f"({name},{duplication})"
    return duplication + ";\n"


def pd_gap_calls(
    alignment: dict[str, str],
    atom: dict[str, object],
    raw_sequences: dict[str, str],
) -> list[dict[str, object]]:
    if "Atha_P" not in alignment or "Atha_D" not in alignment:
        return []
    p, d = alignment["Atha_P"], alignment["Atha_D"]
    raw_position = {"P": 0, "D": 0}
    column_raw = {"P": {}, "D": {}}
    for column in range(len(p)):
        for role, sequence in (("P", p), ("D", d)):
            if sequence[column] != "-":
                column_raw[role][column] = raw_position[role]
                raw_position[role] += 1
    calls = []
    column = 0
    while column < len(p):
        pbase, dbase = p[column], d[column]
        if (pbase == "-") == (dbase == "-"):
            column += 1
            continue
        role = "P" if pbase != "-" else "D"
        start_column = column
        while column < len(p):
            present = p[column] != "-" if role == "P" else d[column] != "-"
            absent = d[column] == "-" if role == "P" else p[column] == "-"
            if not (present and absent):
                break
            column += 1
        end_column = column
        length = sum(
            (p[index] if role == "P" else d[index]).upper() in "ACGT"
            for index in range(start_column, end_column)
        )
        if not 1 <= length <= 10:
            continue
        # Reject reciprocal/complex P-D gaps immediately adjacent to the call.
        complex_nearby = False
        for index in range(max(0, start_column - 5), min(len(p), end_column + 5)):
            if start_column <= index < end_column:
                continue
            if (p[index] == "-") != (d[index] == "-"):
                complex_nearby = True
                break
        if complex_nearby:
            continue
        raw_indices = [
            column_raw[role][index]
            for index in range(start_column, end_column)
            if index in column_raw[role]
        ]
        if len(raw_indices) != length or max(raw_indices) - min(raw_indices) + 1 != length:
            continue
        info = atom[f"atha_{role}_coordinates"]
        if info["strand"] == "+":
            genomic = [int(info["start"]) + index for index in raw_indices]
        else:
            genomic = [int(info["end"]) - 1 - index for index in raw_indices]
        gstart, gend = min(genomic), max(genomic) + 1
        states = {}
        for name, sequence in alignment.items():
            bases = sum(
                sequence[index].upper() in "ACGT"
                for index in range(start_column, end_column)
            )
            states[name] = (
                "PRESENT" if bases == length else "ABSENT" if bases == 0 else "PARTIAL"
            )
        segment = "".join(
            (p[index] if role == "P" else d[index]).upper()
            for index in range(start_column, end_column)
            if (p[index] if role == "P" else d[index]) != "-"
        )
        calls.append(
            {
                "unique_role": role,
                "chrom": info["chrom"],
                "start_0based": gstart,
                "end_0based": gend,
                "fragment_bp": length,
                "segment_sequence": segment,
                "states": states,
            }
        )
    return calls


def call_key(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row["unique_role"],
        row["chrom"],
        int(row["start_0based"]),
        int(row["end_0based"]),
        int(row["fragment_bp"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--pilot",
        type=Path,
        help="Override the legacy core_500 directory.",
    )
    parser.add_argument(
        "--inclusive-queue",
        type=Path,
        help="Custom length-ungated event/P0 queue; bypasses legacy event sets.",
    )
    parser.add_argument("--junction-anchor", type=int, default=60)
    parser.add_argument("--tile-bp", type=int, default=250)
    parser.add_argument("--tile-step", type=int, default=200)
    parser.add_argument(
        "--prank",
        type=Path,
        default=Path("prank"),
        help="PRANK executable (default: resolve 'prank' from PATH)",
    )
    parser.add_argument(
        "--age-bins",
        default="time1,time2,time3",
        help="Comma-separated controlled event age bins to analyze",
    )
    parser.add_argument(
        "--age-free-pd",
        action="store_true",
        help=(
            "Do not require a mapped post-duplication species.  P/D is already "
            "fixed by singleton-state synteny and any mapped P0 can polarize "
            "a locally comparable atom."
        ),
    )
    parser.add_argument(
        "--event-set",
        choices=("controlled55", "mapping_union"),
        default="controlled55",
        help=(
            "Use the formal 55-event endpoint or the nonredundant union of "
            "strict, partial-postduplication, and deeper-P0 mapping queues"
        ),
    )
    parser.add_argument(
        "--dual-p0-span",
        choices=("first", "union", "intersection"),
        default="first",
        help=(
            "How to cut a dual-query P0 target span; union is the symmetric "
            "local-MSA expansion control"
        ),
    )
    args = parser.parse_args()
    project, output = args.project, args.output
    selected_ages = set(args.age_bins.split(","))
    output.mkdir(parents=True, exist_ok=True)
    raw_dir, align_dir = output / "atomic_haplotypes", output / "alignments"
    raw_dir.mkdir(exist_ok=True)
    align_dir.mkdir(exist_ok=True)

    expansion = project / "08_event_inclusion_sensitivity"
    pilot = args.pilot or expansion / "core_500"
    if args.inclusive_queue:
        event_meta = {
            row["event_id"]: row
            for row in read_tsv(args.inclusive_queue)
            if row["age_bin"] in selected_ages
        }
        all_queues = dict(event_meta)
    else:
        event_meta = {}
        all_queues = {}
    queue_paths = (
        pilot / "pilot/strict_primary_event_queue.tsv",
        pilot / "pilot/partial_postdup_event_queue.tsv",
        expansion / "deeper_P0_sensitivity_queue.tsv",
    )
    if not args.inclusive_queue:
        for path in queue_paths:
            for row in read_tsv(path):
                if row["age_bin"] in selected_ages:
                    # The queue files are nested; first occurrence preserves the
                    # strongest tier: strict > partial postdup > deeper P0.
                    all_queues.setdefault(row["event_id"], row)
        if args.event_set == "controlled55":
            event_meta = {
                row["event_id"]: row
                for row in read_tsv(
                    expansion / "controlled_expansion_endpoint_events.tsv"
                )
                if row["age_bin"] in selected_ages
            }
        else:
            event_meta = dict(all_queues)
    events = {
        row["event_id"]: row
        for row in read_tsv(pilot / "inputs/high_priority_core_eligible.tsv")
        if row["event_id"] in event_meta
    }
    queues = {
        event_id: row
        for event_id, row in all_queues.items()
        if event_id in event_meta
    }
    mapping = {
        (row["event_id"], row["species"]): row
        for row in read_tsv(pilot / "outgroup_mapping/event_species_blastn_summary.tsv")
    }
    blocks_by_event: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(pilot / "core/homologous_core_blocks.tsv"):
        if row["event_id"] in events:
            blocks_by_event[row["event_id"]].append(row)
    genome = read_fasta(
        project / "01_reference/prepared_data/TAIR12.Col-CC.annotation_softmasked.fa"
    )
    core_fasta = read_fasta(pilot / "core/TAIR12_PD_homologous_cores.fa")
    candidate_sequences = {
        species: read_fasta(
            pilot / f"outgroup_mapping/candidate_regions/{species}.candidate_regions.fa"
        )
        for species in SPECIES
    }
    hit_index: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for species in SPECIES:
        with (
            pilot / f"outgroup_mapping/blastn_aligned/{species}.aligned_hits.tsv"
        ).open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                event_id = fields[0].split("|")[0]
                if event_id in events and fields[1].split("|")[0] == event_id:
                    hit_index[
                        (event_id, species, fields[0].split("|")[1], fields[1])
                    ].append(parse_hsp(fields))

    repeat_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with (
        project / "01_reference/prepared_data/TAIR12.annotated_repeats.merged.bed"
    ).open() as handle:
        for line in handle:
            chrom, start, end, *_ = line.rstrip().split("\t")
            repeat_intervals[chrom].append((int(start), int(end)))
    features = {kind: defaultdict(list) for kind in ("CDS", "exon", "gene")}
    with (
        project / "01_reference/prepared_data/TAIR12.Col-CC.annotation.gff3"
    ).open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) >= 5 and fields[2] in features:
                features[fields[2]][fields[0]].append(
                    (int(fields[3]) - 1, int(fields[4]))
                )

    atoms = []
    atom_sequences: dict[str, dict[str, str]] = {}
    target_qc = Counter()
    atom_no = 0
    for event_id, event in events.items():
        age = event_meta[event_id]["age_bin"]
        blocks = block_ranges(blocks_by_event[event_id])
        total = int(blocks[-1]["core_end"])
        intervals = set()
        for block in blocks[:-1]:
            boundary = int(block["core_end"])
            start = max(0, boundary - args.junction_anchor)
            end = min(total, boundary + args.junction_anchor)
            if end - start >= 80:
                intervals.add((start, end, "junction"))
        for start in range(0, total, args.tile_step):
            end = min(total, start + args.tile_bp)
            if end - start >= 100:
                intervals.add((start, end, "tile"))
        core_sequences = {
            role: core_fasta[f"{event_id}|{role}|{age}"] for role in ("P", "D")
        }
        for core_start, core_end, atom_type in sorted(intervals):
            p = extract_atha_atom(genome, blocks, "P", core_start, core_end)
            d = extract_atha_atom(genome, blocks, "D", core_start, core_end)
            if not p or not d:
                continue
            sequences = {"Atha_P": p[0], "Atha_D": d[0]}
            boundary_species = queues[event_id]["boundary_P0_species"]
            for species in SPECIES:
                row = mapping.get((event_id, species))
                if not row or row["mapping_status"] != "PASS":
                    continue
                if row["expected_locus_class"] == "preduplication_single_copy":
                    value, qc = target_atom(
                        event_id,
                        species,
                        row["P_to_P_candidate_id"],
                        ("P", "D"),
                        core_start,
                        core_end,
                        core_sequences,
                        hit_index,
                        candidate_sequences,
                        True,
                        dual_span_mode=args.dual_p0_span,
                    )
                    target_qc[("P0", qc)] += 1
                    if value:
                        sequences[f"{species}_P0"] = value
                else:
                    for role in ("P", "D"):
                        value, qc = target_atom(
                            event_id,
                            species,
                            row[f"{role}_to_{role}_candidate_id"],
                            (role,),
                            core_start,
                            core_end,
                            core_sequences,
                            hit_index,
                            candidate_sequences,
                            False,
                        )
                        target_qc[("postduplication", qc)] += 1
                        if value:
                            sequences[f"{species}_{role}"] = value
            has_postdup = any(
                f"{species}_P" in sequences and f"{species}_D" in sequences
                for species in SPECIES
            )
            has_p0 = any(f"{species}_P0" in sequences for species in SPECIES)
            if not has_p0 or (
                not args.age_free_pd and age != "time1" and not has_postdup
            ):
                continue
            atom_no += 1
            atom_id = f"SDATOM{atom_no:05d}"
            tree = fixed_tree(set(sequences))
            atom = {
                "atom_id": atom_id,
                "event_id": event_id,
                "age_bin": age,
                "atom_type": atom_type,
                "core_start_0based": core_start,
                "core_end_0based": core_end,
                "core_bp": core_end - core_start,
                "Atha_P_bp": len(p[0]),
                "Atha_D_bp": len(d[0]),
                "sequence_count": len(sequences),
                "sequence_names": ",".join(sequences),
                "boundary_P0_species": boundary_species,
                "tree_newick": tree.strip(),
                "atha_P_coordinates": p[1],
                "atha_D_coordinates": d[1],
            }
            atoms.append(atom)
            atom_sequences[atom_id] = sequences
            (raw_dir / f"{atom_id}.fa").write_text(fasta_text(sequences))
            (raw_dir / f"{atom_id}.tree.nwk").write_text(tree)

    mafft = shutil.which("mafft")
    muscle = shutil.which("muscle")
    prank = str(args.prank)
    if not mafft or not muscle or not args.prank.exists():
        raise SystemExit("MAFFT, MUSCLE5, and hvnlr PRANK are required")
    atom_lookup = {str(atom["atom_id"]): atom for atom in atoms}
    raw_calls: dict[tuple[str, str], list[dict[str, object]]] = {}
    for mode in ("MAFFT_LINSI", "MUSCLE5"):
        (align_dir / mode).mkdir(exist_ok=True)
    for atom in atoms:
        atom_id = str(atom["atom_id"])
        input_path = raw_dir / f"{atom_id}.fa"
        mafft_result = subprocess.run(
            [mafft, "--localpair", "--maxiterate", "1000", "--quiet", str(input_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        mafft_path = align_dir / "MAFFT_LINSI" / f"{atom_id}.fa"
        mafft_path.write_text(mafft_result.stdout)
        muscle_path = align_dir / "MUSCLE5" / f"{atom_id}.fa"
        subprocess.run(
            [muscle, "-align", str(input_path), "-output", str(muscle_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        for mode, path in (("MAFFT_LINSI", mafft_path), ("MUSCLE5", muscle_path)):
            raw_calls[(atom_id, mode)] = pd_gap_calls(
                read_fasta(path), atom, atom_sequences[atom_id]
            )

    stable_by_atom: dict[str, dict[tuple[object, ...], dict[str, dict[str, object]]]] = {}
    for atom in atoms:
        atom_id = str(atom["atom_id"])
        mafft_calls = {call_key(row): row for row in raw_calls[(atom_id, "MAFFT_LINSI")]}
        muscle_calls = {call_key(row): row for row in raw_calls[(atom_id, "MUSCLE5")]}
        shared = set(mafft_calls) & set(muscle_calls)
        if shared:
            stable_by_atom[atom_id] = {
                key: {"MAFFT_LINSI": mafft_calls[key], "MUSCLE5": muscle_calls[key]}
                for key in shared
            }

    (align_dir / "PRANK_FIXED_TREE").mkdir(exist_ok=True)
    prank_failures = Counter()
    for atom_id, calls in stable_by_atom.items():
        prefix = align_dir / "PRANK_FIXED_TREE" / atom_id
        command = [
            prank,
            f"-d={raw_dir / f'{atom_id}.fa'}",
            f"-t={raw_dir / f'{atom_id}.tree.nwk'}",
            f"-o={prefix}",
            "-DNA",
            "+F",
            "-once",
            "-showanc",
            "-showevents",
            "-quiet",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        best = Path(f"{prefix}.best.fas")
        if result.returncode or not best.exists():
            prank_failures["execution_failure"] += 1
            continue
        prank_calls = {
            call_key(row): row
            for row in pd_gap_calls(
                read_fasta(best), atom_lookup[atom_id], atom_sequences[atom_id]
            )
        }
        for key in list(calls):
            if key in prank_calls:
                calls[key]["PRANK_FIXED_TREE"] = prank_calls[key]

    inferred = []
    seen_keys = {}
    for atom_id, calls in stable_by_atom.items():
        atom = atom_lookup[atom_id]
        for key, by_mode in calls.items():
            if "PRANK_FIXED_TREE" not in by_mode:
                continue
            state_names = set.intersection(
                *(set(row["states"]) for row in by_mode.values())
            )
            consensus = {}
            for name in state_names:
                values = [row["states"][name] for row in by_mode.values()]
                consensus[name] = values[0] if len(set(values)) == 1 else "DISCORDANT"
            base = by_mode["PRANK_FIXED_TREE"]
            unique_role = str(base["unique_role"])
            opposite_role = "D" if unique_role == "P" else "P"
            concordant, discordant, tested = [], [], []
            for species in SPECIES:
                unique = consensus.get(f"{species}_{unique_role}")
                opposite = consensus.get(f"{species}_{opposite_role}")
                if unique in {"PRESENT", "ABSENT"} and opposite in {"PRESENT", "ABSENT"}:
                    tested.append(species)
                    if unique == "PRESENT" and opposite == "ABSENT":
                        concordant.append(species)
                    elif unique != opposite:
                        discordant.append(species)
            p0_states = {
                species: consensus[f"{species}_P0"]
                for species in SPECIES
                if consensus.get(f"{species}_P0") in {"PRESENT", "ABSENT"}
            }
            boundary = str(atom["boundary_P0_species"])
            boundary_state = p0_states.get(boundary, "UNKNOWN")
            deeper = {s: value for s, value in p0_states.items() if s != boundary}
            if boundary_state != "UNKNOWN":
                p0_state, p0_source = boundary_state, boundary
            elif len(deeper) >= 2 and len(set(deeper.values())) == 1:
                p0_state, p0_source = next(iter(deeper.values())), ",".join(sorted(deeper))
            else:
                p0_state, p0_source = "UNKNOWN", "NA"
            p0_discordance = (
                p0_state != "UNKNOWN"
                and any(value != p0_state for value in deeper.values())
            )
            age = str(atom["age_bin"])
            time1_multispecies_P0 = (
                age == "time1"
                and boundary_state != "UNKNOWN"
                and any(value == boundary_state for value in deeper.values())
                and not p0_discordance
            )
            if age == "time1" and not time1_multispecies_P0:
                status, direction = "UNRESOLVED", "time1_multispecies_P0_insufficient"
            elif age != "time1" and discordant:
                status, direction = "UNRESOLVED", "postduplication_state_discordant"
            elif age != "time1" and not concordant:
                status, direction = "UNRESOLVED", "not_replicated_postduplication"
            elif p0_state == "UNKNOWN":
                status, direction = "UNRESOLVED", "P0_state_unresolved"
            elif p0_discordance:
                status, direction = "UNRESOLVED", "P0_state_discordant"
            elif p0_state == "ABSENT":
                status, direction = "POLARIZED", f"{unique_role}_insertion"
            else:
                status, direction = "POLARIZED", f"{opposite_role}_deletion"
            chrom = str(base["chrom"])
            start, end = int(base["start_0based"]), int(base["end_0based"])
            segment_original = genome[chrom][start:end]
            repeat = interval_overlap(repeat_intervals, chrom, start, end)
            context = genome[chrom][max(0, start - 10) : min(len(genome[chrom]), end + 10)]
            low_complexity = longest_homopolymer(context) >= 6 or sequence_entropy(context) < 1.2
            sequence_qc = (
                all(base_value in "ACGT" for base_value in segment_original)
                and not repeat
                and not low_complexity
            )
            if interval_overlap(features["CDS"], chrom, start, end):
                genomic_context = "CDS"
            elif interval_overlap(features["exon"], chrom, start, end):
                genomic_context = "exon_nonCDS"
            elif interval_overlap(features["gene"], chrom, start, end):
                genomic_context = "intron_or_gene_body"
            else:
                genomic_context = "intergenic"
            evidence = (
                "PRIMARY"
                if status == "POLARIZED" and sequence_qc and boundary_state != "UNKNOWN"
                else "SENSITIVITY"
                if status == "POLARIZED"
                else "UNRESOLVED"
            )
            row = {
                "denovo_indel_id": "PENDING",
                "event_id": atom["event_id"],
                "age_bin": atom["age_bin"],
                "atom_id": atom_id,
                "atom_type": atom["atom_type"],
                "unique_role": unique_role,
                "chrom": chrom,
                "start_0based": start,
                "end_0based": end,
                "fragment_bp": base["fragment_bp"],
                "segment_sequence": base["segment_sequence"],
                "MAFFT_MUSCLE_PRANK_exact_coordinate_concordance": "PASS",
                "postduplication_species_tested": ",".join(tested),
                "postduplication_concordant_species": ",".join(concordant),
                "postduplication_discordant_species": ",".join(discordant),
                "time1_multispecies_P0_support": (
                    "PASS" if time1_multispecies_P0 else "NA"
                ),
                "boundary_P0_species": boundary,
                "boundary_P0_state": boundary_state,
                "deeper_P0_states": ",".join(
                    f"{species}:{value}" for species, value in sorted(deeper.items())
                ),
                "P0_state": p0_state,
                "P0_source_species": p0_source,
                "direction_status": status,
                "parsimonious_direction": direction,
                "fully_uppercase": all(value in "ACGT" for value in segment_original),
                "annotated_repeat_overlap": repeat,
                "low_complexity_context": low_complexity,
                "genomic_context": genomic_context,
                "coding_size_class": (
                    "frameshift_candidate"
                    if genomic_context == "CDS" and int(base["fragment_bp"]) % 3
                    else "inframe_candidate"
                    if genomic_context == "CDS"
                    else "noncoding_or_unresolved"
                ),
                "evidence_tier": evidence,
            }
            dedup_key = (
                row["event_id"],
                row["unique_role"],
                row["chrom"],
                row["start_0based"],
                row["end_0based"],
            )
            current = seen_keys.get(dedup_key)
            rank = {"PRIMARY": 2, "SENSITIVITY": 1, "UNRESOLVED": 0}
            if current is None or rank[row["evidence_tier"]] > rank[current["evidence_tier"]]:
                seen_keys[dedup_key] = row

    inferred = sorted(
        seen_keys.values(),
        key=lambda row: (row["event_id"], row["chrom"], int(row["start_0based"])),
    )
    # Overlapping atomic windows can split one contiguous event into adjacent
    # calls. Merge only touching/overlapping calls with the same fully inferred
    # branch direction and evidence tier; never merge across an intervening bp.
    merged = []
    for row in sorted(
        inferred,
        key=lambda value: (
            value["event_id"],
            value["chrom"],
            value["unique_role"],
            value["parsimonious_direction"],
            value["evidence_tier"],
            int(value["start_0based"]),
        ),
    ):
        if (
            merged
            and merged[-1]["event_id"] == row["event_id"]
            and merged[-1]["chrom"] == row["chrom"]
            and merged[-1]["unique_role"] == row["unique_role"]
            and merged[-1]["parsimonious_direction"] == row["parsimonious_direction"]
            and merged[-1]["evidence_tier"] == row["evidence_tier"]
            and int(row["start_0based"]) <= int(merged[-1]["end_0based"])
            and max(int(merged[-1]["end_0based"]), int(row["end_0based"]))
            - min(int(merged[-1]["start_0based"]), int(row["start_0based"]))
            <= 10
        ):
            previous = merged[-1]
            previous["start_0based"] = min(
                int(previous["start_0based"]), int(row["start_0based"])
            )
            previous["end_0based"] = max(
                int(previous["end_0based"]), int(row["end_0based"])
            )
            previous["fragment_bp"] = (
                int(previous["end_0based"]) - int(previous["start_0based"])
            )
            previous["segment_sequence"] = genome[str(previous["chrom"])][
                int(previous["start_0based"]) : int(previous["end_0based"])
            ].upper()
            previous["atom_id"] = ",".join(
                sorted(set(str(previous["atom_id"]).split(",") + [str(row["atom_id"])]))
            )
        else:
            merged.append(dict(row))
    inferred = sorted(
        merged,
        key=lambda row: (row["event_id"], row["chrom"], int(row["start_0based"])),
    )
    for index, row in enumerate(inferred, 1):
        row["denovo_indel_id"] = f"DNIND{index:05d}"
    atom_rows = []
    for atom in atoms:
        atom_rows.append(
            {
                key: value
                for key, value in atom.items()
                if not key.startswith("atha_")
            }
        )
    write_tsv(output / "atomic_region_manifest.tsv", atom_rows)
    write_tsv(
        output / "target_extraction_qc_summary.tsv",
        [
            {"phylogenetic_role": key[0], "status": key[1], "tests": value}
            for key, value in sorted(target_qc.items())
        ],
    )
    write_tsv(output / "denovo_microindel_inference.tsv", inferred)
    flow = [
        {"stage": "atomic_regions_with_required_outgroup_and_P0", "count": len(atoms)},
        {"stage": "atoms_with_MAFFT_MUSCLE_shared_gap", "count": len(stable_by_atom)},
        {
            "stage": "unique_three_aligner_gap_candidates",
            "count": len(inferred),
        },
        {
            "stage": "polarized_all_QC_tiers",
            "count": sum(row["direction_status"] == "POLARIZED" for row in inferred),
        },
        {
            "stage": "primary_polarized",
            "count": sum(row["evidence_tier"] == "PRIMARY" for row in inferred),
        },
    ]
    write_tsv(output / "analysis_flow.tsv", flow)

    endpoint_rows = []
    def event_callable(row: dict[str, str]) -> int:
        return int(
            row.get(
                "primary_bidir_callable_sites",
                row.get("jointly_callable_uppercase_acgt_bp", 0),
            )
        )

    all_callable_bp = sum(event_callable(row) for row in event_meta.values())
    for endpoint, rows in (
        ("primary", [row for row in inferred if row["evidence_tier"] == "PRIMARY"]),
        (
            "primary_plus_sensitivity",
            [row for row in inferred if row["direction_status"] == "POLARIZED"],
        ),
    ):
        p_count = sum(str(row["parsimonious_direction"]).startswith("P_") for row in rows)
        d_count = sum(str(row["parsimonious_direction"]).startswith("D_") for row in rows)
        total = p_count + d_count
        endpoint_rows.append(
            {
                "endpoint": endpoint,
                "P_branch_indels": p_count,
                "D_branch_indels": d_count,
                "total_indels": total,
                "D_to_P_count_ratio": f"{d_count / p_count:.6f}" if p_count else "Inf",
                "exact_binomial_two_sided_p": (
                    f"{exact_two_sided_binomial(d_count, total):.8g}" if total else "NA"
                ),
                "events_with_polarized_indels": len({row["event_id"] for row in rows}),
                "ancestral_callable_bp_all_selected_events": all_callable_bp,
                "P_indels_per_kb_callable": f"{1000 * p_count / all_callable_bp:.8f}",
                "D_indels_per_kb_callable": f"{1000 * d_count / all_callable_bp:.8f}",
                "note": "site_count_test_is_descriptive_events_are_the_independent_units",
            }
        )
    write_tsv(output / "PD_denovo_microindel_statistics.tsv", endpoint_rows)

    primary_for_events = [
        row for row in inferred if row["evidence_tier"] == "PRIMARY"
    ]
    event_rows = []
    d_greater = p_greater = ties = 0
    for event_id, meta in sorted(event_meta.items()):
        rows = [row for row in primary_for_events if row["event_id"] == event_id]
        p_count = sum(
            str(row["parsimonious_direction"]).startswith("P_") for row in rows
        )
        d_count = sum(
            str(row["parsimonious_direction"]).startswith("D_") for row in rows
        )
        callable_bp = event_callable(meta)
        if d_count > p_count:
            d_greater += 1
        elif p_count > d_count:
            p_greater += 1
        else:
            ties += 1
        event_rows.append(
            {
                "event_id": event_id,
                "age_bin": meta["age_bin"],
                "ancestral_callable_bp": callable_bp,
                "P_branch_indels": p_count,
                "D_branch_indels": d_count,
                "P_indels_per_kb": f"{1000 * p_count / callable_bp:.8f}",
                "D_indels_per_kb": f"{1000 * d_count / callable_bp:.8f}",
                "D_minus_P_indels": d_count - p_count,
            }
        )
    write_tsv(output / "event_level_PD_microindel_rates.tsv", event_rows)
    non_ties = d_greater + p_greater
    write_tsv(
        output / "event_level_paired_statistical_summary.tsv",
        [
            {
                "endpoint": f"primary_{len(event_meta)}_event_paired",
                "events_D_greater": d_greater,
                "events_P_greater": p_greater,
                "events_tied": ties,
                "non_tied_events": non_ties,
                "exact_two_sided_sign_test_p": (
                    f"{exact_two_sided_binomial(d_greater, non_ties):.8g}"
                    if non_ties
                    else "NA"
                ),
                "independent_unit": "SD_event",
            }
        ],
    )
    age_rows = []
    for age in sorted({row["age_bin"] for row in event_rows}):
        age_event_rows = [row for row in event_rows if row["age_bin"] == age]
        if not age_event_rows:
            continue
        age_d = sum(
            int(row["D_branch_indels"]) > int(row["P_branch_indels"])
            for row in age_event_rows
        )
        age_p = sum(
            int(row["P_branch_indels"]) > int(row["D_branch_indels"])
            for row in age_event_rows
        )
        age_ties = len(age_event_rows) - age_d - age_p
        age_non_ties = age_d + age_p
        age_primary = [
            row for row in primary_for_events if row["age_bin"] == age
        ]
        age_rows.append(
            {
                "age_bin": age,
                "controlled_events": len(age_event_rows),
                "primary_microindels": len(age_primary),
                "P_branch_microindels": sum(
                    str(row["parsimonious_direction"]).startswith("P_")
                    for row in age_primary
                ),
                "D_branch_microindels": sum(
                    str(row["parsimonious_direction"]).startswith("D_")
                    for row in age_primary
                ),
                "events_D_greater": age_d,
                "events_P_greater": age_p,
                "events_tied": age_ties,
                "exact_two_sided_sign_test_p": (
                    f"{exact_two_sided_binomial(age_d, age_non_ties):.8g}"
                    if age_non_ties
                    else "NA"
                ),
            }
        )
    write_tsv(output / "age_stratified_PD_microindel_summary.tsv", age_rows)

    primary = [row for row in inferred if row["evidence_tier"] == "PRIMARY"]
    direction_counts = Counter(row["parsimonious_direction"] for row in primary)
    context_counts = Counter(
        (row["genomic_context"], row["parsimonious_direction"]) for row in primary
    )
    write_tsv(
        output / "primary_functional_context_summary.tsv",
        [
            {
                "genomic_context": key[0],
                "direction": key[1],
                "microindels": value,
            }
            for key, value in sorted(context_counts.items())
        ],
    )
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    directions = ("P_insertion", "P_deletion", "D_insertion", "D_deletion")
    axes[0].bar(
        range(4),
        [direction_counts[value] for value in directions],
        color=["#3973ac", "#78a6d0", "#c45145", "#db8c83"],
    )
    axes[0].set_xticks(range(4), ("P ins", "P del", "D ins", "D del"), rotation=25)
    axes[0].set_ylabel("Primary de novo micro-indels")
    axes[0].set_title("A  Polarized direction", loc="left")
    for role, color in (("P", "#3973ac"), ("D", "#c45145")):
        values = [
            sum(
                str(row["parsimonious_direction"]).startswith(f"{role}_")
                and int(row["fragment_bp"]) == length
                for row in primary
            )
            for length in range(1, 11)
        ]
        axes[1].plot(range(1, 11), values, marker="o", color=color, label=role)
    axes[1].set_xticks(range(1, 11))
    axes[1].set_xlabel("Length (bp)")
    axes[1].set_ylabel("Primary de novo micro-indels")
    axes[1].set_title("B  Length spectrum", loc="left")
    axes[1].legend(frameon=False)
    figure_dir = output / "figures"
    figure_dir.mkdir(exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(figure_dir / f"denovo_microindel_summary.{suffix}", dpi=300)
    plt.close(fig)
    print(
        f"atoms={len(atoms)} shared_atoms={len(stable_by_atom)} "
        f"three_aligner_candidates={len(inferred)} primary={len(primary)}"
    )


if __name__ == "__main__":
    main()
