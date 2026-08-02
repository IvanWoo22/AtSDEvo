#!/usr/bin/env python3
"""Recover locally comparable P/D/P0 SNP blocks from the mapping-queue union.

The workflow replaces an event-wide P/D mismatch cutoff with a local MSA rule.
MAFFT and MUSCLE must map P, D, and P0 to the same P/D genomic coordinate pair.
P/D gap runs > ``--large-gap-bp`` and their flanks are excluded.  Remaining
sites must lie in a local window with P/D mismatch <= ``--max-local-mismatch``.
Events with sufficient two-aligner sequence are then checked with fixed-tree
PRANK; final sites require exact coordinate, ancestral-base, and class agreement
across all three aligners.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif line:
                records[name].append(line)
    return {name: "".join(parts) for name, parts in records.items()}


def exact_binomial(successes: int, trials: int) -> float:
    if not trials:
        return math.nan
    observed = math.comb(trials, successes)
    return min(
        1.0,
        sum(
            math.comb(trials, value)
            for value in range(trials + 1)
            if math.comb(trials, value) <= observed
        )
        / (2**trials),
    )


def large_gap_mask(p: str, d: str, minimum: int, flank: int) -> set[int]:
    masked: set[int] = set()
    column = 0
    while column < len(p):
        if (p[column] == "-") == (d[column] == "-"):
            column += 1
            continue
        start = column
        present_bases = 0
        while (
            column < len(p)
            and (p[column] == "-") != (d[column] == "-")
        ):
            present_bases += (p[column] != "-") or (d[column] != "-")
            column += 1
        if present_bases >= minimum:
            masked.update(
                range(max(0, start - flank), min(len(p), column + flank))
            )
    return masked


def classify(p: str, d: str, ancestral: str) -> str:
    if p == d == ancestral:
        return "invariant"
    if d == ancestral and p != ancestral:
        return "P_specific"
    if p == ancestral and d != ancestral:
        return "D_specific"
    if p == d and p != ancestral:
        return "shared_PD"
    return "tri_allelic_unresolved"


def genomic_coordinate(info: dict[str, object], raw_index: int) -> int:
    if info["strand"] == "+":
        return int(info["start"]) + raw_index
    return int(info["end"]) - 1 - raw_index


def alignment_sites(
    alignment: dict[str, str],
    raw: dict[str, str],
    atom: dict[str, object],
    large_gap_bp: int,
    gap_flank: int,
    local_radius: int,
    max_local_mismatch: float,
    min_local_pairs: int,
    p0_rule: str,
) -> dict[tuple[object, ...], dict[str, object]]:
    if not {"Atha_P", "Atha_D"} <= set(alignment):
        return {}
    lengths = {len(value) for value in alignment.values()}
    if len(lengths) != 1:
        return {}
    p, d = alignment["Atha_P"], alignment["Atha_D"]
    boundary = str(atom["boundary_P0_species"])
    boundary_name = f"{boundary}_P0"
    if boundary_name not in alignment:
        return {}
    boundary_index = SPECIES.index(boundary)
    deeper_names = [
        f"{species}_P0"
        for species in SPECIES[boundary_index + 1 :]
        if f"{species}_P0" in alignment
    ]
    if p0_rule == "multispecies" and not deeper_names:
        return {}

    raw_positions = {name: 0 for name in alignment}
    column_raw: dict[str, dict[int, int]] = {
        name: {} for name in alignment
    }
    for column in range(len(p)):
        for name, sequence in alignment.items():
            if sequence[column] != "-":
                column_raw[name][column] = raw_positions[name]
                raw_positions[name] += 1
    if any(raw_positions[name] != len(raw[name]) for name in raw):
        return {}

    boundary_alignment = alignment[boundary_name]
    # Remove both present-day P/D structural gaps and large P0-to-copy gaps;
    # SNPs immediately flanking either class are alignment-sensitive.
    excluded = (
        large_gap_mask(p, d, large_gap_bp, gap_flank)
        | large_gap_mask(p, boundary_alignment, large_gap_bp, gap_flank)
        | large_gap_mask(d, boundary_alignment, large_gap_bp, gap_flank)
    )
    sites = {}
    for column in range(len(p)):
        if column in excluded:
            continue
        pbase, dbase = p[column].upper(), d[column].upper()
        abase = alignment[boundary_name][column].upper()
        if pbase not in "ACGT" or dbase not in "ACGT" or abase not in "ACGT":
            continue
        p_raw = column_raw["Atha_P"].get(column)
        d_raw = column_raw["Atha_D"].get(column)
        if p_raw is None or d_raw is None:
            continue
        # Preserve TAIR12 soft-mask as an exclusion criterion.
        if (
            raw["Atha_P"][p_raw] not in "ACGT"
            or raw["Atha_D"][d_raw] not in "ACGT"
        ):
            continue

        left, right = max(0, column - local_radius), min(
            len(p), column + local_radius + 1
        )
        paired = mismatches = 0
        for index in range(left, right):
            if index in excluded:
                continue
            pb, db = p[index].upper(), d[index].upper()
            if pb in "ACGT" and db in "ACGT":
                paired += 1
                mismatches += pb != db
        if paired < min_local_pairs or mismatches / paired > max_local_mismatch:
            continue

        deeper_states = [
            alignment[name][column].upper()
            for name in deeper_names
            if alignment[name][column].upper() in "ACGT"
        ]
        p0_conflict = any(base != abase for base in deeper_states)
        if p0_rule == "multispecies" and (
            not deeper_states or p0_conflict
        ):
            continue
        if p0_rule == "boundary_no_conflict" and p0_conflict:
            continue

        p_info = atom["atha_P_coordinates"]
        d_info = atom["atha_D_coordinates"]
        p_coord = genomic_coordinate(p_info, p_raw)
        d_coord = genomic_coordinate(d_info, d_raw)
        key = (
            p_info["chrom"],
            p_coord,
            d_info["chrom"],
            d_coord,
        )
        sites[key] = {
            "event_id": atom["event_id"],
            "age_bin": atom["age_bin"],
            "atom_id": atom["atom_id"],
            "alignment_column_1based": column + 1,
            "P_chrom": p_info["chrom"],
            "P_position_0based": p_coord,
            "D_chrom": d_info["chrom"],
            "D_position_0based": d_coord,
            "P_base": pbase,
            "D_base": dbase,
            "ancestral_base": abase,
            "polarized_class": classify(pbase, dbase, abase),
            "boundary_P0_species": boundary,
            "deeper_P0_species_support": ",".join(
                name.removesuffix("_P0")
                for name in deeper_names
                if alignment[name][column].upper() == abase
            ),
            "local_paired_bases": paired,
            "local_PD_mismatch_fraction": f"{mismatches / paired:.8f}",
            "large_gap_and_flank_filter": "PASS",
            "uppercase_PD_filter": "PASS",
            "multispecies_P0_filter": (
                "CONFLICT"
                if p0_conflict
                else "PASS"
                if deeper_states
                else "NOT_AVAILABLE"
            ),
            "P0_rule": p0_rule,
            "deeper_P0_conflict": "YES" if p0_conflict else "NO",
        }
    return sites


def reconstruct_atoms(
    project: Path, atom_root: Path, pilot: Path
) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, str]]]:
    workflow = (
        project
        / "09_variant_type_analysis_55"
        / "workflow_scripts"
    )
    sys.path.insert(0, str(workflow))
    from analyze_denovo_msa_microindels import (  # type: ignore
        block_ranges,
        extract_atha_atom,
    )

    blocks: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(pilot / "core/homologous_core_blocks.tsv"):
        blocks[row["event_id"]].append(row)
    genome = read_fasta(
        project / "01_reference/prepared_data/TAIR12.Col-CC.annotation_softmasked.fa"
    )
    atoms = {}
    raw_sequences = {}
    for manifest in read_tsv(atom_root / "atomic_region_manifest.tsv"):
        atom_id, event_id = manifest["atom_id"], manifest["event_id"]
        normalized = block_ranges(blocks[event_id])
        p = extract_atha_atom(
            genome,
            normalized,
            "P",
            int(manifest["core_start_0based"]),
            int(manifest["core_end_0based"]),
        )
        d = extract_atha_atom(
            genome,
            normalized,
            "D",
            int(manifest["core_start_0based"]),
            int(manifest["core_end_0based"]),
        )
        if not p or not d:
            continue
        raw = read_fasta(atom_root / f"atomic_haplotypes/{atom_id}.fa")
        if p[0] != raw["Atha_P"] or d[0] != raw["Atha_D"]:
            raise AssertionError(f"atom reconstruction mismatch: {atom_id}")
        atoms[atom_id] = {
            **manifest,
            "atha_P_coordinates": p[1],
            "atha_D_coordinates": d[1],
        }
        raw_sequences[atom_id] = raw
    return atoms, raw_sequences


def summarize_events(
    event_ids: set[str],
    sites: dict[tuple[object, ...], dict[str, object]],
    queue_meta: dict[str, dict[str, str]],
    threshold: int,
    endpoint: str,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    by_event: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in sites.values():
        by_event[str(row["event_id"])].append(row)
    rows = []
    d_greater = p_greater = ties = 0
    for event_id in sorted(event_ids):
        values = by_event.get(event_id, [])
        if len(values) < threshold:
            continue
        counts = Counter(str(row["polarized_class"]) for row in values)
        p_count, d_count = counts["P_specific"], counts["D_specific"]
        if d_count > p_count:
            d_greater += 1
        elif p_count > d_count:
            p_greater += 1
        else:
            ties += 1
        rows.append(
            {
                "event_id": event_id,
                "age_bin": queue_meta[event_id]["age_bin"],
                "queue_tier": queue_meta[event_id]["analysis_tier"],
                "whole_event_PD_mismatch_pct": queue_meta[event_id][
                    "present_day_PD_mismatch_pct_callable"
                ],
                "local_MSA_callable_sites": len(values),
                "P_specific_SNP": p_count,
                "D_specific_SNP": d_count,
                "D_to_P_ratio": (
                    f"{d_count / p_count:.8f}" if p_count else "Inf"
                ),
                "D_minus_P": d_count - p_count,
            }
        )
    non_ties = d_greater + p_greater
    summary = {
        "endpoint": endpoint,
        "minimum_local_MSA_callable_sites": threshold,
        "events": len(rows),
        "time1": sum(row["age_bin"] == "time1" for row in rows),
        "time2": sum(row["age_bin"] == "time2" for row in rows),
        "time3": sum(row["age_bin"] == "time3" for row in rows),
        "age_free": sum(row["age_bin"] == "age_free" for row in rows),
        "callable_sites": sum(int(row["local_MSA_callable_sites"]) for row in rows),
        "P_specific_SNP": sum(int(row["P_specific_SNP"]) for row in rows),
        "D_specific_SNP": sum(int(row["D_specific_SNP"]) for row in rows),
        "events_D_greater": d_greater,
        "events_P_greater": p_greater,
        "events_tied": ties,
        "event_sign_test_p": (
            f"{exact_binomial(d_greater, non_ties):.8g}"
            if non_ties else "NA"
        ),
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--atom-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--pilot",
        type=Path,
        help="Override the legacy core_500 directory.",
    )
    parser.add_argument(
        "--inclusive-queue",
        type=Path,
        help="Custom length-ungated event/P0 queue.",
    )
    parser.add_argument("--large-gap-bp", type=int, default=11)
    parser.add_argument("--gap-flank", type=int, default=10)
    parser.add_argument("--local-radius", type=int, default=25)
    parser.add_argument("--max-local-mismatch", type=float, default=0.40)
    parser.add_argument("--min-local-pairs", type=int, default=20)
    parser.add_argument(
        "--p0-rule",
        choices=("boundary", "boundary_no_conflict", "multispecies"),
        default="boundary_no_conflict",
        help=(
            "boundary matches the former SNP endpoint; boundary_no_conflict "
            "rejects deeper-P0 conflicts when observable; multispecies also "
            "requires at least one concordant deeper P0"
        ),
    )
    parser.add_argument("--prank-event-threshold", type=int, default=200)
    parser.add_argument(
        "--site-thresholds",
        default="20,50,100,200,500",
        help="Comma-separated local callable-site endpoint thresholds.",
    )
    parser.add_argument(
        "--prank",
        default="prank",
        help="PRANK executable (default: resolve 'prank' from PATH)",
    )
    parser.add_argument(
        "--prank-cache",
        type=Path,
        help="Optional directory containing reusable fixed-tree PRANK alignments",
    )
    args = parser.parse_args()
    prank = shutil.which(args.prank) or (
        args.prank if Path(args.prank).is_file() else None
    )
    if not prank:
        raise SystemExit("PRANK was not found; add it to PATH or pass --prank")
    args.output.mkdir(parents=True, exist_ok=True)

    pilot = (
        args.pilot
        or args.project / "08_event_inclusion_sensitivity/core_500"
    )
    atoms, raw_sequences = reconstruct_atoms(args.project, args.atom_root, pilot)
    queue_meta = {}
    if args.inclusive_queue:
        queue_meta = {
            row["event_id"]: row for row in read_tsv(args.inclusive_queue)
        }
    else:
        for path in (
            args.project
            / "08_event_inclusion_sensitivity/core_500/pilot/strict_primary_event_queue.tsv",
            args.project
            / "08_event_inclusion_sensitivity/core_500/pilot/partial_postdup_event_queue.tsv",
            args.project
            / "08_event_inclusion_sensitivity/deeper_P0_sensitivity_queue.tsv",
        ):
            for row in read_tsv(path):
                queue_meta.setdefault(row["event_id"], row)

    mode_sites: dict[tuple[str, str], dict[tuple[object, ...], dict[str, object]]] = {}
    for atom_id, atom in atoms.items():
        for mode in ("MAFFT_LINSI", "MUSCLE5"):
            path = args.atom_root / f"alignments/{mode}/{atom_id}.fa"
            mode_sites[(atom_id, mode)] = alignment_sites(
                read_fasta(path),
                raw_sequences[atom_id],
                atom,
                args.large_gap_bp,
                args.gap_flank,
                args.local_radius,
                args.max_local_mismatch,
                args.min_local_pairs,
                args.p0_rule,
            )

    two_aligner = {}
    contributing_atoms: dict[str, set[str]] = defaultdict(set)
    for atom_id, atom in atoms.items():
        mafft = mode_sites[(atom_id, "MAFFT_LINSI")]
        muscle = mode_sites[(atom_id, "MUSCLE5")]
        for key in set(mafft) & set(muscle):
            left, right = mafft[key], muscle[key]
            if (
                left["ancestral_base"] == right["ancestral_base"]
                and left["polarized_class"] == right["polarized_class"]
            ):
                dedup = (atom["event_id"], *key)
                two_aligner.setdefault(dedup, left)
                contributing_atoms[str(atom["event_id"])].add(atom_id)

    two_counts = Counter(str(row["event_id"]) for row in two_aligner.values())
    prank_events = {
        event_id
        for event_id, count in two_counts.items()
        if count >= args.prank_event_threshold
    }
    prank_dir = args.prank_cache or args.output / "PRANK_FIXED_TREE"
    prank_dir.mkdir(exist_ok=True)
    prank_sites = {}
    failures = Counter()
    for event_id in sorted(prank_events):
        for atom_id in sorted(contributing_atoms[event_id]):
            prefix = prank_dir / atom_id
            best = Path(f"{prefix}.best.fas")
            if not best.exists():
                command = [
                    prank,
                    f"-d={args.atom_root / f'atomic_haplotypes/{atom_id}.fa'}",
                    f"-t={args.atom_root / f'atomic_haplotypes/{atom_id}.tree.nwk'}",
                    f"-o={prefix}",
                    "-DNA",
                    "+F",
                    "-once",
                    "-showanc",
                    "-showevents",
                    "-quiet",
                ]
                result = subprocess.run(
                    command, capture_output=True, text=True
                )
                if result.returncode or not best.exists():
                    failures["PRANK_execution_failure"] += 1
                    continue
            prank_sites[atom_id] = alignment_sites(
                read_fasta(best),
                raw_sequences[atom_id],
                atoms[atom_id],
                args.large_gap_bp,
                args.gap_flank,
                args.local_radius,
                args.max_local_mismatch,
                args.min_local_pairs,
                args.p0_rule,
            )

    three_aligner = {}
    for dedup, row in two_aligner.items():
        event_id = str(row["event_id"])
        if event_id not in prank_events:
            continue
        key = tuple(dedup[1:])
        # A site may occur in overlapping atoms; accept if any contributing
        # atom reproduces the same coordinate, ancestral base, and class.
        for atom_id in contributing_atoms[event_id]:
            candidate = prank_sites.get(atom_id, {}).get(key)
            if (
                candidate
                and candidate["ancestral_base"] == row["ancestral_base"]
                and candidate["polarized_class"] == row["polarized_class"]
            ):
                accepted = dict(row)
                accepted["ASR_atom_id"] = atom_id
                accepted["ASR_alignment_column_1based"] = candidate[
                    "alignment_column_1based"
                ]
                three_aligner[dedup] = accepted
                break

    write_tsv(args.output / "two_aligner_local_callable_sites.tsv", list(two_aligner.values()))
    write_tsv(
        args.output / "three_aligner_local_callable_sites.tsv",
        list(three_aligner.values()),
    )
    event_ids = set(queue_meta)
    summaries = []
    for label, site_set in (
        ("MAFFT_MUSCLE_local_MSA", two_aligner),
        ("MAFFT_MUSCLE_PRANK_local_MSA", three_aligner),
    ):
        for threshold in (
            int(value) for value in args.site_thresholds.split(",")
        ):
            rows, summary = summarize_events(
                event_ids,
                site_set,
                queue_meta,
                threshold,
                label,
            )
            write_tsv(
                args.output / f"{label}.ge{threshold}.event_metrics.tsv",
                rows,
            )
            summaries.append(summary)
    write_tsv(args.output / "local_MSA_endpoint_summary.tsv", summaries)
    write_tsv(
        args.output / "workflow_qc.tsv",
        [
            {"metric": "mapping_union_events", "value": len(queue_meta)},
            {"metric": "events_with_extracted_atoms", "value": len({str(a["event_id"]) for a in atoms.values()})},
            {"metric": "atomic_regions", "value": len(atoms)},
            {"metric": "two_aligner_callable_sites", "value": len(two_aligner)},
            {"metric": "events_selected_for_PRANK", "value": len(prank_events)},
            {"metric": "PRANK_alignments", "value": len(prank_sites)},
            {"metric": "three_aligner_callable_sites", "value": len(three_aligner)},
            *[
                {"metric": key, "value": value}
                for key, value in sorted(failures.items())
            ],
        ],
    )
    print(
        f"atoms={len(atoms)} two_aligner_sites={len(two_aligner)} "
        f"prank_events={len(prank_events)} "
        f"three_aligner_sites={len(three_aligner)}"
    )


if __name__ == "__main__":
    main()
