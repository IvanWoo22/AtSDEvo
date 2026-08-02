#!/usr/bin/env python3
"""Add historical near-species synteny as auxiliary, non-gating evidence."""

from __future__ import annotations

import argparse
import csv
import importlib.util
from collections import Counter
from pathlib import Path


CORROBORATORS = {
    "Ahalleri": ("Alyrata", "N1"),
    "Rislandica": ("Bstricta", "N2"),
    "Esyriacum": ("Dstrictus", "N3"),
}


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


def load_projection_helpers(project: Path):
    source = (
        project
        / "06_sd_age_tracing_preparation/workflow_scripts/"
        "evaluate_primary_node_sd.py"
    )
    spec = importlib.util.spec_from_file_location("primary_eval", source)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--pilot", required=True, type=Path)
    parser.add_argument("--historical-root", required=True, type=Path)
    parser.add_argument("--flank", type=int, default=2000)
    parser.add_argument("--overlap-threshold", type=float, default=0.5)
    args = parser.parse_args()

    helper = load_projection_helpers(args.project)
    atha_genes = helper.read_atha_genes(
        args.project / "05_mcscanx_synteny/inputs/Atha.gff"
    )
    events = {
        row["event_id"]: row
        for row in read_tsv(
            args.pilot / "inputs/high_priority_core_eligible.tsv"
        )
    }
    queue_ids = {
        row["event_id"]
        for row in read_tsv(
            args.pilot / "pilot/strict_primary_event_queue.tsv"
        )
    }

    projections = {}
    for species in CORROBORATORS:
        collinearity = (
            args.historical_root / f"Atha_{species}.collinearity"
        )
        merged, _, _ = helper.project_blocks(
            collinearity, atha_genes, args.flank
        )
        projections[species] = merged

    output = []
    for event_id in sorted(queue_ids):
        event = events[event_id]
        copy_intervals = [
            (
                event["representative_copy1_chrom"],
                int(event["copy1_envelope_start"]),
                int(event["copy1_envelope_end"]),
            ),
            (
                event["representative_copy2_chrom"],
                int(event["copy2_envelope_start"]),
                int(event["copy2_envelope_end"]),
            ),
        ]
        for species, (primary, node) in CORROBORATORS.items():
            fractions = []
            for chrom, start, end in copy_intervals:
                overlap = helper.interval_overlap(
                    start, end, projections[species].get(chrom, [])
                )
                fractions.append(overlap / (end - start))
            state = (
                (1 if fractions[0] >= args.overlap_threshold else 0)
                + (2 if fractions[1] >= args.overlap_threshold else 0)
            )
            primary_state = int(event[f"{primary}_state"])
            output.append(
                {
                    "event_id": event_id,
                    "age_bin": event["strict_age_bin"],
                    "node": node,
                    "primary_species": primary,
                    "primary_state": primary_state,
                    "corroborator_species": species,
                    "corroborator_state": state,
                    "copy1_overlap_fraction": f"{fractions[0]:.6f}",
                    "copy2_overlap_fraction": f"{fractions[1]:.6f}",
                    "exact_state_agreement": (
                        "PASS" if state == primary_state else "FAIL"
                    ),
                    "opposite_single_copy_conflict": (
                        "YES"
                        if {state, primary_state} == {1, 2}
                        else "NO"
                    ),
                    "evidence_scope": (
                        "AUXILIARY_HISTORICAL_MCScanX_NOT_SEQUENCE_GATING"
                    ),
                }
            )
    write_tsv(
        args.pilot
        / "outgroup_mapping/historical_corroborator_synteny.tsv",
        output,
    )

    summary = []
    for species, (primary, node) in CORROBORATORS.items():
        subset = [row for row in output if row["corroborator_species"] == species]
        states = Counter(int(row["corroborator_state"]) for row in subset)
        summary.append(
            {
                "node": node,
                "primary_species": primary,
                "corroborator_species": species,
                "strict_queue_events": len(subset),
                "exact_state_agreement": sum(
                    row["exact_state_agreement"] == "PASS" for row in subset
                ),
                "opposite_single_copy_conflicts": sum(
                    row["opposite_single_copy_conflict"] == "YES"
                    for row in subset
                ),
                "corroborator_state0": states[0],
                "corroborator_state1": states[1],
                "corroborator_state2": states[2],
                "corroborator_state3": states[3],
            }
        )
    write_tsv(
        args.pilot
        / "outgroup_mapping/historical_corroborator_synteny_summary.tsv",
        summary,
    )
    print(f"Added {len(output)} auxiliary corroborator observations")


if __name__ == "__main__":
    main()
