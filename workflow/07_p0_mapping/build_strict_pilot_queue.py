#!/usr/bin/env python3
"""Build the strict time1-time3 P/D pilot queue and analysis design."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", required=True, type=Path)
    parser.add_argument(
        "--mapping-rule",
        choices=("strict", "relaxed"),
        default="strict",
    )
    parser.add_argument("--queue-output", type=Path)
    parser.add_argument("--summary-output", type=Path)
    args = parser.parse_args()

    qc = {
        row["event_id"]: row
        for row in read_tsv(args.pilot / "core/core_event_qc.tsv")
    }
    readiness = read_tsv(
        args.pilot / "outgroup_mapping/event_mapping_readiness.tsv"
    )
    mappings = read_tsv(
        args.pilot / "outgroup_mapping/event_species_blastn_summary.tsv"
    )
    by_event: dict[str, list[dict[str, str]]] = {}
    for row in mappings:
        by_event.setdefault(row["event_id"], []).append(row)
    corroborator_path = (
        args.pilot
        / "outgroup_mapping/historical_corroborator_synteny.tsv"
    )
    corroborators = {}
    if corroborator_path.exists():
        corroborators = {
            (row["event_id"], row["node"]): row
            for row in read_tsv(corroborator_path)
        }

    queue: list[dict[str, object]] = []
    for ready in readiness:
        readiness_field = f"{args.mapping_rule}_mapping_ready"
        if ready[readiness_field] != "PASS":
            continue
        event_id = ready["event_id"]
        event_qc = qc[event_id]
        mapped = [
            row for row in by_event.get(event_id, [])
            if row["mapping_status"] == "PASS"
        ]
        mapped.sort(key=lambda row: row["species"])
        boundary = ready["boundary_species"]
        postdup = ready["mapped_postduplication_species"]
        divergence = float(event_qc["present_day_PD_mismatch_pct_callable"])
        callable_bp = int(event_qc["jointly_callable_uppercase_acgt_bp"])
        callable_fraction = float(
            event_qc["jointly_callable_fraction_of_M"]
        )
        sensitivity_flags = []
        if divergence > 40:
            sensitivity_flags.append("present_day_PD_mismatch_gt_40pct")
        if callable_fraction < 0.8:
            sensitivity_flags.append("joint_callable_fraction_lt_0.8")
        boundary_node = {
            "time1": "N1", "time2": "N2", "time3": "N3"
        }[ready["age_bin"]]
        corroborator = corroborators.get((event_id, boundary_node), {})
        # Ranking is deterministic and intentionally favors information content:
        # more independently mapped species, then longer callable core, then lower
        # present-day P/D divergence (less alignment uncertainty).
        queue.append(
            {
                "event_id": event_id,
                "age_bin": ready["age_bin"],
                "analysis_tier": (
                    "STRICT_PRIMARY"
                    if args.mapping_rule == "strict"
                    else "PARTIAL_POSTDUP_SENSITIVITY"
                ),
                "boundary_P0_species": boundary,
                "boundary_historical_corroborator": corroborator.get(
                    "corroborator_species", "NOT_AVAILABLE"
                ),
                "boundary_corroborator_exact_state_agreement": (
                    corroborator.get("exact_state_agreement", "NOT_AVAILABLE")
                ),
                "boundary_corroborator_evidence_scope": (
                    "AUXILIARY_NOT_SEQUENCE_GATING"
                    if corroborator
                    else "NOT_AVAILABLE"
                ),
                "postduplication_PD_species": postdup,
                "all_mapping_pass_species": ",".join(
                    row["species"] for row in mapped
                ),
                "mapping_pass_species_count": len(mapped),
                "paired_M_core_bp": event_qc["paired_M_core_bp"],
                "jointly_callable_uppercase_acgt_bp": callable_bp,
                "present_day_PD_mismatch_bp": event_qc[
                    "present_day_PD_mismatch_bp"
                ],
                "present_day_PD_mismatch_pct_callable": (
                    f"{divergence:.6f}"
                ),
                "jointly_callable_fraction_of_M": (
                    f"{callable_fraction:.6f}"
                ),
                "endpoint_sensitivity_flags": (
                    ",".join(sensitivity_flags)
                    if sensitivity_flags
                    else "NONE"
                ),
                "BISER_I_bp": event_qc["mate1_insertion_I_bp"],
                "BISER_D_bp": event_qc["mate2_insertion_D_bp"],
                "BISER_softmasked_S_bp": event_qc[
                    "mate1_softmasked_S_bp"
                ],
                "BISER_softmasked_N_bp": event_qc[
                    "mate2_softmasked_N_bp"
                ],
                "primary_alignment": "PRANK_fixed_topology",
                "alignment_sensitivity": "MAFFT",
                "substitution_ASR_primary": "IQ-TREE3_fixed_topology",
                "substitution_ASR_sensitivity": "RAxML-NG",
                "indel_endpoint": "indelMaP",
                "conversion_filter": "GENECONV",
                "_sort": (
                    {"time1": 1, "time2": 2, "time3": 3}[ready["age_bin"]],
                    0
                    if corroborator.get("exact_state_agreement") == "PASS"
                    else 1,
                    0 if not sensitivity_flags else 1,
                    -len(mapped),
                    -callable_bp,
                    divergence,
                    event_id,
                ),
            }
        )

    queue.sort(key=lambda row: row.pop("_sort"))
    for index, row in enumerate(queue, 1):
        row["run_order"] = index
    # Put run order first without losing deterministic column order.
    queue = [{"run_order": row.pop("run_order"), **row} for row in queue]
    queue_output = (
        args.queue_output
        or args.pilot / "pilot/strict_primary_event_queue.tsv"
    )
    write_tsv(queue_output, queue)

    summary = []
    counts = Counter(row["age_bin"] for row in queue)
    for age in ("time1", "time2", "time3"):
        subset = [row for row in queue if row["age_bin"] == age]
        summary.append(
            {
                "age_bin": age,
                "strict_primary_events": counts[age],
                "median_callable_core_bp": (
                    sorted(
                        int(row["jointly_callable_uppercase_acgt_bp"])
                        for row in subset
                    )[len(subset) // 2]
                    if subset
                    else 0
                ),
                "median_present_day_PD_mismatch_pct": (
                    f"{sorted(float(row['present_day_PD_mismatch_pct_callable']) for row in subset)[len(subset) // 2]:.6f}"
                    if subset
                    else "NA"
                ),
            }
        )
    summary_output = (
        args.summary_output
        or args.pilot / "pilot/strict_primary_summary.tsv"
    )
    write_tsv(summary_output, summary)
    print(f"Queued {len(queue)} events under {args.mapping_rule} mapping rule")


if __name__ == "__main__":
    main()
