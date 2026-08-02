#!/usr/bin/env python3
"""Quantify controlled relaxations of P/D event admission rules."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
BOUNDARY = {"time1": "Alyrata", "time2": "Bstricta", "time3": "Dstrictus"}
POSTDUP = {
    "time1": (),
    "time2": ("Alyrata",),
    "time3": ("Alyrata", "Bstricta"),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def distinct(left: dict[str, str], right: dict[str, str]) -> bool:
    return (
        left["outgroup_scaffold"] != right["outgroup_scaffold"]
        or int(left["mapped_genomic_end"]) <= int(right["mapped_genomic_start"])
        or int(right["mapped_genomic_end"]) <= int(left["mapped_genomic_start"])
    )


def main() -> None:
    root = Path(__file__).resolve().parent
    pilot = root / "core_500"
    events = {
        row["event_id"]: row
        for row in read_tsv(pilot / "inputs/high_priority_core_eligible.tsv")
    }
    mappings = read_tsv(
        pilot / "outgroup_mapping/event_matched_blastn_mappings.tsv"
    )
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in mappings:
        grouped[
            (
                row["event_id"],
                row["species"],
                row["query_role"],
                row["target_copy_role"],
            )
        ].append(row)

    def best(
        event_id: str, species: str, query_role: str, target_role: str
    ) -> dict[str, str] | None:
        return max(
            grouped.get((event_id, species, query_role, target_role), []),
            key=lambda row: (
                row["two_sided_ordered_anchor_status"] == "PASS",
                float(row["query_coverage"]),
                float(row["nonoverlap_weighted_identity"]),
            ),
            default=None,
        )

    def species_pass(
        event_id: str, species: str, state: int, qcov: float, identity: float
    ) -> bool:
        required = [best(event_id, species, "P", "P")]
        if state == 3:
            required.append(best(event_id, species, "D", "D"))
        if any(
            row is None
            or row["two_sided_ordered_anchor_status"] != "PASS"
            or float(row["query_coverage"]) < qcov
            or float(row["nonoverlap_weighted_identity"]) < identity
            for row in required
        ):
            return False
        return state != 3 or distinct(required[0], required[1])  # type: ignore[arg-type]

    rows = []
    event_detail = []
    for rule, qcov, identity in (
        ("strict_mapping", 0.50, 0.60),
        ("moderate_mapping", 0.40, 0.58),
        ("exploratory_mapping", 0.35, 0.55),
    ):
        totals: Counter[tuple[str, str]] = Counter()
        for event_id, event in events.items():
            age = event["strict_age_bin"]
            states = {species: int(event[f"{species}_state"]) for species in SPECIES}
            status = {
                species: (
                    species_pass(event_id, species, states[species], qcov, identity)
                    if states[species] != 0
                    else False
                )
                for species in SPECIES
            }
            boundary = BOUNDARY[age]
            postdup = POSTDUP[age]
            mapped_postdup = sum(status[species] for species in postdup)
            strict_postdup = mapped_postdup == len(postdup)
            partial_postdup = (
                strict_postdup
                if len(postdup) <= 1
                else mapped_postdup >= len(postdup) - 1
            )
            boundary_pass = status[boundary]
            boundary_index = SPECIES.index(boundary)
            deeper_expected = [
                species
                for species in SPECIES[boundary_index + 1 :]
                if states[species] in (1, 2)
            ]
            deeper_pass = [species for species in deeper_expected if status[species]]
            robust_fallback_minimum = 2 if age == "time1" else 1
            robust_fallback = (
                not boundary_pass
                and len(deeper_pass) >= robust_fallback_minimum
            )
            categories = {
                "strict_boundary_strict_postdup": boundary_pass and strict_postdup,
                "strict_boundary_partial_postdup": boundary_pass and partial_postdup,
                "boundary_or_robust_deeper_P0_strict_postdup": (
                    (boundary_pass or robust_fallback) and strict_postdup
                ),
                "boundary_or_robust_deeper_P0_partial_postdup": (
                    (boundary_pass or robust_fallback) and partial_postdup
                ),
            }
            for category, passed in categories.items():
                if passed:
                    totals[(age, category)] += 1
                    totals[("all", category)] += 1
            event_detail.append(
                {
                    "mapping_rule": rule,
                    "event_id": event_id,
                    "age_bin": age,
                    "boundary_species": boundary,
                    "boundary_mapping": "PASS" if boundary_pass else "FAIL",
                    "postdup_expected": len(postdup),
                    "postdup_mapped": mapped_postdup,
                    "deeper_single_copy_expected": ",".join(deeper_expected),
                    "deeper_single_copy_mapped": ",".join(deeper_pass),
                    "robust_deeper_P0_fallback": (
                        "PASS" if robust_fallback else "FAIL"
                    ),
                    **{
                        category: "PASS" if passed else "FAIL"
                        for category, passed in categories.items()
                    },
                }
            )
        for age in ("time1", "time2", "time3", "all"):
            rows.append(
                {
                    "mapping_rule": rule,
                    "minimum_query_coverage": qcov,
                    "minimum_identity": identity,
                    "age_bin": age,
                    **{
                        category: totals[(age, category)]
                        for category in (
                            "strict_boundary_strict_postdup",
                            "strict_boundary_partial_postdup",
                            "boundary_or_robust_deeper_P0_strict_postdup",
                            "boundary_or_robust_deeper_P0_partial_postdup",
                        )
                    },
                }
            )
    write_tsv(root / "mapping_rule_grid.tsv", rows)
    write_tsv(root / "mapping_rule_event_detail.tsv", event_detail)

    strict_queue = read_tsv(pilot / "pilot/strict_primary_event_queue.tsv")
    qc = {
        row["event_id"]: row
        for row in read_tsv(pilot / "core/core_event_qc.tsv")
    }
    readiness = {
        row["event_id"]: row
        for row in read_tsv(
            pilot / "outgroup_mapping/event_mapping_readiness.tsv"
        )
    }
    species_summary = read_tsv(
        pilot / "outgroup_mapping/event_species_blastn_summary.tsv"
    )
    mapped_species: dict[str, list[str]] = defaultdict(list)
    for row in species_summary:
        if row["mapping_status"] == "PASS":
            mapped_species[row["event_id"]].append(row["species"])
    fallback_detail = [
        row
        for row in event_detail
        if row["mapping_rule"] == "strict_mapping"
        and row["robust_deeper_P0_fallback"] == "PASS"
        and row["boundary_or_robust_deeper_P0_strict_postdup"] == "PASS"
        and row["strict_boundary_strict_postdup"] == "FAIL"
    ]
    fallback_rows = []
    for detail in fallback_detail:
        event_id = str(detail["event_id"])
        event_qc = qc[event_id]
        deeper = str(detail["deeper_single_copy_mapped"]).split(",")
        deeper.sort(key=SPECIES.index)
        divergence_value = float(
            event_qc["present_day_PD_mismatch_pct_callable"]
        )
        callable_fraction = float(
            event_qc["jointly_callable_fraction_of_M"]
        )
        flags = ["deeper_P0_fallback"]
        if divergence_value > 40:
            flags.append("present_day_PD_mismatch_gt_40pct")
        if callable_fraction < 0.8:
            flags.append("joint_callable_fraction_lt_0.8")
        fallback_rows.append(
            {
                "run_order": 0,
                "event_id": event_id,
                "age_bin": detail["age_bin"],
                "analysis_tier": "DEEPER_P0_SENSITIVITY",
                "boundary_P0_species": deeper[0],
                "boundary_historical_corroborator": "NOT_AVAILABLE",
                "boundary_corroborator_exact_state_agreement": "NOT_AVAILABLE",
                "boundary_corroborator_evidence_scope": "NOT_AVAILABLE",
                "postduplication_PD_species": readiness[event_id][
                    "mapped_postduplication_species"
                ],
                "all_mapping_pass_species": ",".join(
                    sorted(mapped_species[event_id], key=SPECIES.index)
                ),
                "mapping_pass_species_count": len(mapped_species[event_id]),
                "paired_M_core_bp": event_qc["paired_M_core_bp"],
                "jointly_callable_uppercase_acgt_bp": event_qc[
                    "jointly_callable_uppercase_acgt_bp"
                ],
                "present_day_PD_mismatch_bp": event_qc[
                    "present_day_PD_mismatch_bp"
                ],
                "present_day_PD_mismatch_pct_callable": (
                    f"{divergence_value:.6f}"
                ),
                "jointly_callable_fraction_of_M": (
                    f"{callable_fraction:.6f}"
                ),
                "endpoint_sensitivity_flags": ",".join(flags),
                "BISER_I_bp": event_qc["mate1_insertion_I_bp"],
                "BISER_D_bp": event_qc["mate2_insertion_D_bp"],
                "BISER_softmasked_S_bp": event_qc["mate1_softmasked_S_bp"],
                "BISER_softmasked_N_bp": event_qc["mate2_softmasked_N_bp"],
                "primary_alignment": "PRANK_fixed_topology",
                "alignment_sensitivity": "MAFFT",
                "substitution_ASR_primary": "IQ-TREE3_fixed_topology",
                "substitution_ASR_sensitivity": "RAxML-NG",
                "indel_endpoint": "indelMaP",
                "conversion_filter": "GENECONV",
            }
        )
    combined_queue = strict_queue + fallback_rows
    combined_queue.sort(
        key=lambda row: (
            {"time1": 1, "time2": 2, "time3": 3}[row["age_bin"]],
            row["event_id"],
        )
    )
    for index, row in enumerate(combined_queue, 1):
        row["run_order"] = index
    write_tsv(root / "deeper_P0_sensitivity_queue.tsv", combined_queue)

    metrics = read_tsv(
        pilot / "sequence_variation/event_pd_substitution_metrics.tsv"
    )
    divergence = {
        row["event_id"]: float(row["present_day_PD_mismatch_pct_callable"])
        for row in read_tsv(pilot / "pilot/strict_primary_event_queue.tsv")
    }
    projection_rows = []
    for divergence_limit in (40.0, 50.0, 60.0, None):
        for callable_minimum in (100, 200, 300, 400, 500, 750, 1000):
            retained = [
                row
                for row in metrics
                if int(row["primary_bidir_callable_sites"]) >= callable_minimum
                and (
                    divergence_limit is None
                    or divergence[row["event_id"]] <= divergence_limit
                )
            ]
            projection_rows.append(
                {
                    "maximum_present_day_PD_mismatch_pct": (
                        "none" if divergence_limit is None else divergence_limit
                    ),
                    "minimum_bidirectional_exact_projection_bp": callable_minimum,
                    "retained_events": len(retained),
                    "time1": sum(row["age_bin"] == "time1" for row in retained),
                    "time2": sum(row["age_bin"] == "time2" for row in retained),
                    "time3": sum(row["age_bin"] == "time3" for row in retained),
                }
            )
    write_tsv(root / "bidirectional_projection_grid.tsv", projection_rows)


if __name__ == "__main__":
    main()
