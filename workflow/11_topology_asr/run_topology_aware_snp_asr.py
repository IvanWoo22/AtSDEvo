#!/usr/bin/env python3
"""Topology-aware duplication-node ASR on three-aligner local SNP sites."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
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


def read_fasta_names(path: Path) -> set[str]:
    with path.open() as handle:
        return {
            line[1:].split()[0] for line in handle if line.startswith(">")
        }


def nested(names: list[str]) -> str:
    value = names[0]
    for name in names[1:]:
        value = f"({value},{name})"
    return value


def labeled_tree(names: set[str]) -> str:
    def clade(role: str) -> str:
        ordered = [
            f"Atha_{role}",
            f"Alyrata_{role}",
            f"Bstricta_{role}",
            f"Dstrictus_{role}",
        ]
        return nested([name for name in ordered if name in names])

    value = f"({clade('P')},{clade('D')})DUP"
    for species in SPECIES:
        name = f"{species}_P0"
        if name in names:
            value = f"({name},{value})"
    return value + ";\n"


def parse_state(path: Path) -> dict[int, tuple[str, float]]:
    result = {}
    with path.open() as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            if row["Node"] != "DUP":
                continue
            state = row["State"]
            result[int(row["Site"])] = (state, float(row[f"p_{state}"]))
    return result


def classify(p: str, d: str, ancestor: str) -> str:
    if p == d == ancestor:
        return "invariant"
    if d == ancestor and p != ancestor:
        return "P_specific"
    if p == ancestor and d != ancestor:
        return "D_specific"
    if p == d and p != ancestor:
        return "shared_PD"
    return "tri_allelic_unresolved"


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
        / 2**trials,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument(
        "--iqtree",
        type=Path,
        default=Path("iqtree3"),
        help="IQ-TREE executable (default: resolve 'iqtree3' from PATH)",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--posterior", type=float, default=0.90)
    args = parser.parse_args()
    root = args.root
    candidate_path = (
        root
        / "snp_local_msa_boundary_candidates/"
        "three_aligner_local_callable_sites.tsv"
    )
    candidates = read_tsv(candidate_path)
    atom_ids = sorted({row["ASR_atom_id"] for row in candidates})
    alignment_dir = root / "snp_local_msa/PRANK_FIXED_TREE"
    asr_dir = root / "snp_topology_asr/iqtree"
    asr_dir.mkdir(parents=True, exist_ok=True)

    def run_atom(atom_id: str) -> tuple[str, str]:
        alignment = alignment_dir / f"{atom_id}.best.fas"
        tree = asr_dir / f"{atom_id}.tree.nwk"
        tree.write_text(labeled_tree(read_fasta_names(alignment)))
        for model in ("JC", "K2P"):
            prefix = asr_dir / f"{atom_id}.{model}"
            state = Path(f"{prefix}.state")
            if state.exists():
                continue
            result = subprocess.run(
                [
                    str(args.iqtree),
                    "-s",
                    str(alignment),
                    "-te",
                    str(tree),
                    "-m",
                    model,
                    "-asr",
                    "-redo",
                    "-pre",
                    str(prefix),
                    "-nt",
                    "1",
                    "-quiet",
                ],
                capture_output=True,
                text=True,
            )
            if result.returncode or not state.exists():
                return atom_id, f"FAIL_{model}"
        return atom_id, "PASS"

    statuses = {}
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(run_atom, atom): atom for atom in atom_ids}
        for future in as_completed(futures):
            atom, status = future.result()
            statuses[atom] = status

    states = {}
    for atom in atom_ids:
        if statuses[atom] != "PASS":
            continue
        for model in ("JC", "K2P"):
            states[(atom, model)] = parse_state(
                Path(f"{asr_dir / f'{atom}.{model}'}.state")
            )

    output_rows = []
    for row in candidates:
        atom = row["ASR_atom_id"]
        site = int(row["ASR_alignment_column_1based"])
        jc = states.get((atom, "JC"), {}).get(site)
        k2p = states.get((atom, "K2P"), {}).get(site)
        asr_pass = (
            jc is not None
            and k2p is not None
            and jc[0] == k2p[0]
            and jc[1] >= args.posterior
            and k2p[1] >= args.posterior
        )
        ancestor = jc[0] if asr_pass and jc else "NA"
        direct_pass = (
            bool(row["deeper_P0_species_support"])
            and row["deeper_P0_conflict"] == "NO"
        )
        output_rows.append(
            {
                **row,
                "direct_boundary_plus_deeper_status": (
                    "PASS" if direct_pass else "FAIL"
                ),
                "JC_DUP_state": jc[0] if jc else "NA",
                "JC_DUP_posterior": f"{jc[1]:.6f}" if jc else "NA",
                "K2P_DUP_state": k2p[0] if k2p else "NA",
                "K2P_DUP_posterior": f"{k2p[1]:.6f}" if k2p else "NA",
                "topology_ASR_status": "PASS" if asr_pass else "FAIL",
                "topology_ASR_ancestral_base": ancestor,
                "topology_ASR_class": (
                    classify(row["P_base"], row["D_base"], ancestor)
                    if asr_pass
                    else "NA"
                ),
                "ASR_changes_boundary_P0_base": (
                    "YES"
                    if asr_pass and ancestor != row["ancestral_base"]
                    else "NO"
                    if asr_pass
                    else "NA"
                ),
                "ASR_resolves_deeper_P0_conflict": (
                    "YES"
                    if asr_pass and row["deeper_P0_conflict"] == "YES"
                    else "NO"
                ),
            }
        )
    write_tsv(root / "snp_topology_asr/site_level_ASR.tsv", output_rows)

    summaries = []
    endpoint_specs = (
        (
            "direct_boundary_plus_concordant_deeper",
            lambda row: row["direct_boundary_plus_deeper_status"] == "PASS",
            "polarized_class",
        ),
        (
            "fixed_topology_JC_K2P_posterior_ge_0.90",
            lambda row: row["topology_ASR_status"] == "PASS",
            "topology_ASR_class",
        ),
    )
    for endpoint, predicate, class_field in endpoint_specs:
        selected = [row for row in output_rows if predicate(row)]
        by_event = defaultdict(list)
        for row in selected:
            by_event[row["event_id"]].append(row)
        for threshold in (20, 50, 100, 200, 500):
            retained = {
                event: rows
                for event, rows in by_event.items()
                if len(rows) >= threshold
            }
            counts = Counter(
                row[class_field]
                for rows in retained.values()
                for row in rows
            )
            d_greater = sum(
                sum(r[class_field] == "D_specific" for r in rows)
                > sum(r[class_field] == "P_specific" for r in rows)
                for rows in retained.values()
            )
            p_greater = sum(
                sum(r[class_field] == "P_specific" for r in rows)
                > sum(r[class_field] == "D_specific" for r in rows)
                for rows in retained.values()
            )
            ties = len(retained) - d_greater - p_greater
            summaries.append(
                {
                    "endpoint": endpoint,
                    "minimum_callable_sites": threshold,
                    "events": len(retained),
                    "time1": sum(
                        rows[0]["age_bin"] == "time1"
                        for rows in retained.values()
                    ),
                    "time2": sum(
                        rows[0]["age_bin"] == "time2"
                        for rows in retained.values()
                    ),
                    "time3": sum(
                        rows[0]["age_bin"] == "time3"
                        for rows in retained.values()
                    ),
                    "callable_sites": sum(map(len, retained.values())),
                    "P_specific_SNP": counts["P_specific"],
                    "D_specific_SNP": counts["D_specific"],
                    "D_to_P_ratio": (
                        f"{counts['D_specific'] / counts['P_specific']:.8f}"
                        if counts["P_specific"]
                        else "Inf"
                    ),
                    "events_D_greater": d_greater,
                    "events_P_greater": p_greater,
                    "events_tied": ties,
                    "event_sign_test_p": (
                        f"{exact_binomial(d_greater, d_greater + p_greater):.8g}"
                        if d_greater + p_greater
                        else "NA"
                    ),
                }
            )
    write_tsv(root / "snp_topology_asr/endpoint_summary.tsv", summaries)
    posterior_rows = []
    for posterior in (0.80, 0.90, 0.95, 0.99):
        selected = [
            row
            for row in output_rows
            if row["JC_DUP_state"] in "ACGT"
            and row["JC_DUP_state"] == row["K2P_DUP_state"]
            and float(row["JC_DUP_posterior"]) >= posterior
            and float(row["K2P_DUP_posterior"]) >= posterior
        ]
        by_event = defaultdict(list)
        for row in selected:
            by_event[row["event_id"]].append(row)
        retained = {
            event: rows
            for event, rows in by_event.items()
            if len(rows) >= 200
        }
        classes = Counter(
            classify(row["P_base"], row["D_base"], row["JC_DUP_state"])
            for rows in retained.values()
            for row in rows
        )
        d_greater = p_greater = 0
        for rows in retained.values():
            counts = Counter(
                classify(row["P_base"], row["D_base"], row["JC_DUP_state"])
                for row in rows
            )
            d_greater += counts["D_specific"] > counts["P_specific"]
            p_greater += counts["P_specific"] > counts["D_specific"]
        posterior_rows.append(
            {
                "minimum_JC_and_K2P_DUP_posterior": posterior,
                "minimum_event_callable_sites": 200,
                "events": len(retained),
                "callable_sites": sum(map(len, retained.values())),
                "P_specific_SNP": classes["P_specific"],
                "D_specific_SNP": classes["D_specific"],
                "D_to_P_ratio": (
                    f"{classes['D_specific'] / classes['P_specific']:.8f}"
                    if classes["P_specific"]
                    else "Inf"
                ),
                "events_D_greater": d_greater,
                "events_P_greater": p_greater,
                "event_sign_test_p": (
                    f"{exact_binomial(d_greater, d_greater + p_greater):.8g}"
                    if d_greater + p_greater
                    else "NA"
                ),
            }
        )
    write_tsv(
        root / "snp_topology_asr/posterior_threshold_sensitivity.tsv",
        posterior_rows,
    )
    qc = [
        {"metric": "candidate_sites_boundary_rule", "value": len(candidates)},
        {"metric": "candidate_atoms", "value": len(atom_ids)},
        {
            "metric": "atoms_with_JC_and_K2P_ASR",
            "value": sum(value == "PASS" for value in statuses.values()),
        },
        {
            "metric": "direct_boundary_plus_deeper_sites",
            "value": sum(
                row["direct_boundary_plus_deeper_status"] == "PASS"
                for row in output_rows
            ),
        },
        {
            "metric": "topology_ASR_high_confidence_sites",
            "value": sum(
                row["topology_ASR_status"] == "PASS" for row in output_rows
            ),
        },
        {
            "metric": "ASR_resolved_deeper_P0_conflict_sites",
            "value": sum(
                row["ASR_resolves_deeper_P0_conflict"] == "YES"
                for row in output_rows
            ),
        },
        {
            "metric": "ASR_changed_boundary_P0_base_sites",
            "value": sum(
                row["ASR_changes_boundary_P0_base"] == "YES"
                for row in output_rows
            ),
        },
    ]
    write_tsv(root / "snp_topology_asr/workflow_qc.tsv", qc)
    print(
        f"sites={len(candidates)} atoms={len(atom_ids)} "
        f"ASR_pass={qc[4]['value']} conflict_rescued={qc[5]['value']}"
    )


if __name__ == "__main__":
    main()
