#!/usr/bin/env python3
"""Polarize TAIR12 P/D substitutions with syntenic pre-duplication loci.

The primary ancestral base requires the same anchored outgroup locus to be
aligned independently from both TAIR12 P and D queries, with concordant
projected bases. This prevents P-query ascertainment from automatically
inflating D-specific changes.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


SPECIES = ("Alyrata", "Bstricta", "Dstrictus", "Cviolacea")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


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


def read_fasta(path: Path) -> dict[str, str]:
    result: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                result[name] = []
            else:
                result[name].append(line.strip())
    return {name: "".join(parts) for name, parts in result.items()}


def intervals_coverage(intervals: list[tuple[int, int]]) -> int:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def parse_hsp(fields: list[str]) -> dict[str, object]:
    return {
        "query": fields[0],
        "target": fields[1],
        "qlen": int(fields[2]),
        "length": int(fields[3]),
        "pident": float(fields[4]),
        "qstart": int(fields[7]),
        "qend": int(fields[8]),
        "sstart": int(fields[9]),
        "send": int(fields[10]),
        "evalue": float(fields[11]),
        "bitscore": float(fields[12]),
        # NCBI BLAST+ omits the unsupported qstrand token; nucleotide query
        # coordinates are emitted in forward orientation. Field 14 is sstrand.
        "qstrand": "plus",
        "sstrand": fields[13],
        "qseq": fields[14],
        "sseq": fields[15],
    }


def effective_subject_strand(hsp: dict[str, object]) -> str:
    strand = str(hsp["sstrand"])
    if hsp["qstrand"] == "plus":
        return strand
    return "minus" if strand == "plus" else "plus"


def cluster_hsps(
    hsps: list[dict[str, object]], max_subject_gap: int
) -> list[list[dict[str, object]]]:
    clusters = []
    for strand in ("plus", "minus"):
        subset = [
            hsp for hsp in hsps
            if effective_subject_strand(hsp) == strand
        ]
        subset.sort(
            key=lambda hsp: min(int(hsp["sstart"]), int(hsp["send"]))
        )
        current: list[dict[str, object]] = []
        current_end = -1
        for hsp in subset:
            start = min(int(hsp["sstart"]), int(hsp["send"]))
            end = max(int(hsp["sstart"]), int(hsp["send"]))
            if current and start > current_end + max_subject_gap:
                clusters.append(current)
                current = []
            current.append(hsp)
            current_end = max(current_end, end)
        if current:
            clusters.append(current)
    return clusters


def cluster_score(cluster: list[dict[str, object]]) -> tuple[int, float]:
    coverage = intervals_coverage(
        [
            (
                min(int(hsp["qstart"]), int(hsp["qend"])) - 1,
                max(int(hsp["qstart"]), int(hsp["qend"])),
            )
            for hsp in cluster
        ]
    )
    return coverage, sum(float(hsp["bitscore"]) for hsp in cluster)


def subject_span(cluster: list[dict[str, object]]) -> tuple[int, int]:
    return (
        min(min(int(hsp["sstart"]), int(hsp["send"])) for hsp in cluster),
        max(max(int(hsp["sstart"]), int(hsp["send"])) for hsp in cluster),
    )


def span_distance(a: tuple[int, int], b: tuple[int, int]) -> int:
    if a[1] < b[0]:
        return b[0] - a[1]
    if b[1] < a[0]:
        return a[0] - b[1]
    return 0


def project_cluster(
    cluster: list[dict[str, object]], query_sequence: str
) -> tuple[dict[int, tuple[str, int]], int]:
    projected: dict[int, tuple[float, str, int]] = {}
    query_validation_failures = 0
    for hsp in sorted(
        cluster, key=lambda row: float(row["bitscore"]), reverse=True
    ):
        if hsp["qstrand"] == "plus":
            position = int(hsp["qstart"]) - 1
            step = 1
            normalize = lambda base: base
        else:
            position = int(hsp["qstart"]) - 1
            step = -1
            normalize = lambda base: base.translate(COMPLEMENT)
        subject_position = int(hsp["sstart"]) - 1
        subject_step = (
            1 if int(hsp["send"]) >= int(hsp["sstart"]) else -1
        )
        for query_base, subject_base in zip(
            str(hsp["qseq"]), str(hsp["sseq"])
        ):
            if query_base == "-":
                if subject_base != "-":
                    subject_position += subject_step
                continue
            normalized_query = normalize(query_base).upper()
            if (
                0 <= position < len(query_sequence)
                and normalized_query in "ACGT"
                and query_sequence[position].upper() in "ACGT"
                and normalized_query != query_sequence[position].upper()
            ):
                query_validation_failures += 1
            if subject_base != "-":
                base = normalize(subject_base).upper()
                if base in "ACGT" and position not in projected:
                    projected[position] = (
                        float(hsp["bitscore"]),
                        base,
                        subject_position,
                    )
                subject_position += subject_step
            position += step
    return {
        position: (value[1], value[2])
        for position, value in projected.items()
    }, query_validation_failures


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


def metric_fields(prefix: str, counts: Counter[str]) -> dict[str, object]:
    callable_sites = sum(counts.values())
    polarized = counts["P_specific"] + counts["D_specific"]
    return {
        f"{prefix}_callable_sites": callable_sites,
        f"{prefix}_invariant": counts["invariant"],
        f"{prefix}_P_specific": counts["P_specific"],
        f"{prefix}_D_specific": counts["D_specific"],
        f"{prefix}_shared_PD": counts["shared_PD"],
        f"{prefix}_tri_allelic_unresolved": counts[
            "tri_allelic_unresolved"
        ],
        f"{prefix}_polarized_changes": polarized,
        f"{prefix}_P_specific_rate": (
            f"{counts['P_specific'] / callable_sites:.8f}"
            if callable_sites else "NA"
        ),
        f"{prefix}_D_specific_rate": (
            f"{counts['D_specific'] / callable_sites:.8f}"
            if callable_sites else "NA"
        ),
        f"{prefix}_D_minus_P_rate": (
            f"{(counts['D_specific'] - counts['P_specific']) / callable_sites:.8f}"
            if callable_sites else "NA"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", required=True, type=Path)
    parser.add_argument("--max-subject-gap", type=int, default=5000)
    parser.add_argument("--max-PD-locus-distance", type=int, default=5000)
    parser.add_argument("--queue", type=Path)
    parser.add_argument("--output-directory", type=Path)
    args = parser.parse_args()

    queue_path = (
        args.queue or args.pilot / "pilot/strict_primary_event_queue.tsv"
    )
    output_directory = (
        args.output_directory or args.pilot / "sequence_variation"
    )
    queue = read_tsv(queue_path)
    events = {
        row["event_id"]: row
        for row in read_tsv(
            args.pilot / "inputs/high_priority_core_eligible.tsv"
        )
    }
    mapping = {
        (row["event_id"], row["species"]): row
        for row in read_tsv(
            args.pilot / "outgroup_mapping/event_species_blastn_summary.tsv"
        )
    }
    fasta = read_fasta(
        args.pilot / "core/TAIR12_PD_homologous_cores.fa"
    )
    strict_ids = {row["event_id"] for row in queue}

    hits: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for species in SPECIES:
        path = (
            args.pilot
            / f"outgroup_mapping/blastn_aligned/{species}.aligned_hits.tsv"
        )
        with path.open() as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                event_id = fields[0].split("|", 1)[0]
                if event_id not in strict_ids:
                    continue
                role = fields[0].split("|")[1]
                target_event = fields[1].split("|", 1)[0]
                if event_id != target_event:
                    continue
                hits[(event_id, species, role, fields[1])].append(
                    parse_hsp(fields)
                )

    event_rows = []
    species_rows = []
    site_rows = []
    postdup_rows = []
    for queued in queue:
        event_id = queued["event_id"]
        event = events[event_id]
        age = queued["age_bin"]
        # Preserve FASTA case: lowercase bases are soft-masked and must not
        # enter the jointly callable substitution denominator.
        p_sequence = fasta[f"{event_id}|P|{age}"]
        d_sequence = fasta[f"{event_id}|D|{age}"]
        if len(p_sequence) != len(d_sequence):
            raise AssertionError(f"Unequal P/D core length: {event_id}")

        species_projections: dict[str, dict[str, object]] = {}
        for species in SPECIES:
            species_map = mapping.get((event_id, species))
            if (
                not species_map
                or species_map["expected_locus_class"]
                != "preduplication_single_copy"
                or species_map["mapping_status"] != "PASS"
            ):
                continue
            target = species_map["P_to_P_candidate_id"]
            p_clusters = cluster_hsps(
                hits.get((event_id, species, "P", target), []),
                args.max_subject_gap,
            )
            if not p_clusters:
                continue
            p_cluster = max(p_clusters, key=cluster_score)
            p_span = subject_span(p_cluster)
            p_strand = effective_subject_strand(p_cluster[0])
            d_clusters = [
                cluster
                for cluster in cluster_hsps(
                    hits.get((event_id, species, "D", target), []),
                    args.max_subject_gap,
                )
                if effective_subject_strand(cluster[0]) == p_strand
            ]
            d_cluster = min(
                d_clusters,
                key=lambda cluster: (
                    span_distance(p_span, subject_span(cluster)),
                    -cluster_score(cluster)[0],
                    -cluster_score(cluster)[1],
                ),
                default=None,
            )
            d_distance = (
                span_distance(p_span, subject_span(d_cluster))
                if d_cluster else -1
            )
            p_projection, p_fail = project_cluster(p_cluster, p_sequence)
            d_projection: dict[int, tuple[str, int]] = {}
            d_fail = 0
            if d_cluster and d_distance <= args.max_PD_locus_distance:
                d_projection, d_fail = project_cluster(
                    d_cluster, d_sequence
                )
            bidirectional = {
                position: value[0]
                for position, value in p_projection.items()
                if d_projection.get(position) == value
            }
            species_projections[species] = {
                "P_query": {
                    position: value[0]
                    for position, value in p_projection.items()
                },
                "D_query": {
                    position: value[0]
                    for position, value in d_projection.items()
                },
                "bidirectional": bidirectional,
            }
            species_rows.append(
                {
                    "event_id": event_id,
                    "age_bin": age,
                    "species": species,
                    "candidate_id": target,
                    "P_cluster_query_coverage_bp": cluster_score(p_cluster)[0],
                    "D_cluster_query_coverage_bp": (
                        cluster_score(d_cluster)[0] if d_cluster else 0
                    ),
                    "PD_subject_locus_distance_bp": d_distance,
                    "P_query_projected_bp": len(p_projection),
                    "D_query_projected_bp": len(d_projection),
                    "bidirectional_concordant_projected_bp": len(
                        bidirectional
                    ),
                    "P_query_validation_failures": p_fail,
                    "D_query_validation_failures": d_fail,
                }
            )

        # Independent equal-time validation for post-duplication species:
        # compare A. thaliana P to outgroup P and A. thaliana D to outgroup D.
        # Both loci are selected symmetrically and are already required to be
        # distinct, two-sided-anchor-supported loci by the strict mapping rule.
        for species in SPECIES:
            species_map = mapping.get((event_id, species))
            if (
                not species_map
                or species_map["expected_locus_class"]
                != "postduplication_two_copy"
                or species_map["mapping_status"] != "PASS"
            ):
                continue
            role_projection = {}
            role_cluster_coverage = {}
            for role in ("P", "D"):
                target = species_map[f"{role}_to_{role}_candidate_id"]
                clusters = cluster_hsps(
                    hits.get((event_id, species, role, target), []),
                    args.max_subject_gap,
                )
                if not clusters:
                    role_projection[role] = {}
                    role_cluster_coverage[role] = 0
                    continue
                cluster = max(clusters, key=cluster_score)
                projected, validation_failures = project_cluster(
                    cluster,
                    p_sequence if role == "P" else d_sequence,
                )
                if validation_failures:
                    raise AssertionError(
                        f"Query validation failed: {event_id}/{species}/{role}"
                    )
                role_projection[role] = {
                    position: value[0]
                    for position, value in projected.items()
                }
                role_cluster_coverage[role] = cluster_score(cluster)[0]
            joint_positions = [
                position
                for position in range(len(p_sequence))
                if p_sequence[position] in "ACGT"
                and d_sequence[position] in "ACGT"
                and position in role_projection["P"]
                and position in role_projection["D"]
            ]
            p_mismatches = sum(
                p_sequence[position] != role_projection["P"][position]
                for position in joint_positions
            )
            d_mismatches = sum(
                d_sequence[position] != role_projection["D"][position]
                for position in joint_positions
            )
            postdup_rows.append(
                {
                    "event_id": event_id,
                    "age_bin": age,
                    "postduplication_species": species,
                    "joint_callable_sites": len(joint_positions),
                    "P_terminal_mismatches": p_mismatches,
                    "D_terminal_mismatches": d_mismatches,
                    "P_terminal_mismatch_rate": (
                        f"{p_mismatches / len(joint_positions):.8f}"
                        if joint_positions else "NA"
                    ),
                    "D_terminal_mismatch_rate": (
                        f"{d_mismatches / len(joint_positions):.8f}"
                        if joint_positions else "NA"
                    ),
                    "D_minus_P_terminal_rate": (
                        f"{(d_mismatches - p_mismatches) / len(joint_positions):.8f}"
                        if joint_positions else "NA"
                    ),
                    "P_cluster_query_coverage_bp": role_cluster_coverage["P"],
                    "D_cluster_query_coverage_bp": role_cluster_coverage["D"],
                    "endpoint_sensitivity_flags": queued[
                        "endpoint_sensitivity_flags"
                    ],
                }
            )

        boundary = queued["boundary_P0_species"]
        boundary_projection = species_projections.get(boundary, {})
        boundary_bidir: dict[int, str] = boundary_projection.get(
            "bidirectional", {}
        )
        boundary_pquery: dict[int, str] = boundary_projection.get(
            "P_query", {}
        )
        primary_counts: Counter[str] = Counter()
        pquery_counts: Counter[str] = Counter()
        multi_counts: Counter[str] = Counter()
        single_species = sorted(species_projections)
        for position, (p_base, d_base) in enumerate(
            zip(p_sequence, d_sequence)
        ):
            if p_base not in "ACGT" or d_base not in "ACGT":
                continue
            primary_class = "NA"
            pquery_class = "NA"
            multi_class = "NA"
            if position in boundary_bidir:
                primary_class = classify(
                    p_base, d_base, boundary_bidir[position]
                )
                primary_counts[primary_class] += 1
            if position in boundary_pquery:
                pquery_class = classify(
                    p_base, d_base, boundary_pquery[position]
                )
                pquery_counts[pquery_class] += 1
            bases = [
                species_projections[species]["bidirectional"].get(position)
                for species in single_species
            ]
            bases = [
                base
                for base in bases
                if isinstance(base, str) and base in "ACGT"
            ]
            if len(bases) >= 2 and len(set(bases)) == 1:
                multi_class = classify(p_base, d_base, bases[0])
                multi_counts[multi_class] += 1
            if primary_class != "NA":
                site_rows.append(
                    {
                        "event_id": event_id,
                        "age_bin": age,
                        "core_position_1based": position + 1,
                        "P_base": p_base,
                        "D_base": d_base,
                        "boundary_species": boundary,
                        "boundary_bidirectional_ancestral_base": (
                            boundary_bidir[position]
                        ),
                        "primary_class": primary_class,
                        "boundary_Pquery_class": pquery_class,
                        "multi_outgroup_concordant_class": multi_class,
                    }
                )

        event_rows.append(
            {
                "run_order": queued["run_order"],
                "event_id": event_id,
                "age_bin": age,
                "boundary_species": boundary,
                "strict_core_bp": len(p_sequence),
                "single_copy_outgroups_projected": ",".join(single_species),
                "single_copy_outgroup_count": len(single_species),
                "endpoint_sensitivity_flags": queued[
                    "endpoint_sensitivity_flags"
                ],
                "boundary_corroborator_exact_state_agreement": queued[
                    "boundary_corroborator_exact_state_agreement"
                ],
                **metric_fields("primary_bidir", primary_counts),
                **metric_fields("boundary_Pquery", pquery_counts),
                **metric_fields("multi_outgroup", multi_counts),
            }
        )

    write_tsv(
        output_directory / "outgroup_projection_qc.tsv",
        species_rows,
    )
    write_tsv(
        output_directory / "event_pd_substitution_metrics.tsv",
        event_rows,
    )
    write_tsv(
        output_directory / "polarized_sites.tsv",
        site_rows,
    )
    write_tsv(
        output_directory / "postduplication_terminal_branch_metrics.tsv",
        postdup_rows,
    )
    print(
        f"Polarized {len(event_rows)} events and wrote "
        f"{len(site_rows)} primary callable sites"
    )


if __name__ == "__main__":
    main()
