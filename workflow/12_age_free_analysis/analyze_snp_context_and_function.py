#!/usr/bin/env python3
"""Build SBS96-like, substitution-subtype, and paired functional-context audits."""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import Counter
from pathlib import Path


ROOT = Path(os.environ.get("SD_PROJECT_ROOT", Path(__file__).resolve().parents[2]))
OUT = Path(
    os.environ.get(
        "SD_REASSESSMENT_ROOT", ROOT / "14_pd_polarization_reassessment"
    )
)
SNP_ROOT = ROOT / "12_inclusive_pd_sequence_variation"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
SIX = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
CONTEXTS = ("CDS", "exon_nonCDS", "intron_or_gene_body", "promoter_2kb", "intergenic")


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, object]]) -> None:
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


def fasta(path: Path) -> dict[str, str]:
    result: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                result[name] = []
            elif line:
                result[name].append(line.upper())
    return {key: "".join(value) for key, value in result.items()}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def normalize_context(trimer: str, ancestral: str, derived: str) -> tuple[str, str]:
    if ancestral in "AG":
        trimer = reverse_complement(trimer)
        ancestral = ancestral.translate(COMPLEMENT)
        derived = derived.translate(COMPLEMENT)
    change = f"{ancestral}>{derived}"
    return f"{trimer[0]}[{change}]{trimer[2]}", change


def cosine(left: list[float], right: list[float]) -> float:
    dot = sum(a * b for a, b in zip(left, right))
    norm = math.sqrt(sum(a * a for a in left) * sum(b * b for b in right))
    return dot / norm if norm else math.nan


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", type=Path, default=SNP_ROOT)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--minimum-sites", type=int, default=200)
    args = parser.parse_args()
    analysis_root = args.analysis_root
    output = args.output
    minimum_sites = args.minimum_sites
    events = {
        row["event_id"]
        for row in read(
            analysis_root
            / f"snp_local_msa/MAFFT_MUSCLE_PRANK_local_MSA.ge{minimum_sites}.event_metrics.tsv"
        )
    }
    direct = [
        row
        for row in read(
            analysis_root / "snp_local_msa/three_aligner_local_callable_sites.tsv"
        )
        if row["event_id"] in events
    ]
    topology_path = analysis_root / "snp_topology_asr/site_level_ASR.tsv"
    candidates = read(topology_path) if topology_path.exists() else []
    candidate_index = {}
    for row in candidates:
        key = (
            row["event_id"],
            row["P_chrom"],
            row["P_position_0based"],
            row["D_chrom"],
            row["D_position_0based"],
        )
        candidate_index[key] = row

    align_cache: dict[str, dict[str, str]] = {}
    opportunities = Counter()
    spectra = Counter()
    usable = missing_join = missing_context = 0
    site_rows = []
    for row in direct:
        key = (
            row["event_id"],
            row["P_chrom"],
            row["P_position_0based"],
            row["D_chrom"],
            row["D_position_0based"],
        )
        source = candidate_index.get(key, row if row.get("ASR_atom_id") else None)
        if source is None:
            missing_join += 1
            continue
        atom = source["ASR_atom_id"]
        if atom not in align_cache:
            align_cache[atom] = fasta(
                (
                    analysis_root
                    / "microindel_local_msa/alignments/PRANK_FIXED_TREE"
                    / f"{atom}.best.fas"
                )
                if (
                    analysis_root
                    / "microindel_local_msa/alignments/PRANK_FIXED_TREE"
                    / f"{atom}.best.fas"
                ).exists()
                else analysis_root
                / "snp_local_msa/PRANK_FIXED_TREE"
                / f"{atom}.best.fas"
            )
        alignment = align_cache[atom]
        boundary_name = f"{row['boundary_P0_species']}_P0"
        sequence = alignment.get(boundary_name)
        if sequence is None:
            missing_context += 1
            continue
        column = int(source["ASR_alignment_column_1based"]) - 1
        left = column - 1
        while left >= 0 and sequence[left] == "-":
            left -= 1
        right = column + 1
        while right < len(sequence) and sequence[right] == "-":
            right += 1
        if left < 0 or right >= len(sequence):
            missing_context += 1
            continue
        trimer = sequence[left] + sequence[column] + sequence[right]
        if any(base not in "ACGT" for base in trimer):
            missing_context += 1
            continue
        ancestral = row["ancestral_base"]
        if sequence[column] != ancestral:
            missing_context += 1
            continue
        normalized_trimer = (
            reverse_complement(trimer) if ancestral in "AG" else trimer
        )
        opportunities[normalized_trimer] += 1
        usable += 1
        if row["polarized_class"] not in {"P_specific", "D_specific"}:
            continue
        role = row["polarized_class"][0]
        derived = row[f"{role}_base"]
        label, change = normalize_context(trimer, ancestral, derived)
        spectra[(role, label)] += 1
        site_rows.append(
            {
                "event_id": row["event_id"],
                "copy_role": role,
                "ancestral_trinucleotide": normalized_trimer,
                "six_class_change": change,
                "SBS96_like_channel": label,
                "P_chrom": row["P_chrom"],
                "P_position_0based": row["P_position_0based"],
                "D_chrom": row["D_chrom"],
                "D_position_0based": row["D_position_0based"],
            }
        )
    write(output / "snp_context/SBS96_like_polarized_sites.tsv", site_rows)

    labels = [
        f"{left}[{change}]{right}"
        for change in SIX
        for left in "ACGT"
        for right in "ACGT"
    ]
    totals = {
        role: sum(spectra[(role, label)] for label in labels)
        for role in ("P", "D")
    }
    spectrum_rows = []
    for role in ("P", "D"):
        for label in labels:
            trimer = label[0] + label[2] + label[-1]
            count = spectra[(role, label)]
            opportunity = opportunities[trimer]
            spectrum_rows.append(
                {
                    "copy_role": role,
                    "SBS96_like_channel": label,
                    "SNP_count": count,
                    "fraction_within_role": (
                        f"{count / totals[role]:.8f}" if totals[role] else "NA"
                    ),
                    "ancestral_context_callable_sites": opportunity,
                    "SNP_per_context_callable_site": (
                        f"{count / opportunity:.8f}" if opportunity else "NA"
                    ),
                }
            )
    write(output / "snp_context/SBS96_like_spectrum.tsv", spectrum_rows)

    raw_p = [spectra[("P", label)] for label in labels]
    raw_d = [spectra[("D", label)] for label in labels]
    rate_p = [
        spectra[("P", label)] / opportunities[label[0] + label[2] + label[-1]]
        if opportunities[label[0] + label[2] + label[-1]]
        else 0
        for label in labels
    ]
    rate_d = [
        spectra[("D", label)] / opportunities[label[0] + label[2] + label[-1]]
        if opportunities[label[0] + label[2] + label[-1]]
        else 0
        for label in labels
    ]
    write(
        output / "statistics/SBS96_like_summary.tsv",
        [
            {
                "direct_callable_sites": len(direct),
                "callable_sites_with_ancestral_trinucleotide": usable,
                "missing_candidate_join": missing_join,
                "missing_ancestral_trinucleotide": missing_context,
                "P_SNP_with_trinucleotide": totals["P"],
                "D_SNP_with_trinucleotide": totals["D"],
                "P_D_raw_spectrum_cosine_similarity": f"{cosine(raw_p, raw_d):.8f}",
                "P_D_opportunity_normalized_cosine_similarity": f"{cosine(rate_p, rate_d):.8f}",
            }
        ],
    )

    subtype = Counter()
    for row in site_rows:
        subtype[(row["copy_role"], row["six_class_change"])] += 1
    subtype_rows = []
    for change in SIX:
        p = subtype[("P", change)]
        d = subtype[("D", change)]
        mutation_group = "transition" if change in {"C>T", "T>C"} else "transversion"
        subtype_rows.append(
            {
                "mutation_group": mutation_group,
                "six_class_change": change,
                "P_count": p,
                "D_count": d,
                "D_to_P_ratio": f"{d / p:.6f}" if p else "NA",
                "P_fraction_within_group": (
                    f"{p / sum(subtype[('P', c)] for c in SIX if ('transition' if c in {'C>T', 'T>C'} else 'transversion') == mutation_group):.8f}"
                ),
                "D_fraction_within_group": (
                    f"{d / sum(subtype[('D', c)] for c in SIX if ('transition' if c in {'C>T', 'T>C'} else 'transversion') == mutation_group):.8f}"
                ),
            }
        )
    write(output / "snp_context/transition_transversion_subtypes.tsv", subtype_rows)
    group_rows = []
    for role in ("P", "D"):
        ti = sum(subtype[(role, change)] for change in ("C>T", "T>C"))
        tv = sum(subtype[(role, change)] for change in ("C>A", "C>G", "T>A", "T>G"))
        group_rows.append(
            {
                "copy_role": role,
                "transitions": ti,
                "transversions": tv,
                "Ti_to_Tv_ratio": f"{ti / tv:.8f}" if tv else "NA",
                "C_to_T_over_T_to_C": (
                    f"{subtype[(role, 'C>T')] / subtype[(role, 'T>C')]:.8f}"
                    if subtype[(role, "T>C")]
                    else "NA"
                ),
            }
        )
    write(output / "snp_context/transition_transversion_role_summary.tsv", group_rows)

    positions = read(
        analysis_root
        / f"snp_types_ge{minimum_sites}/local_MSA_callable_genomic_positions.tsv"
    )
    if len(positions) % 2:
        raise AssertionError("Callable context rows are not P/D paired")
    context_matrix = Counter()
    same_context = Counter()
    p_intergenic_d_context = Counter()
    paired_rows = []
    for index in range(0, len(positions), 2):
        p, d = positions[index], positions[index + 1]
        if p["copy_role"] != "P" or d["copy_role"] != "D" or p["event_id"] != d["event_id"]:
            raise AssertionError(f"Broken P/D context pair at rows {index + 1}-{index + 2}")
        p_context, d_context = p["context"], d["context"]
        p_change = p["is_role_specific_change"] == "True"
        d_change = d["is_role_specific_change"] == "True"
        context_matrix[(p_context, d_context, "callable")] += 1
        context_matrix[(p_context, d_context, "P_change")] += p_change
        context_matrix[(p_context, d_context, "D_change")] += d_change
        if p_context == d_context:
            same_context[(p_context, "callable")] += 1
            same_context[(p_context, "P_change")] += p_change
            same_context[(p_context, "D_change")] += d_change
        if p_context == "intergenic":
            p_intergenic_d_context[(d_context, "callable")] += 1
            p_intergenic_d_context[(d_context, "P_change")] += p_change
            p_intergenic_d_context[(d_context, "D_change")] += d_change
        paired_rows.append(
            {
                "event_id": p["event_id"],
                "P_context": p_context,
                "D_context": d_context,
                "P_role_specific_change": p_change,
                "D_role_specific_change": d_change,
            }
        )
    write(output / "functional_context/paired_callable_contexts.tsv", paired_rows)
    write(
        output / "functional_context/PD_context_transition_matrix.tsv",
        [
            {
                "P_context": p_context,
                "D_context": d_context,
                "callable_pairs": context_matrix[(p_context, d_context, "callable")],
                "P_specific_SNP": context_matrix[(p_context, d_context, "P_change")],
                "D_specific_SNP": context_matrix[(p_context, d_context, "D_change")],
            }
            for p_context in CONTEXTS
            for d_context in CONTEXTS
            if context_matrix[(p_context, d_context, "callable")]
        ],
    )
    write(
        output / "functional_context/same_context_PD_rates.tsv",
        [
            {
                "shared_context": context,
                "callable_pairs": same_context[(context, "callable")],
                "P_specific_SNP": same_context[(context, "P_change")],
                "D_specific_SNP": same_context[(context, "D_change")],
                "P_rate": (
                    f"{same_context[(context, 'P_change')] / same_context[(context, 'callable')]:.8f}"
                    if same_context[(context, "callable")]
                    else "NA"
                ),
                "D_rate": (
                    f"{same_context[(context, 'D_change')] / same_context[(context, 'callable')]:.8f}"
                    if same_context[(context, "callable")]
                    else "NA"
                ),
            }
            for context in CONTEXTS
        ],
    )
    write(
        output / "functional_context/P_intergenic_corresponding_D_context.tsv",
        [
            {
                "D_context_at_P_intergenic_pair": context,
                "callable_pairs": p_intergenic_d_context[(context, "callable")],
                "P_specific_SNP": p_intergenic_d_context[(context, "P_change")],
                "D_specific_SNP": p_intergenic_d_context[(context, "D_change")],
            }
            for context in CONTEXTS
            if p_intergenic_d_context[(context, "callable")]
        ],
    )


if __name__ == "__main__":
    main()
