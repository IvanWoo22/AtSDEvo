#!/usr/bin/env python3
"""Classify polarized SNP spectra and genomic contexts for the controlled 55 events."""

from __future__ import annotations

import argparse
import bisect
import csv
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats


COMPLEMENT = str.maketrans("ACGT", "TGCA")
GENETIC_CODE = {
    codon: aa
    for aa, codons in {
        "F": ("TTT", "TTC"), "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
        "I": ("ATT", "ATC", "ATA"), "M": ("ATG",), "V": ("GTT", "GTC", "GTA", "GTG"),
        "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
        "P": ("CCT", "CCC", "CCA", "CCG"), "T": ("ACT", "ACC", "ACA", "ACG"),
        "A": ("GCT", "GCC", "GCA", "GCG"), "Y": ("TAT", "TAC"),
        "*": ("TAA", "TAG", "TGA"), "H": ("CAT", "CAC"), "Q": ("CAA", "CAG"),
        "N": ("AAT", "AAC"), "K": ("AAA", "AAG"), "D": ("GAT", "GAC"),
        "E": ("GAA", "GAG"), "C": ("TGT", "TGC"), "W": ("TGG",),
        "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
        "G": ("GGT", "GGC", "GGA", "GGG"),
    }.items()
    for codon in codons
}
SIX_CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
EFFECT_PRIORITY = {
    "stop_gain": 5,
    "stop_loss": 5,
    "start_loss": 4,
    "missense": 3,
    "synonymous": 2,
    "CDS_unresolved": 1,
    "noncoding": 0,
}


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
    sequences: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            else:
                sequences[name].append(line.strip())
    return {name: "".join(parts) for name, parts in sequences.items()}


def parse_attributes(text: str) -> dict[str, str]:
    result = {}
    for field in text.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
    return result


def normalize_change(ancestral: str, derived: str) -> str:
    if ancestral in "AG":
        ancestral = ancestral.translate(COMPLEMENT)
        derived = derived.translate(COMPLEMENT)
    return f"{ancestral}>{derived}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    project = args.project
    expansion = project / "08_event_inclusion_sensitivity"
    pilot = expansion / "core_500"
    controlled = read_tsv(expansion / "controlled_expansion_endpoint_events.tsv")
    event_tier = {row["event_id"]: row["admission_tier"] for row in controlled}
    event_callable = {
        row["event_id"]: int(row["primary_bidir_callable_sites"])
        for row in controlled
    }
    sources = {
        "core_500_strict": pilot / "sequence_variation/polarized_sites.tsv",
        "partial_postduplication": (
            pilot / "sequence_variation_partial_postdup/polarized_sites.tsv"
        ),
        "deeper_P0_fallback": (
            pilot / "sequence_variation_deeper_P0/polarized_sites.tsv"
        ),
    }
    polarized = []
    for tier, path in sources.items():
        wanted = {event for event, value in event_tier.items() if value == tier}
        for row in read_tsv(path):
            if row["event_id"] in wanted and row["primary_class"] in {
                "P_specific", "D_specific"
            }:
                polarized.append(row)

    blocks: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(pilot / "core/homologous_core_blocks.tsv"):
        if row["event_id"] in event_tier:
            blocks[row["event_id"]].append(row)
    for event_rows in blocks.values():
        event_rows.sort(key=lambda row: int(row["core_block_index"]))

    sites: list[dict[str, object]] = []
    for number, row in enumerate(polarized, 1):
        role = row["primary_class"][0]
        position = int(row["core_position_1based"]) - 1
        cumulative = 0
        selected = None
        offset = -1
        for block in blocks[row["event_id"]]:
            length = int(block["core_block_bp"])
            if cumulative <= position < cumulative + length:
                selected = block
                offset = position - cumulative
                break
            cumulative += length
        if selected is None:
            raise AssertionError(f"Core coordinate not found: {row['event_id']}/{position}")
        copy = "copy1" if selected["copy1_role"] == role else "copy2"
        strand = selected[f"{copy}_strand"]
        start = int(selected[f"{copy}_start"])
        end = int(selected[f"{copy}_end"])
        genomic_position = start + offset if strand == "+" else end - 1 - offset
        derived = row[f"{role}_base"]
        ancestral = row["boundary_bidirectional_ancestral_base"]
        derived_genomic = (
            derived if strand == "+" else derived.translate(COMPLEMENT)
        )
        ancestral_genomic = (
            ancestral if strand == "+" else ancestral.translate(COMPLEMENT)
        )
        change = normalize_change(ancestral, derived)
        sites.append(
            {
                "site_id": f"SDVS{number:07d}",
                "event_id": row["event_id"],
                "age_bin": row["age_bin"],
                "admission_tier": event_tier[row["event_id"]],
                "copy_role": role,
                "core_position_1based": position + 1,
                "chrom": selected[f"{copy}_chrom"],
                "genomic_position_1based": genomic_position + 1,
                "alignment_strand": strand,
                "ancestral_base_alignment": ancestral,
                "derived_base_alignment": derived,
                "ancestral_base_genomic": ancestral_genomic,
                "derived_base_genomic": derived_genomic,
                "six_class_change": change,
                "transition_transversion": (
                    "transition"
                    if {ancestral, derived} in ({"A", "G"}, {"C", "T"})
                    else "transversion"
                ),
                "multi_outgroup_concordant_support": (
                    "PASS"
                    if row["multi_outgroup_concordant_class"]
                    == row["primary_class"]
                    else "FAIL"
                ),
                "_position0": genomic_position,
                "_flags": set(),
                "_transcript_effects": [],
            }
        )

    by_chrom: dict[str, list[dict[str, object]]] = defaultdict(list)
    for site in sites:
        by_chrom[str(site["chrom"])].append(site)
    chrom_positions = {}
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda row: int(row["_position0"]))
        chrom_positions[chrom] = [int(row["_position0"]) for row in rows]

    gff = project / "01_reference/prepared_data/TAIR12.Col-CC.annotation.gff3"
    cds_by_transcript: dict[str, list[tuple[str, int, int, str, int]]] = defaultdict(list)
    with gff.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, feature = fields[0], fields[2]
            if chrom not in by_chrom or feature not in {"gene", "exon", "CDS"}:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            positions = chrom_positions[chrom]
            left, right = bisect.bisect_left(positions, start), bisect.bisect_left(positions, end)
            flag = {"gene": "gene", "exon": "exon", "CDS": "CDS"}[feature]
            for site in by_chrom[chrom][left:right]:
                site["_flags"].add(flag)  # type: ignore[union-attr]
            if feature == "gene":
                strand = fields[6]
                promoter_start, promoter_end = (
                    (max(0, start - 2000), start)
                    if strand == "+"
                    else (end, end + 2000)
                )
                pleft = bisect.bisect_left(positions, promoter_start)
                pright = bisect.bisect_left(positions, promoter_end)
                for site in by_chrom[chrom][pleft:pright]:
                    site["_flags"].add("promoter_2kb")  # type: ignore[union-attr]
            if feature == "CDS":
                attributes = parse_attributes(fields[8])
                phase = int(fields[7]) if fields[7] in {"0", "1", "2"} else 0
                for transcript in attributes.get("Parent", "").split(","):
                    if transcript:
                        cds_by_transcript[transcript].append(
                            (chrom, start, end, fields[6], phase)
                        )

    repeats = project / "01_reference/prepared_data/TAIR12.annotated_repeats.merged.bed"
    with repeats.open() as handle:
        for line in handle:
            chrom, start_text, end_text, *_ = line.rstrip().split("\t")
            if chrom not in by_chrom:
                continue
            start, end = int(start_text), int(end_text)
            positions = chrom_positions[chrom]
            left, right = bisect.bisect_left(positions, start), bisect.bisect_left(positions, end)
            for site in by_chrom[chrom][left:right]:
                site["_flags"].add("repeat")  # type: ignore[union-attr]

    genome = read_fasta(
        project / "01_reference/prepared_data/TAIR12.Col-CC.annotation_softmasked.fa"
    )
    site_lookup: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for site in sites:
        site_lookup[(str(site["chrom"]), int(site["_position0"]))].append(site)
    for transcript, segments in cds_by_transcript.items():
        strand = segments[0][3]
        if any(segment[3] != strand for segment in segments):
            continue
        segments.sort(key=lambda value: value[1], reverse=strand == "-")
        sequence_parts = []
        mappings: list[tuple[dict[str, object], int]] = []
        offset = 0
        for chrom, start, end, _, _phase in segments:
            part = genome[chrom][start:end].upper()
            if strand == "-":
                part = part.translate(COMPLEMENT)[::-1]
            sequence_parts.append(part)
            positions = chrom_positions.get(chrom, [])
            left, right = bisect.bisect_left(positions, start), bisect.bisect_left(positions, end)
            for site in by_chrom.get(chrom, [])[left:right]:
                position = int(site["_position0"])
                index = offset + (position - start if strand == "+" else end - 1 - position)
                mappings.append((site, index))
            offset += end - start
        sequence = "".join(sequence_parts)
        if not mappings or len(sequence) % 3:
            continue
        for site, index in mappings:
            if index < 0 or index >= len(sequence):
                continue
            derived_genomic = str(site["derived_base_genomic"])
            ancestral_genomic = str(site["ancestral_base_genomic"])
            derived_tx = (
                derived_genomic
                if strand == "+"
                else derived_genomic.translate(COMPLEMENT)
            )
            ancestral_tx = (
                ancestral_genomic
                if strand == "+"
                else ancestral_genomic.translate(COMPLEMENT)
            )
            if sequence[index] != derived_tx:
                continue
            codon_start = index - index % 3
            codon = sequence[codon_start : codon_start + 3]
            ancestral_codon = (
                codon[: index % 3] + ancestral_tx + codon[index % 3 + 1 :]
            )
            current_aa = GENETIC_CODE.get(codon, "X")
            ancestral_aa = GENETIC_CODE.get(ancestral_codon, "X")
            if "X" in (current_aa, ancestral_aa):
                effect = "CDS_unresolved"
            elif current_aa == ancestral_aa:
                effect = "synonymous"
            elif ancestral_aa == "*" and current_aa != "*":
                effect = "stop_loss"
            elif ancestral_aa != "*" and current_aa == "*":
                effect = "stop_gain"
            elif codon_start == 0 and ancestral_aa == "M" and current_aa != "M":
                effect = "start_loss"
            else:
                effect = "missense"
            site["_transcript_effects"].append(
                (effect, transcript, ancestral_codon, codon, ancestral_aa, current_aa)
            )  # type: ignore[union-attr]

    output_rows = []
    for site in sites:
        flags: set[str] = site.pop("_flags")  # type: ignore[assignment]
        transcript_effects: list[tuple[str, str, str, str, str, str]] = site.pop(
            "_transcript_effects"
        )  # type: ignore[assignment]
        site.pop("_position0")
        if "CDS" in flags:
            gene_context = "CDS"
        elif "exon" in flags:
            gene_context = "exon_nonCDS"
        elif "gene" in flags:
            gene_context = "intron_or_gene_body"
        elif "promoter_2kb" in flags:
            gene_context = "promoter_2kb"
        else:
            gene_context = "intergenic"
        best = (
            max(transcript_effects, key=lambda item: EFFECT_PRIORITY[item[0]])
            if transcript_effects
            else None
        )
        output_rows.append(
            {
                **site,
                "gene_context": gene_context,
                "repeat_overlap": "PASS" if "repeat" in flags else "FAIL",
                "coding_effect": (
                    best[0] if best else "CDS_unresolved" if "CDS" in flags else "noncoding"
                ),
                "representative_transcript": best[1] if best else "NA",
                "ancestral_codon": best[2] if best else "NA",
                "derived_codon": best[3] if best else "NA",
                "ancestral_amino_acid": best[4] if best else "NA",
                "derived_amino_acid": best[5] if best else "NA",
                "all_transcript_effects": (
                    ",".join(sorted({item[0] for item in transcript_effects}))
                    if transcript_effects
                    else "NA"
                ),
            }
        )
    write_tsv(args.output / "polarized_snp_sites_annotated.tsv", output_rows)

    spectrum = Counter((row["copy_role"], row["six_class_change"]) for row in output_rows)
    write_tsv(
        args.output / "snp_six_class_spectrum.tsv",
        [
            {
                "copy_role": role,
                "six_class_change": change,
                "polarized_snp_count": spectrum[(role, change)],
                "rate_per_joint_callable_site": (
                    f"{spectrum[(role, change)] / sum(event_callable.values()):.8f}"
                ),
            }
            for role in ("P", "D")
            for change in SIX_CLASSES
        ],
    )
    multi_spectrum = Counter(
        (row["copy_role"], row["six_class_change"])
        for row in output_rows
        if row["multi_outgroup_concordant_support"] == "PASS"
    )
    write_tsv(
        args.output / "snp_six_class_multi_outgroup_sensitivity.tsv",
        [
            {
                "copy_role": role,
                "six_class_change": change,
                "multi_outgroup_concordant_snp_count": multi_spectrum[
                    (role, change)
                ],
            }
            for role in ("P", "D")
            for change in SIX_CLASSES
        ],
    )
    event_counts: Counter[tuple[str, str, str]] = Counter()
    event_titv: Counter[tuple[str, str, str]] = Counter()
    for row in output_rows:
        event_counts[
            (str(row["event_id"]), str(row["copy_role"]), str(row["six_class_change"]))
        ] += 1
        event_titv[
            (
                str(row["event_id"]),
                str(row["copy_role"]),
                str(row["transition_transversion"]),
            )
        ] += 1
    event_rows = []
    for event in controlled:
        event_id = event["event_id"]
        row: dict[str, object] = {
            "event_id": event_id,
            "age_bin": event["age_bin"],
            "admission_tier": event["admission_tier"],
            "joint_callable_sites": event_callable[event_id],
        }
        for role in ("P", "D"):
            for change in SIX_CLASSES:
                row[f"{role}_{change}"] = event_counts[(event_id, role, change)]
            for kind in ("transition", "transversion"):
                row[f"{role}_{kind}"] = event_titv[(event_id, role, kind)]
        event_rows.append(row)
    write_tsv(args.output / "event_snp_type_metrics.tsv", event_rows)

    rng = np.random.default_rng(20260724)
    statistical = []
    for mutation in (*SIX_CLASSES, "transition", "transversion"):
        p = np.array(
            [
                int(row[f"P_{mutation}"]) / int(row["joint_callable_sites"])
                for row in event_rows
            ]
        )
        d = np.array(
            [
                int(row[f"D_{mutation}"]) / int(row["joint_callable_sites"])
                for row in event_rows
            ]
        )
        differences = d - p
        samples = rng.integers(0, len(event_rows), size=(20_000, len(event_rows)))
        medians = np.median(differences[samples], axis=1)
        nonzero = differences[differences != 0]
        sign_p = (
            stats.binomtest(int(np.sum(nonzero > 0)), len(nonzero), 0.5).pvalue
            if len(nonzero)
            else 1.0
        )
        wilcoxon_p = (
            stats.wilcoxon(d, p, zero_method="wilcox").pvalue
            if np.any(differences)
            else 1.0
        )
        statistical.append(
            {
                "mutation_class": mutation,
                "events": len(event_rows),
                "P_count": sum(int(row[f"P_{mutation}"]) for row in event_rows),
                "D_count": sum(int(row[f"D_{mutation}"]) for row in event_rows),
                "D_to_P_count_ratio": (
                    f"{sum(int(row[f'D_{mutation}']) for row in event_rows) / sum(int(row[f'P_{mutation}']) for row in event_rows):.6f}"
                    if sum(int(row[f"P_{mutation}"]) for row in event_rows)
                    else "NA"
                ),
                "median_D_minus_P_rate": f"{np.median(differences):.8f}",
                "median_difference_CI95_low": f"{np.quantile(medians, 0.025):.8f}",
                "median_difference_CI95_high": f"{np.quantile(medians, 0.975):.8f}",
                "events_D_gt_P": int(np.sum(differences > 0)),
                "events_P_gt_D": int(np.sum(differences < 0)),
                "events_tied": int(np.sum(differences == 0)),
                "event_sign_test_p": f"{sign_p:.12g}",
                "paired_wilcoxon_p": f"{wilcoxon_p:.12g}",
            }
        )
    write_tsv(args.output / "snp_type_statistical_summary.tsv", statistical)

    for column, path in (
        ("gene_context", "snp_gene_context_summary.tsv"),
        ("coding_effect", "snp_coding_effect_summary.tsv"),
        ("repeat_overlap", "snp_repeat_overlap_summary.tsv"),
    ):
        counts = Counter((row["copy_role"], row[column]) for row in output_rows)
        categories = sorted({str(row[column]) for row in output_rows})
        write_tsv(
            args.output / path,
            [
                {
                    "copy_role": role,
                    column: category,
                    "polarized_snp_count": counts[(role, category)],
                    "fraction_within_copy_role": (
                        f"{counts[(role, category)] / sum(counts[(role, value)] for value in categories):.8f}"
                    ),
                }
                for role in ("P", "D")
                for category in categories
            ],
        )
    print(f"Annotated {len(output_rows)} polarized SNPs from {len(controlled)} events")


if __name__ == "__main__":
    main()
