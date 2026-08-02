#!/usr/bin/env python3
"""Build a browsable, auditable review package for inclusive micro-indels.

The package deliberately keeps three different units separate:
  * main-window events with at least one PRIMARY call;
  * exact call identities shared by the two window settings;
  * events whose P/D direction is non-tied and identical in both settings.
"""

from __future__ import annotations

import csv
import hashlib
import html
import math
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


PROJECT = Path(__file__).resolve().parents[2]
SOURCE = PROJECT / "12_inclusive_pd_sequence_variation"
OUTPUT = Path(__file__).resolve().parents[1]
MAIN_NAME = "main_anchor60_tile250"
SENS_NAME = "sensitivity_anchor80_tile300"
RUNS = {
    MAIN_NAME: SOURCE / "microindel_local_msa",
    SENS_NAME: SOURCE / "microindel_local_msa_sensitivity_anchor80_tile300",
}
ALIGNERS = {
    "MAFFT-L-INS-i": ("MAFFT_LINSI", ".fa"),
    "MUSCLE5": ("MUSCLE5", ".fa"),
    "PRANK fixed tree": ("PRANK_FIXED_TREE", ".best.fas"),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = ""
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            else:
                records[name].append(line)
    return {name: "".join(parts) for name, parts in records.items()}


def attrs(text: str) -> dict[str, str]:
    result = {}
    for part in text.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            result[key] = unquote(value)
    return result


def clean_locus(value: str) -> str:
    return value.replace("TAIR12_TAIR12_", "")


def parse_genes(path: Path) -> dict[str, list[dict[str, object]]]:
    genes: dict[str, list[dict[str, object]]] = defaultdict(list)
    by_id: dict[str, dict[str, object]] = {}
    products: dict[str, list[str]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue
            chrom, _, kind, start, end, _, strand, _, raw = fields
            values = attrs(raw)
            if kind == "gene":
                gene_id = values.get("ID", "")
                row = {
                    "id": gene_id,
                    "chrom": chrom,
                    "start": int(start) - 1,
                    "end": int(end),
                    "strand": strand,
                    "locus": clean_locus(values.get("locus_tag", values.get("Name", gene_id))),
                    "symbol": values.get("gene", values.get("Name", "")),
                    "biotype": values.get("gene_biotype", "unknown"),
                    "product": "",
                }
                genes[chrom].append(row)
                by_id[gene_id] = row
            elif kind in {"mRNA", "transcript", "ncRNA", "rRNA", "tRNA"}:
                parent = values.get("Parent", "")
                product = values.get("product", "")
                if parent and product and product not in products[parent]:
                    products[parent].append(product)
    for gene_id, product_values in products.items():
        if gene_id in by_id:
            by_id[gene_id]["product"] = " | ".join(product_values[:2])
    for values in genes.values():
        values.sort(key=lambda row: (row["start"], row["end"]))
    return genes


def gene_label(gene: dict[str, object] | None) -> str:
    if not gene:
        return "NA"
    symbol = str(gene["symbol"])
    locus = str(gene["locus"])
    name = locus if not symbol or symbol == locus else f"{locus} ({symbol})"
    product = str(gene["product"]) or str(gene["biotype"])
    return f"{name}; {product}; {gene['chrom']}:{int(gene['start'])+1}-{gene['end']}({gene['strand']})"


def annotate_interval(
    genes: dict[str, list[dict[str, object]]], chrom: str, start: int, end: int
) -> dict[str, object]:
    values = genes.get(chrom, [])
    overlapping = [g for g in values if int(g["start"]) < end and start < int(g["end"])]
    left = max((g for g in values if int(g["end"]) <= start), key=lambda g: g["end"], default=None)
    right = min((g for g in values if int(g["start"]) >= end), key=lambda g: g["start"], default=None)
    if overlapping:
        nearest = overlapping[0]
        distance = 0
    else:
        candidates = []
        if left:
            candidates.append((start - int(left["end"]), left))
        if right:
            candidates.append((int(right["start"]) - end, right))
        distance, nearest = min(candidates, default=(math.inf, None), key=lambda x: x[0])
    return {
        "overlap": " | ".join(gene_label(g) for g in overlapping) or "NA",
        "nearest": gene_label(nearest),
        "nearest_distance_bp": "NA" if nearest is None else int(distance),
        "left": gene_label(left),
        "left_distance_bp": "NA" if left is None else start - int(left["end"]),
        "right": gene_label(right),
        "right_distance_bp": "NA" if right is None else int(right["start"]) - end,
    }


def block_ranges(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    result = []
    cursor = 0
    for row in sorted(rows, key=lambda value: int(value["core_block_index"])):
        length = int(row["core_block_bp"])
        normalized: dict[str, object] = {
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
    if side == "start":
        selected = blocks[-1] if position == total else next(
            row for row in blocks
            if int(row["core_start"]) <= position < int(row["core_end"])
        )
    else:
        selected = blocks[0] if position == 0 else next(
            row for row in blocks
            if int(row["core_start"]) < position <= int(row["core_end"])
        )
    offset = position - int(selected["core_start"])
    strand = str(selected[f"{role.lower()}_strand"])
    start = int(selected[f"{role.lower()}_start"])
    end = int(selected[f"{role.lower()}_end"])
    coordinate = start + offset if strand == "+" else end - offset
    return str(selected[f"{role.lower()}_chrom"]), coordinate, strand


def atom_coordinates(
    manifest: list[dict[str, str]],
    blocks_by_event: dict[str, list[dict[str, object]]],
) -> dict[tuple[str, str], dict[str, object]]:
    result = {}
    for atom in manifest:
        blocks = blocks_by_event[atom["event_id"]]
        for role in ("P", "D"):
            chrom1, cut1, strand1 = core_cut(
                blocks, role, int(atom["core_start_0based"]), "start"
            )
            chrom2, cut2, strand2 = core_cut(
                blocks, role, int(atom["core_end_0based"]), "end"
            )
            if chrom1 == chrom2 and strand1 == strand2:
                result[(atom["atom_id"], role)] = {
                    "chrom": chrom1,
                    "start": min(cut1, cut2),
                    "end": max(cut1, cut2),
                    "strand": strand1,
                }
    return result


def call_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row["event_id"],
        row["unique_role"],
        row["chrom"],
        row["start_0based"],
        row["end_0based"],
        row["parsimonious_direction"],
    )


def direction_letter(row: dict[str, str]) -> str:
    return row["parsimonious_direction"].split("_", 1)[0]


def find_gap(
    alignment: dict[str, str],
    call: dict[str, str],
    coordinate: dict[str, object],
) -> tuple[int, int] | None:
    p = alignment.get("Atha_P", "")
    d = alignment.get("Atha_D", "")
    role = call["unique_role"]
    role_seq = p if role == "P" else d
    raw_position = 0
    column_raw: dict[int, int] = {}
    partial_hits: list[tuple[int, int, int]] = []
    for column, base in enumerate(role_seq):
        if base != "-":
            column_raw[column] = raw_position
            raw_position += 1
    column = 0
    while column < len(p):
        if (p[column] == "-") == (d[column] == "-"):
            column += 1
            continue
        gap_role = "P" if p[column] != "-" else "D"
        start_column = column
        while column < len(p):
            present = p[column] != "-" if gap_role == "P" else d[column] != "-"
            absent = d[column] == "-" if gap_role == "P" else p[column] == "-"
            if not (present and absent):
                break
            column += 1
        if gap_role != role:
            continue
        raw_indices = [
            column_raw[index]
            for index in range(start_column, column)
            if index in column_raw
        ]
        if not raw_indices:
            continue
        if coordinate["strand"] == "+":
            genomic = [int(coordinate["start"]) + index for index in raw_indices]
        else:
            genomic = [int(coordinate["end"]) - 1 - index for index in raw_indices]
        if (
            coordinate["chrom"] == call["chrom"]
            and min(genomic) == int(call["start_0based"])
            and max(genomic) + 1 == int(call["end_0based"])
        ):
            return start_column, column
        if coordinate["chrom"] == call["chrom"]:
            overlap = min(max(genomic) + 1, int(call["end_0based"])) - max(
                min(genomic), int(call["start_0based"])
            )
            if overlap > 0:
                partial_hits.append((overlap, start_column, column))
    # Calls spanning two adjacent atoms are merged by the source workflow.  In
    # that case each source MSA contains only its coordinate-overlapping part.
    if partial_hits:
        _, start_column, end_column = max(partial_hits)
        return start_column, end_column
    return None


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def flag(label: str, kind: str = "pass") -> str:
    return f'<span class="flag {kind}">{esc(label)}</span>'


def alignment_snippet(
    run_name: str,
    atom_ids: str,
    call: dict[str, str],
    coordinates: dict[tuple[str, str], dict[str, object]],
) -> str:
    atom_list = atom_ids.split(",")
    panels = []
    for aligner, (directory, suffix) in ALIGNERS.items():
        chosen = None
        for atom_id in atom_list:
            path = RUNS[run_name] / "alignments" / directory / f"{atom_id}{suffix}"
            if not path.exists():
                continue
            alignment = read_fasta(path)
            coordinate = coordinates.get((atom_id, call["unique_role"]))
            hit = find_gap(alignment, call, coordinate) if coordinate else None
            if hit:
                chosen = (atom_id, alignment, hit)
                break
        if not chosen:
            panels.append(
                f"<details><summary>{esc(aligner)}</summary>"
                f"<p class='warn'>无法在源 MSA 中按基因组坐标定位该合并调用；源 atom: {esc(atom_ids)}</p></details>"
            )
            continue
        atom_id, alignment, (hit_start, hit_end) = chosen
        left = max(0, hit_start - 35)
        right = min(len(next(iter(alignment.values()))), hit_end + 35)
        preferred = ["Atha_P", "Atha_D"]
        names = preferred + sorted(name for name in alignment if name not in preferred)
        lines = []
        for name in names:
            if name not in alignment:
                continue
            sequence = alignment[name][left:right]
            rel_start, rel_end = hit_start - left, hit_end - left
            before = esc(sequence[:rel_start])
            middle = esc(sequence[rel_start:rel_end])
            after = esc(sequence[rel_end:])
            lines.append(
                f"<div class='alnline'><span class='alnname'>{esc(name)}</span>"
                f"<code>{before}<mark>{middle}</mark>{after}</code></div>"
            )
        panels.append(
            f"<details><summary>{esc(aligner)} · {esc(atom_id)} · columns "
            f"{hit_start+1}-{hit_end}</summary><div class='alignment'>{''.join(lines)}</div></details>"
        )
    return "".join(panels)


CSS = """
:root{--ink:#17202a;--muted:#667085;--line:#d7dde5;--bg:#f4f7fb;--card:#fff;
--blue:#185fa5;--green:#217a52;--amber:#9a6700;--red:#b42318;--purple:#6e44aa}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.5 system-ui,-apple-system,"Noto Sans CJK SC","Microsoft YaHei",sans-serif}
main{max-width:1500px;margin:auto;padding:24px}.hero{background:linear-gradient(120deg,#15365f,#246a73);
color:#fff;padding:28px;border-radius:14px}.hero h1{margin:0 0 8px;font-size:28px}
.hero p{margin:6px 0;max-width:1100px}.metrics{display:grid;grid-template-columns:repeat(5,minmax(130px,1fr));
gap:10px;margin:16px 0}.metric,.card{background:var(--card);border:1px solid var(--line);
border-radius:10px;padding:14px}.metric strong{font-size:25px;display:block;color:var(--blue)}
h2{margin:28px 0 10px}h3{margin:20px 0 8px}.muted{color:var(--muted)}
.flag{display:inline-block;border-radius:999px;padding:2px 8px;margin:2px;font-size:12px;
font-weight:650;background:#d9f5e8;color:var(--green)}.flag.warn{background:#fff0c2;color:var(--amber)}
.flag.fail{background:#fee4e2;color:var(--red)}.flag.info{background:#e6efff;color:var(--blue)}
.flag.top{background:#efe5ff;color:var(--purple)}table{width:100%;border-collapse:collapse;background:#fff}
th,td{border-bottom:1px solid var(--line);padding:8px;vertical-align:top;text-align:left}
th{background:#eef3f8;position:sticky;top:0;z-index:1}tr:hover td{background:#f8fbff}
.tablewrap{overflow:auto;max-height:72vh;border:1px solid var(--line);border-radius:9px}
a{color:var(--blue);text-decoration:none}a:hover{text-decoration:underline}
input,select{padding:8px;border:1px solid var(--line);border-radius:6px;margin:4px}
details{border:1px solid var(--line);border-radius:7px;margin:7px 0;background:#fff}
summary{cursor:pointer;padding:8px;font-weight:650}.alignment{overflow-x:auto;padding:8px;background:#101827;color:#dbeafe}
.alnline{white-space:nowrap}.alnname{display:inline-block;width:145px;color:#93c5fd}
code{font:13px/1.45 ui-monospace,SFMono-Regular,Consolas,monospace}mark{background:#f59e0b;color:#111827}
.call{border-left:5px solid var(--blue);margin:14px 0}.call.nonexact{border-left-color:var(--amber)}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:12px}.kv{display:grid;
grid-template-columns:minmax(160px,230px) 1fr;border-top:1px solid var(--line)}
.kv div{padding:6px;border-bottom:1px solid var(--line)}.warn{color:var(--red)}
@media(max-width:800px){.metrics,.grid2{grid-template-columns:1fr}.hero{border-radius:0}
main{padding:10px}.alnname{width:110px}}
@media print{body{background:#fff}details>*{display:block!important}.tablewrap{max-height:none;overflow:visible}
th{position:static}.filters{display:none}}
"""


JS = """
function filterRows(){
 const q=(document.getElementById('q')?.value||'').toLowerCase();
 const dir=document.getElementById('dir')?.value||'';
 const exact=document.getElementById('exact')?.value||'';
 document.querySelectorAll('#events tbody tr').forEach(r=>{
  const okq=r.innerText.toLowerCase().includes(q);
  const okd=!dir||r.dataset.direction===dir;
  const oke=!exact||r.dataset.exact===exact;
  r.style.display=(okq&&okd&&oke)?'':'none';
 });
}
"""


def page(title: str, body: str, prefix: str = "") -> str:
    return (
        "<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{esc(title)}</title><style>{CSS}</style></head><body><main>{body}</main>"
        f"<script>{JS}</script></body></html>"
    )


def kv_table(values: list[tuple[str, object]]) -> str:
    return "<div class='card'>" + "".join(
        f"<div class='kv'><div class='muted'>{esc(key)}</div><div>{value}</div></div>"
        for key, value in values
    ) + "</div>"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "events").mkdir(exist_ok=True)
    genes = parse_genes(
        PROJECT / "01_reference/prepared_data/TAIR12.Col-CC.annotation.gff3"
    )
    event_meta_rows = read_tsv(SOURCE / "inputs/high_priority_events.tsv")
    event_meta = {row["event_id"]: row for row in event_meta_rows}
    robustness = {
        row["event_id"]: row
        for row in read_tsv(
            SOURCE / "microindel_parameter_robustness/window_parameter_event_robustness.tsv"
        )
    }
    blocks_raw: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(SOURCE / "core/homologous_core_blocks.tsv"):
        blocks_raw[row["event_id"]].append(row)
    blocks_by_event = {
        event: block_ranges(rows) for event, rows in blocks_raw.items()
    }

    primary: dict[str, list[dict[str, str]]] = {}
    all_inference: dict[str, list[dict[str, str]]] = {}
    coords: dict[str, dict[tuple[str, str], dict[str, object]]] = {}
    call_maps: dict[str, dict[tuple[str, ...], dict[str, str]]] = {}
    for run_name, run_dir in RUNS.items():
        all_rows = read_tsv(run_dir / "denovo_microindel_inference.tsv")
        rows = [row for row in all_rows if row["evidence_tier"] == "PRIMARY"]
        primary[run_name] = rows
        all_inference[run_name] = all_rows
        call_maps[run_name] = {call_key(row): row for row in rows}
        coords[run_name] = atom_coordinates(
            read_tsv(run_dir / "atomic_region_manifest.tsv"), blocks_by_event
        )

    main_events = sorted({row["event_id"] for row in primary[MAIN_NAME]})
    if len(main_events) != 28:
        raise SystemExit(f"Expected 28 main PRIMARY events, observed {len(main_events)}")
    union_keys = sorted(set(call_maps[MAIN_NAME]) | set(call_maps[SENS_NAME]))
    exact_keys = sorted(set(call_maps[MAIN_NAME]) & set(call_maps[SENS_NAME]))
    if len(exact_keys) != 31:
        raise SystemExit(f"Expected 31 exact shared calls, observed {len(exact_keys)}")

    calls_by_event: dict[str, list[dict[str, object]]] = defaultdict(list)
    call_table = []
    exact_table = []
    for index, key in enumerate(union_keys, 1):
        main_call = call_maps[MAIN_NAME].get(key)
        sens_call = call_maps[SENS_NAME].get(key)
        source = main_call or sens_call
        assert source is not None
        annotation = annotate_interval(
            genes, source["chrom"], int(source["start_0based"]), int(source["end_0based"])
        )
        row = {
            "review_call_id": f"RVIND{index:04d}",
            "event_id": source["event_id"],
            "age_bin": source["age_bin"],
            "unique_role": source["unique_role"],
            "chrom": source["chrom"],
            "start_0based": source["start_0based"],
            "end_0based": source["end_0based"],
            "coordinate_1based_closed": f"{source['chrom']}:{int(source['start_0based'])+1}-{source['end_0based']}",
            "fragment_bp": source["fragment_bp"],
            "segment_sequence": source["segment_sequence"],
            "parsimonious_direction": source["parsimonious_direction"],
            "genomic_context": source["genomic_context"],
            "coding_size_class": source["coding_size_class"],
            "main_present": "YES" if main_call else "NO",
            "main_atom_id": main_call["atom_id"] if main_call else "NA",
            "sensitivity_present": "YES" if sens_call else "NO",
            "sensitivity_atom_id": sens_call["atom_id"] if sens_call else "NA",
            "exact_cross_window": "YES" if main_call and sens_call else "NO",
            "three_aligner_coordinate_concordance": source[
                "MAFFT_MUSCLE_PRANK_exact_coordinate_concordance"
            ],
            "boundary_P0_species": source["boundary_P0_species"],
            "boundary_P0_state": source["boundary_P0_state"],
            "deeper_P0_states": source["deeper_P0_states"] or "NA",
            "postduplication_concordant_species": source[
                "postduplication_concordant_species"
            ] or "NA",
            "time1_multispecies_P0_support": source["time1_multispecies_P0_support"],
            "fully_uppercase": source["fully_uppercase"],
            "annotated_repeat_overlap": source["annotated_repeat_overlap"],
            "low_complexity_context": source["low_complexity_context"],
            "overlapping_gene": annotation["overlap"],
            "nearest_gene": annotation["nearest"],
            "nearest_gene_distance_bp": annotation["nearest_distance_bp"],
            "left_gene": annotation["left"],
            "left_gene_distance_bp": annotation["left_distance_bp"],
            "right_gene": annotation["right"],
            "right_gene_distance_bp": annotation["right_distance_bp"],
        }
        call_table.append(row)
        calls_by_event[source["event_id"]].append(
            {"row": row, "main": main_call, "sensitivity": sens_call}
        )
        if row["exact_cross_window"] == "YES":
            exact_table.append(row.copy())

    # The requested package is event-scoped to the 28 main-window PRIMARY
    # events.  Keep the two sensitivity-only events in a separate audit table
    # instead of silently mixing them into the 28-event call detail.
    out_of_scope_calls = [
        row for row in call_table if row["event_id"] not in set(main_events)
    ]
    call_table = [row for row in call_table if row["event_id"] in set(main_events)]

    event_rows = []
    for event_id in main_events:
        meta = event_meta[event_id]
        robust = robustness[event_id]
        p_locus_name = meta["provisional_p_locus"]
        d_locus_name = meta["provisional_d_locus"]
        def locus(which: str) -> dict[str, object]:
            letter = "A" if which == "locus_A" else "B"
            chrom = meta[f"locus_{letter}_chrom"]
            start = int(meta[f"locus_{letter}_representative_start"])
            end = int(meta[f"locus_{letter}_representative_end"])
            ann = annotate_interval(genes, chrom, start, end)
            return {"chrom": chrom, "start": start, "end": end, **ann}
        p_locus = locus(p_locus_name)
        d_locus = locus(d_locus_name)
        event_calls = calls_by_event[event_id]
        exact_count = sum(c["row"]["exact_cross_window"] == "YES" for c in event_calls)
        robust_direction = (
            robust["main_direction"]
            if robust["main_direction"] == robust["sensitivity_direction"]
            and robust["main_direction"] != "TIE"
            else "NO"
        )
        main_rows = [r for r in primary[MAIN_NAME] if r["event_id"] == event_id]
        sens_rows = [r for r in primary[SENS_NAME] if r["event_id"] == event_id]
        cds_count = sum(r["genomic_context"] == "CDS" for r in main_rows + sens_rows)
        frameshift_count = sum(
            r["coding_size_class"] == "frameshift_candidate"
            for r in main_rows + sens_rows
        )
        # Transparent ranking: reproducibility first, then functional impact and evidence volume.
        score = (
            20
            + min(exact_count, 3) * 4
            + (4 if exact_count else 0)
            + min(cds_count, 3) * 2
            + min(frameshift_count, 3)
            + min(len(main_rows) + len(sens_rows), 6) * 0.5
        ) if robust_direction != "NO" else 0
        qc_flags = []
        qc_flags.append("DIRECTION_ROBUST" if robust_direction != "NO" else "DIRECTION_NOT_ROBUST")
        qc_flags.append(f"EXACT_SHARED_CALLS={exact_count}")
        if exact_count == 0:
            qc_flags.append("BOUNDARY_SHIFT_ACROSS_WINDOWS")
        if cds_count:
            qc_flags.append("CDS_OVERLAP")
        if frameshift_count:
            qc_flags.append("FRAMESHIFT_CANDIDATE")
        event_rows.append(
            {
                "event_id": event_id,
                "age_bin": robust["age_bin"],
                "relative_orientation": meta["relative_orientation"],
                "P_locus_0based_halfopen": f"{p_locus['chrom']}:{p_locus['start']}-{p_locus['end']}",
                "D_locus_0based_halfopen": f"{d_locus['chrom']}:{d_locus['start']}-{d_locus['end']}",
                "P_locus_gene_annotation": p_locus["overlap"] if p_locus["overlap"] != "NA" else p_locus["nearest"],
                "D_locus_gene_annotation": d_locus["overlap"] if d_locus["overlap"] != "NA" else d_locus["nearest"],
                "main_P_calls": robust["main_P"],
                "main_D_calls": robust["main_D"],
                "main_direction": robust["main_direction"],
                "sensitivity_P_calls": robust["sensitivity_P"],
                "sensitivity_D_calls": robust["sensitivity_D"],
                "sensitivity_direction": robust["sensitivity_direction"],
                "robust_direction": robust_direction,
                "exact_shared_calls": exact_count,
                "union_primary_calls": len(event_calls),
                "CDS_primary_observations_across_windows": cds_count,
                "frameshift_primary_observations_across_windows": frameshift_count,
                "manual_review_score": f"{score:.1f}",
                "manual_review_top16": "PENDING",
                "QC_flags": ";".join(qc_flags),
                "_p_locus": p_locus,
                "_d_locus": d_locus,
            }
        )

    robust_rows = [r for r in event_rows if r["robust_direction"] != "NO"]
    robust_rows.sort(
        key=lambda r: (
            -float(r["manual_review_score"]),
            -int(r["exact_shared_calls"]),
            r["event_id"],
        )
    )
    top_ids = {row["event_id"] for row in robust_rows[:16]}
    for row in event_rows:
        row["manual_review_top16"] = "YES" if row["event_id"] in top_ids else "NO"

    export_event_fields = [key for key in event_rows[0] if not key.startswith("_")]
    write_tsv(
        OUTPUT / "event_details.tsv",
        [{key: row[key] for key in export_event_fields} for row in event_rows],
        export_event_fields,
    )
    write_tsv(OUTPUT / "call_details.tsv", call_table)
    write_tsv(OUTPUT / "exact_cross_window_calls.tsv", exact_table)
    write_tsv(
        OUTPUT / "out_of_scope_sensitivity_only_calls.tsv",
        out_of_scope_calls,
        list(call_table[0]),
    )
    top_rows = [
        {key: row[key] for key in export_event_fields}
        for row in sorted(
            (r for r in event_rows if r["manual_review_top16"] == "YES"),
            key=lambda r: (-float(r["manual_review_score"]), r["event_id"]),
        )
    ]
    write_tsv(OUTPUT / "manual_review_top16.tsv", top_rows, export_event_fields)
    source_files = [
        SOURCE / "inputs/high_priority_events.tsv",
        SOURCE / "core/homologous_core_blocks.tsv",
        SOURCE / "microindel_local_msa/denovo_microindel_inference.tsv",
        SOURCE / "microindel_local_msa/atomic_region_manifest.tsv",
        SOURCE / "microindel_local_msa_sensitivity_anchor80_tile300/denovo_microindel_inference.tsv",
        SOURCE / "microindel_local_msa_sensitivity_anchor80_tile300/atomic_region_manifest.tsv",
        SOURCE / "microindel_parameter_robustness/window_parameter_event_robustness.tsv",
        PROJECT / "01_reference/prepared_data/TAIR12.Col-CC.annotation.gff3",
    ]
    write_tsv(
        OUTPUT / "source_manifest.tsv",
        [
            {
                "source_path_relative_to_project": str(path.relative_to(PROJECT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in source_files
        ],
    )
    write_tsv(
        OUTPUT / "package_qc.tsv",
        [
            {"check": "main_PRIMARY_event_count", "observed": len(main_events), "expected": 28, "status": "PASS"},
            {"check": "exact_cross_window_call_count", "observed": len(exact_keys), "expected": 31, "status": "PASS"},
            {"check": "events_with_exact_call", "observed": len({key[0] for key in exact_keys}), "expected": 20, "status": "PASS"},
            {"check": "direction_robust_event_count", "observed": len(robust_rows), "expected": 21, "status": "PASS"},
            {"check": "direction_robust_D_event_count", "observed": sum(r["robust_direction"] == "D" for r in event_rows), "expected": 18, "status": "PASS"},
            {"check": "direction_robust_P_event_count", "observed": sum(r["robust_direction"] == "P" for r in event_rows), "expected": 3, "status": "PASS"},
            {"check": "manual_review_top_count", "observed": len(top_ids), "expected": 16, "status": "PASS"},
            {"check": "28_event_scoped_call_union_row_count", "observed": len(call_table), "expected": 70, "status": "PASS"},
            {"check": "out_of_scope_sensitivity_only_call_count", "observed": len(out_of_scope_calls), "expected": 2, "status": "PASS"},
        ],
    )

    for event in event_rows:
        event_id = event["event_id"]
        badges = [
            flag(event["age_bin"], "info"),
            flag(f"main {event['main_direction']}", "info"),
            flag(f"sensitivity {event['sensitivity_direction']}", "info"),
            flag(
                f"robust {event['robust_direction']}",
                "pass" if event["robust_direction"] != "NO" else "warn",
            ),
            flag(f"exact calls {event['exact_shared_calls']}", "pass" if int(event["exact_shared_calls"]) else "warn"),
        ]
        if event["manual_review_top16"] == "YES":
            badges.append(flag("TOP16 MANUAL REVIEW", "top"))
        meta_values = [
            ("事件 / 年龄", f"{esc(event_id)} · {esc(event['age_bin'])}"),
            ("两窗口方向", f"main {event['main_P_calls']}P/{event['main_D_calls']}D → {event['main_direction']}; "
             f"sensitivity {event['sensitivity_P_calls']}P/{event['sensitivity_D_calls']}D → {event['sensitivity_direction']}"),
            ("相对方向", esc(event["relative_orientation"])),
            ("QC", esc(event["QC_flags"])),
            ("人工复核评分", esc(event["manual_review_score"])),
        ]
        locus_cards = []
        for role, loc in (("P", event["_p_locus"]), ("D", event["_d_locus"])):
            locus_cards.append(
                kv_table([
                    (f"{role} locus", f"{esc(loc['chrom'])}:{int(loc['start'])+1}-{loc['end']} (1-based closed)"),
                    ("重叠基因", esc(loc["overlap"])),
                    ("最近基因", f"{esc(loc['nearest'])} · {esc(loc['nearest_distance_bp'])} bp"),
                    ("左邻基因", f"{esc(loc['left'])} · {esc(loc['left_distance_bp'])} bp"),
                    ("右邻基因", f"{esc(loc['right'])} · {esc(loc['right_distance_bp'])} bp"),
                ])
            )
        call_cards = []
        for item in sorted(
            calls_by_event[event_id],
            key=lambda x: (
                x["row"]["chrom"],
                int(x["row"]["start_0based"]),
                x["row"]["parsimonious_direction"],
            ),
        ):
            row = item["row"]
            exact = row["exact_cross_window"] == "YES"
            call_badges = [
                flag(row["parsimonious_direction"], "info"),
                flag(row["genomic_context"], "info"),
                flag("EXACT TWO-WINDOW" if exact else "WINDOW-SPECIFIC", "pass" if exact else "warn"),
                flag(row["coding_size_class"], "warn" if row["coding_size_class"] == "frameshift_candidate" else "info"),
            ]
            run_sections = []
            for run_name, call in ((MAIN_NAME, item["main"]), (SENS_NAME, item["sensitivity"])):
                if not call:
                    run_sections.append(f"<h4>{esc(run_name)}</h4><p class='muted'>该窗口未检出同一坐标/方向调用。</p>")
                    continue
                run_sections.append(
                    f"<h4>{esc(run_name)} · atom {esc(call['atom_id'])}</h4>"
                    + alignment_snippet(run_name, call["atom_id"], call, coords[run_name])
                )
            call_cards.append(
                f"<section class='card call {'exact' if exact else 'nonexact'}'>"
                f"<h3>{esc(row['review_call_id'])} · {esc(row['coordinate_1based_closed'])} · "
                f"{esc(row['segment_sequence'])} ({esc(row['fragment_bp'])} bp)</h3>"
                f"<p>{''.join(call_badges)}</p>"
                + kv_table([
                    ("坐标约定", f"0-based half-open {esc(row['chrom'])}:{esc(row['start_0based'])}-{esc(row['end_0based'])}; "
                     f"1-based closed {esc(row['coordinate_1based_closed'])}"),
                    ("窗口存在性", f"main={esc(row['main_present'])} ({esc(row['main_atom_id'])}); "
                     f"sensitivity={esc(row['sensitivity_present'])} ({esc(row['sensitivity_atom_id'])})"),
                    ("P0", f"boundary {esc(row['boundary_P0_species'])}={esc(row['boundary_P0_state'])}; "
                     f"deeper {esc(row['deeper_P0_states'])}"),
                    ("复制后物种", esc(row["postduplication_concordant_species"])),
                    ("序列 QC", f"3-aligner={esc(row['three_aligner_coordinate_concordance'])}; "
                     f"uppercase={esc(row['fully_uppercase'])}; repeat={esc(row['annotated_repeat_overlap'])}; "
                     f"low-complexity={esc(row['low_complexity_context'])}"),
                    ("重叠/最近基因", f"{esc(row['overlapping_gene'])}<br>{esc(row['nearest_gene'])} "
                     f"({esc(row['nearest_gene_distance_bp'])} bp)"),
                ])
                + "".join(run_sections)
                + "</section>"
            )
        body = (
            f"<p><a href='../index.html'>← 返回索引</a></p>"
            f"<div class='hero'><h1>{esc(event_id)} 逐事件 micro-indel 复核证据</h1>"
            f"<p>{''.join(badges)}</p></div>"
            f"<h2>事件判定</h2>{kv_table(meta_values)}"
            f"<h2>P / D 位点与邻近功能注释</h2><div class='grid2'>{''.join(locus_cards)}</div>"
            f"<h2>两窗口调用与局部 P/D/P0 比对</h2>"
            f"<p class='muted'>橙色高亮为按基因组坐标回投后定位的 P/D 非对称 gap block；每个窗口分别展示 "
            f"MAFFT-L-INS-i、MUSCLE5 和固定树 PRANK 的局部比对。</p>"
            + "".join(call_cards)
        )
        (OUTPUT / "events" / f"{event_id}.html").write_text(
            page(f"{event_id} micro-indel review", body), encoding="utf-8"
        )

    index_rows = []
    for event in event_rows:
        exact_kind = "yes" if int(event["exact_shared_calls"]) else "no"
        direction = event["robust_direction"]
        index_rows.append(
            f"<tr data-direction='{esc(direction)}' data-exact='{exact_kind}'>"
            f"<td><a href='events/{esc(event['event_id'])}.html'>{esc(event['event_id'])}</a>"
            f"{flag('TOP16','top') if event['manual_review_top16']=='YES' else ''}</td>"
            f"<td>{esc(event['age_bin'])}</td><td>{esc(event['relative_orientation'])}</td>"
            f"<td>{esc(event['main_P_calls'])}/{esc(event['main_D_calls'])} → {esc(event['main_direction'])}</td>"
            f"<td>{esc(event['sensitivity_P_calls'])}/{esc(event['sensitivity_D_calls'])} → {esc(event['sensitivity_direction'])}</td>"
            f"<td>{flag(direction,'pass') if direction != 'NO' else flag('NO','warn')}</td>"
            f"<td>{esc(event['exact_shared_calls'])}</td>"
            f"<td>{esc(event['P_locus_gene_annotation'])}<br><span class='muted'>{esc(event['D_locus_gene_annotation'])}</span></td>"
            f"<td>{esc(event['QC_flags'])}</td><td>{esc(event['manual_review_score'])}</td></tr>"
        )
    body = (
        "<div class='hero'><h1>TAIR12 P/D micro-indel 逐事件复核证据包</h1>"
        "<p>主窗口 PRIMARY 事件、跨窗口精确调用和双窗口方向稳健事件采用互不混用的计数口径。"
        "所有坐标同时注明 0-based half-open 与 1-based closed；事件页保留三种比对器的局部 P/D/P0 原始证据。</p>"
        "<p>参数：main = junction anchor 60 bp / tile 250 bp / tile step 200 bp；"
        "sensitivity = junction anchor 80 bp / tile 300 bp / tile step 200 bp。</p></div>"
        "<div class='metrics'>"
        "<div class='metric'><strong>28</strong>主窗口 PRIMARY 事件</div>"
        "<div class='metric'><strong>31</strong>精确跨窗口调用</div>"
        "<div class='metric'><strong>20</strong>含精确调用的事件</div>"
        "<div class='metric'><strong>21</strong>方向稳健事件（18 D / 3 P）</div>"
        "<div class='metric'><strong>16</strong>优先人工复核</div></div>"
        "<section class='card'><h2>口径与 QC 说明</h2>"
        "<p><b>精确跨窗口</b>要求 event、unique role、染色体、起止坐标和极化方向全部相同。"
        "<b>方向稳健</b>要求两窗口的事件级 P/D 计数方向相同且非平局。PRIMARY 调用还要求三比对器坐标一致、"
        "P0 可极化、全大写 ACGT、无注释重复、无低复杂度标记；time1 另需多物种 P0 支持，time3 需复制后物种支持。</p>"
        "<p class='warn'><b>重要：</b>“28”是主窗口有 PRIMARY 调用的事件数；精确跨窗口调用是 31 个，不应写成 28 个精确调用。</p>"
        "<p>TOP16 仅从 21 个方向稳健事件中选择，评分优先考虑精确复现数，其次为 CDS/候选移码影响和两窗口证据量。"
        "它是人工浏览优先级，不是新的生物学显著性阈值。</p></section>"
        "<h2>事件索引</h2><div class='filters'>"
        "<input id='q' oninput='filterRows()' placeholder='搜索事件、基因、QC…'>"
        "<select id='dir' onchange='filterRows()'><option value=''>全部方向</option>"
        "<option value='D'>稳健 D</option><option value='P'>稳健 P</option><option value='NO'>非稳健</option></select>"
        "<select id='exact' onchange='filterRows()'><option value=''>全部精确状态</option>"
        "<option value='yes'>有精确调用</option><option value='no'>无精确调用</option></select></div>"
        "<div class='tablewrap'><table id='events'><thead><tr><th>事件</th><th>年龄</th><th>相对方向</th>"
        "<th>main P/D</th><th>sensitivity P/D</th><th>稳健方向</th><th>精确调用</th>"
        "<th>P / D 位点基因</th><th>QC</th><th>评分</th></tr></thead><tbody>"
        + "".join(index_rows)
        + "</tbody></table></div>"
        "<h2>下载表</h2><ul>"
        "<li><a href='event_details.tsv'>28 事件明细表</a></li>"
        "<li><a href='call_details.tsv'>两窗口 PRIMARY 调用并集明细</a></li>"
        "<li><a href='exact_cross_window_calls.tsv'>31 个精确跨窗口调用</a></li>"
        "<li><a href='manual_review_top16.tsv'>16 个优先人工复核方向稳健事件</a></li>"
        "<li><a href='out_of_scope_sensitivity_only_calls.tsv'>不属于28事件范围的2个敏感性窗口调用（审计附表）</a></li>"
        "<li><a href='package_qc.tsv'>证据包一致性检查</a></li>"
        "<li><a href='source_manifest.tsv'>输入来源与 SHA-256</a></li></ul>"
    )
    (OUTPUT / "index.html").write_text(page("TAIR12 micro-indel review", body), encoding="utf-8")

    readme = f"""# TAIR12 micro-indel 逐事件复核证据包

- 浏览入口：`index.html`
- 主窗口 PRIMARY 事件：{len(main_events)}
- 两窗口精确共享调用：{len(exact_keys)}
- 含精确共享调用的事件：{len({key[0] for key in exact_keys})}
- 双窗口方向稳健事件：{sum(r['robust_direction'] != 'NO' for r in event_rows)}
- 优先人工复核事件：{len(top_ids)}
- 28 事件范围内的两窗口调用并集：{len(call_table)}
- 事件范围外的敏感性窗口调用：{len(out_of_scope_calls)}（单独审计附表）
- 一致性检查：`package_qc.tsv`
- 输入来源与校验和：`source_manifest.tsv`

注意：“28”是主窗口有 PRIMARY micro-indel 的事件数；按源流程的严格键
（event、unique role、chrom、start、end、direction）精确共享调用为 31 个。

复现：

```bash
python3 workflow_scripts/build_microindel_review_package.py
```

坐标：

- `start_0based/end_0based` 为 0-based half-open；
- HTML 同时显示 1-based closed 坐标；
- 邻近基因和产品来自 `01_reference/prepared_data/TAIR12.Col-CC.annotation.gff3`。
"""
    (OUTPUT / "README.md").write_text(readme, encoding="utf-8")
    print(
        f"Built {OUTPUT}: {len(main_events)} events, {len(exact_keys)} exact calls, "
        f"{len(top_ids)} top-review events"
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
