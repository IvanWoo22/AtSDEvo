#!/usr/bin/env python3
"""Build publication-style and inline visualizations of representative SD events."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


EXAMPLES = [
    ("SDJGI00300", "典型D偏高"),
    ("SDJGI00368", "强D偏高"),
    ("SDJGI00274", "部分近缘P/D缺失"),
    ("SDJGI00762", "更深P0替代"),
    ("SDJGI00352", "P偏高反例"),
]
EXAMPLE_EN = {
    "SDJGI00300": "typical D excess",
    "SDJGI00368": "strong D excess",
    "SDJGI00274": "partial post-duplication mapping",
    "SDJGI00762": "deeper P0 fallback",
    "SDJGI00352": "P-excess counterexample",
}
AGE_COLORS = {"time1": "#3973ac", "time2": "#e17c30", "time3": "#3d9970"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def event_sites(root: Path, tiers: dict[str, str]) -> dict[str, list[dict[str, str]]]:
    pilot = root / "core_500"
    paths = {
        "core_500_strict": pilot / "sequence_variation/polarized_sites.tsv",
        "partial_postduplication": (
            pilot / "sequence_variation_partial_postdup/polarized_sites.tsv"
        ),
        "deeper_P0_fallback": (
            pilot / "sequence_variation_deeper_P0/polarized_sites.tsv"
        ),
    }
    wanted = set(tiers)
    result: dict[str, list[dict[str, str]]] = {event: [] for event in wanted}
    for tier, path in paths.items():
        tier_events = {event for event in wanted if tiers[event] == tier}
        if not tier_events:
            continue
        for row in read_tsv(path):
            if row["event_id"] in tier_events:
                result[row["event_id"]].append(row)
    return result


def locus_label(event: dict[str, str], role: str) -> str:
    p_copy = event["provisional_p_copy"]
    copy = p_copy if role == "P" else ("copy2" if p_copy == "copy1" else "copy1")
    return (
        f"{event[f'representative_{copy}_chrom']}:"
        f"{int(event[f'representative_{copy}_start']) + 1:,}–"
        f"{int(event[f'representative_{copy}_end']):,}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", required=True, type=Path)
    parser.add_argument("--html-output", required=True, type=Path)
    args = parser.parse_args()
    root = args.analysis_root
    project = root.parent
    pilot = root / "core_500"
    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    events = read_tsv(root / "controlled_expansion_endpoint_events.tsv")
    by_id = {row["event_id"]: row for row in events}
    tiers = {row["event_id"]: row["admission_tier"] for row in events}
    core_meta = {
        row["event_id"]: row
        for row in read_tsv(pilot / "inputs/high_priority_core_eligible.tsv")
    }
    strict_queue = read_tsv(pilot / "pilot/strict_primary_event_queue.tsv")
    partial_queue = read_tsv(pilot / "pilot/partial_postdup_event_queue.tsv")
    fallback_queue = read_tsv(root / "deeper_P0_sensitivity_queue.tsv")
    queue = {}
    for row in strict_queue + partial_queue + fallback_queue:
        queue[row["event_id"]] = row
    sites = event_sites(root, tiers)

    example_data = []
    for event_id, interpretation in EXAMPLES:
        row = by_id[event_id]
        positions = [
            int(site["core_position_1based"])
            for site in sites[event_id]
            if site["primary_class"] in {"P_specific", "D_specific"}
        ]
        core_length = max(
            int(row["primary_bidir_callable_sites"]),
            max(positions, default=1),
        )
        edges = np.linspace(1, core_length + 1, 41)
        p_positions = [
            int(site["core_position_1based"])
            for site in sites[event_id]
            if site["primary_class"] == "P_specific"
        ]
        d_positions = [
            int(site["core_position_1based"])
            for site in sites[event_id]
            if site["primary_class"] == "D_specific"
        ]
        p_hist, _ = np.histogram(p_positions, bins=edges)
        d_hist, _ = np.histogram(d_positions, bins=edges)
        meta = core_meta[event_id]
        q = queue[event_id]
        example_data.append(
            {
                "event_id": event_id,
                "interpretation": interpretation,
                "age": row["age_bin"],
                "tier": row["admission_tier"],
                "callable": int(row["primary_bidir_callable_sites"]),
                "p_count": int(row["P_specific_changes"]),
                "d_count": int(row["D_specific_changes"]),
                "p_rate": float(row["P_specific_rate"]),
                "d_rate": float(row["D_specific_rate"]),
                "core_length": core_length,
                "p_locus": locus_label(meta, "P"),
                "d_locus": locus_label(meta, "D"),
                "p0_species": q["boundary_P0_species"],
                "bin_edges": [round(value) for value in edges.tolist()],
                "p_bins": p_hist.tolist(),
                "d_bins": d_hist.tolist(),
            }
        )

    scatter_data = [
        {
            "id": row["event_id"],
            "age": row["age_bin"],
            "tier": row["admission_tier"],
            "p": float(row["P_specific_rate"]),
            "d": float(row["D_specific_rate"]),
            "callable": int(row["primary_bidir_callable_sites"]),
            "selected": row["event_id"] in {event for event, _ in EXAMPLES},
        }
        for row in events
    ]
    payload = {"scatter": scatter_data, "examples": example_data}
    (root / "sd_example_visualization_data.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig = plt.figure(figsize=(10.5, 12.5), constrained_layout=True)
    grid = fig.add_gridspec(7, 2, height_ratios=[2.6, 0.1, 1, 1, 1, 1, 1])
    ax_scatter = fig.add_subplot(grid[0, 0])
    ax_bar = fig.add_subplot(grid[0, 1])
    selected_ids = {event for event, _ in EXAMPLES}
    for age in ("time1", "time2", "time3"):
        subset = [row for row in scatter_data if row["age"] == age]
        ax_scatter.scatter(
            [row["p"] for row in subset],
            [row["d"] for row in subset],
            s=[18 + np.sqrt(row["callable"]) for row in subset],
            color=AGE_COLORS[age],
            alpha=0.72,
            edgecolor="none",
            label=age,
        )
    maximum = max(max(row["p"], row["d"]) for row in scatter_data) * 1.07
    ax_scatter.plot([0, maximum], [0, maximum], color="#666666", lw=1, ls="--")
    for row in scatter_data:
        if row["id"] in selected_ids:
            ax_scatter.scatter(
                row["p"], row["d"], s=78, facecolor="none", edgecolor="#202020", lw=1.2
            )
            ax_scatter.annotate(
                row["id"].replace("SDJGI", ""),
                (row["p"], row["d"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7,
            )
    ax_scatter.set(xlim=(0, maximum), ylim=(0, maximum))
    ax_scatter.set_xlabel("P-specific substitution rate")
    ax_scatter.set_ylabel("D-specific substitution rate")
    ax_scatter.set_title("A  Controlled 55-event set", loc="left", fontweight="normal")
    ax_scatter.legend(frameon=False, ncol=3, fontsize=8)

    labels = [item["event_id"].replace("SDJGI", "") for item in example_data]
    y = np.arange(len(example_data))
    p_counts = [item["p_count"] for item in example_data]
    d_counts = [item["d_count"] for item in example_data]
    ax_bar.barh(y + 0.18, p_counts, height=0.34, color="#3973ac", label="P-specific")
    ax_bar.barh(y - 0.18, d_counts, height=0.34, color="#c45145", label="D-specific")
    ax_bar.set_yticks(y, labels)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Polarized substitutions")
    ax_bar.set_title("B  Representative events", loc="left", fontweight="normal")
    ax_bar.legend(frameon=False, fontsize=8)

    for index, item in enumerate(example_data):
        ax = fig.add_subplot(grid[index + 2, :])
        centers = (
            np.array(item["bin_edges"][:-1]) + np.array(item["bin_edges"][1:])
        ) / 2
        widths = np.diff(item["bin_edges"]) * 0.86
        ax.bar(centers, item["d_bins"], width=widths, color="#c45145", label="D-specific")
        ax.bar(
            centers,
            -np.array(item["p_bins"]),
            width=widths,
            color="#3973ac",
            label="P-specific",
        )
        ax.axhline(0, color="#666666", lw=0.7)
        limit = max(max(item["d_bins"]), max(item["p_bins"]), 1)
        ax.set_ylim(-limit * 1.35, limit * 1.35)
        ax.set_xlim(1, item["core_length"])
        ax.set_ylabel("sites/bin")
        ax.set_xlabel("BISER homologous-core coordinate (bp)")
        ax.set_title(
            f"{chr(67 + index)}  {item['event_id']} · "
            f"{EXAMPLE_EN[item['event_id']]} · "
            f"{item['age']} · P={item['p_count']}, D={item['d_count']} · "
            f"P0={item['p0_species']}",
            loc="left",
            fontsize=8.5,
            fontweight="normal",
        )
        if index == 0:
            ax.legend(frameon=False, ncol=2, fontsize=8, loc="upper right")
    for suffix in ("png", "pdf"):
        fig.savefig(
            figures / f"representative_sd_variation_examples.{suffix}",
            dpi=300 if suffix == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)

    data_json = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    fragment = f"""
<div id="sd-variation-examples" style="width:100%;color:var(--foreground);">
  <style>
    #sd-variation-examples .panel-label {{ margin: 0 0 6px; font-weight: 500; }}
    #sd-variation-examples .legend {{ display:flex; flex-wrap:wrap; gap:12px; margin:4px 0 10px; color:var(--muted-foreground); }}
    #sd-variation-examples .legend span {{ display:inline-flex; align-items:center; gap:5px; }}
    #sd-variation-examples .swatch {{ width:10px; height:10px; display:inline-block; }}
    #sd-variation-examples svg {{ width:100%; height:auto; display:block; overflow:visible; }}
    #sd-variation-examples .example {{ margin-top:12px; }}
    #sd-variation-examples .meta {{ color:var(--muted-foreground); margin-bottom:2px; }}
    #sd-variation-examples .axis {{ stroke:var(--border); stroke-width:1; }}
    #sd-variation-examples .grid {{ stroke:var(--border); stroke-width:1; opacity:.55; }}
    #sd-variation-examples text {{ fill:var(--foreground); font:12px system-ui,sans-serif; }}
    #sd-variation-examples .minor {{ fill:var(--muted-foreground); font-size:11px; }}
    @media (max-width:520px) {{
      #sd-variation-examples text {{ font-size:10px; }}
      #sd-variation-examples .minor {{ font-size:9px; }}
    }}
  </style>
  <div class="panel-label">A　55个受控事件的P/D特异替换率</div>
  <div id="sd-scatter"></div>
  <div class="legend">
    <span><i class="swatch" style="background:var(--viz-series-3)"></i>time1</span>
    <span><i class="swatch" style="background:var(--viz-series-4)"></i>time2</span>
    <span><i class="swatch" style="background:var(--viz-series-5)"></i>time3</span>
    <span>虚线：P=D；圈出：下方实例</span>
  </div>
  <div class="panel-label">B　代表性SD事件的同源核心变异轨迹（40等宽窗口）</div>
  <div class="legend">
    <span><i class="swatch" style="background:var(--viz-series-1)"></i>P-specific（向下）</span>
    <span><i class="swatch" style="background:var(--viz-series-2)"></i>D-specific（向上）</span>
  </div>
  <div id="sd-tracks"></div>
  <script>
  (() => {{
    const data={data_json};
    const root=document.getElementById("sd-variation-examples");
    const NS="http://www.w3.org/2000/svg";
    const make=(tag,attrs={{}},text="")=>{{
      const el=document.createElementNS(NS,tag);
      Object.entries(attrs).forEach(([k,v])=>el.setAttribute(k,v));
      if(text) el.textContent=text;
      return el;
    }};
    const scatter=root.querySelector("#sd-scatter");
    const W=700,H=345,m={{l:58,r:22,t:18,b:48}},iw=W-m.l-m.r,ih=H-m.t-m.b;
    const max=Math.max(...data.scatter.flatMap(d=>[d.p,d.d]))*1.08;
    const sx=x=>m.l+x/max*iw, sy=y=>m.t+ih-y/max*ih;
    const svg=make("svg",{{viewBox:`0 0 ${{W}} ${{H}}`,"aria-label":"P和D特异替换率散点图"}});
    [0,.25,.5,.75,1].forEach(f=>{{
      const val=max*f, x=sx(val),y=sy(val);
      svg.append(make("line",{{x1:x,y1:m.t,x2:x,y2:m.t+ih,class:"grid"}}));
      svg.append(make("line",{{x1:m.l,y1:y,x2:m.l+iw,y2:y,class:"grid"}}));
      svg.append(make("text",{{x:x,y:m.t+ih+18,"text-anchor":"middle",class:"minor"}},val.toFixed(2)));
      svg.append(make("text",{{x:m.l-8,y:y+4,"text-anchor":"end",class:"minor"}},val.toFixed(2)));
    }});
    svg.append(make("line",{{x1:m.l,y1:m.t+ih,x2:m.l+iw,y2:m.t,class:"axis","stroke-dasharray":"5 4"}}));
    const ageColor={{time1:"var(--viz-series-3)",time2:"var(--viz-series-4)",time3:"var(--viz-series-5)"}};
    data.scatter.forEach(d=>{{
      const c=make("circle",{{cx:sx(d.p),cy:sy(d.d),r:d.selected?6:3.5,fill:ageColor[d.age],opacity:d.selected?1:.72,stroke:d.selected?"var(--foreground)":"none","stroke-width":d.selected?1.5:0}});
      c.append(make("title",{{}},`${{d.id}} · ${{d.age}} · P=${{d.p.toFixed(4)}} · D=${{d.d.toFixed(4)}} · callable=${{d.callable}}`));
      svg.append(c);
      if(d.selected) svg.append(make("text",{{x:sx(d.p)+7,y:sy(d.d)-6,class:"minor"}},d.id.replace("SDJGI","")));
    }});
    svg.append(make("text",{{x:m.l+iw/2,y:H-6,"text-anchor":"middle"}},"P-specific rate"));
    const yl=make("text",{{x:14,y:m.t+ih/2,"text-anchor":"middle",transform:`rotate(-90 14 ${{m.t+ih/2}})`}},"D-specific rate");
    svg.append(yl); scatter.append(svg);
    const tracks=root.querySelector("#sd-tracks");
    data.examples.forEach((d,i)=>{{
      const wrap=document.createElement("div"); wrap.className="example";
      const meta=document.createElement("div"); meta.className="meta";
      meta.textContent=`${{d.event_id}} · ${{d.interpretation}} · ${{d.age}} · ${{d.tier}} · P=${{d.p_count}}, D=${{d.d_count}} · P0=${{d.p0_species}}`;
      wrap.append(meta);
      const loc=document.createElement("div"); loc.className="meta";
      loc.textContent=`P: ${{d.p_locus}}　D: ${{d.d_locus}}　exact-callable: ${{d.callable}} bp`;
      wrap.append(loc);
      const tw=700,th=116,tm={{l:42,r:12,t:10,b:28}},tiw=tw-tm.l-tm.r,mid=tm.t+(th-tm.t-tm.b)/2;
      const ts=make("svg",{{viewBox:`0 0 ${{tw}} ${{th}}`,"aria-label":`${{d.event_id}}的P和D特异变异窗口轨迹`}});
      const ymax=Math.max(1,...d.p_bins,...d.d_bins),scale=(mid-tm.t-4)/ymax,bw=tiw/d.p_bins.length*.82;
      ts.append(make("line",{{x1:tm.l,y1:mid,x2:tm.l+tiw,y2:mid,class:"axis"}}));
      d.p_bins.forEach((p,j)=>{{
        const x=tm.l+(j+.5)*tiw/d.p_bins.length;
        const dh=d.d_bins[j]*scale,ph=p*scale;
        if(dh) ts.append(make("rect",{{x:x-bw/2,y:mid-dh,width:bw,height:dh,fill:"var(--viz-series-2)"}}));
        if(ph) ts.append(make("rect",{{x:x-bw/2,y:mid,width:bw,height:ph,fill:"var(--viz-series-1)"}}));
        const hit=make("rect",{{x:x-tiw/d.p_bins.length/2,y:tm.t,width:tiw/d.p_bins.length,height:th-tm.t-tm.b,fill:"transparent"}});
        hit.append(make("title",{{}},`${{d.bin_edges[j]}}–${{d.bin_edges[j+1]-1}} bp · P=${{p}} · D=${{d.d_bins[j]}}`));
        ts.append(hit);
      }});
      [0,.25,.5,.75,1].forEach(f=>{{
        const x=tm.l+f*tiw;
        ts.append(make("text",{{x:x,y:th-7,"text-anchor":"middle",class:"minor"}},Math.round(1+f*(d.core_length-1)).toString()));
      }});
      ts.append(make("text",{{x:tm.l-7,y:mid-ymax*scale+4,"text-anchor":"end",class:"minor"}},ymax.toString()));
      ts.append(make("text",{{x:tm.l-7,y:mid+ymax*scale+4,"text-anchor":"end",class:"minor"}},(-ymax).toString()));
      wrap.append(ts); tracks.append(wrap);
    }});
  }})();
  </script>
</div>
""".strip() + "\n"
    args.html_output.parent.mkdir(parents=True, exist_ok=True)
    args.html_output.write_text(fragment)


if __name__ == "__main__":
    main()
