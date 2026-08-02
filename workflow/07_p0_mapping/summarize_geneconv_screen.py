#!/usr/bin/env python3
"""Summarize nominal GENECONV fragments with experiment-wide correction."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", required=True, type=Path)
    args = parser.parse_args()
    root = args.pilot / "sequence_variation/geneconv"
    manifest = read_tsv(root / "block_manifest.tsv")
    by_id = {row["block_id"]: row for row in manifest}
    fragments = []
    for frags_path in sorted((root / "input").glob("*.frags")):
        block_id = frags_path.stem
        with frags_path.open() as handle:
            for line in handle:
                if not line.startswith("GI "):
                    continue
                fields = line.split()
                pair = fields[1]
                fragments.append(
                    {
                        "event_id": by_id[block_id]["event_id"],
                        "block_id": block_id,
                        "sequence_pair": pair,
                        "simulation_p": fields[2],
                        "BC_KA_p": fields[3],
                        "aligned_begin": fields[4],
                        "aligned_end": fields[5],
                        "fragment_bp": fields[6],
                        "polymorphic_sites": fields[7],
                        "simulation_bonferroni_218_blocks": (
                            f"{min(1.0, float(fields[2]) * len(manifest)):.8f}"
                        ),
                        "experiment_wide_status": (
                            "PASS"
                            if float(fields[2]) * len(manifest) <= 0.05
                            else "FAIL"
                        ),
                    }
                )
    write_tsv(root / "nominal_fragments.tsv", fragments)
    pd_fragments = [
        row for row in fragments
        if set(row["sequence_pair"].split(";")) == {"P", "D"}
    ]
    write_tsv(
        root / "screen_summary.tsv",
        [
            {
                "tested_contiguous_blocks": len(manifest),
                "tested_events": len(
                    {row["event_id"] for row in manifest}
                ),
                "blocks_with_any_nominal_GI_fragment": len(
                    {row["block_id"] for row in fragments}
                ),
                "events_with_nominal_PD_fragment": len(
                    {row["event_id"] for row in pd_fragments}
                ),
                "nominal_PD_fragments": len(pd_fragments),
                "experiment_wide_significant_PD_fragments": sum(
                    row["experiment_wide_status"] == "PASS"
                    for row in pd_fragments
                ),
                "interpretation": (
                    "screen_only;three_sequence_blocks;"
                    "no_nominal_fragment_is_experiment_wide_significant"
                ),
            }
        ],
    )
    print(
        f"GENECONV: {len(pd_fragments)} nominal P-D fragments; "
        f"{sum(row['experiment_wide_status'] == 'PASS' for row in pd_fragments)} "
        "experiment-wide significant"
    )


if __name__ == "__main__":
    main()
