#!/usr/bin/env python3
"""Precompute fixed-tree PRANK alignments for atom FASTAs in parallel."""

from __future__ import annotations

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--prank", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=10)
    args = parser.parse_args()
    raw = args.root / "atomic_haplotypes"
    out = args.root / "alignments/PRANK_FIXED_TREE"
    out.mkdir(parents=True, exist_ok=True)
    inputs = sorted(raw.glob("*.fa"))

    def align(path: Path) -> bool:
        atom = path.stem
        prefix = out / atom
        best = Path(f"{prefix}.best.fas")
        if best.exists() and best.stat().st_size:
            return True
        command = [
            str(args.prank), f"-d={path}", f"-t={raw / f'{atom}.tree.nwk'}",
            f"-o={prefix}", "-DNA", "+F", "-once", "-showanc", "-showevents", "-quiet",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        return result.returncode == 0 and best.exists()

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        status = list(pool.map(align, inputs))
    print(f"prank_ok={sum(status)} prank_failed={len(status) - sum(status)}")


if __name__ == "__main__":
    main()
