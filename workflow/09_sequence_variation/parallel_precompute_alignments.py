#!/usr/bin/env python3
"""Precompute MAFFT L-INS-i and MUSCLE5 atom alignments in parallel."""

from __future__ import annotations

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--mafft", required=True, type=Path)
    parser.add_argument("--muscle", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    raw = args.root / "atomic_haplotypes"
    mafft_dir = args.root / "alignments/MAFFT_LINSI"
    muscle_dir = args.root / "alignments/MUSCLE5"
    mafft_dir.mkdir(parents=True, exist_ok=True)
    muscle_dir.mkdir(parents=True, exist_ok=True)
    inputs = sorted(path for path in raw.glob("*.fa") if not path.name.endswith(".tree.nwk"))

    def align(path: Path) -> str:
        mafft_out = mafft_dir / path.name
        muscle_out = muscle_dir / path.name
        if not mafft_out.exists() or not mafft_out.stat().st_size:
            result = subprocess.run(
                [str(args.mafft), "--localpair", "--maxiterate", "1000", "--quiet", str(path)],
                check=True, capture_output=True, text=True,
            )
            mafft_out.write_text(result.stdout)
        if not muscle_out.exists() or not muscle_out.stat().st_size:
            subprocess.run(
                [str(args.muscle), "-align", str(path), "-output", str(muscle_out)],
                check=True, capture_output=True, text=True,
            )
        return path.name

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        completed = sum(1 for _ in pool.map(align, inputs))
    print(f"aligned_atoms={completed}")


if __name__ == "__main__":
    main()
