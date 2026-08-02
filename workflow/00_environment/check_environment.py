#!/usr/bin/env python3
"""检查 AtSDEvo 的 Python 包和命令行工具是否可用。"""

from __future__ import annotations

import argparse
import importlib
import shutil
import subprocess
import sys


PYTHON_PACKAGES = ("numpy", "scipy", "matplotlib")
TOOLS = (
    ("mafft", "09：局部序列比对", "current"),
    ("muscle", "09：局部序列比对（必须为 MUSCLE5）", "current"),
    ("prank", "09：固定树局部序列比对", "current"),
    ("iqtree3", "11：拓扑祖先状态敏感性分析", "downstream"),
    ("blastn", "07：P0 候选映射", "downstream"),
    ("blastp", "05：外群共线性", "downstream"),
    ("makeblastdb", "05/07：BLAST 数据库", "downstream"),
    ("gffread", "05：蛋白序列提取", "downstream"),
    ("bedtools", "01–03：区间与掩蔽处理", "downstream"),
    ("samtools", "01：FASTA 索引", "downstream"),
    ("MCScanX", "05：外群共线性", "downstream"),
    ("biser", "03：片段重复发现", "full"),
)
SCOPE_LEVEL = {"current": 0, "downstream": 1, "full": 2}


def first_line(command: list[str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True, timeout=15)
    output = (result.stdout + result.stderr).strip()
    return output.splitlines()[0] if output else "无法读取版本"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scope",
        choices=tuple(SCOPE_LEVEL),
        default="current",
        help="current=当前主分析；downstream=阶段04–12；full=阶段01–12",
    )
    args = parser.parse_args()
    failures: list[str] = []
    print(f"Python\t{sys.version.split()[0]}\t{sys.executable}")
    if sys.version_info[:2] != (3, 11):
        failures.append("Python 必须使用 3.11.x")

    for package in PYTHON_PACKAGES:
        try:
            module = importlib.import_module(package)
            print(f"Python包\t{package}\t{module.__version__}")
        except (ImportError, AttributeError) as error:
            failures.append(f"缺少 Python 包 {package}: {error}")

    resolved: dict[str, str] = {}
    selected_tools = [
        (command, stage)
        for command, stage, scope in TOOLS
        if SCOPE_LEVEL[scope] <= SCOPE_LEVEL[args.scope]
    ]
    for command, stage in selected_tools:
        path = shutil.which(command)
        if path:
            resolved[command] = path
            print(f"命令\t{command}\t{path}\t{stage}")
        else:
            failures.append(f"缺少 {command}（{stage}）")

    muscle = resolved.get("muscle")
    if muscle:
        version = first_line([muscle, "-version"])
        print(f"版本\tmuscle\t{version}")
        if "muscle 5" not in version.lower():
            failures.append(f"需要 MUSCLE5，当前为：{version}")

    for command, flag in (("mafft", "--version"), ("blastn", "-version")):
        if command in resolved:
            print(f"版本\t{command}\t{first_line([resolved[command], flag])}")

    if failures:
        print("\n环境检查失败：", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print("\n环境检查通过。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
