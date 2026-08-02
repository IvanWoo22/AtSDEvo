#!/usr/bin/env bash
set -euo pipefail

report_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$report_dir"

pandoc REPORT_FOR_ADVISOR_PDF.md \
  --from markdown+raw_tex \
  --pdf-engine=xelatex \
  --lua-filter=cjk_linebreak.lua \
  --include-in-header=pdf_unicode_font_header.tex \
  --resource-path="$report_dir:$report_dir/.." \
  --output REPORT_FOR_ADVISOR.pdf

echo "Created: $report_dir/REPORT_FOR_ADVISOR.pdf"
