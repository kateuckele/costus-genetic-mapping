#!/usr/bin/env bash
# compare_all_columns.sh  FILE1  FILE2
# Prints, for every column, the set-difference of values between FILE1 and FILE2
# Assumes:
#   • Comma-separated, single-line header.
#   • No embedded commas inside fields (for fully RFC-4180–compliant CSVs, swap
#     cut -d, for csvcut –c).

set -euo pipefail

f1="$1"
f2="$2"

# How many columns?  (based on the header of the first file)
ncols=$(awk -F, 'NR==1 {print NF}' "$f1")

echo "Comparing $ncols columns between:"
echo "  • $f1"
echo "  • $f2"
echo

for ((col=1; col<=ncols; col++)); do
  echo "================  Column $col  ================"

  comm -3 \
    <(tail -n +2 "$f1" | cut -d, -f"$col" | sort -u) \
    <(tail -n +2 "$f2" | cut -d, -f"$col" | sort -u) \
    | sed 's/^/\t/'          # indent for readability

  echo
done

