#!/usr/bin/env python3
"""Section 6.1 of the note: "Among the 582 sums 19 <= S <= 600, 507 carry a primitive class."

Counts, from theory/diophantine/data/per_S.csv (written by enumerate_fast.py), the sums S with
19 <= S <= 600 and those among them with at least one primitive class. Exits nonzero on a mismatch.

    python3 checks/note/count_sums.py
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
rows = [r for r in csv.DictReader(open(ROOT / "theory" / "diophantine" / "data" / "per_S.csv", encoding="utf-8"))
        if 19 <= int(r["S"]) <= 600]
n_sums = len(rows)
n_prim = sum(1 for r in rows if int(r["prim_classes"]) > 0)
ok = (n_sums, n_prim) == (582, 507)
print(f"{'PASS' if ok else 'FAIL'} sums 19 <= S <= 600: {n_sums}; with a primitive class: {n_prim} (stated: 582, 507)")
sys.exit(0 if ok else 1)
