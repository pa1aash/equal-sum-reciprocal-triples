#!/usr/bin/env python3
"""Regenerate Tables 1 and 2 of the note from the committed data.

    python3 tools/make_tables.py           # write tables/fibres.tex and tables/density.tex
    python3 tools/make_tables.py --check   # exit 1 unless both files are current

Table 1 (the first class of each size k): theory/diophantine/data/groups.csv (written by
enumerate_fast.py) and the ranks and torsion groups of theory/diophantine/data/ranks.txt (ranks.py).

Table 2 (cumulative counts): pairs and classes from the S <= 6000 enumeration
(checks/enumeration6000/check_enum.txt, written by check_enum.py), cross-checked against
theory/diophantine/data/per_S.csv for S <= 4800; the local exponent log(N(S)/N(S/2))/log 2 from the
exact counts, rounded once (and checked against the transcript); the primitive classes and the largest class from
per_S.csv, and at S = 6000 from checks/note/check_note.txt (check_note.py --enum).
"""
import csv
import math
import re
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "theory" / "diophantine" / "data"
ENUM6000 = ROOT / "checks" / "enumeration6000" / "check_enum.txt"
NOTECHECK = ROOT / "checks" / "note" / "check_note.txt"
OUT = ROOT / "tables"


def tex_triple(t):
    return "(" + ",".join(str(x) for x in t) + ")"


def ranks():
    """lambda -> (proven rank, torsion group in TeX) from ranks.txt."""
    out = {}
    txt = (DATA / "ranks.txt").read_text(encoding="utf-8")
    for blk in txt.split("\nlambda = ")[1:]:
        lam = blk.split()[0]
        rk = re.search(r"rank:\s+(\d+) <= rank <= (\d+)\s+-> PROVEN", blk)
        assert rk and rk.group(1) == rk.group(2), lam
        tors = re.search(r"elltors:\s+\[(\d+), \[([\d, ]+)\]", blk)
        inv = [int(x) for x in tors.group(2).split(",")]
        assert math.prod(inv) == int(tors.group(1)), lam
        group = "\\times".join(f"\\Z/{n}" for n in sorted(inv))
        out[Fraction(lam)] = (int(rk.group(1)), f"${group}$")
    return out


def table_fibres():
    first = {}
    for r in csv.DictReader(open(DATA / "groups.csv", encoding="utf-8")):
        first.setdefault(int(r["multiplicity"]), r)
    assert sorted(first) == [2, 3, 4, 5, 6], sorted(first)
    rk = ranks()
    rows = []
    for k in sorted(first):
        r = first[k]
        S = int(r["S"])
        lam = S * Fraction(int(r["R_num"]), int(r["R_den"]))
        lam_t = str(lam.numerator) if lam.denominator == 1 else f"{lam.numerator}/{lam.denominator}"
        rank, tors = str(rk[lam][0]), rk[lam][1]
        triples = [tuple(int(x) for x in t.strip().strip("()").split(",")) for t in r["triples"].split(";")]
        mem = "; ".join("$" + tex_triple(t) + "$" for t in triples)
        rows.append(f"{k} & {S} & ${lam_t}$ & {rank} & {tors} & {mem} \\\\")
    return ("\\begin{tabular}{@{}rrrrl>{\\raggedright\\arraybackslash}p{0.44\\textwidth}@{}}\n\\toprule\n"
            "$k$ & $S$ & $\\Lambda$ & rank & torsion & members\\\\\n\\midrule\n"
            + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")


def table_density():
    per = {int(r["S"]): r for r in csv.DictReader(open(DATA / "per_S.csv", encoding="utf-8"))}
    txt = ENUM6000.read_text(encoding="utf-8")
    sec = txt.split("test to S=6000", 1)[1]
    enum = {}
    for line in sec.splitlines()[2:]:
        f = line.split()
        if len(f) != 7 or not f[0].isdigit():
            break
        enum[int(f[0])] = (int(f[1]), int(f[3]), f[5])
    note = NOTECHECK.read_text(encoding="utf-8")
    prim6000 = int(re.search(r"primitive classes with S <= 6000: (\d+)", note).group(1))
    largest6000 = int(re.search(r"no class of size 7 with S <= 6000: largest (\d+)", note).group(1))

    maxfib, cur = {}, 0
    for S in sorted(per):
        cur = max(cur, int(per[S]["max_fibre"]))
        maxfib[S] = cur
    r18 = per[18]
    rows = [f"18 & {r18['cum_pairs']} & {r18['cum_classes']} & {r18['cum_prim_classes']} & {maxfib[18]} & -- \\\\"]
    cum = {S: int(r["cum_pairs"]) for S, r in per.items()}
    for S in (100, 200, 300, 400, 500, 600, 1000, 2000, 3000, 4000, 4800, 6000):
        pairs, classes, loc4 = enum[S]
        # The local exponent log(N(S)/N(S/2))/log 2 from the exact counts, rounded once to three decimals
        # (rounding the transcript's four-decimal value again would round twice). The transcript must agree.
        exact = math.log(pairs / cum[S // 2]) / math.log(2)
        assert abs(exact - float(loc4)) <= 5e-5, (S, exact, loc4)
        loc = f"{exact:.3f}"
        if S in per:
            assert (int(per[S]["cum_pairs"]), int(per[S]["cum_classes"])) == (pairs, classes), S
            prim, mf = per[S]["cum_prim_classes"], maxfib[S]
        else:
            prim, mf = prim6000, largest6000
        rows.append(f"{S} & {pairs} & {classes} & {prim} & {mf} & {loc} \\\\")
    return ("\\begin{tabular}{@{}rrrrrr@{}}\n\\toprule\n"
            "$X$ & pairs $\\N(X)$ & classes $\\Ncl(X)$ & primitive classes & largest size & local exponent\\\\\n\\midrule\n"
            + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")


def main():
    check = "--check" in sys.argv
    stale = []
    for name, fn in (("fibres", table_fibres), ("density", table_density)):
        text, path = fn(), OUT / f"{name}.tex"
        if check:
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                stale.append(name)
        else:
            OUT.mkdir(exist_ok=True)
            path.write_text(text, encoding="utf-8")
    if check:
        if stale:
            print("out of date: " + ", ".join(f"tables/{n}.tex" for n in stale))
            return 1
        print("tables/fibres.tex and tables/density.tex are current")
    else:
        print("wrote tables/fibres.tex and tables/density.tex")
    return 0


if __name__ == "__main__":
    sys.exit(main())
