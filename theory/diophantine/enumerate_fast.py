#!/usr/bin/env python3
"""
Fast exact enumeration of two-coefficient degeneracies, 10 <= S <= S_MAX.

A degeneracy class at sum S is a set of k >= 2 distinct hyperbolic triples
2 <= p <= q <= r with p + q + r = S sharing R = 1/p + 1/q + 1/r. The heavy loop
lives in enumerate_core.c (integer arithmetic, reduced-fraction keys); this
driver compiles it, runs it in parallel over blocks of S, and writes

    data/per_S.csv    one row per S: triads, distinct R-values, pairs, classes,
                      maximum fibre, primitive / scaled split, fibre-size
                      histogram, cumulative counts under both conventions
    data/groups.csv   every degeneracy class, with R, content and triples

Conventions: N_pairs counts C(k,2) per class,
N_classes counts 1 per class. A class is primitive when the gcd of every entry
of every triple in it is 1; a pair is primitive when the gcd of its six entries
is 1.

Before anything is written the output is asserted against the existing
harness, code/enumerate_degeneracies.py, for every S <= 600: per-S triad,
pair, class and primitive-class counts, and the full list of classes (R,
triples, content) must match exactly, and the committed data/degeneracy-
groups.csv must be reproduced row for row.

Usage: enumerate_fast.py [S_MAX]   (default 4800; the kernel's key packing is
exact up to S = 4800 and refuses anything larger)
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
DATA = HERE / "data"
CORE_SRC = HERE / "enumerate_core.c"

S_MIN = 10
S_MAX_DEFAULT = 4800
CHECK_MAX = 600


def build_core(workdir: Path) -> Path:
    exe = workdir / "enumerate_core"
    subprocess.run(
        ["cc", "-O3", "-Wall", "-Wextra", "-Werror", "-o", str(exe), str(CORE_SRC)],
        check=True,
    )
    return exe


def blocks(lo: int, hi: int, n: int) -> list[tuple[int, int]]:
    """Split [lo, hi] into about n blocks of roughly equal cost (cost ~ S^2)."""
    total = sum(s * s for s in range(lo, hi + 1))
    target = total / n
    out, start, acc = [], lo, 0
    for s in range(lo, hi + 1):
        acc += s * s
        if acc >= target and s < hi:
            out.append((start, s))
            start, acc = s + 1, 0
    out.append((start, hi))
    return out


def run_core(exe: Path, lo: int, hi: int) -> str:
    res = subprocess.run([str(exe), "range", str(lo), str(hi)],
                         check=True, capture_output=True, text=True)
    return res.stdout


def parse(text: str, per_s: dict, groups: list) -> None:
    for line in text.splitlines():
        f = line.split("\t")
        if f[0] == "S":
            S, triads, distinct, pairs, classes, maxk, pc, pp = map(int, f[1:])
            per_s[S] = dict(S=S, triads=triads, distinct_R=distinct, pairs=pairs,
                            classes=classes, max_fibre=maxk,
                            prim_classes=pc, prim_pairs=pp)
        elif f[0] == "G":
            S, num, den, k, content = map(int, f[1:6])
            triples = tuple(tuple(int(x) for x in t.split(",")) for t in f[6].split(";"))
            assert len(triples) == k
            groups.append((S, num, den, k, content, triples))
        else:
            raise RuntimeError(f"unparsed kernel output: {line!r}")


def verify_groups(groups: list) -> None:
    """Independent exact re-check of every class: Fractions, not the kernel."""
    for S, num, den, k, content, triples in groups:
        R = Fraction(num, den)
        assert len(set(triples)) == k >= 2
        g = 0
        for t in triples:
            p, q, r = t
            assert 2 <= p <= q <= r and p + q + r == S, (S, t)
            assert Fraction(1, p) + Fraction(1, q) + Fraction(1, r) == R, (S, t, R)
            assert R < 1
            for m in t:
                g = gcd(g, m)
        assert g == content, (S, triples, g, content)


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


def check_against_harness(per_s: dict, groups: list) -> None:
    """Assert exact agreement with code/enumerate_degeneracies.py for S <= 600."""
    sys.path.insert(0, str(REPO / "code"))
    import enumerate_degeneracies as H  # noqa: E402

    records, detail = H.sweep(S_MIN, CHECK_MAX)
    mine: dict[int, list] = {}
    for g in groups:
        if g[0] <= CHECK_MAX:
            mine.setdefault(g[0], []).append(g)
    for rec in records:
        S = rec.S
        row = per_s[S]
        assert row["triads"] == rec.n_triads, (S, "triads", row["triads"], rec.n_triads)
        assert row["pairs"] == rec.n_pairs, (S, "pairs")
        assert row["classes"] == rec.n_classes, (S, "classes")
        assert row["prim_classes"] == rec.n_primitive_classes, (S, "primitive classes")
        theirs = sorted((g.R, g.triples, g.content) for g in detail[S])
        ours = sorted((Fraction(n, d), t, c) for (_, n, d, _, c, t) in mine.get(S, []))
        assert ours == theirs, (S, ours, theirs)
    # cumulative checkpoints printed in the paper, both conventions
    for S, (want_p, want_c) in H.CHECKPOINTS.items():
        assert sum(per_s[s]["pairs"] for s in range(S_MIN, S + 1)) == want_p, S
        assert sum(per_s[s]["classes"] for s in range(S_MIN, S + 1)) == want_c, S

    # and the committed groups file, row for row
    committed = []
    with (REPO / "data" / "degeneracy-groups.csv").open() as fh:
        for r in csv.DictReader(fh):
            trip = tuple(tuple(int(x) for x in t.strip(" ()").split(","))
                         for t in r["triples"].split(";"))
            committed.append((int(r["S"]), Fraction(int(r["R_num"]), int(r["R_den"])),
                              trip, int(r["content"])))
    ours_all = [(S, Fraction(n, d), t, c) for (S, n, d, _, c, t) in groups if S <= CHECK_MAX]
    assert sorted(committed) == sorted(ours_all), "data/degeneracy-groups.csv not reproduced"
    print(f"  PASS  harness agreement for every S <= {CHECK_MAX}: "
          f"{len(records)} sums, {len(ours_all)} classes, checkpoints in both conventions")


def finish_rows(per_s: dict, groups: list, s_max: int) -> list[dict]:
    hist: dict[int, dict[int, int]] = {}
    for S, _, _, k, _, _ in groups:
        hist.setdefault(S, {}).setdefault(k, 0)
        hist[S][k] += 1
    rows = []
    cum = dict(pairs=0, classes=0, prim_classes=0, prim_pairs=0)
    for S in range(S_MIN, s_max + 1):
        row = dict(per_s[S])
        h = hist.get(S, {})
        assert sum(h.values()) == row["classes"]
        assert sum(k * (k - 1) // 2 * c for k, c in h.items()) == row["pairs"]
        row["scaled_classes"] = row["classes"] - row["prim_classes"]
        row["scaled_pairs"] = row["pairs"] - row["prim_pairs"]
        row["fibres_2"] = h.get(2, 0)
        row["fibres_3"] = h.get(3, 0)
        row["fibres_4"] = h.get(4, 0)
        row["fibres_5plus"] = sum(c for k, c in h.items() if k >= 5)
        for key in cum:
            cum[key] += row[key]
            row["cum_" + key] = cum[key]
        rows.append(row)
    return rows


COLUMNS = ["S", "triads", "distinct_R", "pairs", "classes", "max_fibre",
           "prim_pairs", "scaled_pairs", "prim_classes", "scaled_classes",
           "fibres_2", "fibres_3", "fibres_4", "fibres_5plus",
           "cum_pairs", "cum_classes", "cum_prim_pairs", "cum_prim_classes"]


def write_outputs(rows: list[dict], groups: list) -> None:
    DATA.mkdir(exist_ok=True)
    with (DATA / "per_S.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r[c] for c in COLUMNS})
    with (DATA / "groups.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["S", "R_num", "R_den", "multiplicity", "content", "triples"])
        for S, n, d, k, c, t in sorted(groups):
            w.writerow([S, n, d, k, c, " ; ".join("(%d,%d,%d)" % x for x in t)])


def load_groups(path: Path = DATA / "groups.csv") -> list:
    """Read data/groups.csv back as (S, num, den, k, content, triples) tuples."""
    out = []
    with path.open() as fh:
        for r in csv.DictReader(fh):
            trip = tuple(tuple(int(x) for x in t.strip(" ()").split(","))
                         for t in r["triples"].split(";"))
            out.append((int(r["S"]), int(r["R_num"]), int(r["R_den"]),
                        int(r["multiplicity"]), int(r["content"]), trip))
    return out


def load_per_s(path: Path = DATA / "per_S.csv") -> list[dict]:
    with path.open() as fh:
        return [{k: int(v) for k, v in r.items()} for r in csv.DictReader(fh)]


def main(argv: list[str]) -> int:
    s_max = int(argv[0]) if argv else S_MAX_DEFAULT
    if s_max < CHECK_MAX:
        raise SystemExit(f"S_MAX must be >= {CHECK_MAX} so the harness check runs")
    with tempfile.TemporaryDirectory() as tmp:
        exe = build_core(Path(tmp))
        jobs = blocks(S_MIN, s_max, 4 * (os.cpu_count() or 4))
        with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as pool:
            outputs = list(pool.map(lambda b: run_core(exe, *b), jobs))
    per_s: dict = {}
    groups: list = []
    for text in outputs:
        parse(text, per_s, groups)
    assert sorted(per_s) == list(range(S_MIN, s_max + 1)), "missing sums"

    print(f"enumerated {S_MIN} <= S <= {s_max}: {len(groups)} classes")
    verify_groups(groups)
    print("  PASS  every class re-verified in Fraction arithmetic")
    check_against_harness(per_s, groups)
    rows = finish_rows(per_s, groups, s_max)
    write_outputs(rows, groups)

    by = {r["S"]: r for r in rows}
    print("\n     S   cum_pairs  cum_classes  cum_prim_cls  N/S^2     max_fibre<=S")
    running_max = 0
    for r in rows:
        running_max = max(running_max, r["max_fibre"])
        r["running_max"] = running_max
    for S in [100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000,
              2500, 3000, 3500, 4000, 4500, 4800]:
        if S in by:
            r = by[S]
            print(f"  {S:5d}  {r['cum_pairs']:9d}  {r['cum_classes']:11d}  "
                  f"{r['cum_prim_classes']:12d}  {r['cum_pairs'] / S**2:.5f}  {r['running_max']:4d}")
    for k in range(3, 10):
        first = next((r["S"] for r in rows if r["max_fibre"] >= k), None)
        if first is None:
            print(f"  no fibre of size {k} up to S = {s_max}")
            break
        print(f"  first fibre of size >= {k}: S = {first}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
