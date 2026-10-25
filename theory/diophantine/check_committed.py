#!/usr/bin/env python3
"""Consistency of the committed degeneracy data.  Every check is an assert; nothing is
enumerated at scale and nothing is written.

The 4800-sum enumeration (enumerate_fast.py) takes hours, so the suite does not repeat it. It
checks instead that the three committed files are consistent with one another and with an
independent recount on a sample of sums:

  (a) data/groups.csv, row by row: every class has >= 2 distinct hyperbolic triples
      2 <= p <= q <= r with p + q + r = S and the stated R = 1/p + 1/q + 1/r (reduced),
      the stated multiplicity and content (gcd of all entries);
  (b) data/per_S.csv against groups.csv: classes, pairs (C(k,2) per class), maximum fibre,
      fibre histogram, primitive classes and the cumulative columns, for every S in 10..4800;
  (c) data/per_S.csv against the committed harness table data/degeneracies.csv for S <= 600
      (triads, pairs, classes, primitive classes, cumulative counts);
  (d) an independent recount from the triads themselves (triads, distinct R, classes, pairs,
      max fibre) for every S <= 150 and for S in {300, 600, 1200, 2400, 4800}.

usage: python check_committed.py
"""
import csv
import sys
from collections import Counter, defaultdict
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
S_MIN, S_MAX = 10, 4800


def triples_of(cell):
    return [tuple(int(x) for x in t.strip(" ()").split(",")) for t in cell.split(";")]


def load_groups():
    groups = defaultdict(list)          # S -> [(k, content, R)]
    n = 0
    with (HERE / "data" / "groups.csv").open() as fh:
        for r in csv.DictReader(fh):
            S, rn, rd, k, c = (int(r[x]) for x in ("S", "R_num", "R_den", "multiplicity", "content"))
            ts = triples_of(r["triples"])
            assert len(ts) == k >= 2 and len(set(ts)) == k, (S, ts)
            g = 0
            for p, q, rr in ts:
                assert 2 <= p <= q <= rr and p + q + rr == S, (S, (p, q, rr))
                num, den = q * rr + p * rr + p * q, p * q * rr
                assert den > num, ("not hyperbolic", (p, q, rr))
                d = gcd(num, den)
                assert (num // d, den // d) == (rn, rd), (S, (p, q, rr), rn, rd)
                g = gcd(g, gcd(gcd(p, q), rr))
            assert g == c, (S, ts, g, c)
            groups[S].append((k, c, (rn, rd)))
            n += 1
    return groups, n


def check_per_S(groups):
    rows = list(csv.DictReader((HERE / "data" / "per_S.csv").open()))
    assert [int(r["S"]) for r in rows] == list(range(S_MIN, S_MAX + 1))
    cp = cc = cpc = 0
    for r in rows:
        S = int(r["S"])
        gs = groups.get(S, [])
        classes = len(gs)
        pairs = sum(k * (k - 1) // 2 for k, _, _ in gs)
        prim_classes = sum(1 for _, c, _ in gs if c == 1)
        hist = Counter(k for k, _, _ in gs)
        assert int(r["classes"]) == classes and int(r["pairs"]) == pairs, S
        assert int(r["max_fibre"]) == max([k for k, _, _ in gs] + [1]), S
        assert int(r["prim_classes"]) == prim_classes and int(r["scaled_classes"]) == classes - prim_classes, S
        assert (int(r["fibres_2"]), int(r["fibres_3"]), int(r["fibres_4"])) == (hist[2], hist[3], hist[4]), S
        assert int(r["fibres_5plus"]) == sum(v for k, v in hist.items() if k >= 5), S
        cp, cc, cpc = cp + pairs, cc + classes, cpc + prim_classes
        assert (int(r["cum_pairs"]), int(r["cum_classes"]), int(r["cum_prim_classes"])) == (cp, cc, cpc), S
    assert set(groups) <= set(range(S_MIN, S_MAX + 1))
    return {int(r["S"]): r for r in rows}


def check_harness(per_s):
    rows = list(csv.DictReader((REPO / "data" / "degeneracies.csv").open()))
    assert rows and int(rows[-1]["S"]) == 600
    for r in rows:
        S = int(r["S"])
        p = per_s[S]
        assert int(r["n_triads"]) == int(p["triads"]), S
        assert int(r["N_pairs"]) == int(p["pairs"]) and int(r["N_classes"]) == int(p["classes"]), S
        assert int(r["N_primitive_classes"]) == int(p["prim_classes"]), S
        assert int(r["cum_pairs"]) == int(p["cum_pairs"]) and int(r["cum_classes"]) == int(p["cum_classes"]), S
        assert int(r["cum_primitive_classes"]) == int(p["cum_prim_classes"]), S
    return len(rows)


def recount(S):
    by_R = Counter()
    triads = 0
    for p in range(2, S // 3 + 1):
        for q in range(p, (S - p) // 2 + 1):
            r = S - p - q
            num, den = q * r + p * r + p * q, p * q * r
            if num >= den:
                continue
            triads += 1
            d = gcd(num, den)
            by_R[(num // d, den // d)] += 1
    return triads, len(by_R), sum(1 for v in by_R.values() if v > 1), \
        sum(v * (v - 1) // 2 for v in by_R.values()), max(by_R.values(), default=1)


def check_recount(per_s):
    sums = list(range(S_MIN, 151)) + [300, 600, 1200, 2400, 4800]
    for S in sums:
        t, dr, cl, pr, mf = recount(S)
        p = per_s[S]
        got = (int(p["triads"]), int(p["distinct_R"]), int(p["classes"]), int(p["pairs"]), int(p["max_fibre"]))
        assert got == (t, dr, cl, pr, mf), (S, got, (t, dr, cl, pr, mf))
    return len(sums)


def main():
    groups, n = load_groups()
    print(f"(a) groups.csv: {n} classes, every triple hyperbolic with the stated R, sum, multiplicity and content")
    per_s = check_per_S(groups)
    print(f"(b) per_S.csv: {len(per_s)} sums agree with groups.csv (classes, pairs, fibres, primitive, cumulative)")
    m = check_harness(per_s)
    print(f"(c) per_S.csv agrees with data/degeneracies.csv for all {m} sums S <= 600")
    k = check_recount(per_s)
    print(f"(d) independent recount of triads, distinct R, classes, pairs, max fibre at {k} sums: all agree")
    print("DIOPHANTINE COMMITTED-DATA CHECK PASSED")


if __name__ == "__main__":
    sys.exit(main())
