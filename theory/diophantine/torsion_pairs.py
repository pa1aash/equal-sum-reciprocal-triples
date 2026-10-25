#!/usr/bin/env python3
"""
Which primitive pairs differ by a torsion point (the note, Section 3).

    python3 torsion_pairs.py          (a few minutes; reads data/groups.csv, S <= 4800)

For a positive point P of C_lambda and a torsion point T, the point P + T is positive exactly when T lies on the
identity component of C_lambda(R) (the positive points form the other component, Lemma 2 of the note), that is,
exactly when T is not positive. The base points are never positive. In case (ii) of the torsion theorem the coset
P_0 + Z/6 consists of positive or mixed-sign points; in case (iii) the coset G + Z/6, G = (1 : r : r^2), is
positive for r > 0 and consists of points with mixed signs for r < 0 (lambda > 9 needs r + 1/r <= -4). So, besides
the base points, only the curves with torsion Z/12 and r < 0 carry torsion translations that preserve positivity.

Every pair {x, x'} of a primitive class with S <= 4800 is classified by the points P of x (one ordering) and P' of
x' (all six orderings, which realise -P' and the translates by the points of order 3):
  dual       P' - P or P' + P is T_2, (1:0:0) or (0:1:0), the coset T_2 + Z/3 (reciprocation, Proposition 3, then a
             permutation);
  torsion    P' - P or P' + P is a torsion point outside the base group;
  other      neither P' - P nor P' + P is torsion (membership in the torsion group listed by Theorem 1).
(A difference in the subgroup Z/3 of the base points would make x' a permutation of x; that is asserted not to
happen.) Asserted: every "torsion" pair lies on a curve with torsion Z/12 and r < 0; the pairs
{(165,240,440), (182,208,455)} (S = 845) and {(208,455,520), (231,336,616)} (S = 1183) on C_{169/16} (r = -4) are of
type "torsion"; for S <= 600 there are 1753 primitive pairs, 423 dual and none of type "torsion" (the note, Section 3).
Writes data/torsion_pairs.txt.
"""
import csv
import os
from fractions import Fraction
from itertools import permutations
from math import gcd, isqrt

from cubic_group import Cubic, O, lam, normalize

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def rat_sqrt(q):
    q = Fraction(q)
    if q < 0:
        return None
    a, b = isqrt(q.numerator), isqrt(q.denominator)
    return Fraction(a, b) if a * a == q.numerator and b * b == q.denominator else None


def z12_r(L):
    """r with lambda = (r + 1 + 1/r)^2, r rational, or None (two roots r, 1/r; return the one of modulus >= 1)."""
    w = rat_sqrt(L)
    if w is None:
        return None
    for s in (w, -w):                      # r + 1/r = s - 1
        d = rat_sqrt((s - 1) ** 2 - 4)
        if d is not None:
            r = ((s - 1) + d) / 2
            return r if abs(r) >= 1 else 1 / r
    return None


_TORS = {}


def torsion_set(L):
    """C_lambda(Q)_tors as an explicit set (Theorem 1 of the note): the base points, and the coset of a point
    (1 : 1 : t) or (1 : r : r^2); checked against Cubic.order."""
    if L in _TORS:
        return _TORS[L]
    C = Cubic(L)
    T = {O, (0, 1, -1), (1, 0, -1), (0, 0, 1), (1, 0, 0), (0, 1, 0)}
    gen = None
    r = z12_r(L)
    if r is not None:
        gen = normalize((r.denominator ** 2, r.numerator * r.denominator, r.numerator ** 2))
    else:
        d = rat_sqrt((L - 1) * (L - 9))
        if d is not None:
            tt = ((L - 5) + d) / 4                           # 2t^2 + (5 - L) t + 2 = 0
            gen = normalize((tt.denominator, tt.denominator, tt.numerator))
    if gen is not None:
        assert C.on(gen) and C.order(gen) in (2, 12), (L, gen)
        T |= {C.add(gen, Q) for Q in list(T)}
    assert all(C.order(Q) is not None for Q in T) and len(T) in (6, 12)
    _TORS[L] = T
    return T


def classify(x, xp):
    C = Cubic(lam(x))
    tors = torsion_set(lam(x))
    P = normalize(x)
    best = "other"
    for perm in set(permutations(xp)):
        Pp = normalize(perm)
        for D in (C.add(Pp, C.neg(P)), C.add(Pp, P)):
            if D in (O, (0, 1, -1), (1, 0, -1)):          # O and the points of order 3: a permutation
                raise AssertionError(f"{x}, {xp}: a permutation")
            if D in ((0, 0, 1), (1, 0, 0), (0, 1, 0)):     # T_2 + (order 3): a permutation of the dual
                return "dual"
            if D in tors:
                best = "torsion"
    return best


def main():
    counts = {}
    tors = []
    for row in csv.DictReader(open(os.path.join(DATA, "groups.csv"))):
        S = int(row["S"])
        members = [tuple(int(v) for v in t.strip().strip("()").split(",")) for t in row["triples"].split(";")]
        for i in range(len(members)):
            for j in range(i):
                x, xp = members[i], members[j]
                if gcd(*x, *xp) != 1:
                    continue
                kind = classify(x, xp)
                for X in (600, 4800):
                    if S <= X:
                        counts[(X, kind)] = counts.get((X, kind), 0) + 1
                if kind == "torsion":
                    r = z12_r(lam(x))
                    assert r is not None and r < 0, (x, xp, lam(x), r)
                    tors.append((S, x, xp, lam(x), r))
    tot = {X: sum(v for (Y, _), v in counts.items() if Y == X) for X in (600, 4800)}
    assert tot[600] == 1753 and counts[(600, "dual")] == 423 and counts.get((600, "torsion"), 0) == 0, counts
    ex = {(845, (165, 240, 440), (182, 208, 455)), (1183, (208, 455, 520), (231, 336, 616))}
    found = {(S, tuple(sorted(min(x, xp))), tuple(sorted(max(x, xp)))) for S, x, xp, _, _ in tors}
    assert ex <= found, (ex, found)
    out = [f"primitive pairs with S <= {X}: {tot[X]}; dual {counts.get((X, 'dual'), 0)}, torsion (outside the base group) "
           f"{counts.get((X, 'torsion'), 0)}, other {counts.get((X, 'other'), 0)}" for X in (600, 4800)]
    out.append("pairs of torsion type (all on curves with torsion Z/12 and r < 0):")
    out += [f"  S = {S}: {x} ; {xp}  lambda = {L}, r = {r}" for S, x, xp, L, r in sorted(tors)]
    open(os.path.join(DATA, "torsion_pairs.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
