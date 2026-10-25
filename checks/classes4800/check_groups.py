"""Exact group-law checks on the plane cubics C_lam (base point O = (1:-1:0)).

Orders of the six base points; P + T2 = iota(P) and the 12-point orbit;
P = (4,9,18) on C_{155/12}: infinite order via Mazur, 3P, odd multiples positive);
the 12 rational torsion points of C_{27/2}; positive ones), isosceles points are torsion.

Run from the repository root:
    python3 checks/classes4800/check_groups.py
"""
import os
import sys
from fractions import Fraction as Fr
from itertools import permutations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dio_common import (O, F, add, neg, mul, order, normalize, is_positive, lam_of, check, finish,
                        Fail, third)  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "check_groups.txt")
log = []
BASE = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -1, 0), (0, 1, -1), (1, 0, -1)]


def perm_points(P):
    return {normalize(tuple(P[i] for i in p)) for p in permutations(range(3))}


def recip(P):
    x, y, z = P
    return normalize((y * z, x * z, x * y))


try:
    Oo = normalize(O)
    # ---------------- base points --------------------------------------------
    for lam in [Fr(27, 2), Fr(68, 5), Fr(155, 12), Fr(11), Fr(23, 2), Fr(101, 7)]:
        ords = {}
        for P in BASE:
            check(F(lam, P) == 0, f"lam={lam}: base point {P} on C", log)
            ords[P] = order(lam, P)
        log.append(f"     lam={lam}: orders of base points " + ", ".join(f"{P}:{ords[P]}" for P in BASE))
        check(sorted(ords.values()) == [1, 2, 3, 3, 6, 6], f"lam={lam}: base point orders are 1,2,3,3,6,6", log)
        check(ords[(1, -1, 0)] == 1 and ords[(0, 0, 1)] == 2, f"lam={lam}: O=(1:-1:0) order 1, T2=(0:0:1) order 2", log)
        # closure: base points form a group
        G = {normalize(P) for P in BASE}
        check(all(add(lam, P, Q) in G for P in G for Q in G), f"lam={lam}: the six base points form a subgroup (cyclic of order 6)", log)

    # ---------------- reciprocation and orbit ----------------------------------
    lam = Fr(155, 12)
    P = (4, 9, 18)
    T2 = normalize((0, 0, 1))
    pts = [mul(lam, n, P) for n in range(1, 5)] + [add(lam, P, (1, 0, -1))]
    for Q in pts:
        check(add(lam, Q, T2) == recip(Q), f"P + T2 = iota(P) for P = {Q if max(map(abs,Q))<10**6 else '(large)'}", log)
    # orbit of a point under S3 and iota = {+-P + T}
    G = [normalize(B) for B in BASE]
    for Q in pts[:2]:
        orbit = perm_points(Q) | perm_points(recip(Q))
        pmT = {add(lam, Q, T) for T in G} | {add(lam, neg(lam, Q), T) for T in G}
        check(orbit == pmT and len(orbit) == 12, "S3 x <iota> orbit of a non-torsion point = {+-P+T : T in Z/6} (12 points)", log)
    # transposition y<->z is P -> T - P with T = sigma(O) = (1:0:-1)
    Tsig = normalize((1, 0, -1))
    for Q in pts[:3]:
        sw = normalize((Q[0], Q[2], Q[1]))
        check(sw == add(lam, Tsig, neg(lam, Q)), "transposition y<->z acts as P -> (1:0:-1) - P", log)
    log.append("     => a fixed point (u:v:v) of y<->z satisfies 2P = (1:0:-1), a torsion section: every isosceles point is torsion")
    o_T = order(lam, Tsig)
    log.append(f"     order of (1:0:-1) is {o_T}")
    # check on many isosceles points
    bad = 0
    orders_seen = {}
    cnt = 0
    for v in range(1, 25):
        for u in range(1, 25):
            if u == v:
                continue
            from math import gcd
            if gcd(u, v) != 1:
                continue
            lm = lam_of((u, v, v))
            Pi = normalize((u, v, v))
            o = order(lm, Pi, bound=12)
            cnt += 1
            orders_seen[o] = orders_seen.get(o, 0) + 1
            if o is None:
                bad += 1
    check(bad == 0, f"all {cnt} primitive isosceles points (u,v,v), u!=v<=24, are torsion; orders seen {orders_seen}", log)

    # ---------------- large fibres -----------------------------------------------
    lam = Fr(155, 12)
    P = (4, 9, 18)
    check(lam_of(P) == lam and F(lam, P) == 0, "(4,9,18) lies on C_{155/12}", log)
    Q = normalize(P)
    mults = {}
    for n in range(1, 13):
        mults[n] = mul(lam, n, P)
        check(mults[n] != Oo, f"{n}P != O", log)
    log.append("     Mazur: a rational torsion point has order in {1..10,12}; nP != O for all n <= 12 => P has infinite order")
    # Compare as an ordered projective point: the curve is symmetric, so a permutation of 3P is
    # another point on it ((162833463 : 287876366 : 723926268) is (1:0:-1) - 3P), not 3P.
    check(mults[3] == normalize((162833463, 723926268, 287876366)),
          f"3P = {mults[3]} = (162833463 : 723926268 : 287876366) as an ordered point", log)
    check(mults[3] != normalize((162833463, 287876366, 723926268)),
          "the permutation (162833463 : 287876366 : 723926268) printed after Thm 8.9 is not 3P", log)
    for n in range(1, 10):
        pos = is_positive(mults[n])
        check(pos == (n % 2 == 1), f"{n}P positive: {pos} (expected {n % 2 == 1}); height ~ 10^{len(str(max(abs(c) for c in mults[n])))-1}", log)
    # odd multiples are distinct triples up to permutation
    trip = [tuple(sorted(mults[n])) for n in (1, 3, 5, 7, 9)]
    check(len(set(trip)) == 5, "P,3P,5P,7P,9P give 5 distinct triples (up to order)", log)

    # ---------------- C_{27/2} torsion ----------------------------------------------
    lam = Fr(27, 2)
    gens = [normalize(B) for B in BASE] + [normalize((1, 4, 4))]
    grp = set(gens)
    while True:
        new = {add(lam, a, b) for a in grp for b in grp} | {neg(lam, a) for a in grp}
        if new <= grp:
            break
        grp |= new
    check(len(grp) == 12, f"group generated by base points and (1:4:4) on C_27/2 has 12 elements", log)
    ords = {Pp: order(lam, Pp) for Pp in grp}
    hist = {}
    for o in ords.values():
        hist[o] = hist.get(o, 0) + 1
    log.append("     orders: " + ", ".join(f"{k}:{hist[k]}" for k in sorted(hist)))
    check(hist == {1: 1, 2: 3, 3: 2, 6: 6}, "order histogram {1:1, 2:3, 3:2, 6:6} = Z/2 x Z/6", log)
    check(ords[normalize((1, 4, 4))] == 6, f"(1:4:4) has order {ords[normalize((1,4,4))]}", log)
    pos = sorted(Pp for Pp in grp if is_positive(Pp))
    log.append("     12 points: " + ", ".join(str(Pp) for Pp in sorted(grp)))
    log.append("     positive: " + ", ".join(str(Pp) for Pp in pos))
    expected = perm_points((1, 4, 4)) | perm_points((1, 1, 4))
    check(set(pos) == expected and len(pos) == 6, "positive points = permutations of (1,4,4) and (1,1,4) (6 points)", log)
    trans = {add(lam, normalize(B), normalize((1, 4, 4))) for B in BASE}
    check(grp == {normalize(B) for B in BASE} | trans, "12 points = six base points and their translates by (1:4:4)", log)

    # ---------------- positive torsion points are isosceles or geometric progressions --------
    from math import gcd as _g
    res, n = {}, 0
    for S in range(4, 121):
        for a in range(1, S // 3 + 1):
            for b in range(a, (S - a) // 2 + 1):
                c = S - a - b
                if _g(_g(a, b), c) != 1 or a == c:
                    continue
                n += 1
                o = order(lam_of((a, b, c)), (a, b, c), 12)
                if o:
                    res.setdefault(o, []).append((a, b, c))
    other = [t for v in res.values() for t in v if not (t[0] == t[1] or t[1] == t[2] or t[1] ** 2 == t[0] * t[2])]
    log.append(f"     {n} primitive positive triples with sum <= 120 (excluding (1,1,1)); torsion ones by order: "
               + ", ".join(f"{k}:{len(v)}" for k, v in sorted(res.items())))
    check(not other, "every positive torsion point with sum <= 120 is isosceles or a geometric progression", log)
    check(all(t[1] ** 2 == t[0] * t[2] for t in res.get(12, [])), "order-12 positive points are exactly the GP triples", log)
    log.append("ALL GROUP CHECKS PASSED")
    finish(log, OUT)
except Fail:
    finish(log, OUT)
    sys.exit(1)
