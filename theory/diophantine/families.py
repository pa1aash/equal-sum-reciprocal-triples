#!/usr/bin/env python3
"""
Parametric families of primitive degeneracies: identities, the lower bound
they give, how much of the data they explain, and a systematic search for
further low-degree families.

Every identity is checked symbolically (sympy, for all parameter values);
every count is exact integer arithmetic; floating point appears only in the
asymptotic constants printed next to exact counts.

  1. Identities, verified for all parameters:
       Vieta     lambda(qr/p, q, r) = lambda(p, q, r),  lambda = e1 e2 / e3
       Dual      {e2(t) t, e1(t) (bc, ca, ab)} : sum e1 e2, R = 1/e3, t = (a,b,c)
       Isosceles {(2u+v)(u,v,v), (u+2v)(v,u,u)} : S = (2u+v)(u+2v), R = 1/(uv)
       Bundle    {(u(q+r)+vqr)(u, vq, vr), (u+v(q+r))(vqr, ur, uq)}
  2. Isosceles family against the data: every hyperbolic member with
     S <= S_MAX is present in data/groups.csv, and its primitive count
     tracks c_iso S with c_iso = 3 log 2 / (2 pi^2).
  3. The lower bound N(X) >= sum over isosceles D of #{k >= 4 : k S_D <= X},
     compared with floor(X/18), c_iso X log X and the true N(X).
  4. Decomposition of the primitive pairs in data/groups.csv: dual
     (reciprocal) pairs, isosceles pairs, and the rest; growth of each.
  5. Search for degree-2 families. A one-parameter family whose two triples
     are linear in the parameter before the sums are equalised is a pair of
     lines (G1, G2) in P^2 and a Moebius map phi with
     lambda o G1 = lambda o G2 o phi. Lines with coefficients up to B are
     grouped by the branch-value discriminant of lambda restricted to them,
     candidates are tested exactly, and pairs where phi is induced by a
     permutation of coordinates are discarded as trivial.
  6. Relations Q = +-m P + T (T torsion, m = 2, 3) between the two triples
     of each non-dual primitive pair, on the curve C_lambda.

Usage: families.py [S_MAX]   (default 4800; reads data/groups.csv, data/per_S.csv;
                              writes data/families.txt, data/line_families.csv)
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd, log, pi
from multiprocessing import get_context
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cubic_group import Cubic, dual, lam, normalize  # noqa: E402
from enumerate_fast import load_groups, load_per_s  # noqa: E402

DATA = HERE / "data"
C_ISO = 3 * log(2) / (2 * pi ** 2)
OUT: list[str] = []


def say(s: str = "") -> None:
    print(s, flush=True)
    OUT.append(s)


# --------------------------------------------------------------------------
# 1. identities
# --------------------------------------------------------------------------

def R_of(t):
    return sum(1 / sp.Integer(1) / m for m in t)


def check_identities() -> None:
    a, b, c, p, q, r, u, v = sp.symbols("a b c p q r u v", positive=True)
    lam_s = lambda t: (t[0] + t[1] + t[2]) * (t[0]*t[1] + t[1]*t[2] + t[2]*t[0]) / (t[0]*t[1]*t[2])
    assert sp.simplify(lam_s((q * r / p, q, r)) - lam_s((p, q, r))) == 0
    say("  PASS  Vieta: lambda(qr/p, q, r) = lambda(p, q, r)")

    e1, e2, e3 = a + b + c, a*b + b*c + c*a, a*b*c
    A = [e2 * a, e2 * b, e2 * c]
    B = [e1 * b * c, e1 * c * a, e1 * a * b]
    assert sp.expand(sum(A) - e1 * e2) == 0 and sp.expand(sum(B) - e1 * e2) == 0
    assert sp.simplify(sum(1 / x for x in A) - 1 / e3) == 0
    assert sp.simplify(sum(1 / x for x in B) - 1 / e3) == 0
    say("  PASS  Dual: {e2 t, e1 (bc,ca,ab)} has S = e1 e2 and R = 1/e3 for all (a,b,c)")

    A = [(2*u + v) * u, (2*u + v) * v, (2*u + v) * v]
    B = [(u + 2*v) * v, (u + 2*v) * u, (u + 2*v) * u]
    assert sp.expand(sum(A) - (2*u + v) * (u + 2*v)) == 0
    assert sp.expand(sum(B) - (2*u + v) * (u + 2*v)) == 0
    assert sp.simplify(sum(1 / x for x in A) - 1 / (u * v)) == 0
    assert sp.simplify(sum(1 / x for x in B) - 1 / (u * v)) == 0
    say("  PASS  Isosceles: {(2u+v)(u,v,v), (u+2v)(v,u,u)} has S = (2u+v)(u+2v), R = 1/(uv)")

    al, be = q + r, q * r
    A = [(u*al + v*be) * u, (u*al + v*be) * v * q, (u*al + v*be) * v * r]
    B = [(u + v*al) * v * be, (u + v*al) * u * r, (u + v*al) * u * q]
    S = (u*al + v*be) * (u + v*al)
    assert sp.expand(sum(A) - S) == 0 and sp.expand(sum(B) - S) == 0
    assert sp.simplify(sum(1 / x for x in A) - sum(1 / x for x in B)) == 0
    say("  PASS  Conic bundle: {(u(q+r)+vqr)(u,vq,vr), (u+v(q+r))(vqr,ur,uq)} has equal S and R")


# --------------------------------------------------------------------------
# 2-3. isosceles family
# --------------------------------------------------------------------------

def isosceles_primitive(s_max: int) -> list[tuple[int, tuple, tuple]]:
    """Primitive isosceles dual pairs with S <= s_max (u < v or u > v, unordered)."""
    out = []
    u = 1
    while (2 * u + 1) * (u + 2) <= 3 * s_max:
        v = 1
        while (2 * u + v) * (u + 2 * v) <= 3 * s_max:
            if u < v and gcd(u, v) == 1:
                A = ((2*u + v) * u, (2*u + v) * v, (2*u + v) * v)
                B = ((u + 2*v) * v, (u + 2*v) * u, (u + 2*v) * u)
                g = gcd(gcd(*A), gcd(*B))
                A = tuple(sorted(x // g for x in A))
                B = tuple(sorted(x // g for x in B))
                S = sum(A)
                if S <= s_max:
                    out.append((S, A, B))
            v += 1
        u += 1
    return sorted(out)


def hyperbolic(t) -> bool:
    return min(t) >= 2 and Fraction(1, t[0]) + Fraction(1, t[1]) + Fraction(1, t[2]) < 1


def check_isosceles(groups, per_s, s_max: int) -> None:
    pairs_in_data = set()
    for S, _, _, _, _, T in groups:
        for t1, t2 in combinations(T, 2):
            pairs_in_data.add((S, t1, t2))
            pairs_in_data.add((S, t2, t1))
    prim = isosceles_primitive(s_max)
    for S, A, B in prim:
        assert sum(A) == sum(B) == S
        assert Fraction(1, A[0]) + Fraction(1, A[1]) + Fraction(1, A[2]) == \
               Fraction(1, B[0]) + Fraction(1, B[1]) + Fraction(1, B[2])
    # every hyperbolic scaled member must be in the enumeration
    members = 0
    for S0, A, B in prim:
        k = 1
        while k * S0 <= s_max:
            Ak = tuple(k * x for x in A)
            Bk = tuple(k * x for x in B)
            if hyperbolic(Ak) and hyperbolic(Bk):
                assert (k * S0, Ak, Bk) in pairs_in_data, (k * S0, Ak, Bk)
                members += 1
            k += 1
    say(f"  PASS  all {members} hyperbolic isosceles pairs with S <= {s_max} are in data/groups.csv")

    by = {r["S"]: r for r in per_s}
    say(f"  {'X':>5s} {'prim iso':>9s} {'c_iso X':>8s} {'iso pairs k>=4':>15s} {'c_iso XlogX':>12s} "
        f"{'floor(X/18)':>11s} {'N(X)':>7s}")
    for X in [600, 1000, 2000, 3000, 4000, 4800]:
        if X > s_max:
            continue
        p = sum(1 for S, _, _ in prim if S <= X)
        lb = sum(max(0, X // S - 3) for S, _, _ in prim)   # k = 4 .. floor(X/S)
        say(f"  {X:5d} {p:9d} {C_ISO * X:8.1f} {lb:15d} {C_ISO * X * log(X):12.1f} "
            f"{X // 18:11d} {by[X]['cum_pairs']:7d}")
        assert lb <= by[X]["cum_pairs"]
        assert lb > X // 18 or X < 200


# --------------------------------------------------------------------------
# 4. decomposition of the data
# --------------------------------------------------------------------------

def is_isosceles(t) -> bool:
    return t[0] == t[1] or t[1] == t[2]


def decompose(groups, s_max: int) -> list[tuple]:
    cnt = defaultdict(lambda: [0, 0, 0, 0])   # X -> prim pairs, dual, iso dual, non-dual
    nondual = []
    events = []
    for S, _, _, _, _, T in groups:
        for t1, t2 in combinations(T, 2):
            if gcd(gcd(*t1), gcd(*t2)) != 1:
                continue
            d = dual(t1) == tuple(sorted(normalize(t2)))
            iso = d and is_isosceles(t1)
            events.append((S, d, iso))
            if not d:
                nondual.append((S, t1, t2))
    say(f"  {'X':>5s} {'prim pairs':>10s} {'dual':>7s} {'iso dual':>8s} {'non-dual':>8s} {'dual %':>7s}")
    rows = []
    for X in [200, 600, 1000, 2000, 3000, 4000, 4800]:
        if X > s_max:
            continue
        tot = sum(1 for S, _, _ in events if S <= X)
        du = sum(1 for S, d, _ in events if S <= X and d)
        iso = sum(1 for S, _, i in events if S <= X and i)
        rows.append((X, tot, du, iso, tot - du))
        say(f"  {X:5d} {tot:10d} {du:7d} {iso:8d} {tot - du:8d} {100 * du / tot:7.1f}")
    for label, col in [("dual", 2), ("non-dual", 4), ("all primitive", 1)]:
        a, b = rows[-3], rows[-1]
        e = log(b[col] / a[col]) / log(b[0] / a[0])
        say(f"  local exponent of primitive {label} pairs over {a[0]}..{b[0]}: {e:.3f}")
    return nondual


# --------------------------------------------------------------------------
# 5. degree-2 families: pairs of lines
# --------------------------------------------------------------------------

s_, t_, L_ = sp.symbols("s t L")


def line_points(ell):
    """Two integer points spanning the line a x + b y + c z = 0."""
    M = sp.Matrix([list(ell)])
    ns = M.nullspace()
    basis = []
    for vec in ns:
        den = sp.ilcm(*[sp.fraction(x)[1] for x in vec])
        basis.append([int(x * den) for x in vec])
    return basis[0], basis[1]


def restricted_map(ell):
    """lambda restricted to the line, as reduced binary forms (N, D) in (s, t)."""
    P0, P1 = line_points(ell)
    X = [P0[i] * s_ + P1[i] * t_ for i in range(3)]
    e1 = X[0] + X[1] + X[2]
    e2 = X[0] * X[1] + X[1] * X[2] + X[2] * X[0]
    e3 = X[0] * X[1] * X[2]
    N, D = sp.expand(e1 * e2), sp.expand(e3)
    g = sp.gcd(N, D)
    N, D = sp.cancel(N / g), sp.cancel(D / g)
    return (P0, P1), sp.expand(N), sp.expand(D)


def branch_signature(N, D):
    """Canonical primitive integer polynomial whose roots are the branch values."""
    F = sp.Poly(sp.expand(N - L_ * D), s_, t_)
    deg = F.total_degree()
    coeffs = [F.coeff_monomial(s_ ** (deg - i) * t_ ** i) for i in range(deg + 1)]
    if deg == 3:
        a, b, c, d = coeffs
        disc = b*b*c*c - 4*a*c**3 - 4*b**3*d - 27*a*a*d*d + 18*a*b*c*d
    elif deg == 2:
        a, b, c = coeffs
        disc = b*b - 4*a*c
    else:
        return None
    P = sp.Poly(sp.expand(disc), L_)
    if P.is_zero:
        return None
    P = P.primitive()[1]
    if P.LC() < 0:
        P = -P
    return deg, tuple(P.all_coeffs())


def admissible(ell) -> bool:
    """The line meets the open positive triangle."""
    return any(x > 0 for x in ell) and any(x < 0 for x in ell)


def moebius_equivalences(m1, m2):
    """All Moebius phi (over Q) with lambda o G1 = lambda o G2 o phi, found exactly."""
    (_, N1, D1), (_, N2, D2) = m1, m2
    f1 = lambda S, T: Fraction(int(N1.subs({s_: S, t_: T})), 1) / int(D1.subs({s_: S, t_: T})) \
        if D1.subs({s_: S, t_: T}) != 0 else None
    # three rational source points and the rational points of their target fibres
    srcs, fibres = [], []
    for S0, T0 in [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (1, 2)]:
        val = f1(S0, T0)
        if val is None:
            continue
        poly = sp.Poly(sp.expand(N2 * val.denominator - D2 * val.numerator), s_, t_)
        roots = []
        # rational roots of the binary form: t/s = x, plus s = 0
        if poly.as_expr().subs(s_, 0) == 0:
            roots.append((0, 1))
        uni = sp.Poly(poly.as_expr().subs(s_, 1), t_)
        if not uni.is_zero:
            for x in sp.roots(uni, filter="Q").keys():
                x = sp.Rational(x)
                roots.append((x.q, x.p))
        if roots:
            srcs.append((S0, T0))
            fibres.append(roots)
        if len(srcs) == 3:
            break
    if len(srcs) < 3:
        return []
    found = []
    for images in product(*fibres):
        # Moebius (s,t) -> M (s,t) with M src_i ~ img_i
        a, b, c, d = sp.symbols("a b c d")
        k1, k2, k3 = sp.symbols("k1 k2 k3")
        eqs = []
        for (S0, T0), (S1, T1), k in zip(srcs, images, (1, k2, k3)):
            eqs += [a * S0 + b * T0 - k * S1, c * S0 + d * T0 - k * T1]
        sol = sp.solve(eqs, [a, b, c, d, k2, k3], dict=True)
        if not sol:
            continue
        sol = sol[0]
        M = [sol.get(a, a), sol.get(b, b), sol.get(c, c), sol.get(d, d)]
        if any(x.free_symbols for x in M):
            continue
        if M[0] * M[3] - M[1] * M[2] == 0:
            continue
        Ss, Ts = M[0] * s_ + M[1] * t_, M[2] * s_ + M[3] * t_
        lhs = sp.expand(N1 * D2.subs({s_: Ss, t_: Ts}, simultaneous=True)
                        - N2.subs({s_: Ss, t_: Ts}, simultaneous=True) * D1)
        if lhs == 0:
            found.append(tuple(M))
    return found


def is_trivial(m1, m2, M) -> bool:
    """phi induced by a coordinate permutation: G2(phi(.)) is a permutation of G1(.)."""
    (P0, P1), _, _ = m1
    (Q0, Q1), _, _ = m2
    a, b, c, d = M
    G1 = [P0[i] * s_ + P1[i] * t_ for i in range(3)]
    G2 = [Q0[i] * (a * s_ + b * t_) + Q1[i] * (c * s_ + d * t_) for i in range(3)]
    for sig in permutations(range(3)):
        H = [G2[sig[i]] for i in range(3)]
        if all(sp.expand(G1[i] * H[j] - G1[j] * H[i]) == 0 for i in range(3) for j in range(3)):
            return True
    return False


def _map_job(ell):
    m = restricted_map(ell)
    return ell, m, branch_signature(m[1], m[2])


def search_line_pairs(B: int) -> list[dict]:
    lines = set()
    for ell in product(range(-B, B + 1), repeat=3):
        if gcd(gcd(*ell), 0) != 1 or not admissible(ell):
            continue
        ell = normalize(ell)
        # one representative per S3-orbit (permuted lines give permuted families)
        lines.add(min(normalize(tuple(ell[i] for i in sg)) for sg in permutations(range(3))))
    lines = sorted(lines)
    with get_context("spawn").Pool() as pool:
        results = pool.map(_map_job, lines, chunksize=16)
    buckets = defaultdict(list)
    for ell, m, sig in results:
        if sig is not None:
            buckets[sig].append((ell, m))
    say(f"  {len(lines)} lines up to S3 with |coefficients| <= {B}; "
        f"{len(buckets)} distinct branch signatures")
    families = []
    for sig, members in buckets.items():
        # every line against every permutation of every line in the bucket (itself included)
        for i, (l1, m1) in enumerate(members):
            for l2, _ in members[i:]:
                pass_pair(l1, m1, l2, families)
    return dedupe(families)


def pass_pair(l1, m1, l2, families) -> None:
    for sg in set(permutations(range(3))):
        l2p = normalize(tuple(l2[i] for i in sg))
        m2 = restricted_map(l2p)
        for M in moebius_equivalences(m1, m2):
            if is_trivial(m1, m2, M):
                continue
            families.append(dict(line1=l1, line2=l2p, phi=M,
                                 degree=sp.Poly(m1[2], s_, t_).total_degree(),
                                 G1=m1[0], G2=m2[0]))


def dedupe(families: list[dict]) -> list[dict]:
    """One family per unordered pair of S3-orbits of lines."""
    canon = lambda ell: min(normalize(tuple(ell[i] for i in sg)) for sg in permutations(range(3)))
    uniq, seen = [], set()
    for f in families:
        a, b = canon(f["line1"]), canon(f["line2"])
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)
        uniq.append(f)
    return uniq


def family_pairs(f, s_max: int):
    """Primitive hyperbolic degeneracies produced by a line-pair family, S <= s_max."""
    (P0, P1), (Q0, Q1) = f["G1"], f["G2"]
    a, b, c, d = [sp.Rational(x) for x in f["phi"]]
    den = sp.ilcm(*[x.q for x in (a, b, c, d)])
    a, b, c, d = [int(x * den) for x in (a, b, c, d)]
    out = set()
    lim = int(s_max ** 0.5) * 4 + 4
    for S0 in range(-lim, lim + 1):
        for T0 in range(0, lim + 1):
            if gcd(S0, T0) != 1:
                continue
            x = [P0[i] * S0 + P1[i] * T0 for i in range(3)]
            y = [Q0[i] * (a * S0 + b * T0) + Q1[i] * (c * S0 + d * T0) for i in range(3)]
            for sx in (1, -1):
                for sy in (1, -1):
                    X = [sx * e for e in x]
                    Y = [sy * e for e in y]
                    if min(X) <= 0 or min(Y) <= 0:
                        continue
                    X, Y = tuple(sorted(normalize(X))), tuple(sorted(normalize(Y)))
                    if X == Y or lam(X) != lam(Y):
                        continue
                    sx_, sy_ = sum(X), sum(Y)
                    Lc = sx_ * sy_ // gcd(sx_, sy_)
                    if Lc > s_max:
                        continue
                    A = tuple(m * (Lc // sx_) for m in X)
                    Bt = tuple(m * (Lc // sy_) for m in Y)
                    if hyperbolic(A) and hyperbolic(Bt):
                        out.add((Lc,) + tuple(sorted((A, Bt))))
    return out


# --------------------------------------------------------------------------
# 6. multiplication relations on non-dual pairs
# --------------------------------------------------------------------------

def _mult_job(item):
    S, t1, t2 = item
    C = Cubic(lam(t1))
    P = tuple(t1)
    hits = []
    for Q in set(permutations(t2)):
        for m in (2, 3):
            for base, other in ((P, Q), (Q, P)):
                mB = C.mul(m, base)
                for sign in (1, -1):
                    Dv = C.add(other, C.neg(mB)) if sign == 1 else C.add(other, mB)
                    if C.order(Dv) is not None:
                        hits.append(m)
    return S, sorted(set(hits))


def main(argv: list[str]) -> int:
    s_max = int(argv[0]) if argv else 4800
    groups = [g for g in load_groups() if g[0] <= s_max]
    per_s = load_per_s()

    say("== 1. identities (symbolic, all parameter values) ==")
    check_identities()

    say("\n== 2-3. isosceles family and the lower bound ==")
    check_isosceles(groups, per_s, s_max)

    say("\n== 4. decomposition of primitive pairs ==")
    nondual = decompose(groups, s_max)

    say("\n== 5. degree-2 families from pairs of lines ==")
    fams = search_line_pairs(B=int(argv[1]) if len(argv) > 1 else 6)
    data_pairs = set()
    for S, _, _, _, _, T in groups:
        for t1, t2 in combinations(T, 2):
            data_pairs.add((S,) + tuple(sorted((t1, t2))))
    with (DATA / "line_families.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["line1", "line2", "phi", "map_degree", "pairs_S_le_max", "in_data", "dual_family"])
        for f in fams:
            pairs = family_pairs(f, s_max)
            assert pairs <= data_pairs, "a family member is missing from the enumeration"
            prim = [p for p in pairs if gcd(gcd(*p[1]), gcd(*p[2])) == 1]
            is_dual = all(dual(p[1]) == tuple(sorted(normalize(p[2]))) for p in prim) if prim else None
            w.writerow([f["line1"], f["line2"], f["phi"], f["degree"], len(pairs), len(pairs), is_dual])
            say(f"  family {f['line1']} ~ {f['line2']}  phi={f['phi']}  deg {f['degree']}: "
                f"{len(pairs)} pairs with S <= {s_max} (all in data), "
                f"{'reciprocal (dual) pairs' if is_dual else 'NOT the dual family'}")
    say(f"  {len(fams)} non-trivial line-pair families found")

    say("\n== 6. multiplication relations Q = +-mP + T on non-dual pairs (S <= 400) ==")
    items = [x for x in nondual if x[0] <= 400]
    with get_context("spawn").Pool() as pool:
        res = pool.map(_mult_job, items, chunksize=8)
    hit = sum(1 for _, h in res if h)
    say(f"  {len(items)} non-dual primitive pairs; {hit} satisfy Q = +-mP + T or P = +-mQ + T "
        f"with m in (2, 3) and T torsion")

    (DATA / "families.txt").write_text("\n".join(OUT) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
