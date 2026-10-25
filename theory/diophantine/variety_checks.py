#!/usr/bin/env python3
"""
Exact checks of Propositions 2-4 of the note (the pencil, its maps and base points).

  1. Scale invariance: a set of triples is a degeneracy class (equal S, equal R)
     after rescaling iff the triples share lambda = e1 e2 / e3; lambda = S R.
  2. The pencil e1 e2 = lambda e3 is Beauville's Gamma_1(6) pencil
     (X+Y)(Y+Z)(Z+X) + t XYZ = 0 with t = 1 - lambda (polynomial identity).
  3. A Weierstrass model over Q(lambda), derived here, its discriminant and
     c4, the Kodaira types I6 (inf), I3 (1), I2 (0), I1 (9), and agreement of
     its j-invariant with the Bremner-Guy-Nowakowski model
         tau^2 = sigma (sigma^2 + (n^2 - 6n - 3) sigma + 16 n).
  4. The six base points have orders 1, 2, 3, 3, 6, 6 under the chord-tangent
     law with O = (1:-1:0), at many values of lambda.
  5. Reciprocation is translation by the 2-torsion point (0:0:1), proved
     symbolically for every point of every C_lambda.
  6. The dual surface: the six cubics e2 x, e2 y, e2 z, e1 yz, e1 zx, e1 xy
     span a 5-dimensional space (one linear relation, equal sums) and their
     base locus is the three coordinate points and the two points
     e1 = e2 = 0; no three of the five are collinear. Hence the image is an
     anticanonical quartic del Pezzo surface in P^4.
  3c. The explicit Q-birational map from C_lambda to the BGN model and its
     inverse, and the integral model handed to PARI (ranks are in ranks.py).
  7. The explicit birational map between the degeneracy cone and the fibre
     square of the pencil over the lambda-line.

Usage: variety_checks.py      (writes data/variety_checks.txt)
"""

from __future__ import annotations

import random
import sys
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cubic_group import O, TRIVIAL, Cubic  # noqa: E402

OUT: list[str] = []


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def main() -> int:
    x, y, z, L, t = sp.symbols("x y z lambda t")
    e1, e2, e3 = x + y + z, x*y + y*z + z*x, x*y*z

    # 1. scale invariance
    k = sp.symbols("k", positive=True)
    lam_s = lambda a, b, c: (a + b + c) * (a*b + b*c + c*a) / (a*b*c)
    assert sp.simplify(lam_s(k*x, k*y, k*z) - lam_s(x, y, z)) == 0
    say("PASS 1  lambda = e1 e2/e3 is scale invariant; lambda = S * R")

    # 2. Beauville
    assert sp.expand((x + y) * (y + z) * (z + x) - (e1 * e2 - e3)) == 0
    say("PASS 2  (x+y)(y+z)(z+x) = e1 e2 - e3, so e1 e2 - lambda e3 = Beauville's "
        "Gamma_1(6) pencil with t = 1 - lambda")

    # 3. Weierstrass model, derived by s = (x+y)/z, W = ((x-y)/z)(lambda-1-s), s = 1/X
    X = sp.symbols("X")
    a_, b_ = L - 1, L - 5
    cub = sp.expand((1 - a_ * X) * (4 * X**2 - b_ * X + 1))
    c3, c2, c1, c0 = sp.Poly(cub, X).all_coeffs()
    A2, A4, A6 = sp.factor(c2), sp.factor(c1 * c3), sp.factor(c0 * c3**2)
    assert sp.expand(A2 - (L - 3)**2) == 0
    assert sp.expand(A4 - 8 * (L - 3) * (L - 1)) == 0
    assert sp.expand(A6 - 16 * (L - 1)**2) == 0
    b2, b4, b6, b8 = 4*A2, 2*A4, 4*A6, 4*A2*A6 - A4**2
    c4 = sp.factor(b2**2 - 24*b4)
    Delta = sp.factor(-b2**2*b8 - 8*b4**3 - 27*b6**2 + 9*b2*b4*b6)
    assert sp.expand(Delta - 4096 * L**2 * (L - 9) * (L - 1)**3) == 0
    say("PASS 3a model y^2 = x^3 + (l-3)^2 x^2 + 8(l-3)(l-1) x + 16(l-1)^2, "
        "Delta = 2^12 l^2 (l-9)(l-1)^3")
    # Kodaira types: multiplicative where ord(Delta) = n > 0 and c4 is a unit
    for place, n in [(0, 2), (1, 3), (9, 1)]:
        assert c4.subs(L, place) != 0
        assert sp.Poly(Delta, L).as_expr().subs(L, place) == 0
        say(f"        lambda = {place}: ord Delta = {n}, c4 = {c4.subs(L, place)} != 0  ->  I{n}")
    # at infinity: weights deg a_i <= i, so ord_inf Delta = 12 - deg Delta, ord_inf c4 = 4 - deg c4
    dD, dc4 = sp.degree(Delta, L), sp.degree(sp.expand(c4), L)
    assert (12 - dD, 4 - dc4) == (6, 0)
    say("        lambda = inf: ord Delta = 6, ord c4 = 0  ->  I6;  Euler numbers 6+3+2+1 = 12 "
        "(rational elliptic surface)")
    # BGN model and j-invariants
    j_ours = sp.cancel(c4**3 / Delta)
    aB, bB = L**2 - 6*L - 3, 16*L
    c4B = 16 * (aB**2 - 3 * bB)                  # for y^2 = x^3 + a x^2 + b x
    DB = 16 * bB**2 * (aB**2 - 4 * bB)
    assert sp.cancel(j_ours - c4B**3 / DB) == 0
    assert sp.factor(aB**2 - 4*bB) == sp.factor((L - 1)**3 * (L - 9))
    say("PASS 3b j-invariant equals that of the Bremner-Guy-Nowakowski model; "
        "their discriminant (n-1)^3 (n-9) reproduced")

    # 3c. explicit Q-birational map C_lambda -> BGN model, its inverse, and the PARI model
    F = e1 * e2 - L * e3
    aB_, bB_ = L**2 - 6*L - 3, 16*L
    sig = -4 * e2 / z**2
    tau = 4 * L * (L - 1) * (x - y) / (x + y - (L - 1) * z)
    W = sp.together(tau**2 - sig * (sig**2 + aB_ * sig + bB_))
    num = sp.numer(W).subs(z, 1)
    F1 = F.subs(z, 1)
    rem = sp.Poly(num, x, domain=sp.QQ.frac_field(y, L)).rem(sp.Poly(F1, x, domain=sp.QQ.frac_field(y, L)))
    assert rem.is_zero
    # inverse: s = (x+y)/z, d = (x-y)/z from (sigma, tau)
    sg, tu = sp.symbols("sigma tau")
    s_inv = sg * (L - 1) / (sg - 4 * L)
    d_inv = tu * (s_inv - L + 1) / (4 * L * (L - 1))
    # inverse o forward = identity on C_lambda
    for target, expr in [((x + y), s_inv), ((x - y), d_inv)]:
        diff = sp.together(expr.subs({sg: sig, tu: tau}).subs(z, 1) - target)
        r = sp.Poly(sp.numer(diff), x, domain=sp.QQ.frac_field(y, L)).rem(sp.Poly(F1, x, domain=sp.QQ.frac_field(y, L)))
        assert r.is_zero, target
    # forward o inverse lands on C_lambda: F(x(s,d), y(s,d), 1) vanishes on the Weierstrass curve
    xi, yi = (s_inv + d_inv) / 2, (s_inv - d_inv) / 2
    Fi = sp.together(F.subs({x: xi, y: yi, z: 1}))
    Fn = sp.Poly(sp.expand(sp.numer(Fi)), tu)
    Fn = Fn.rem(sp.Poly(tu**2 - sg * (sg**2 + aB_ * sg + bB_), tu))
    assert sp.simplify(Fn.as_expr()) == 0
    # integral model: lambda = A/B, X = B^2 sigma, Y = B^3 tau
    A_, B_ = sp.symbols("A B", positive=True)
    Xs, Ys = sp.symbols("Xs Ys")
    lhs = (B_**3 * tu)**2 - (B_**2 * sg) * ((B_**2 * sg)**2 + (A_**2 - 6*A_*B_ - 3*B_**2) * (B_**2 * sg)
                                         + 16 * A_ * B_**3)
    rhs = B_**6 * (tu**2 - sg * (sg**2 + aB_ * sg + bB_)).subs(L, A_ / B_)
    assert sp.expand(lhs - rhs) == 0
    say("PASS 3c C_lambda -> tau^2 = sigma(sigma^2 + (l^2-6l-3) sigma + 16 l) via sigma = -4 e2/z^2, "
        "tau = 4 l (l-1)(x-y)/(x+y-(l-1)z); inverse s = sigma(l-1)/(sigma-4l), "
        "x-y = tau (s-l+1)/(4l(l-1)); both compositions are the identity on the curves; "
        "X = B^2 sigma, Y = B^3 tau gives the integral model [0, A^2-6AB-3B^2, 0, 16AB^3, 0] used by ranks.py")

    # 4. torsion of the base points
    rng = random.Random(1)
    lams = [Fraction(rng.randint(10, 400), rng.randint(1, 30)) for _ in range(40)]
    lams = [l for l in lams if l not in (0, 1, 9)]
    for l in lams:
        C = Cubic(l)
        orders = [C.order(T) for T in TRIVIAL]
        assert sorted(orders) == [1, 2, 3, 3, 6, 6], (l, orders)
    say(f"PASS 4  base points have orders {{1,2,3,3,6,6}} (a Z/6 subgroup) on {len(lams)} curves")

    # 5. reciprocation IS translation by the 2-torsion point T2 = (0:0:1) (symbolic, no permutation)
    T2 = (0, 0, 1)
    Fp = lambda P: sp.expand((P[0] + P[1] + P[2]) * (P[0]*P[1] + P[1]*P[2] + P[2]*P[0]) - L * P[0]*P[1]*P[2])
    P = (x, y, z)
    Q1 = (x * z, y * z, x * y)                 # claimed third point of the line P T2
    Q2 = (y * z, x * z, x * y)                 # claimed third point of the line O Q1 = P + T2
    det = lambda U, V, Wv: sp.expand(sp.Matrix([U, V, Wv]).det())
    assert det(P, T2, Q1) == 0 and det(O, Q1, Q2) == 0
    for Q in (Q1, Q2):
        q, r_ = sp.div(Fp(Q), Fp(P), x, y, z)
        assert r_ == 0                         # Q lies on C_lambda whenever P does
    # Q1 differs generically from P and T2, Q2 from O and Q1, so they are the third intersections
    assert sp.expand(Q1[0] * P[1] - Q1[1] * P[0]) == 0 and sp.expand(Q1[0] * P[2] - Q1[2] * P[0]) != 0
    assert sp.expand(Q2[0] * Q1[1] - Q2[1] * Q1[0]) != 0
    # P + T2 = O * (P * T2) = Q2 = (yz : xz : xy) = (1/x : 1/y : 1/z); applying twice returns P,
    # so 2 T2 = O.
    say("PASS 5  P * T2 = (xz:yz:xy) and O * (P * T2) = (yz:xz:xy): with base point O = (1:-1:0), "
        "P + T2 = (1/x : 1/y : 1/z) exactly, for every P on every C_lambda; hence 2 T2 = O")

    # 6. dual surface
    cubics = [e2*x, e2*y, e2*z, e1*y*z, e1*z*x, e1*x*y]
    mons = sorted(sp.Poly(sum(cubics), x, y, z).monoms())
    allm = sorted(set(m for c in cubics for m in sp.Poly(c, x, y, z).monoms()))
    M = sp.Matrix([[sp.Poly(c, x, y, z).coeff_monomial(m) for m in allm] for c in cubics])
    assert M.rank() == 5
    assert sp.expand(sum(cubics[:3]) - sum(cubics[3:])) == 0
    w = (-1 + sp.sqrt(3) * sp.I) / 2
    pts = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, w, w**2), (1, w**2, w)]
    for P in pts:
        assert all(sp.expand(c.subs({x: P[0], y: P[1], z: P[2]})) == 0 for c in cubics)
    for P, Q, R in combinations(pts, 3):
        assert sp.expand(sp.Matrix([P, Q, R]).det()) != 0
    say("PASS 6  six cubics span dimension 5 (single relation: equal sums); base points "
        "= 3 coordinate points + {e1=e2=0}, no three collinear -> anticanonical dP4 in P^4, "
        "Pic over Q of rank 1 + 3 + 1 = 5")

    # 7. the degeneracy cone and the fibre square of the pencil, explicitly
    xp, yp, zp = sp.symbols("xp yp zp")
    f1, f2, f3 = e1, e2, e3
    g1, g2, g3 = xp + yp + zp, xp*yp + yp*zp + zp*xp, xp*yp*zp
    # Psi: (P, P') with lambda(P) = lambda(P') -> (g1 P, f1 P'): equal sums and equal R
    T1 = [g1 * x, g1 * y, g1 * z]
    T2p = [f1 * xp, f1 * yp, f1 * zp]
    E1 = lambda t: t[0] + t[1] + t[2]
    E2 = lambda t: t[0]*t[1] + t[1]*t[2] + t[2]*t[0]
    E3 = lambda t: t[0]*t[1]*t[2]
    assert sp.expand(E1(T1) - E1(T2p)) == 0
    # e2 e3' - e2' e3 of the rescaled pair equals f1^2 g1^2 (f1 f2 g3 - g1 g2 f3),
    # which vanishes exactly when lambda(P) = lambda(P')
    lhs = sp.expand(E2(T1) * E3(T2p) - E2(T2p) * E3(T1))
    assert sp.expand(lhs - f1**2 * g1**2 * (f1 * f2 * g3 - g1 * g2 * f3)) == 0
    # Phi: (t, t') on the cone -> ([t], [t']) is inverse to Psi up to the overall scale
    say("PASS 7  (P, P') on a common C_lambda  <->  (e1(P') P, e1(P) P') on the degeneracy cone: "
        "e2 e3' - e2' e3 = e1^2 e1'^2 (e1 e2 e3' - e1' e2' e3); so the cone is birational to "
        "P^2 x_{P^1} P^2 over lambda, i.e. to the fibre square of the pencil's elliptic surface")

    (HERE / "data" / "variety_checks.txt").write_text("\n".join(OUT) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
