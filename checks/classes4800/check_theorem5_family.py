"""Check of Lemma 3 and Theorem 4 (the conic bundle of dual pairs) against the family of the proof.

Family: coprime q <= r, coprime u, v >= 1, t = (u, vq, vr), alpha = q+r, beta = qr.
Dual pair / v = {(alpha u + beta v)(u, vq, vr), (u + alpha v)(beta v, ur, uq)}, common sum
Q(u,v) = (alpha u + beta v)(u + alpha v).

Checks: the identity; the multiplicity (tuples per reduced pair D) on an exhaustive range;
the box inequality; the constant, symbolically; the phi-sum asymptotics; the uniformity
exponent; and the family count against c_A y log y.

Run from the repository root:
    python3 checks/classes4800/check_theorem5_family.py
"""
import os
import sys
from fractions import Fraction as Fr
from math import gcd, log, pi, sqrt

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dio_common import check, finish, Fail  # noqa: E402

OUT = os.path.join(HERE, "check_theorem5_family.txt")
log_ = []


def reduced_pair(q, r, u, v):
    al, be = q + r, q * r
    A = [(al * u + be * v) * x for x in (u, v * q, v * r)]
    B = [(u + al * v) * x for x in (be * v, u * r, u * q)]
    g = 0
    for x in A + B:
        g = gcd(g, x)
    A = tuple(sorted(x // g for x in A))
    B = tuple(sorted(x // g for x in B))
    return (min(A, B), max(A, B)), sum(A)


def is_gp(t):
    a, b, c = sorted(t)
    return b * b == a * c


try:
    # 1. identity
    q, r, u, v = sp.symbols("q r u v", positive=True)
    al, be = q + r, q * r
    A = [(al * u + be * v) * x for x in (u, v * q, v * r)]
    B = [(u + al * v) * x for x in (be * v, u * r, u * q)]
    Qf = (al * u + be * v) * (u + al * v)
    check(sp.expand(sum(A) - Qf) == 0 and sp.expand(sum(B) - Qf) == 0, "bundle pair: both sums = (au+bv)(u+av)", log_)
    check(sp.simplify(sum(1 / x for x in A) - sum(1 / x for x in B)) == 0, "bundle pair: equal reciprocal sums", log_)
    t = (u, v * q, v * r)
    e1, e2 = sum(t), t[0] * t[1] + t[1] * t[2] + t[2] * t[0]
    check(sp.expand(e2 - v * (al * u + be * v)) == 0 and sp.expand(e1 - (u + al * v)) == 0,
          "it is the dual pair of t=(u,vq,vr) divided by v (e2 = v(au+bv), e1 = u+av)", log_)

    # 2. multiplicity, exhaustive for Q <= Y
    Y = 60000
    seen = {}
    ntup = 0
    for qq in range(1, Y):
        if (2 * qq) * qq * qq > Y and qq > 1:  # alpha*beta lower bound for r >= q
            break
        for rr in range(qq, Y):
            if gcd(qq, rr) != 1:
                continue
            a_, b_ = qq + rr, qq * rr
            if (a_ + b_) * (1 + a_) > Y:
                break
            for vv in range(1, Y):
                if (a_ + b_ * vv) * (1 + a_ * vv) > Y:
                    break
                for uu in range(1, Y):
                    Qv = (a_ * uu + b_ * vv) * (uu + a_ * vv)
                    if Qv > Y:
                        break
                    if gcd(uu, vv) != 1:
                        continue
                    if is_gp((uu, vv * qq, vv * rr)):
                        continue
                    D, S = reduced_pair(qq, rr, uu, vv)
                    seen.setdefault(D, []).append((qq, rr, uu, vv))
                    ntup += 1
    mult = {}
    for D, L in seen.items():
        mult[len(L)] = mult.get(len(L), 0) + 1
    log_.append(f"     tuples (q,r,u,v) with Q <= {Y}, non-GP: {ntup}; distinct reduced pairs D: {len(seen)}; "
                f"histogram of tuples per D: {dict(sorted(mult.items()))}")
    check(max(mult) <= 6, "every reduced pair D arises from at most 6 tuples (multiplicity 6 is correct)", log_)
    # tuples with both members producing the same D: the second member is again of the form (u', v'q', v'r')
    both = sum(1 for L in seen.values() if len({(x[0], x[1]) for x in L}) > 1)
    log_.append(f"     pairs D reached from two different slopes (q,r): {both}")

    # 3. box inequality: beta/alpha <= alpha and Q <= 4 alpha beta V^2 on [1,U]x[1,V]
    for (qq, rr) in [(1, 1), (1, 2), (2, 3), (3, 7), (5, 8)]:
        a_, b_ = qq + rr, qq * rr
        V = Fr(10 ** 3)
        U = b_ * V / a_
        corner = (a_ * U + b_ * V) * (U + a_ * V)
        check(Fr(b_, a_) <= a_ and corner <= 4 * a_ * b_ * V * V, f"box corner bound for (q,r)=({qq},{rr})", [])
    log_.append("PASS box [1,U]x[1,V], U = beta V/alpha, V = sqrt(y/(4 alpha beta)): Q <= y (corner check, beta/alpha <= alpha)")

    # 4. constant, symbolically
    P = sp.pi
    box_density = sp.Rational(6) / P ** 2 * sp.Rational(1, 4)          # (6/pi^2) U V = (6/pi^2) y / (4 alpha^2)
    phi_sum = sp.Rational(3) / P ** 2                                  # sum_{n<=M} phi(n)/(2 n^2) ~ (3/pi^2) log M
    c_A = box_density * phi_sum * sp.Rational(1, 6) * sp.Rational(1, 8)   # / multiplicity, log M = log y / 8
    check(sp.simplify(c_A - sp.Rational(3, 32) / P ** 4) == 0, f"c_A = (3/(2pi^2))(3/pi^2)(1/6)(1/8) = {c_A}", log_)
    c_N = c_A * sp.Rational(1, 2) * sp.Rational(1, 2)                  # X/(2 S_D) and (c_A/2)(log)^2
    check(sp.simplify(c_N - sp.Rational(3, 128) / P ** 4) == 0, f"c_N = c_A * (1/2) * (1/2) = {c_N} (as stated)", log_)
    theta = sp.Rational(1, 6)
    best = box_density * phi_sum * sp.Rational(1, 6) * theta * sp.Rational(1, 2)
    log_.append(f"     with the same argument, log M = theta log y for any theta < 1/6, and #{{k>=4}} ~ X/S_D on S_D <= X/log X:"
                f" constant -> {sp.simplify(best)} (sup over the method); 3/(128 pi^4) is valid, not optimal")

    # 5. phi-sum asymptotics and the lower bound over q <= r <= M
    for M in [100, 1000, 3000]:
        phi = list(range(M + 1))
        for p in range(2, M + 1):
            if phi[p] == p:
                for k in range(p, M + 1, p):
                    phi[k] -= phi[k] // p
        s_phi = sum(phi[n] / (2 * n * n) for n in range(3, M + 1)) + 1 / 4  # n = 2: the single pair q=r=1
        s_pairs = sum(1 / (qq + rr) ** 2 for rr in range(1, M + 1) for qq in range(1, rr + 1) if gcd(qq, rr) == 1) \
            if M <= 1000 else None
        msg = f"     M={M}: sum_(n<=M) phi(n)/(2n^2) = {s_phi:.5f}, (3/pi^2) log M = {3 / pi ** 2 * log(M):.5f}"
        if s_pairs is not None:
            msg += f", sum over coprime q<=r<=M of 1/(q+r)^2 = {s_pairs:.5f}"
            check(s_pairs >= s_phi - 1e-12, f"M={M}: sum over q<=r<=M dominates the n<=M phi-sum", [])
        log_.append(msg)
    log_.append("PASS sum over coprime q<=r<=M of 1/(q+r)^2 >= sum_(n<=M) phi(n)/(2n^2) ~ (3/pi^2) log M")

    # 6. uniformity: error O(M^3 sqrt(y) log y) = o(y) iff 3 theta + 1/2 < 1
    check(3 * Fr(1, 8) + Fr(1, 2) < 1, "M = y^(1/8): total error y^(7/8) log y = o(y) (needs theta < 1/6)", log_)

    # 7. family count against c_A y log y
    cA = 3 / (32 * pi ** 4)
    Svals = sorted(S for D, S in [(D, sum(D[0])) for D in seen])
    for yy in [10000, 30000, 60000]:
        n = sum(1 for s in Svals if s <= yy)
        log_.append(f"     y={yy}: distinct family pairs with S_D <= y: {n}; c_A y log y = {cA * yy * log(yy):.0f}")
        check(n >= cA * yy * log(yy), f"y={yy}: family count exceeds c_A y log y", [])
    log_.append("PASS family counts exceed c_A y log y on the computed range")
    log_.append("ALL THEOREM-5 FAMILY CHECKS PASSED")
    finish(log_, OUT)
except Fail:
    finish(log_, OUT)
    sys.exit(1)
