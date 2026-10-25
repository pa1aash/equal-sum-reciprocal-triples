"""Shared exact-arithmetic helpers for the checks in this folder.

Everything here is independent of theory/diophantine/.  The curve is

    C_lam :  (x+y+z)(xy+yz+zx) = lam * x*y*z            (projective plane cubic)

with group law (chord and tangent) taken with base point O = (1:-1:0).
The Weierstrass model is

    E_lam :  W^2 = s (s^2 + (lam^2 - 6 lam - 3) s + 16 lam),
    s = -4 e2 / z^2 .

All arithmetic is in Fraction / int.
"""
from fractions import Fraction as Fr
from math import gcd
import sys

O = (1, -1, 0)


def F(lam, P):
    x, y, z = P
    return (x + y + z) * (x * y + y * z + z * x) - lam * x * y * z


def grad(lam, P):
    x, y, z = P
    e1 = x + y + z
    e2 = x * y + y * z + z * x
    return (e2 + e1 * (y + z) - lam * y * z,
            e2 + e1 * (x + z) - lam * x * z,
            e2 + e1 * (x + y) - lam * x * y)


def normalize(P):
    """Projective normal form: primitive integers, first nonzero coordinate > 0."""
    P = [Fr(c) for c in P]
    if all(c == 0 for c in P):
        raise ValueError("zero vector")
    den = 1
    for c in P:
        den = den * c.denominator // gcd(den, c.denominator)
    Q = [int(c * den) for c in P]
    g = 0
    for c in Q:
        g = gcd(g, abs(c))
    Q = [c // g for c in Q]
    for c in Q:
        if c != 0:
            if c < 0:
                Q = [-d for d in Q]
            break
    return tuple(Q)


def _cubic_coeffs(lam, P, Q):
    """Coefficients (c3,c2,c1,c0) of F(r P + s Q) = c3 r^3 + c2 r^2 s + c1 r s^2 + c0 s^3."""
    def val(r, s):
        return F(lam, tuple(r * p + s * q for p, q in zip(P, Q)))
    c3 = val(1, 0)
    c0 = val(0, 1)
    f1 = val(1, 1)    # c3+c2+c1+c0
    fm = val(1, -1)   # c3-c2+c1-c0
    # c2 + c1 = f1 - c3 - c0 ;  -c2 + c1 = fm - c3 + c0
    s1 = f1 - c3 - c0
    s2 = fm - c3 + c0
    c1 = Fr(s1 + s2, 2)
    c2 = Fr(s1 - s2, 2)
    return c3, c2, c1, c0


def third(lam, P, Q):
    """Third intersection of the line PQ (tangent if P == Q) with C_lam."""
    P = normalize(P)
    Q = normalize(Q)
    assert F(lam, P) == 0 and F(lam, Q) == 0, (P, Q)
    if P != Q:
        c3, c2, c1, c0 = _cubic_coeffs(lam, P, Q)
        assert c3 == 0 and c0 == 0
        # F = r s (c2 r + c1 s): third root (r:s) = (c1 : -c2)
        R = tuple(c1 * p - c2 * q for p, q in zip(P, Q))
        if all(c == 0 for c in R):
            raise ArithmeticError("line contained in curve?")
        return normalize(R)
    g = grad(lam, P)
    assert any(c != 0 for c in g), "singular point"
    # a point D != P on the tangent line g . X = 0
    cands = [(g[1], -g[0], 0), (g[2], 0, -g[0]), (0, g[2], -g[1])]
    D = None
    for cand in cands:
        if any(c != 0 for c in cand):
            Dn = normalize(cand)
            if Dn != P:
                D = Dn
                break
    assert D is not None
    c3, c2, c1, c0 = _cubic_coeffs(lam, P, D)
    assert c3 == 0 and c2 == 0
    # F = s^2 (c1 r + c0 s)
    if c1 == 0:
        return P  # flex
    R = tuple(c0 * p - c1 * d for p, d in zip(P, D))
    return normalize(R)


def add(lam, P, Q):
    return third(lam, O, third(lam, P, Q))


def neg(lam, P):
    OO = third(lam, O, O)
    return third(lam, P, OO)


def mul(lam, n, P):
    if n < 0:
        return mul(lam, -n, neg(lam, P))
    R = normalize(O)
    A = normalize(P)
    while n:
        if n & 1:
            R = add(lam, R, A)
        A = add(lam, A, A)
        n >>= 1
    return R


def order(lam, P, bound=24):
    Q = normalize(P)
    for n in range(1, bound + 1):
        if Q == normalize(O):
            return n
        Q = add(lam, Q, P)
    return None


def is_positive(P):
    P = normalize(P)
    return all(c > 0 for c in P) or all(c < 0 for c in P)


def lam_of(t):
    a, b, c = t
    e1 = a + b + c
    e2 = a * b + b * c + c * a
    e3 = a * b * c
    return Fr(e1 * e2, e3)


def esym(t):
    a, b, c = t
    return a + b + c, a * b + b * c + c * a, a * b * c


def to_weier(lam, P):
    """(x:y:z) on C_lam with z != 0 and generic -> (s, W) on E_lam (affine), or None for O."""
    x, y, z = [Fr(c) for c in P]
    if z == 0:
        raise ValueError("z = 0: use group law instead")
    x, y = x / z, y / z
    s = -4 * (x * y + x + y)
    W = -4 * (-lam * x * y - lam * x + lam * y + 2 * x * y * y + 3 * x * y + x + 2 * y * y + y)
    return s, W


def on_weier(lam, s, W):
    return W * W == s * (s * s + (lam * lam - 6 * lam - 3) * s + 16 * lam)


class Fail(Exception):
    pass


def check(cond, msg, log):
    log.append(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        raise Fail(msg)


def finish(log, path):
    with open(path, "w") as fh:
        fh.write("\n".join(log) + "\n")
    print("\n".join(log))
