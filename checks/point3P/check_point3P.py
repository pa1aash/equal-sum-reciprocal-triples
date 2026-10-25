"""The coordinates of 3P on C_{155/12}, P = (4:9:18), base point O = (1:-1:0).

Independent of checks/classes4800/dio_common.py: the group law is re-implemented here from
the definition (third intersection of a line with the cubic), in exact integer arithmetic.

  (a) P lies on C_L: (X+Y+Z)(XY+YZ+ZX) = L XYZ, L = 155/12
  (b) 3P = (162833463 : 723926268 : 287876366) as an ORDERED projective point
  (c) the point (162833463 : 287876366 : 723926268), with Y and Z exchanged, lies on C_L but is a
      different point: it is the image of 3P under the transposition Y <-> Z, i.e. (1:0:-1) - 3P
  (d) the same 3P from both association orders, P + (P + P) and (P + P) + P

Run from the repository root:
    python3 checks/point3P/check_point3P.py
Writes checks/point3P/check_point3P.txt.  Exits nonzero on any failure.
"""
import os
import sys
from fractions import Fraction as Fr
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "check_point3P.txt")
log = []
fails = 0


def check(cond, msg):
    global fails
    log.append(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        fails += 1


L = Fr(155, 12)


def F(P):
    x, y, z = P
    return 12 * (x + y + z) * (x * y + y * z + z * x) - 155 * x * y * z      # 12 * equation, integer


def norm(P):
    g = 0
    for c in P:
        g = gcd(g, c)
    P = tuple(c // g for c in P)
    for c in P:                     # first nonzero coordinate positive
        if c != 0:
            return P if c > 0 else tuple(-x for x in P)


def cubic_in_uv(P, Q):
    """Coefficients (a, b, c, d) of F(uP + vQ) = a u^3 + b u^2 v + c u v^2 + d v^3, exactly."""
    # evaluate at four (u, v) and solve the Vandermonde system in integers via finite differences
    vals = {}
    for (u, v) in [(1, 0), (0, 1), (1, 1), (1, -1), (1, 2)]:
        vals[(u, v)] = F(tuple(u * p + v * q for p, q in zip(P, Q)))
    a, d = vals[(1, 0)], vals[(0, 1)]
    s1 = vals[(1, 1)] - a - d          # b + c
    s2 = vals[(1, -1)] - a + d         # -b + c   (F(P - Q) = a - b + c - d)
    b = (s1 - s2) // 2
    c = (s1 + s2) // 2
    assert vals[(1, 2)] == a + 2 * b + 4 * c + 8 * d
    return a, b, c, d


def third(P, Q):
    """Third intersection with C of the line PQ (tangent line if P = Q)."""
    if norm(P) != norm(Q):
        a, b, c, d = cubic_in_uv(P, Q)
        assert a == 0 and d == 0
        return norm(tuple(c * p - b * q for p, q in zip(P, Q)))
    # tangent: gradient of F at P, pick R on the tangent line, R != P
    x, y, z = P
    h = 1
    gx = F((x + h, y, z)) - F((x - h, y, z))
    gy = F((x, y + h, z)) - F((x, y - h, z))
    gz = F((x, y, z + h)) - F((x, y, z - h))
    # central differences of a cubic: (F(x+1)-F(x-1)) = 2 F_x + (1/3) F_xxx; remove the cubic part
    def third_diff(i):
        e = [0, 0, 0]; e[i] = 1
        f = lambda k: F(tuple(p + k * ei for p, ei in zip(P, e)))
        return f(2) - 2 * f(1) + 2 * f(-1) - f(-2)          # = 2 F_iii
    g = [ (d1 - Fr(td, 6)) / 2 for d1, td in zip([gx, gy, gz], [third_diff(0), third_diff(1), third_diff(2)])]
    g = [int(v) if v.denominator == 1 else v for v in g]
    # R = cross(g, P) lies on the tangent line g.X = 0 and is not P (g.P = 0 by Euler, P x g != 0)
    R = (g[1] * z - g[2] * y, g[2] * x - g[0] * z, g[0] * y - g[1] * x)
    den = 1
    for c in R:
        if isinstance(c, Fr):
            den = den * c.denominator // gcd(den, c.denominator)
    R = tuple(int(c * den) for c in R)
    a, b, c, d = cubic_in_uv(P, R)
    assert a == 0 and b == 0
    return norm(tuple(d * p - c * r for p, r in zip(P, R)))


O = norm((1, -1, 0))


def add(P, Q):
    return third(O, third(P, Q))


P = norm((4, 9, 18))
check(F(P) == 0, "(a) P = (4:9:18) lies on C_{155/12}")
check(F(O) == 0, "(a) O = (1:-1:0) lies on C_{155/12}")
P2 = add(P, P)
P3a = add(P, P2)
P3b = add(P2, P)
check(F(P2) == 0 and F(P3a) == 0, f"(b) 2P = {P2} and 3P lie on the curve")
check(P3a == P3b, "(d) P + 2P = 2P + P")
check(P3a == (162833463, 723926268, 287876366), f"(b) 3P = {P3a} as an ordered point")
printed = (162833463, 287876366, 723926268)
check(F(printed) == 0, "(c) the printed triple lies on the curve (the curve is symmetric)")
check(norm(printed) != P3a, "(c) the printed point is not 3P")
T = norm((1, 0, -1))
OO = third(O, O)
neg = lambda X: third(OO, X)                 # -X is the third point on the line through X and O*O
check(add(P3a, neg(P3a)) == O, "(c) negation: 3P + (-3P) = O")
check(norm(printed) == add(T, neg(P3a)), "(c) the printed point is (1:0:-1) - 3P, the image of 3P under Y <-> Z")
check(norm((P3a[0], P3a[2], P3a[1])) == norm(printed), "(c) printed = 3P with Y and Z exchanged")

with open(OUT, "w") as fh:
    fh.write("\n".join(log) + f"\n\n{len(log)} checks, {fails} failures\n")
print("\n".join(log))
print(f"{len(log)} checks, {fails} failures")
sys.exit(1 if fails else 0)
