"""Numerical checks of Theorem 3 (the isosceles family) and of the count of its copies.

A(y) = #{coprime 1 <= u < v : S(D_uv) = (2u+v)(u+2v)/g <= y}, claimed c_iso y + O(sqrt(y) log y),
c_iso = 3 log 2 / (2 pi^2).  Then L(X) = sum_D #{k >= 4 : k S_D <= X} = (c_iso + o(1)) X log X.

Run from the repository root:
    python3 checks/classes4800/check_isosceles.py
"""
import os
import sys
from math import gcd, log, pi, sqrt, isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dio_common import check, finish, Fail  # noqa: E402

OUT = os.path.join(HERE, "check_isosceles.txt")
log_ = []
c_iso = 3 * log(2) / (2 * pi ** 2)
YMAX = 4 * 10 ** 7

try:
    Svals = []
    vmax = isqrt(YMAX // 2) + 2          # (2u+v)(u+2v) >= 2v^2 and g <= 3  =>  v <= sqrt(3Y/2)
    vmax = isqrt(3 * YMAX // 2) + 2
    for v in range(2, vmax + 1):
        for u in range(1, v):
            if gcd(u, v) != 1:
                continue
            g = 3 if (u - v) % 3 == 0 else 1
            S = (2 * u + v) * (u + 2 * v) // g
            if S <= YMAX:
                Svals.append(S)
    Svals.sort()
    check(sum(1 for s in Svals if s <= 4800) == 506, "A(4800) = 506", log_)
    log_.append(f"     c_iso = {c_iso:.8f};  c_iso * 4800 = {c_iso*4800:.2f}")
    log_.append("     y            A(y)       c_iso*y      (A-c y)/(sqrt(y) log y)")
    import bisect
    worst = 0.0
    for y in [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, YMAX]:
        A = bisect.bisect_right(Svals, y)
        err = (A - c_iso * y) / (sqrt(y) * log(y))
        worst = max(worst, abs(err))
        log_.append(f"     {y:<12d} {A:<10d} {c_iso*y:<12.1f} {err:+.4f}")
    check(worst < 1.0, "normalised error |A(y) - c_iso y| / (sqrt(y) log y) stays < 1", log_)
    # partial summation: L(X) = sum_{S_D <= X/4} (floor(X/S_D) - 3)
    log_.append("     X            L(X)          L/(X log X)   (c_iso = %.5f)" % c_iso)
    ratios = []
    for X in [10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, YMAX]:
        L = sum(X // s - 3 for s in Svals if 4 * s <= X)
        ratios.append(L / (X * log(X)))
        log_.append(f"     {X:<12d} {L:<13d} {L/(X*log(X)):.5f}   (L - c X log X)/X = {(L - c_iso*X*log(X))/X:+.4f}")
    check(abs((ratios[-1] * YMAX * log(YMAX) - c_iso * YMAX * log(YMAX)) / YMAX) < 1.0,
          "L(X) - c_iso X log X = O(X) numerically (secondary term bounded)", log_)
    log_.append("ALL ISOSCELES CHECKS PASSED")
    finish(log_, OUT)
except Fail:
    finish(log_, OUT)
    sys.exit(1)
