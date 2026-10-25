#!/usr/bin/env python3
"""
Exact enumeration of hyperbolic n-cone spherical orbifolds.

An n-cone pillow O(m_1, ..., m_n) is the 2-sphere carrying n cone points of
orders 2 <= m_1 <= ... <= m_n. It is hyperbolic when

    sum_i (1 - 1/m_i) > 2,

equivalently when its orbifold Euler characteristic is negative. For n = 3 this
reduces exactly to the paper's condition 1/p + 1/q + 1/r < 1.

Everything here is exact: integers and fractions.Fraction only, never a float.
The module is stdlib-only by design -- it is the arithmetic floor the rest of
the harness stands on, so it must not inherit a third-party dependency.

Written to be reused unchanged at n = 4 and n = 5. Nothing is specialised to
triples: no triple is hardcoded anywhere in this file.

AREA CONVENTION
---------------
Gauss-Bonnet for a closed 2-orbifold gives

    Area(O) = -2*pi*chi_orb(O),        chi_orb(O) = 2 - sum_i (1 - 1/m_i).

For n = 3 this is Area = 2*pi*(1 - R) with R = sum_i 1/m_i, which is the
paper's area formula. The factor 2 matters: Area = pi*(1 - R) would be wrong by a factor of two,
and `selftest()` below asserts against that error directly.

Area is returned as an exact rational multiple of pi -- the coefficient of pi,
not a float -- so `area_over_pi` is the honest name and no rounding enters.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import Iterable, Iterator, Sequence

__all__ = [
    "Pillow",
    "is_hyperbolic",
    "reciprocal_sum",
    "power_sum",
    "euler_characteristic",
    "area_over_pi",
    "enumerate_pillows",
    "enumerate_by_sum",
    "signature",
    "signature_exact",
    "selftest",
]

Orders = tuple[int, ...]


# --------------------------------------------------------------------------
# Core invariants. Each takes a plain tuple of cone orders and is valid for
# every n >= 1.
# --------------------------------------------------------------------------

def reciprocal_sum(orders: Sequence[int]) -> Fraction:
    """R = sum_i 1/m_i, exact."""
    total = Fraction(0)
    for m in orders:
        total += Fraction(1, m)
    return total


def order_sum(orders: Sequence[int]) -> int:
    """S_1 = sum_i m_i."""
    return sum(orders)


def power_sum(orders: Sequence[int], k: int) -> int:
    """P_k = sum_i m_i**k, for arbitrary k >= 0 (k = 3, 5, 7 are used later)."""
    if k < 0:
        raise ValueError(f"power_sum expects k >= 0, got {k}")
    return sum(m ** k for m in orders)


def euler_characteristic(orders: Sequence[int]) -> Fraction:
    """
    chi_orb = 2 - sum_i (1 - 1/m_i) = 2 - n + R.

    For n = 3 this is R - 1, the paper's chi(O).
    """
    n = len(orders)
    return Fraction(2 - n) + reciprocal_sum(orders)


def area_over_pi(orders: Sequence[int]) -> Fraction:
    """
    Area / pi, exact.

    Area = -2*pi*chi_orb, so this returns -2*chi_orb. For n = 3 that is
    2*(1 - R), i.e. Area = 2*pi*(1 - R), the paper's eq:area.
    """
    return -2 * euler_characteristic(orders)


def is_hyperbolic(orders: Sequence[int]) -> bool:
    """
    sum_i (1 - 1/m_i) > 2, i.e. chi_orb < 0, equivalently R < n - 2.

    Integer-only and gcd-free: clearing denominators by prod(m_i) turns the
    test into a comparison of two integers. This sits in the innermost loop of
    every enumeration, so it avoids building Fraction objects; it is exactly as
    rigorous, since no division is ever performed.
    """
    n = len(orders)
    if n <= 2:
        return False  # no closed hyperbolic sphere with fewer than 3 cones
    prod = 1
    for m in orders:
        prod *= m
    lhs = 0
    for m in orders:
        lhs += prod // m
    return lhs < (n - 2) * prod


def signature_exact(orders: Sequence[int]) -> tuple[int, Fraction]:
    """
    The two-coefficient signature sigma = (S_1, R) as exact objects.

    This is the invariant the first two heat coefficients resolve. Keying a
    degeneracy search on anything else -- in particular on any invariant that
    omits S_1 -- does not search for what the paper defines.
    """
    return order_sum(orders), reciprocal_sum(orders)


def signature(orders: Sequence[int]) -> tuple[int, int, int]:
    """
    sigma as a hashable integer triple (S_1, num, den) with num/den = R in
    lowest terms.

    Same information as `signature_exact`, but built from machine integers
    rather than Fraction objects. Used in the hot enumeration path, where it is
    several times faster and every bit as exact.
    """
    s1 = order_sum(orders)
    num, den = _reciprocal_sum_reduced(orders)
    return s1, num, den


def _reciprocal_sum_reduced(orders: Sequence[int]) -> tuple[int, int]:
    """R as a reduced (numerator, denominator) pair, using integers only."""
    num, den = 0, 1
    for m in orders:
        num, den = num * m + den, den * m
        g = gcd(num, den)
        if g > 1:
            num //= g
            den //= g
    return num, den


# --------------------------------------------------------------------------
# Pillow record
# --------------------------------------------------------------------------

class Pillow:
    """One hyperbolic n-cone pillow, with its invariants computed on demand."""

    __slots__ = ("orders",)

    def __init__(self, orders: Iterable[int]) -> None:
        t = tuple(orders)
        if len(t) < 1:
            raise ValueError("a pillow needs at least one cone order")
        if any(m < 2 for m in t):
            raise ValueError(f"cone orders must be >= 2, got {t}")
        if list(t) != sorted(t):
            raise ValueError(f"cone orders must be non-decreasing, got {t}")
        self.orders = t

    # -- invariants -------------------------------------------------------
    @property
    def n(self) -> int:
        return len(self.orders)

    @property
    def S1(self) -> int:
        return order_sum(self.orders)

    @property
    def R(self) -> Fraction:
        return reciprocal_sum(self.orders)

    @property
    def chi(self) -> Fraction:
        return euler_characteristic(self.orders)

    @property
    def area_over_pi(self) -> Fraction:
        return area_over_pi(self.orders)

    def P(self, k: int) -> int:
        return power_sum(self.orders, k)

    @property
    def signature(self) -> tuple[int, int, int]:
        return signature(self.orders)

    @property
    def signature_exact(self) -> tuple[int, Fraction]:
        return signature_exact(self.orders)

    def is_hyperbolic(self) -> bool:
        return is_hyperbolic(self.orders)

    def scaled(self, k: int) -> "Pillow":
        """The pillow with every cone order multiplied by k."""
        if k < 1:
            raise ValueError(f"scale factor must be >= 1, got {k}")
        return Pillow(m * k for m in self.orders)

    @property
    def content(self) -> int:
        """gcd of the cone orders. content > 1 means this is a scaled pillow."""
        g = 0
        for m in self.orders:
            g = gcd(g, m)
        return g

    # -- plumbing ---------------------------------------------------------
    def __eq__(self, other: object) -> bool:
        return isinstance(other, Pillow) and self.orders == other.orders

    def __hash__(self) -> int:
        return hash(self.orders)

    def __lt__(self, other: "Pillow") -> bool:
        return self.orders < other.orders

    def __repr__(self) -> str:
        return f"O{self.orders}"


# --------------------------------------------------------------------------
# Enumeration
# --------------------------------------------------------------------------

def enumerate_by_sum(n: int, S: int) -> Iterator[Orders]:
    """
    Yield every multiset 2 <= m_1 <= ... <= m_n with sum S that is hyperbolic,
    in lexicographic order, as plain tuples.

    General in n. Recursive with an interval bound on each coordinate, so no
    non-hyperbolic branch is walked further than it has to be.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if S < 2 * n:
        return

    orders: list[int] = []

    def walk(depth: int, low: int, remaining: int) -> Iterator[Orders]:
        slots = n - depth
        if slots == 1:
            m = remaining
            if m >= low:
                orders.append(m)
                cand = tuple(orders)
                if is_hyperbolic(cand):
                    yield cand
                orders.pop()
            return
        # m_depth ranges over [low, floor(remaining / slots)]: it may not
        # exceed the average of what is left, or the sequence could not stay
        # non-decreasing.
        high = remaining // slots
        for m in range(low, high + 1):
            orders.append(m)
            yield from walk(depth + 1, m, remaining - m)
            orders.pop()

    yield from walk(0, 2, S)


def enumerate_pillows(n: int, S: int) -> list[Pillow]:
    """`enumerate_by_sum` wrapped as Pillow records."""
    return [Pillow(t) for t in enumerate_by_sum(n, S)]


def min_hyperbolic_sum(n: int, search_ceiling: int = 10_000) -> int:
    """
    Smallest cone-order sum admitting a hyperbolic n-cone pillow.

    Computed, not tabulated. For n = 3 this returns 10, realised only by
    O(3,3,4) -- the paper's claim at the head of Section 3.2.
    """
    for S in range(2 * n, search_ceiling + 1):
        for _ in enumerate_by_sum(n, S):
            return S
    raise RuntimeError(
        f"no hyperbolic {n}-cone pillow with sum <= {search_ceiling}; "
        "raise search_ceiling"
    )


# --------------------------------------------------------------------------
# Self-test. Run this file directly to execute it.
# --------------------------------------------------------------------------

def selftest() -> None:
    """
    Regression tests for this module. Every check is an assert: this function
    either returns having proved its claims, or raises.
    """
    checks = 0

    def note(label: str) -> None:
        nonlocal checks
        checks += 1
        print(f"  PASS  {label}")

    # -- area convention, the factor-of-2 regression --------------------
    # eq:area -- Area = 2*pi*(1 - R) for n = 3. The superseded upstream script
    # used pi*(1 - R). Assert the correct value AND assert we are not
    # reproducing the wrong one.
    for orders in [(2, 8, 8), (3, 3, 12), (3, 3, 4), (2, 3, 7), (7, 11, 13)]:
        R = reciprocal_sum(orders)
        assert area_over_pi(orders) == 2 * (1 - R), (
            f"area convention broken at {orders}: "
            f"{area_over_pi(orders)} != {2 * (1 - R)}"
        )
        assert area_over_pi(orders) != (1 - R), (
            f"area at {orders} equals the superseded pi*(1-R) form; "
            "the factor-of-2 regression is back"
        )
    note("eq:area -- Area/pi = 2*(1 - R) for n = 3, and not the pi*(1-R) form")

    # Area is positive exactly on hyperbolic pillows, for every n tested.
    for n in (3, 4, 5):
        for S in range(2 * n, 2 * n + 25):
            for orders in enumerate_by_sum(n, S):
                assert area_over_pi(orders) > 0, f"non-positive area at {orders}"
                assert euler_characteristic(orders) < 0, f"chi >= 0 at {orders}"
    note("Area > 0 and chi_orb < 0 on every enumerated pillow, n = 3, 4, 5")

    # -- chi and the n = 3 reduction ------------------------------------
    for orders in [(2, 8, 8), (3, 3, 12), (5, 6, 7), (4, 4, 10)]:
        assert euler_characteristic(orders) == reciprocal_sum(orders) - 1, (
            f"chi != R - 1 at {orders}"
        )
    note("chi_orb = R - 1 at n = 3")

    # -- hyperbolicity reduces to the paper's condition at n = 3 ----
    for S in range(6, 60):
        general = {t for t in enumerate_by_sum(3, S)}
        manuscript = set()
        for p in range(2, S // 3 + 1):
            for q in range(p, (S - p) // 2 + 1):
                r = S - p - q
                if r >= q and Fraction(1, p) + Fraction(1, q) + Fraction(1, r) < 1:
                    manuscript.add((p, q, r))
        assert general == manuscript, (
            f"n=3 reduction mismatch at S={S}: "
            f"only-general={general - manuscript}, only-manuscript={manuscript - general}"
        )
    note("n = 3 hyperbolicity reduces exactly to 1/p + 1/q + 1/r < 1, S <= 59")

    # -- enumeration is ordered, in range, and complete ------------------
    for n in (3, 4, 5):
        for S in range(2 * n, 2 * n + 20):
            seen = list(enumerate_by_sum(n, S))
            assert len(seen) == len(set(seen)), f"duplicates at n={n}, S={S}"
            assert seen == sorted(seen), f"not lexicographic at n={n}, S={S}"
            for t in seen:
                assert len(t) == n, f"wrong arity {t}"
                assert sum(t) == S, f"wrong sum {t}"
                assert list(t) == sorted(t), f"unsorted {t}"
                assert t[0] >= 2, f"order < 2 in {t}"
    note("enumeration: no duplicates, sorted, correct arity and sum, n = 3, 4, 5")

    # Brute force agreement, n = 4, small S: the recursive bound must not skip.
    from itertools import combinations_with_replacement
    for S in range(8, 30):
        brute = {
            t for t in combinations_with_replacement(range(2, S + 1), 4)
            if sum(t) == S and is_hyperbolic(t)
        }
        assert set(enumerate_by_sum(4, S)) == brute, f"n=4 brute-force mismatch at S={S}"
    note("n = 4 agrees with brute force over 8 <= S <= 29")

    # -- least hyperbolic sum, computed not tabulated --------------------
    assert min_hyperbolic_sum(3) == 10, "least hyperbolic triple sum is not 10"
    assert [t for t in enumerate_by_sum(3, 10)] == [(3, 3, 4)], (
        "S = 10 is not realised uniquely by O(3,3,4)"
    )
    for S in range(6, 10):
        assert not list(enumerate_by_sum(3, S)), f"unexpected hyperbolic triple at S={S}"
    note("least hyperbolic n=3 sum is 10, realised only by O(3,3,4)")

    # -- power sums ------------------------------------------------------
    assert power_sum((2, 8, 8), 3) == 1032, "P_3(2,8,8) != 1032"
    assert power_sum((3, 3, 12), 3) == 1782, "P_3(3,3,12) != 1782"
    assert power_sum((2, 8, 8), 1) == 18 == power_sum((3, 3, 12), 1)
    assert power_sum((2, 8, 8), 0) == 3, "P_0 should count the cone points"
    for k in (5, 7):
        assert power_sum((2, 8, 8), k) == 2 ** k + 2 * 8 ** k, f"P_{k} wrong"
    note("P_k for k = 0, 1, 3, 5, 7, including P_3 = 1032 and 1782")

    # -- signature agreement between the fast and exact forms -----------
    for n in (3, 4):
        for S in range(2 * n, 2 * n + 30):
            for orders in enumerate_by_sum(n, S):
                s1, num, den = signature(orders)
                s1e, Re = signature_exact(orders)
                assert s1 == s1e, f"S1 disagreement at {orders}"
                assert Fraction(num, den) == Re, f"R disagreement at {orders}"
                assert gcd(num, den) == 1, f"R not reduced at {orders}"
    note("fast integer signature agrees with the exact Fraction signature")

    # -- scaling ---------------------------------------------------------
    base = Pillow((2, 8, 8))
    for k in range(1, 15):
        sc = base.scaled(k)
        assert sc.S1 == k * base.S1, f"scaled S1 wrong at k={k}"
        assert sc.R == base.R / k, f"scaled R wrong at k={k}"
        assert sc.is_hyperbolic(), f"scaling broke hyperbolicity at k={k}"
        assert sc.content == k * base.content, f"content wrong at k={k}"
    note("scaling by k multiplies S1 by k, divides R by k, preserves hyperbolicity")

    # -- Pillow rejects malformed input ----------------------------------
    for bad in [(1, 8, 8), (8, 2, 8), ()]:
        try:
            Pillow(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Pillow accepted malformed orders {bad}")
    note("Pillow rejects orders < 2, unsorted orders, and the empty multiset")

    print(f"\norbifold_enum selftest: {checks} checks, all passed")


if __name__ == "__main__":
    selftest()
