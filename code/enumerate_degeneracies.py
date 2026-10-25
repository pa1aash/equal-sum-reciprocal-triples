#!/usr/bin/env python3
"""
Two-coefficient spectral degeneracies among hyperbolic triangular pillows.

A two-coefficient degeneracy is a set of distinct hyperbolic triples sharing
the paper's signature

    sigma = (S_1, R),    S_1 = p + q + r,    R = 1/p + 1/q + 1/r.

Both components matter. Keying on any invariant that omits S_1 does not search
for the paper's object at all.

This module sweeps 10 <= S <= 600 exactly -- integers and Fraction only -- and
reports both counting conventions side by side, because both are used:

    N_pairs(S)    = sum over sigma-groups of C(k, 2)
    N_classes(S)  = number of sigma-groups with k > 1

They agree up to S = 100 and diverge from S = 200 onward, once some sigma-value
is shared by three or more triples.

INTEGRITY
---------
There is no fallback data anywhere in this file. If a search that must succeed
returns nothing, that is a failure and raises; it is never papered over with a
canned result. Every checkpoint is an assert, and `main` exits nonzero if any
of them fails.

Floating point appears in exactly one function, `power_law_exponent`, which is
marked as such. Everything else is exact.
"""

from __future__ import annotations

import csv
import sys
from fractions import Fraction
from math import gcd
from pathlib import Path
from typing import Iterator, NamedTuple, Sequence

from orbifold_enum import Pillow, enumerate_by_sum, signature

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "data"

S_MIN = 10
S_MAX = 600

# Independently computed cumulative counts. These are the hard self-check: the
# sweep must reproduce every one of them or the run fails.
CHECKPOINTS: dict[int, tuple[int, int]] = {
    100: (92, 92),
    200: (386, 380),
    300: (840, 822),
    400: (1496, 1468),
    500: (2210, 2158),
    600: (3067, 2977),
}

# Fit window for the empirical growth exponent.
FIT_LO, FIT_HI = 50, 600


class Degeneracy(NamedTuple):
    """One sigma-group of two or more distinct triples."""

    S: int
    R: Fraction
    triples: tuple[tuple[int, ...], ...]
    content: int          # gcd of every cone order in every triple
    primitive: bool       # content == 1
    base_sum: int         # S // content -- the sum this is a scaling of

    @property
    def multiplicity(self) -> int:
        return len(self.triples)

    @property
    def pair_count(self) -> int:
        k = len(self.triples)
        return k * (k - 1) // 2


class SumRecord(NamedTuple):
    """Degeneracy counts at one cone-order sum."""

    S: int
    n_triads: int
    n_pairs: int
    n_classes: int
    n_primitive_classes: int
    cum_pairs: int
    cum_classes: int
    cum_primitive_classes: int


# --------------------------------------------------------------------------
# Exact enumeration
# --------------------------------------------------------------------------

def degeneracy_groups(S: int) -> list[Degeneracy]:
    """
    Every sigma-group of size >= 2 at cone-order sum S, exactly.

    Within a fixed S the signature reduces to R, so grouping is by reduced R.
    Returns [] when S genuinely admits no degeneracy -- that is a real result,
    not a failure, and callers that require a non-empty answer assert it
    themselves.
    """
    buckets: dict[tuple[int, int], list[tuple[int, ...]]] = {}
    for orders in enumerate_by_sum(3, S):
        _, num, den = signature(orders)
        buckets.setdefault((num, den), []).append(orders)

    groups: list[Degeneracy] = []
    for (num, den), triples in buckets.items():
        if len(triples) < 2:
            continue
        content = 0
        for t in triples:
            for m in t:
                content = gcd(content, m)
        groups.append(
            Degeneracy(
                S=S,
                R=Fraction(num, den),
                triples=tuple(sorted(triples)),
                content=content,
                primitive=(content == 1),
                base_sum=S // content,
            )
        )
    groups.sort(key=lambda g: g.R)
    return groups


def count_triads(S: int) -> int:
    """Number of hyperbolic triads of cone-order sum S."""
    return sum(1 for _ in enumerate_by_sum(3, S))


def sweep(s_min: int = S_MIN, s_max: int = S_MAX) -> tuple[list[SumRecord], dict[int, list[Degeneracy]]]:
    """
    Enumerate degeneracies for every sum in [s_min, s_max].

    Returns the per-sum records (with running cumulatives) and the full group
    detail keyed by sum.
    """
    records: list[SumRecord] = []
    detail: dict[int, list[Degeneracy]] = {}
    cum_p = cum_c = cum_prim = 0

    for S in range(s_min, s_max + 1):
        groups = degeneracy_groups(S)
        detail[S] = groups
        n_pairs = sum(g.pair_count for g in groups)
        n_classes = len(groups)
        n_prim = sum(1 for g in groups if g.primitive)
        cum_p += n_pairs
        cum_c += n_classes
        cum_prim += n_prim
        records.append(
            SumRecord(
                S=S,
                n_triads=count_triads(S),
                n_pairs=n_pairs,
                n_classes=n_classes,
                n_primitive_classes=n_prim,
                cum_pairs=cum_p,
                cum_classes=cum_c,
                cum_primitive_classes=cum_prim,
            )
        )
    return records, detail


# --------------------------------------------------------------------------
# The one place floating point is permitted
# --------------------------------------------------------------------------

def power_law_exponent(
    records: Sequence[SumRecord],
    field: str,
    lo: int = FIT_LO,
    hi: int = FIT_HI,
) -> float:
    """
    ###################################################################
    ## THE ONLY FLOATING-POINT COMPUTATION IN THIS HARNESS.          ##
    ## Ordinary least squares slope of log(cumulative count) against ##
    ## log(S) over [lo, hi]. This is an empirical fit reported as    ##
    ## evidence for a conjecture -- it is not, and must not be used  ##
    ## as, a verification of any exact claim.                        ##
    ###################################################################
    """
    from math import log

    xs, ys = [], []
    for rec in records:
        if lo <= rec.S <= hi:
            value = getattr(rec, field)
            if value > 0:
                xs.append(log(rec.S))
                ys.append(log(value))
    if len(xs) < 2:
        raise RuntimeError(
            f"power_law_exponent: fewer than two usable points for {field} "
            f"in [{lo}, {hi}]; refusing to report a slope"
        )
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        raise RuntimeError("power_law_exponent: degenerate abscissa")
    return num / den


# --------------------------------------------------------------------------
# Self-checks. Every one is an assert.
# --------------------------------------------------------------------------

def check_base_case(detail: dict[int, list[Degeneracy]]) -> None:
    """S = 18: exactly one degeneracy, exactly {O(2,8,8), O(3,3,12)}, R = 3/4."""
    for S in range(S_MIN, 18):
        groups = detail[S]
        assert not groups, (
            f"Theorem A violated: degeneracy found at S = {S} <= 17: {groups}"
        )

    groups = detail[18]
    assert len(groups) == 1, (
        f"S = 18 must carry exactly one degeneracy, found {len(groups)}: {groups}"
    )
    g = groups[0]
    assert g.R == Fraction(3, 4), f"S = 18 degeneracy at R = {g.R}, expected 3/4"
    assert g.triples == ((2, 8, 8), (3, 3, 12)), (
        f"S = 18 pair is {g.triples}, expected ((2,8,8), (3,3,12))"
    )
    assert g.primitive, "the S = 18 base degeneracy must be primitive"
    assert g.content == 1, f"S = 18 content is {g.content}, expected 1"
    # The third heat coefficient must separate them.
    p3 = [Pillow(t).P(3) for t in g.triples]
    assert p3 == [1032, 1782], f"P_3 values {p3}, expected [1032, 1782]"
    assert p3[0] != p3[1], "P_3 fails to separate the base degeneracy"


def check_s36(detail: dict[int, list[Degeneracy]]) -> None:
    """S = 36: the primitive pair {O(6,15,15), O(8,8,20)} at R = 3/10 exists."""
    groups = detail[36]
    if not groups:
        raise RuntimeError(
            "S = 36 returned no degeneracies; the sweep is broken. "
            "A search that must succeed returning nothing is a failure, "
            "never a reason to substitute canned data."
        )
    by_R = {g.R: g for g in groups}

    target = Fraction(3, 10)
    assert target in by_R, (
        f"the primitive S = 36 degeneracy at R = 3/10 was not found; "
        f"R-values present: {sorted(by_R)}"
    )
    prim = by_R[target]
    assert prim.triples == ((6, 15, 15), (8, 8, 20)), (
        f"S = 36 primitive pair is {prim.triples}, expected ((6,15,15), (8,8,20))"
    )
    assert prim.primitive, (
        f"the S = 36 pair at R = 3/10 must be primitive, but content = {prim.content}"
    )
    assert prim.content == 1, f"content is {prim.content}, expected 1"

    # And the scaled copy of the base pair must also be there, flagged as such.
    scaled_R = Fraction(3, 8)
    assert scaled_R in by_R, f"the scaled copy at R = 3/8 is missing from S = 36"
    sc = by_R[scaled_R]
    assert sc.triples == ((4, 16, 16), (6, 6, 24)), (
        f"S = 36 scaled pair is {sc.triples}, expected ((4,16,16), (6,6,24))"
    )
    assert not sc.primitive, "the S = 36 pair at R = 3/8 must be flagged as scaled"
    assert sc.content == 2, f"scaled content is {sc.content}, expected 2"
    assert sc.base_sum == 18, f"scaled base_sum is {sc.base_sum}, expected 18"


def check_primitivity_consistency(detail: dict[int, list[Degeneracy]]) -> None:
    """
    Every non-primitive degeneracy must reduce, on dividing out its content, to
    a genuine degeneracy at the smaller sum. This makes the primitive/scaled
    flag falsifiable rather than decorative.
    """
    for S, groups in detail.items():
        for g in groups:
            if g.primitive:
                continue
            k = g.content
            reduced = tuple(sorted(tuple(m // k for m in t) for t in g.triples))
            base = g.base_sum
            assert base in detail, (
                f"degeneracy at S = {S} claims base sum {base}, outside the swept range"
            )
            matches = [b for b in detail[base] if b.triples == reduced]
            assert matches, (
                f"degeneracy at S = {S} with content {k} does not reduce to a "
                f"degeneracy at S = {base}: reduced = {reduced}"
            )
            assert matches[0].R == g.R * k, (
                f"scaling law broken: R at S={S} is {g.R}, base R is {matches[0].R}"
            )


def check_lower_bound(records: Sequence[SumRecord]) -> None:
    """Theorem 5.2: cumulative count >= floor(S / 18) for every S >= 18."""
    for rec in records:
        if rec.S >= 18:
            assert rec.cum_pairs >= rec.S // 18, (
                f"lower bound violated at S = {rec.S}: "
                f"cum_pairs = {rec.cum_pairs} < {rec.S // 18}"
            )
            assert rec.cum_classes >= rec.S // 18, (
                f"lower bound violated at S = {rec.S}: "
                f"cum_classes = {rec.cum_classes} < {rec.S // 18}"
            )


def check_checkpoints(records: Sequence[SumRecord]) -> list[str]:
    """Compare the sweep against the independently computed cumulative counts."""
    by_S = {r.S: r for r in records}
    failures: list[str] = []
    for S, (want_pairs, want_classes) in sorted(CHECKPOINTS.items()):
        if S not in by_S:
            failures.append(f"S = {S}: not present in sweep")
            continue
        rec = by_S[S]
        ok_p = rec.cum_pairs == want_pairs
        ok_c = rec.cum_classes == want_classes
        status = "PASS" if (ok_p and ok_c) else "FAIL"
        print(
            f"  {status}  S = {S:3d}   pairs {rec.cum_pairs:5d} (want {want_pairs:5d})"
            f"   classes {rec.cum_classes:5d} (want {want_classes:5d})"
        )
        if not ok_p:
            failures.append(
                f"S = {S}: cum_pairs = {rec.cum_pairs}, expected {want_pairs}"
            )
        if not ok_c:
            failures.append(
                f"S = {S}: cum_classes = {rec.cum_classes}, expected {want_classes}"
            )
    return failures


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def write_csvs(records: Sequence[SumRecord], detail: dict[int, list[Degeneracy]]) -> tuple[Path, Path]:
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    per_sum = DATA_DIR / "degeneracies.csv"
    with per_sum.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "S", "n_triads",
            "N_pairs", "N_classes", "N_primitive_classes",
            "cum_pairs", "cum_classes", "cum_primitive_classes",
        ])
        for r in records:
            w.writerow([
                r.S, r.n_triads,
                r.n_pairs, r.n_classes, r.n_primitive_classes,
                r.cum_pairs, r.cum_classes, r.cum_primitive_classes,
            ])

    groups_csv = DATA_DIR / "degeneracy-groups.csv"
    with groups_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "S", "R_num", "R_den", "multiplicity", "pair_count",
            "primitive", "content", "base_sum", "triples",
        ])
        for S in sorted(detail):
            for g in detail[S]:
                w.writerow([
                    g.S, g.R.numerator, g.R.denominator,
                    g.multiplicity, g.pair_count,
                    "primitive" if g.primitive else "scaled",
                    g.content, g.base_sum,
                    " ; ".join("(" + ",".join(str(m) for m in t) + ")" for t in g.triples),
                ])

    return per_sum, groups_csv


# --------------------------------------------------------------------------

def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    s_max = S_MAX
    if argv:
        try:
            s_max = int(argv[0])
        except ValueError:
            print(f"usage: {Path(__file__).name} [S_MAX]", file=sys.stderr)
            return 2

    print(f"enumerating two-coefficient degeneracies, {S_MIN} <= S <= {s_max}")
    print("signature: sigma = (S_1, R), both components. exact rational arithmetic.\n")

    records, detail = sweep(S_MIN, s_max)

    failures: list[str] = []

    print("structural self-checks")
    for label, fn in [
        ("Theorem A / B base case at S = 18", lambda: check_base_case(detail)),
        ("primitive and scaled pairs at S = 36", lambda: check_s36(detail)),
        ("primitivity flags reduce correctly", lambda: check_primitivity_consistency(detail)),
        ("Theorem 5.2 lower bound floor(S/18)", lambda: check_lower_bound(records)),
    ]:
        try:
            fn()
        except (AssertionError, RuntimeError) as exc:
            print(f"  FAIL  {label}\n        {exc}")
            failures.append(f"{label}: {exc}")
        else:
            print(f"  PASS  {label}")

    print("\ncumulative-count checkpoints")
    if s_max >= max(CHECKPOINTS):
        failures.extend(check_checkpoints(records))
    else:
        partial = {S: v for S, v in CHECKPOINTS.items() if S <= s_max}
        if partial:
            failures.extend(check_checkpoints(records))
        print(f"  NOTE  sweep stopped at S = {s_max}; "
              f"{len([S for S in CHECKPOINTS if S > s_max])} checkpoint(s) not exercised")

    per_sum, groups_csv = write_csvs(records, detail)
    print(f"\nwrote {per_sum.relative_to(REPO)}")
    print(f"wrote {groups_csv.relative_to(REPO)}")

    last = records[-1]
    n_groups = sum(len(g) for g in detail.values())
    n_prim = sum(1 for gs in detail.values() for g in gs if g.primitive)
    print(
        f"\nat S = {last.S}: cum_pairs = {last.cum_pairs}, "
        f"cum_classes = {last.cum_classes}, "
        f"primitive classes = {last.cum_primitive_classes} "
        f"({n_prim}/{n_groups} of all groups)"
    )

    try:
        exp_pairs = power_law_exponent(records, "cum_pairs")
        exp_classes = power_law_exponent(records, "cum_classes")
        print(
            f"OLS exponent over {FIT_LO} <= S <= {min(FIT_HI, s_max)} "
            f"(floating point, empirical): "
            f"pairs {exp_pairs:.3f}, classes {exp_classes:.3f}"
        )
    except RuntimeError as exc:
        print(f"  NOTE  exponent not computed: {exc}")

    if failures:
        print(f"\nFAIL -- {len(failures)} check(s) failed:")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("\nPASS -- all degeneracy checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
