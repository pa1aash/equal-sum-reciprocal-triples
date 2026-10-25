"""F9 (Sections 5.4 and 8.4): the curve C_{27/2} and the growth of the degeneracies.

(a) The real locus of C_{27/2}: (x+y+z)(xy+yz+zx) = (27/2) xyz in the affine chart x+y+z = 1,
    drawn in barycentric coordinates (the positive triangle x, y, z > 0 is outlined; equal aspect).
    Its only positive rational points are the permutations of (1,4,4)/9 and (1,1,4)/6, i.e.
    (2,8,8) and (3,3,12) up to scale (rank C_{27/2} = 0, theory/diophantine/ranks.py), drawn in
    the pillow hues.
(b) calN(S)/S, with calN(S) the number of unordered pairs of hyperbolic triads with equal sum S' <= S
    and equal R (the cumulative count, calligraphic N in the paper, where N(S) counts the pairs
    at sum exactly S; theory/diophantine/data/per_S.csv, cum_pairs), ink step curve; the empirical
    law calN(S) = calN(599) + c sum_{s=600}^{S} (log s)^q fitted on 600 <= S <= 4800 with q = 4.5 (dashed),
    drawn heavier, on top, over its window (the fits with q = 4.0 and 5.0 are asserted to stay
    within 15% of the data there too, so a band would have no visible width); the exact count of the
    isosceles family {(2u+v)(u,v,v), (u+2v)(v,u,u)}/g and all its hyperbolic multiples (the base
    pair D_{1,4} is (2,8,8) ~ (3,3,12)), a rigorous lower bound at every S, whose growth is
    (c_iso + o(1)) S log S (thin grey step curve; SPEC rule 3).
    The split disc at S = 18 is the first pair, (2,8,8) ~ (3,3,12).

Asserted: the six points lie exactly on the curve; every row of groups.csv with S R = 27/2 is a
multiple of {(1,4,4), (1,1,4)}; per_S.csv cum_pairs is the cumulative sum of its pairs column;
the free cumulative fit (exponent_fits.fit_cumulative, window 600..4800) reproduces k = 4.56 of
data/exponent_fits.txt; the k >= 4 isosceles count equals families.txt at X = 600, ..., 4800; all
hyperbolic isosceles pairs with S <= 4800 number 2917 (families.py) and each is a pair of some group
in groups.csv; the first is (2,8,8) ~ (3,3,12).
"""
import re
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd, log

import numpy as np

from figlib import ROOT, check_only, finish, fs, rows
from common_import import import_from

LAMBDA = F(27, 2)


def curve_points():
    pts = {"2,8,8": sorted(set(permutations((F(1, 9), F(4, 9), F(4, 9))))),
           "3,3,12": sorted(set(permutations((F(1, 6), F(1, 6), F(4, 6)))))}
    for v in pts.values():
        assert len(v) == 3
        for x, y, z in v:
            assert x + y + z == 1 and (x + y + z) * (x * y + y * z + z * x) == LAMBDA * x * y * z
    return pts


def hyperbolic(t):
    p, q, r = t
    return q * r + p * r + p * q < p * q * r


def isosceles_pairs(X, kmin=1):
    """(S, triad1, triad2): the primitive pairs D_{u,v} = {(2u+v)(u,v,v), (u+2v)(v,u,u)}/g of the
    isosceles family (coprime u < v, g = gcd(2u+v, u+2v)) and their multiples kD, k >= kmin,
    with S <= X and both triads hyperbolic. D_{1,4} = {(2,8,8), (3,3,12)}."""
    out = []
    v = 2
    while (2 + v) * (1 + 2 * v) // 3 <= X:
        for u in range(1, v):
            if gcd(u, v) != 1:
                continue
            a, b = 2 * u + v, u + 2 * v
            g = gcd(a, b)
            SD = a * b // g
            t1 = tuple(sorted(x * a // g for x in (u, v, v)))
            t2 = tuple(sorted(x * b // g for x in (v, u, u)))
            k = kmin
            while k * SD <= X:
                k1, k2 = tuple(k * x for x in t1), tuple(k * x for x in t2)
                if hyperbolic(k1) and hyperbolic(k2):
                    out.append((k * SD, k1, k2))
                k += 1
        v += 1
    return out


def data():
    pts = curve_points()
    groups = rows("theory/diophantine/data/groups.csv")
    tri = lambda s: [tuple(int(x) for x in t.strip().strip("()").split(",")) for t in s.split(";")]
    on_curve = [g for g in groups if int(g["S"]) * F(int(g["R_num"]), int(g["R_den"])) == LAMBDA]
    assert on_curve and on_curve[0]["S"] == "18"
    for g in on_curve:
        k = int(g["S"]) // 18
        assert int(g["S"]) == 18 * k and sorted(tri(g["triples"])) == [(2 * k, 8 * k, 8 * k), (3 * k, 3 * k, 12 * k)]

    per = rows("theory/diophantine/data/per_S.csv")
    S = np.array([int(r["S"]) for r in per])
    pairs = np.array([int(r["pairs"]) for r in per])
    cum = np.array([int(r["cum_pairs"]) for r in per])
    assert np.array_equal(np.cumsum(pairs), cum) and S[-1] == 4800 and cum[-1] == 102719

    ef, _ = import_from("theory/diophantine", "exponent_fits")
    a, b = 600, 4800
    sel = (S >= a) & (S <= b)
    m = (sel.sum() // ef.BLOCK) * ef.BLOCK
    k_free = ef.fit_cumulative(S[sel][:m].astype(float), pairs[sel][:m].astype(float))
    txt = (ROOT / "theory/diophantine/data/exponent_fits.txt").read_text()
    k_txt = float(re.search(r"^all\s+600-4800\s+\S+ \+- \S+\s+(\S+) \+-", txt, re.M).group(1))
    assert abs(k_free - k_txt) < 0.005, (k_free, k_txt)
    base = cum[S == a - 1][0]
    fits = {}
    for q in (4.0, 4.5, 5.0):                 # c by least squares on log N, as fit_cumulative does
        w = np.cumsum(np.log(S[sel].astype(float)) ** q)
        y = cum[sel] - base
        c = np.exp(np.mean(np.log(y) - np.log(w)))
        fits[q] = (S[sel], base + c * w)
        assert np.max(np.abs(np.log(base + c * w) - np.log(cum[sel]))) < 0.15, q

    iso4 = isosceles_pairs(4800, kmin=4)                 # the k >= 4 count of families.txt
    iso = isosceles_pairs(4800)                          # every hyperbolic multiple: the curve drawn
    assert len(iso) == 2917 and min(iso)[:3] == (18, (2, 8, 8), (3, 3, 12))
    assert set(iso4) <= set(iso)
    known = {(int(g["S"]), t) for g in groups for t in combinations(sorted(tri(g["triples"])), 2)}
    for Sx, t1, t2 in iso:
        assert (Sx, tuple(sorted((t1, t2)))) in known, (Sx, t1, t2)
    fam = (ROOT / "theory/diophantine/data/families.txt").read_text()
    table = re.findall(r"^\s+(\d+)\s+\d+\s+[\d.]+\s+(\d+)\s+[\d.]+\s+\d+\s+\d+$", fam, re.M)
    assert len(table) == 6, table
    for X, n in table:
        assert sum(1 for s, *_ in iso4 if s <= int(X)) == int(n), (X, n)
    iso_S = np.sort([s for s, *_ in iso])
    iso_cum = np.searchsorted(iso_S, S, side="right")
    assert np.all(iso_cum <= cum)
    c_iso = 3 * log(2) / (2 * np.pi ** 2)
    assert abs(c_iso - 0.10535) < 1e-5
    return pts, (S, cum), fits, iso_cum, k_txt


def bary(x, y, z):
    """Barycentric (x+y+z = 1) to the plane: vertices (0,0), (1,0), (1/2, sqrt3/2)."""
    return y + z / 2, z * np.sqrt(3) / 2


def draw(pts, Ncum, fits, iso_cum, _k):
    fs.use()
    fig = fs.figure(62)
    # (a) the curve
    ax = fig.add_axes([0.03, 0.05, 0.40, 0.86 + 0.0])
    X, Y = np.meshgrid(np.linspace(-0.55, 1.55, 1400), np.linspace(-0.55, 1.42, 1300))
    z = 2 * Y / np.sqrt(3)
    y = X - z / 2
    x = 1 - y - z
    f = x * y + y * z + z * x - float(LAMBDA) * x * y * z
    ax.contour(X, Y, f, levels=[0], colors=[fs.GREY["ink"]], linewidths=fs.LW["regular"])
    tri = np.array([bary(1, 0, 0), bary(0, 1, 0), bary(0, 0, 1), bary(1, 0, 0)])
    ax.plot(tri[:, 0], tri[:, 1], color=fs.GREY["light"], lw=fs.LW["thin"], zorder=0)
    for sig, marker in (("2,8,8", "o"), ("3,3,12", "o")):
        for p in pts[sig]:
            u, v = bary(*(float(c) for c in p))
            ax.plot([u], [v], ls="none", marker=marker, ms=fs.MARKER_PT["large"], mfc=fs.PILLOW[sig],
                    mec=fs.GREY["ink"], mew=fs.LW["hair"], zorder=5)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(-0.55, 1.42)
    ax.set_aspect("equal")
    ax.set_axis_off()
    fs.letter_at(fig, 0.01, 0.985, "a")
    # (b) growth
    bx = fig.add_axes([0.58, 0.18, 0.40, 0.73])
    S, cum = Ncum
    m = S >= 18
    bx.step(S[m], cum[m] / S[m], where="post", color=fs.GREY["ink"], lw=fs.LW["regular"])
    bx.plot(fits[4.5][0], fits[4.5][1] / fits[4.5][0], color=fs.GREY["mid"], lw=fs.LW["heavy"],
            ls=(0, (3, 2)), zorder=4)
    mi = iso_cum > 0
    bx.step(S[mi], iso_cum[mi] / S[mi], where="post", color=fs.GREY["mid"], lw=fs.LW["thin"])
    bx.plot([18], [1 / 18], ls="none", marker="o", ms=fs.MARKER_PT["large"], fillstyle="left", mfc=fs.PILLOW["2,8,8"],
            mfcalt=fs.PILLOW["3,3,12"], mec=fs.GREY["ink"], mew=fs.LW["hair"], zorder=5)
    bx.set_xscale("log")
    bx.set_yscale("log")
    bx.set_xlim(15, 5000)
    bx.set_ylim(8e-3, 40)
    fs.decimal_log_ticks(bx.xaxis, [100, 1000])
    fs.decimal_log_ticks(bx.yaxis, [0.01, 0.1, 1, 10])
    bx.set_xlabel(r"$S$")
    bx.set_ylabel(r"$\mathcal{N}(S)/S$")
    fs.letter_at(fig, 0.48, 0.985, "b")
    return fig


if __name__ == "__main__":
    d = data()
    if not check_only():
        finish(draw(*d), "F9")
    print("F9: assertions passed")
