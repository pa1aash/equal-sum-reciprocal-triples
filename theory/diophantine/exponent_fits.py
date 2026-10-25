#!/usr/bin/env python3
"""
The log-power kappa, fitted the same way on the per-S and the cumulative series.

Model. n(s) = c (log s)^k for the number of degenerate pairs at sum s.
  * per-S:      Poisson maximum likelihood for n(s), s in a window [a, b];
  * cumulative: least squares on log of N(S) - N(a-1) against the exact sum
                sum_{s=a}^{S} c (log s)^k, S in [a, b].
Both fit the same two parameters on the same window, so the two estimates of
k are comparable. The per-S fit is the exact Poisson maximum-likelihood estimate
(Newton's method on the log-linear model). Uncertainties: residual block
bootstrap, applied identically to both fits: the per-S model is fitted, the
multiplicative residuals y/mu are resampled in blocks of 50 consecutive sums and
multiplied back onto the fitted trend, and both fits are repeated (400
resamples). (Resampling blocks of the counts themselves onto fixed abscissae
would destroy the trend and measure the spread of a trend-free null instead.)
For comparison the pure power law n(s) = c s^beta is fitted by Poisson maximum
likelihood on the same window, with the log-likelihood difference and the
Pearson dispersion of the log-power fit.

For contrast it also refits the earlier closed form N(S) = c S (log S)^K from
S = 100, and shows that K exceeds k by the secondary term of
sum_{s<=S} (log s)^k = S (log S)^k (1 - k/log S + ...).

Series fitted: all pairs; primitive pairs; dual (reciprocal) pairs with
scalings; primitive dual pairs; non-dual pairs (all, with scalings); primitive
non-dual pairs. Point estimates use the whole window; the bootstrap uses the
largest whole number of blocks in it. Manin's conjecture on the
dual quartic del Pezzo predicts k = 4 for primitive dual pairs and k = 5 for
dual pairs with scalings (per-S densities of the cumulative X (log X)^4 and
X (log X)^5).

holdouts() below: A(X) and N(X) to 6000 from the third enumeration; the family c s^b0 (log s)^q with b0
fixed (asymptotic exponent 1 + b0 for N), its deviance and the errors of both forward holdouts; the free law
c s^beta (log s)^q; and q, beta on sliding windows of width 2400. All with asserts.

Writes data/exponent_fits.txt.
"""

from __future__ import annotations

import sys
from itertools import combinations
from math import gcd
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cubic_group import dual, normalize  # noqa: E402
from enumerate_fast import load_groups, load_per_s  # noqa: E402

S_MAX = 4800
BLOCK, NBOOT = 50, 400
OUT: list[str] = []


def say(s=""):
    print(s, flush=True)
    OUT.append(s)


def series() -> dict[str, np.ndarray]:
    n = {k: np.zeros(S_MAX + 1) for k in ("all", "primitive", "dual", "primitive dual", "non-dual",
                                          "primitive non-dual")}
    for S, _, _, _, _, T in load_groups():
        for t1, t2 in combinations(T, 2):
            prim = gcd(gcd(*t1), gcd(*t2)) == 1
            d = dual(t1) == tuple(sorted(normalize(t2)))
            n["all"][S] += 1
            n["primitive"][S] += prim
            n["dual"][S] += d
            n["primitive dual"][S] += prim and d
            n["non-dual"][S] += not d
            n["primitive non-dual"][S] += prim and not d
    check = np.array([r["pairs"] for r in load_per_s()])
    assert np.array_equal(n["all"][10:], check), "group file and per-S file disagree"
    return n


def poisson_fit(x, y):
    """Exact Poisson MLE of log mu = a + b x by Newton's method; returns (a, b, mu)."""
    X = np.column_stack([np.ones_like(x), x])
    th = np.array([np.log(max(y.mean(), 1e-9)) - 4.0 * x.mean(), 4.0])
    for _ in range(100):
        mu = np.exp(X @ th)
        g = X.T @ (y - mu)
        H = X.T @ (X * mu[:, None])
        step = np.linalg.solve(H, g)
        th = th + step
        if np.max(np.abs(step)) < 1e-13:
            break
    mu = np.exp(X @ th)
    return th[0], th[1], mu


def fit_per_s(s, y):
    return poisson_fit(np.log(np.log(s)), y)[1]


def nll(y, mu):
    return float(np.sum(mu - y * np.log(mu)))


def fit_cumulative(s, y):
    cum = np.cumsum(y)
    ok = cum > 0
    ll = np.log(np.log(s))

    def resid(th):
        model = np.cumsum(np.exp(th[0] + th[1] * ll))
        return np.log(model[ok]) - np.log(cum[ok])
    return least_squares(resid, [0.0, 4.0]).x[1]


def bootstrap(fn, s, y, rng):
    """Residual block bootstrap around the fitted per-S trend (the trend is kept)."""
    nb = len(s) // BLOCK
    mu = poisson_fit(np.log(np.log(s)), y)[2]
    r = y / mu
    est = []
    for _ in range(NBOOT):
        idx = rng.integers(0, nb, nb)
        rb = np.concatenate([r[i * BLOCK:(i + 1) * BLOCK] for i in idx])
        est.append(fn(s, mu * rb))
    return np.std(est)


def closed_form_K(S_all, N_all, lo=100):
    m = S_all >= lo
    S, N = S_all[m], N_all[m]
    r = least_squares(lambda th: th[0] + np.log(S) + th[1] * np.log(np.log(S)) - np.log(N), [-8, 5])
    return r.x[1]


def main() -> int:
    rng = np.random.default_rng(20261001)
    n = series()
    s_all = np.arange(S_MAX + 1, dtype=float)
    windows = [(600, 4800), (1200, 4800), (600, 2400), (2400, 4800)]
    say("k in n(s) = c (log s)^k; same model, same window, both series; +- block-bootstrap sd")
    say(f"{'series':16s} {'window':>11s} {'per-S k':>14s} {'cumulative k':>14s}")
    table = {}
    for name in n:
        for a, b in windows:
            s = s_all[a:b + 1]
            y = n[name][a:b + 1]
            kp = fit_per_s(s, y)          # point estimates on the whole window [a, b]
            kc = fit_cumulative(s, y)
            # the bootstrap uses a whole number of blocks, so that resamples have equal length
            m = (len(s) // BLOCK) * BLOCK
            sp = bootstrap(fit_per_s, s[:m], y[:m], rng)
            sc = bootstrap(fit_cumulative, s[:m], y[:m], rng)
            table[(name, a, b)] = (kp, sp, kc, sc)
            say(f"{name:16s} {a:5d}-{b:<5d} {kp:7.2f} +- {sp:4.2f} {kc:7.2f} +- {sc:4.2f}")
        say()

    say("pure power law n(s) = c s^beta, Poisson MLE, window 600-4800 (for comparison with the log-power model)")
    for name in ("all", "primitive"):
        s = s_all[600:4801]
        y = n[name][600:4801]
        _, beta, mu_p = poisson_fit(np.log(s), y)
        _, k, mu_l = poisson_fit(np.log(np.log(s)), y)
        disp = float(np.sum((y - mu_l) ** 2 / mu_l) / (len(y) - 2))
        d = nll(y, mu_p) - nll(y, mu_l)
        say(f"{name:16s} beta = {beta:.3f} (cumulative exponent about {1 + beta:.2f}); log-power k = {k:.3f}; "
            f"NLL(power) - NLL(log-power) = {d:.2f}; Pearson dispersion {disp:.2f}; scaled difference {d / disp:.2f}")
    say()

    # the sums 4801..6000 as a holdout. Both laws are fitted per sum on 600..4800 (above); their
    # predictions of N(6000) - N(4800) are compared with the exact count, 145679 - 102719 (Table 2 of the note;
    # N(6000) from the third enumeration, N(4800) from per_S.csv).
    N6000, N4800 = 145679, int(np.cumsum(n["all"])[4800])
    assert N4800 == 102719, N4800
    s = s_all[600:4801]
    y = n["all"][600:4801]
    a_p, beta, mu_p = poisson_fit(np.log(s), y)
    a_l, k, mu_l = poisson_fit(np.log(np.log(s)), y)
    disp = float(np.sum((y - mu_l) ** 2 / mu_l) / (len(y) - 2))
    sh = np.arange(4801, 6001, dtype=float)
    pred_l = float(np.sum(np.exp(a_l + k * np.log(np.log(sh)))))
    pred_p = float(np.sum(np.exp(a_p + beta * np.log(sh))))
    act = N6000 - N4800
    say(f"holdout 4801-6000 (all pairs): exact N(6000) - N(4800) = {act}; log-power law predicts {pred_l:.0f} "
        f"({100 * (pred_l / act - 1):+.1f}%), power law {pred_p:.0f} ({100 * (pred_p / act - 1):+.1f}%); "
        f"misses in dispersion-scaled Poisson sd: {abs(pred_l - act) / np.sqrt(disp * act):.1f} and "
        f"{abs(pred_p - act) / np.sqrt(disp * act):.1f}")
    say()

    # how much the per-sum exponent depends on the arithmetic of s (Lemma 5: n(s) is a divisor sum), and
    # where the two laws are tangent: the log-slope k/log s of c (log s)^k equals beta at s = exp(k/beta).
    s = s_all[600:4801]
    y = n["all"][600:4801]
    si = s.astype(int)
    classes = {"gcd(s,6)=1": np.gcd(si, 6) == 1, "s odd": si % 2 == 1, "s even": si % 2 == 0, "6 | s": si % 6 == 0}
    kres = {}
    for lab, mask in classes.items():
        kres[lab] = fit_per_s(s[mask], y[mask])
    say("per-sum k for all pairs on 600-4800, by residue class of s: "
        + "; ".join(f"{lab}: {v:.2f}" for lab, v in kres.items()))
    assert 4.0 < kres["gcd(s,6)=1"] < 4.1 and 4.25 < kres["s odd"] < 4.35 and 4.55 < kres["s even"] < 4.65 \
        and 4.9 < kres["6 | s"] < 5.0, kres
    s_tan = float(np.exp(k / beta))
    say(f"tangency of the two fitted laws: k/log s = beta at s = {s_tan:.0f}")
    assert 2000 < s_tan < 2200, s_tan
    assert abs(100 * (pred_l / act - 1) - 0.25) < 0.01 and abs(100 * (pred_p / act - 1) - 2.33) < 0.01
    say()

    N_all = np.cumsum(n["all"])
    K = closed_form_K(s_all[10:], N_all[10:])
    kp = table[("all", 600, 4800)][0]
    L = np.log(3000)
    say(f"closed form c S (log S)^K fitted on 100..4800 (the earlier cumulative fit): K = {K:.2f}")
    say(f"secondary-term prediction from per-S k = {kp:.2f}: K ~ k (1 + 1/log S) = "
        f"{kp * (1 + 1 / L):.2f} at S ~ 3000")
    # (a) a second, disjoint holdout: both laws fitted per sum on 600..3600 predict the pairs with
    # 3601 <= s <= 4800; (b) the number A(X) of primitive pairs with sum at most X, the quantity of Question E,
    # and its local exponent log2(A(X)/A(X/2)).
    s = s_all[600:3601]
    y = n["all"][600:3601]
    a_p2, beta2, _ = poisson_fit(np.log(s), y)
    a_l2, k2, _ = poisson_fit(np.log(np.log(s)), y)
    sh = np.arange(3601, 4801, dtype=float)
    act2 = float(np.sum(n["all"][3601:4801]))
    pl2 = float(np.sum(np.exp(a_l2 + k2 * np.log(np.log(sh)))))
    pp2 = float(np.sum(np.exp(a_p2 + beta2 * np.log(sh))))
    el2, ep2 = 100 * (pl2 / act2 - 1), 100 * (pp2 / act2 - 1)
    say(f"holdout 3601-4800 (fits on 600-3600): exact {act2:.0f}; log-power law {el2:+.2f}%, power law {ep2:+.2f}%")
    assert abs(el2 - 0.21) < 0.02 and abs(ep2 - 2.32) < 0.02, (el2, ep2)
    A = np.cumsum(n["primitive"])
    vals = {X: int(A[X]) for X in (500, 1000, 2000, 2400, 4000, 4800)}
    loc = {X: float(np.log2(A[X] / A[X // 2])) for X in (1000, 2000, 4000, 4800)}
    say("primitive pairs A(X): " + ", ".join(f"A({X}) = {v}" for X, v in vals.items()) + "; local exponents "
        + ", ".join(f"{X}: {v:.3f}" for X, v in loc.items()))
    assert vals[1000] == 4141 and vals[2000] == 12512 and vals[4000] == 35567 and vals[4800] == 46254, vals
    say()
    holdouts(n, s_all)
    (HERE / "data" / "exponent_fits.txt").write_text("\n".join(OUT) + "\n")
    return 0


# what the holdouts can and cannot show ---------------------------------------------------
SUMMARY6000 = HERE.parents[1] / "checks" / "enumeration6000" / "data" / "summary_10_6000.txt"


def per_sum_6000(n):
    """Per-sum counts of all pairs and of primitive pairs for 10 <= s <= 6000, from the third enumeration
    (checks/enumeration6000/enum.c; columns S ntriads nfibres npairs maxfibre adj nonadj primfibres primpairs),
    asserted against data/per_S.csv for s <= 4800 and against N(6000) = 145679."""
    al, pr = np.zeros(6001), np.zeros(6001)
    for line in SUMMARY6000.read_text().split("\n"):
        if line.strip():
            v = [int(t) for t in line.split()]
            al[v[0]], pr[v[0]] = v[3], v[8]
    assert np.array_equal(al[10:4801], n["all"][10:4801]) and np.array_equal(pr[10:4801], n["primitive"][10:4801])
    assert int(al.sum()) == 145679, al.sum()
    return al, pr


def poisson_fit_offset(x, y, off):
    """Poisson MLE of log mu = a + b x + off (off fixed), by Newton's method; returns (a, b, mu)."""
    X = np.column_stack([np.ones_like(x), x])
    th = np.array([np.log(max(y.mean(), 1e-9)) - float(np.mean(off)) - 3.0 * x.mean(), 3.0])
    for _ in range(200):
        mu = np.exp(X @ th + off)
        step = np.linalg.solve(X.T @ (X * mu[:, None]), X.T @ (y - mu))
        th = th + step
        if np.max(np.abs(step)) < 1e-12:
            break
    return th[0], th[1], np.exp(X @ th + off)


def deviance(y, mu):
    t = np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0)
    return float(2 * np.sum(t - (y - mu)))


def holdouts(n, s_all):
    al, pr = per_sum_6000(n)
    A = np.cumsum(pr)
    a_vals = {X: int(A[X]) for X in (500, 1000, 2000, 2400, 2700, 3000, 4000, 4800, 5400, 6000)}
    say("primitive pairs to 6000 (third enumeration): "
        + ", ".join(f"A({X}) = {v}" for X, v in a_vals.items())
        + f"; local exponents 5400: {np.log2(A[5400] / A[2700]):.3f}, 6000: {np.log2(A[6000] / A[3000]):.3f}")
    assert a_vals[500] == 1328 and a_vals[2400] == 16565 and a_vals[6000] == 63876, a_vals
    N = np.cumsum(al)
    say(f"N(2700) = {int(N[2700])}, N(5400) = {int(N[5400])}, local exponent at 5400 {np.log2(N[5400] / N[2700]):.3f}")
    assert int(N[5400]) == 123559

    # (a) the family n(s) = c s^b0 (log s)^q, b0 fixed, (c, q) fitted per sum on 600..4800 (Poisson MLE); under this
    # law N(X) has asymptotic exponent 1 + b0. Deviance on the window and the errors of both forward holdouts.
    s, y = s_all[600:4801], al[600:4801]
    s2, y2 = s_all[600:3601], al[600:3601]
    h1, h2 = np.arange(4801, 6001, dtype=float), np.arange(3601, 4801, dtype=float)
    act1, act2 = float(al[4801:6001].sum()), float(al[3601:4801].sum())
    say("family c s^b0 (log s)^q, (c, q) fitted per sum (Poisson MLE); deviance on 600..4800; "
        "holdout errors 4801..6000 (fit 600..4800) and 3601..4800 (fit 600..3600)")
    fam = {}
    for b0 in (0.0, 0.1, 0.2, 0.3, 0.5):
        c0, q, mu = poisson_fit_offset(np.log(np.log(s)), y, b0 * np.log(s))
        dev = deviance(y, mu)
        p1 = float(np.sum(np.exp(c0 + q * np.log(np.log(h1)) + b0 * np.log(h1))))
        c2, q2, _ = poisson_fit_offset(np.log(np.log(s2)), y2, b0 * np.log(s2))
        p2 = float(np.sum(np.exp(c2 + q2 * np.log(np.log(h2)) + b0 * np.log(h2))))
        fam[b0] = (q, dev, 100 * (p1 / act1 - 1), 100 * (p2 / act2 - 1))
        say(f"  b0 = {b0:.1f}: q = {q:.2f}, deviance {dev:.1f} on {len(y) - 2} d.f., holdouts {fam[b0][2]:+.2f}% and "
            f"{fam[b0][3]:+.2f}%")
    _, beta, mu_p = poisson_fit(np.log(s), y)
    dev_p = deviance(y, mu_p)
    say(f"  power law (q = 0): beta = {beta:.3f}, deviance {dev_p:.1f}; deviance / d.f. of the b0 = 0 law: "
        f"{fam[0.0][1] / (len(y) - 2):.2f}")
    # b0 = 0 is the law of the note; laws with asymptotic exponent 1.1 to 1.3 fit almost as well
    assert abs(fam[0.0][0] - 4.49) < 0.01 and abs(fam[0.0][2] - 0.25) < 0.01 and abs(fam[0.0][3] - 0.21) < 0.01
    assert all(fam[b][1] - fam[0.0][1] < 10 for b in (0.1, 0.2, 0.3)), fam
    assert all(abs(fam[b][2]) < 1.5 and abs(fam[b][3]) < 1.5 for b in (0.1, 0.2, 0.3)), fam
    assert fam[0.0][2] < fam[0.1][2] < fam[0.2][2] < fam[0.3][2] < fam[0.5][2], fam
    assert 5.5 < fam[0.0][1] / (len(y) - 2) < 7.0
    assert dev_p - fam[0.0][1] > 10

    # (b) the free three-parameter law c s^beta (log s)^q on 600..4800 and on 600..6000
    for lo, hi in ((600, 4800), (600, 6000)):
        ss, yy = np.arange(lo, hi + 1, dtype=float), al[lo:hi + 1]
        X = np.column_stack([np.ones_like(ss), np.log(ss), np.log(np.log(ss))])
        th = np.array([np.log(yy.mean()) - 4.0 * np.log(np.log(ss)).mean(), 0.0, 4.0])
        for _ in range(200):
            mu = np.exp(X @ th)
            step = np.linalg.solve(X.T @ (X * mu[:, None]), X.T @ (yy - mu))
            th = th + step
            if np.max(np.abs(step)) < 1e-12:
                break
        say(f"  free law c s^beta (log s)^q on {lo}..{hi}: beta = {th[1]:.2f}, q = {th[2]:.2f}")
        assert th[1] < 0, th

    # (c) sliding windows of width 2400: per-sum q (all, primitive) and power-law beta (all)
    say("sliding windows of width 2400, per sum: q (all pairs), q (primitive pairs), beta (all pairs)")
    qs, qp, bs = [], [], []
    for lo in range(600, 3601, 600):
        hi = lo + 2400
        ss = np.arange(lo, hi + 1, dtype=float)
        qs.append(fit_per_s(ss, al[lo:hi + 1]))
        qp.append(fit_per_s(ss, pr[lo:hi + 1]))
        bs.append(poisson_fit(np.log(ss), al[lo:hi + 1])[1])
        say(f"  [{lo}, {hi}]: q = {qs[-1]:.2f}, q_prim = {qp[-1]:.2f}, beta = {bs[-1]:.3f}")
    assert qs[0] > qs[-1] + 0.25 and qp[0] > qp[-1] + 0.5 and bs[0] > bs[-1] + 0.1, (qs, qp, bs)
    assert 4.5 < qs[0] < 4.6 and 4.2 < qs[-1] < 4.3 and 3.7 < qp[0] < 3.8 and 3.0 < qp[-1] < 3.2, (qs, qp)
    assert 0.6 < bs[0] < 0.65 and 0.48 < bs[-1] < 0.52, bs


if __name__ == "__main__":
    sys.exit(main())
