#!/usr/bin/env python3
"""
The growth law of the cumulative degeneracy count, tested out of sample, and
an empirical check of the birthday heuristic of the paper's section 5.3.

Reads data/per_S.csv (written by enumerate_fast.py). Floating point is used
for fitting and for the heuristic's predicted counts only; every count it is
compared against is exact.

Part 1  Growth-law fits. Models for the cumulative pair count N(S):
          A  c S^2                      (Conjecture 5.3 as stated)
          B  c S^2 / log S
          C  c S^a
          D  c S^2 (1 + d / log S)
          E  c S (log S)^k
          F  c S^a (log S)^k
        fitted by least squares on log N over a training window and scored
        on a disjoint later window: fit 100..1000 -> predict 1001..2000 (the
        protocol asked for), and fit 100..2400 -> predict 2401..4800.
        Also the local exponent d log N / d log S over doubling windows, and
        the same analysis for cumulative primitive classes.

Part 2  Birthday heuristic check, on sampled sums S:
          T(S)   triads, D(S) distinct R-values, P(S) pairs (exact)
          M_eff  = T^2 / (2 P), the number of "possible R-values" that a
                   uniform birthday model would need to reproduce P(S)
          H(S)   number of reduced fractions in (0,1) with denominator
                   <= (S/3)^3 -- the honest count of candidate R-values
          P_null = sum_d C(n_d, 2) / (phi(d) * w): the pair count expected if,
                   conditional on its reduced denominator d, each R were
                   uniform among the phi(d) w admissible numerators (w the
                   width of the R-range). Comparing P_null with P tests the
                   independence assumption directly.

Usage: growth_fits.py            (writes data/growth_fits.txt, data/birthday.csv)
"""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from collections import Counter
from math import comb, log, pi
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from enumerate_fast import build_core, load_per_s  # noqa: E402

DATA = HERE / "data"

MODELS = {
    "A  c S^2":            (lambda th, S: th[0] + 2 * np.log(S), [-5.0]),
    "B  c S^2/log S":      (lambda th, S: th[0] + 2 * np.log(S) - np.log(np.log(S)), [-3.0]),
    "C  c S^a":            (lambda th, S: th[0] + th[1] * np.log(S), [-5.0, 2.0]),
    "D  c S^2(1+d/log S)": (lambda th, S: th[0] + 2 * np.log(S) + np.log(np.abs(1 + th[1] / np.log(S))), [-5.0, -1.0]),
    "E  c S (log S)^k":    (lambda th, S: th[0] + np.log(S) + th[1] * np.log(np.log(S)), [-8.0, 4.0]),
    "F  c S^a (log S)^k":  (lambda th, S: th[0] + th[1] * np.log(S) + th[2] * np.log(np.log(S)), [-8.0, 1.0, 4.0]),
}


def fit(S, y, model):
    f, th0 = model
    res = least_squares(lambda th: f(th, S) - y, th0)
    return res.x


def score(S, N, model, train, test):
    S = np.asarray(S, float)
    y = np.log(np.asarray(N, float))
    tr = (S >= train[0]) & (S <= train[1])
    te = (S >= test[0]) & (S <= test[1])
    th = fit(S[tr], y[tr], model)
    f = model[0]
    rel_in = np.exp(f(th, S[tr]) - y[tr]) - 1
    rel_out = np.exp(f(th, S[te]) - y[te]) - 1
    return th, np.sqrt(np.mean(rel_in ** 2)), np.sqrt(np.mean(rel_out ** 2)), rel_out[-1], np.max(np.abs(rel_out))


def local_exponents(S, N, windows):
    idx = {s: i for i, s in enumerate(S)}
    out = []
    for a, b in windows:
        out.append((a, b, log(N[idx[b]] / N[idx[a]]) / log(b / a)))
    return out


# --------------------------------------------------------------------------
# birthday check
# --------------------------------------------------------------------------

def spf_table(n: int) -> list[int]:
    spf = list(range(n + 1))
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def primes_of(m: int, spf: list[int]) -> set[int]:
    out = set()
    while m > 1:
        p = spf[m]
        out.add(p)
        while m % p == 0:
            m //= p
    return out


def birthday_row(exe: Path, S: int, spf: list[int]) -> dict:
    text = subprocess.run([str(exe), "dump", str(S)], check=True,
                          capture_output=True, text=True).stdout
    keys = Counter()
    den_count = Counter()
    phi_of = {}
    rmin, rmax = 1.0, 0.0
    for line in text.splitlines():
        num, den, p, q, r = map(int, line.split())
        keys[(num, den)] += 1
        den_count[den] += 1
        if den not in phi_of:
            ph = den
            for pr in primes_of(p, spf) | primes_of(q, spf) | primes_of(r, spf):
                if den % pr == 0:
                    ph = ph // pr * (pr - 1)
            phi_of[den] = ph
        x = num / den
        rmin, rmax = min(rmin, x), max(rmax, x)
    T = sum(keys.values())
    D = len(keys)
    P = sum(comb(k, 2) for k in keys.values())
    w = rmax - rmin
    # pairs sharing a denominator, and the uniform-given-denominator expectation
    same_den_pairs = sum(comb(n, 2) for n in den_count.values())
    p_null = sum(comb(n, 2) / (phi_of[d] * w) for d, n in den_count.items() if n >= 2)
    dens = np.array(sorted(den_count.elements()), float)
    H = (S / 3) ** 3
    return dict(
        S=S, triads=T, distinct_R=D, pairs=P,
        M_eff=T * T / (2 * P) if P else float("inf"),
        candidate_R=3 / pi ** 2 * H * H,          # reduced fractions in (0,1), den <= (S/3)^3
        same_den_pairs=same_den_pairs,
        P_null=p_null,
        median_den_over_S3=float(np.median(dens)) / S ** 3,
        frac_den_below_S2=float(np.mean(dens <= S * S)),
    )


def main() -> int:
    rows = load_per_s()
    S = np.array([r["S"] for r in rows])
    out_lines = []

    def say(s=""):
        print(s)
        out_lines.append(s)

    for label, field in [("cumulative pairs N(S)", "cum_pairs"),
                         ("cumulative primitive classes", "cum_prim_classes")]:
        N = np.array([r[field] for r in rows])
        say(f"== {label} ==")
        say("local exponent dlogN/dlogS:")
        for a, b, e in local_exponents(list(S), list(N),
                                       [(100, 200), (200, 400), (400, 800), (800, 1600),
                                        (1600, 3200), (2400, 4800)]):
            say(f"   {a:5d}-{b:<5d} {e:.3f}")
        for train, test in [((100, 1000), (1001, 2000)), ((100, 2400), (2401, 4800))]:
            say(f"fit on {train[0]}..{train[1]}, predict {test[0]}..{test[1]}")
            say(f"   {'model':22s} {'params':>26s} {'rms in':>8s} {'rms out':>8s} {'max out':>8s} {'end err':>8s}")
            for name, model in MODELS.items():
                th, rin, rout, end, mx = score(S, N, model, train, test)
                ps = ", ".join(f"{v:.4g}" for v in ([np.exp(th[0])] + list(th[1:])))
                say(f"   {name:22s} {ps:>26s} {rin:8.4f} {rout:8.4f} {mx:8.4f} {end:+8.4f}")
        say()

    # per-S averages: is N(S) = Theta(S)?
    per = np.array([r["pairs"] for r in rows], float)
    say("== per-S pair count, window means ==")
    say(f"   {'window':>11s} {'mean N(S)':>10s} {'mean N(S)/S':>12s} {'mean/(log S)^4':>15s}")
    for a in [200, 400, 800, 1600, 3200, 4400]:
        m = (S >= a) & (S < a + 400)
        mid = a + 200
        say(f"   {a:5d}-{a+399:<5d} {per[m].mean():10.2f} {per[m].mean()/mid:12.5f} {per[m].mean()/log(mid)**4:15.5f}")
    say()

    # birthday check
    spf = spf_table(5000)
    brows = []
    with tempfile.TemporaryDirectory() as tmp:
        exe = build_core(Path(tmp))
        for s in [200, 400, 600, 800, 1200, 1600, 2000, 3000, 4000, 4800]:
            brows.append(birthday_row(exe, s, spf))
    cols = list(brows[0].keys())
    with (DATA / "birthday.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in brows:
            w.writerow(r)
    say("== birthday heuristic check ==")
    say(f"   {'S':>5s} {'T':>8s} {'D':>8s} {'T-D':>5s} {'P':>4s} {'M_eff':>10s} {'M_eff/S^3':>9s} "
        f"{'cand.R':>9s} {'P_null':>8s} {'P/P_null':>8s} {'med d/S^3':>9s}")
    for r in brows:
        say(f"   {r['S']:5d} {r['triads']:8d} {r['distinct_R']:8d} {r['triads']-r['distinct_R']:5d} "
            f"{r['pairs']:4d} {r['M_eff']:10.3e} {r['M_eff']/r['S']**3:9.4f} {r['candidate_R']:9.2e} "
            f"{r['P_null']:8.3f} {r['pairs']/r['P_null']:8.1f} {r['median_den_over_S3']:9.4f}")
    Ss = np.array([r["S"] for r in brows], float)
    Me = np.array([r["M_eff"] for r in brows])
    Pn = np.array([r["P_null"] for r in brows])
    slope_M = np.polyfit(np.log(Ss), np.log(Me), 1)[0]
    slope_null = np.polyfit(np.log(Ss), np.log(Pn), 1)[0]
    say(f"   fitted exponent of M_eff: {slope_M:.3f}   (heuristic asserts 3)")
    say(f"   fitted exponent of P_null: {slope_null:.3f}  (uniform-given-denominator expectation)")
    (DATA / "growth_fits.txt").write_text("\n".join(out_lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
