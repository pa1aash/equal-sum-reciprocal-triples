#!/usr/bin/env python3
"""Third enumeration of two-coefficient collisions, to S = 6000.  Exits nonzero on failure.

1. Pure-Python exact enumeration (reduced integer fractions) of every hyperbolic triad
   with S <= 600; collision fibres compared line by line with the C enumerator (enum.c).
2. C enumerator over 10 <= S <= SMAX_C (default 6000), run in parallel chunks; outputs in
   checks/enumeration6000/data/.
3. Claims: Cor 2 / Thm A,B; tab:enum listing; S=36 classes; tab:density in both counting
   conventions; power-law exponent; growth to SMAX_C; first collisions.
"""
import os
import sys
import math
import subprocess
from fractions import Fraction as Fr
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SMAX_C = int(os.environ.get("SMAX_C", "6000"))
DUMPMAX = 2000
FAIL = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAIL.append(msg)


# ---------------------------------------------------------------- C enumerator
def run_c():
    os.makedirs(DATA, exist_ok=True)
    summ = os.path.join(DATA, f"summary_10_{SMAX_C}.txt")
    fib = os.path.join(DATA, f"fibres_10_{DUMPMAX}.txt")
    if os.path.exists(summ) and os.path.exists(fib):
        return summ, fib
    exe = os.path.join(DATA, "enum")
    subprocess.run(["cc", "-O2", "-o", exe, os.path.join(HERE, "enum.c")], check=True)
    nch = 10
    b = [10] + [round(SMAX_C * (k / nch) ** (1 / 3)) for k in range(1, nch + 1)]
    procs, parts = [], []
    for i in range(nch):
        lo = b[i] if i == 0 else b[i] + 1
        hi = b[i + 1]
        ps, pf = os.path.join(DATA, f"_s{i}.txt"), os.path.join(DATA, f"_f{i}.txt")
        parts.append((ps, pf))
        procs.append(subprocess.Popen([exe, str(lo), str(hi), ps, pf, str(DUMPMAX)]))
    for pr in procs:
        assert pr.wait() == 0
    with open(summ, "w") as fs, open(fib, "w") as ff:
        for ps, pf in parts:
            fs.write(open(ps).read()); ff.write(open(pf).read())
            os.remove(ps); os.remove(pf)
    return summ, fib


summ_path, fib_path = run_c()
summary = {}
for line in open(summ_path):
    v = list(map(int, line.split()))
    summary[v[0]] = dict(ntri=v[1], nfib=v[2], npairs=v[3], maxk=v[4], adj=v[5], nonadj=v[6], primfib=v[7], primpairs=v[8])
check(sorted(summary) == list(range(10, SMAX_C + 1)), f"C summary covers 10..{SMAX_C}")
cfib = {}
for line in open(fib_path):
    parts = line.split("|")
    S, num, den = map(int, parts[0].split())
    mem = tuple(sorted(tuple(map(int, x.split())) for x in parts[1:]))
    cfib.setdefault(S, set()).add(((num, den), mem))

# ---------------------------------------------------------------- pure Python S <= 600
print("== Python exact enumeration S <= 600")
PY = 600
pyfib = {}
ntri_py = {}
for S in range(10, PY + 1):
    groups = {}
    n = 0
    for p in range(2, S // 3 + 1):
        for q in range(p, (S - p) // 2 + 1):
            r = S - p - q
            e2, e3 = p * q + q * r + r * p, p * q * r
            if e2 < e3:
                g = gcd(e2, e3)
                groups.setdefault((e2 // g, e3 // g), []).append((p, q, r))
                n += 1
    ntri_py[S] = n
    pyfib[S] = {(k, tuple(sorted(v))) for k, v in groups.items() if len(v) > 1}
check(all(ntri_py[S] == summary[S]["ntri"] for S in range(10, PY + 1)), "triad counts agree (Python vs C), S<=600")
check(all(pyfib[S] == cfib.get(S, set()) for S in range(10, PY + 1)), "collision fibres agree exactly (Python vs C), S<=600")
# cross-check C pair counts from fibres
check(all(summary[S]["npairs"] == sum(len(m) * (len(m) - 1) // 2 for _, m in cfib.get(S, ())) for S in range(10, DUMPMAX + 1)),
      "C pair counts consistent with dumped fibres, S<=2000")

# ---------------------------------------------------------------- Cor 2 / Thm A,B / tab:enum
print("== Theorem A/B, Corollary 2, tab:enum")
check(all(not pyfib[S] for S in range(10, 18)), "no collision for S<=17")
check(pyfib[18] == {((3, 4), ((2, 8, 8), (3, 3, 12)))}, "S=18: exactly one collision, {(2,8,8),(3,3,12)}, R=3/4")
with open(os.path.join(HERE, "tab_enum_reference.txt"), "w") as fo:
    fo.write("# every hyperbolic triad with 10<=S<=18, R exact (listing of the triads of small sum)\n")
    tot = 0
    for S in range(10, 19):
        trs = [(p, q, S - p - q) for p in range(2, S // 3 + 1) for q in range(p, (S - p) // 2 + 1)
               if p * q + q * (S - p - q) + (S - p - q) * p < p * q * (S - p - q)]
        Rs = [Fr(p * q + q * r + r * p, p * q * r) for p, q, r in trs]
        tot += len(trs)
        fo.write(f"S={S} ({len(trs)} triads): " + "; ".join(f"({p},{q},{r}) {R}" for (p, q, r), R in zip(trs, Rs)) + "\n")
        if S <= 17:
            assert len(set(Rs)) == len(Rs)
    fo.write(f"total {tot}\n")
print(f"   tab:enum reference listing written (total triads 10<=S<=18: {tot})")
counts = {S: summary[S]['ntri'] for S in range(10, 19)}
print("   triads per S:", counts)

# ---------------------------------------------------------------- S = 36
print("== S=36")
f36 = sorted(cfib[36])
for k_, m_ in f36:
    print("   ", k_, m_)
check(f36 == sorted({((3, 10), ((6, 15, 15), (8, 8, 20))), ((3, 8), ((4, 16, 16), (6, 6, 24)))}),
      "S=36: exactly two fibres {(4,16,16),(6,6,24)} R=3/8 and {(6,15,15),(8,8,20)} R=3/10")
print("   note: (6,15,15) has least order 6 and (8,8,20) least order 8: a NON-adjacent collision")

# ---------------------------------------------------------------- density table
print("== tab:density (cumulative)")
claimed = {18: 1, 100: 92, 200: 386, 300: 840, 400: 1496, 500: 2210, 600: 3067}
cumF = {}; cumP = {}; cumPF = {}; cumPP = {}; cumA = {}
cf = cp = cpf = cpp = ca = 0
for S in range(10, SMAX_C + 1):
    d = summary[S]
    cf += d["nfib"]; cp += d["npairs"]; cpf += d["primfib"]; cpp += d["primpairs"]; ca += d["adj"]
    cumF[S], cumP[S], cumPF[S], cumPP[S], cumA[S] = cf, cp, cpf, cpp, ca
print("   S  | claimed | pairs | fibres | primitive fibres | adjacent-strata pairs")
for S in claimed:
    print(f"   {S:4d} | {claimed[S]:5d} | {cumP[S]:5d} | {cumF[S]:5d} | {cumPF[S]:5d} | {cumA[S]:5d}")
check(all(cumP[S] == claimed[S] for S in claimed), "tab:density matches the PAIR convention (C(k,2) pairs per fibre)")
check(cumF[600] == 2977, "fibre (class) count to 600 = 2977 ")
check(all(cumF[S] <= cumP[S] for S in claimed), "fibres <= pairs")
# ratios printed in tab:density
ok = True
for S in claimed:
    r1 = round(cumP[S] / S, 2 if S > 18 else 3); r2 = round(cumP[S] / S ** 2, 4)
    print(f"   {S}: N/S={cumP[S]/S:.4f} N/S^2={cumP[S]/S**2:.5f}")
# scaled copies at 600
check(600 // 18 == 33, "floor(600/18)=33")
# total adjacency statistics
tot_adj = sum(summary[S]["adj"] for S in range(10, 601))
tot_non = sum(summary[S]["nonadj"] for S in range(10, 601))
print(f"   S<=600: adjacent-strata pairs {tot_adj}, non-adjacent pairs {tot_non}")
check(tot_non > 0, "non-adjacent collisions exist (contradicts 'localizes every degeneracy to adjacent strata')")

# ---------------------------------------------------------------- power law
print("== power-law fit, 50<=S<=600")


def ols(xs, ys):
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


for name, cum in (("pairs", cumP), ("fibres", cumF)):
    xs = [math.log(S) for S in range(50, 601)]
    ys = [math.log(cum[S]) for S in range(50, 601)]
    sl = ols(xs, ys)
    xs2 = [math.log(S) for S in claimed if S >= 50]
    ys2 = [math.log(cum[S]) for S in claimed if S >= 50]
    sl2 = ols(xs2, ys2)
    print(f"   {name}: slope over all integer S = {sl:.4f}; over table checkpoints = {sl2:.4f}")
    if name == "pairs":
        slope_pairs = sl
check(abs(slope_pairs - 2.03) < 0.02, f"exponent 2.03 reproduced (pairs, all S in [50,600]): {slope_pairs:.4f}")

# ---------------------------------------------------------------- conjecture conj:density
print(f"== conj:density test to S={SMAX_C}")
print("   S     N_pairs   N/S^2     N_fibres  N/S^2   local exponent (pairs, [S/2,S])  N_pairs/(S^2 log S)")
growth = []
for S in [100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000, 2500, 3000, 3600, 4000, 4800, 5400, 6000]:
    if S > SMAX_C:
        continue
    le = math.log(cumP[S] / cumP[S // 2]) / math.log(2)
    print(f"   {S:5d} {cumP[S]:9d} {cumP[S]/S**2:.5f}  {cumF[S]:9d} {cumF[S]/S**2:.5f}   {le:.4f}   {cumP[S]/(S*S*math.log(S)):.6f}")
    growth.append((S, cumP[S] / S ** 2))
# OLS exponent on [S/2, S] windows at top of range
for lo, hi in ((600, 1200), (1200, 2400), (2400, 4800), (3000, SMAX_C)):
    if hi <= SMAX_C:
        xs = [math.log(S) for S in range(lo, hi + 1)]
        print(f"   OLS exponent (pairs) on [{lo},{hi}]: {ols(xs, [math.log(cumP[S]) for S in range(lo, hi + 1)]):.4f};"
              f" fibres {ols(xs, [math.log(cumF[S]) for S in range(lo, hi + 1)]):.4f}")
# per-sum N(S)/S averaged over blocks
for lo, hi in ((500, 600), (1000, 1100), (2000, 2100), (3000, 3100), (4000, 4100), (5900, 6000)):
    if hi <= SMAX_C:
        av = sum(summary[S]["npairs"] for S in range(lo, hi + 1)) / (hi - lo + 1)
        print(f"   mean N(S) on [{lo},{hi}] = {av:.2f};  /S = {av/((lo+hi)/2):.5f};  /(S log S) = {av/(((lo+hi)/2)*math.log((lo+hi)/2)):.6f}")

# ---------------------------------------------------------------- first collisions
print("== first adjacent collision table")


def Sstar(pp):
    return 18 if pp == 2 else 19 if pp == 3 else 3 * pp + 8 if pp <= 8 else 3 * pp + 7


def stratum(S, pp):
    out = []
    for q in range(pp, (S - pp) // 2 + 1):
        r = S - pp - q
        if pp * q + q * r + r * pp < pp * q * r:
            out.append((Fr(pp * q + q * r + r * pp, pp * q * r), (pp, q, r)))
    return out


def window(S, pp):
    A, B = stratum(S, pp), stratum(S, pp + 1)
    lo = min(x[0] for x in A); hi = max(x[0] for x in B)
    return [x for x in A if lo <= x[0] <= hi], [x for x in B if lo <= x[0] <= hi]


first = {}
for S in range(10, DUMPMAX + 1):
    for (_, mem) in cfib.get(S, ()):
        for i in range(len(mem)):
            for j in range(i + 1, len(mem)):
                a, b = mem[i][0], mem[j][0]
                if abs(a - b) == 1:
                    pp = min(a, b)
                    if pp not in first:
                        x, y = (mem[i], mem[j]) if a < b else (mem[j], mem[i])
                        first[pp] = (S, x, y)
claimedTH5 = {2: (18, 18, 0, "1+1", "1+1", (2, 8, 8), (3, 3, 12)), 3: (19, 38, 19, "2+1", "11+3", (3, 14, 21), (4, 6, 28)),
              4: (20, 20, 0, "1+1", "1+1", (4, 8, 8), (5, 5, 10)), 5: (23, 117, 94, "1+1", "49+12", (5, 32, 80), (6, 15, 96)),
              6: (26, 34, 8, "2+1", "6+3", (6, 14, 14), (7, 9, 18)), 7: (29, 62, 33, "2+1", "18+8", (7, 20, 35), (8, 14, 40)),
              8: (32, 64, 32, "2+1", "18+8", (8, 20, 36), (9, 15, 40)), 9: (34, 42, 8, "1+1", "5+3", (9, 15, 18), (10, 12, 20)),
              10: (37, 109, 72, "1+1", "37+18", (10, 44, 55), (11, 28, 70)), 11: (40, 66, 26, "1+1", "14+8", (11, 22, 33), (12, 18, 36)),
              12: (43, 188, 145, "1+1", "74+34", (12, 72, 104), (13, 45, 130)), 13: (46, 94, 48, "1+1", "25+15", (13, 39, 42), (14, 28, 52)),
              14: (49, 69, 20, "1+1", "11+7", (14, 20, 35), (15, 18, 36)), 22: (73, 422, 349, "1+1", "176+96", (22, 92, 308), (23, 77, 322))}
allrows = True
for pp, (ss, fc, gp, w1, w2, t1, t2) in claimedTH5.items():
    S0, x, y = first[pp]
    A, B = window(ss, pp); A2, B2 = window(S0, pp)
    mine = (Sstar(pp), S0, S0 - Sstar(pp), f"{len(A)}+{len(B)}", f"{len(A2)}+{len(B2)}", x, y)
    row_ok = mine == (ss, fc, gp, w1, w2, t1, t2)
    allrows &= row_ok
    print(f"   p={pp:2d}: mine {mine}  {'OK' if row_ok else 'MISMATCH claimed ' + str((ss, fc, gp, w1, w2, t1, t2))}")
    # is the listed colliding pair the unique adjacent pair at the first-collision sum?
    others = [m for (_, m) in cfib[S0] if any(abs(u[0] - v[0]) == 1 and min(u[0], v[0]) == pp for u in m for v in m)]
    if len(others) != 1 or len(others[0]) != 2:
        print(f"      note: at S={S0} adjacent fibres for p={pp}: {others}")
check(allrows, "every row reproduced (S*, first collision, gap, windows, pair)")
p_ok = [pp for pp in first if pp not in (2, 4) and first[pp][0] <= Sstar(pp)]
print(f"   adjacent pairs p with a first collision found (S<={DUMPMAX}): {len(first)} values, p up to {max(first)}")
check(not p_ok, f"Prop 3(2): first adjacent collision strictly after S*(p) for every p not in {{2,4}} found (S<={DUMPMAX})")
check(first[2][0] == 18 and first[4][0] == 20, "p=2,4: first collision at S*")

print()
if FAIL:
    print(f"{len(FAIL)} FAILURE(S):", *FAIL, sep="\n  ")
    sys.exit(1)
print("ALL ENUMERATION CHECKS PASSED")
