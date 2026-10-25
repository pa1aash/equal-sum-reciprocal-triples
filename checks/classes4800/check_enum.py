"""Checks over all degeneracy classes to S = 4800: Sections 3 and 6.1 of the note.

Uses the C enumerator enum_classes.c (compiled and run if its output is missing), after
cross-checking it against an independent pure-Python enumeration (Fraction keys) for S <= 300.

Run from the repository root:
    python3 checks/classes4800/check_enum.py
"""
import os
import subprocess
import sys
from fractions import Fraction as Fr
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from dio_common import (add, neg, order, normalize, lam_of, check, finish, Fail)  # noqa: E402

OUT = os.path.join(HERE, "check_enum.txt")
SMAX = 4800
ENUM = os.path.join(HERE, "enum_4800.txt")
log = []


def run_c(smin, smax, listmax, path=None):
    exe = os.path.join(HERE, "enum_classes")
    if not os.path.exists(exe):
        subprocess.check_call(["cc", "-O2", "-o", exe, os.path.join(HERE, "enum_classes.c")])
    out = subprocess.check_output([exe, str(smin), str(smax), str(listmax)], text=True)
    if path:
        with open(path, "w") as fh:
            fh.write(out)
    return out


def parse(text):
    stats, classes = {}, []
    for line in text.splitlines():
        f = line.split()
        if f[0] == "S":
            stats[int(f[1])] = dict(ntr=int(f[3]), pairs=int(f[5]), classes=int(f[7]), maxsize=int(f[9]))
        elif f[0] == "C":
            S, size, prim = int(f[1]), int(f[2]), int(f[3])
            mem = [tuple(int(v) for v in m.split(",")) for m in f[6:]]
            classes.append(dict(S=S, size=size, prim=prim, R=Fr(int(f[4]), int(f[5])), mem=mem))
    return stats, classes


def python_classes(S):
    groups = {}
    for p in range(2, S // 3 + 1):
        for q in range(p, (S - p) // 2 + 1):
            r = S - p - q
            R = Fr(1, p) + Fr(1, q) + Fr(1, r)
            if R < 1:
                groups.setdefault(R, []).append((p, q, r))
    return {R: sorted(v) for R, v in groups.items() if len(v) >= 2}


def gcd_all(trips):
    g = 0
    for t in trips:
        for x in t:
            g = gcd(g, x)
    return g


def is_dual(t, u):
    a, b, c = t
    d = sorted((b * c, c * a, a * b))
    g = gcd_all([d])
    d = [x // g for x in d]
    gu = gcd_all([u])
    return d == sorted(x // gu for x in u)


def main():
    try:
        # ---- 1. C enumerator vs Python on S <= 300 -------------------------------------
        st_small, cl_small = parse(run_c(10, 300, 300))
        for S in range(10, 301):
            py = python_classes(S)
            cc = {c["R"]: sorted(c["mem"]) for c in cl_small if c["S"] == S}
            if py != cc:
                raise Fail(f"C/Python mismatch at S={S}")
        log.append("PASS C enumerator and independent Python enumeration agree class-by-class for 10 <= S <= 300")

        if not os.path.exists(ENUM):
            run_c(10, SMAX, SMAX, ENUM)
        with open(ENUM) as fh:
            stats, classes = parse(fh.read())
        check(max(stats) == SMAX and len(stats) == SMAX - 9, f"full enumeration present for 10 <= S <= {SMAX}", log)
        tot_tr = sum(v["ntr"] for v in stats.values())
        log.append(f"     {tot_tr} hyperbolic triples, {len(classes)} classes of size >= 2, "
                   f"{sum(v['pairs'] for v in stats.values())} pairs")

        # ---- 2. small sums -----------------------------------------------
        check(all(stats[S]["pairs"] == 0 for S in range(10, 18)) and stats[18]["pairs"] == 1,
              "no degeneracy for S <= 17, exactly one at S = 18", log)
        c18 = [c for c in classes if c["S"] == 18][0]
        check(sorted(c18["mem"]) == [(2, 8, 8), (3, 3, 12)], "S=18 class is {(2,8,8),(3,3,12)}", log)
        c36 = sorted(sorted(c["mem"]) for c in classes if c["S"] == 36)
        check(c36 == [[(4, 16, 16), (6, 6, 24)], [(6, 15, 15), (8, 8, 20)]], "S=36: the two classes", log)
        cum, run = {}, 0
        for S in range(10, SMAX + 1):
            run += stats[S]["pairs"]
            cum[S] = run
        table = {18: 1, 100: 92, 200: 386, 300: 840, 400: 1496, 500: 2210, 600: 3067}
        got = {S: cum[S] for S in table}
        log.append(f"     cumulative pair counts N(S): {got}")
        check(got == table, "table of N(S) reproduced exactly", log)

        # ---- 3. first fibres of each size ----------------------------------------
        first = {}
        for c in classes:
            k = c["size"]
            if k not in first or c["S"] < first[k]["S"]:
                first[k] = c
        firstS = {k: first[k]["S"] for k in sorted(first)}
        log.append(f"     first S with a class of size k: {firstS}")
        for k in range(3, 7):
            log.append(f"     size {k}: S={first[k]['S']}, R={first[k]['R']}, lambda={first[k]['R'] * first[k]['S']}, "
                       f"primitive={bool(first[k]['prim'])}, members {sorted(first[k]['mem'])}")
        check(firstS[3] == 136 and firstS[4] == 408 and firstS[5] == 1849 and firstS[6] == 4600,
              "first sizes 3,4,5,6 occur at S = 136, 408, 1849, 4600", log)
        check(max(firstS) == 6, f"no class of size >= 7 for S <= {SMAX}", log)
        claimed = {
            3: [(15, 55, 66), (16, 40, 80), (17, 34, 85)],
            5: [(168, 820, 861), (172, 645, 1032), (185, 480, 1184), (215, 344, 1290), (253, 276, 1320)],
            6: [(750, 1750, 2100), (756, 1674, 2170), (800, 1400, 2400), (805, 1380, 2415), (882, 1170, 2548),
                (920, 1104, 2576)],
        }
        for k, mem in claimed.items():
            check(sorted(first[k]["mem"]) == mem, f"size-{k} example matches the class found", log)
        c408 = first[4]
        check(c408["R"] * 408 == Fr(68, 5), "S=408 size-4 class lies on lambda = 68/5 (the S=136 curve)", log)
        scaled = sorted(tuple(3 * x for x in t) for t in claimed[3])
        check(all(t in c408["mem"] for t in scaled), "S=408 class = 3 x (S=136 class) plus one more point: "
              + str([t for t in sorted(c408['mem']) if t not in scaled]), log)
        sizes = {}
        for c in classes:
            sizes[c["size"]] = sizes.get(c["size"], 0) + 1
        log.append(f"     class-size histogram S <= {SMAX}: {dict(sorted(sizes.items()))}")
        psizes = {}
        for c in classes:
            if c["prim"]:
                psizes[c["size"]] = psizes.get(c["size"], 0) + 1
        log.append(f"     primitive-class size histogram S <= {SMAX}: {dict(sorted(psizes.items()))}")

        # ---- 4. statistics: primitive pairs with S <= 600 ---------------------------
        prim_pairs, dual_pairs, other = 0, 0, []
        for c in classes:
            if c["S"] > 600:
                continue
            m = c["mem"]
            for i in range(len(m)):
                for j in range(i + 1, len(m)):
                    if gcd_all([m[i], m[j]]) == 1:
                        prim_pairs += 1
                        if is_dual(m[i], m[j]):
                            dual_pairs += 1
                        else:
                            other.append((c["S"], m[i], m[j]))
        log.append(f"     primitive pairs (gcd of the six entries = 1) with S <= 600: {prim_pairs}; dual: {dual_pairs}; other: {len(other)}")
        check((prim_pairs, dual_pairs, len(other)) == (1753, 423, 1330), "counts 1753 = 423 + 1330 reproduced", log)
        tors_rel = 0
        for S, t, u in other:
            lam = lam_of(t)
            P, Q = normalize(t), normalize(u)
            if order(lam, add(lam, Q, neg(lam, P)), 12) is not None or order(lam, add(lam, Q, P), 12) is not None:
                tors_rel += 1
        check(tors_rel == 0, "none of the 1330 non-dual primitive pairs has P'-P or P'+P torsion (Mazur bound 12)", log)

        # ---- 5. isosceles family ----------------------------------------------------
        pairset = {}
        for c in classes:
            pairset.setdefault(c["S"], []).append(set(c["mem"]))
        D = []
        for v in range(2, 200):
            for u in range(1, v):
                if gcd(u, v) != 1:
                    continue
                g = gcd(2 * u + v, u + 2 * v)
                S = (2 * u + v) * (u + 2 * v) // g
                if S <= SMAX:
                    t1 = tuple(sorted(x * (2 * u + v) // g for x in (u, v, v)))
                    t2 = tuple(sorted(x * (u + 2 * v) // g for x in (v, u, u)))
                    D.append((u, v, g, S, t1, t2))
        check(len(D) == 506, f"number of primitive isosceles pairs D_(u,v) with S <= 4800: {len(D)} (claimed 506)", log)
        check(all(g in (1, 3) for (_, _, g, _, _, _) in D), "g in {1,3} for all of them", log)
        check(all(gcd_all([t1, t2]) == 1 for (*_, t1, t2) in D), "all D_(u,v) primitive", log)
        missing, copies = 0, 0
        for (u, v, g, S, t1, t2) in D:
            for k in range(1, SMAX // S + 1):
                a = tuple(k * x for x in t1)
                b = tuple(k * x for x in t2)
                if Fr(1, a[0]) + Fr(1, a[1]) + Fr(1, a[2]) >= 1:
                    continue
                copies += 1
                if not any(a in cs and b in cs for cs in pairset.get(k * S, [])):
                    missing += 1
        check(missing == 0, f"every hyperbolic copy kD_(u,v) with kS <= 4800 ({copies} pairs) is found by the enumeration", log)
        # a class contains at most one isosceles pair
        multi = 0
        for c in classes:
            iso = [t for t in c["mem"] if t[0] == t[1] or t[1] == t[2]]
            if len(iso) > 2:
                multi += 1
        check(multi == 0, "no class contains more than two isosceles triples (at most one isosceles pair per class)", log)
        log.append("ALL ENUMERATION CHECKS PASSED")
        finish(log, OUT)
    except Fail:
        finish(log, OUT)
        sys.exit(1)


if __name__ == "__main__":
    main()
