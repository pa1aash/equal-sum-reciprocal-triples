# Trace: every checkable item of the paper

*Triples with equal sum and equal reciprocal sum*. Numbering is that of the submitted version of the paper. Each row names the item, the script that proves, checks or computes it, and the make target that runs it. A step name in the last column is the name printed by `make` (see `python3 tools/steps.py --list`).

| item | statement | script | make target (step) |
|---|---|---|---|
| Abstract, Section 1 | smallest coincidence $\{(2,8,8),(3,3,12)\}$, $S=18$, $R=3/4$ | `checks/note/check_note.py`, `theory/diophantine/check_committed.py` | `make verify` (note checks, with the enumeration 4801..6000; committed data consistent) |
| Abstract, Section 6.1 | $3.07\times10^9$ admissible triples with $S\le4800$ | `checks/classes4800/check_enum.py` | `make verify` (classes to 4800, second program) |
| Abstract, Sections 1.3, 6.3 | local growth exponent from $2.0$ at $500$ to $1.58$ at $6000$ | `checks/enumeration6000/check_enum.py` | `make verify` (enumeration to 6000, third program) |
| Lemma 1 | every coincidence is admissible | `checks/note/check_note.py` | `make verify` (note checks, with the enumeration 4801..6000) |
| Proposition 1 | classes are points on a cubic; $\Lambda=SR$ | `theory/diophantine/variety_checks.py` | `make verify` (curves of the pencil) |
| Proposition 2 | pencil identity, fibre types, maps to $E_\Lambda$, discriminant, base points of orders 2, 3, 6 | `theory/diophantine/variety_checks.py` | `make verify` (curves of the pencil) |
| Proposition 2 | orders of the base points from five lines through the flex $O$ | `checks/note/check_note.py` | `make verify` (note checks, with the enumeration 4801..6000) |
| Proposition 2 | $(1:4:4)\mapsto(-6,45)$ at $\Lambda=27/2$ | `checks/note/check_note.py` | `make verify` (note checks, with the enumeration 4801..6000) |
| Lemma 2 | odd multiples of a positive point are positive | `checks/classes4800/check_groups.py` | `make verify` (group law on the curves) |
| Theorem 1 | $(\Lambda^2-6\Lambda-3)^2-64\Lambda=(\Lambda-1)^3(\Lambda-9)$; isosceles criterion | `checks/note/check_note.py` | `make verify` (note checks, with the enumeration 4801..6000) |
| Theorem 1 | geometric progressions of order 12 (order 4 when the middle term is last); isosceles orientations of orders 6, 6, 2; positive torsion points with $S\le60$ | `checks/note/check_note.py` | `make verify` (note checks, with the enumeration 4801..6000) |
| Theorem 1 | torsion $\mathbb Z/2\times\mathbb Z/6$ at $\Lambda=27/2$ | `theory/diophantine/ranks.py` | `make verify` (ranks (PARI ellrank)) |
| Proposition 3, Corollary 1 | reciprocation is translation by $T_2$; the dual family | `theory/diophantine/variety_checks.py`, `theory/diophantine/families.py` | `make verify` (curves of the pencil; families and dual pairs) |
| Section 3 | 1753 primitive pairs with $S\le600$, 423 dual | `theory/diophantine/families.py` | `make verify` (families and dual pairs) |
| Section 3 | for the other 1330, neither $P'-P$ nor $P'+P$ has order at most 12 | `checks/classes4800/check_enum.py` | `make verify` (classes to 4800, second program) |
| Lemma 3 | $\mathcal N(X)=\sum_D\lfloor X/S_D\rfloor$ over the primitive pairs, for every $X\le4800$ | `checks/note/check_note.py` | `make verify` (note checks, with the enumeration 4801..6000) |
| Theorem 2 | the isosceles family; 506 primitive isosceles pairs with $S\le4800$ | `checks/note/check_note.py`, `theory/diophantine/families.py` | `make verify` (note checks, with the enumeration 4801..6000; families and dual pairs) |
| Theorem 2 | copies $kD_{u,v}$, $k\ge1$: $0.608X$ at $X=4800$ and $c_{\rm iso}X\log X-0.286X$; threshold $\log X>437.8$ for Theorem 3 | `checks/note/check_note.py` | `make verify` (note checks, with the enumeration 4801..6000) |
| Lemma 4, Theorem 3 | conic bundle of dual pairs; the three geometric-progression ratios; multiplicity; constants $3/(32\pi^4)$, $3/(128\pi^4)$ | `checks/classes4800/check_theorem5_family.py`, `theory/diophantine/families.py` | `make verify` (conic-bundle family; families and dual pairs) |
| Theorem 4 | $P=(4:9:18)$ on $C_{155/12}$ has $nP\ne O$ for $n\le12$ | `checks/classes4800/check_groups.py` | `make verify` (group law on the curves) |
| Section 5 | rank 0 of $C_{27/2}$ and its twelve rational points | `theory/diophantine/ranks.py` | `make verify` (ranks (PARI ellrank)) |
| Section 5 | rank of $C_{155/12}$ is 2 (`ellrank`, $r_1=r_2$) | `theory/diophantine/ranks.py` | `make verify` (ranks (PARI ellrank)) |
| Section 5 | independent 2-isogeny descents give the same ranks | `checks/classes4800/check_descent.py` | `make verify` (2-isogeny descents (PARI points)) |
| Section 5 | $3P=(162833463:723926268:287876366)$ | `checks/point3P/check_point3P.py` | `make verify` (the point 3P) |
| Section 5 | first classes of sizes 3-6 at $S=136,408,1849,4600$; no class of size 7 to 6000 | `checks/note/check_note.py`, `checks/classes4800/check_enum.py` | `make verify` (note checks, with the enumeration 4801..6000; classes to 4800, second program) |
| Table 1 | first class of each size, $\Lambda$, ranks, torsion groups | `tools/make_tables.py` (from `groups.csv`, `ranks.txt`) | `make tables` (tables 1 and 2) |
| Table 2 | cumulative counts and local exponents | `tools/make_tables.py` (from `check_enum.txt`, `per_S.csv`, `check_note.txt`) | `make tables` (tables 1 and 2) |
| Section 6.1 | two independent enumerations agree to $S=4800$; class by class against exact fractions for $S\le300$ | `theory/diophantine/enumerate_fast.py`, `checks/classes4800/check_enum.py` | `make verify`, `make full` (classes to 4800, second program; enumeration to 4800 from scratch (enumerate_fast)) |
| Section 6.1 | a third enumeration to $S=6000$ | `checks/enumeration6000/check_enum.py` | `make verify`, `make full` (enumeration to 6000, third program; third program from scratch to 6000) |
| Section 6.1 | the second program on $4801\le S\le6000$: counts at 5400 and 6000, 62401 primitive classes | `checks/note/check_note.py --enum` | `make verify` (note checks, with the enumeration 4801..6000) |
| Section 6.1 | 582 sums $19\le S\le600$, 507 with a primitive class | `checks/note/count_sums.py` | `make verify` (sums with a primitive class) |
| Section 6.1 | class-size histogram to 4800 | `checks/classes4800/check_enum.py` | `make verify` (classes to 4800, second program) |
| Section 6.1 | dual pairs: 32.9% to 200, 16.5% of 46 254 to 4800; 669 non-dual pairs to 400 | `theory/diophantine/families.py` | `make verify` (families and dual pairs) |
| Section 1.4 | $R/r=(\Lambda-1)/4$ under the Ravi substitution | `checks/note/check_note.py` | `make verify` (note checks, with the enumeration 4801..6000) |
| Proposition 4 | factorisation identity; recount of $n(S)$ for $S\le160$ | `checks/note/check_note.py` | `make verify` (note checks, with the enumeration 4801..6000) |
| Section 6.3 | fits of the exponent $q$ and their errors | `theory/diophantine/exponent_fits.py` | `make verify` (exponent fits) |
| Section 6.3 | a quadratic law is not supported | `theory/diophantine/growth_fits.py` | `make verify` (growth fits) |
| Section 6.3 | $\mathcal N(S)/S^2$ from $0.0094$ to $0.0040$ | `checks/enumeration6000/check_enum.py` | `make verify` (enumeration to 6000, third program) |
| Section 6.4 | 189 lines; 84% of primitive pairs to 4800 are not dual | `theory/diophantine/families.py` | `make verify` (families and dual pairs) |
| Section 6.4 | pairs of points and the fibre square | `theory/diophantine/variety_checks.py` | `make verify` (curves of the pencil) |
| Figure 1 | the curve $C_{27/2}$ and the growth of the counts | `figures/src/F9.py` | `make figures` (figure 1 (F9)) |
