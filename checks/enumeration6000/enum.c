/* Exhaustive enumeration of two-coefficient collisions among hyperbolic triads.
 *
 * For every S in [Smin, Smax], list all 2 <= p <= q <= r with p+q+r = S and
 * pq+qr+rp < pqr (i.e. R < 1), and group them by the exact value of
 * R = (pq+qr+rp)/(pqr).  Grouping: sort by the double value of e2/e3 (IEEE division
 * is correctly rounded, so equal rationals give bit-identical doubles), then split each
 * run of equal doubles by exact cross-multiplication in 128-bit integers.  Distinct
 * rationals that happen to share a double are therefore never merged.
 *
 * Usage: enum Smin Smax summary.txt fibres.txt dumpmax
 *   summary line per S:
 *     S ntriads nfibres npairs maxfibre adjpairs nonadjpairs primfibres primpairs
 *   adjpairs     = pairs whose least orders differ by exactly 1
 *   primfibres   = fibres whose members have overall gcd 1
 *   primpairs    = pairs {T,T'} NOT of the form k*(U,U') with k>1 and U,U' hyperbolic
 *   fibres.txt (only for S <= dumpmax):  S num den p q r | p q r | ...
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct { double v; uint64_t e2, e3; int p, q, r; } T;

static int cmp(const void *a, const void *b) {
    double x = ((const T *)a)->v, y = ((const T *)b)->v;
    return (x > y) - (x < y);
}
static uint64_t g64(uint64_t a, uint64_t b) { while (b) { uint64_t t = a % b; a = b; b = t; } return a; }
static int eqR(const T *a, const T *b) {
    return (unsigned __int128)a->e2 * b->e3 == (unsigned __int128)b->e2 * a->e3;
}
static int hyp(uint64_t p, uint64_t q, uint64_t r) { return p*q + q*r + r*p < p*q*r; }

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage\n"); return 2; }
    int Smin = atoi(argv[1]), Smax = atoi(argv[2]);
    int dumpmax = atoi(argv[5]);
    FILE *fs = fopen(argv[3], "w"), *ff = fopen(argv[4], "w");
    if (!fs || !ff) return 2;
    size_t cap = (size_t)Smax * Smax / 12 + 1000;
    T *a = malloc(cap * sizeof(T));
    int *used = malloc(4096 * sizeof(int));
    T *grp = malloc(4096 * sizeof(T));
    for (int S = Smin; S <= Smax; S++) {
        size_t n = 0;
        for (int p = 2; 3 * p <= S; p++)
            for (int q = p; 2 * q <= S - p; q++) {
                int r = S - p - q;
                uint64_t e2 = (uint64_t)p*q + (uint64_t)q*r + (uint64_t)r*p, e3 = (uint64_t)p*q*r;
                if (e2 >= e3) continue;
                if (n >= cap) { fprintf(stderr, "cap\n"); return 3; }
                a[n].v = (double)e2 / (double)e3; a[n].e2 = e2; a[n].e3 = e3;
                a[n].p = p; a[n].q = q; a[n].r = r; n++;
            }
        qsort(a, n, sizeof(T), cmp);
        long long nfib = 0, npairs = 0, adj = 0, nonadj = 0, primfib = 0, primpairs = 0; int maxk = 1;
        size_t i = 0;
        while (i < n) {
            size_t j = i + 1;
            while (j < n && a[j].v == a[i].v) j++;
            size_t run = j - i;
            if (run > 1) {
                if (run > 4096) { fprintf(stderr, "run too long\n"); return 4; }
                for (size_t u = 0; u < run; u++) used[u] = 0;
                for (size_t u = 0; u < run; u++) {
                    if (used[u]) continue;
                    int kk = 0;
                    for (size_t w = u; w < run; w++)
                        if (!used[w] && eqR(&a[i + u], &a[i + w])) { used[w] = 1; grp[kk++] = a[i + w]; }
                    if (kk < 2) continue;
                    nfib++; npairs += (long long)kk * (kk - 1) / 2; if (kk > maxk) maxk = kk;
                    uint64_t g = 0;
                    for (int x = 0; x < kk; x++) { g = g64(g, grp[x].p); g = g64(g, grp[x].q); g = g64(g, grp[x].r); }
                    if (g == 1) primfib++;
                    for (int x = 0; x < kk; x++) for (int y = x + 1; y < kk; y++) {
                        int d = grp[y].p - grp[x].p; if (d < 0) d = -d;
                        if (d == 1) adj++; else nonadj++;
                        uint64_t h = g64(g64(g64(grp[x].p, grp[x].q), g64(grp[x].r, grp[y].p)), g64(grp[y].q, grp[y].r));
                        int scaled = 0;
                        for (uint64_t k = 2; k <= h; k++) if (h % k == 0) {
                            if (hyp(grp[x].p / k, grp[x].q / k, grp[x].r / k) && hyp(grp[y].p / k, grp[y].q / k, grp[y].r / k)) { scaled = 1; break; }
                        }
                        if (!scaled) primpairs++;
                    }
                    if (S <= dumpmax) {
                        uint64_t gg = g64(grp[0].e2, grp[0].e3);
                        fprintf(ff, "%d %llu %llu", S, (unsigned long long)(grp[0].e2 / gg), (unsigned long long)(grp[0].e3 / gg));
                        for (int x = 0; x < kk; x++) fprintf(ff, " | %d %d %d", grp[x].p, grp[x].q, grp[x].r);
                        fprintf(ff, "\n");
                    }
                }
            }
            i = j;
        }
        fprintf(fs, "%d %zu %lld %lld %d %lld %lld %lld %lld\n", S, n, nfib, npairs, maxk, adj, nonadj, primfib, primpairs);
    }
    fclose(fs); fclose(ff);
    return 0;
}
