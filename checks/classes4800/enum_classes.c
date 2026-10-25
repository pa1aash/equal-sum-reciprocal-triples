/* Exact enumeration of two-coefficient degeneracy classes.
 *
 * For every S in [SMIN, SMAX] enumerate all hyperbolic triples 2 <= p <= q <= r, p+q+r = S,
 * R = 1/p+1/q+1/r < 1, and group them by the exact rational R = e2/e3.
 * Grouping: the double e2/e3 is the correctly rounded value of the exact rational
 * (e2, e3 < 2^53), so equal rationals give bit-identical doubles; a hash on the double bits
 * finds all candidates, and every candidate match is confirmed by exact 128-bit
 * cross-multiplication e2*e3' == e2'*e3.  Distinct rationals with the same double are kept
 * apart (probing continues).  No pair can be missed.
 *
 * Output (stdout), one line per S:
 *   S <S> ntriples <n> pairs <N(S)> classes <c> maxsize <m>
 * and one line per class of size >= 2 (all S if LIST=1, else only size >= 3):
 *   C <S> <size> <prim> <e2/g> <e3/g> p,q,r p,q,r ...
 * where prim = 1 iff the gcd of all entries of the class is 1.
 *
 * Usage: enum_classes SMIN SMAX LISTMAX   (classes of size 2 are listed for S <= LISTMAX)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef unsigned __int128 u128;
typedef struct { uint64_t key; uint32_t p, q; uint32_t stamp; int32_t cls; } Ent;

static uint64_t g64(uint64_t a, uint64_t b) { while (b) { uint64_t t = a % b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s SMIN SMAX LISTMAX\n", argv[0]); return 2; }
    long SMIN = atol(argv[1]), SMAX = atol(argv[2]), LISTMAX = atol(argv[3]);
    size_t TB = 1; while (TB < (size_t)(SMAX * SMAX / 4 + 16)) TB <<= 1;  /* > 3 x max triples */
    Ent *T = calloc(TB, sizeof(Ent));
    /* class storage */
    size_t CAP = 1 << 20;
    uint32_t *cp = malloc(CAP * sizeof(uint32_t)), *cq = malloc(CAP * sizeof(uint32_t));
    int32_t *cnext = malloc(CAP * sizeof(int32_t));
    int32_t *chead = malloc(CAP * sizeof(int32_t)), *csize = malloc(CAP * sizeof(int32_t));
    uint32_t *crp = malloc(CAP * sizeof(uint32_t)), *crq = malloc(CAP * sizeof(uint32_t));
    if (!T || !cp || !cq || !cnext || !chead || !csize || !crp || !crq) { fprintf(stderr, "oom\n"); return 3; }
    uint32_t stamp = 0;
    for (long S = SMIN; S <= SMAX; S++) {
        stamp++;
        long ntr = 0; int ncls = 0; size_t nmem = 0;
        for (long p = 2; 3 * p <= S; p++) {
            for (long q = p; 2 * q <= S - p; q++) {
                long r = S - p - q;
                uint64_t e2 = (uint64_t)(p * q + q * r + r * p), e3 = (uint64_t)p * q * r;
                if (e2 >= e3) continue;               /* not hyperbolic */
                ntr++;
                double R = (double)e2 / (double)e3;
                uint64_t key; memcpy(&key, &R, 8);
                size_t h = (size_t)((key * 0x9E3779B97F4A7C15ULL) >> 20) & (TB - 1);
                for (;;) {
                    Ent *E = &T[h];
                    if (E->stamp != stamp) {          /* empty slot: insert */
                        E->stamp = stamp; E->key = key; E->p = (uint32_t)p; E->q = (uint32_t)q; E->cls = -1;
                        break;
                    }
                    if (E->key == key) {
                        uint64_t P2 = E->p, Q2 = E->q, R2 = (uint64_t)S - P2 - Q2;
                        uint64_t f2 = P2 * Q2 + Q2 * R2 + R2 * P2, f3 = P2 * Q2 * R2;
                        if ((u128)e2 * f3 == (u128)f2 * e3) {   /* exact equality */
                            if (E->cls < 0) {
                                if ((size_t)ncls + 2 >= CAP) { fprintf(stderr, "class overflow\n"); return 4; }
                                E->cls = ncls; chead[ncls] = -1; csize[ncls] = 1;
                                crp[ncls] = E->p; crq[ncls] = E->q; ncls++;
                            }
                            if (nmem + 2 >= CAP) { fprintf(stderr, "member overflow\n"); return 4; }
                            cp[nmem] = (uint32_t)p; cq[nmem] = (uint32_t)q;
                            cnext[nmem] = chead[E->cls]; chead[E->cls] = (int32_t)nmem; nmem++;
                            csize[E->cls]++;
                            break;
                        }
                    }
                    h = (h + 1) & (TB - 1);
                }
            }
        }
        long pairs = 0; int maxs = 0;
        for (int c = 0; c < ncls; c++) { pairs += (long)csize[c] * (csize[c] - 1) / 2; if (csize[c] > maxs) maxs = csize[c]; }
        printf("S %ld ntriples %ld pairs %ld classes %d maxsize %d\n", S, ntr, pairs, ncls, maxs);
        for (int c = 0; c < ncls; c++) {
            if (csize[c] < 3 && S > LISTMAX) continue;
            uint64_t P2 = crp[c], Q2 = crq[c], R2 = (uint64_t)S - P2 - Q2;
            uint64_t f2 = P2 * Q2 + Q2 * R2 + R2 * P2, f3 = P2 * Q2 * R2, gg = g64(f2, f3);
            uint64_t G = g64(g64(P2, Q2), R2);
            for (int m = chead[c]; m >= 0; m = cnext[m]) G = g64(G, g64(g64(cp[m], cq[m]), (uint64_t)S - cp[m] - cq[m]));
            printf("C %ld %d %d %llu %llu %u,%u,%llu", S, csize[c], G == 1, (unsigned long long)(f2 / gg),
                   (unsigned long long)(f3 / gg), crp[c], crq[c], (unsigned long long)R2);
            for (int m = chead[c]; m >= 0; m = cnext[m]) printf(" %u,%u,%lu", cp[m], cq[m], (unsigned long)(S - cp[m] - cq[m]));
            printf("\n");
        }
        fflush(stdout);
    }
    return 0;
}
