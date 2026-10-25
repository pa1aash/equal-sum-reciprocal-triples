/*
 * Exact kernel for the two-coefficient degeneracy enumeration.
 *
 * For each cone-order sum S in [lo, hi] it lists every hyperbolic triple
 * 2 <= p <= q <= r, p + q + r = S, keys it by the reduced fraction
 * R = e2/e3 = (pq + qr + rp)/(pqr), sorts by key, and groups equal keys.
 * Integer arithmetic only. Every bound that keeps the arithmetic exact is
 * asserted at run time rather than assumed.
 *
 * Output, one line per record, tab separated:
 *   S  S  triads  distinct_R  pairs  classes  max_fibre  prim_classes  prim_pairs
 *   G  S  num  den  k  content  p,q,r;p,q,r;...
 * where prim_pairs counts unordered pairs whose six entries have gcd 1.
 *
 * Mode "dump": a single S, one line "num den p q r" per hyperbolic triple.
 *
 * Usage: enumerate_core range LO HI
 *        enumerate_core dump S
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    uint64_t key;   /* (num << 32) | den, both asserted < 2^32 */
    uint32_t p, q;
} Item;

static uint64_t gcd64(uint64_t a, uint64_t b) {
    while (b) { uint64_t t = a % b; a = b; b = t; }
    return a;
}

static void die(const char *msg, long S) {
    fprintf(stderr, "enumerate_core: %s (S = %ld)\n", msg, S);
    exit(1);
}

/* LSD radix sort on the 64-bit key, 8 bits per pass. */
static void radix_sort(Item *a, Item *tmp, size_t n) {
    for (int shift = 0; shift < 64; shift += 8) {
        size_t count[257];
        memset(count, 0, sizeof count);
        for (size_t i = 0; i < n; i++) count[((a[i].key >> shift) & 0xff) + 1]++;
        int trivial = 0;
        for (int b = 0; b < 256; b++) if (count[b + 1] == n) trivial = 1;
        if (trivial) continue;
        for (int b = 0; b < 256; b++) count[b + 1] += count[b];
        for (size_t i = 0; i < n; i++) tmp[count[(a[i].key >> shift) & 0xff]++] = a[i];
        memcpy(a, tmp, n * sizeof(Item));
    }
}

static size_t fill(long S, Item *buf) {
    size_t n = 0;
    for (uint64_t p = 2; 3 * p <= (uint64_t)S; p++) {
        for (uint64_t q = p; 2 * q <= (uint64_t)S - p; q++) {
            uint64_t r = (uint64_t)S - p - q;
            uint64_t e2 = p * q + q * r + r * p;
            uint64_t e3 = p * q * r;
            if (e2 >= e3) continue;               /* not hyperbolic */
            uint64_t g = gcd64(e2, e3);
            uint64_t num = e2 / g, den = e3 / g;
            if (num >> 32 || den >> 32) die("key overflow", S);
            buf[n].key = (num << 32) | den;
            buf[n].p = (uint32_t)p;
            buf[n].q = (uint32_t)q;
            n++;
        }
    }
    return n;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: enumerate_core range LO HI | dump S\n"); return 2; }
    int dump = strcmp(argv[1], "dump") == 0;
    long lo = atol(argv[2]);
    long hi = dump ? lo : atol(argv[3]);
    if (lo < 6 || hi < lo || hi > 4800) die("S out of supported range", hi);

    size_t cap = (size_t)hi * (size_t)hi / 12 + (size_t)hi + 16;
    Item *buf = malloc(cap * sizeof(Item));
    Item *tmp = malloc(cap * sizeof(Item));
    if (!buf || !tmp) die("out of memory", hi);

    for (long S = lo; S <= hi; S++) {
        size_t n = fill(S, buf);
        if (n > cap) die("buffer overrun", S);
        if (dump) {
            for (size_t i = 0; i < n; i++) {
                uint64_t p = buf[i].p, q = buf[i].q, r = (uint64_t)S - p - q;
                printf("%llu %llu %llu %llu %llu\n",
                       (unsigned long long)(buf[i].key >> 32),
                       (unsigned long long)(buf[i].key & 0xffffffffULL),
                       (unsigned long long)p, (unsigned long long)q, (unsigned long long)r);
            }
            continue;
        }
        radix_sort(buf, tmp, n);

        uint64_t distinct = 0, pairs = 0, classes = 0, maxk = n ? 1 : 0;
        uint64_t prim_classes = 0, prim_pairs = 0;
        size_t i = 0;
        while (i < n) {
            size_t j = i + 1;
            while (j < n && buf[j].key == buf[i].key) j++;
            distinct++;
            uint64_t k = j - i;
            if (k > maxk) maxk = k;
            if (k >= 2) {
                classes++;
                pairs += k * (k - 1) / 2;
                uint64_t content = 0;
                for (size_t a = i; a < j; a++) {
                    content = gcd64(content, buf[a].p);
                    content = gcd64(content, buf[a].q);
                }
                content = gcd64(content, (uint64_t)S);
                if (content == 1) prim_classes++;
                for (size_t a = i; a < j; a++)
                    for (size_t b = a + 1; b < j; b++) {
                        uint64_t g = gcd64(gcd64(buf[a].p, buf[a].q), gcd64(buf[b].p, buf[b].q));
                        if (gcd64(g, (uint64_t)S) == 1) prim_pairs++;
                    }
                /* triples in lexicographic order: sort the small slice by (p, q) */
                for (size_t a = i + 1; a < j; a++) {
                    Item x = buf[a]; size_t b = a;
                    while (b > i && (buf[b - 1].p > x.p || (buf[b - 1].p == x.p && buf[b - 1].q > x.q))) {
                        buf[b] = buf[b - 1]; b--;
                    }
                    buf[b] = x;
                }
                printf("G\t%ld\t%llu\t%llu\t%llu\t%llu\t", S,
                       (unsigned long long)(buf[i].key >> 32),
                       (unsigned long long)(buf[i].key & 0xffffffffULL),
                       (unsigned long long)k, (unsigned long long)content);
                for (size_t a = i; a < j; a++)
                    printf("%s%u,%u,%llu", a == i ? "" : ";", buf[a].p, buf[a].q,
                           (unsigned long long)((uint64_t)S - buf[a].p - buf[a].q));
                printf("\n");
            }
            i = j;
        }
        printf("S\t%ld\t%zu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n", S, n,
               (unsigned long long)distinct, (unsigned long long)pairs,
               (unsigned long long)classes, (unsigned long long)maxk,
               (unsigned long long)prim_classes, (unsigned long long)prim_pairs);
    }
    free(buf); free(tmp);
    return 0;
}
