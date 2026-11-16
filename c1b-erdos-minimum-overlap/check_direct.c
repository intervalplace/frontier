/* Independent exact check for a C_1b certificate (shares no code with verify.py): reads n, Q and heights_int from the
   JSON file, checks 0 <= k_i <= Q and sum k_i = n Q / 2, and for every shift m sums X_m = sum_i k_i k_{i+m} directly in
   128-bit integers (O(n^2)); the overlap numerator is Q * S_m - X_m with S_m the sum of k_i over the overlapping range.
   Prints the shift with the largest overlap and the unreduced fraction 2 * best / (n Q^2).
   cc -O2 -o check_direct check_direct.c && ./check_direct C1b_15360_0.3808564298.json */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef __int128 i128;
static void p128(i128 x) { char b[64]; int i = 63; b[i] = 0; if (x == 0) { printf("0"); return; } while (x > 0) { b[--i] = '0' + (int)(x % 10); x /= 10; } printf("%s", b + i); }
static long long num_after(const char *s, const char *key) { const char *p = strstr(s, key); if (!p) { fprintf(stderr, "missing %s\n", key); exit(1); } return atoll(p + strlen(key)); }
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb"); fseek(f, 0, SEEK_END); long sz = ftell(f); rewind(f);
    char *s = malloc(sz + 1); if (fread(s, 1, sz, f) != (size_t)sz) return 1; s[sz] = 0;
    long n = num_after(s, "\"n\":"); long long Q = num_after(s, "\"Q\":");
    char *p = strstr(s, "\"heights_int\":["); if (!p) return 1; p += strlen("\"heights_int\":[");
    long long *k = malloc(n * sizeof(long long)); i128 tot = 0;
    for (long i = 0; i < n; i++) { k[i] = strtoll(p, &p, 10); if (k[i] < 0 || k[i] > Q) { printf("height %ld out of [0, Q]\n", i); return 1; } tot += k[i]; while (*p == ',' || *p == ' ') p++; }
    if (*p != ']') { printf("expected exactly n heights\n"); return 1; }
    if (2 * tot != (i128)n * Q) { printf("mass is not exactly 1\n"); return 1; }
    i128 *P = malloc((n + 1) * sizeof(i128)); P[0] = 0; for (long i = 0; i < n; i++) P[i + 1] = P[i] + k[i];
    i128 best = -1; long arg = 0;
    for (long m = -(n - 1); m <= n - 1; m++) {
        long lo = m < 0 ? -m : 0, hi = m > 0 ? n - 1 - m : n - 1;
        i128 X = 0; for (long i = lo; i <= hi; i++) X += (i128)k[i] * k[i + m];
        i128 o = (i128)Q * (P[hi + 1] - P[lo]) - X;
        if (o > best) { best = o; arg = m; }
    }
    printf("n = %ld, maximum at shift m = %ld\nC_1b <= ", n, arg); p128(2 * best); printf(" / "); p128((i128)n * Q * Q); printf("  (unreduced)\n");
    return 0;
}
