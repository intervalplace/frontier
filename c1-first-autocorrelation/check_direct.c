/* Independent exact check for a C1 (1a) certificate with nonnegative integer heights a_i <= 2^50: every
   c_k = sum_{i+j=k} a_i a_j is summed directly in 128-bit integers (O(n^2)); prints max c_k, its index and S = sum a_i.
   (c_k <= n * 2^100 < 2^119 and S < 2^69, so nothing overflows.)  The bound is 2n max c_k / S^2. */
#include <stdio.h>
#include <stdlib.h>
typedef __int128 i128;
static void p128(i128 x) { char b[64]; int i = 63; b[i] = 0; if (x == 0) { printf("0"); return; } while (x > 0) { b[--i] = '0' + (int)(x % 10); x /= 10; } printf("%s", b + i); }
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); size_t cap = 1 << 20, n = 0; long long *a = malloc(cap * sizeof(long long)), x;
    while (fscanf(f, "%lld", &x) == 1) { if (x < 0 || x > (1LL << 50)) { printf("height out of range\n"); return 1; } a[n++] = x; }
    i128 S = 0; for (size_t i = 0; i < n; i++) S += a[i];
    i128 best = -1; size_t arg = 0;
    for (size_t k = 0; k <= 2 * n - 2; k++) {
        size_t lo = k >= n ? k - n + 1 : 0, hi = k < n ? k : n - 1; i128 c = 0;
        for (size_t i = lo, j = k - lo; i < j; i++, j--) c += 2 * (i128)a[i] * a[j];
        if ((k % 2) == 0 && k / 2 >= lo && k / 2 <= hi) c += (i128)a[k / 2] * a[k / 2];
        if (c > best) { best = c; arg = k; }
    }
    printf("%zu ", n); p128(best); printf(" %zu ", arg); p128(S); printf("\n");
    return 0;
}
