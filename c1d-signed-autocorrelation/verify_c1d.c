/* Exact check of a C1d certificate, independent of verify.py: every knot value
   c_m = sum_i v_i v_{m-i} is summed directly in 128-bit integers (O(n^2)).
   Prints the unreduced fraction 2n max|c_m| / S^2, the reduced fraction and the deciding knot.
   Refuses inputs where a 128-bit overflow would be possible.   cc -O2 -o verify_c1d verify_c1d.c */
#include <stdio.h>
#include <stdlib.h>
typedef __int128 i128;
static void pr(i128 x){ char b[64]; int k=0; if(x<0){putchar('-'); x=-x;} do{ b[k++]='0'+(int)(x%10); x/=10; }while(x); while(k) putchar(b[--k]); }
static i128 gcd(i128 a, i128 b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ i128 t=a%b; a=b; b=t; } return a; }
int main(int argc, char **argv){
  FILE *f = fopen(argv[1], "r"); if(!f){ perror("open"); return 2; }
  size_t cap = 1<<16, n = 0; long long *v = malloc(cap*sizeof *v), x;
  while(fscanf(f, "%lld", &x) == 1){ if(n == cap){ cap *= 2; v = realloc(v, cap*sizeof *v); } v[n++] = x; }
  fclose(f);
  long long mx = 0; i128 S = 0;
  for(size_t i=0;i<n;i++){ long long a = v[i] < 0 ? -v[i] : v[i]; if(a > mx) mx = a; S += v[i]; }
  /* |c_m| <= n mx^2 and S^2 <= (n mx)^2: require both, and 2n max|c|, below 2^125 */
  long double bound = (long double)n * mx * mx * 2 * n;
  if(bound >= 4.0e37L || (long double)n*mx*(long double)n*mx >= 4.0e37L){ fprintf(stderr, "refused: possible overflow\n"); return 2; }
  i128 best = -1; long long arg = -1, nmax = 0;
  for(long long m = 0; m < 2*(long long)n - 1; m++){
    long long lo = m - (long long)n + 1 > 0 ? m - (long long)n + 1 : 0, hi = m < (long long)n - 1 ? m : (long long)n - 1;
    i128 c = 0; long long i = lo, j = hi;
    for(; i < j; i++, j--) c += (i128)v[i] * v[j];
    c *= 2; if(i == j) c += (i128)v[i] * v[i];
    if(c < 0) c = -c;
    if(c > best){ best = c; arg = m; nmax = 1; } else if(c == best) nmax++;
  }
  i128 num = (i128)2*n*best, den = S*S, g = gcd(num, den);
  printf("n = %zu, max|v| = %lld, deciding knot m = %lld (%lld knot(s) at the maximum)\n", n, mx, arg, nmax);
  printf("unreduced: "); pr(num); printf(" / "); pr(den); printf("\n");
  printf("reduced:   "); pr(num/g); printf(" / "); pr(den/g); printf("\n");
  printf("value ~ %.20Lf\n", (long double)num/(long double)den);
  return 0;
}
