# An improved upper bound for the signed third autocorrelation constant C_1d

C_1d is the largest constant with  max_t |(f*f)(t)| >= C_1d (∫ f)^2  for every real f supported on [-1/4, 1/4]
(constant 1d of https://github.com/teorth/optimizationproblems). Any f with ∫ f != 0 gives an upper bound. This
directory holds a certificate for

    C_1d  <=  279676395642374208108178964480 / 192135142203433701871000999249  =  1.45562333071922556161...

which improves the best certificate we know of, 279676466330630562685427712000/192135175720687825901095191121
= 1.45562344470022455281... (Kerem Özdemir with Claude, 26 Oct 2026, https://github.com/Keremozdemirra/tao-1c-1d-certificates,
itself refined from the 10240-cell certificate of echul-seal, https://github.com/echul-seal/c1d-upper-bound), by 1.140e-7.

## Certificate

`C1d_40960_1.4556233308.json`: f has n = 40960 equal cells on [-1/4, 1/4] with signed integer values v_i (key `v_int`,
max |v_i| = 10^12). f*f is piecewise linear with knots at multiples of 1/(2n), so

    max_t |f*f| / (∫ f)^2  =  2n * max_m |c_m| / (sum v)^2,   c = v*v (full discrete autoconvolution),

attained at m = 15313. Keys `score_num`, `score_den`, `argmax_m` record the claim.

## Verify

    sha256sum -c SHA256SUMS
    python3 verify.py C1d_40960_1.4556233308.json --beat 279676466330630562685427712000/192135175720687825901095191121
    cc -O2 -o verify_direct verify_direct.c && python3 -c "import json;print('\n'.join(map(str,json.load(open('C1d_40960_1.4556233308.json'))['v_int'])))" | ./verify_direct

`verify.py` (standard library, ~3 s) gets every c_m from three no-carry big-integer products (v = p - q with p, q >= 0)
and re-sums 30 of them term by term. `verify_direct.c` is a plain O(n^2) loop over all 81919 values of m in 128-bit
integers, sharing no code with `verify.py`; it prints the same numerator and denominator.

## How it was found

Starting point: the Özdemir 40960-cell certificate (1.4556234447). Polished with L-BFGS on a log-sum-exp smoothing of
max_m |c_m| (soft-max over +c and -c), beta raised from 1e8 by a factor 1.6 per stage to 2e10, ftol 1e-16, then scaled
to max |v| = 10^12 and rounded to integers; the certified number is that of the rounded function. Small-kick restarts
(0.2-0.5% multiplicative noise, re-polish from beta 1e9) find only ~1e-10 more, so this basin is nearly converged.

## Search code

`search_polish.py` (numpy, scipy): `python3 search_polish.py <start.txt> <out.txt> <seconds> <beta0> [noise seed]`.

## Draft row for constants/1d.md (Known upper bounds)

| $1.4556233308$* | [Grad2026] | **Construction:** a step function with $n = 40960$ equal cells on $[-1/4,1/4]$ and signed integer values ($\max \lvert v\_i \rvert = 10^{12}$). **Exact value:** $279676395642374208108178964480/192135142203433701871000999249 = 1.45562333071\ldots$, the maximum of $\lvert v * v \rvert$ at index $m = 15313$; $1.14 \cdot 10^{-7}$ below the Özdemir certificate it was polished from. **Certificate:** SHA-256 of `C1d_40960_1.4556233308.json` in SHA256SUMS; a standard-library exact verifier and an independent $O(n^2)$ C check are included. |
