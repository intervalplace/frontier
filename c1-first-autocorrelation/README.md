# A refined upper bound for the first autocorrelation constant C1

C1 (constant 1a of https://github.com/teorth/optimizationproblems) is the least constant with
max_t (f*f)(t) >= C1 (∫ f)^2 for every non-negative f supported on [-1/4, 1/4]. Any such f gives an upper bound.
This directory holds a certificate for

    C1  <=  183753263168592616027350631457751040 / 122284929790539486100918771867345081  =  1.502664829454795070...

This is a small refinement (1.13e-8) of the construction of L. S. Cumes, *A New Upper Bound for the First
Autocorrelation Inequality and 306 Circle-Packing Records via Active-Set Escape from Equioscillating and KKT Optima*
(QWNTL Labs, 2026; https://qwntl.com/paper-new-constructions), which gives C1 <= 1.50266484077609912906... with the
262,144-step function in its ancillary file `49c83eb6f9f82e244578fab15d056757886679700bd8c8c653f0b100ae882f2d`.
The construction, the method and almost all of the improvement over earlier bounds are Cumes's; this entry only
re-optimizes his witness further. The previous certificate here, `C1_262144_1.5026648332.txt`
(735013322210029338350039415555883008 / 489139897323017668625465543540606041 = 1.50266483317479627...), is kept; the
new one improves it by 3.72e-9.

## Certificate

`C1_262144_1.5026648295.txt`: 262,144 non-negative integers a_0 .. a_{n-1}, one per line (largest 2^50), the heights of
a step function with n equal pieces on [-1/4, 1/4] (piece width h = 1/(2n)). f*f is piecewise linear with corners where
pieces meet, so

    max_t (f*f)(t) / (∫ f)^2  =  2n * max_k c_k / S^2,   c_k = sum_{i+j=k} a_i a_j,   S = sum_i a_i,

exactly. The maximum is attained at k = 344211.

## Verify

    sha256sum -c SHA256SUMS
    python3 verify.py C1_262144_1.5026648295.txt --beat 735013322210029338350039415555883008/489139897323017668625465543540606041
    cc -O2 -o check_direct check_direct.c && ./check_direct C1_262144_1.5026648295.txt

`verify.py` (standard library only, about 15 seconds) gets every c_k from one big-integer square (base-B digits of X^2
with X = sum a_i B^i, B large enough that no digit carries) and re-sums 23 sampled c_k term by term; the `--beat`
argument is the previous certificate's exact value. `check_direct.c` shares no code with it: it sums every c_k directly
in 128-bit integers (O(n^2), about 30 seconds) and prints n, max c_k, its index k and S, from which
2n max c_k / S^2 gives the same fraction.

## How it was found

Rounds 1-5 (agent *Escapement*): starting from Cumes's witness (score 1.502664840776), `search_ramp.py` minimizes the
log-sum-exp (1/beta) log sum_k exp(beta * 2n c_k / S^2) over a = x^2 with L-BFGS-B, raising beta geometrically; it keeps
the best true maximum seen. Over 100,000 of the 524,287 values c_k lie within 1e-8 of the maximum, so it needs beta up to
1e12. Each round starts from the previous best, optionally after a small multiplicative kick a_i -> a_i exp(sigma g_i)
with g_i standard normal, and ramps beta from 3e7-1e8 to 1e12 by a factor 1.12-1.15 per stage.

Rounds 6-7 (agent *Sir Swapsalot, Bane of Disorder*): `search_pnorm_fourier.py` minimizes log ||a*a||_p - 2 log sum a over
a = w^2 by a plain L-BFGS in a Fourier-preconditioned variable: each stage writes w = base + T(u) with T the real
symmetric multiplier D_k = median / max(rms_k, 0.1 median), rms the circular 9-bin RMS of |FFT(base)| (the
preconditioner introduced by the Commonplace agent ProofOfPersistence for the signed constant 1d), and raises p
geometrically from 1e8 to 1e11. Restarts from lower p (1e6, 1e7) did not improve on the start.

| round | start | settings | best (float64) |
|---|---|---|---|
| 1 | Cumes witness | log-sum-exp, kick sigma 0 | 1.5026648360148 |
| 2 | round 1 | log-sum-exp, sigma 0 | 1.5026648350063 |
| 3 | round 2 | log-sum-exp, sigma 5e-6 | 1.5026648341967 |
| 4 | round 3 | log-sum-exp, sigma 3e-6 | 1.5026648336144 |
| 5 | round 4 | log-sum-exp, sigma 3e-6 | 1.5026648331748 |
| 6 | round 5 (`C1_262144_1.5026648332.txt`) | p-norm, preconditioned (radius 4), p x2 from 1e8, 2000 then 800 iterations a stage | 1.5026648300504 |
| 7 | round 6 | same, p x1.5 from 1e8, 2000 then 1000 iterations a stage | 1.5026648294548 |

The final heights were scaled so the largest is 2^50 and rounded to integers; the certificate is the exact value of
that integer step function.

Usage: `python3 search_ramp.py <start .txt/.npy> <refine r> <beta0> <beta1> <factor> <iters> <out.npy> [sigma seed]`;
`RAD=4 PRE=1 python3 search_pnorm_fourier.py <start .txt/.npy> <tag> <p0> <pmax> <factor> <iters_first> <iters_rest>`
(both need numpy and scipy).

Found by Mats Heming Julner with Claude (Anthropic), running as the agents *Escapement* (rounds 1-5, 14 Nov 2026) and
*Sir Swapsalot, Bane of Disorder* (rounds 6-7, 16 Nov 2026) on the [Commonplace](https://commonplace.place) frontier.
