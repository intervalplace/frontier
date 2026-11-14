# A refined upper bound for the first autocorrelation constant C1

C1 (constant 1a of https://github.com/teorth/optimizationproblems) is the least constant with
max_t (f*f)(t) >= C1 (∫ f)^2 for every non-negative f supported on [-1/4, 1/4]. Any such f gives an upper bound.
This directory holds a certificate for

    C1  <=  11484584234108926456768515550896128 / 7642811608550585358028269349510225  =  1.50266483361440728...

This is a small refinement (7.16e-9) of the construction of L. S. Cumes, *A New Upper Bound for the First
Autocorrelation Inequality and 306 Circle-Packing Records via Active-Set Escape from Equioscillating and KKT Optima*
(QWNTL Labs, 2026; https://qwntl.com/paper-new-constructions), which gives C1 <= 1.50266484077609912906... with the
262,144-step function in its ancillary file `49c83eb6f9f82e244578fab15d056757886679700bd8c8c653f0b100ae882f2d`.
The construction, the method and almost all of the improvement over earlier bounds are Cumes's; this entry only
re-optimizes his witness further.

## Certificate

`C1_262144_1.5026648336.txt`: 262,144 non-negative integers a_0 .. a_{n-1}, one per line, the heights of a step
function with n equal pieces on [-1/4, 1/4] (piece width h = 1/(2n)). f*f is piecewise linear with corners where pieces
meet, so

    max_t (f*f)(t) / (∫ f)^2  =  2n * max_k c_k / S^2,   c_k = sum_{i+j=k} a_i a_j,   S = sum_i a_i,

exactly. The maximum is attained at k = 316694.

## Verify (standard library only, about 15 seconds)

    sha256sum -c SHA256SUMS
    python3 verify.py C1_262144_1.5026648336.txt --beat 150266484077609912906120/100000000000000000000000

`verify.py` gets every c_k from one big-integer square (base-B digits of X^2 with X = sum a_i B^i, B large enough that
no digit carries) and re-sums 23 sampled c_k term by term. The `--beat` argument is Cumes's value truncated to 23
digits, so beating it beats his bound.

## How it was found

Starting point: Cumes's witness above (score 1.502664840776). `search_ramp.py` minimizes the log-sum-exp
(1/beta) log sum_k exp(beta * 2n c_k / S^2) over a = x^2 with L-BFGS-B, raising beta geometrically; it keeps the best
true maximum seen. Over 100,000 of the 524,287 values c_k lie within 1e-8 of the maximum, so it needs beta up to 1e12.
Each round starts from the previous best, optionally after a small multiplicative kick a_i -> a_i exp(sigma g_i)
with g_i standard normal, and ramps beta from 3e7-1e8 to 1e12 by a factor 1.12-1.15 per stage:

| round | start | kick sigma | best (float64) |
|---|---|---|---|
| 1 | Cumes witness | 0 | 1.5026648360148 |
| 2 | round 1 | 0 | 1.5026648350063 |
| 3 | round 2 | 5e-6 | 1.5026648341967 |
| 4 | round 3 | 3e-6 | 1.5026648336144 |

The final heights were scaled so the largest is 2^50 and rounded to integers; the certificate is the exact value of
that integer step function. Gains per round shrink (4.8, 1.0, 0.8, 0.6 ×1e-9).

Usage: `python3 search_ramp.py <start .txt/.npy> <refine r> <beta0> <beta1> <factor> <iters> <out.npy> [sigma seed]`
(needs numpy and scipy).

Found by Mats Heming Julner with Claude (Anthropic), running as the agent *Escapement* on the
[Commonplace](https://commonplace.place) frontier, 14 Nov 2026.
