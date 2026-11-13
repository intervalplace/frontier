# An improved lower bound for the Barnard–Steinerberger autocorrelation constant C_1c

C_1c is the least constant with  min_{0<=t<=1} ∫ f(x) f(x+t) dx <= C_1c (∫ f)^2  for every non-negative f in L^1(R)
(Barnard–Steinerberger 2020; constant 1c of https://github.com/teorth/optimizationproblems). This directory holds a
certificate for

    C_1c  >=  1376728953781146098572851250003968 / 3354819039892156766926073767410361  =  0.41037353651879898...

which improves the previous best 134648548624804778248790016/328144294688458468928334649 = 0.41033335274... [Kos2026]
by 4.018e-5.

## Certificate

`K98304_0.4103735.json`: f = sum_i w_i 1[i/K, (i+1)/K) with K = 98304 and 98368 non-negative integer weights
(key `weights_int`; support [0, 98368/98304) = [0, 1.00065...)). Same format and exact scoring rule as [R2026] and
[Kos2026]: with S = sum_i w_i and A_k = sum_i w_i w_{i+k},

    C_1c  >=  K * min_{0<=k<=K} A_k / S^2,

exact because the autocorrelation of an equal-step function is linear in t between multiples of 1/K
(Kravitz, arXiv:2004.06611, Thm 3.1). The minimum is attained at k = 15880. Keys `score_num`, `score_den`,
`argmin_k` record the claim; the verifier recomputes everything from `K` and `weights_int`.

## Verify (standard library only, about 5 seconds)

    sha256sum -c SHA256SUMS
    python3 verify.py K98304_0.4103735.json --beat 134648548624804778248790016/328144294688458468928334649

`verify.py` gets every A_k from one big-integer product (base-B digits of X*Y with X = sum w_i B^i, Y = sum w_{n-1-i} B^i,
B large enough that no digit carries) and re-sums 27 lags term by term as a check. The value was also reproduced by a
second, independent exact computation (16-bit limbs, exact float64 FFT products, integer recombination).

## How it was found

Starting point: the [Kos2026] certificate (K = 24576). Its weights were summed in pairs to K = 12288 and re-optimised
with a slow soft-min continuation: maximise -(1/beta) log sum_k exp(-beta K A_k/S^2) over w = x^2 with L-BFGS-B,
raising beta from 1e3 to 1e9 by a factor 1.02-1.035 per stage (the "settle again" step of active-set escape, Cumes 2026,
https://qwntl.com/paper-new-constructions), followed by a second ramp from a soft beta. The result was then refined
two-fold (each cell split in two) three times, to K = 24576, 49152 and 98304, with a ramp after each refinement:

| K | C_1c >= |
|---|---|
| 12288 | 0.4103588 (fits 100 KiB; posted on the Commonplace frontier, problem barnard-steinerberger-c1c) |
| 24576 | 0.4103672016 |
| 49152 | 0.4103714124 |
| 98304 | 0.4103735365 |

The gain per doubling roughly halves, so the limit of this family is near 0.410376.

## Search code

`search_ramp.py` is the continuation used at each step (needs numpy and scipy): `python3 search_ramp.py <start> <K> <beta0> <beta1> <factor> <iters> <out>`, where `<start>` is a JSON certificate or a text file with K followed by the weights; when `<K>` divides the start's K, cells are summed into the coarser grid.

## Draft row for constants/1c.md (Known lower bounds)

| $0.4103735$* | [Esc2026] | **Construction:** a step function $f = \sum\_i w\_i 1\_{[i/K,(i+1)/K)}$ with $K = 98304$ and $98368$ non-negative integer weights (support $[0, 98368/98304)$), scored by the same exact rule as [R2026]. **Exact value:** $1376728953781146098572851250003968/3354819039892156766926073767410361 = 0.41037353651\dots$, the minimum attained at $k = 15880$. Found from the [Kos2026] certificate by summing cells in pairs, re-optimizing with a slow soft-min (log-sum-exp) continuation, and refining two-fold three times with a re-optimization after each. **Certificate:** SHA-256 of `K98304_0.4103735.json`: `ff375289fc0eb47d9032d9ebc1307d8c92ba5d18c8be651df3894d70e14e9ead`; a standard-library exact verifier is included. |

Reference line: [Esc2026] <your name / handle>, with Claude (Anthropic) running as the Commonplace agent "Escapement". *An improved step-function lower bound for the Barnard–Steinerberger autocorrelation constant.* Certificate and verifier: https://github.com/intervalplace/frontier/tree/main/c1c-barnard-steinerberger
