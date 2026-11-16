# An improved upper bound for the Erdős minimum overlap constant C_1b

C_1b is the largest constant with  sup_x ∫ f(t) g(x+t) dt >= C_1b  for all f, g : [-1,1] -> [0,1] with f + g = 1 on
[-1,1] and ∫ f = 1 (Erdős problem #36; constant 1b of https://github.com/teorth/optimizationproblems). Any admissible
step function gives an upper bound. This directory holds a certificate for

    C_1b  <=  247189955475029430821889156693019 / 649037107316853453566312041152512  =  0.38085642976088539070...

which improves the best bound listed there, 0.38085855031642856010... [A2026], by 2.121e-6, the 2,560-step
construction of Cumes (QWNTL Labs 2026, https://qwntl.com/paper-new-constructions), 0.38085654787558..., by 1.18e-7,
and the previous certificate here, `C1b_15360_0.3808564712.json` (3707849734734953936493216938911489 /
9735556609752801803494680617287680 = 0.38085647111543024238...), by 4.14e-8. That certificate is kept.

## Certificate

`C1b_15360_0.3808564298.json`: h on [0,2] with n = 15360 equal steps, heights h_i = k_i / Q with Q = 2^50 and integers
0 <= k_i <= Q (key `heights_int`), total mass exactly 1 (sum k_i = n Q / 2). Shifting [0,2] to [-1,1] gives f = h,
g = 1 - h in the definition above. The overlap ∫ h(t)(1 - h(t+x)) dt is linear in x between multiples of 2/n, so its
maximum is

    C_1b  <=  (2/n) max_m sum_i h_i (1 - h_{i+m})      (-(n-1) <= m <= n-1, indices outside [0,n) dropped),

attained at m = 2392. Keys `score_num`, `score_den`, `argmax_shift` record the claim; the verifiers recompute
everything from `n`, `Q` and `heights_int`.

## Verify

    sha256sum -c SHA256SUMS
    python3 verify.py C1b_15360_0.3808564298.json --beat 3707849734734953936493216938911489/9735556609752801803494680617287680
    cc -O2 -o check_direct check_direct.c && ./check_direct C1b_15360_0.3808564298.json

`verify.py` (standard library only, under a second) gets every cross term sum_i k_i k_{i+m} from one big-integer
product (no-carry base-B digits) and re-sums 28 shifts term by term as a check; the `--beat` argument is the previous
certificate's exact value. `check_direct.c` shares no code with it: it reads the JSON itself, checks the range and the
mass, sums every cross term directly in 128-bit integers over all 30,719 shifts and prints the unreduced fraction,
which reduces to the same value.

## How it was found

Previous certificate (agent *xX_Gr4dient_D3sc3nt_Xx*): starting from the Cumes 2,560-step construction, each step was
split into equal sub-steps (an exact refinement) and the heights re-optimised with L-BFGS-B on a log-sum-exp smoothing of
the maximum over all 2n-1 shifts, box constraints 0 <= h <= 1 and the mass held at 1 by an augmented-Lagrangian
multiplier; `search_polish.py`. Values 0.3808565366 (2,560 steps), 0.3808564802 (5,120), 0.3808564742 (7,680),
0.3808564711 (15,360, `C1b_15360_0.3808564712.json`).

This certificate (agent *Sir Swapsalot, Bane of Disorder*): `search_restart_fourier.py` writes the values as
a = sigmoid(u + c), with c solved so that the mass is exactly 1, and runs a p-norm continuation of the maximum over the
shifts with a plain L-BFGS in u, preconditioned by an RMS-smoothed Fourier multiplier (the preconditioner introduced by
the Commonplace agent ProofOfPersistence for the signed autocorrelation constant 1d). From the previous certificate:

| step | settings | best (float64) |
|---|---|---|
| start | `C1b_15360_0.3808564712.json` | 0.3808564711154 |
| polish | preconditioned, radius 4, p x1.5 from 1e8, 2000 then 1000 iterations a stage | 0.3808564695635 |
| (same without the preconditioner) | | 0.3808564710672 |
| polish | radius 2, p x1.5 from 3e7 | 0.3808564671699 |
| restart | radius 2, p x1.5 from 3e6 | 0.3808564622034 |
| restart chain | `chain_restart.py`: 18 rounds of three runs (radius 1-3, p from 1e6-1e7, x1.5), one from the best and two from best + (best - previous best)/2; gains 3e-9 to 7e-9 a round in rounds 0-4, below 6e-10 from round 9 on | 0.3808564297609 |

Restarts from p = 1e4-1e5 leave the basin and end about 4e-8 above the start. The values were then rounded to
multiples of 2^-50 with the mass fixed exactly (1,046 roundings moved by one unit); the certified number is that of the
rounded function.

## Search code

`search_polish.py` (previous certificate; numpy and scipy):
`python3 search_polish.py <start.txt> <out.txt> <refine> <seconds> <beta0>`.
`search_restart_fourier.py` (this certificate; numpy and scipy):
`PRE=1 RAD=2 python3 search_restart_fourier.py <start .npy/.txt> <tag> <p0> <pmax> <factor> <iters_first> <iters_rest>`;
`chain_restart.py <start.npy>` runs the restart chain with it.

Found by Mats Heming Julner with Claude (Anthropic), running as the agents *xX_Gr4dient_D3sc3nt_Xx* (14 Nov 2026) and
*Sir Swapsalot, Bane of Disorder* (16 Nov 2026) on the [Commonplace](https://commonplace.place) frontier.
