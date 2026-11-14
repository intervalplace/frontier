# An improved upper bound for the Erdős minimum overlap constant C_1b

C_1b is the largest constant with  sup_x ∫ f(t) g(x+t) dt >= C_1b  for all f, g : [-1,1] -> [0,1] with f + g = 1 on
[-1,1] and ∫ f = 1 (Erdős problem #36; constant 1b of https://github.com/teorth/optimizationproblems). Any admissible
step function gives an upper bound. This directory holds a certificate for

    C_1b  <=  3707849734734953936493216938911489 / 9735556609752801803494680617287680  =  0.38085647111543024238...

which improves the best bound listed there, 0.38085855031642856010... [A2026], by 2.079e-6, and the 2,560-step
construction of Cumes (QWNTL Labs 2026, https://qwntl.com/paper-new-constructions), 0.38085654787558..., by 7.68e-8.

## Certificate

`C1b_15360_0.3808564712.json`: h on [0,2] with n = 15360 equal steps, heights h_i = k_i / Q with Q = 2^50 and integers
0 <= k_i <= Q (key `heights_int`), total mass exactly 1 (sum k_i = n Q / 2). Shifting [0,2] to [-1,1] gives f = h,
g = 1 - h in the definition above. The overlap ∫ h(t)(1 - h(t+x)) dt is linear in x between multiples of 2/n, so its
maximum is

    C_1b  <=  (2/n) max_m sum_i h_i (1 - h_{i+m})      (-(n-1) <= m <= n-1, indices outside [0,n) dropped),

attained at m = 6384. Keys `score_num`, `score_den`, `argmax_shift` record the claim; the verifier recomputes
everything from `n`, `Q` and `heights_int`.

## Verify (standard library only, under a second)

    sha256sum -c SHA256SUMS
    python3 verify.py C1b_15360_0.3808564712.json --beat 38085855031642856010/100000000000000000000

`verify.py` gets every cross term sum_i k_i k_{i+m} from one big-integer product (no-carry base-B digits) and re-sums
28 shifts term by term as a check. The value was also reproduced by a second, independent exact computation: a plain
O(n^2) loop over all 30,719 shifts in Python integers, sharing no code with `verify.py`.

## How it was found

Starting point: the Cumes (QWNTL Labs) 2,560-step construction, 0.38085654787558... Each step was split into equal
sub-steps (an exact refinement, which leaves the bound unchanged) and the heights re-optimised with L-BFGS-B on a
log-sum-exp smoothing of the maximum over all 2n-1 shifts, with box constraints 0 <= h <= 1 and the mass held at 1 by
an augmented-Lagrangian multiplier (a plain quadratic penalty lets the optimiser trade mass for overlap, and
re-normalising then loses the gain). The smoothing parameter beta starts at 1e8-1e9 (lower values erase the active set
of the start) and rises by a factor 1.6 per stage to 3e10, with ftol 1e-16.

| steps | C_1b <= |
|---|---|
| 2,560 (native) | 0.3808565366 |
| 5,120 | 0.3808564802 |
| 7,680 | 0.3808564742 (fits 100 KiB; posted on the Commonplace frontier, problem erdos-minimum-overlap) |
| 10,240 | 0.3808564771 |
| 15,360 | 0.3808564711 (this certificate) |
| 23,040 | 0.3808564718 |

The heights were then rounded to multiples of 2^-50 and the mass fixed exactly; the certified number is that of the
rounded function. Gains are now below 1e-9 per refinement, so this family is close to its limit near 0.38085647.

## Search code

`search_polish.py` (numpy and scipy): `python3 search_polish.py <start.txt> <out.txt> <refine> <seconds> <beta0>`,
where `<start.txt>` lists the step heights and each step is split into `<refine>` equal sub-steps before polishing.

## Draft row for constants/1b.md (Known upper bounds)

| $0.3808564712$* | [Grad2026] | **Construction:** a step function $h$ on $[0,2]$ with $n = 15360$ equal steps, heights $k\_i/2^{50}$ with integers $0 \le k\_i \le 2^{50}$ and $\sum\_i k\_i = n \cdot 2^{49}$ (mass exactly $1$). **Exact value:** $3707849734734953936493216938911489/9735556609752801803494680617287680 = 0.38085647111\dots$, the maximum attained at shift $m = 6384$. Found from the 2,560-step construction of [Cum2026] by refining each step (an exact refinement) and re-optimizing with a log-sum-exp continuation. **Certificate:** SHA-256 of `C1b_15360_0.3808564712.json`: `7bfe24a36566ee3544d22f4064e9cb64495b21431a31f521de390dc5f0865b7b`; a standard-library exact verifier is included. |
