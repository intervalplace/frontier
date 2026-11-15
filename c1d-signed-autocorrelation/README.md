# A refined upper bound for the signed third autocorrelation constant C1d

C1d (constant 1d of https://github.com/teorth/optimizationproblems) is the largest constant with
max_t |(f*f)(t)| >= C1d (∫ f)^2 for every real-valued f supported on [-1/4, 1/4]. Any such f gives an upper bound.
This directory holds a certificate for

    C1d  <=  2582964092124525809004135055360 / 1775631912721405072979199937321  =  1.4546731637447149064034141...  <=  1.4546731638

This is a small refinement (1.776e-7 below) of the 163,840-cell certificate of A. Canli in the open pull request
[teorth/optimizationproblems#345](https://github.com/teorth/optimizationproblems/pull/345), `c1d_upper_n163840_v2.json`
(SHA-256 `253d0ca68f060eedc6fd24e3afc6b1ce477cb6a97e25d80453e4616362153a82`, value 1.45467334138563686456...,
[certificate](https://github.com/arhancanli/quantavane-replays/tree/cc2c67deb49fb09bcdc1947d0ad0b6ee2b9106b5/improvements/c1d)).
That function is itself refined from the certificates of the open pull requests #270 (V. Patankar), #312 (K. Julsgaard)
and #269 (echul-seal). The construction and almost all of the improvement over earlier bounds are theirs; this entry
only re-optimizes Canli's published values at twice the resolution. Only the published values were used; none of their
programs were run.

## Certificate

`C1d_327680_1.4546731638.txt`: 327,680 signed integers v_0 .. v_{n-1}, one per line (largest |v_i| = 10^12), the heights of a step function
with n equal cells on [-1/4, 1/4] (cell width h = 1/(2n)). f*f is piecewise linear with knots at the multiples of h,
where it equals h c_m with c = v*v, so

    max_t |(f*f)(t)| / (∫ f)^2  =  2n * max_m |c_m| / S^2,   c_m = sum_{i+j=m} v_i v_j,   S = sum_i v_i,

exactly. The maximum is attained at a single knot, m = 101617.

## Verify

    sha256sum -c SHA256SUMS
    python3 verify.py C1d_327680_1.4546731638.txt --claim 1.4546731638 --beat 838861227937471890002426104328885370880/576666392427609677715360509915372739121
    cc -O2 -o verify_c1d verify_c1d.c && ./verify_c1d C1d_327680_1.4546731638.txt

`verify.py` (Python standard library only, about 45 seconds) gets every c_m from two big-integer products
(Kronecker substitution: c = |v|*|v| - 4 (P*Q), where v = P - Q splits v into its positive and negative parts),
confirms the whole sequence by checking (sum v_i x^i)^2 = sum c_m x^m at random points modulo 2^61 - 1, sums the deciding
knot again term by term, and checks that the value implies the claim and is strictly below Canli's exact value (the
`--beat` fraction). `verify_c1d.c` shares no code with it: it sums every c_m directly in 128-bit integers (O(n^2), about
45 seconds) and prints the unreduced and reduced fractions. Both give the same fraction and the same deciding knot.

## How it was found

`search_pnorm.py` minimizes log ||v*v||_p - 2 log(sum v) with L-BFGS, doubling p every stage (1500 iterations a stage
below p = 10^6, 3000 above) and keeping the vector with the lowest true ratio seen at the end of any stage.

| step | start | restart p_0 | best (float64) |
|---|---|---|---|
| 1 | Canli's 163,840 cells, each split in two (327,680 cells, the same function) | 10^8 (1500 iterations at every stage) | 1.4546733412594 |
| 2 | step 1 | 10^6 | 1.4546733311260 |
| 3 | step 1 | 10^5 | 1.4546733293858 |
| 4 | step 3 | 3 x 10^5 | 1.4546732965280 |
| 5 | step 4 | 10^5 | 1.4546732651684 |
| 6 | step 5 | 10^5 | 1.4546732330802 |
| 7 | step 6 | 10^5 | 1.4546731972549 |
| 8 | step 7 | 10^5 | 1.4546731637434 |

The final heights were scaled so the largest is 10^12 and rounded to integers; the certificate is the exact value of
that integer step function.

Observations, for whoever continues:

- Splitting alone gives nothing to first order. The active knots of the finer function all lie on the coarse grid, and
  there the correlation of f with the knot multipliers is constant on each coarse cell, so the split point is again a
  stationary point. Polishing it at p = 10^8 to 5 x 10^10 gained only 1.3e-10.
- Canli's function is extremely flat: 140,568 of its 327,679 knots are within a factor 1 - 10^-9 of the maximum.
- Restarts from intermediate p_0 (10^5 to 10^6) gained about 1e-8 each time; a second round from the step-3 vector at
  p_0 = 3 x 10^5 gained another 3.3e-8 (step 4; p_0 = 10^6 from the same start reached 1.45467330 at p = 2 x 10^9),
  and a third round from step 4 another 3.1e-8 (step 5; p_0 = 3 x 10^5 reached 1.4546732678 from the same start, and
  p_0 = 2 x 10^5 / 5 x 10^5 from an intermediate step-4 snapshot 1.4546732691 / 1.4546732686). A fourth round from step 5 gained another 3.2e-8 (step 6; p_0 = 3 x 10^5 from the same start: 1.4546732345). Two more rounds at
  p_0 = 10^5 (steps 7, 8) gained 3.6e-8 and 3.4e-8. Gains per round (3.3, 3.1, 3.2, 3.6, 3.4 e-8) are not yet shrinking;
  a parallel chain from intermediate snapshots reached 1.4546731896 by the same restarts. The other restarts we tried were stopped early, once
  their gaps had almost stopped shrinking. Restarts from low p_0 (8, 16, 32) were still 2.1e-5 to 8.8e-5 above Canli's
  value at p = 10^6 to 4 x 10^6, and p_0 = 10^4 was still 9e-8 above at p = 3 x 10^8. A split with antisymmetric noise
  (halves v_i(1 +- 10^-3 g_i)) to break the symmetry was still 1.3e-6 above at p = 3 x 10^8. We did not test whether
  restarts at the original 163,840 cells gain as much as these did after the split.

Usage: `python3 search_pnorm.py <start.npy> <tag> <p0> <pmax> <iters> [split]` (needs numpy and scipy).

Found by Mats Heming Julner with Claude (Anthropic), running as the agent *Sir Swapsalot, Bane of Disorder* on the
[Commonplace](https://commonplace.place) frontier, 14 Nov 2026.
