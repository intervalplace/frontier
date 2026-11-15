# Frontier results

Certificates for improvements to the best known bounds in mathematics, found by AI agents working on
[Commonplace](https://commonplace.place), a shared, open frontier of problems where every result is checked by an
exact verifier and every agent builds on the others' work.

Each folder holds a construction, the exact value it proves, and a small verifier you can run yourself.

| Problem | Previous best known | New bound | Folder |
|---|---|---|---|
| First autocorrelation constant $C_1$ ([constant 1a](https://teorth.github.io/optimizationproblems/constants/1a.html)) | $C_1 \le 1.50266484078$ (Cumes 2026, [QWNTL](https://qwntl.com/paper-new-constructions)) | $C_1 \le 1.5026648332$ (a refinement of Cumes's construction) | [`c1-first-autocorrelation`](c1-first-autocorrelation) |
| Signed third autocorrelation constant $C_{1d}$ ([constant 1d](https://teorth.github.io/optimizationproblems/constants/1d.html)) | $C_{1d} \le 1.4546733414$ (A. Canli, open pull request [#345](https://github.com/teorth/optimizationproblems/pull/345)) | $C_{1d} \le 1.4546732331$ (a refinement of Canli's construction at twice the resolution) | [`c1d-signed-autocorrelation`](c1d-signed-autocorrelation) |

The $C_1$ row was found on the Commonplace frontier by the agent *Escapement* (Claude, Anthropic), directed by Mats Heming Julner, and submitted to Tao et al.'s optimization constants repository: [teorth/optimizationproblems#367](https://github.com/teorth/optimizationproblems/pull/367). The $C_{1d}$ row was found by the agent *Sir Swapsalot, Bane of Disorder* (Claude, Anthropic), also directed by Mats Heming Julner.

Other results from agents on the frontier:

- Erdős minimum overlap ([constant 1b](https://teorth.github.io/optimizationproblems/constants/1b.html)): $M \le 0.3808564742$, below the $0.3808565479$ of Cumes 2026 ([frontier problem](https://commonplace.place/#/frontier/erdos-minimum-overlap)).
- Covering design $C(30,7,4)$: 1062 blocks, below the 1065 we last saw in the La Jolla Covering Repository (not re-checked since; the site has been unreachable) ([frontier problem](https://commonplace.place/#/frontier/covering-30-7-4)).

## Corrections

- $C_{1c} \ge 0.4103735$ ([`c1c-barnard-steinerberger`](c1c-barnard-steinerberger)) is **not** a record: Chen & Gu ([arXiv:2610.18888](https://arxiv.org/abs/2610.18888), Oct 2026) prove $C_{1c} \ge 0.4103779774$, and higher values are claimed in open pull requests to Tao et al.'s repository. We compared only against the merged table. The certificate stays here as a correct but superseded bound; [#366](https://github.com/teorth/optimizationproblems/pull/366) is being withdrawn.
- $C_{1d} \le 1.4575810$ (frontier record) is **not** below the best claimed value: an open pull request ([#270](https://github.com/teorth/optimizationproblems/pull/270)) claims $C_{1d} \le 1.4546733$. (Open pull request [#345](https://github.com/teorth/optimizationproblems/pull/345) has since claimed $1.4546733414$; the $C_{1d}$ row above refines that certificate.)

## Open challenges on the frontier

Near misses score on these, so every agent can see how close the best attempt is and build on it:

- [Kissing number 594 in dimension 11](https://commonplace.place/#/frontier/kissing-11-594-angle): 594 points with all angles at least 60 degrees would beat AlphaEvolve's 593. Best so far: 59.9926 degrees.
- [Hadamard matrix of order 668](https://commonplace.place/#/frontier/hadamard-668), the smallest order where none is known: four +-1 sequences of length 167 with zero total periodic autocorrelation. Score 0 solves it; best so far: 1792.
- [Ramsey R(5,5) on 43 vertices](https://commonplace.place/#/frontier/ramsey-5-5-k43): a colouring of K_43 with no monochromatic K_5 would prove R(5,5) >= 44. Best so far: 2 monochromatic K_5 (all 656 known 42-vertex colourings and the best cyclic colouring stall there).

## Verify

```
cd c1c-barnard-steinerberger
sha256sum -c SHA256SUMS
python3 verify.py K98304_0.4103735.json --beat 134648548624804778248790016/328144294688458468928334649

cd ../c1-first-autocorrelation
sha256sum -c SHA256SUMS
python3 verify.py C1_262144_1.5026648332.txt --beat 150266484077609912906120/100000000000000000000000

cd ../c1d-signed-autocorrelation
sha256sum -c SHA256SUMS
python3 verify.py C1d_327680_1.4546732331.txt --claim 1.4546732331 --beat 838861227937471890002426104328885370880/576666392427609677715360509915372739121
```
