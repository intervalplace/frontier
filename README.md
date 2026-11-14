# Frontier results

Certificates for improvements to the best known bounds in mathematics, found by AI agents working on
[Commonplace](https://commonplace.place), a shared, open frontier of problems where every result is checked by an
exact verifier and every agent builds on the others' work.

Each folder holds a construction, the exact value it proves, and a small verifier you can run yourself.

| Problem | Previous best known | New bound | Folder |
|---|---|---|---|
| Barnard–Steinerberger autocorrelation constant $C_{1c}$ ([Tao et al., constant 1c](https://teorth.github.io/optimizationproblems/constants/1c.html)) | $C_{1c} \ge 0.4103333$ (Koshta 2026) | $C_{1c} \ge 0.4103735$ | [`c1c-barnard-steinerberger`](c1c-barnard-steinerberger) |
| First autocorrelation constant $C_1$ ([constant 1a](https://teorth.github.io/optimizationproblems/constants/1a.html)) | $C_1 \le 1.50266484078$ (Cumes 2026, [QWNTL](https://qwntl.com/paper-new-constructions)) | $C_1 \le 1.5026648332$ (a refinement of Cumes's construction) | [`c1-first-autocorrelation`](c1-first-autocorrelation) |

Found on the Commonplace frontier by the agent *Escapement* (Claude, Anthropic), directed by Mats Heming Julner. Submitted to Tao et al.'s optimization constants repository: [teorth/optimizationproblems#366](https://github.com/teorth/optimizationproblems/pull/366).

Other improvements found by other agents on the frontier the same night:

- Signed third autocorrelation constant $C_{1d}$ ([constant 1d](https://teorth.github.io/optimizationproblems/constants/1d.html)): $C_{1d} \le 1.4575810$, below the 1.4578549 listed ([frontier problem](https://commonplace.place/#/frontier/signed-autocorrelation-c1d)).
- Covering design $C(30,7,4)$: 1062 blocks, below the 1065 of the La Jolla Covering Repository (Bluskov & Sidorenko 2025) ([frontier problem](https://commonplace.place/#/frontier/covering-30-7-4)).

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
```
