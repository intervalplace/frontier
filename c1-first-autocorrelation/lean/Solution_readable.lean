-- READING COPY ONLY: the hexadecimal data of the step function is elided, so this file does not compile.
-- The real file is Solution.lean.
/-
Copyright 2025 The Formal Conjectures Authors.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-/
module
public import FormalConjecturesUtil

/-!
# Tao's Optimization constant 1a / An autocorrelation constant related to Sidon sets

*References:*
- [Tao's optimization constant 1a](https://teorth.github.io/optimizationproblems/constants/1a.html)
- [M2010] Matolcsi, Máté, and Carlos Vinuesa. "Improved bounds on the supremum of autoconvolutions."
  Journal of mathematical analysis and applications 372.2 (2010): 439-447. [arXiv:0907.1379](https://arxiv.org/abs/0907.1379)
- [Y2026] Yuksekgonul, Mert et al., "Learning to Discover at Test Time," 2026, [arXiv:2601.16175](https://arxiv.org/abs/2601.16175)
- [T2026] Together AI, "EinsteinArena-new-SOTA: State-of-the-art results on open math problems,"
  2026, [GitHub](https://github.com/togethercomputer/EinsteinArena-new-SOTA)
- [PBV2026] Piterbarg, Andrei, Jai Bajaj, and Derrick Vincent. "A multi-scale arcsine lower bound
  for the Sidon autocorrelation constant $C_{1a}$," 2026,
  [GitHub](https://github.com/AndreiPiterbarg/sidon-autocorrelation)
-/

@[expose] public section

open MeasureTheory Set

namespace Constant1a

/-- **Tao's Optimization constant 1a / An autocorrelation constant related to Sidon sets**:
The largest constant $C$ such that
$$\max_{-1/2 \le t \le 1/2} \int_{\mathbb{R}} f(t - x) f(x)\,dx
  \ge C \left(\int_{-1/4}^{1/4} f(x)\,dx\right)^2$$
for all nonnegative $f \colon \mathbb{R} \to \mathbb{R}$.

The inequality is only required for square-integrable $f$. This does not change the constant
(see [M2010, §2]), and it guarantees that every autoconvolution integral exists and that the
supremum over $t$ is finite.
This number is related to the maximal size of Sidon sets in additive combinatorics. -/
@[optimization_constant "1a"]
noncomputable def C1a : ℝ :=
  sSup {C : ℝ | ∀ ⦃f : ℝ → ℝ⦄, 0 ≤ f → MemLp f 2 →
    C * (∫ x in (- 1 / 4)..(1 / 4), f x) ^ 2
      ≤ sSup {∫ x, f (t - x) * f x | t ∈ Icc (- 1 / 2 : ℝ) (1 / 2)}}

/-- A lower bound, proven by Matolcsi-Vinuesa in [M2010]. -/
@[category research solved, AMS 5 11 26]
theorem c1a_lower_bound : 1.2748 ≤ C1a := by
  sorry

/-- An upper bound, proven by Yuksekgonul et al. in [Y2026]. -/
@[category research solved, AMS 5 11 26]
theorem c1a_upper_bound : C1a ≤ 1.5029 := by
  sorry

def mem_Ico_c1a.answer : ℝ := sorry

/-- How can the best known upper bound $1.5028628587\ldots$ of [T2026] be improved?

The step function published with [T2026] gives $1.50286285870\ldots$. The figure $1.502862$ in
the table of the first reference is this value truncated, so it is not a known upper bound. -/
@[category research open, question value, AMS 5 11 26]
theorem mem_Ico_c1a : mem_Ico_c1a.answer ∈ Set.Ico C1a 1.5028628587 := by
  sorry

def mem_Ioc_c1a.answer : ℝ := sorry

/-- How can the best known lower bound $1.292$ of [PBV2026] be improved? -/
@[category research open, question value, AMS 5 11 26]
theorem mem_Ioc_c1a : mem_Ioc_c1a.answer ∈ Set.Ioc 1.292 C1a := by
  sorry

def c1a_eq.answer : ℝ := sorry

/-- What is the exact value of the constant? -/
@[category research open, question value, AMS 5 11 26]
theorem c1a_eq : C1a = c1a_eq.answer := by
  sorry

namespace C1aCert

set_option maxRecDepth 100000
set_option maxHeartbeats 0

noncomputable section


/-- slot width in bits -/
def BW : ℕ := 60

/-- base-`2^BW` digit -/
def dig (x i : ℕ) : ℕ := x / (2 ^ BW) ^ i % 2 ^ BW

def sumDigits : ℕ → ℕ → ℕ
  | 0, x => x
  | d + 1, x => sumDigits d (x % (1 <<< (BW * 2 ^ d))) + sumDigits d (x >>> (BW * 2 ^ d))

def checkAll (T : ℕ) : ℕ → ℕ → Bool
  | 0, x => decide (x ≤ T)
  | d + 1, x => checkAll T d (x % (1 <<< (BW * 2 ^ d))) && checkAll T d (x >>> (BW * 2 ^ d))

/-- the `j`-th block of `2^e` digits -/
def slice (e x j : ℕ) : ℕ := (x >>> (BW * 2 ^ e * j)) % (1 <<< (BW * 2 ^ e))

theorem pow_split (d : ℕ) : (2 : ℕ) ^ (BW * 2 ^ d) = ((2 : ℕ) ^ BW) ^ (2 ^ d) := by
  rw [pow_mul]

theorem dig_mod (x i k : ℕ) (hik : i < k) : dig (x % ((2 ^ BW) ^ k)) i = dig x i := by
  unfold dig
  obtain ⟨r, rfl⟩ : ∃ r, k = i + (r + 1) := ⟨k - i - 1, by omega⟩
  rw [pow_add, Nat.mod_mul_right_div_self, pow_succ]
  exact Nat.mod_mul_left_mod _ _ _

theorem dig_div (x i k : ℕ) : dig (x / ((2 ^ BW) ^ k)) i = dig x (k + i) := by
  unfold dig
  rw [Nat.div_div_eq_div_mul, pow_add]

theorem checkAll_spec (T : ℕ) : ∀ (d x : ℕ), checkAll T d x = true → ∀ i < 2 ^ d, dig x i ≤ T := by
  intro d
  induction d with
  | zero =>
    intro x h i hi
    have hi0 : i = 0 := by omega
    subst hi0
    simp only [checkAll, decide_eq_true_eq] at h
    unfold dig
    simp only [pow_zero, Nat.div_one]
    exact (Nat.mod_le _ _).trans h
  | succ d ih =>
    intro x h i hi
    simp only [checkAll, Bool.and_eq_true] at h
    obtain ⟨h1, h2⟩ := h
    rw [Nat.one_shiftLeft, pow_split] at h1
    rw [Nat.shiftRight_eq_div_pow, pow_split] at h2
    by_cases hlt : i < 2 ^ d
    · have := ih _ h1 i hlt
      rwa [dig_mod _ _ _ hlt] at this
    · have hi2 : i - 2 ^ d < 2 ^ d := by rw [pow_succ] at hi; omega
      have := ih _ h2 (i - 2 ^ d) hi2
      rw [dig_div] at this
      have e : 2 ^ d + (i - 2 ^ d) = i := by omega
      rwa [e] at this

theorem sumDigits_spec : ∀ (d x : ℕ), x < ((2 : ℕ) ^ BW) ^ (2 ^ d) →
    sumDigits d x = ∑ i ∈ Finset.range (2 ^ d), dig x i := by
  intro d
  induction d with
  | zero =>
    intro x hx
    simp only [pow_zero, pow_one] at hx
    simp [sumDigits, dig, Nat.mod_eq_of_lt hx]
  | succ d ih =>
    intro x hx
    simp only [sumDigits]
    rw [Nat.one_shiftLeft, Nat.shiftRight_eq_div_pow, pow_split]
    have hpos : 0 < ((2 : ℕ) ^ BW) ^ (2 ^ d) := by positivity
    have hlo : x % ((2 : ℕ) ^ BW) ^ (2 ^ d) < ((2 : ℕ) ^ BW) ^ (2 ^ d) := Nat.mod_lt _ hpos
    have hhi : x / ((2 : ℕ) ^ BW) ^ (2 ^ d) < ((2 : ℕ) ^ BW) ^ (2 ^ d) := by
      apply Nat.div_lt_of_lt_mul
      rw [← pow_add]
      have : 2 ^ d + 2 ^ d = 2 ^ (d + 1) := by rw [pow_succ]; omega
      rwa [this]
    rw [ih _ hlo, ih _ hhi]
    have e : 2 ^ (d + 1) = 2 ^ d + 2 ^ d := by rw [pow_succ]; omega
    rw [e, Finset.sum_range_add]
    congr 1
    · apply Finset.sum_congr rfl
      intro i hi
      exact dig_mod _ _ _ (Finset.mem_range.mp hi)
    · apply Finset.sum_congr rfl
      intro i _
      exact dig_div _ _ _

theorem slice_lt (e x j : ℕ) : slice e x j < ((2 : ℕ) ^ BW) ^ (2 ^ e) := by
  unfold slice
  rw [Nat.one_shiftLeft, pow_split]
  exact Nat.mod_lt _ (by positivity)

theorem dig_slice (e x j i : ℕ) (hi : i < 2 ^ e) : dig (slice e x j) i = dig x (2 ^ e * j + i) := by
  unfold slice
  rw [Nat.one_shiftLeft, Nat.shiftRight_eq_div_pow, pow_split, dig_mod _ _ _ hi]
  have : (2 : ℕ) ^ (BW * 2 ^ e * j) = ((2 : ℕ) ^ BW) ^ (2 ^ e * j) := by
    rw [mul_assoc, pow_mul]
  rw [this, dig_div]

/-- digit expansion -/
theorem digits_expansion : ∀ (n x : ℕ), x < ((2 : ℕ) ^ BW) ^ n →
    x = ∑ i ∈ Finset.range n, dig x i * ((2 : ℕ) ^ BW) ^ i := by
  intro n
  induction n with
  | zero => intro x hx; simp at hx; simp [hx]
  | succ n ih =>
    intro x hx
    have hB : 0 < (2 : ℕ) ^ BW := by positivity
    have hq : x / 2 ^ BW < ((2 : ℕ) ^ BW) ^ n := by
      apply Nat.div_lt_of_lt_mul
      rwa [pow_succ, mul_comm] at hx
    have h1 := ih _ hq
    rw [Finset.sum_range_succ']
    have h2 : ∀ i, dig (x / 2 ^ BW) i = dig x (i + 1) := by
      intro i
      have := dig_div x i 1
      rw [pow_one] at this
      rw [this, add_comm]
    have h3 : ∑ i ∈ Finset.range n, dig x (i + 1) * ((2 : ℕ) ^ BW) ^ (i + 1)
        = (x / 2 ^ BW) * 2 ^ BW := by
      conv_rhs => rw [h1]
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro i _
      rw [h2, pow_succ, mul_assoc]
    rw [h3]
    have h4 : dig x 0 * ((2 : ℕ) ^ BW) ^ 0 = x % 2 ^ BW := by simp [dig]
    rw [h4]
    exact (Nat.div_add_mod' x (2 ^ BW)).symm

/-- digits of a sum `∑ c m * B^m` with all `c m < B` -/
theorem dig_of_sum (c : ℕ → ℕ) (hc : ∀ m, c m < 2 ^ BW) : ∀ (j K : ℕ), j < K →
    dig (∑ m ∈ Finset.range K, c m * ((2 : ℕ) ^ BW) ^ m) j = c j := by
  intro j
  induction j generalizing c with
  | zero =>
    intro K hK
    obtain ⟨K, rfl⟩ : ∃ K', K = K' + 1 := ⟨K - 1, by omega⟩
    rw [Finset.sum_range_succ']
    unfold dig
    simp only [pow_zero, Nat.div_one, mul_one]
    have : ∑ m ∈ Finset.range K, c (m + 1) * ((2 : ℕ) ^ BW) ^ (m + 1)
        = 2 ^ BW * ∑ m ∈ Finset.range K, c (m + 1) * ((2 : ℕ) ^ BW) ^ m := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro i _
      rw [pow_succ]; ring
    rw [this, Nat.mul_add_mod]
    exact Nat.mod_eq_of_lt (hc 0)
  | succ j ih =>
    intro K hK
    obtain ⟨K, rfl⟩ : ∃ K', K = K' + 1 := ⟨K - 1, by omega⟩
    rw [Finset.sum_range_succ']
    have hsum : ∑ m ∈ Finset.range K, c (m + 1) * ((2 : ℕ) ^ BW) ^ (m + 1)
        = 2 ^ BW * ∑ m ∈ Finset.range K, c (m + 1) * ((2 : ℕ) ^ BW) ^ m := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro i _
      rw [pow_succ]; ring
    have hdiv : (∑ m ∈ Finset.range K, c (m + 1) * ((2 : ℕ) ^ BW) ^ (m + 1) + c 0 * ((2 : ℕ) ^ BW) ^ 0)
        / 2 ^ BW = ∑ m ∈ Finset.range K, c (m + 1) * ((2 : ℕ) ^ BW) ^ m := by
      rw [hsum, pow_zero, mul_one, Nat.mul_add_div (by positivity), Nat.div_eq_of_lt (hc 0), add_zero]
    have key := dig_div (∑ m ∈ Finset.range K, c (m + 1) * ((2 : ℕ) ^ BW) ^ (m + 1) + c 0 * ((2 : ℕ) ^ BW) ^ 0) j 1
    rw [pow_one, hdiv, add_comm 1 j] at key
    rw [← key]
    exact ih (fun m => c (m + 1)) (fun m => hc (m + 1)) K (by omega)



open Finset

/-- the coefficient `∑_{i+j=m, i,j<n} v i * v j` -/
def conv (n : ℕ) (v : ℕ → ℕ) (m : ℕ) : ℕ :=
  ∑ p ∈ (range n ×ˢ range n).filter (fun p => p.1 + p.2 = m), v p.1 * v p.2

theorem conv_le_sq (n : ℕ) (v : ℕ → ℕ) (m : ℕ) : conv n v m ≤ (∑ i ∈ range n, v i) ^ 2 := by
  unfold conv
  calc ∑ p ∈ (range n ×ˢ range n).filter (fun p => p.1 + p.2 = m), v p.1 * v p.2
      ≤ ∑ p ∈ range n ×ˢ range n, v p.1 * v p.2 :=
        Finset.sum_le_sum_of_subset (Finset.filter_subset _ _)
    _ = (∑ i ∈ range n, v i) ^ 2 := by
        rw [sq, Finset.sum_mul_sum, Finset.sum_product]

theorem sq_eq_sum_conv (n : ℕ) (v : ℕ → ℕ) (B : ℕ) :
    (∑ i ∈ range n, v i * B ^ i) * (∑ i ∈ range n, v i * B ^ i)
      = ∑ m ∈ range (2 * n), conv n v m * B ^ m := by
  rw [Finset.sum_mul_sum, ← Finset.sum_product']
  have hmaps : ∀ p ∈ range n ×ˢ range n, p.1 + p.2 ∈ range (2 * n) := by
    intro p hp
    rw [Finset.mem_product, Finset.mem_range, Finset.mem_range] at hp
    rw [Finset.mem_range]; omega
  rw [← Finset.sum_fiberwise_of_maps_to hmaps]
  apply Finset.sum_congr rfl
  intro m _
  unfold conv
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro p hp
  rw [Finset.mem_filter] at hp
  rw [← hp.2, pow_add]; ring

/-- regrouping a weighted double sum by `i + j` (real weights) -/
theorem sum_weight_conv (n : ℕ) (v : ℕ → ℕ) (w : ℕ → ℝ) :
    ∑ i ∈ range n, ∑ j ∈ range n, ((v i : ℝ) * (v j : ℝ)) * w (i + j)
      = ∑ m ∈ range (2 * n), (conv n v m : ℝ) * w m := by
  rw [← Finset.sum_product']
  have hmaps : ∀ p ∈ range n ×ˢ range n, p.1 + p.2 ∈ range (2 * n) := by
    intro p hp
    rw [Finset.mem_product, Finset.mem_range, Finset.mem_range] at hp
    rw [Finset.mem_range]; omega
  rw [← Finset.sum_fiberwise_of_maps_to hmaps]
  apply Finset.sum_congr rfl
  intro m _
  unfold conv
  rw [Nat.cast_sum, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro p hp
  rw [Finset.mem_filter] at hp
  rw [← hp.2]; push_cast; ring


/-- overlap length of two intervals of length `h` whose left ends differ by `q - h` -/
noncomputable def Wh (h q : ℝ) : ℝ := max 0 (min h q - max 0 (q - h))

theorem Wh_nonneg (h q : ℝ) : 0 ≤ Wh h q := le_max_left _ _

theorem sum_Wh_le (h : ℝ) (hh : 0 ≤ h) (q : ℝ) : ∀ K : ℕ,
    ∑ m ∈ range K, Wh h (q - m * h) ≤ max 0 (min h (((K : ℝ) + 1) * h - q)) := by
  intro K
  induction K with
  | zero => simp
  | succ K ih =>
    rw [Finset.sum_range_succ]
    have key : max 0 (min h (((K : ℝ) + 1) * h - q)) + Wh h (q - K * h)
        ≤ max 0 (min h ((((K + 1 : ℕ) : ℝ) + 1) * h - q)) := by
      unfold Wh
      push_cast
      simp only [max_def, min_def]
      split_ifs <;> nlinarith
    linarith

theorem sum_Wh_le_h (h : ℝ) (hh : 0 ≤ h) (q : ℝ) (K : ℕ) : ∑ m ∈ range K, Wh h (q - m * h) ≤ h := by
  refine (sum_Wh_le h hh q K).trans ?_
  simp only [max_def, min_def]
  split_ifs <;> linarith

open MeasureTheory Set

/-- left end of the `i`-th cell; the cells have width `1/(2n)` and tile `(-1/4, 1/4]` -/
noncomputable def a (n : ℕ) (i : ℕ) : ℝ := (i : ℝ) * (1 / (2 * (n : ℝ))) - 1 / 4

def cell (n i : ℕ) : Set ℝ := Set.Ioc (a n i) (a n (i + 1))

/-- the step function with heights `v 0, …, v (n-1)` -/
noncomputable def f0 (n : ℕ) (v : ℕ → ℕ) (x : ℝ) : ℝ :=
  ∑ i ∈ range n, (v i : ℝ) * Set.indicator (cell n i) 1 x

theorem a_succ (n i : ℕ) : a n (i + 1) = a n i + 1 / (2 * (n : ℝ)) := by
  unfold a; push_cast; ring

theorem a_mono (n : ℕ) {i j : ℕ} (hij : i ≤ j) : a n i ≤ a n j := by
  unfold a
  have h : (0 : ℝ) ≤ 1 / (2 * (n : ℝ)) := by positivity
  have : (i : ℝ) ≤ j := by exact_mod_cast hij
  nlinarith

theorem a_zero (n : ℕ) : a n 0 = -(1 / 4) := by unfold a; simp

theorem a_self (n : ℕ) (hn : 0 < n) : a n n = 1 / 4 := by
  unfold a
  have : (n : ℝ) ≠ 0 := by positivity
  field_simp; ring

theorem f0_nonneg (n : ℕ) (v : ℕ → ℕ) : 0 ≤ f0 n v := by
  intro x
  unfold f0
  apply Finset.sum_nonneg
  intro i _
  apply mul_nonneg (Nat.cast_nonneg _)
  exact Set.indicator_nonneg (fun _ _ => zero_le_one) _

theorem f0_memLp (n : ℕ) (v : ℕ → ℕ) : MemLp (f0 n v) 2 := by
  unfold f0
  apply memLp_finsetSum
  intro i _
  apply MemLp.const_mul
  have : (Set.indicator (cell n i) (1 : ℝ → ℝ)) = Set.indicator (cell n i) (fun _ => (1 : ℝ)) := rfl
  rw [this]
  exact memLp_indicator_const 2 measurableSet_Ioc 1 (Or.inr measure_Ioc_lt_top.ne)

theorem intervalIntegral_indicator_Ioc (c d l r : ℝ) (hcl : c ≤ l) (hlr : l ≤ r) (hrd : r ≤ d) :
    ∫ x in c..d, Set.indicator (Set.Ioc l r) (1 : ℝ → ℝ) x = r - l := by
  rw [intervalIntegral.integral_of_le (hcl.trans (hlr.trans hrd))]
  rw [MeasureTheory.setIntegral_indicator measurableSet_Ioc]
  have : Set.Ioc c d ∩ Set.Ioc l r = Set.Ioc l r :=
    Set.inter_eq_right.mpr (Set.Ioc_subset_Ioc hcl hrd)
  rw [this, ← intervalIntegral.integral_of_le hlr]
  simp

theorem f0_integral (n : ℕ) (hn : 0 < n) (v : ℕ → ℕ) :
    ∫ x in (-1 / 4 : ℝ)..(1 / 4), f0 n v x = (∑ i ∈ range n, (v i : ℝ)) * (1 / (2 * (n : ℝ))) := by
  unfold f0
  rw [intervalIntegral.integral_finsetSum]
  · rw [Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro i hi
    rw [intervalIntegral.integral_const_mul]
    congr 1
    have hi' : i + 1 ≤ n := Finset.mem_range.mp hi
    have h1 : (-1 / 4 : ℝ) ≤ a n i := by
      have := a_mono n (Nat.zero_le i); rw [a_zero] at this; linarith
    have h2 : a n i ≤ a n (i + 1) := a_mono n (Nat.le_succ i)
    have h3 : a n (i + 1) ≤ 1 / 4 := by
      have := a_mono n hi'; rwa [a_self n hn] at this
    unfold cell
    rw [intervalIntegral_indicator_Ioc _ _ _ _ h1 h2 h3, a_succ]; ring
  · intro i _
    apply IntervalIntegrable.const_mul
    rw [intervalIntegrable_iff]
    exact (integrableOn_const (by simp)).indicator measurableSet_Ioc

/-- the set where both indicators are non-zero -/
def E (n : ℕ) (t : ℝ) (i j : ℕ) : Set ℝ := ((fun x : ℝ => t - x) ⁻¹' cell n i) ∩ cell n j

theorem measurableSet_E (n : ℕ) (t : ℝ) (i j : ℕ) : MeasurableSet (E n t i j) :=
  (measurableSet_preimage (measurable_const.sub measurable_id) measurableSet_Ioc).inter measurableSet_Ioc

theorem E_subset (n : ℕ) (t : ℝ) (i j : ℕ) :
    E n t i j ⊆ Set.Icc (max (a n j) (t - a n (i + 1))) (min (a n (j + 1)) (t - a n i)) := by
  intro x hx
  obtain ⟨h1, h2⟩ := hx
  simp only [Set.mem_preimage, cell, Set.mem_Ioc] at h1 h2
  refine ⟨max_le h2.1.le (by linarith), le_min h2.2 (by linarith)⟩

theorem overlap_eq (A h q : ℝ) :
    max (min (A + h) (A + q) - max A (A + (q - h))) 0 = Wh h q := by
  unfold Wh
  have h1 : min (A + h) (A + q) = A + min h q := (min_add_add_left A h q)
  have h2 : max A (A + (q - h)) = A + max 0 (q - h) := by
    have := max_add_add_left A 0 (q - h)
    rw [add_zero] at this
    exact this
  rw [h1, h2, max_comm]
  congr 1
  ring

theorem vol_E_le (n : ℕ) (t : ℝ) (i j : ℕ) :
    volume.real (E n t i j) ≤ Wh (1 / (2 * (n : ℝ))) (t + 1 / 2 - ((i + j : ℕ) : ℝ) * (1 / (2 * (n : ℝ)))) := by
  have hmono := measureReal_mono (μ := volume) (E_subset n t i j) (by
    rw [Real.volume_Icc]; exact ENNReal.ofReal_ne_top)
  refine hmono.trans (le_of_eq ?_)
  rw [Measure.real, Real.volume_Icc, ENNReal.toReal_ofReal']
  obtain ⟨q, hq⟩ : ∃ q : ℝ, q = t - a n i - a n j := ⟨_, rfl⟩
  have e1 : t - a n (i + 1) = a n j + (q - 1 / (2 * (n : ℝ))) := by
    rw [a_succ, hq]; ring
  have e2 : t - a n i = a n j + q := by
    rw [hq]; ring
  rw [a_succ n j, e1, e2, overlap_eq, hq]
  congr 1
  unfold a
  push_cast
  ring

theorem f0_mul_expand (n : ℕ) (v : ℕ → ℕ) (t x : ℝ) :
    f0 n v (t - x) * f0 n v x
      = ∑ i ∈ range n, ∑ j ∈ range n, ((v i : ℝ) * (v j : ℝ)) * Set.indicator (E n t i j) 1 x := by
  unfold f0
  rw [Finset.sum_mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  apply Finset.sum_congr rfl
  intro j _
  have : Set.indicator (cell n i) (1 : ℝ → ℝ) (t - x) * Set.indicator (cell n j) (1 : ℝ → ℝ) x
      = Set.indicator (E n t i j) (1 : ℝ → ℝ) x := by
    unfold E
    by_cases h1 : t - x ∈ cell n i <;> by_cases h2 : x ∈ cell n j <;>
      simp [Set.indicator, h1, h2]
  rw [← this]; ring

theorem conv_integral_le (n : ℕ) (v : ℕ → ℕ) (T : ℕ) (hT : ∀ m, conv n v m ≤ T) (t : ℝ) :
    ∫ x, f0 n v (t - x) * f0 n v x ≤ (T : ℝ) * (1 / (2 * (n : ℝ))) := by
  have hh : (0 : ℝ) ≤ 1 / (2 * (n : ℝ)) := by positivity
  have hint : ∀ i j : ℕ, Integrable (fun x : ℝ => ((v i : ℝ) * (v j : ℝ)) * Set.indicator (E n t i j) 1 x) := by
    intro i j
    apply Integrable.const_mul
    have hfin : volume (E n t i j) ≠ ⊤ := by
      refine ne_top_of_le_ne_top ?_ (measure_mono (E_subset n t i j))
      rw [Real.volume_Icc]; exact ENNReal.ofReal_ne_top
    exact (integrable_indicator_iff (measurableSet_E n t i j)).mpr (integrableOn_const hfin)
  have h1 : ∫ x, f0 n v (t - x) * f0 n v x
      = ∑ i ∈ range n, ∑ j ∈ range n, ((v i : ℝ) * (v j : ℝ)) * volume.real (E n t i j) := by
    simp_rw [f0_mul_expand]
    rw [integral_finsetSum _ (fun i _ => integrable_finsetSum _ (fun j _ => hint i j))]
    apply Finset.sum_congr rfl
    intro i _
    rw [integral_finsetSum _ (fun j _ => hint i j)]
    apply Finset.sum_congr rfl
    intro j _
    rw [integral_const_mul, integral_indicator_one (measurableSet_E n t i j)]
  rw [h1]
  set h : ℝ := 1 / (2 * (n : ℝ)) with hdef
  calc ∑ i ∈ range n, ∑ j ∈ range n, ((v i : ℝ) * (v j : ℝ)) * volume.real (E n t i j)
      ≤ ∑ i ∈ range n, ∑ j ∈ range n, ((v i : ℝ) * (v j : ℝ)) * Wh h (t + 1 / 2 - ((i + j : ℕ) : ℝ) * h) := by
        apply Finset.sum_le_sum; intro i _
        apply Finset.sum_le_sum; intro j _
        exact mul_le_mul_of_nonneg_left (vol_E_le n t i j) (by positivity)
    _ = ∑ m ∈ range (2 * n), (conv n v m : ℝ) * Wh h (t + 1 / 2 - (m : ℝ) * h) :=
        sum_weight_conv n v (fun m => Wh h (t + 1 / 2 - (m : ℝ) * h))
    _ ≤ ∑ m ∈ range (2 * n), (T : ℝ) * Wh h (t + 1 / 2 - (m : ℝ) * h) := by
        apply Finset.sum_le_sum; intro m _
        exact mul_le_mul_of_nonneg_right (by exact_mod_cast hT m) (Wh_nonneg _ _)
    _ = (T : ℝ) * ∑ m ∈ range (2 * n), Wh h (t + 1 / 2 - (m : ℝ) * h) := by rw [Finset.mul_sum]
    _ ≤ (T : ℝ) * h := mul_le_mul_of_nonneg_left (sum_Wh_le_h h hh _ _) (Nat.cast_nonneg _)


theorem sum_range_mul' (f : ℕ → ℕ) (K : ℕ) : ∀ J : ℕ,
    ∑ k ∈ range (K * J), f k = ∑ j ∈ range J, ∑ i ∈ range K, f (K * j + i) := by
  intro J
  induction J with
  | zero => simp
  | succ J ih =>
    rw [Nat.mul_succ, Finset.sum_range_add, ih, Finset.sum_range_succ]

/-- From the kernel-checked facts about the packed integer `P` to the bound on the convolution coefficients. -/
theorem cert_core (e JP P T : ℕ) (s : ℕ → ℕ)
    (hP : P < ((2 : ℕ) ^ BW) ^ (2 ^ e * JP))
    (hS : ∀ j < JP, sumDigits e (slice e P j) = s j)
    (hQ : ∀ j < 2 * JP, checkAll T e (slice e (P * P) j) = true)
    (hB : (∑ j ∈ range JP, s j) ^ 2 < 2 ^ BW) :
    (∑ i ∈ range (2 ^ e * JP), dig P i = ∑ j ∈ range JP, s j) ∧
      ∀ m, conv (2 ^ e * JP) (fun i => dig P i) m ≤ T := by
  have hsum : ∑ i ∈ range (2 ^ e * JP), dig P i = ∑ j ∈ range JP, s j := by
    rw [sum_range_mul']
    apply Finset.sum_congr rfl
    intro j hj
    rw [← hS j (Finset.mem_range.mp hj), sumDigits_spec e _ (slice_lt e P j)]
    apply Finset.sum_congr rfl
    intro i hi
    exact (dig_slice e P j i (Finset.mem_range.mp hi)).symm
  refine ⟨hsum, ?_⟩
  intro m
  set n := 2 ^ e * JP with hn
  have hlt : ∀ m, conv n (fun i => dig P i) m < 2 ^ BW := by
    intro m
    refine lt_of_le_of_lt (conv_le_sq n _ m) ?_
    rw [hsum]; exact hB
  by_cases hm : m < 2 * n
  · have hPP : P * P = ∑ m ∈ range (2 * n), conv n (fun i => dig P i) m * ((2 : ℕ) ^ BW) ^ m := by
      rw [← sq_eq_sum_conv n (fun i => dig P i) (2 ^ BW)]
      rw [← digits_expansion n P hP]
    have hd : dig (P * P) m = conv n (fun i => dig P i) m := by
      rw [hPP]
      exact dig_of_sum _ hlt m (2 * n) hm
    rw [← hd]
    have hj : m / 2 ^ e < 2 * JP := by
      apply Nat.div_lt_of_lt_mul
      rw [hn] at hm
      calc m < 2 * (2 ^ e * JP) := hm
        _ = 2 ^ e * (2 * JP) := by ring
    have hi : m % 2 ^ e < 2 ^ e := Nat.mod_lt _ (by positivity)
    have := checkAll_spec T e _ (hQ _ hj) _ hi
    rw [dig_slice e (P * P) _ _ hi, Nat.div_add_mod] at this
    exact this
  · have : conv n (fun i => dig P i) m = 0 := by
      unfold conv
      apply Finset.sum_eq_zero
      intro p hp
      rw [Finset.mem_filter, Finset.mem_product, Finset.mem_range, Finset.mem_range] at hp
      omega
    rw [this]; exact Nat.zero_le _

/-- The analytic step: a step function whose convolution coefficients are at most `T` bounds `C1a`. -/
theorem C1a_le_of_cert (n : ℕ) (hn : 0 < n) (v : ℕ → ℕ) (T : ℕ) (hT : ∀ m, conv n v m ≤ T)
    (hS : 0 < ∑ i ∈ range n, v i) :
    C1a ≤ (2 * (n : ℝ) * (T : ℝ)) / ((∑ i ∈ range n, v i : ℕ) : ℝ) ^ 2 := by
  have hSr : (0 : ℝ) < ((∑ i ∈ range n, v i : ℕ) : ℝ) := by exact_mod_cast hS
  have hnr : (0 : ℝ) < (n : ℝ) := by exact_mod_cast hn
  unfold C1a
  apply csSup_le
  · refine ⟨0, ?_⟩
    intro f hf _
    rw [zero_mul]
    apply Real.sSup_nonneg
    rintro y ⟨t, _, rfl⟩
    exact integral_nonneg (fun x => mul_nonneg (hf _) (hf _))
  · intro C hC
    have h := hC (f0_nonneg n v) (f0_memLp n v)
    have hsup := Real.sSup_le (s := {x | ∃ t ∈ Icc (-1 / 2 : ℝ) (1 / 2), ∫ x, f0 n v (t - x) * f0 n v x = x})
      (a := (T : ℝ) * (1 / (2 * (n : ℝ))))
      (by rintro y ⟨t, _, rfl⟩; exact conv_integral_le n v T hT t) (by positivity)
    have h2 := h.trans hsup
    rw [f0_integral n hn v, ← Nat.cast_sum] at h2
    set Sr : ℝ := ((∑ i ∈ range n, v i : ℕ) : ℝ) with hSrdef
    rw [le_div_iff₀ (by positivity)]
    have hh : (0 : ℝ) < 1 / (2 * (n : ℝ)) := by positivity
    have e1 : (C * Sr ^ 2 * (1 / (2 * (n : ℝ)))) * (1 / (2 * (n : ℝ))) ≤ (T : ℝ) * (1 / (2 * (n : ℝ))) := by
      calc (C * Sr ^ 2 * (1 / (2 * (n : ℝ)))) * (1 / (2 * (n : ℝ)))
          = C * (Sr * (1 / (2 * (n : ℝ)))) ^ 2 := by ring
        _ ≤ (T : ℝ) * (1 / (2 * (n : ℝ))) := h2
    have e2 : C * Sr ^ 2 * (1 / (2 * (n : ℝ))) ≤ (T : ℝ) := le_of_mul_le_mul_right e1 hh
    have e3 : C * Sr ^ 2 = (C * Sr ^ 2 * (1 / (2 * (n : ℝ)))) * (2 * (n : ℝ)) := by field_simp
    rw [e3]
    calc (C * Sr ^ 2 * (1 / (2 * (n : ℝ)))) * (2 * (n : ℝ)) ≤ (T : ℝ) * (2 * (n : ℝ)) :=
          mul_le_mul_of_nonneg_right e2 (by positivity)
      _ = 2 * (n : ℝ) * (T : ℝ) := by ring

noncomputable def blk0 : ℕ := (0x034c01a4…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0xb4747252…(16384 hex digits)… + ((0xe712f051…(16384 hex digits)… + ((0x00020d01…(16384 hex digits)… + ((0xfd54d656…(16384 hex digits)… + ((0x468581b5…(16384 hex digits)… + ((0x489c2db8…(16384 hex digits)… + ((0xf58f14ae…(16384 hex digits)… + ((0x09025371…(16384 hex digits)… + ((0x9319a468…(16384 hex digits)… + ((0x72d84c51…(16384 hex digits)… + ((0x21d3fbce…(16384 hex digits)… + ((0x26c08d62…(16384 hex digits)… + ((0x89612378…(16384 hex digits)… + ((0x632f484f…(16384 hex digits)… + (0) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)
noncomputable def blk1 : ℕ := (0x3904a49b…(16384 hex digits)… + ((0xcd431930…(16384 hex digits)… + ((0xfc0b072c…(16384 hex digits)… + ((0xcda3988d…(16384 hex digits)… + ((0x13ac2f50…(16384 hex digits)… + ((0x05851070…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0xfd207000…(16384 hex digits)… + ((0x8f0d743c…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x9343320d…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0xe92da4e9…(16384 hex digits)… + (0) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)
noncomputable def blk2 : ℕ := (0x7319d265…(16384 hex digits)… + ((0xe943a8ce…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x55e955e9…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x6e01643f…(16384 hex digits)… + ((0xeae3398e…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x14e44e91…(16384 hex digits)… + ((0xcbb32d6c…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + (0) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)
noncomputable def blk3 : ℕ := (0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0x00000000…(16384 hex digits)… + ((0xf86f88c3…(16384 hex digits)… + ((0x7195d157…(16384 hex digits)… + ((0x3de54ed6…(16384 hex digits)… + (0) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)) <<< 65536)
/-- low 14 bits of the heights, in 14-bit slots -/
noncomputable def L : ℕ := (blk0 + ((blk1 + ((blk2 + ((blk3 + (0) <<< 1048576)) <<< 1048576)) <<< 1048576)) <<< 1048576)
/-- the heights of 15 bits or more: 24-bit entries `gap * 64 + high part`, in increasing cell order -/
noncomputable def X : ℕ := 0x43000043…(14126 hex digits)…
/-- one doubling step of the unpacking: in every block of `2 * 2^k` slots the upper half moves up -/
noncomputable def stage (k x : ℕ) : ℕ :=
  (x &&& (((1 <<< (2 ^ k * 14)) - 1) * (((1 <<< (2 ^ 18 * 60)) - 1) / ((1 <<< (2 ^ (k + 1) * 60)) - 1)))) +
    ((x &&& ((((1 <<< (2 ^ k * 14)) - 1) * (((1 <<< (2 ^ 18 * 60)) - 1) / ((1 <<< (2 ^ (k + 1) * 60)) - 1)))
      <<< (2 ^ k * 14))) <<< (2 ^ k * 46))

noncomputable def unpack : ℕ → ℕ → ℕ
  | 0, x => x
  | k + 1, x => unpack k (stage k x)

/-- sum of the gap fields of the `2^d` lowest entries -/
noncomputable def gapSum : ℕ → ℕ → ℕ
  | 0, lst => (lst >>> 6) % 262144
  | d + 1, lst => gapSum d (lst % (1 <<< (24 * 2 ^ d))) + gapSum d (lst >>> (24 * 2 ^ d))

/-- the high parts of the `2^d` lowest entries, placed in 60-bit slots relative to the cell before the first -/
noncomputable def excB : ℕ → ℕ → ℕ
  | 0, lst => (lst % 64) <<< (60 * ((lst >>> 6) % 262144))
  | d + 1, lst => excB d (lst % (1 <<< (24 * 2 ^ d))) +
      (excB d (lst >>> (24 * 2 ^ d))) <<< (60 * gapSum d (lst % (1 <<< (24 * 2 ^ d))))

/-- The packed step function: `P = ∑ i < 2^18, v i * 2^(60 i)`. -/
noncomputable def P : ℕ := unpack 18 L + (excB 12 X) <<< 14
/-- The largest coefficient of the autoconvolution, `max_m ∑_(i+j=m) v i * v j`. -/
def T : ℕ := 1215979245790
def sv (j : ℕ) : ℕ := if j = 0 then 113451105 else if j = 1 then 111486517 else if j = 2 then 132797775 else if j = 3 then 23191925 else if j = 4 then 87809943 else if j = 5 then 50586821 else if j = 6 then 21541757 else 110486688

theorem kP : P < 1 <<< (BW * 2 ^ 15 * 8) := by decide +kernel
theorem kS0 : sumDigits 15 (slice 15 P 0) = 113451105 := by decide +kernel
theorem kS1 : sumDigits 15 (slice 15 P 1) = 111486517 := by decide +kernel
theorem kS2 : sumDigits 15 (slice 15 P 2) = 132797775 := by decide +kernel
theorem kS3 : sumDigits 15 (slice 15 P 3) = 23191925 := by decide +kernel
theorem kS4 : sumDigits 15 (slice 15 P 4) = 87809943 := by decide +kernel
theorem kS5 : sumDigits 15 (slice 15 P 5) = 50586821 := by decide +kernel
theorem kS6 : sumDigits 15 (slice 15 P 6) = 21541757 := by decide +kernel
theorem kS7 : sumDigits 15 (slice 15 P 7) = 110486688 := by decide +kernel
theorem kQ0 : checkAll T 15 (slice 15 (P * P) 0) = true := by decide +kernel
theorem kQ1 : checkAll T 15 (slice 15 (P * P) 1) = true := by decide +kernel
theorem kQ2 : checkAll T 15 (slice 15 (P * P) 2) = true := by decide +kernel
theorem kQ3 : checkAll T 15 (slice 15 (P * P) 3) = true := by decide +kernel
theorem kQ4 : checkAll T 15 (slice 15 (P * P) 4) = true := by decide +kernel
theorem kQ5 : checkAll T 15 (slice 15 (P * P) 5) = true := by decide +kernel
theorem kQ6 : checkAll T 15 (slice 15 (P * P) 6) = true := by decide +kernel
theorem kQ7 : checkAll T 15 (slice 15 (P * P) 7) = true := by decide +kernel
theorem kQ8 : checkAll T 15 (slice 15 (P * P) 8) = true := by decide +kernel
theorem kQ9 : checkAll T 15 (slice 15 (P * P) 9) = true := by decide +kernel
theorem kQ10 : checkAll T 15 (slice 15 (P * P) 10) = true := by decide +kernel
theorem kQ11 : checkAll T 15 (slice 15 (P * P) 11) = true := by decide +kernel
theorem kQ12 : checkAll T 15 (slice 15 (P * P) 12) = true := by decide +kernel
theorem kQ13 : checkAll T 15 (slice 15 (P * P) 13) = true := by decide +kernel
theorem kQ14 : checkAll T 15 (slice 15 (P * P) 14) = true := by decide +kernel
theorem kQ15 : checkAll T 15 (slice 15 (P * P) 15) = true := by decide +kernel

theorem cert_facts :
    (∑ i ∈ range (2 ^ 15 * 8), dig P i = 651352531) ∧ ∀ m, conv (2 ^ 15 * 8) (fun i => dig P i) m ≤ T := by
  have hP : P < ((2 : ℕ) ^ BW) ^ (2 ^ 15 * 8) := by
    have h := kP
    rwa [Nat.one_shiftLeft, mul_assoc, pow_mul] at h
  have hS : ∀ j < 8, sumDigits 15 (slice 15 P j) = sv j := by
    intro j hj
    interval_cases j
    · exact kS0
    · exact kS1
    · exact kS2
    · exact kS3
    · exact kS4
    · exact kS5
    · exact kS6
    · exact kS7
  have hQ : ∀ j < 2 * 8, checkAll T 15 (slice 15 (P * P) j) = true := by
    intro j hj
    interval_cases j
    · exact kQ0
    · exact kQ1
    · exact kQ2
    · exact kQ3
    · exact kQ4
    · exact kQ5
    · exact kQ6
    · exact kQ7
    · exact kQ8
    · exact kQ9
    · exact kQ10
    · exact kQ11
    · exact kQ12
    · exact kQ13
    · exact kQ14
    · exact kQ15
  have hsv : ∑ j ∈ range 8, sv j = 651352531 := by simp [Finset.sum_range_succ, sv]
  have hB : (∑ j ∈ range 8, sv j) ^ 2 < 2 ^ BW := by
    rw [hsv]; norm_num [BW]
  obtain ⟨hsum, hconv⟩ := cert_core 15 8 P T sv hP hS hQ hB
  exact ⟨hsum.trans hsv, hconv⟩

end

end C1aCert

/-- An upper bound from an explicit step function with $2^{18}$ equal steps on $[-1/4,1/4]$:
the record of the Commonplace frontier board `autocorrelation-c1-large` (exact value $1.50266482948\ldots$),
with its heights rounded to $20$ bits, which gives $1.50267087879\ldots$.
This is below the bound $1.5028628587$ asked for in `mem_Ico_c1a`. -/
@[category research solved, AMS 5 11 26]
theorem c1a_le_frontier : C1a ≤ 1.5026709 := by
  obtain ⟨hsum, hconv⟩ := C1aCert.cert_facts
  have hpos : 0 < ∑ i ∈ Finset.range (2 ^ 15 * 8), C1aCert.dig C1aCert.P i := by rw [hsum]; norm_num
  have h := C1aCert.C1a_le_of_cert (2 ^ 15 * 8) (by positivity) (fun i => C1aCert.dig C1aCert.P i)
    C1aCert.T hconv hpos
  refine h.trans ?_
  rw [hsum]
  rw [div_le_iff₀ (by positivity)]
  norm_num [C1aCert.T]

-- TODO: Formalise relationship to Sidon sets.
-- TODO: State results from https://arxiv.org/abs/2602.07292

end Constant1a
