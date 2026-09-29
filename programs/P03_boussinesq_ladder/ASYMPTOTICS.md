# The λ → 1⁺ limit of the Boussinesq branch: what is proved, what is derived, what is only measured

*Status of every statement: **[P]** proved (elementary), **[F]** formal matched asymptotics (derived, not
rigorous), **[N]** numerical measurement with its error bar, **[O]** open. Nothing here is a theorem about the
Boussinesq equations.*

## 0. Summary
- **[P] Lemma.** Write the smoothness defect along the branch as F(z) = R(z)[cos Θ(z) + η(z)]. Suppose R > 0, Θ is
  increasing and unbounded, and |η| < 1. Then F has infinitely many zeros, at least one per increase of Θ by π. If,
  in addition, Θ(z) = Cz + Θ₀ + o(1) and η → 0, then the zeros satisfy λ_n = 1 + C/(nπ + c + o(1)).
- **[F] WKB representation.** Matched asymptotics of the quasi-stagnant boundary layer give F in this form:
  - Θ = Re Φ + arg K, where Φ(λ) = ε⁻¹∫κ ds is the phase of the leading WKB wave;
  - κ is the continued root of the local (boundary-layer) eigenproblem;
  - C = 2 Re a₀ with a₀ = ∫κ₀ ds;
  - η = O(ε) = O(1/z).

  The O(1) remainder Θ₀ + o(1) requires the stagnant-layer profile to have a regular expansion in ε up to the
  front. This hypothesis (U) is **not** satisfied uniformly near the front, where the dip closes like ε^{0.18}.
  Without (U) the derivation gives only Θ = Cz + o(z). That is still enough for infinitely many zeros and for the
  limiting spacing π/C.
- **[P given F] C > 0.** In the Hou–Luo model, Re κ₀ > 0 follows in closed form from Ω > 0 in the layer. In 2D it is
  checked pointwise along every computed profile.
- **[N] Measurements.**
  - Re Φ(λ_n) − nπ = δ_n converges (2.27, 2.11, 2.00, 1.92, 1.85, 1.83, 1.81 for n = 1…7; increments shrink like
    z⁻²).
  - Re Φ(z) = Cz + Φ₀ + E with |E| ≤ 0.03 on z ∈ [4, 11.3].
  - C ∈ [2.09, 2.13], i.e. an asymptotic spacing Δz_∞ = π/C ∈ [1.476, 1.50]. The central WKB value is 1.496.
- **[O] Open.**
  - The exact value of C. Δz_∞ = 3/2 (C = 2π/3) is a candidate consistent with the data. It is **not derived**.
  - Uniqueness of the zero in each half period (observed, not shown).
  - A rigorous version of any [F] step.

## 1. Setting
- **Branch.** For each λ > 1 the least-singular self-similar profile has Θ ≈ −|y₁|^m at the stagnation point.
  - m = (λ−1)/ε, with ε = 1+λ−A and A = −∂₁U₁(0).
  - Write z = 1/(λ−1). Then 1/ε = m z exactly.
  - The smoothness defect is F(z) = m − 2; smooth profiles are its zeros.
  - Along the computed branch |F| ≤ 7.4e-2 for z ≥ 1.1, so 1/ε = 2z(1 + O(F)).
- **Radial transport speed.** On the wall (x = y₁ = e^s) it is D(s) = V₁/x = (1+λ) + U₁/x, with D(−∞) = ε.

## 2. An elementary lemma [P]
**Lemma.** Let F ∈ C([z₀, ∞)) satisfy F(z) = R(z)[cos Θ(z) + η(z)], with:
- R, Θ continuous and R > 0;
- Θ strictly increasing with Θ → ∞;
- sup|η| = η̄ < 1.

Then:
1. For every integer k with kπ ≥ Θ(z₀), F(Θ⁻¹(kπ)) and F(Θ⁻¹((k+1)π)) have opposite signs, (−1)^k and (−1)^{k+1}.
   So F has a zero in every interval Θ⁻¹([kπ, (k+1)π]), infinitely many in total, and they accumulate only at z = ∞.
2. If η(z) → 0, every zero z* satisfies Θ(z*) = (k+½)π + O(|η(z*)|) for some integer k.
3. If in addition Θ, η ∈ C¹ with Θ' ≥ c > 0 and η' → 0, each such interval contains exactly one zero for large z,
   and the zero is simple. This is because (cos Θ + η)' = −Θ' sin Θ + η' ≠ 0 where |sin Θ| ≈ 1.
4. If Θ(z) = Cz + Θ₀ + o(1), the zeros, numbered consecutively, satisfy z_n = (nπ + c')/C + o(1). Equivalently
   λ_n = 1 + C/(nπ + c + o(1)).
5. If only Θ = Cz + o(z) is known, then still z_n ~ nπ/C and z_{n+1} − z_n → π/C, provided Θ'(z) → C.

*Proof.* Items 1 and 2 follow from the intermediate value theorem and |cos Θ| = |η| at a zero. Item 3 is the
implicit function theorem. Items 4 and 5 follow by inverting Θ. ∎

In the reviewer's notation, Φ(λ_n) = nπ + δ + o(1) is item 2 with Φ = Θ − arg K, and Φ(λ) = C/(λ−1) + Φ₀ + o(1) is
the hypothesis of item 4.

## 3. The WKB representation of F [F]
### 3.1 Regions (verified on the computed profiles, THEORY.md §2)
| region | extent along the wall | D | numerical check |
|---|---|---|---|
| I (inner) | x ≲ ε (s ≲ ln ε) | ε(1 + O(x)) | exact local expansion; s_start −12 → −20 changes λ₁ by 2e-6 |
| II (quasi-stagnant) | ε ≪ x < x_c | εD̂(s), D̂ = O(1) | D̂ ∈ [D̂_min, 1] for z = 2.5–11.3 |
| III (front) | x − x_c = O(ε) in s | rises from O(ε) to O(1) | max d ln D/ds ∝ 1/ε (Hou–Luo 1.1/ε) |
| IV (outer) | x > x_c | O(1), D ≈ k(x−x_c)^{1/2} near the front | x_c ≈ 0.72, k ≈ 1.2 (2D) |

The dip of D̂ just inside the front deepens like D̂_min ∝ ε^{0.18} in 2D and ε^{0.35–0.38} in Hou–Luo. Region II
therefore has a regular ε-expansion only away from the front. This is the non-uniformity mentioned in §0.

### 3.2 Why F is exponentially small
Suppose m ≠ 2. The leading vorticity near the origin ∝ r^{m−1} then induces a wall velocity ∝ (m−2)x^{m−1}. For
m = 2 the induced wall velocity vanishes (Ψ = (C/2)y₁y₂²). In region II the total wall velocity must equal
−(1+λ)x up to O(ε), so m − 2 = O(ε). A formal expansion of the profile in powers of ε in region II repeats this
argument at every order. Hence the formal series is that of a smooth profile at every λ, and F is beyond all orders.

The defect is carried by an exponentially small solution of the linearised steady (μ = 0) self-similar operator
around the formal profile. In region II such solutions have WKB form.

### 3.3 Eikonal equation: the local eigenproblem
**Ansatz.** Near a wall point (x, 0) of region II use the vertical variable Y = y/(εx). Take
θ' = θ̃(Y) e^{iϕ}, ω' = ω̃(Y) e^{iϕ}/(εx), ψ' = εx ψ̃(Y) e^{iϕ}, with ϕ = ε⁻¹∫^s κ ds.

**Local base flow.** Expanded in Y:
- V₁ = εx(D̂ + ĉY) + O(ε²Y²), with ĉ = −Ω_b;
- V₂ = μ_v y + O(y²), with μ_v = 2(1+λ) − D − xD';
- Θ_x = G and Θ_y taken at the wall.

**Leading order O(1)** (this is `bq_local_eig.py`):
```
[iκ(D̂+ĉY) + μ_v Y∂_Y] θ̃ = −G ∂_Yψ̃ + iκ Θ_y ψ̃
[1 + iκ(D̂+ĉY) + μ_v Y∂_Y] ω̃ = iκ θ̃
−∂_Y²ψ̃ + κ²ψ̃ = ω̃,   ψ̃(0) = 0,  ψ̃ bounded (the growing mode e^{κY} excluded)
```

**Terms dropped are all O(ε) relative:**
- (1−λ)θ' = −mεθ';
- u'·∇Ω in the vorticity equation;
- the Y² terms of V₁;
- the slow x-variation of the amplitude and of the base data over one wavelength (relative size ε/κ);
- the log-polar metric.

So κ = κ₀(s) + εκ₁(s) + O(ε²) wherever the base data have a regular expansion and κ is bounded away from 0.

**Hou–Luo analogue** (1D, Mellin symbol ~ −1/k):
- The eigenproblem reduces to i D̂ κ² + κ − Ω = 0.
- The root is κ_HL = (−1 + √(1+4iΩ))/(2iD̂).

### 3.4 Phase and amplitude
Solvability at O(ε) gives a transport equation for the amplitude A(s), and the full phase is
Φ = ε⁻¹∫κ₀ ds + ∫κ₁ ds + arg A + O(ε).

The WKB form fails in two places, each giving an O(1) connection phase:
1. **Region I.** As s → −∞, κ₀ → 0 like e^{s} (Γ = e^{−s}Θ_s/D̂ ∝ x). The approximation needs ε/κ ≪ 1, so it fails
   for x ≲ ε. There the wave matches the linearised stagnation-point problem, which converts the wave amplitude
   into a strain perturbation δA and hence F = (2/ε)δA. The phase integral converges at −∞, and this region adds
   an O(1) constant plus O(ε).
2. **Region III.** κ → 0 as D̂ → ∞. The front has width O(ε) in s and contributes an O(1) constant, provided the
   front structure has a regular inner limit. Because the dip closes (D̂_min → 0), this contribution is only
   known to be o(1/ε). The phase error there is E(z) = O(z^{q}) with an unknown q < 1, possibly q = 0.

**Result [F].**
- F(z) = |K| e^{−Im Φ}[cos(Re Φ + arg K) + O(ε)].
- Re Φ(z) = C z + Φ₀ + E(z), with C = 2 Re a₀ and a₀ = ∫_{−∞}^{s_c} κ₀ ds.
- E = o(1) under (U); E = o(z) in general.
- With the Lemma this gives infinitely many zeros and the spacing π/C. Under (U) it also gives λ_n = 1 + C/(nπ + c + o(1)).

## 4. Positivity of C [P given the WKB form; N in 2D]
**Hou–Luo.** Write √(1+4iΩ) = p + iq (principal branch).
- p² − q² = 1 and 2pq = 4Ω, so q = 2Ω/p > 0 and p > 1 when Ω > 0.
- κ_HL = [q − i(p−1)]/(2D̂), so Re κ > 0 and Im κ < 0 wherever Ω > 0 and D̂ > 0.
- Ω > 0 holds throughout the computed layers, so Re a₀ > 0 and C > 0. The sign Im κ < 0 is the exponential decay
  of the wave from the front towards the stagnation point, i.e. the beyond-all-orders smallness of F.

**2D.** The root is continued from κ_HL at every point. Along all computed profiles min Re κ > 0 and max Im κ < 0
(`wkb_sens.py`, table in §5).

## 5. Numerical verification [N]
### 5.1 Phase at the smooth profiles
- WKB phase cut at the front (D/ε > 3); λ₄–λ₇ at hs = 0.0125, the others at hs = 0.025.
- δ_n = Re Φ(λ_n) − nπ.

| n | z_n | Re Φ(λ_n) | δ_n | Re εΦ |
|---|---|---|---|---|
| 1 | 2.50566 | 5.409 | 2.267 | 1.0793 |
| 2 | 3.96277 | 8.390 | 2.107 | 1.0586 |
| 3 | 5.42733 | 11.420 | 1.995 | 1.0521 |
| 4 | 6.89723 | 14.490 | 1.924 | 1.0505 |
| 5 | 8.37004 | 17.562 | 1.854 | 1.0491 |
| 6 | 9.84429 | 20.677 | 1.827 | 1.0502 |
| 7 | 11.32010 | 23.797 | 1.806 | 1.0511 |

**Quantisation.**
- The δ_n decrease by 0.16, 0.11, 0.071, 0.070, 0.027, 0.021, roughly like z⁻².
- A fit δ_n = δ + d/z (n = 5–7) gives δ = 1.67. A fit with d/z² gives 1.74.
- So Φ(λ_n) = nπ + δ + o(1) holds with δ = 1.7 ± 0.05 (this δ refers to the cut integral; the cut-off moves Φ₀ and
  δ together).

**Linear law.** Least squares for Re Φ = Cz + Φ₀:
| data | C | Φ₀ | max residual | π/C |
|---|---|---|---|---|
| n = 2–7 | 2.094 | 0.06 | 0.029 | 1.500 |
| n = 3–7 | 2.100 | 0.007 | 0.024 | 1.496 |
| n = 4–7 | 2.105 | −0.04 | 0.016 | 1.492 |

Adding a 1/z term gives C = 2.12–2.15. That fit is less stable (four points, three parameters) and is not used for
the central value.

**Directly measured spacings.**
- z_{n+1} − z_n = 1.4646, 1.4699, 1.4728, 1.4743, 1.4758 (n = 3…7), increasing monotonically.
- The local slopes π/Δz_n = 2.1451, 2.1373, 2.1331, 2.1310, 2.1287 decrease monotonically, which gives C < 2.129.
- The λ₇ entry carries an uncertainty of about 5e-3 in z from the hs-sensitivity of the deepest crossing.

**Estimate.** Combining the WKB slope and the monotone bound, C ∈ [2.09, 2.13], i.e. Δz_∞ ∈ [1.476, 1.50].

### 5.2 Hou–Luo (sharper test of the same statements)
- There are 12 crossings (z up to 15).
- ΔRe Φ per rung is 3.02–3.09, tending to π.
- The observed spacings rise from 1.2538 to 1.2674, against the WKB limit π/(2 Re a₀) ≈ 1.27 with a₀ ≈ 1.24 − 0.50i.
  The deepest spacing is within 0.3 % of the prediction.

## 6. The value 3/2 [O]
- Δz_∞ = 3/2 is equivalent to Re a₀ = π/3 = 1.0472. The measured Re εΦ_cut is 1.049–1.051 for n = 4–7.
- Two things make the difference hard to settle:
  - the cut-off contribution (§5.3, `wkb_sens.py`);
  - the unknown sign of the finite-ε drift. In Hou–Luo, εΦ_cut approaches a₀ from below.
- So Re a₀ = π/3 is neither excluded nor supported beyond this. Nothing in the leading-order problem suggests a
  closed form for a₀: it is an integral of a numerically determined local root over a numerically determined layer
  profile.
- Recommended wording: "λ_n = 1 + C/(nπ + c + o(1)) with C = 2.10 ± 0.02, i.e. an asymptotic spacing
  1/(λ_{n+1}−1) − 1/(λ_n−1) → π/C ≈ 1.49 ± 0.01; the value 3/2 is a candidate."
