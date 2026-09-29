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
  Without (U) the derivation gives only Θ = Cz + o(z); the measured dip scaling indicates Θ = Cz + O(z^{0.6}).
  That is still enough for infinitely many zeros and for the limiting spacing π/C, but convergence to that
  spacing is slow.
- **[P given F] C > 0.** In the Hou–Luo model, Re κ₀ > 0 follows in closed form from Ω > 0 in the layer. In 2D it is
  checked pointwise along every computed profile.
- **[N] Measurements.**
  - **Quantization.** The WKB phase integrated through most of the front (cut at D/ε = 5 or 8) satisfies
    Re Φ(λ_n) − nπ = δ + O(0.03) for n = 4–7, with δ ≈ 2.0 (cut 5) or 2.3 (cut 8). Cut at D/ε = 2–3, δ_n still
    drifts by −0.04 per rung.
  - **Linear law.** Re Φ(z) = Cz + Φ₀ + E with |E| ≤ 0.024 on z ∈ [6.9, 11.3] for every cut-off.
  - **Asymptotic spacing.** π/C lies between 1.476 and about 1.51. Its value depends on the form of the corrections
    (§5.3, §5.5).
    - The measured spacings increase monotonically; if that trend continues, 1.4758 (n = 7) is a lower bound.
    - Analytic 1/z phase corrections extrapolate to 1.478.
    - The measured scaling of the dip in front of the front suggests a non-analytic phase correction ∝ z^{0.6}.
      That extrapolates to 1.49–1.50.
    - All these models fit the seven rungs to ≤ 1e-3; the data cannot tell them apart.
- **[O] Open.**
  - The exact value of C.
  - Δz_∞ = 3/2 (C = 2π/3) is a candidate. It is consistent with the non-analytic extrapolation suggested by the dip
    scaling, larger than the analytic-correction extrapolation (1.478), and **not derived**.
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
- The eigenproblem reduces to i D̂²κ² + D̂κ − Ω = 0, i.e. v(1 + v) = iΩ with v = iD̂κ.
- The root is κ_HL = (−1 + √(1+4iΩ))/(2iD̂).

### 3.4 Phase and amplitude
Solvability at O(ε) gives a transport equation for the amplitude A(s), and the full phase is
Φ = ε⁻¹∫κ₀ ds + ∫κ₁ ds + arg A + O(ε).

The WKB form fails in two places, each giving an O(1) connection phase:
1. **Region I.** As s → −∞, κ₀ → 0 like e^{s} (Γ = e^{−s}Θ_s/D̂ ∝ x). The approximation needs ε/κ ≪ 1, so it fails
   for x ≲ ε. There the wave matches the linearised stagnation-point problem, which converts the wave amplitude
   into a strain perturbation δA and hence F = (2/ε)δA. The phase integral converges at −∞, and this region adds
   an O(1) constant plus O(ε).
2. **Region III.** κ → 0 as D̂ → ∞. The front has width O(ε) in s and would contribute an O(1) constant if the
   front structure had a regular inner limit. It does not: the dip just inside the front closes (D̂_min ∝ ε^{0.18})
   over a width ∝ ε^{0.69}, and the local root there grows like D̂^{−3/2} (§5.5). The scaling estimate for its phase
   contribution is E(z) = O(z^{q}) with q ≈ 0.58 < 1. So the formal result is Θ = Cz + O(z^{0.6}), not
   Cz + Θ₀ + o(1).

**Result [F].**
- F(z) = |K| e^{−Im Φ}[cos(Re Φ + arg K) + O(ε)].
- Re Φ(z) = C z + Φ₀ + E(z), with C = 2 Re a₀ and a₀ = ∫_{−∞}^{s_c} κ₀ ds.
- E = o(1) under (U); E = o(z) in general. The dip scaling indicates E ~ z^{0.58} (§5.5).
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
### 5.1 Phase at the smooth profiles (front cut-off D/ε = 3)
- λ₄–λ₇ at hs = 0.0125, the others at hs = 0.025.
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

### 5.2 Dependence on the front cut-off (`wkb_sens.py`, `wkb_sens.out`)
The phase integral stops where D/ε first exceeds a cut-off after the dip. The table fits Re Φ = Cz + Φ₀ to
n = 4–7 (hs = 0.0125).

| cut D/ε | C | π/C | max residual | ΔRe Φ (4→5, 5→6, 6→7) | δ_n (n = 4…7) |
|---|---|---|---|---|---|
| 2 | 2.1039 | 1.4932 | 0.001 | 3.096, 3.102, 3.107 | 1.78, 1.73, 1.69, 1.66 |
| 3 | 2.1051 | 1.4924 | 0.016 | 3.071, 3.115, 3.121 | 1.92, 1.85, 1.83, 1.81 |
| 5 | 2.1165 | 1.4843 | 0.024 | 3.079, 3.153, 3.118 | 2.05, 1.99, 2.00, 1.98 |
| 8 | 2.1264 | 1.4775 | 0.022 | 3.098, 3.128, 3.181 | 2.35, 2.30, 2.29, 2.33 |

**Resolution checks.**
- Halving the quadrature step changes εΦ by ≤ 4e-4.
- Halving hs of the profile changes εΦ by 8e-5 at cut 3, but by 3e-3 at cut 5 (the front is hs-sensitive).
- Root-tracking failures occur only for s < −6.7, where κ < 0.005. They contribute 0.0017 to εΦ.

**Why the front matters.** The region between D/ε = 2 and 8 has width 0.41, 0.33, 0.26, 0.225 in s at
z = 6.9, 8.4, 9.8, 11.3, i.e. ≈ 2.6/z ∝ ε. So its phase contribution is O(1), consistent with the formal
derivation. It still changes with z in the computed range, which biases the finite-z slopes. Including more of the
front has two effects:
- it makes δ_n constant (the quantization is then satisfied to ±0.03);
- it lowers π/C towards the directly measured spacing.

### 5.3 Direct fits of the resonance positions
The fits are of nπ = C z_n + a (+ b/z_n or b/z_n²) to the exact zeros.

| data | linear | + b/z | + b/z² |
|---|---|---|---|
| n = 2–7 | 1.4717 | 1.4780 | 1.4759 |
| n = 3–7 | 1.4733 | 1.4784 | 1.4767 |
| n = 4–7 | 1.4743 | 1.4788 | 1.4773 |

**Uncertainty from λ₇.** λ₇ is uncertain by about 1e-5; the two solvers differ by 9.6e-6, which gives δz₇ ≈ 1.3e-3.
That moves the b/z value by about ±0.002. A 5e-5 uncertainty would move it by ±0.006.

**Directly measured spacings.**
- z_{n+1} − z_n = 1.4646, 1.4699, 1.4728, 1.4743, 1.4758 (n = 3…7): monotone, with increments 5.3e-3, 2.9e-3, 1.5e-3,
  1.5e-3.
- A power-law extrapolation of the increments (exponent 2.0–2.3) gives a limit of 1.483–1.485.

**Extrapolation depends on the correction exponent.** Fit the local phase slope π/Δz_n = C + B z^{−β} (interval
midpoints z = 4.7–10.6):

| β | meaning | π/C (n ≥ 2) | π/C (n ≥ 3) | max residual |
|---|---|---|---|---|
| 2 | analytic 1/z phase corrections | 1.4781 | 1.4785 | 4e-4 |
| 1 | log z phase correction | 1.4850 | 1.4839 | 3e-4 |
| 0.75 | | 1.4896 | 1.4875 | 4e-4 |
| 0.5 | phase correction ∝ z^{1/2} | 1.4989 | 1.4949 | 4e-4 |
| 0.42 | dip scaling, §5.5 | 1.5043 | 1.4991 | 5e-4 |
| 0.3 | | 1.5180 | 1.5097 | 5e-4 |

Fitting the phase itself, nπ = Cz + a + b z^{1−β}, gives the same: 1.4963–1.4986 for β = 0.42, 1.4834–1.4837 for
β = 1.

**Estimate.** π/C lies between 1.476 (the monotone lower bound) and about 1.51. The analytic-correction value is
1.478, and the value suggested by the dip scaling is 1.49–1.50.

### 5.4 Where hypothesis (U) holds (`wkb_uniform.py`, `wkb_uniform.out`)
The base data D̂(s) and the local root κ(s) are compared at fixed s on the five crossing profiles n = 3–7
(ε = 0.092 … 0.044), fitting a + bε.

| s (x = e^s) | D̂: ε → 0 limit, max residual | κ: ε → 0 limit, max residual |
|---|---|---|
| −3 (0.05) | 0.9995, 1.3e-4 | 0.0972 − 0.0122i, 5e-5 |
| −2 (0.14) | 0.9949, 2.4e-4 | 0.2380 − 0.0724i, 2.4e-4 |
| −1.5 (0.22) | 0.9817, 1.1e-3 | 0.3514 − 0.1568i, 5.5e-4 |
| −1 (0.37) | 0.9545, 3.6e-3 | 0.4948 − 0.3166i, 3.4e-3 |
| −0.8 (0.45) | 0.971, 1e-2 (slope −1.4) | 0.51 − 0.37i, 1.4e-2 |
| −0.6 (0.55) | 0.85, 3e-2 | 0.69 − 0.56i, 6e-2 |
| −0.5 (0.61) | 0.23 (dip closing), 4e-2 | 1.6 − 1.9i, 6e-2 |

**Result.**
- Hypothesis (U) holds, with O(ε) corrections, in the inner part of the layer (x ≲ 0.4).
- It fails in the dip region x ≈ 0.45–0.72, where the base flow still changes at O(1) as ε decreases: the dip deepens
  like ε^{0.18}.
- That region carries about a third of Re a₀. Its contribution to Φ need not have an analytic expansion in ε. This is
  the reason the remainder E(z) is only known to be o(z), and the reason for the upper widening of the range for π/C.
- This is the main open point of the asymptotic analysis. It needs an inner problem for the dip/front region as
  ε → 0, which is not attempted here.

### 5.5 The dip: a source of non-analytic corrections (`dip_scaling.py`)
Measured on the crossing profiles n = 3–7:

| quantity | fit on the four deepest |
|---|---|
| distance from the dip minimum to the front | ∝ ε^{1.00} |
| dip depth D̂_min | ∝ ε^{0.18} |
| dip half-width | ∝ ε^{0.69} |
| local root in the dip | ∝ D̂^{−1.56} (the Hou–Luo scaling κ ~ D̂^{−3/2}) |

So the finite-ε dip changes the phase by roughly ε⁻¹·w·κ_dip ~ ε^{−1+0.69−0.27} ≈ z^{0.58}. That is a
non-analytic correction to Φ = Cz + …, i.e. β ≈ 0.42 in §5.3.

This is a scaling estimate, not a matched inner solution. It shows that analytic extrapolations may underestimate
π/C.

### 5.5b Hou–Luo to ε = 0.012 (z ≈ 40): the phase coefficient converges (`hl_dip_scaling.py`, `hl_wkb2.out`)
**Dip scaling.** Over ε < 0.03 (ten states):
- depth D̂_min ∝ ε^{0.31};
- half-width ∝ ε^{0.34};
- distance from dip to front ∝ ε^{1.23};
- the dip and the front converge to x ≈ 0.61.

At moderate ε (0.04–0.09) the same quantities give different local exponents; the half-width, for example, scales
like ε^{0.96} there. The 2D exponents of §5.5, measured only at ε ≥ 0.044, are therefore pre-asymptotic.

**Phase coefficient.** Re εΦ_cut rises monotonically from 1.125 (ε = 0.116) to 1.211 (ε = 0.012). Fits over z > 12:

| model | a₀ | max residual |
|---|---|---|
| a + bε^{2/3} | 1.236 | 7.2e-4 |
| a + bε | 1.223 | 7.7e-4 |
| a + bε^{−1/6} (divergent) | — | 2.3e-3 |

So the data favour convergence, a₀ = 1.22–1.25. In Hou–Luo this means Φ = C/(λ−1) + O((λ−1)^{−1/3}): the leading
term dominates, and the correction is non-analytic with exponent 1/3 (β = 2/3 in §5.3).

**Spacing.** π/(2a₀) = 1.26–1.29, against measured spacings of 1.2583–1.2652 (z ≤ 13.6).

**Consequence for 2D.** If the 2D correction has the Hou–Luo exponent (β = 2/3), the §5.3 extrapolation gives
π/C ≈ 1.49. The 2D range 1.476–1.51 is kept.

### 5.6 Hou–Luo (sharper test of the same statements)
- There are 11 crossings resolved above the noise (z ≤ 13.56).
- ΔRe Φ per rung is 3.02–3.09, tending to π.
- The spacings are 1.2583, 1.2612, 1.2637, 1.2652, 1.2638, against the WKB limit π/(2 Re a₀) ≈ 1.267 with
  a₀ ≈ 1.24 − 0.50i (0.2 % agreement).
- Re κ > 0 and Im κ ≤ 0 hold at every point of the layer for z = 10.6–40.6, because Ω > 0 there.

## 6. The value 3/2 [O]
- Δz_∞ = 3/2 is equivalent to C = 2π/3 = 2.0944.
- The data (seven rungs) cannot fix the correction exponent β:
  - analytic corrections give 1.478;
  - the non-analytic corrections suggested by the dip scaling (β ≈ 0.4–0.5) give 1.49–1.50, i.e. 3/2 within the
    uncertainty.
- Deciding it needs either the inner (dip/front) problem as ε → 0, or several more rungs at higher resolution.
  The amplitude |m−2| ~ 1e-8 at n = 8 makes the latter expensive but feasible.
- Nothing in the leading-order problem suggests a closed form for a₀.
- Recommended wording: "λ_n = 1 + C/(nπ + c + o(1)); the asymptotic spacing π/C lies between 1.476 and ≈ 1.51, and
  3/2 is a candidate consistent with the non-analytic extrapolation but not derived."

## 7. The unstable spectrum: one eigen-condition from the same phase [F; checks N]
The numerical facts are in INSTABILITY_LADDER.md §3 and §6:
- the lower unstable eigenvalues form a ladder μ ≈ (λ − 1)(k + c);
- along the continuous branch the ladder moves up by one spacing per rung interval;
- a new member enters at the bottom once per interval.

What follows is a formal account in the language of §3.

### 7.1 The WKB root does not depend on the growth rate
With perturbations ∝ e^{μτ} the transport symbols of §3.3 become μ + iD̂κ (temperature) and 1 + μ + iD̂κ
(vorticity). All other dropped terms stay O(ε).
- **Hou–Luo.** The local relation becomes (μ + iD̂κ)(1 + μ + iD̂κ) = iΩ. It depends on κ and μ only through
  v = μ + iD̂κ. Hence

      κ(μ) = κ(0) + iμ/D̂   for every μ,

  so Re κ, and with it Re Φ = ε⁻¹ Re∫κ ds, is independent of μ.
- **2D.** The local eigenproblem of §3.3 with μ added was solved by shooting (sample wall points x = 0.1–0.5,
  with and without Θ_y). The result is κ(μ) − κ(0) = iμ/D̂ to 1e-10 for μ = 0.02–0.3.
- **Meaning.** The imaginary shift is the WKB form of the exact identity L_μ(Θ̄φ) = Θ̄ L_{μ+λ−1}φ:
  exp(−μ∫ds/(εD̂)) = |Θ̄|^{−μ/(λ−1)}. The growth rate changes only the real amplitude of the waves, not their phase.

### 7.2 The stagnation point supplies a phase linear in μ
**Region I (x ≲ ε), Hou–Luo.**
- Use X = x/ε, ν = μ/ε, and Θ̄ = x² + …; for m = 2, ε = (λ−1)/2.
- To leading order ω' = θ'_X/ε and u' = H(θ') + (outer strain), with H the Hilbert transform in X.
- The temperature equation becomes

      (ν − 2)θ' + Xθ'_X = −2X H(θ') + (strain forcing).

**Fourier solution** (k conjugate to X, H ↔ −i sgn k).
- For k ≠ 0 the homogeneous equation gives (ν − 3)θ̂ = (k − 2 sgn k)θ̂_k. So θ̂ ∝ |k₀ − |k||^{ν−3}, with k₀ = 2 the
  wavenumber of the WKB wave at the stagnation point (κ ≈ Γ ∝ x in §3.4).
- Regularity at the stagnation point forbids a power-law tail at large |k|. The regular solution is therefore
  supported on |k| < k₀ (plus terms at k = 0 from the strain).
- Its endpoint singularity at k → k₀⁻ gives, at large X,

      θ' ∼ A Γ(ν − 2) X^{2−ν} cos(k₀X − π(ν − 2)/2),

  a standing wave with phase −π(ν − 2)/2 = −π(μ/(λ−1) − 1).
- This is the phase of the Mellin transform of the wave, ∫X^{s−1}e^{ik₀X}dX = Γ(s)(−ik₀)^{−s}. It depends only on
  the transport near the stagnation point being the Euler operator (ν − 2) + X∂_X. The same phase is therefore
  expected in 2D, where the wall transport near the stagnation point has the same form (the 2D region-I problem
  itself is not solved here).

### 7.3 Eigen-condition
**The structure of an eigenfunction.**
- It consists of a smooth part: θ' ≈ Θ̄·const in the layer, and localized at the dip/front (INSTABILITY_LADDER.md §4).
- On top of it sits a WKB standing wave, generated at the front and decaying toward the stagnation point.
- The smooth part is regular at x = 0. The wave's continuation through region I must be the regular solution of §7.2.

**Matching.** Equating the phase of the wave arriving from the front with the phase of that regular solution gives

    Re Φ(λ) − π μ/(λ−1) + θ₀ = kπ + o(1),   k ∈ ℤ,

where θ₀ collects the O(1) connection phases at the front and at the stagnation point.
- Re Φ(λ) is the phase that quantizes the profiles (§3.4, Re Φ(λ_n) = nπ + δ).
- The sign of the μ-term is fixed by §7.2. It predicts that the ladder moves up as λ decreases, which is what is
  observed.

### 7.4 Consequences and checks
| consequence | check |
|---|---|
| spacing of the lower eigenvalues → λ − 1 | extrapolated spacing 1.00–1.02 (HL; deep states running), 1.003 (2D) (`fig_spectral_flow.png` d) |
| the offset moves with Re Φ: +1 spacing per rung interval, linear in z between rungs; new members enter at μ = 0 | HL: 9 branch points, within 0.011; entering members at their predicted positions (INSTABILITY_LADDER.md §6) |
| one new unstable eigenvalue per half-turn of Re Φ, so index(n) = n + const, with const fixed by one rung | counts n for n = 0–7 (2D), 0–10 (HL) |
| the uniform ladder has about Re Φ/π members and ends near C/π = 1/Δz_∞ | computed uniform-ladder tops 0.73–0.74 (HL, n = 8–10) and 0.65 (2D, n = 7), against C/π ≈ 0.79 and 0.66–0.68. The highest eigenvalue is displaced a further 0.06–0.10 upward, which is outside this leading-order account |

### 7.5 What remains open
- The constants θ₀ and hence the offsets c: 0.65–0.67 in Hou–Luo and 0.70–0.72 in 2D, in units of the local
  spacing. They need the front connection and the region-I problem with the strain forcing.
- The upper end, where μ = O(1). There the top member lies 0.06–0.10 above the uniform ladder, and in Hou–Luo
  adjacent members collide into weakly complex pairs.
