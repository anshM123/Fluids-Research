# Why the Boussinesq ladder is infinite: WKB quantisation in a quasi-stagnant boundary region

*Working notes. The numbers are updated as the scans progress; see README.md for the data.*

## 1. The branch and the smoothness function
For every λ > 1 the least-singular self-similar profile has the local exponent m(λ) = (λ−1)/ε(λ) at the
stagnation point, with ε = 1+λ−A and A = −∂₁U₁(0). Smooth profiles are the zeros of

  F(λ) := m(λ) − 2.

Along one connected branch, F oscillates about 0 as λ ↓ 1. The oscillation is regular in z = 1/(λ−1) and its
amplitude decays roughly geometrically. Each zero is a smooth profile: λ₀ (stable), λ₁, λ₂, … . The 1D Hou–Luo
model, u_x = Hω on the boundary, shows exactly the same structure (hl_solver.py).

## 2. Structure of the profiles as λ → 1
The data give three regions along the boundary (x = y₁, s = ln x). Here D(x) = V₁/x = 1+λ+U₁/x is the radial
speed of the self-similar flow along the boundary.

1. **Quasi-stagnant region** x < x_c: D = ε D̂(s) with D̂ = O(1). D̂ → 1 at the origin and has a dip of depth
   D̂_min ≈ 0.4–0.8 before the front. Θ_b ≈ −x² and Ω_b ≈ ∂ₓΘ_b (vorticity slaved to the temperature gradient).
2. **Front** at x_c ≈ 0.5–0.6 (fixed as λ → 1). Here D rises from O(ε) to O(1) over a width that shrinks like ε
   (Hou–Luo: max d ln D/dη ≈ 1.1/ε).
3. **Outer region** x > x_c: an O(1) flow, where Θ saturates (Θ grows only like x^{(λ−1)/(1+λ)}).

**The front becomes a square-root cusp as λ → 1.** On the outer side the radial speed approaches
D ≈ k (x − x_c)^{1/2} in both models:
- Hou–Luo: x_c ≈ 0.61 and k ≈ 1.1 (`hl_front.png`); 2D: x_c ≈ 0.72 and k ≈ 1.2 (`bq_front.png`).
- The curves collapse onto the square-root law outside a regularising inner layer that shrinks with ε.

This is the terminal structure of the CCF branch (P02, den ∝ |x−x_s|^{1/2}). Here, however, the stagnant
interval x < x_c can exist only when ε → 0, so the cusp is reached only in the limit λ → 1. A leading-order
reason: on the stagnant interval the velocity must equal −(1+λ)x to O(ε). That is a finite-Hilbert-transform
("airfoil") condition on the vorticity, whose generic solution has an inverse-square-root edge singularity.

**The dip does not close at finite ε.** The dip just inside the front deepens and moves towards x_c, but slowly:
- Hou–Luo: D̂_min = 0.58, 0.52, 0.47, 0.44, 0.41 at ε = 0.116, 0.085, 0.067, 0.055, 0.047, i.e. D̂_min ∝ ε^{0.38}.
- 2D: D̂_min = 0.83, 0.73, 0.665, 0.618 at ε = 0.20, 0.10, 0.06, 0.04, i.e. D̂_min ∝ ε^{0.18}.

So D_min = εD̂_min → 0 only as λ → 1. A sonic point, which is how the CCF ladder terminates, is not formed at
any λ > 1 on the computed range (scans continue to λ = 1.02 in 2D and λ = 1.01 in Hou–Luo).

**Why m → 2.** A least-singular profile with m ≠ 2 carries a non-analytic velocity component c(m) x^{m−1} along
the boundary, with c(m) ∝ (m−2). The particular solution of −ΔΨ = C y₁^{m−1} in the corner has
Ψ_β(β = 0) ∝ (m−2); it vanishes for m = 2, where Ψ = (C/2) y₁y₂². In the quasi-stagnant region, D − ε itself
must stay O(ε), so (m−2) = O(ε). The same holds at every order of the formal expansion in ε. The smoothness
defect is therefore beyond all orders in ε, and set by exponentially small (WKB) terms.

## 3. WKB modes of the quasi-stagnant region
In the quasi-stagnant region the radial transport speed is O(ε) while all other rates are O(1). Consider
perturbations ∝ exp(i∫K dx) with K = κ/(εx) and κ = O(1).

**Hou–Luo (closed form).** The Mellin symbol of Ω ↦ U/ξ behaves as −1/k for |k| → ∞. The linearised equations
D θ'_η + q'Θ_η = δθ', D ω'_η + q'Ω_η + ω' = e^{−η}θ'_η and q' = M ∗ ω' then give, at leading order,

  i D̂ κ² + κ − Γ = 0,   Γ = e^{−η}Θ_η/D̂ ≈ Ω/D̂   ⇒   κ = (−1 + √(1 + 4iΩ)) / (2iD̂).

The root is complex: the WKB modes oscillate on the scale ε and decay exponentially.

**2D Boussinesq (vertical eigenproblem).** Near a boundary point, use the vertical coordinate y = εx·Y, which is
comparable to the radial wavelength. The local base flow is V₁ = εx(D̂ + ĉY) with ĉ = −Ω_b, and V₂ = μy with
μ = 2(1+λ) − D − xD'. The perturbation equations at leading order are

  [iκ(D̂+ĉY) + μY∂_Y] θ̃ = −G ∂_Yψ̃ + iκΘ_y ψ̃,
  [1 + iκ(D̂+ĉY) + μY∂_Y] ω̃ = iκ θ̃,
  −∂_Y²ψ̃ + κ²ψ̃ = ω̃,   ψ̃(0) = 0,  ψ̃ → 0 as Y → ∞,

with G = ∂ₓΘ and Θ_y = ∂_yΘ at the boundary. This is a nonlinear eigenvalue problem for κ. It is solved by
shooting (the regular solution at Y = 0, with the growing mode e^{κY} suppressed); see `bq_local_eig.py`. Its root
continues the Hou–Luo root. For the same local data, Re κ is smaller and |Im κ| larger than in Hou–Luo.

## 4. Quantisation: the ladder law
The smoothness defect is carried by the WKB wave between the front and the stagnation point:

  F(λ) ≈ Re[ K(λ) e^{iΦ(λ)} ],   Φ = (1/ε)∫ κ(s) ds = m z ∫ κ ds   (1/ε = m z, z = 1/(λ−1)).

With εΦ → a (complex) as ε → 0:
- zeros are spaced asymptotically by **Δz = π/(m Re a)**, which gives the empirical law 1/(λ_n − 1) ≈ z_* + nΔz;
- the amplitude decays like e^{−m|Im a|z}, a factor e^{π|Im a|/Re a} per half period, times a slowly varying
  algebraic factor (fits suggest ∝ z^{−p} with p ≈ 1.5).

## 5. Numerical tests
The WKB phase is computed with the root tracked from the Hou–Luo closed form at every point, integrated from the
stagnation point to the front (first point beyond the dip with D/ε > 3).
- 2D: `bq_wkb2d.py`, evaluated on the refined smooth profiles.
- Hou–Luo: `hl_wkb2.py`, along the branch and interpolated to the crossings.

**2D Boussinesq.** Refined smooth profiles (Nb = 32, hs = 0.025):

| n | λ_n | z_n = 1/(λ_n−1) | spacing | εΦ (to the front) | Re Φ = 2z εΦ | ΔRe Φ |
|---|---|---|---|---|---|---|
| 0 | 1.9205593 | 1.08630 | — | 1.2271 − 0.4984i | 2.666 | — |
| 1 | 1.3990960 | 2.50566 | 1.41937 | 1.0793 − 0.5443i | 5.409 | 2.743 |
| 2 | 1.2523486 | 3.96278 | 1.45711 | 1.0586 − 0.5612i | 8.390 | 2.981 |
| 3 | 1.1842533 | 5.42733 | 1.46455 | 1.0521 − 0.5722i | 11.420 | 3.030 |
| 4 | 1.1449864 | 6.89718 | 1.46985 | 1.0505 − 0.5809i | 14.492 | 3.072 |
| 5 | 1.1194818 | 8.37040 | 1.47322 | 1.0516 − 0.5883i | 17.602 | 3.110 |
| 6 | 1.1015235 | 9.84999 | 1.47959 | | | |

- The phase gained between consecutive smooth profiles, 2.98, 3.03, 3.07, 3.11, tends to π. This is the
  quantisation condition Re Φ(λ_n) = Φ₀ + nπ.
- εΦ levels off at a ≈ 1.05 − 0.6i, so the asymptotic spacing is Δz = π/(2 Re a) ≈ 1.50. The observed spacings
  (1.457 → 1.480) increase towards it.
- The two-point law of Wang et al. (slope 1.4187, the line through λ₀ and λ₁) therefore underestimates the
  asymptotic slope by about 5 %.

**Hou–Luo.**
- Crossings at z = 4.742, 6.000, 7.261, 8.525, 9.790, 11.054, 12.333, 13.56. ΔRe Φ = 3.025, 3.052, 3.074, 3.069,
  3.092, 3.094, tending to π.
- εΦ_cut increases towards the same limit as the full integral, a ≈ 1.24 − 0.50i (1.2006 − 0.4820i at z = 23),
  so Δz_∞ = π/(2 Re a) ≈ 1.27. Observed spacings 1.2538 → 1.2674.

**Amplitude.** The per-interval WKB factor e^{|ΔIm Φ|} is ≈ 6 in 2D (|ΔIm Φ| = 1.65 → 1.84) and ≈ 3.5 in Hou–Luo
(≈ 1.24). The observed extremum ratios (2D: 10.9, 9.4, 8.7, …; HL: 6.7 → 4.55) are larger. After dividing the
extrema by the computed WKB factor e^{Im Φ(z)}, what remains is a clean power law:
- Hou–Luo, 8 extrema (z = 5.2–14.1, |m−2| from 5.5e-4 to 1.8e-8): |m − 2| e^{−Im Φ} ∝ z^{−p}, with p = 1.58 overall
  (local values 1.48–1.76).
- 2D, extrema at z = 1.5–7.4: local p = 1.01, 1.23, 1.33, 1.29, not yet asymptotic.

Hence, empirically, m(λ) − 2 ≈ C (λ−1)^{p} Re[e^{iΦ(λ) + iφ}] with p ≈ 3/2.

**Consequence.** Both models have an infinite ladder of smooth self-similar profiles accumulating at λ = 1, with
1/(λ_n − 1) growing linearly in n. The empirical law of Wang et al. is the leading-order WKB quantisation of the
quasi-stagnant region. The asymptotic 2D slope is π/(2 Re a) ≈ 1.50, larger than the slope 1.4187 of the line through λ₀ and λ₁.
