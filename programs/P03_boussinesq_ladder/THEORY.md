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
**Hou–Luo** (εΦ from `hl_wkb.py`, extrapolated to ε → 0):
- a ≈ 1.238 − 0.506i, which predicts Δz = 1.269. Observed spacings 1.2634–1.2653 at z = 5–8 (0.4 %).
- The phase accumulated between consecutive zeros is 3.08–3.11 (π = 3.14) using the finite-ε εΦ.
- Decay per half period: predicted e^{1.28} = 3.6, observed 6.7 → 4.55 (z = 2 → 8). The ratio of the two is
  consistent with a z^{−1.5} prefactor for z = 4–8.

**2D Boussinesq** (εΦ from `bq_wkb2d.py` on the computed profiles):
- εΦ = 1.139 − 0.568i (z = 2.5), 1.087 − 0.577i (z = 4.0), 1.054 − 0.588i (z = 6.7), extrapolating to
  a ≈ 1.0 − 0.61i.
- Local spacing predictions are 1.38, 1.45 and 1.49; the observed spacings are 1.419, 1.457 and 1.467.
  The asymptotic spacing should be approached slowly, like O(1/z).
- Phase between consecutive zeros ≈ 2.9 (7 % below π), consistent with O(1/z) corrections to the leading-order
  local problem.

**Consequence.** Both models have an infinite ladder of smooth self-similar profiles accumulating at λ = 1, with
1/(λ_n − 1) growing linearly in n. The empirical law of Wang et al. is the leading-order WKB quantisation of the
quasi-stagnant region. The asymptotic 2D slope is π/(2 Re a) ≈ 1.5–1.56, larger than the slope 1.4187 fitted to
n ≤ 3.
