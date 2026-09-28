# SI S1 — Theory of the terminal cusp of the CCF self-similar branch

## S1.1 Setting
Profile equation for θ = (T−t)^λ Θ(ξ), ξ = x/(T−t)^{1+λ}, with Θ even:

  ξ Θ'(ξ) · den(ξ) = λ Θ(ξ),    den(ξ) = 1 + λ + HΘ(ξ)/ξ,    Hf(x) = (1/π) p.v.∫ f(y)/(x−y) dy.       (1)

We use η = ln ξ, so d ln Θ/dη = λ/den.

Two quantities describe a solution:
- **Sonic depth** δ = min_ξ den(ξ), attained at ξ_s = e^{η_s}.
- **Local exponent** at the origin, p = λ/(1 + λ + h₁), with h₁ = lim_{ξ→0} HΘ/ξ = −(2/π)∫₀^∞ Θ/ξ² dξ.
  Then Θ ~ ξ^p as ξ → 0, and Θ is analytic at 0 iff p ∈ 2ℕ.

Hilbert transforms of power functions, used throughout (for −1 < Re σ < 1, and by continuation):
  H[|x|^σ] = −tan(πσ/2) sgn(x)|x|^σ,    H[sgn(x)|x|^σ] = cot(πσ/2) |x|^σ.                            (2)
In particular H[sgn(x)|x|^{1/2}] = |x|^{1/2}.

## S1.2 The square-root cusp (δ = 0)
Suppose den vanishes at ξ_s and den ≥ 0 nearby. Try Θ = Θ_s + B sgn(ξ−ξ_s)|ξ−ξ_s|^{1/2} + (regular terms).

- **Step 1: den.** By (2), the singular part of HΘ is B|ξ−ξ_s|^{1/2}, which is even and non-negative for
  B > 0. With the sonic condition 1 + λ + [regular part of HΘ](ξ_s)/ξ_s = 0,
    den(ξ) = (B/ξ_s)|ξ−ξ_s|^{1/2} + O(ξ−ξ_s).
- **Step 2: the ODE.** Insert this into (1):
    Θ'(ξ) = λΘ_s/(ξ_s den) ≈ λΘ_s/(B|ξ−ξ_s|^{1/2}).
  Integrating gives Θ − Θ_s = (2λΘ_s/B) sgn(ξ−ξ_s)|ξ−ξ_s|^{1/2}.
- **Step 3: self-consistency.** The coefficients must agree, so
    **B² = 2λΘ_s.**                                                                                        (3)

In η = ln ξ, write x = η − η_s. Then den ≈ k|x|^{1/2} with k = Bξ_s^{−1/2}, so
    k² = 2λΘ_s/ξ_s                                                                                          (4)
and ln Θ − ln Θ_s ≈ (2λ/k) sgn(x)|x|^{1/2}.

- The jump of ln Θ across the cusp is finite, so the terminal profile is continuous with a vertical tangent.
- (3) and (4) are checked to 1 % on the computed branch at δ = 1.7e−4, where the deviation is dominated by
  O(δ/√|x|) corrections.

## S1.3 Inner layer (δ > 0) and the width law
For small δ > 0 the cusp is regularised on the scale ℓ at which k|x|^{1/2} = δ, that is ℓ = δ²/k².

Write den = δ f(X), with X = x/ℓ, and ln Θ − ln Θ_s = (λℓ/δ) F(X), where F' = 1/f. Keep only the local
(singular) part of H and use (4). All constants cancel and the problem becomes

  f(X) − 1 = ½[HF(X) − HF(0)] = (1/π) p.v.∫₀^∞ F(y) X²/(y(X²−y²)) dy,    F odd,  F' = 1/f,  f(0) = 1.     (5)

Problem (5) contains no parameter. Its solution satisfies f → |X|^{1/2} automatically, because
p.v.∫₀^∞ t^{−1/2}/(1−t²) dt = π/2. Numerically (`inner_layer.py`, residual 7.5e−15), f''(0) = 0.08900.

The curvature width w = [den_min/(den''/2)]^{1/2} used in the code is therefore

  w = ℓ √(2/f''(0)) = 4.7405 δ²/k².                                                                         (6)

- With k² from (4): w/δ² = 9.00 using the cusp value, 9.17 using the value measured at δ = 1.7e−4.
- Computed: 9.098, constant to 0.5 % for δ ∈ [1e−5, 1e−2].

## S1.4 Linearisation about the cusp: universal operator
Perturb the terminal profile: Θ = Θ₀(1 + u), den = den₀ + den₁.

From (1), u' = −λ den₁/den₀². Near the cusp:
- den₁ ≈ (Θ_s/ξ_s) H[u], using (2) and the scale invariance of H;
- den₀² = k²|x| = (2λΘ_s/ξ_s)|x|, by (4).

Hence

  **u'(x) = −H[u](x) / (2|x|)**,                                                                             (7)

a universal nonlocal Euler-type operator: every profile-dependent constant has cancelled.

## S1.5 Local exponents
Insert power laws into (7) and use (2):
- **Even modes**, u = |x|^σ: σ = ½ tan(πσ/2). Roots: σ = 0 (scaling/amplitude); σ = ±½ (−½ is the
  translation mode −Θ₀'/Θ₀ ∝ |x|^{−1/2}); ±2.8910, ±4.9357, …
- **Odd modes**, u = sgn(x)|x|^σ: σ = −½ cot(πσ/2). Roots: **σ = ±iτ**, where
    τ tanh(πτ/2) = ½,   τ = 0.6494242963517879,                                                              (8)
  and the real roots ±1.8302, ±3.9192, …

With σ = iτ, the odd equation becomes iτ = −½ cot(iπτ/2) = ½ i coth(πτ/2), which is (8). The equation
y tanh y = c has no non-zero real root in y for c ≤ 0, and exactly one positive root for c > 0; here c = π/4.

## S1.6 Matching and the log-periodic spiral
The inner layer (5) perturbs ln Θ by O(λℓ/δ) = O(δ) at |x| ~ ℓ. Seen from the outer region, it acts as a
source localised at the cusp.

- Singular outer modes (Re σ < 0) are excited with amplitude O(δ) ℓ^{−σ}. The even translation mode σ = −½
  gives O(δ²).
- The marginal odd pair σ = ±iτ is excited with amplitude O(δ) and phase τ ln ℓ = 2τ ln δ + const.

Global functionals of the solution, in particular p and λ, therefore satisfy

  p(δ) = p* + δ[B_p + C_p cos(2τ ln δ) + D_p sin(2τ ln δ)] + o(δ),                                          (9)

and likewise for λ(δ). Consequences:
- The branch spirals into (λ*, p*) in the (λ, p) plane.
- Consecutive extrema of p(δ) are separated by the factor e^{π/(2τ)} = 11.2318 in δ.
- Their deviations from p* decrease by the same factor.

No δ ln δ term arises at this order because no exponent is degenerate; fits confirm a coefficient
consistent with zero.

Numerical confirmation:
- With 2τ free, fits over δ ≤ 3e−3 give 1.302 (from p) and 1.286–1.292 (from λ), against 2τ = 1.2988.
- Extrema of p: maximum 2.013428 at δ ≈ 1.5e−2; minimum 2.00495 at δ ≈ 1.6e−3; maximum 2.00584 at
  δ ≈ 1.5e−4; minimum ≈ 2.005766 at δ ≈ 1e−5.
- p* = 2.005772 ± 1e−6 and λ* = 0.4535845 ± 5e−7.

## S1.7 Finiteness of the ladder
- By (9), |p(δ) − p*| = O(δ). Since p* − 2 = 5.77e−3 ≠ 0, p(δ) = 2 has no solution for
  δ < δ_c ≈ 5.77e−3/max|B + C cos + D sin|.
- For δ above that threshold, the computed branch shows p ∈ (2.0049, 2.0134) all the way back to λ₂.
- On the other end (λ → ∞), p decreases monotonically from 2 to 1⁺.
- Hence the branch has exactly three crossings of p = 2: λ₀, λ₁ and λ₂.

## S1.8 Fractional generalisation (prediction)
For θ_t + (HΛ^sθ)θ_x = 0, 0 ≤ s < 1, the same argument gives:
- Cusp: Θ − Θ_s ∝ sgn(x)|x|^a with a = (1+s)/2, and den ∝ |x|^{(1−s)/2}.
- Universal operator: u' = −(a/c_s(a)) HΛ^s[u]/|x|^{1−s}, where HΛ^s[sgn|x|^σ] = c_s(σ)|x|^{σ−s} and
    c_s(σ) = cot(π(σ−s)/2) · 2^s Γ(1+σ/2)Γ((1+s−σ)/2)/[Γ((1−σ)/2)Γ(1+(σ−s)/2)].
- Odd roots: σ = s/2 ± iτ(s), with τ = 0.6494, 0.7046, 0.7594, 0.7872 at s = 0, 0.2, 0.5, 0.8.
- Spiral: layer width ∝ δ^{2/(1−s)}, oscillation amplitude ∝ δ^{1/(1−s)}, log-frequency 2τ(s)/(1−s).
- As s → 1, the operator becomes local (HΛ = −∂_x) and the spiral disappears, consistent with the monotone
  vanishing-order ladder of Burgers.
