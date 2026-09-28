# P02 — Numerical method for self-similar blow-up profiles of 1D nonlocal transport equations

## Problem
CCF equation θ_t + (Hθ)θ_x = 0, H f(x) = (1/π) p.v.∫ f(y)/(x−y) dy.
Self-similar ansatz θ = (T−t)^λ Θ(ξ), ξ = x/(T−t)^{1+λ}; Θ even, Θ(0)=0:

        −λΘ + [(1+λ)ξ + HΘ(ξ)] Θ'(ξ) = 0.                                        (1)

Far field Θ ~ C|ξ|^β, β = λ/(1+λ) (matching to a time-independent outer solution).

## Key structural facts used
1. **First-order ODE with nonlocal coefficient.** (1) ⇔ d lnΘ / d lnξ = λ / d(ξ), with the *sonic factor*
   d(ξ) = 1 + λ + HΘ(ξ)/ξ > 0.
2. **Local exponent.** As ξ→0, HΘ/ξ → h1 = −(2/π)∫_0^∞ Θ/y² dy, so Θ ~ a ξ^p with
   p(λ) = λ / (1 + λ + h1).  For every λ in a range there is a solution (a continuous family); it is
   analytic at the origin iff p ∈ {2,4,…}. **Smooth self-similar profiles = level crossings p(λ)=2.**
3. **Scaling symmetry** Θ → κΘ(ξ/κ) (continuous family at fixed λ) must be fixed by a normalisation that
   is monotone along the orbit; we fix the far-field amplitude C (value of lnΘ − β lnξ at ξ=e^{50}).
   (A point-value normalisation Θ(1)=const meets each orbit twice and creates spurious folds.)

## Discretisation
- Logarithmic variable η = ln ξ, Θ = e^{cη}Ψ(η) with β < c < 1 so that Ψ decays exponentially at both ends.
- **Exact Hilbert transform in log variables.** For even Θ:
  HΘ(e^η) = e^{cη} (K_c * Ψ)(η), K_c(x) = e^{−cx}/(π sinh x) (p.v.), with Fourier multiplier
  K̂_c(k) = −i tanh(π(k − ic)/2)  (validated to 1e−12).
- **Integral (shooting) form** lnΨ(η) − lnΨ(η₀) − ∫_{η₀}^{η} [λ/d − c] dη' = 0, which avoids the spurious
  constraint that square spectral discretisations of the index-1 differential operator introduce.
- 8th-order cumulative quadrature.
- **Uniform grid** (FFT convolution, N up to 2^18) or **adaptive mapped grid** η = s − (1−ε)σ tanh((s−η_d)/σ)
  clustering points in the internal sonic layer; Hilbert convolution by the Sidi–Israeli alternating-point
  trapezoidal rule (spectrally accurate for p.v. kernels): (K*Ψ)_i = 2h_s Σ_{j−i odd} K_c(η_i−η_j) g'_j Ψ_j.
- Matrix-free Newton–Krylov (GMRES converges in 2–30 iterations: the Jacobian is identity + Fredholm-type
  operator); bordered systems for (i) smoothness constraint p=2 with λ unknown, (ii) pseudo-arclength
  continuation in (Ψ, λ); on-the-fly detection and refinement of every p=2 crossing; automatic re-gridding
  when the sonic layer thins.

## Linear stability (instability order)
Self-similar time τ=−ln(T−t); perturbations δ e^{μτ}:  μδ = λδ − d δ_η − (Hδ/ξ) Θ̄_η.
For Re μ > 0 the solution regular at ξ=0 is the particular solution
δ = T_μ δ,  (T_μ δ)(η) = −∫_{−∞}^{η} exp((λ−μ)∫_{η'}^{η} ds/d) (Hδ/ξ) Θ̄_η/d dη'
(the homogeneous solution ~ξ^{2(1−μ/λ)} is non-smooth and excluded). Eigenvalues μ ⇔ 1 ∈ spec(T_μ),
located by Arnoldi on T_μ along real μ. Checks: trivial time-translation mode μ=1 reproduced exactly;
λ₁ profile unstable eigenvalue μ≈0.366 (literature 0.36525).

## Cost
Each profile: 2–30 s on one CPU core (N≈5000 mapped points, 12-digit convergence).
