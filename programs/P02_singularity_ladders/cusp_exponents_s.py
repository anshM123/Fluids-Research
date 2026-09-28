"""Cusp exponents for the fractional family θ_t + (HΛ^sθ)θ_x = 0 (prediction, IDEA118).
At a sonic point the profile has Θ − Θ_s ∝ sgn(x)|x|^a with a = (1+s)/2, and the linearised profile equation
reduces to the universal operator  u' = −(a/c_s(a)) HΛ^s[u] / |x|^{1−s}.  Power modes:
   even u = |x|^σ :        σ = −(a/c_s(a)) e_s(σ),     HΛ^s|x|^σ      = e_s(σ) sgn(x)|x|^{σ−s}
   odd  u = sgn(x)|x|^σ :  σ = −(a/c_s(a)) c_s(σ),     HΛ^s[sgn|x|^σ] = c_s(σ) |x|^{σ−s}
e_s(σ) = 2^s Γ((1+σ)/2)Γ((1+s−σ)/2)/[Γ(−σ/2)Γ(1+(σ−s)/2)]  (validated in fccf_nk.py),
c_s(σ) = cot(π(σ−s)/2) · 2^s Γ(1+σ/2)Γ((1+s−σ)/2)/[Γ((1−σ)/2)Γ(1+(σ−s)/2)].
s = 0 reproduces σ = ½tan(πσ/2), σ = −½cot(πσ/2) (τ tanh(πτ/2) = ½).  Tracks the complex odd root in s."""
import numpy as np
from scipy.special import loggamma


def c_s(s, z):
    return (np.cos(np.pi * (z - s) / 2) / np.sin(np.pi * (z - s) / 2)) * np.exp(
        s * np.log(2.0) + loggamma(1 + z / 2) + loggamma((1 + s - z) / 2) - loggamma((1 - z) / 2) - loggamma(1 + (z - s) / 2))


def odd_eq(s, z):
    a = (1 + s) / 2
    return z + (a / c_s(s, a)) * c_s(s, z)


def newton(fun, z, tol=1e-14):
    for _ in range(100):
        f = fun(z); h = 1e-7
        df = (fun(z + h) - fun(z - h)) / (2 * h)
        dz = -f / df
        z = z + dz
        if abs(dz) < tol:
            return z, True
    return z, False


z = 0.6494242963517879j
print(" s      a=(1+s)/2   odd root σ = ρ + iτ        spiral: amplitude δ-power, log-frequency in ln δ")
for s in np.round(np.arange(0.0, 0.951, 0.05), 3):
    z, ok = newton(lambda q: odd_eq(s, q), z)
    a = (1 + s) / 2
    # den ≈ δ + k|x|^{a−s}: layer width w ∝ δ^{2/(1−s)}; the layer perturbs ln Θ by O(λw/δ) = O(δ^{(1+s)/(1−s)});
    # mode amplitude δ^{(1+s)/(1−s)}·w^{−σ} ∝ δ^{(1+s−2ρ)/(1−s)}, phase (2τ/(1−s)) ln δ
    print(f"{s:4.2f}   {a:.3f}      {z.real:+.6f} {z.imag:+.6f}i   {'ok' if ok else 'NO CONV'}   "
          f"amp ∝ δ^{(1+s-2*z.real)/(1-s):.4f}, ω = {2*abs(z.imag)/(1-s):.4f}")
