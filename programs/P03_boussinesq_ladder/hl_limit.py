"""The λ → 1 limit of the Hou–Luo self-similar profile, solved directly (no finite-ε data).

At λ = 1 the profile splits into two regions (front at ξ = 1 after the scaling ξ → ξ/x_c):
  * stalled layer |ξ| < 1: D = 0, i.e. U = −2ξ, i.e. H[Ω] = −2 on (−1, 1)   (a finite-Hilbert / airfoil equation),
                           Θ(ξ) = ∫_0^ξ Ω   (vorticity slaved to the buoyancy gradient);
  * outer region  ξ > 1:   Θ = Θ(1) constant, Ω transported:  Ω + ξ D Ω_ξ = 0,  D = 2 + U/ξ,  U' = H[Ω],
                           Ω(1⁺) = Ω_f   (one closure number, set by the front inner problem).
Ω is odd. Layer: Ω = Σ_j b_j T_j(x)/√(1−x²) (j odd), which carries the square-root endpoint singularity exactly;
T[Ω] = −Σ b_j U_{j−1} on the layer, Σ b_j (x − √(x²−1))^j/√(x²−1) outside. Outer: y = 1 + σ² grid (log grid far out).

The WKB phase constant of the profile ladder is a property of this limit solution:
    a₀ = ∫_0^1 κ₀ dξ/ξ,   κ₀ = (−1 + √(1 + 4iΩ))/(2i D̂),   D̂ = 2Θ/(ξΩ)   (the layer relation DΘ_η = (λ−1)Θ),
    C = 2 Re a₀,   predicted asymptotic spacing of z_n = 1/(λ_n − 1):  a = π/C.
usage: python3 hl_limit.py OMEGA_F [OMEGA_F ...]"""
import numpy as np, sys

J = 160                                    # odd Chebyshev modes j = 1, 3, ..., 2J−1
jj = np.arange(1, 2 * J, 2)
M = 1200                                   # Gauss–Chebyshev (2nd kind) nodes on the layer
th_k = np.arange(1, M + 1) * np.pi / (M + 1)
x_k = np.cos(th_k); w_k = np.pi / (M + 1) * np.sin(th_k) ** 2          # ∫ f √(1−x²) dx ≈ Σ w_k f(x_k)
Uk = np.sin(np.outer(jj, th_k)) / np.sin(th_k)                          # U_{j−1}(x_k), shape (J, M)

# outer grid: y = 1 + σ², σ ∈ (0, 3], then logarithmic to 1e8
sig = np.linspace(0, 3.0, 6001)[1:]
y1 = 1 + sig ** 2; dy1 = 2 * sig * (sig[1] - sig[0])
ylog = np.exp(np.linspace(np.log(10.0), np.log(1e8), 3000))[1:]
y = np.concatenate([y1, ylog])
dy = np.concatenate([dy1, np.gradient(ylog)])


def layer_H_outside(b, x):
    """H of the layer vorticity at x > 1"""
    s = np.sqrt(x ** 2 - 1)
    r = x - s
    return np.sum(b[:, None] * r[None, :] ** jj[:, None], axis=0) / s


def _kernel_on_layer(x):
    return 2 * y[None, :] / (x[:, None] ** 2 - y[None, :] ** 2) * dy[None, :] / np.pi


def _kernel_outside():
    """matrix of the principal-value operator Ω_O ↦ (1/π) p.v.∫_1^∞ Ω_O(y) 2y/(x² − y²) dy at the outer grid points,
    with the local subtraction on |y − x| < d (d = min(x − 1, x/2)) and the analytic p.v. of the kernel there"""
    n = len(y)
    K = np.empty((n, n))
    for i in range(n):
        x = y[i]; d = min(x - 1, x / 2)
        with np.errstate(divide='ignore'):
            row = 2 * y / (x ** 2 - y ** 2) * dy
        row[i] = 0.0
        near = np.abs(y - x) < d
        row_i = -row[near].sum() + np.log((2 * x - d) / (2 * x + d))
        K[i] = row; K[i, i] += row_i
    return K / np.pi


K_LAYER = _kernel_on_layer(x_k)
K_OUT = _kernel_outside()


def outer_H_on_layer(Om_o, x=None):
    return K_LAYER @ Om_o


def outer_H_outside(Om_o, xs=None):
    return K_OUT @ Om_o


def solve(Om_f, iters=400, tol=1e-9, verbose=False, init=None, relax=0.3):
    """fixed-point iteration for given Ω_f, started from `init` (a previous solution) or from the Ω_f = 0 solution"""
    if init is None:
        Om_o = np.zeros_like(y); b = np.zeros(J); b[0] = 2.0      # Ω_f = 0: Ω = 2ξ/√(1−ξ²), D = 2√(1 − 1/ξ²)
    else:
        Om_o = init['Om_o'] * (Om_f / init['Om_f'] if init['Om_f'] > 0 else 1.0); b = init['b'].copy()
        if init['Om_f'] == 0:
            Om_o = Om_f * init['shape']
    err = np.inf
    for it in range(iters):
        g = -2.0 - outer_H_on_layer(Om_o, x_k)
        b_new = -(2 / np.pi) * (Uk * (w_k * g)[None, :]).sum(axis=1)
        H = layer_H_outside(b_new, y) + outer_H_outside(Om_o, y)
        U = -2.0 + np.cumsum(H * dy)
        D = 2.0 + U / y
        Dp = np.maximum(D, 1e-8)
        Om_new = Om_f * np.exp(-np.cumsum(dy / (y * Dp)))
        err = max(np.abs(b_new - b).max(), np.abs(Om_new - Om_o).max() / max(Om_f, 1e-12))
        b = (1 - relax) * b + relax * b_new; Om_o = (1 - relax) * Om_o + relax * Om_new
        if verbose and it % 20 == 0:
            print(f"   it {it}: change {err:.2e}  Σb = {b.sum():.6f}  min D = {D.min():.4f}", flush=True)
        if err < tol:
            break
    return dict(ok=bool(D.min() > 0 and err < 1e-6), it=it, b=b, Om_o=Om_o, D=D, U=U, H=H, err=err, Om_f=Om_f,
                Dmin=float(D.min()), shape=np.exp(-np.cumsum(dy / (y * np.maximum(D, 1e-8)))))


def layer_fields(b, th):
    """Ω, Θ on the layer at x = cos θ"""
    x = np.cos(th)
    Om = (b[:, None] * np.cos(np.outer(jj, th))).sum(axis=0) / np.sin(th)
    Th = (b[:, None] * ((np.sin(jj * np.pi / 2) / jj)[:, None] - np.sin(np.outer(jj, th)) / jj[:, None])).sum(axis=0)
    return x, Om, Th


def phase(b, n=400000):
    """a₀ = ∫_0^1 κ₀ dx/x; in θ (x = cos θ): dx/x = tan θ dθ, integrand ~ θ^{−1/2} at the front → θ = u², dθ = 2u du"""
    u = (np.arange(n) + 0.5) / n * np.sqrt(np.pi / 2)
    th = u ** 2; dth = 2 * u * (u[1] - u[0])
    x, Om, Th = layer_fields(b, th)
    Dh = 2 * Th / (x * Om)
    kap = (-1 + np.sqrt(1 + 4j * Om)) / (2j * Dh)
    return np.sum(kap * np.tan(th) * dth)


if __name__ == "__main__":
    targets = [float(v) for v in sys.argv[1:]]
    prev = dict(Om_f=0.0, b=None, Om_o=np.zeros_like(y), shape=np.exp(-np.cumsum(dy / (y * 2 * np.sqrt(1 - 1 / y ** 2)))))
    prev['b'] = np.zeros(J); prev['b'][0] = 2.0
    Om_f = 0.0
    for tgt in targets:
        while Om_f < tgt - 1e-12:
            Om_f = min(tgt, Om_f + 0.1)
            r = solve(Om_f, init=prev)
            prev = r
        b = r['b']
        Theta_c = float(np.sum(b * np.sin(jj * np.pi / 2) / jj))
        th0 = np.array([np.pi / 2 - 1e-4]); _, Om0, _ = layer_fields(b, th0)
        c0 = float(Om0[0] / np.cos(th0[0]) / 2)        # Θ ≈ c0 ξ² at the origin (units x_c = 1)
        alpha = b.sum() / np.sqrt(2)
        a0 = phase(b); C = 2 * a0.real
        print(f"Ω_f = {Om_f:.3f}: ok={r['ok']} ({r['it']} it, change {r['err']:.1e}, min D out {r['Dmin']:.4f}); "
              f"Θ(1) = {Theta_c:.5f}, c0 = {c0:.5f} (x_c in the Θ≈ξ² normalisation), Θ(1)/c0... ; α = {alpha:.4f}; "
              f"a₀ = {a0.real:.5f} {a0.imag:+.5f}i, C = {C:.5f}, π/C = {np.pi / C:.5f}", flush=True)
        np.savez(f"hl_limit_Omf{Om_f:.3f}.npz", b=b, y=y, Om_o=r['Om_o'], D=r['D'], U=r['U'])
