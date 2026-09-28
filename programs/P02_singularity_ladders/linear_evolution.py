"""Independent check of instability: evolve the linearised CCF equation in self-similar time,
        δ_τ = λ δ − d̄ δ_η − (Hδ/ξ) Θ̄_η ,
on the mapped grid (3rd-order upwind in s, RK4), inflow boundary δ=0 at η=−L1, from a random smooth
perturbation vanishing like ξ² at the origin. The asymptotic growth rate equals max Re μ over the
admissible spectrum (time-translation mode gives μ=1)."""
import numpy as np, sys
from ccf_mapped import CCFMapped


def growth(M, phi, lam, T=25.0, seed=0, project_T=False):
    Th = np.exp(phi + M.c * M.eta)
    d = 1 + lam + M.G(phi)
    Th_eta = lam * Th / d
    pre, post = np.exp(-M.c * M.eta), np.exp((M.c - 1) * M.eta)
    v = d / (M.gp * M.hs)               # transport speed in s-units
    dt = 0.4 / v.max()
    rng = np.random.default_rng(seed)
    # smooth random perturbation ~ Θ̄ × (random combination of bumps)
    r = np.zeros(M.N)
    for _ in range(12):
        c0, w0 = rng.uniform(-5, 5), rng.uniform(0.3, 2)
        r += rng.standard_normal() * np.exp(-((M.eta - c0) / w0) ** 2)
    delta = Th * (1 + r) * np.exp(-np.maximum(M.eta, 0) * 0.0)

    def rhs(dl):
        ds = np.zeros_like(dl)
        # 3rd-order upwind (speed > 0): (2 f_{i+1} + 3 f_i − 6 f_{i−1} + f_{i−2}) / 6
        ds[2:-1] = (2 * dl[3:] + 3 * dl[2:-1] - 6 * dl[1:-2] + dl[:-3]) / 6
        ds[1] = dl[1] - dl[0]
        ds[-1] = dl[-1] - dl[-2]
        Hxi = post * (M.Hm @ (pre * dl))
        out = lam * dl - v * ds - Hxi * Th_eta
        out[0] = 0.0
        return out
    t, hist = 0.0, []
    nrm0 = np.linalg.norm(delta)
    while t < T:
        k1 = rhs(delta); k2 = rhs(delta + 0.5 * dt * k1); k3 = rhs(delta + 0.5 * dt * k2); k4 = rhs(delta + dt * k3)
        delta = delta + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        n = np.linalg.norm(delta[np.abs(M.eta) < 15])
        hist.append((t, n))
        if n > 1e200:
            break
    h = np.array(hist)
    lg = np.log(h[:, 1])
    rates = np.gradient(lg, h[:, 0])
    return h, rates, delta


if __name__ == "__main__":
    import time
    from ccf_nk import CCFNK
    S0 = CCFNK(30, 120, 16384, 0.7)
    which = sys.argv[1]
    if which == 'l0':
        d = np.load('ladder_F16k_up_cross2.npy'); lam, eta, phi, etad, eps = d[0], S0.eta, d[1:], -1.0, 1.0
    elif which == 'l1':
        d = np.load('ladder_F16k_up_cross1.npy'); lam, eta, phi, etad, eps = d[0], S0.eta, d[1:], -1.45, 0.3
    else:
        d = np.load('mapped_l2_hs0.02_eps0.025.npy'); lam, eta, phi, etad, eps = 0.4713242277712, d[0], d[1], -0.9977, 0.025
    M = CCFMapped(30, 120, hs=0.03, eta_d=etad, eps=eps, sigma=2.0, c=0.7)
    ph, lam, ok = M.solve_p(M.from_other(eta, phi), lam, 2.0)
    t0 = time.time()
    h, rates, delta = growth(M, ph, lam)
    n = len(rates)
    print(f"{which}: lam={lam:.10f}  growth rate at τ≈{h[n//2,0]:.1f}: {rates[n//2]:.4f}, τ≈{h[3*n//4,0]:.1f}: {rates[3*n//4]:.4f}, "
          f"τ≈{h[-1,0]:.1f}: {rates[-5]:.4f}  ({time.time()-t0:.0f}s)")
    np.save(f"linevo_{which}.npy", h)
