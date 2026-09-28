"""Leading eigenvalues of the linearised self-similar CCF operator by subspace iteration on the linear evolution
(RK4 + QR re-orthonormalisation every Δτ): returns the k largest Re μ (Lyapunov-type spectrum of the linear flow)."""
import numpy as np, sys, time
from ccf_mapped import CCFMapped
from ccf_nk import CCFNK

def setup(which):
    S0 = CCFNK(30, 120, 16384, 0.7)
    if which == 'l0':
        d = np.load('ladder_F16k_up_cross2.npy'); lam, eta, phi, etad, eps = d[0], S0.eta, d[1:], -1.0, 1.0
    elif which == 'l1':
        d = np.load('ladder_F16k_up_cross1.npy'); lam, eta, phi, etad, eps = d[0], S0.eta, d[1:], -1.45, 0.3
    else:
        d = np.load('mapped_l2_hs0.02_eps0.025.npy'); lam, eta, phi, etad, eps = 0.4713242277712, d[0], d[1], -0.9977, 0.025
    M = CCFMapped(30, 120, hs=0.03, eta_d=etad, eps=eps, sigma=2.0, c=0.7)
    ph, lam, ok = M.solve_p(M.from_other(eta, phi), lam, 2.0)
    return M, ph, lam

def run(which, k=4, T=30.0, dT=0.5):
    M, phi, lam = setup(which)
    Th = np.exp(phi + M.c * M.eta); d = 1 + lam + M.G(phi); Th_eta = lam * Th / d
    pre, post = np.exp(-M.c * M.eta), np.exp((M.c - 1) * M.eta)
    v = d / (M.gp * M.hs); dt = 0.4 / v.max()
    mask = np.abs(M.eta) < 15
    def rhs(D):   # D: (N, k)
        ds = np.zeros_like(D)
        ds[2:-1] = (2 * D[3:] + 3 * D[2:-1] - 6 * D[1:-2] + D[:-3]) / 6
        ds[1] = D[1] - D[0]; ds[-1] = D[-1] - D[-2]
        Hxi = post[:, None] * (M.Hm @ (pre[:, None] * D))
        out = lam * D - v[:, None] * ds - Hxi * Th_eta[:, None]
        out[0] = 0.0
        return out
    rng = np.random.default_rng(1)
    D = np.zeros((M.N, k))
    for j in range(k):
        r = np.zeros(M.N)
        for _ in range(12):
            c0, w0 = rng.uniform(-5, 5), rng.uniform(0.3, 2)
            r += rng.standard_normal() * np.exp(-((M.eta - c0) / w0) ** 2)
        D[:, j] = Th * r
    W = np.where(mask, 1.0, 0.0)[:, None]
    Q, R = np.linalg.qr(D * W); D = Q
    sums = np.zeros(k); t = 0.0; nre = 0; hist = []
    while t < T:
        tn = t + dT
        while t < tn - 1e-12:
            h = min(dt, tn - t)
            k1 = rhs(D); k2 = rhs(D + 0.5*h*k1); k3 = rhs(D + 0.5*h*k2); k4 = rhs(D + h*k3)
            D = D + h/6*(k1 + 2*k2 + 2*k3 + k4); t += h
        Q, R = np.linalg.qr(D * W)
        D = D @ np.linalg.inv(R)          # keep full vectors, normalised on the core
        if t > T / 3:
            sums += np.log(np.abs(np.diag(R))); nre += 1
            hist.append(np.log(np.abs(np.diag(R))) / dT)
    rates = sums / (nre * dT)
    return lam, rates, np.array(hist)

if __name__ == "__main__":
    which = sys.argv[1]; t0 = time.time()
    lam, rates, hist = run(which)
    print(f"{which}: lam={lam:.10f} leading Re(mu) (subspace iteration): {np.round(rates, 4)}  "
          f"(last-window: {np.round(hist[-5:].mean(0), 4)})  {time.time()-t0:.0f}s", flush=True)
