"""Track smooth self-similar profiles (p=2) of the fractional CCF family θ_t + (HΛ^s θ)θ_x = 0 from s=0
(CCF: λ0=1.1808 stable, λ1=0.6057 1st unstable) towards s→1 (Burgers limit: explicit ladder λ_i = 1 + 1/i).
The weight c is kept at the centre of the admissible strip (β, 1+s)."""
import numpy as np, sys, time
from fccf_nk import FCCF

which = sys.argv[1]                       # 'l0' or 'l1'
smax = float(sys.argv[2]) if len(sys.argv) > 2 else 0.95
N = int(sys.argv[3]) if len(sys.argv) > 3 else 16384
L1, L2 = 30.0, 120.0
src = {"l0": "ladder_F16k_up_cross2.npy", "l1": "ladder_F16k_up_cross1.npy"}[which]
d = np.load(src)
lam, phi_c07 = d[0], d[1:]
eta_src = -30.0 + (150.0 / 16384) * np.arange(16384)            # grid of the source file (L1=30, L2=120, c=0.7)
lnTheta_src = phi_c07 + 0.7 * eta_src


def setup(s, lam):
    beta = (1 + s) * lam / (1 + lam)
    c = min(max(0.75, beta + 0.25), 1 + s - 0.1)
    return FCCF(s=s, L1=L1, L2=L2, N=N, c=c)


def to_grid(S, eta_old, lnTheta_old):
    # ln Θ is linear in η in both tails -> linear extrapolation outside the old range
    lt = np.interp(S.eta, eta_old, lnTheta_old)
    left = S.eta < eta_old[0]
    if left.any():
        sl = (lnTheta_old[50] - lnTheta_old[0]) / (eta_old[50] - eta_old[0])
        lt[left] = lnTheta_old[0] + sl * (S.eta[left] - eta_old[0])
    return lt - S.c * S.eta


S = setup(0.0, lam)
phi = to_grid(S, eta_src, lnTheta_src)
phi, lam, ok = S.solve_p(phi, lam, 2.0, tol=1e-10)
print(f"{which} s=0: lam={lam:.12f} ok={ok}", flush=True)
rows = [(0.0, lam, S.c, float((S.b(lam) + S.G(phi)).min()))]
s, ds = 0.0, 0.02
eta_prev, lt_prev = S.eta.copy(), phi + S.c * S.eta
t0 = time.time()
while s < smax - 1e-12:
    s_new = min(s + ds, smax)
    Sn = setup(s_new, lam)
    ph = to_grid(Sn, eta_prev, lt_prev)
    ph, ln, ok = Sn.solve_p(ph, lam, 2.0, tol=1e-10)
    if ok and abs(ln - lam) > 0.05 + 3 * ds:
        ok = False                       # jumped to another branch
    if not ok:
        ds *= 0.5
        print(f"   fail at s={s_new:.4f}; ds -> {ds:.4f}", flush=True)
        if ds < 1e-4:
            break
        continue
    s, lam, S, phi = s_new, ln, Sn, ph
    den = S.b(lam) + S.G(phi)
    rows.append((s, lam, S.c, float(den.min())))
    print(f"{which} s={s:.4f}: lam={lam:.10f} c={S.c:.3f} min_den={den.min():.4e} at eta={S.eta[np.argmin(den)]:.2f} "
          f"(Burgers-limit ladder 2, 1.5, 1.333…)  t={time.time()-t0:.0f}s", flush=True)
    eta_prev, lt_prev = S.eta.copy(), phi + S.c * S.eta
    np.save(f"fccf_track_{which}.npy", np.array(rows))
    np.save(f"fccf_{which}_s{s:.3f}.npy", np.concatenate([[lam, S.c], phi]))
    ds = min(ds * 1.2, 0.05)
