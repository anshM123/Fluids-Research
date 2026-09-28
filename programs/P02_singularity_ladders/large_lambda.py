"""The other end of the CCF branch: natural continuation in λ from λ0 up to large λ (uniform FFT log grid;
no internal layers there).  The weight is kept inside the admissible strip, c = β + 0.4(1−β), and the domain
is lengthened as the far-field decay rate c−β of Ψ shrinks.  Checks that p(λ) stays below 2."""
import numpy as np, time
from ccf_nk import CCFNK

d = np.load("ladder_F16k_up_cross2.npy")               # λ0 profile, grid L1=30, L2=120, N=16384, c=0.7
lam = d[0]
eta_old = -30.0 + (150.0 / 16384) * np.arange(16384)
lt_old = d[1:] + 0.7 * eta_old                          # ln Θ


def make(lam):
    beta = lam / (1 + lam)
    c = beta + 0.4 * (1 - beta)
    L2 = float(min(max(120.0, 36.0 / (c - beta)), 3000.0))
    N = int(2 ** np.ceil(np.log2((30.0 + L2) / 0.009)))
    return CCFNK(30.0, L2, N, c)


def transfer(S, eta_old, lt_old, beta):
    lt = np.interp(S.eta, eta_old, lt_old)
    right = S.eta > eta_old[-1]
    if right.any():
        lt[right] = lt_old[-1] + beta * (S.eta[right] - eta_old[-1])
    return lt - S.c * S.eta


rows = []
t0 = time.time()
lams = list(np.arange(1.2, 2.0, 0.1)) + list(np.arange(2.0, 6.0, 0.25)) + list(np.arange(6.0, 20.01, 1.0))
S = None
for L in lams:
    S = make(L)
    phi = transfer(S, eta_old, lt_old, lam / (1 + lam))
    ph, ok = S.solve_fixed(phi, L)
    if not ok:
        print(f"fail at λ={L} (N={S.N}, L2={S.L2:.0f}, c={S.c:.4f})", flush=True)
        break
    lam = L
    eta_old, lt_old = S.eta.copy(), ph + S.c * S.eta
    F, den = S.F(ph, lam, S.normval(ph, lam))
    p = S.p_of(ph, lam)
    rows.append((lam, p, den.min(), S.h1(ph)))
    print(f"λ={lam:7.3f} p={p:.8f} min_den={den.min():.4f} h1={S.h1(ph):.6f} β={lam/(1+lam):.5f} "
          f"c={S.c:.4f} L2={S.L2:.0f} N={S.N} t={time.time()-t0:.0f}s", flush=True)
    np.save("large_lambda_branch.npy", np.array(rows))
