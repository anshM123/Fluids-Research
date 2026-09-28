import numpy as np, sys
from fccf_nk import FCCF
d = np.load("ladder_F16k_up_cross2.npy"); lam = d[0]
eta_src = -30.0 + (150.0 / 16384) * np.arange(16384); lt = d[1:] + 0.7 * eta_src
for s in [0.0, 0.01, 0.02, 0.03, 0.04]:
    for N, c in [(16384, 0.75), (32768, 0.75)]:
        S = FCCF(s=s, L1=30, L2=120, N=N, c=c)
        phi = np.interp(S.eta, eta_src, lt) - c * S.eta
        phi, lam2, ok = S.solve_p(phi, lam, 2.0, tol=1e-10)
        F, den = S.F(phi, lam2, S.normval(phi, lam2))
        print(f"s={s} N={N} c={c}: lam={lam2:.10f} ok={ok} |F|={np.abs(F).max():.1e} minden={den.min():.4f} h_s={S.h1(phi):.6f} p={S.p_of(phi, lam2):.8f}", flush=True)
        if N == 16384:
            lam_keep, lt_keep = lam2, phi + c * S.eta
    lam, lt, eta_src = lam_keep, lt_keep, FCCF(s=0, L1=30, L2=120, N=16384, c=0.75).eta
