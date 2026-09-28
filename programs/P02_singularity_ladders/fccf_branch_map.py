"""Branch map p(λ; s) of the fractional CCF family at fixed s (uniform FFT log grid, natural continuation
in λ downward from the stable-profile continuation λ₀(s)).  Counts p=2 crossings (smooth profiles) that are
resolvable on the uniform grid (sonic margin ≳ 0.02).  Usage: fccf_branch_map.py S [N]"""
import numpy as np, sys, glob, time
from fccf_nk import FCCF

s = float(sys.argv[1])
N = int(sys.argv[2]) if len(sys.argv) > 2 else 65536
files = sorted(glob.glob("fccf_l0_s*.npy"), key=lambda f: abs(float(f[9:14]) - s))
d = np.load(files[0]); s_file = float(files[0][9:14])
lam, c = d[0], d[1]
phi16 = d[2:]
eta16 = -30.0 + (150.0 / phi16.size) * np.arange(phi16.size)
S = FCCF(s=s, L1=30.0, L2=120.0, N=N, c=c)
phi = np.interp(S.eta, eta16, phi16)
phi, lam, ok = S.solve_p(phi, lam, 2.0, tol=1e-10)
print(f"s={s} (start file s={s_file}): λ0(s)={lam:.10f} ok={ok} c={c:.3f}", flush=True)
rows, cross = [], []
t0 = time.time()
dl = -0.005
phi0 = S.normval(phi, lam)
p_prev = 2.0
L = lam
while L > 0.05:
    Ln = L + dl
    beta = Ln / S.b(Ln)
    if beta >= S.c - 0.02:
        break
    ph, okf = S.solve_fixed(phi.copy(), Ln, phi0=S.normval(phi, L) + (S.c - beta) * 0 , maxit=40)
    if not okf:
        dl *= 0.5
        if abs(dl) < 2e-5:
            print(f"   stop at λ={L:.5f} (natural continuation fails: fold or unresolved layer)", flush=True)
            break
        continue
    F, den = S.F(ph, Ln, S.normval(ph, Ln))
    p = S.p_of(ph, Ln)
    rows.append((Ln, p, den.min(), S.eta[np.argmin(den)]))
    if (p_prev - 2) * (p - 2) < 0:
        cross.append(0.5 * (L + Ln))
        print(f"   p=2 crossing near λ={0.5*(L+Ln):.5f}", flush=True)
    p_prev, L, phi = p, Ln, ph
    if len(rows) % 10 == 0:
        print(f"   λ={L:.4f} p={p:.6f} min_den={den.min():.4f} at η={S.eta[np.argmin(den)]:.2f} t={time.time()-t0:.0f}s", flush=True)
    dl = max(dl * 1.2, -0.01)
np.save(f"fccf_bmap_s{s:.3f}.npy", np.array(rows))
print(f"s={s}: crossings (resolved part): λ0={lam:.6f} + {['%.5f' % x for x in cross]}", flush=True)
