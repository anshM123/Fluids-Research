"""Test of hypothesis (U) of ASYMPTOTICS.md: at fixed s inside the stalled layer, do the base data D̂(s) and the local
root κ(s) converge as ε → 0 with O(ε) corrections? Uses the crossing profiles n = 3…7."""
import numpy as np, warnings
warnings.filterwarnings('ignore')
from bq_wkb2d import wkb_phase, boundary_params
files = ["Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.184253.npy"] + \
        [f"Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam{l}.npy" for l in ("1.144986", "1.119474", "1.101582", "1.088338")]
svals = np.array([-3.0, -2.0, -1.5, -1.0, -0.8, -0.6, -0.5])
E, K, DH = [], [], []
for f in files:
    lam, m, eps, s, D, mu, G, Thy, Omb = boundary_params(f)
    r = wkb_phase(f, return_kappa=True)
    grid, kap = r[6], r[7]
    k_at = np.interp(svals, grid, kap.real) + 1j * np.interp(svals, grid, kap.imag)
    k_at[svals > grid[-1]] = np.nan
    E.append(eps); K.append(k_at); DH.append(np.interp(svals, s, D / eps))
    print(f"ε={eps:.5f}: D̂(s)={np.round(DH[-1], 4).tolist()}  κ(s)={np.round(k_at, 4).tolist()}", flush=True)
E = np.array(E); K = np.array(K); DH = np.array(DH)
print("\nlinear-in-ε fit at fixed s (value at ε=0, slope, max residual):")
for j, sv in enumerate(svals):
    ok = np.isfinite(K[:, j])
    for name, Y in (("D̂", DH[:, j]), ("Re κ", K[:, j].real), ("Im κ", K[:, j].imag)):
        A = np.vstack([np.ones(ok.sum()), E[ok]]).T
        c, *_ = np.linalg.lstsq(A, Y[ok], rcond=None)
        res = np.abs(Y[ok] - A @ c).max()
        print(f"  s={sv:5.2f} {name:5s}: ε→0 limit {c[0]:+.4f}, slope {c[1]:+.3f}, max residual {res:.1e}")
