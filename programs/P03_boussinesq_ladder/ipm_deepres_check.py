"""Resolution check of the deep IPM branch: the same z solved on h_s = 0.0125 (continuation e5) and on h_s = 0.00625
(ipm_shift.py from the e5 state).  For each grid: m − 2, the wall-dip geometry and the resolution indicator of
ipm_deep.py (|Ω_b + ∂ₓR| / max ∂ₓR over the dip-to-front window: the IPM identity Ω_b = −∂ₓR holds exactly for the
continuous problem), the travel time T = ∫ds/D̂ and the repaired phase (wkb_phase3).
usage: ipm_deepres_check.py Z STATE_0125 STATE_00625 [Z_0125]  (writes one line per grid; deepres_z{Z}_hs0.00625.out is
the solver log).  Z_0125 is the exact z of the continuation state (e5 steps from z₀ = 1/0.1339405 = 7.4660017, so
its z differ from the rounded values by 1.7e-6): m is recomputed by marching the state, and a march at a λ the
state was not solved for gives a meaningless m (the O(1) geometry and phase are unaffected)."""
import numpy as np, sys, os
from ipm_solver import IPM
from bq_newton import full
from ipm_wkb import wkb_phase3, wall_data
from scipy.integrate import cumulative_trapezoid

z = float(sys.argv[1])
z125 = float(sys.argv[4]) if len(sys.argv) > 4 else z
for f, hs, lam in ((sys.argv[2], 0.0125, 1 / z125), (sys.argv[3], 0.00625, 1 / z)):
    B = IPM(lam, hs=hs, Nb=32, s_sw=12.0, s_start=-20.0); Y = np.load(f)
    r = B.march(full(B, Y) / B.ea2[:, None], return_all=True)
    D = (1 + lam) + r['Ur'][:, 0]; D0 = 1 + lam - r['A']; Dh = D / D0
    sel = (B.s > -3) & (B.s < 2); i = np.where(sel)[0][np.argmin(Dh[sel])]
    a, b, c = Dh[i - 1], Dh[i], Dh[i + 1]; den = a - 2 * b + c
    Dmin_p = b - 0.125 * (a - c) ** 2 / den if den > 0 else b
    kc = i + np.argmax(Dh[i:] > 3.0)
    G = np.exp(-B.s) * np.gradient(r['Th'][:, 0], B.s); Om = r['Omega'][:, 0]
    w = slice(max(i - 20, 0), kc + 1)
    mis = float(np.max(np.abs(G[w] + Om[w])) / np.max(np.abs(G[w])))
    icut = i + np.argmax(Dh[i:] > 2.0); s_cut = np.interp(2.0, Dh[icut - 1:icut + 1], B.s[icut - 1:icut + 1])
    i0 = np.searchsorted(B.s, -10.0)
    g = np.concatenate([B.s[i0:icut], [s_cut]]); yv = np.concatenate([1 / Dh[i0:icut], [0.5]])
    T = float(cumulative_trapezoid(yv, g)[-1])
    o = wkb_phase3(f, lam=lam, Nb=32, hs=hs, Dcut=2.0)
    print(f"z={1 / lam:.7f} hs={hs}: m-2={r['m'] - 2:+.3e} D̂min={Dh[i]:.4f} (parab {Dmin_p:.4f}) x_dip={np.exp(B.s[i]):.4f} "
          f"x_cut={np.exp(s_cut):.4f} max∂ₓR={np.max(G[w]):.2f} mis={mis:.2e} T={T:.4f} "
          f"ReΦ={o['Phi'].real:.4f} I={o['I'].real:.5f} front_fail={o['nfail_front']}", flush=True)
