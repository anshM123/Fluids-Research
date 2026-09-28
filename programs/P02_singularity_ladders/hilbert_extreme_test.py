"""Accuracy of the alternating-point Hilbert matrix on geometric grids with extreme centre refinement:
compare Hm @ Ψ with the exact log-variable transform (FFT on a fine uniform grid) for a smooth Ψ."""
import numpy as np
from scipy.interpolate import CubicSpline
from ccf_sinhgrid import CCFSinh
from ccf_nk import CCFNK

S = CCFNK(30.0, 120.0, 2 ** 17, 0.7)
f = lambda e: np.exp(-(e + 0.9) ** 2 / 2.0) + 0.3 * np.exp(-(e - 2.0) ** 2)
exact = CubicSpline(S.eta, S.conv(f(S.eta)))
for hc in [1e-4, 1e-6, 1e-8, 1e-10, 1e-11]:
    M = CCFSinh(30.0, 120.0, hs=0.03, eta_d=-0.96, hc=hc, B=1.0, c=0.7)
    approx = M.Hm @ f(M.eta)
    sel = np.abs(M.x) < 0.5
    err = np.abs(approx[sel] - exact(M.eta[sel]))
    ic = int(np.argmin(np.abs(M.x)))
    print(f"hc={hc:.0e} N={M.N}: max |Hm·Ψ − exact| for |η−η_d|<0.5: {err.max():.2e}; at centre {abs(approx[ic]-exact(M.eta[ic])):.2e}",
          flush=True)
    del M
