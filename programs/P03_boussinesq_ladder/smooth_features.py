"""Local features of a smooth profile (m = 2): boundary velocity V_r/r − ε ≈ b3 r² + ..., dip, knee, Θ saturation."""
import numpy as np, sys
from bq_solver import BQ
from bq_newton import full
import re
for f in sys.argv[1:]:
    lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    B = BQ(lam); Y = np.load(f); X = full(B, Y)
    r = B.march(X / B.ea2[:, None], return_all=True)
    lam, A, m = B.lam, r['A'], r['m']; eps = 1 + lam - A
    vr = (1 + lam) + r['Ur']; s = B.s
    ss = (s > -6) & (s < -2.5)
    rr = np.exp(2 * s[ss]); co = np.polyfit(rr, (vr[ss, 0] - eps) / rr, 2)
    k = np.argmin(vr[:, 0])
    cpow = np.where(B.cb > 0, B.cb, 0.0)
    Th = r['Th'] * cpow[None, :] ** m
    ratio = Th[:, 0] / (-np.exp(m * s)); kk = np.argmax(ratio < 0.5)
    print(f"λ={lam:.8f} m={m:.8f} ε={eps:.6f} b3={co[-1]:+.5f} b5={co[-2]:+.4f} b3/ε={co[-1]/eps:+.4f}  dip: vr={vr[k,0]:.5f} "
          f"ratio={vr[k,0]/eps:.4f} at r={np.exp(s[k]):.4f};  half-saturation r={np.exp(s[kk]):.4f};  Θ(r=10,β=0)={np.interp(np.log(10), s, Th[:,0]):.4f}")
