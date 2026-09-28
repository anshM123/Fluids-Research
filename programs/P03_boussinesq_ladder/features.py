"""Profile features along the branch from saved states: A, m, boundary Θ(y1,0), saturation, Ω extremum, cubic coeff."""
import numpy as np, sys, glob
from bq_solver import BQ
from bq_newton import full
files = sys.argv[1:]
print(f"{'lam':>8} {'m':>10} {'eps':>8} {'Th(1,0)':>10} {'Th(10,0)':>10} {'Th(1e3,0)':>10} {'s_half':>7} {'minOm':>9} {'r_minOm':>8} {'b_minOm':>8} {'b3':>9}")
for f in files:
    lam = float(f.split('lam')[1][:6]) if 'lam' in f else None
    B = BQ(lam)
    Y = np.load(f); X = full(B, Y)
    r = B.march(X / B.ea2[:, None], return_all=True)
    m, A = r['m'], r['A']; eps = 1 + lam - A
    cpow = np.where(B.cb > 0, B.cb, 0.0)
    Th = r['Th'] * cpow[None, :] ** m; Om = r['Omega']
    s = B.s
    thb = Th[:, 0]                      # boundary β = 0
    def at(sv): return np.interp(sv, s, thb)
    # s where boundary Θ/(-r^m) drops to 1/2 (saturation scale)
    ratio = thb / (-np.exp(m * s))
    k = np.argmax(ratio < 0.5); s_half = s[k] if ratio[k] < 0.5 else np.nan
    sel = (s > -8) & (s < 8)
    i, j = np.unravel_index(np.argmin(Om[sel]), Om[sel].shape)
    rmin = np.exp(s[sel][i]); bmin = B.beta[j]
    # cubic coefficient of V1 on the boundary: V_r/r(β=0) = ε + b3 r² + ...
    vr0 = (1 + lam) + r['Ur'][:, 0]
    ss = (s > -6) & (s < -3)
    b3 = np.polyfit(np.exp(2 * s[ss]), vr0[ss], 1)[0]
    print(f"{lam:8.4f} {m:10.6f} {eps:8.5f} {at(0):10.4e} {at(np.log(10)):10.4e} {at(np.log(1e3)):10.4e} {s_half:7.3f} {Om[sel].min():9.4f} {rmin:8.4f} {bmin:8.4f} {b3:9.4f}")
