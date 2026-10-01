"""Deep continuation of the IPM least-singular branch in z = 1/λ (endpoint study, task: does the dip of the wall
speed close at finite λ, does the branch fold, or does it continue?).  Saves every state (ipm_deep_{tag}_lam*.npy,
not committed) for the WKB phase (ipm_wkb.py) and records m, A, min V_r/r and the wall dip (position, depth).
usage: ipm_deep.py LAM_START LAM_END DZ TAG START.npy [Nb hs NB_IN HS_IN tol]"""
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import newton, full
from bq_regrid import transfer
from bq_logpolar import BQLogPolar

lam0, lam1, dz, tag, start = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], sys.argv[5]
Nb = int(sys.argv[6]) if len(sys.argv) > 6 else 32
hs = float(sys.argv[7]) if len(sys.argv) > 7 else 0.025
nbi = int(sys.argv[8]) if len(sys.argv) > 8 else Nb
hsi = float(sys.argv[9]) if len(sys.argv) > 9 else hs
tol = float(sys.argv[10]) if len(sys.argv) > 10 else 1e-10
import os
SST = float(os.environ.get('IPM_SSTART', -20.0))          # origin truncation (the neglected local corrections are O(e^{s_start}))
SST_IN = float(os.environ.get('IPM_SSTART_IN', -20.0))
log = open(f"ipm_deep_{tag}.out", "a")


def out(msg):
    print(msg, flush=True); log.write(msg + "\n"); log.flush()


def dip(B, Y):
    r = B.march(full(B, Y) / B.ea2[:, None], return_all=True)
    D = (1 + B.lam) + r['Ur'][:, 0]; D0 = 1 + B.lam - r['A']
    sel = (B.s > -3) & (B.s < 2); i = np.where(sel)[0][np.argmin(D[sel])]
    kc = i + np.argmax(D[i:] / D0 > 3.0)
    G = np.exp(-B.s) * np.gradient(r['Th'][:, 0], B.s); Om = r['Omega'][:, 0]
    w = slice(max(i - 20, 0), kc + 1)                      # dip-to-front window
    mis = float(np.max(np.abs(G[w] + Om[w])) / np.max(np.abs(G[w])))   # IPM identity Ω_b = −∂ₓR (resolution check)
    return B.s[i], D[i] / D0, D[i], B.s[kc], mis, float(np.max(G[w]))


t0 = time.time()
z, zmax, dzc = 1 / lam0, 1 / lam1, dz
out(f"# ipm_deep {tag}: λ {lam0} -> {lam1}, Δz = {dz}, Nb={Nb} hs={hs} s_start={SST} tol={tol} start={start} (grid in: Nb={nbi} hs={hsi} s_start={SST_IN})")
acc, rows = [], []
import glob as _glob
if os.environ.get('RESUME', '1') == '1' and os.path.exists(f"ipm_deep_{tag}.npy"):
    _done = np.load(f"ipm_deep_{tag}.npy")
    _st = sorted(((1 / float(f.split('_lam')[1][:-4]), f) for f in _glob.glob(f"ipm_deep_{tag}_lam*.npy")))
    if len(_done) and _st:
        rows = [tuple(r) for r in _done]
        acc = [(zz, np.load(ff)) for zz, ff in _st[-2:]]
        z = acc[-1][0] + dz; dzc = dz
        out(f"# resumed after z = {acc[-1][0]:.4f} ({len(rows)} points done)")
Y0 = np.load(start)
if not acc and (nbi, hsi, SST_IN) != (Nb, hs, SST):
    Y0 = transfer(BQLogPolar(lam0, hs=hsi, Nb=nbi, s_start=SST_IN), Y0, IPM(lam0, hs=hs, Nb=Nb, s_sw=12.0, s_start=SST))
while z <= zmax + 1e-12:
    lam = 1 / z
    B = IPM(lam, hs=hs, Nb=Nb, s_sw=12.0, s_start=SST)
    if len(acc) >= 2:
        (za, Ya), (zb, Yb) = acc[-2], acc[-1]
        Yg = Yb + (Yb - Ya) * (z - zb) / (zb - za)
    else:
        Yg = acc[-1][1] if acc else Y0
    ok = False
    try:
        Y, info, ok = newton(B, Yg, tol=tol, maxit=16, verbose=False, t0=t0, fd='central', pert=1e-6)
    except RuntimeError as e:
        out(f"λ={lam:.7f}: invalid predictor ({e})")
    if not ok:
        if not acc:
            out("failed at the start"); break
        dzc *= 0.5
        out(f"λ={lam:.7f}: no convergence (vrmin of last iterate {info['vrmin'] if 'info' in dir() else float('nan'):.4f}); Δz -> {dzc:.2e}")
        if abs(dzc) < 2e-3:
            out("STOP: step below 2e-3"); break
        z = acc[-1][0] + dzc
        continue
    acc.append((z, Y)); acc = acc[-2:]
    sd, dh, dd, sc, mis, gmax = dip(B, Y)
    rows.append((lam, info['A'], info['m'], info['vrmin'], sd, dh, dd, sc, mis, gmax))
    out(f"λ={lam:.7f} z={z:.4f}: m−2={info['m'] - 2:+.4e} D0={1+lam-info['A']:.5f} vrmin={info['vrmin']:.5f} "
        f"dip s={sd:.3f} D̂_dip={dh:.4f} front s={sc:.3f} max∂ₓR={gmax:.2f} |Ω_b+∂ₓR|/max={mis:.1e} t={time.time() - t0:.0f}s")
    np.save(f"ipm_deep_{tag}.npy", np.array(rows))
    np.save(f"ipm_deep_{tag}_lam{lam:.7f}.npy", Y)
    dzc = min(dz, 1.5 * dzc)
    z = z + dzc
