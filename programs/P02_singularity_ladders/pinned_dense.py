"""Dense sampling of the pinned branch in δ (fixed ratio per step, no step growth) for a model-free location
of the extrema of p(δ) and λ(δ).  Every step is re-gridded to exactly `pts` points per layer width, so all
recorded values are fully resolved.  Usage: pinned_dense.py START.npy TAG FAC DMIN [pts] [hs]"""
import numpy as np, sys, time
import scipy.linalg as sla
from scipy.interpolate import CubicSpline
from mapped_continuation import layer_info, make_grid
import pinned_delta as P

start, tag, fac, dmin = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
pts = int(sys.argv[5]) if len(sys.argv) > 5 else 24
P.HS = float(sys.argv[6]) if len(sys.argv) > 6 else 0.03
d = np.load(start)
eta0, phi0, lam = d[0], d[1], float(d[2][0])
sp = CubicSpline(eta0, phi0)
Mt = make_grid(-0.96, 0.0015, hs=0.03, pts=24)
Ft, dent = Mt.F(sp(Mt.eta), lam, Mt.normval(sp(Mt.eta), lam))
e, md, w = layer_info(Mt, dent)
M, ic = P.new_grid(e, w, pts)
phi = sp(M.eta)
F, den = P.Fpin(M, phi, lam, ic)
delta = den[ic]
phi, lam, ok, hist, _ = P.solve(M, phi, lam, ic, delta)
print(f"start δ={delta:.6e} λ={lam:.12f} p={M.p_of(phi, lam):.12f} ok={ok}", flush=True)
rows = [(delta, lam, M.p_of(phi, lam), w, M.N)]
t0 = time.time()
while (fac < 1 and delta > dmin) or (fac > 1 and delta < dmin):
    F, den = P.Fpin(M, phi, lam, ic)
    A = P.bordered(M, phi, lam, den, ic)
    t = sla.solve(A, np.concatenate([np.zeros(M.N), [1.0]]), check_finite=False)
    dn = delta * fac
    ph, lm = phi + (dn - delta) * t[:-1], lam + (dn - delta) * t[-1]
    ph, lm, ok, hist, _ = P.solve(M, ph, lm, ic, dn)
    if not ok:
        print(f"   step failed at δ={dn:.4e}", flush=True)
        break
    phi, lam, delta = ph, lm, dn
    F, den = P.Fpin(M, phi, lam, ic)
    e, md, w = layer_info(M, den)
    # re-grid to exactly pts points per width at every step (values below are on the fresh grid)
    spn = CubicSpline(M.eta, phi)
    M, ic = P.new_grid(M.eta[ic], w, pts)
    phi = spn(M.eta)
    phi, lam, okr, hist, _ = P.solve(M, phi, lam, ic, delta)
    if not okr:
        print("   regrid solve failed", flush=True)
        break
    p = M.p_of(phi, lam)
    rows.append((delta, lam, p, w, M.N))
    print(f"δ={delta:.6e} λ={lam:.13f} p={p:.13f} w/δ²={w/delta**2:.4f} N={M.N} t={time.time()-t0:.0f}s", flush=True)
    np.save(f"pdense_{tag}.npy", np.array(rows))
    np.save(f"pdense_{tag}_last.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
    np.save(f"pdense_states/{tag}_{len(rows):03d}.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
