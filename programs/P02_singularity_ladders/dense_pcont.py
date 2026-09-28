"""Continuation in the local exponent p (λ unknown) with dense bordered Newton on the mapped grid."""
import numpy as np, sys, time
import scipy.linalg as sla
from ccf_mapped import CCFMapped
from mapped_continuation import layer_info, make_grid
from dense_arclength import jac

d = np.load('pdiag_start.npy'); eta0, phi0, lam = d[0], d[1], float(d[2][0])
M = make_grid(-0.963, 0.002, hs=0.03)
phi = M.from_other(eta0, phi0); target = M.normval(phi, lam)
for it in range(20):
    F, den = M.F(phi, lam, target)
    if np.max(np.abs(F)) < 1e-11: break
    J, Fl = jac(M, phi, lam, den); phi = phi - sla.solve(J, F)
p = M.p_of(phi, lam)
print(f"start lam={lam:.10f} p={p:.10f}", flush=True)
dp = -2e-4
t0 = time.time()
rows = []
while p > 1.95:
    pt = p + dp
    ph, lm, ok = phi.copy(), lam, False
    for it in range(15):
        F, den = M.F(ph, lm, target)
        if den.min() <= 0: break
        g = M.h1(ph) - (lm/pt - 1 - lm)
        res = max(np.abs(F).max(), abs(g))
        if res < 1e-10: ok = True; break
        J, Fl = jac(M, ph, lm, den)
        N = M.N
        A = np.zeros((N+1, N+1)); A[:N,:N] = J; A[:N,N] = Fl
        A[N,:N] = M.dh1(ph); A[N,N] = -(1/pt - 1)
        dd = sla.solve(A, -np.concatenate([F, [g]]))
        tt = 1.0
        while tt > 1e-3:
            Fn, dn = M.F(ph + tt*dd[:N], lm + tt*dd[N], target)
            gn = M.h1(ph + tt*dd[:N]) - ((lm+tt*dd[N])/pt - 1 - (lm+tt*dd[N]))
            if dn.min() > 0 and max(np.abs(Fn).max(), abs(gn)) < 2*res + 1e-9: break
            tt *= 0.5
        ph, lm = ph + tt*dd[:N], lm + tt*dd[N]
    if not ok:
        dp *= 0.5
        print(f"   p-step failed at target {pt:.8f}; dp -> {dp:.2e}", flush=True)
        if abs(dp) < 1e-8: print("STOP: cannot continue in p (turning point in p?)", flush=True); break
        continue
    phi, lam, p = ph, lm, pt
    F, den = M.F(phi, lam, target); e, md, w = layer_info(M, den)
    rows.append((p, lam, md, e, w))
    print(f"p={p:.8f} lam={lam:.10f} min_den={md:.4e} layer={e:.4f} w={w:.2e} pts/w={w/(M.hs*M.eps):.1f} N={M.N} it={it} t={time.time()-t0:.0f}s", flush=True)
    np.save("dpc_last.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
    if (p - 2.0) * (p - dp - 2.0) <= 0 or abs(p - 2.0) < 1e-12:
        print(f"*** p=2 reached: lam={lam:.12f}", flush=True)
    if w/(M.hs*M.eps) < 12:
        Mn = make_grid(e, w, hs=0.03); phi = Mn.from_other(M.eta, phi); M = Mn
        for it2 in range(15):
            F, den = M.F(phi, lam, target)
            if np.abs(F).max() < 1e-11: break
            J, Fl = jac(M, phi, lam, den); phi = phi - sla.solve(J, F)
        print(f"   regrid eps={M.eps:.2e} |F|={np.abs(F).max():.1e}", flush=True)
    if it <= 3: dp = max(dp*1.5, -2e-3)
np.save("dpc_branch.npy", np.array(rows))
