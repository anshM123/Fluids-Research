"""Continuation of the CCF self-similar branch toward the sonic limit, parametrised by the sonic depth
δ = (1+λ+HΘ/ξ)(η_p) at a probe point η_p placed at the internal layer (updated on re-gridding).
Unknowns (φ, λ); equations F(φ,λ)=0 and den(η_p) = δ. Dense bordered Newton on the adaptive mapped grid.
Records p(δ), λ(δ): decides whether p returns to 2 (a further smooth profile) before δ→0."""
import numpy as np, sys, time
import scipy.linalg as sla
from ccf_mapped import CCFMapped
from mapped_continuation import layer_info, make_grid
from dense_arclength import jac

start = sys.argv[1] if len(sys.argv) > 1 else 'pdiag_start.npy'
d = np.load(start); eta0, phi0, lam = d[0], d[1], float(d[2][0])
M = make_grid(-0.963, 0.0015, hs=0.03, pts=24)
phi = M.from_other(eta0, phi0); target = M.normval(phi, lam)
for it in range(20):
    F, den = M.F(phi, lam, target)
    if np.max(np.abs(F)) < 1e-11:
        break
    J, Fl = jac(M, phi, lam, den); phi = phi - sla.solve(J, F)
F, den = M.F(phi, lam, target)
ip = int(np.argmin(den)); delta = den[ip]
print(f"start lam={lam:.10f} p={M.p_of(phi, lam):.10f} delta={delta:.5e} probe eta={M.eta[ip]:.4f} N={M.N}", flush=True)
fac = 0.93
rows = []
t0 = time.time()
while delta > 1e-7:
    dt = delta * fac
    ph, lm, ok = phi.copy(), lam, False
    for it in range(15):
        F, den = M.F(ph, lm, target)
        if den.min() <= 0:
            break
        gcon = den[ip] - dt
        res = max(np.abs(F).max(), abs(gcon) / max(dt, 1e-12) * 1e-3)
        if np.abs(F).max() < 1e-10 and abs(gcon) < 1e-9 * max(dt, 1e-6):
            ok = True
            break
        J, Fl = jac(M, ph, lm, den)
        N = M.N
        A = np.zeros((N + 1, N + 1)); A[:N, :N] = J; A[:N, N] = Fl
        # d den(η_p)/dφ_j = E_p Hm[p, j] Ψ_j ;  d den(η_p)/dλ = 1
        A[N, :N] = M.E[ip] * M.Hm[ip, :] * np.exp(ph); A[N, N] = 1.0
        dd = sla.solve(A, -np.concatenate([F, [gcon]]))
        tt = 1.0
        while tt > 1e-3:
            Fn, dn = M.F(ph + tt * dd[:N], lm + tt * dd[N], target)
            if dn.min() > 0:
                break
            tt *= 0.5
        ph, lm = ph + tt * dd[:N], lm + tt * dd[N]
    if not ok:
        fac = 1 - (1 - fac) * 0.5
        print(f"   step failed (target δ={dt:.3e}); fac -> {fac:.4f}", flush=True)
        if fac > 0.9999:
            print("STOP", flush=True)
            break
        continue
    phi, lam, delta = ph, lm, dt
    F, den = M.F(phi, lam, target)
    e, md, w = layer_info(M, den)
    p = M.p_of(phi, lam)
    rows.append((delta, lam, p, md, e, w))
    print(f"δ={delta:.4e} lam={lam:.10f} p={p:.10f} min_den={md:.4e} layer={e:.4f} w={w:.2e} "
          f"pts/w={w/(M.hs*M.eps):.1f} N={M.N} it={it} t={time.time()-t0:.0f}s", flush=True)
    np.save("delta_branch.npy", np.array(rows))
    np.save("delta_last.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
    if len(rows) > 1 and (rows[-2][2] - 2) * (p - 2) < 0:
        print(f"*** p=2 crossing between δ={rows[-2][0]:.4e} (λ={rows[-2][1]:.10f}) and δ={delta:.4e} (λ={lam:.10f})", flush=True)
    # re-grid when the layer gets under-resolved or the probe drifts off the minimum
    if w / (M.hs * M.eps) < 14 or abs(int(np.argmin(den)) - ip) > 3:
        Mn = make_grid(e, w, hs=0.03, pts=24)
        phi = Mn.from_other(M.eta, phi); M = Mn
        for it2 in range(20):
            F, den = M.F(phi, lam, target)
            if np.abs(F).max() < 1e-11:
                break
            J, Fl = jac(M, phi, lam, den); phi = phi - sla.solve(J, F)
        ip = int(np.argmin(den)); delta = den[ip]
        print(f"   regrid: eps={M.eps:.2e} N={M.N} |F|={np.abs(F).max():.1e} new probe δ={delta:.4e} at {M.eta[ip]:.4f}", flush=True)
    if it <= 3:
        fac = max(fac * fac, 0.8) if fac < 0.95 else fac * 0.99
