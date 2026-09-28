"""Continuation of the CCF self-similar branch to the sonic (cusp) limit δ → 0 with the layer PINNED.

The scaling symmetry Θ → κΘ(ξ/κ) translates solutions in η.  Instead of fixing the far-field amplitude we use
it to pin the minimum of the sonic factor den = 1+λ+HΘ/ξ at the centre η_c of a geometrically graded grid:
        F_i = φ_i − φ_{i0} − ∫_{η_{i0}}^{η_i}(λ/den − c) dη  (i ≠ i0),   F_{i0} = den_{c+1} − den_{c−1} = 0,
plus the continuation condition den_c = δ.  This removes the layer-translation mode (which made λ, δ and
arclength continuations stall) and δ becomes a regular parameter; the grid is refined geometrically as the
layer width w (∝ δ² near the cusp limit) shrinks.  Records λ(δ), p(δ).
"""
import numpy as np, sys, time
import scipy.linalg as sla
from scipy.interpolate import CubicSpline
from mapped_continuation import layer_info
from dense_arclength import jac
from ccf_sinhgrid import CCFSinh


def Fpin(M, phi, lam, ic):
    den = 1 + lam + M.G(phi)
    F = phi - phi[M.i0] - M.cumint(lam / den - M.c)
    F[M.i0] = den[ic + 1] - den[ic - 1]
    return F, den


def bordered(M, phi, lam, den, ic):
    J, Fl = jac(M, phi, lam, den)
    Psi = np.exp(phi)
    J[M.i0, :] = (M.E[ic + 1] * M.Hm[ic + 1, :] - M.E[ic - 1] * M.Hm[ic - 1, :]) * Psi
    Fl[M.i0] = 0.0
    N = M.N
    A = np.empty((N + 1, N + 1))
    A[:N, :N] = J; A[:N, N] = Fl
    A[N, :N] = M.E[ic] * M.Hm[ic, :] * Psi; A[N, N] = 1.0
    return A


def solve(M, phi, lam, ic, delta, tol=1e-11, maxit=12, verbose=False):
    hist = []
    lu = None
    for it in range(maxit):
        F, den = Fpin(M, phi, lam, ic)
        if den.min() <= 0:
            return phi, lam, False, hist, None
        g = den[ic] - delta
        res = max(np.abs(F).max(), abs(g) / max(delta, 1e-14) * 1e-3)
        hist.append(res)
        if verbose:
            print(f"      it {it}: |F|={np.abs(F).max():.2e} g/δ={g/delta:.2e} lam={lam:.12f} argmin-ic={np.argmin(den)-ic}", flush=True)
        if np.abs(F).max() < tol and abs(g) < 1e-9 * delta:
            return phi, lam, True, hist, lu
        if it >= 3 and res > 0.5 * hist[-2]:
            return phi, lam, False, hist, None
        A = bordered(M, phi, lam, den, ic)
        lu = sla.lu_factor(A, check_finite=False)
        dd = sla.lu_solve(lu, -np.concatenate([F, [g]]), check_finite=False)
        tt = 1.0
        while tt > 1e-4:
            dn = 1 + lam + tt * dd[-1] + M.G(phi + tt * dd[:-1])
            if dn.min() > 0:
                break
            tt *= 0.5
        phi, lam = phi + tt * dd[:-1], lam + tt * dd[-1]
    return phi, lam, False, hist, None


def new_grid(etac, w, pts, hs=0.03):
    M = CCFSinh(30.0, 120.0, hs=hs, eta_d=etac, hc=min(w / pts, hs / 2), B=1.0, c=0.7)
    ic = int(np.argmin(np.abs(M.eta - etac)))
    return M, ic


def main(start, tag, pts=24, fac=0.8, dmin=1e-6):
    d = np.load(start)
    eta0, phi0, lam = d[0], d[1], float(d[2][0])
    # locate the layer on the source grid (fine enough: saved from the sinh grid)
    sp = CubicSpline(eta0, phi0)
    from mapped_continuation import make_grid
    Mt = make_grid(-0.96, 0.0015, hs=0.03, pts=24)
    Ft, dent = Mt.F(sp(Mt.eta), lam, Mt.normval(sp(Mt.eta), lam))
    e, md, w = layer_info(Mt, dent)
    M, ic = new_grid(e, w, pts)
    phi = sp(M.eta)
    F, den = Fpin(M, phi, lam, ic)
    delta = den[ic]
    print(f"start: lam={lam:.10f} layer {e:.6f} w={w:.3e} δ(grid centre)={delta:.5e} N={M.N}", flush=True)
    phi, lam, ok, hist, lu = solve(M, phi, lam, ic, delta, verbose=True)
    print(f"pinned solve ok={ok}: lam={lam:.12f} p={M.p_of(phi, lam):.12f}", flush=True)
    if not ok:
        return
    rows = []
    t0 = time.time()
    tphi, tl = None, None
    while delta > dmin:
        # tangent d(φ,λ)/dδ from the converged factorisation (bordered row e_N)
        F, den = Fpin(M, phi, lam, ic)
        A = bordered(M, phi, lam, den, ic)
        t = sla.solve(A, np.concatenate([np.zeros(M.N), [1.0]]), check_finite=False)
        dnew = delta * fac
        ph, lm = phi + (dnew - delta) * t[:-1], lam + (dnew - delta) * t[-1]
        ph, lm, ok, hist, _ = solve(M, ph, lm, ic, dnew)
        if not ok:
            fac = 1 - (1 - fac) * 0.5
            print(f"   step failed (hist {', '.join(f'{h:.1e}' for h in hist)}); fac -> {fac:.4f}", flush=True)
            if fac > 0.999:
                print("STOP", flush=True)
                break
            continue
        phi, lam, delta = ph, lm, dnew
        F, den = Fpin(M, phi, lam, ic)
        e, md, w = layer_info(M, den)
        p = M.p_of(phi, lam)
        rows.append((delta, lam, p, e, w, M.N, M.hc))
        print(f"δ={delta:.6e} lam={lam:.12f} p={p:.12f} w={w:.3e} w/δ²={w/delta**2:.3f} hc={M.hc:.2e} N={M.N} "
              f"nit={len(hist)-1} t={time.time()-t0:.0f}s", flush=True)
        np.save(f"pin_{tag}_branch.npy", np.array(rows))
        np.save(f"pin_{tag}_last.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
        if len(rows) > 1 and (rows[-2][2] - 2) * (p - 2) < 0:
            print(f"*** p=2 CROSSING between δ={rows[-2][0]:.4e} (λ={rows[-2][1]:.10f}) and δ={delta:.4e} (λ={lam:.10f})", flush=True)
        # regrid when the layer has thinned: keep ~pts points per width at the centre
        if w / M.hc < 0.75 * pts:
            spn = CubicSpline(M.eta, phi)
            M, ic = new_grid(M.eta[ic], w, pts)
            phi = spn(M.eta)
            phi, lam, okr, hist, _ = solve(M, phi, lam, ic, delta)
            print(f"   regrid: hc={M.hc:.2e} N={M.N} ok={okr} hist={', '.join(f'{h:.1e}' for h in hist)} "
                  f"lam={lam:.12f} p={M.p_of(phi, lam):.12f}", flush=True)
            if not okr:
                print("STOP (regrid)", flush=True)
                break
        if len(hist) - 1 <= 3:
            fac = max(fac * fac, 0.5) if fac < 0.97 else fac - 0.05


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
