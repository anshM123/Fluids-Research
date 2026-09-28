"""Pseudo-arclength continuation of the CCF self-similar family through the near-sonic fold region, on the
geometrically graded grid of ccf_sinhgrid.py, with dense LU solves.

Why this parametrisation: near the fold beyond λ₂ the null vector of the fixed-λ Jacobian is (almost) a pure
translation of the thin sonic layer.  Neither λ nor the sonic depth δ=min den is a good parameter there (δ is
translation invariant to first order: the δ-bordered Jacobian has σ_min ≈ 2e-4), and an arclength norm
weighted by dη (the earlier dense_arclength.py) gives the layer almost zero weight.  Here the arclength norm
weights every grid point in |η|<10 equally (weight h_s per point, i.e. uniform in the computational
variable s), so the layer translation is fully visible.

Records (λ, p, δ, layer position, width) along the branch, detects p=2 crossings and refines them with a
dense bordered Newton solve for (φ, λ) at p=2.
"""
import numpy as np, sys, time, json
import scipy.linalg as sla
from scipy.interpolate import CubicSpline
from mapped_continuation import layer_info
from dense_arclength import jac
from ccf_sinhgrid import make_sinh_grid


def weights(M):
    return (np.abs(M.eta) < 10.0) * M.hs / 20.0


def solve_bordered(M, ph, lm, target, row_phi, row_lam, rhs_fun, tol=1e-11, maxit=10, verbose=False):
    """Newton for F(φ,λ)=0 plus one scalar constraint g(φ,λ)=rhs_fun(φ,λ) with constant derivative row.
    Returns (φ, λ, ok, lu) with lu the last factorisation (for tangent solves)."""
    N = M.N
    hist = []
    lu = None
    for it in range(maxit):
        F, den = M.F(ph, lm, target)
        if den.min() <= 0:
            return ph, lm, False, None, hist
        g = rhs_fun(ph, lm)
        res = max(np.abs(F).max(), abs(g))
        hist.append(res)
        if verbose:
            print(f"      it {it}: |F|={np.abs(F).max():.2e} g={g:.2e} lam={lm:.12f} minden={den.min():.4e}", flush=True)
        if res < tol:
            return ph, lm, True, lu, hist
        if it >= 2 and res > 0.5 * hist[-2]:       # not contracting: give up (caller reduces step)
            return ph, lm, False, None, hist
        J, Fl = jac(M, ph, lm, den)
        A = np.empty((N + 1, N + 1))
        A[:N, :N] = J; A[:N, N] = Fl; A[N, :N] = row_phi; A[N, N] = row_lam
        lu = sla.lu_factor(A, check_finite=False)
        dd = sla.lu_solve(lu, -np.concatenate([F, [g]]), check_finite=False)
        tt = 1.0
        while tt > 1e-3:
            dn = 1 + (lm + tt * dd[N]) + M.G(ph + tt * dd[:N])
            if dn.min() > 0:
                break
            tt *= 0.5
        ph, lm = ph + tt * dd[:N], lm + tt * dd[N]
    return ph, lm, False, None, hist


def tangent(M, ph, lm, target, tphi_prev, tl_prev, W):
    F, den = M.F(ph, lm, target)
    J, Fl = jac(M, ph, lm, den)
    N = M.N
    A = np.empty((N + 1, N + 1))
    A[:N, :N] = J; A[:N, N] = Fl; A[N, :N] = W * tphi_prev; A[N, N] = tl_prev
    t = sla.solve(A, np.concatenate([np.zeros(N), [1.0]]), check_finite=False)
    nrm = np.sqrt((W * t[:N]) @ t[:N] + t[N] ** 2)
    return t[:N] / nrm, t[N] / nrm


def refine_p2(M, ph, lm, target, p=2.0):
    """bordered Newton for (φ,λ) with p(φ,λ)=p, i.e. h1(φ) = λ/p - 1 - λ"""
    N = M.N
    for it in range(15):
        F, den = M.F(ph, lm, target)
        g = M.h1(ph) - (lm / p - 1 - lm)
        if max(np.abs(F).max(), abs(g)) < 1e-12:
            return ph, lm, True
        J, Fl = jac(M, ph, lm, den)
        A = np.empty((N + 1, N + 1))
        A[:N, :N] = J; A[:N, N] = Fl; A[N, :N] = M.dh1(ph); A[N, N] = -(1 / p - 1)
        dd = sla.solve(A, -np.concatenate([F, [g]]), check_finite=False)
        ph, lm = ph + dd[:N], lm + dd[N]
    return ph, lm, False


def regrid(M, phi, tphi, lam, target, pts, hs):
    F, den = M.F(phi, lam, target)
    e, md, w = layer_info(M, den)
    Mn = make_sinh_grid(e, w, hs=hs, pts=pts)
    sp = CubicSpline(M.eta, phi); st = CubicSpline(M.eta, tphi)
    return Mn, sp(Mn.eta), st(Mn.eta)


def run(start, tag, direction=-1.0, ds=2e-3, nsteps=400, pts=24, hs=0.03, dsmax=0.05):
    d = np.load(start)
    eta0, phi0, lam = d[0], d[1], float(d[2][0])
    # provisional grid -> layer -> proper grid
    from mapped_continuation import make_grid
    M0 = make_grid(-0.963, 0.0015, hs=hs, pts=24)
    target = M0.normval(M0.from_other(eta0, phi0), lam)
    F, den = M0.F(M0.from_other(eta0, phi0), lam, target)
    e, md, w = layer_info(M0, den)
    M = make_sinh_grid(e, w, hs=hs, pts=pts)
    phi = CubicSpline(eta0, phi0)(M.eta)
    for it in range(25):
        F, den = M.F(phi, lam, target)
        if np.abs(F).max() < 1e-11:
            break
        J, Fl = jac(M, phi, lam, den)
        phi = phi - sla.solve(J, F, check_finite=False)
    F, den = M.F(phi, lam, target)
    e, md, w = layer_info(M, den)
    print(f"start lam={lam:.10f} p={M.p_of(phi, lam):.10f} |F|={np.abs(F).max():.1e} δ={md:.5e} layer={e:.5f} "
          f"w={w:.3e} N={M.N}", flush=True)
    W = weights(M)
    # initial tangent: d/dλ direction (fixed-λ Jacobian is nonsingular at the start)
    J, Fl = jac(M, phi, lam, den)
    tphi = -sla.solve(J, Fl, check_finite=False); tl = 1.0
    nrm = np.sqrt((W * tphi) @ tphi + tl ** 2)
    tphi, tl = direction * tphi / nrm, direction * tl / nrm
    rows, crossings = [], []
    p_prev = M.p_of(phi, lam)
    t0 = time.time()
    step = 0
    while step < nsteps:
        ph_pred, lm_pred = phi + ds * tphi, lam + ds * tl
        row_phi, row_lam = W * tphi, tl
        cons = lambda ph, lm: (W * tphi) @ (ph - ph_pred) + tl * (lm - lm_pred)
        ph, lm, ok, lu, hist = solve_bordered(M, ph_pred, lm_pred, target, row_phi, row_lam, cons)
        if not ok:
            ds *= 0.5
            print(f"   corrector failed (hist {', '.join(f'{h:.1e}' for h in hist)}); ds -> {ds:.2e}", flush=True)
            if ds < 1e-8:
                print("STOP: ds too small", flush=True)
                break
            continue
        nit = len(hist) - 1
        tphi_new, tl_new = tangent(M, ph, lm, target, tphi, tl, W)
        phi, lam, tphi, tl = ph, lm, tphi_new, tl_new
        F, den = M.F(phi, lam, target)
        e, md, w = layer_info(M, den)
        p = M.p_of(phi, lam)
        rows.append((lam, p, md, e, w, tl, M.N))
        print(f"step {step}: lam={lam:.10f} p={p:.10f} δ={md:.5e} layer={e:.5f} w={w:.3e} dλ/ds={tl:+.4f} ds={ds:.2e} "
              f"nit={nit} N={M.N} t={time.time()-t0:.0f}s", flush=True)
        np.save(f"sac_{tag}_branch.npy", np.array(rows))
        np.save(f"sac_{tag}_last.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
        if (p_prev - 2) * (p - 2) < 0:
            print(f"*** p=2 CROSSING between λ={rows[-2][0]:.10f} and λ={lam:.10f}; refining", flush=True)
            # interpolate linearly in arclength between the two states is unnecessary: start from current state
            ph2, lm2, ok2 = refine_p2(M, phi.copy(), lam, target)
            F2, den2 = M.F(ph2, lm2, target)
            e2, md2, w2 = layer_info(M, den2)
            crossings.append(dict(step=step, lam=float(lm2), ok=bool(ok2), p=float(M.p_of(ph2, lm2)),
                                  delta=float(md2), layer=float(e2), w=float(w2), N=int(M.N)))
            np.save(f"sac_{tag}_cross{len(crossings)}.npy", np.vstack([M.eta, ph2, np.full(M.N, lm2)]))
            print(f"*** refined smooth profile: λ={lm2:.12f} ok={ok2} δ={md2:.4e} layer={e2:.5f}", flush=True)
            json.dump(crossings, open(f"sac_{tag}_crossings.json", "w"), indent=1)
        p_prev = p
        # re-grid when the layer drifts by > w/3 from the grid centre or thins/thickens by > 25 %
        wc = M.hc * pts
        if abs(e - M.eta_d) > w / 3 or not (0.8 < w / wc < 1.25):
            M, phi, tphi = regrid(M, phi, tphi, lam, target, pts, hs)
            W = weights(M)
            nrm = np.sqrt((W * tphi) @ tphi + tl ** 2); tphi, tl = tphi / nrm, tl / nrm
            # re-converge on the new grid, keeping the arclength position (orthogonal projection)
            ph0, lm0 = phi.copy(), lam
            cons0 = lambda ph, lm: (W * tphi) @ (ph - ph0) + tl * (lm - lm0)
            phi, lam, okr, lu, hist = solve_bordered(M, phi, lam, target, W * tphi, tl, cons0, maxit=15)
            print(f"   regrid: centre {M.eta_d:.5f}, hc={M.hc:.2e}, N={M.N}, ok={okr}, hist={', '.join(f'{h:.1e}' for h in hist)}",
                  flush=True)
            if not okr:
                print("STOP: regrid re-solve failed", flush=True)
                break
            tphi, tl = tangent(M, phi, lam, target, tphi, tl, W)
        if nit <= 4:
            ds = min(ds * 1.3, dsmax)
        elif nit >= 7:
            ds *= 0.7
        step += 1
    json.dump(crossings, open(f"sac_{tag}_crossings.json", "w"), indent=1)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else -1.0,
        ds=float(sys.argv[4]) if len(sys.argv) > 4 else 2e-3)
