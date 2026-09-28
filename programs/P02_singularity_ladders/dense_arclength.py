"""Robust pseudo-arclength continuation of the CCF profile family on the adaptive mapped grid with DENSE
Jacobians and direct (LU) solves — for the fold region beyond λ₂ where Krylov solvers stall.
Tracks p(λ) and detects/refines p=2 crossings (smooth profiles)."""
import numpy as np, sys, time, json
import scipy.linalg as sla
from ccf_mapped import CCFMapped
from mapped_continuation import layer_info, make_grid


def cumint_mat(M, X):
    """apply M.cumint column-wise to matrix X (N×K)"""
    Y = X * M.gp[:, None]
    n, h, w = M.N, M.hs, M.w8
    mid = np.zeros((n - 1, X.shape[1]))
    for j, off in enumerate(range(-3, 5)):
        mid[3:n - 4] += w[j] * Y[3 + off:n - 4 + off]
    for i in list(range(0, 3)) + list(range(n - 4, n - 1)):
        s = i if i < 3 else 7 - (n - 1 - i)
        mid[i] = M.wl[s] @ Y[i - s:i - s + 8]
    C = np.zeros_like(X)
    C[1:] = np.cumsum(mid * h, axis=0)
    return C - C[M.i0]


def jac(M, phi, lam, den):
    N = M.N
    Psi = np.exp(phi)
    dG = (M.E[:, None] * M.Hm) * Psi[None, :]
    J = np.eye(N) - cumint_mat(M, (-lam / den**2)[:, None] * dG)
    J[:, M.i0] -= 1.0
    J[M.i0, :] = 0.0
    J[M.i0, M.iR] = 1.0
    Fl = M.Fl(lam, den)
    return J, Fl


def run(start_file, nsteps, tag, direction=-1.0, ds=0.004, dsmax=0.02):
    d = np.load(start_file)
    eta0, phi0_, lam = d[0], d[1], float(d[2][0])
    F_dummy = None
    M0 = make_grid(-0.963, 0.002, hs=0.03)
    phi = M0.from_other(eta0, phi0_)
    target = M0.normval(phi, lam)
    # re-solve at fixed λ on this grid (dense Newton)
    M = M0
    for it in range(20):
        F, den = M.F(phi, lam, target)
        if np.max(np.abs(F)) < 1e-11:
            break
        J, Fl = jac(M, phi, lam, den)
        phi = phi - sla.solve(J, F)
    F, den = M.F(phi, lam, target)
    e, md, w = layer_info(M, den)
    print(f"start lam={lam:.8f} p={M.p_of(phi, lam):.8f} |F|={np.abs(F).max():.1e} min_den={md:.3e} w={w:.2e} N={M.N}", flush=True)
    W = (np.abs(M.eta) < 10.0) * M.hs * M.gp / 20.0
    J, Fl = jac(M, phi, lam, den)
    tphi = -sla.solve(J, Fl); tl = 1.0
    nrm = np.sqrt((W * tphi) @ tphi + tl**2)
    tphi, tl = direction * tphi / nrm, direction * tl / nrm
    p_prev = M.p_of(phi, lam)
    rows, crossings = [], []
    t0 = time.time()
    step = 0
    while step < nsteps:
        N = M.N
        ph, lm = phi + ds * tphi, lam + ds * tl
        ok = False
        for it in range(14):
            F, den = M.F(ph, lm, target)
            if den.min() <= 0:
                break
            g = (W * tphi) @ (ph - phi) + tl * (lm - lam) - ds
            res = max(np.max(np.abs(F)), abs(g))
            if res < 1e-10:
                ok = True
                break
            J, Fl = jac(M, ph, lm, den)
            A = np.zeros((N + 1, N + 1))
            A[:N, :N] = J; A[:N, N] = Fl; A[N, :N] = W * tphi; A[N, N] = tl
            dd = sla.solve(A, -np.concatenate([F, [g]]))
            tt = 1.0
            while tt > 1e-3:          # damped update keeping the sonic factor positive
                Fn, dn = M.F(ph + tt * dd[:N], lm + tt * dd[N], target)
                if dn.min() > 0 and np.max(np.abs(Fn)) < 2 * res + 1e-9:
                    break
                tt *= 0.5
            ph, lm = ph + tt * dd[:N], lm + tt * dd[N]
            if step < 3:
                print(f"      corr it {it}: res={res:.2e} tt={tt:.3f} minden={dn.min():.3e}", flush=True)
        if not ok:
            ds *= 0.5
            print(f"   corrector failed; ds -> {ds:.2e}", flush=True)
            if ds < 1e-7:
                break
            continue
        J, Fl = jac(M, ph, lm, den)
        A = np.zeros((N + 1, N + 1))
        A[:N, :N] = J; A[:N, N] = Fl; A[N, :N] = W * tphi; A[N, N] = tl
        t = sla.solve(A, np.concatenate([np.zeros(N), [1.0]]))
        nrm = np.sqrt((W * t[:N]) @ t[:N] + t[N] ** 2)
        tphi, tl = t[:N] / nrm, t[N] / nrm
        phi, lam = ph, lm
        p_new = M.p_of(phi, lam)
        e, md, w = layer_info(M, den)
        rows.append((lam, p_new, md, e, w, tl))
        print(f"step {step}: lam={lam:.10f} p={p_new:.10f} min_den={md:.3e} layer={e:.4f} w={w:.2e} pts/w={w/(M.hs*M.eps):.1f} "
              f"dlam/ds={tl:+.4f} ds={ds:.2e} N={N} t={time.time()-t0:.0f}s", flush=True)
        if (p_prev - 2) * (p_new - 2) < 0:
            crossings.append(dict(step=step, lam_left=float(rows[-2][0]) if len(rows) > 1 else None, lam=float(lam), p=float(p_new)))
            np.save(f"dac_{tag}_cross{len(crossings)}.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
            print(f"*** CROSSING p=2 between previous and λ={lam:.10f}", flush=True)
        p_prev = p_new
        np.save(f"dac_{tag}_last.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
        # regrid if the layer becomes under-resolved
        if w / (M.hs * M.eps) < 12:
            Mn = make_grid(e, w, hs=0.03)
            phi = Mn.from_other(M.eta, phi); tphi = Mn.from_other(M.eta, tphi); M = Mn
            W = (np.abs(M.eta) < 10.0) * M.hs * M.gp / 20.0
            for it in range(15):
                F, den = M.F(phi, lam, target)
                if np.max(np.abs(F)) < 1e-11:
                    break
                J, Fl = jac(M, phi, lam, den)
                phi = phi - sla.solve(J, F)
            print(f"   regrid -> eps={M.eps:.2e}, |F|={np.abs(F).max():.1e}", flush=True)
        if it < 4:
            ds = min(ds * 1.3, dsmax)
        step += 1
    np.save(f"dac_{tag}_branch.npy", np.array(rows))
    json.dump(crossings, open(f"dac_{tag}_crossings.json", "w"), indent=1)


if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]), sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else -1.0)
