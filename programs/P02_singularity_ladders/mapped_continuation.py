"""Adaptive-grid pseudo-arclength continuation of the CCF profile family (mapped grid, thin sonic layers).
Re-grids when the internal layer (min of 1+λ+HΘ/ξ) becomes thinner than ~15 local grid spacings or moves.
Detects and refines every p=2 crossing (smooth profiles)."""
import numpy as np, sys, time, json
from ccf_mapped import CCFMapped


def layer_info(M, den):
    j = int(np.argmin(den))
    j = min(max(j, 1), M.N - 2)
    # curvature in η using non-uniform 3-point formula
    e0, e1, e2 = M.eta[j - 1], M.eta[j], M.eta[j + 1]
    d0, d1, d2 = den[j - 1], den[j], den[j + 1]
    kap = 2 * ((d2 - d1) / (e2 - e1) - (d1 - d0) / (e1 - e0)) / (e2 - e0) / 2
    w = np.sqrt(max(d1, 1e-14) / max(kap, 1e-12))
    return M.eta[j], d1, w


def make_grid(etad, w, hs=0.03, pts=20, sigma=2.0, c=0.7, L1=30.0, L2=120.0):
    eps = min(1.0, (w / pts) / hs)
    return CCFMapped(L1, L2, hs=hs, eta_d=etad, eps=eps, sigma=sigma, c=c)


def run(start_file, lam, nsteps, tag, direction=-1.0, ds0=0.01, dsmax=0.05, hs=0.03):
    d = np.load(start_file)
    eta_old, phi_old = d[0], d[1]
    M0 = make_grid(-1.0, 0.05, hs=hs)
    phi0grid = M0.from_other(eta_old, phi_old)
    phi0grid, lam, ok = M0.solve_p(phi0grid, lam, 2.0)
    F, den = M0.F(phi0grid, lam, M0.normval(phi0grid, lam))
    etad, md, w = layer_info(M0, den)
    M = make_grid(etad, w, hs=hs)
    phi = M.from_other(M0.eta, phi0grid)
    phi, lam, ok = M.solve_p(phi, lam, 2.0)
    print(f"start: lam={lam:.12f} ok={ok} N={M.N} eps={M.eps:.3e} layer at {etad:.3f} w={w:.2e} min_den={md:.3e}", flush=True)
    target = M.normval(phi, lam)
    wv = lambda M: (np.abs(M.eta) < 10.0) * M.hs * M.gp / 20.0
    W = wv(M)
    F, den = M.F(phi, lam, target)
    Fl = M.Fl(lam, den)
    tphi, _ = M._gmres(lambda v: M.Jv(v, phi, lam, den), -Fl, M.N)
    tl = 1.0
    nrm = np.sqrt((W * tphi) @ tphi + tl**2)
    tphi, tl = direction * tphi / nrm, direction * tl / nrm
    ds = ds0
    p_prev = M.p_of(phi, lam)
    branch, crossings = [], []
    t0 = time.time()
    step = 0
    while step < nsteps:
        N = M.N
        ph, lm = phi + ds * tphi, lam + ds * tl
        ok = False
        for it in range(15):
            F, den = M.F(ph, lm, target)
            if den.min() <= 0:
                break
            g = (W * tphi) @ (ph - phi) + tl * (lm - lam) - ds
            if max(np.max(np.abs(F)), abs(g)) < 1e-11:
                ok = True
                break
            Fl = M.Fl(lm, den)
            sc = 1.0 / max(np.linalg.norm(W * tphi), abs(tl))

            def mv(x, ph=ph, lm=lm, den=den, Fl=Fl, sc=sc):
                v, s = x[:N], x[N]
                return np.concatenate([M.Jv(v, ph, lm, den) + s * Fl, [sc * ((W * tphi) @ v + tl * s)]])
            dd, info = M._gmres(mv, -np.concatenate([F, [sc * g]]), N + 1, tol=1e-10)
            ph, lm = ph + dd[:N], lm + dd[N]
        if not ok:
            ds *= 0.5
            if ds < 1e-8:
                print("   ds underflow; stop", flush=True)
                break
            continue
        Fl = M.Fl(lm, den)
        sc = 1.0 / max(np.linalg.norm(W * tphi), abs(tl))

        def mvt(x, ph=ph, lm=lm, den=den, Fl=Fl, sc=sc):
            v, s = x[:N], x[N]
            return np.concatenate([M.Jv(v, ph, lm, den) + s * Fl, [sc * ((W * tphi) @ v + tl * s)]])
        t, info = M._gmres(mvt, np.concatenate([np.zeros(N), [sc]]), N + 1, tol=1e-10)
        nrm = np.sqrt((W * t[:N]) @ t[:N] + t[N] ** 2)
        tphi_n, tl_n = t[:N] / nrm, t[N] / nrm
        p_new = M.p_of(ph, lm)
        if (p_prev - 2.0) * (p_new - 2.0) < 0:
            sp, sl = (phi, lam) if abs(p_prev - 2) < abs(p_new - 2) else (ph, lm)
            rphi, rlam, rok = M.solve_p(sp.copy(), sl, 2.0, phi0=target)
            _, rden = M.F(rphi, rlam, target)
            e_, m_, w_ = layer_info(M, rden)
            crossings.append(dict(lam=float(rlam), ok=bool(rok), min_den=float(m_), layer_eta=float(e_),
                                  layer_w=float(w_), N=int(M.N), eps=float(M.eps), dlam_ds=float(tl_n)))
            np.save(f"mc_{tag}_cross{len(crossings)}.npy", np.vstack([M.eta, rphi, np.full(M.N, rlam)]))
            print(f"*** CROSSING {len(crossings)}: lam={rlam:.12f} ok={rok} min_den={m_:.3e} w={w_:.2e} "
                  f"(λ {lam:.6f}→{lm:.6f}, p {p_prev:.6f}→{p_new:.6f}) N={M.N} eps={M.eps:.2e}", flush=True)
        phi, lam, tphi, tl, p_prev = ph, lm, tphi_n, tl_n, p_new
        etad, md, w = layer_info(M, den)
        branch.append((lam, p_new, md, etad, w, tl))
        if step % 5 == 0:
            print(f"step {step}: lam={lam:.10f} p={p_new:.10f} min_den={md:.3e} layer={etad:.4f} w={w:.2e} "
                  f"h_loc={M.hs*M.eps:.1e} dlam/ds={tl:+.4f} ds={ds:.2e} N={M.N} t={time.time()-t0:.0f}s", flush=True)
        # re-grid if layer under-resolved (<12 pts) or over-resolved (>60 pts) or moved
        hloc = M.hs * M.eps
        if (w / hloc < 12 and M.eps > 1e-4) or (w / hloc > 60 and M.eps < 1) or abs(etad - M.eta_d) > 0.3 * M.sigma * max(M.eps, 0.05):
            Mn = make_grid(etad, w, hs=M.hs)
            phi = Mn.from_other(M.eta, phi)
            tphi = Mn.from_other(M.eta, tphi)
            M = Mn
            W = wv(M)
            phi, okr = M.solve_fixed(phi, lam, phi0=target)
            print(f"   regrid: layer {etad:.4f} w={w:.2e} -> eps={M.eps:.2e} N={M.N} resolve_ok={okr}", flush=True)
            if not okr:
                ds *= 0.5
        if it < 5:
            ds = min(ds * 1.3, dsmax)
        step += 1
        if md < 1e-9:
            print("   sonic limit reached; stop", flush=True)
            break
    np.save(f"mc_{tag}_branch.npy", np.array(branch))
    json.dump(crossings, open(f"mc_{tag}_crossings.json", "w"), indent=1)


if __name__ == "__main__":
    run(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], float(sys.argv[5]) if len(sys.argv) > 5 else -1.0)
