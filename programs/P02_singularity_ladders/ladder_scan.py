"""Arclength continuation of the CCF profile family with on-the-fly detection and Newton refinement
of every crossing p(λ)=p_target (smooth self-similar profiles). Saves refined profiles."""
import numpy as np, sys, time, json
from ccf_nk import CCFNK
from scipy.sparse.linalg import LinearOperator, gmres


def run(N=16384, c=0.7, L1=30.0, L2=120.0, nsteps=1500, ptarget=2.0, dsmax=0.1, tag="run", lam_start=0.6, direction=-1.0):
    S = CCFNK(L1, L2, N, c)
    if N > 16384:
        # initialise by interpolating a converged coarse solution (Picard is slow at large N)
        S0 = CCFNK(L1, L2, 16384, c)
        p0 = S0.guess(lam_start, 2.0)
        p0, ok = S0.picard(p0, lam_start, tol=1e-6)
        p0, ok = S0.solve_fixed(p0, lam_start)
        phi = np.interp(S.eta, S0.eta, p0)
    else:
        phi = S.guess(lam_start, 2.0)
        phi, ok = S.picard(phi, lam_start, tol=1e-6)
    phi, ok = S.solve_fixed(phi, lam_start)
    print("init ok", ok, flush=True)
    phi0 = S.normval(phi, lam_start)
    lam = lam_start
    w = (np.abs(S.eta) < 10.0) * S.h / 20.0
    F, den = S.F(phi, lam, phi0)
    Fl = S.Fl(lam, den)
    tphi, _ = S._gmres(lambda v: S.Jv(v, phi, lam, den), -Fl, N)
    tl = 1.0
    nrm = np.sqrt((w * tphi) @ tphi + tl**2)
    tphi, tl = direction * tphi / nrm, direction * tl / nrm   # direction=-1: decreasing λ first
    ds = 0.02
    branch, crossings = [], []
    p_prev = S.p_of(phi, lam)
    t0 = time.time()
    step = 0
    while step < nsteps:
        ph, lm = phi + ds * tphi, lam + ds * tl
        ok = False
        for it in range(20):
            F, den = S.F(ph, lm, phi0)
            if den.min() <= 0:
                break
            g = (w * tphi) @ (ph - phi) + tl * (lm - lam) - ds
            if max(np.max(np.abs(F)), abs(g)) < 1e-11:
                ok = True
                break
            Fl = S.Fl(lm, den)
            sc = 1.0 / max(np.linalg.norm(w * tphi), abs(tl))

            def mv(x, ph=ph, lm=lm, den=den, Fl=Fl, sc=sc):
                v, s = x[:N], x[N]
                return np.concatenate([S.Jv(v, ph, lm, den) + s * Fl, [sc * ((w * tphi) @ v + tl * s)]])
            d, info = S._gmres(mv, -np.concatenate([F, [sc * g]]), N + 1, tol=1e-10)
            ph, lm = ph + d[:N], lm + d[N]
        if not ok:
            ds *= 0.5
            if ds < 1e-7:
                print("   step size underflow; stop", flush=True)
                break
            continue
        Fl = S.Fl(lm, den)
        sc = 1.0 / max(np.linalg.norm(w * tphi), abs(tl))

        def mvt(x, ph=ph, lm=lm, den=den, Fl=Fl, sc=sc):
            v, s = x[:N], x[N]
            return np.concatenate([S.Jv(v, ph, lm, den) + s * Fl, [sc * ((w * tphi) @ v + tl * s)]])
        t, info = S._gmres(mvt, np.concatenate([np.zeros(N), [sc]]), N + 1, tol=1e-10)
        nrm = np.sqrt((w * t[:N]) @ t[:N] + t[N] ** 2)
        tphi_new, tl_new = t[:N] / nrm, t[N] / nrm
        p_new = S.p_of(ph, lm)
        # crossing detection
        if (p_prev - ptarget) * (p_new - ptarget) < 0:
            # refine with λ unknown and p fixed, starting from the endpoint closer to target
            start_phi, start_lam = (phi, lam) if abs(p_prev - ptarget) < abs(p_new - ptarget) else (ph, lm)
            rphi, rlam, rok = S.solve_p(start_phi.copy(), start_lam, ptarget, phi0=phi0)
            _, rden = S.F(rphi, rlam, phi0)
            info_c = dict(step=step, lam=float(rlam), ok=bool(rok), min_den=float(rden.min()),
                          eta_min=float(S.eta[np.argmin(rden)]), dlam_ds=float(tl_new))
            crossings.append(info_c)
            np.save(f"ladder_{tag}_cross{len(crossings)}.npy", np.concatenate([[rlam], rphi]))
            print(f"*** CROSSING {len(crossings)}: lam={rlam:.13f} ok={rok} min_den={rden.min():.5f} "
                  f"(between λ={lam:.6f} p={p_prev:.6f} and λ={lm:.6f} p={p_new:.6f})", flush=True)
        phi, lam, tphi, tl, p_prev = ph, lm, tphi_new, tl_new, p_new
        branch.append((lam, p_new, float(den.min()), float(S.eta[np.argmin(den)]), float(tl)))
        if step % 10 == 0:
            print(f"step {step}: lam={lam:.10f} p={p_new:.10f} min_den={den.min():.3e} "
                  f"eta_min={S.eta[np.argmin(den)]:.3f} dlam/ds={tl:+.4f} ds={ds:.2e} t={time.time()-t0:.0f}s", flush=True)
        if it < 5:
            ds = min(ds * 1.4, dsmax)
        step += 1
        if den.min() < 1e-6:
            print("   approaching sonic limit; stop", flush=True)
            break
    np.save(f"ladder_{tag}_branch.npy", np.array(branch))
    json.dump(crossings, open(f"ladder_{tag}_crossings.json", "w"), indent=1)
    return branch, crossings


if __name__ == "__main__":
    N = int(sys.argv[1]); tag = sys.argv[2]
    nsteps = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
    L2 = float(sys.argv[4]) if len(sys.argv) > 4 else 120.0
    c = float(sys.argv[5]) if len(sys.argv) > 5 else 0.7
    direction = float(sys.argv[6]) if len(sys.argv) > 6 else -1.0
    run(N=N, tag=tag, nsteps=nsteps, L2=L2, c=c, direction=direction)
