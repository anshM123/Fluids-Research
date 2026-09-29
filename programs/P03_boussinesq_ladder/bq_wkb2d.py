"""2D WKB phase Φ = (1/ε)∫κ(s) ds with the local 2D dispersion relation (bq_local_eig) evaluated on the boundary data
of a computed Boussinesq profile. Prediction: crossings spaced by Δz = π/(m Re εΦ); decay per half period exp(π Im/Re)."""
import numpy as np, sys, re
from bq_solver import BQ
from bq_newton import full
from bq_local_eig import local_root, hl_root


def boundary_params(f):
    lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    mNb = re.search(r'Nb(\d+)_hs([0-9.]+?)_', f)
    Nb, hs = (int(mNb.group(1)), float(mNb.group(2))) if mNb else (32, 0.025)
    B = BQ(lam, Nb=Nb, hs=hs); Y = np.load(f); X = full(B, Y)
    r = B.march(X / B.ea2[:, None], return_all=True)
    m, A = r['m'], r['A']; eps = 1 + lam - A
    s = B.s
    D = (1 + lam) + r['Ur'][:, 0]
    dD = np.gradient(D, s)
    cpow = np.where(B.cb > 0, B.cb, 0.0)
    Th = r['Th'] * cpow[None, :] ** m
    Thb = Th[:, 0]
    G = np.exp(-s) * np.gradient(Thb, s)                    # ∂_xΘ on the boundary
    Thy = np.exp(-s) * (r['Th'] @ B.Db.T)[:, 0]              # ∂_yΘ = (1/x) ∂_βΘ̂ at β = 0
    Omb = r['Omega'][:, 0]
    mu = 2 * (1 + lam) - D - dD
    return lam, m, eps, s, D, mu, G, Thy, Omb


def wkb_phase(f, s_lo=-10.0, ds=0.025, Dcut=3.0, verbose=False, return_kappa=False):
    """ε∫κ ds over the quasi-stagnant region, from s_lo up to the front (first s beyond the dip with D/ε > Dcut).
    The root is tracked from the Hou–Luo closed form at every point (robust against branch jumping)."""
    lam, m, eps, s, D, mu, G, Thy, Omb = boundary_params(f)
    Dh_all = D / eps
    sel = (s > -3) & (s < 2)
    kd = np.argmax(sel) + np.argmin(Dh_all[sel])                     # dip
    kc = kd + np.argmax(Dh_all[kd:] > Dcut)                          # front
    s_cut = s[kc]
    grid = np.arange(s_lo, s_cut + 1e-12, ds)
    kap = np.full(len(grid), np.nan, dtype=complex)
    prev = None
    for i, sv in enumerate(grid):
        Dh = np.interp(sv, s, D) / eps; ch = -np.interp(sv, s, Omb)
        p = (Dh, ch, np.interp(sv, s, mu), np.interp(sv, s, G), np.interp(sv, s, Thy))
        hl = hl_root(Dh, ch)
        cands = [hl] if prev is None else [prev, hl]
        for g in cands:
            try:
                k, ok = local_root(*p, g)
            except Exception:
                ok = False
            if ok and np.isfinite(k) and abs(k - hl) < 0.6 * abs(hl) + 1e-6:
                kap[i] = k; prev = k
                break
    good = np.isfinite(kap)
    # fill isolated failures by interpolation
    if (~good).any() and good.sum() > 2:
        kap[~good] = np.interp(grid[~good], grid[good], kap[good].real) + 1j * np.interp(grid[~good], grid[good], kap[good].imag)
    if return_kappa:
        return lam, m, eps, s_cut, np.trapezoid(kap, grid), (~good).sum(), grid, kap
    return lam, m, eps, s_cut, np.trapezoid(kap, grid), (~good).sum()


if __name__ == "__main__":
    for f in sys.argv[1:]:
        lam, m, eps, s_cut, Phi, nf = wkb_phase(f)
        print(f"λ={lam:.8f} z={1/(lam-1):.4f} ε={eps:.5f} m={m:.6f}: s_front={s_cut:.3f}  εΦ_cut = {Phi.real:.5f}{Phi.imag:+.5f}i  "
              f"Φ = m z εΦ = {m/(lam-1)*Phi.real:.5f}{m/(lam-1)*Phi.imag:+.5f}i   (filled {nf})", flush=True)
