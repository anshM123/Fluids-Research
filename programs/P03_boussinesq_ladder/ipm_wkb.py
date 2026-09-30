"""IPM WKB phase Φ = (1/D₀)∫κ(s) ds along the wall of a computed IPM profile (D₀ = 1+λ−A = λ/m, the
stagnation-point value of D = V₁/x, so D̂ = D/D₀ = 1 at the origin), with the local dispersion relation of
ipm_local_eig.  Prediction: m(λ) − 2 ∝ Re exp(iΦ) along the branch, i.e. smooth profiles one per half-turn of ReΦ
and a decay factor exp(ΔImΦ) per half period."""
import numpy as np, sys, re
from ipm_solver import IPM
from bq_newton import full
from ipm_local_eig import local_root, LocalEigIPM


def wall_data(f, lam=None, Nb=None, hs=None):
    if lam is None:
        lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    if Nb is None or hs is None:
        mNb = re.search(r'Nb(\d+)_hs([0-9.]+?)_', f)
        Nb, hs = (int(mNb.group(1)), float(mNb.group(2))) if mNb else (32, 0.025)
    B = IPM(lam, Nb=Nb, hs=hs, s_sw=12.0); Y = np.load(f); X = full(B, Y)
    r = B.march(X / B.ea2[:, None], return_all=True)
    m, A = r['m'], r['A']; D0 = 1 + lam - A
    s = B.s
    D = (1 + lam) + r['Ur'][:, 0]
    dD = np.gradient(D, s)
    Rw = r['Th'][:, 0]                                      # R at the wall (cos β = 1)
    G = np.exp(-s) * np.gradient(Rw, s)                     # ∂_x R
    Ry = np.exp(-s) * (r['Th'] @ B.Db.T)[:, 0]              # ∂_y R = (1/x) ∂_β R̂ at β = 0
    Omb = r['Omega'][:, 0]
    mu = 2 * (1 + lam) - D - dD
    return dict(lam=lam, m=m, A=A, D0=D0, s=s, D=D, mu=mu, G=G, Ry=Ry, Omb=Omb)


def wkb_phase(f, lam=None, Nb=None, hs=None, s_lo=-10.0, ds=0.025, Dcut=3.0, return_kappa=False):
    """(1/D₀)∫κ ds from s_lo to the front (first s beyond the dip where D/D₀ > Dcut); the root is tracked from the
    Hou–Luo-type guess (−1+√(1+4iĉ))/(2iD̂) with continuation from the previous point."""
    from bq_local_eig import hl_root
    W = wall_data(f, lam, Nb, hs)
    s, D0 = W['s'], W['D0']
    Dh_all = W['D'] / D0
    sel = (s > -3) & (s < 2)
    kd = np.argmax(sel) + np.argmin(np.where(sel, Dh_all, np.inf)[np.argmax(sel):np.argmax(sel) + sel.sum()])
    kc = kd + np.argmax(Dh_all[kd:] > Dcut)
    s_cut = s[kc]
    grid = np.arange(s_lo, s_cut + 1e-12, ds)
    kap = np.full(len(grid), np.nan, dtype=complex)
    prev = None
    for i, sv in enumerate(grid):
        Dh = np.interp(sv, s, W['D']) / D0; ch = -np.interp(sv, s, W['Omb'])
        p = (Dh, ch, np.interp(sv, s, W['mu']), np.interp(sv, s, W['G']), np.interp(sv, s, W['Ry']))
        hl = hl_root(Dh, ch)
        for g in ([hl] if prev is None else [prev, hl]):
            try:
                k, ok = local_root(*p, g)
            except Exception:
                ok = False
            if ok and np.isfinite(k) and abs(k) > 1e-8 and abs(k - hl) < 0.8 * abs(hl) + 1e-6:
                kap[i] = k; prev = k
                break
    good = np.isfinite(kap)
    # far tail (κ < 0.005, shooting ill-conditioned): the small-ĉ asymptote κ ≈ κ_HL (1 % accurate at s = −3)
    tail = (~good) & (grid < -6.0)
    for i in np.where(tail)[0]:
        kap[i] = hl_root(np.interp(grid[i], s, W['D']) / D0, -np.interp(grid[i], s, W['Omb']))
    rest = ~np.isfinite(kap)
    if rest.any() and (~rest).sum() > 2:
        ok_ = ~rest
        kap[rest] = np.interp(grid[rest], grid[ok_], kap[ok_].real) + 1j * np.interp(grid[rest], grid[ok_], kap[ok_].imag)
    Phi = np.trapezoid(kap, grid) / D0
    out = dict(lam=W['lam'], m=W['m'], D0=D0, s_cut=s_cut, s_dip=s[kd], Dh_dip=Dh_all[kd], Phi=Phi,
               nfail=int((~good & (grid >= -6.0)).sum()), I=np.trapezoid(kap, grid))
    if return_kappa:
        out.update(grid=grid, kap=kap, good=good)
    return out


if __name__ == "__main__":
    import os
    Dcut = float(sys.argv[1])
    Nb_env = int(os.environ.get('WKB_NB', 0)) or None
    hs_env = float(os.environ.get('WKB_HS', 0)) or None
    for f in sys.argv[2:]:
        lam = None
        if not f.startswith("ipm_rung_"):
            mm = re.search(r'lam([0-9.]+?)(?:_hs[0-9.]+)?\.npy', f)
            lam = float(mm.group(1)) if mm else None
        o = wkb_phase(f, lam=lam, Dcut=Dcut, Nb=Nb_env, hs=hs_env)
        print(f"{f[:58]:58s} λ={o['lam']:.7f} m-2={o['m']-2:+.2e} D0={o['D0']:.5f} dip s={o['s_dip']:.3f} D̂={o['Dh_dip']:.4f} "
              f"s_cut={o['s_cut']:.3f} Φ={o['Phi'].real:.4f}{o['Phi'].imag:+.4f}i I={o['I'].real:.5f}{o['I'].imag:+.5f}i fail={o['nfail']}", flush=True)
