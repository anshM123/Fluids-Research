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
    import os
    ss = float(os.environ.get('WKB_SS', -20.0))
    B = IPM(lam, Nb=Nb, hs=hs, s_sw=12.0, s_start=ss); Y = np.load(f); X = full(B, Y)
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
    # continuous cut-off: linear interpolation of D̂ = Dcut between the bracketing grid points (removes the
    # grid-step jitter of the phase, ≈ κ h_s/D₀); fixed number of quadrature points
    s_cut = s[kc - 1] + (Dcut - Dh_all[kc - 1]) / (Dh_all[kc] - Dh_all[kc - 1]) * (s[kc] - s[kc - 1])
    grid = np.linspace(s_lo, s_cut, int(round((s_cut - s_lo) / ds)) + 1)
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
            mz = re.search(r'_z([0-9.]+?)\.npy', f)
            if mz:
                lam = 1 / float(mz.group(1))
        o = wkb_phase(f, lam=lam, Dcut=Dcut, Nb=Nb_env, hs=hs_env)
        print(f"{f[:58]:58s} λ={o['lam']:.7f} m-2={o['m']-2:+.2e} D0={o['D0']:.5f} dip s={o['s_dip']:.3f} D̂={o['Dh_dip']:.4f} "
              f"s_cut={o['s_cut']:.3f} Φ={o['Phi'].real:.4f}{o['Phi'].imag:+.4f}i I={o['I'].real:.5f}{o['I'].imag:+.5f}i fail={o['nfail']}", flush=True)



def wkb_phase3(f, lam=None, Nb=None, hs=None, s_lo=-10.0, ds=0.025, Dcut=2.0, return_kappa=False, sub=4):
    """Phase with the original tracking up to the dip and continuation in K = κD̂ beyond it (front side), where D̂
    rises steeply and the original κ-tracking can fail on deep profiles.  Beyond the dip the step is ds/sub; the
    predictor is K_prev/D̂ (exact IPM scaling), the acceptance is continuity in K (25 %), and the HL root and the
    previous κ are fallbacks."""
    from bq_local_eig import hl_root
    o = wkb_phase(f, lam=lam, Nb=Nb, hs=hs, s_lo=s_lo, ds=ds, Dcut=Dcut, return_kappa=True)
    W = wall_data(f, lam, Nb, hs); s, D0 = W['s'], W['D0']
    g, k = o['grid'], o['kap'].copy()
    idip = np.argmin(np.abs(g - o['s_dip']))
    par = lambda sv: (np.interp(sv, s, W['D']) / D0, -np.interp(sv, s, W['Omb']), np.interp(sv, s, W['mu']),
                      np.interp(sv, s, W['G']), np.interp(sv, s, W['Ry']))
    gf = np.linspace(g[idip], o['s_cut'], max(2, int(round((o['s_cut'] - g[idip]) / (ds / sub))) + 1))
    kf = np.empty(len(gf), dtype=complex); kf[0] = k[idip]
    Kp = k[idip] * par(gf[0])[0]; prev = k[idip]; nfail = 0
    for j in range(1, len(gf)):
        p = par(gf[j]); Dh = p[0]; got = None
        for guess in (Kp / Dh, prev, hl_root(Dh, p[1])):
            try:
                r, ok = local_root(*p, guess)
            except Exception:
                ok = False
            if ok and np.isfinite(r) and abs(r * Dh - Kp) <= 0.25 * abs(Kp) + 1e-9:
                got = r; break
        if got is None:
            nfail += 1; got = Kp / Dh                 # hold K (not κ) across a failed point
        kf[j] = got; prev = got; Kp = got * Dh
    I = np.trapezoid(k[:idip + 1], g[:idip + 1]) + np.trapezoid(kf, gf)
    out = dict(lam=o['lam'], m=o['m'], D0=D0, s_cut=o['s_cut'], s_dip=o['s_dip'], Dh_dip=o['Dh_dip'],
               I=I, Phi=I / D0, nfail_front=nfail, nfail=o['nfail'])
    if return_kappa:
        out.update(grid=np.concatenate([g[:idip + 1], gf[1:]]), kap=np.concatenate([k[:idip + 1], kf[1:]]))
    return out
