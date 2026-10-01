"""Dip geometry and the anatomy of the IPM phase on a resolved profile (task: divergence class of Φ₀).
Uses the exact scaling of the IPM local problem, κ = K(G, μ, R_y)/D̂, so that Φ₀ = ∫K ds / D.
Reports: D̂_min, the dip position s_d, the curvature D̂'' at the dip (local quadratic fit), the front position
(D̂ = 2 beyond the dip), K at the dip, and the split of Re Φ₀ into layer (s < s_d − w) and dip (the rest, up to the
cut-off) contributions, with the local-quadratic estimate π K_d / √(D̂_min D̂''/2) / D₀ of the dip part.
usage: ipm_dip_geometry.py FILE [FILE …]   (env WKB_NB, WKB_HS, WKB_SS for files without grid tags)"""
import numpy as np, sys, re, os
from ipm_wkb import wkb_phase, wall_data

def analyse(f):
    lam = None
    if not f.startswith("ipm_rung_"):
        mm = re.search(r'lam([0-9.]+?)(?:_hs[0-9.]+)?\.npy', f); lam = float(mm.group(1)) if mm else None
        mz = re.search(r'_z([0-9.]+?)\.npy', f)
        if mz: lam = 1 / float(mz.group(1))
    Nb = int(os.environ.get('WKB_NB', 0)) or None; hs = float(os.environ.get('WKB_HS', 0)) or None
    o = wkb_phase(f, lam=lam, Nb=Nb, hs=hs, Dcut=2.0, return_kappa=True)
    W = wall_data(f, lam=lam, Nb=Nb, hs=hs)
    s, Dh = W['s'], W['D'] / W['D0']
    sel = (s > -3) & (s < 2); i = np.where(sel)[0][np.argmin(Dh[sel])]
    j = np.arange(i - 3, i + 4); c = np.polyfit(s[j] - s[i], Dh[j], 2)        # local quadratic at the dip
    sd = s[i] - c[1] / (2 * c[0]); Dmin = np.polyval(c, sd - s[i]); D2 = 2 * c[0]
    g, kap = o['grid'], o['kap']
    Kd = np.interp(sd, g, (kap * np.interp(g, s, Dh)).real) + 1j * np.interp(sd, g, (kap * np.interp(g, s, Dh)).imag)
    w = np.sqrt(2 * Dmin / D2)                                                   # half-width where D̂ = 2 D̂_min
    lay = g < sd - 3 * w
    P_lay = np.trapezoid(kap[lay], g[lay]) / o['D0']; P_dip = o['Phi'] - P_lay
    est = np.pi * Kd / np.sqrt(Dmin * D2 / 2) / o['D0']
    return dict(lam=o['lam'], z=1 / o['lam'], D0=o['D0'], Dmin=Dmin, D2=D2, sd=sd, w=w, s_cut=o['s_cut'], Kd=Kd,
                Phi=o['Phi'], P_lay=P_lay, P_dip=P_dip, est=est)

if __name__ == "__main__":
    print(f"{'file':40s} {'z':>7s} {'D̂min':>7s} {'D̂″':>7s} {'w':>6s} {'s_d':>6s} {'s_cut':>6s} {'K_d':>16s} "
          f"{'ReΦ':>8s} {'layer':>8s} {'dip':>8s} {'dip est':>8s}")
    for f in sys.argv[1:]:
        r = analyse(f)
        print(f"{f[-40:]:40s} {r['z']:7.4f} {r['Dmin']:7.4f} {r['D2']:7.2f} {r['w']:6.3f} {r['sd']:6.3f} {r['s_cut']:6.3f} "
              f"{r['Kd'].real:7.3f}{r['Kd'].imag:+7.3f}i {r['Phi'].real:8.3f} {r['P_lay'].real:8.3f} {r['P_dip'].real:8.3f} "
              f"{r['est'].real:8.3f}", flush=True)
