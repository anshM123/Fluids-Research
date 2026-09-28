"""Independent check of a computed profile: evaluate the self-similar Boussinesq equations in polar form with
finite differences in s (6th order) and Chebyshev in β, on the fields reconstructed from the solver state.
    R_θ = (1−λ)Θ + V·∇Θ,   R_ω = Ω + V·∇Ω − ∂₁Θ,   R_ψ = ΔΨ + Ω,    V = (1+λ)y + ∇^⊥Ψ."""
import numpy as np, sys, re
from bq_solver import BQ
from bq_newton import full

def d_s(F, h):
    """6th-order centred first derivative along axis 0 (interior only)"""
    c = np.array([-1, 9, -45, 0, 45, -9, 1]) / 60.0
    out = np.full_like(F, np.nan)
    out[3:-3] = sum(c[k] * F[k:len(F) - 6 + k] for k in range(7)) / h
    return out

for f in sys.argv[1:]:
    lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    mNb = re.search(r'Nb(\d+)_hs([0-9.]+?)_', f)
    Nb, hs = (int(mNb.group(1)), float(mNb.group(2))) if mNb else (32, 0.025)
    B = BQ(lam, Nb=Nb, hs=hs); Y = np.load(f); X = full(B, Y)
    print(f"file {f}: Nb={Nb} hs={hs}")
    r = B.march(X / B.ea2[:, None], return_all=True)
    m = r['m']; s = B.s; h = B.hs
    cpow = np.where(B.cb > 0, B.cb, 0.0)
    Th = r['Th'] * cpow[None, :] ** m; Om = r['Omega']
    R = np.exp(s)[:, None]
    Psi = R ** 2 * X
    Db = B.Db
    # derivatives
    Th_s = d_s(Th, h); Om_s = d_s(Om, h); Psi_s = d_s(Psi, h); Psi_ss = d_s(Psi_s, h)
    Th_b = Th @ Db.T; Om_b = Om @ Db.T; Psi_b = Psi @ Db.T; Psi_bb = Psi @ (Db @ Db).T
    Th_r = Th_s / R; Om_r = Om_s / R; Psi_r = Psi_s / R
    Psi_rr = (Psi_ss - Psi_s) / R ** 2
    Ur = Psi_b / R; Ub = -Psi_r
    Vr = (1 + lam) * R + Ur; Vb = Ub
    cb, sb = B.cb[None, :], B.sb[None, :]
    VgTh = Vr * Th_r + Vb * Th_b / R
    VgOm = Vr * Om_r + Vb * Om_b / R
    d1Th = cb * Th_r - sb * Th_b / R
    Rth = (1 - lam) * Th + VgTh
    Rom = Om + VgOm - d1Th
    Rps = Psi_rr + Psi_r / R + Psi_bb / R ** 2 + Om
    sel = (s > -6) & (s < 6)
    def rel(Rs, *terms):
        scale = sum(np.abs(t[sel]) for t in terms) + 1e-300
        return np.nanmax(np.abs(Rs[sel]) / scale.max()), np.nanmax(np.abs(Rs[sel]))
    print(f"λ={lam:.7f} m={m:.10f}: max|R_θ|/max|terms| = {rel(Rth, (1-lam)*Th, VgTh)[0]:.2e}, "
          f"|R_ω| = {rel(Rom, Om, VgOm, d1Th)[0]:.2e}, |R_ψ| = {rel(Rps, Psi_rr, Om)[0]:.2e}   (r ∈ [e^-6, e^6])")
    for name, Rs in (("R_θ", Rth), ("R_ω", Rom), ("R_ψ", Rps)):
        A = np.abs(Rs[sel]); i, j = np.unravel_index(np.nanargmax(A), A.shape)
        print(f"   {name}: max {np.nanmax(A):.2e} at s={s[sel][i]:.3f}, β={B.beta[j]:.4f} (j={j});  max over j<Nb-1 excluding axis: "
              f"{np.nanmax(A[:, :-1]):.2e}; excluding 2 near axis: {np.nanmax(A[:, :-2]):.2e}; interior β only: {np.nanmax(A[:, 1:-1]):.2e}")
    for lo, hi in ((-8, -2), (-2, 1), (1, 3), (3, 6), (6, 20)):
        ss = (s > lo) & (s < hi)
        def prel(Rs, *terms):
            sc = sum(np.abs(t[ss]) for t in terms).max(axis=1)        # row scale: max over β at each s
            return np.nanmax(np.abs(Rs[ss]).max(axis=1) / sc)
        print(f"   s∈[{lo},{hi}]: rel. residual (row-normalised): θ {prel(Rth, (1-lam)*Th, VgTh):.2e}  ω {prel(Rom, Om, VgOm, d1Th):.2e}  ψ {prel(Rps, Psi_rr, Om):.2e}")
    ss = (s > 3) & (s < 6)
    sc = np.abs((1-lam)*Th[ss]) + np.abs(VgTh[ss])
    Q = np.abs(Rth[ss]) / (sc + 1e-12 * sc.max())
    print("   R_θ relative, s∈[3,6], max over s for each β index:", np.array2string(np.nanmax(Q, axis=0), precision=1, max_line_width=250))
    sc = np.abs(Om[ss]) + np.abs(VgOm[ss]) + np.abs(d1Th[ss])
    Q = np.abs(Rom[ss]) / (sc + 1e-12 * sc.max())
    print("   R_ω relative, s∈[3,6], per β:", np.array2string(np.nanmax(Q, axis=0), precision=1, max_line_width=250))
