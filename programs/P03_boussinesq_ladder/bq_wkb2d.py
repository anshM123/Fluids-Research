"""2D WKB phase Φ = (1/ε)∫κ(s) ds with the local 2D dispersion relation (bq_local_eig) evaluated on the boundary data
of a computed Boussinesq profile. Prediction: crossings spaced by Δz = π/(m Re εΦ); decay per half period exp(π Im/Re)."""
import numpy as np, sys, re
from bq_solver import BQ
from bq_newton import full
from bq_local_eig import local_root, hl_root


def boundary_params(f):
    lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    B = BQ(lam); Y = np.load(f); X = full(B, Y)
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


if __name__ == "__main__":
    s_lo, s_hi, ds = -10.0, 4.0, 0.05
    for f in sys.argv[1:]:
        lam, m, eps, s, D, mu, G, Thy, Omb = boundary_params(f)
        grid = np.arange(s_lo, s_hi + 1e-9, ds)
        kap, kap_hl = [], []
        guess = None
        for sv in grid:
            Dh = np.interp(sv, s, D) / eps; ch = -np.interp(sv, s, Omb)
            p = (Dh, ch, np.interp(sv, s, mu), np.interp(sv, s, G), np.interp(sv, s, Thy))
            g = hl_root(Dh, ch) if guess is None else guess
            k, ok = local_root(*p, g)
            if not ok or not np.isfinite(k):
                k = np.nan
            else:
                guess = k
            kap.append(k); kap_hl.append(hl_root(Dh, ch))
        kap = np.array(kap); good = np.isfinite(kap)
        Phi = np.trapezoid(kap[good], grid[good]); Phih = np.trapezoid(np.array(kap_hl), grid)
        print(f"λ={lam:.5f} z={1/(lam-1):.3f} ε={eps:.4f} m={m:.5f}: εΦ_2D = {Phi.real:.5f}{Phi.imag:+.5f}i "
              f"(HL formula {Phih.real:.5f}{Phih.imag:+.5f}i);  Δz_pred = {np.pi/(m*Phi.real):.4f}, "
              f"half-period decay ratio = {np.exp(np.pi*abs(Phi.imag)/Phi.real):.3f}; failed points: {np.sum(~good)}", flush=True)
        np.save(f"wkb2d_lam{lam:.4f}.npy", np.array([grid, kap]))
