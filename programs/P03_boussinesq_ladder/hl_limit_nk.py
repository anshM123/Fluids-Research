"""Newton–Krylov continuation of the Hou–Luo λ → 1 limit family (hl_limit.py) in the front vorticity Ω_f, for the
large-Ω_f end where the damped fixed-point iteration fails. Unknowns: the layer coefficients b and u = log Ω_O on the
outer grid; residual = (fixed-point map) − (unknowns). Reports along the family: α = Σb/√2 (inverse-square-root
coefficient of the layer vorticity at the front), min over the outer region of D/√(ξ−1) and its location, and the
phase constant C = 2 Re a₀ with the implied spacing π/C.
usage: python3 hl_limit_nk.py START_NPZ OMF_START OMF_END STEP"""
import numpy as np, sys, time, warnings
from scipy.optimize import newton_krylov
import hl_limit as L
warnings.filterwarnings('ignore')
import os
TAG = os.environ.get('HL_TAG', '')
y, dy, t = L.y, L.dy, L.y - 1
J = L.J


def fp(b, Om_o, Om_f):
    g = -2.0 - L.K_LAYER @ Om_o
    b_new = -(2 / np.pi) * (L.Uk * (L.w_k * g)[None, :]).sum(axis=1)
    H = L.layer_H_outside(b_new, y) + L.K_OUT @ Om_o
    U = -2.0 + np.cumsum(H * dy)
    D = 2.0 + U / y
    q = dy / (y * np.maximum(D, 1e-10))
    u_new = np.log(Om_f) - np.cumsum(q)
    return b_new, u_new, D


def solve(Om_f, x0):
    def F(x):
        b, u = x[:J], x[J:]
        b_new, u_new, _ = fp(b, np.exp(u), Om_f)
        return np.concatenate([b_new - b, u_new - u])
    x = newton_krylov(F, x0, method='lgmres', f_tol=1e-9, maxiter=60, inner_maxiter=60, verbose=False)
    b, u = x[:J], x[J:]
    _, _, D = fp(b, np.exp(u), Om_f)
    return x, D


if __name__ == "__main__":
    d = np.load(sys.argv[1])
    Om_f, Om_end, step = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    b0 = np.zeros(J); nb = min(J, len(d['b'])); b0[:nb] = d['b'][:nb]
    u0 = np.log(np.maximum(d['Om_o'], 1e-300))
    if len(d['y']) != len(y) or not np.allclose(d['y'], y):       # start state from another grid: interpolate log Ω_O
        u0 = np.interp(np.log(y - 1), np.log(d['y'] - 1), u0)
    x = np.concatenate([b0, u0])
    t0 = time.time()
    while Om_f <= Om_end + 1e-12:
        try:
            x, D = solve(Om_f, x)
        except Exception as e:
            print(f"Ω_f = {Om_f:.3f}: Newton–Krylov failed ({type(e).__name__}: {str(e)[:80]})", flush=True)
            break
        b = x[:J]
        sel = t < 0.5
        r = D[sel] / np.sqrt(t[sel]); i = int(np.argmin(r))
        a0 = L.phase(b); C = 2 * a0.real
        print(f"Ω_f = {Om_f:7.3f}: α = {b.sum() / np.sqrt(2):.5f}; D/√t: front {r[0]:.4f}, min {r[i]:.5f} at t = {t[sel][i]:.2e}; "
              f"min D = {D.min():.2e}; C = {C:.5f}, π/C = {np.pi / C:.5f} ({time.time() - t0:.0f}s)", flush=True)
        np.savez(f"hl_limit{TAG}_Omf{Om_f:.3f}.npz", b=b, y=y, Om_o=np.exp(x[J:]), D=D)
        Om_f += step
