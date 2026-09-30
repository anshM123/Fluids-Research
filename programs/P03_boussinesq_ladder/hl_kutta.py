"""The Hou–Luo λ → 1 limit problem closed by the Kutta condition at the front.

hl_limit.py solves the limit problem (stalled layer |ξ| < 1 with H[Ω] = −2, outer transport Ω + ξDΩ_ξ = 0) as a
one-parameter family in the front vorticity Ω_f = Ω(1⁺). Along the family the inverse-square-root coefficient of
the layer vorticity, α = Σb/√2 (Ω ≈ α(1−ξ)^{−1/2}), decreases to 0 as Ω_f grows, and the phase constant C(Ω_f)
levels off. The limit member is characterised directly by α = 0 (bounded-edge, Kutta-type condition): the outer
vorticity is then unbounded at the front, so its amplitude is fixed at an interior reference point y_* and
determined by α = 0, instead of by Ω_f.
    b = b⁰ + L[Ω_O]  (layer inversion, affine in the outer vorticity),   Ω_O = A·S,   A = −Σb⁰ / Σ L[S],
    S(y) = exp(−∫_{y_*}^{y} dy'/(y' D)).
usage: python3 hl_kutta.py [START_NPZ]"""
import numpy as np, sys
import hl_limit as L

y, dy = L.y, L.dy
YSTAR = 2.0
istar = int(np.argmin(np.abs(y - YSTAR)))
b0 = -(2 / np.pi) * (L.Uk * (L.w_k * (-2.0))[None, :]).sum(axis=1)          # Ω_O = 0 layer solution (b0 = [2, 0, …])


def layer_b(Om_o):
    return b0 + (2 / np.pi) * (L.Uk * (L.w_k * (L.K_LAYER @ Om_o))[None, :]).sum(axis=1)


def outer_from(b, Om_o):
    H = L.layer_H_outside(b, y) + L.K_OUT @ Om_o
    U = -2.0 + np.cumsum(H * dy)
    D = 2.0 + U / y
    return D, U, H


def shape(D):
    q = dy / (y * np.maximum(D, 1e-12))
    cum = np.cumsum(q)
    return np.exp(-(cum - cum[istar]))


def solve_kutta(S0, iters=4000, relax=0.05, tol=1e-10, verbose=True):
    S = S0 / S0[istar]
    for it in range(iters):
        LS = layer_b(S) - b0
        A = -b0.sum() / LS.sum()
        Om_o = A * S
        b = b0 + A * LS
        D, U, H = outer_from(b, Om_o)
        Snew = shape(D)
        err = np.abs(np.log(Snew[: istar * 3]) - np.log(S[: istar * 3])).max()
        S = np.exp((1 - relax) * np.log(S) + relax * np.log(Snew))
        if verbose and it % 200 == 0:
            print(f"  it {it}: change {err:.2e}  A = {A:.6f}  Σb = {b.sum():.2e}  min D = {D.min():.3e}  "
                  f"Ω_O(t=1e-6, 1e-3, 1e-1) = {np.interp([1e-6, 1e-3, 1e-1], y - 1, Om_o).round(3).tolist()}", flush=True)
        if err < tol:
            break
    return dict(b=b, Om_o=Om_o, D=D, U=U, H=H, A=A, err=err, it=it)


if __name__ == "__main__":
    start = sys.argv[1] if len(sys.argv) > 1 else "hl_limit_Omf10.000.npz"
    d = np.load(start)
    r = solve_kutta(d['Om_o'])
    b = r['b']
    a0 = L.phase(b); C = 2 * a0.real
    Theta_c = float(np.sum(b * np.sin(L.jj * np.pi / 2) / L.jj))
    th0 = np.array([np.pi / 2 - 1e-4]); _, Om0, _ = L.layer_fields(b, th0)
    c0 = float(Om0[0] / np.cos(th0[0]) / 2)
    print(f"Kutta limit: change {r['err']:.1e} after {r['it']} it; Σb = {b.sum():.2e}; c0 = {c0:.5f}, Θ(1) = {Theta_c:.5f}; "
          f"a₀ = {a0.real:.5f}{a0.imag:+.5f}i, C = {C:.5f}, π/C = {np.pi / C:.5f}", flush=True)
    np.savez("hl_limit_kutta.npz", b=b, y=y, Om_o=r['Om_o'], D=r['D'], U=r['U'])
