"""Continuation of the Hou–Luo λ→1 limit family (hl_limit.py) to large front vorticity Ω_f, with stronger damping,
reporting the phase constant C(Ω_f) = 2 Re a₀ and the implied spacing π/C, plus the outer front coefficient
k = D/√(ξ − 1) at ξ − 1 = 1e-4 and the layer data used to compare with finite-ε profiles."""
import numpy as np, sys, warnings; warnings.filterwarnings('ignore')
import hl_limit as L
targets = [float(v) for v in sys.argv[1:]]
d = np.load("hl_limit_Omf4.000.npz")
prev = dict(Om_f=4.0, b=d['b'], Om_o=d['Om_o'])
Om_f = 4.0
for tgt in targets:
    while Om_f < tgt - 1e-12:
        step = 0.25 if Om_f < 10 else (1.0 if Om_f < 40 else 5.0)
        Om_f = min(tgt, Om_f + step)
        r = L.solve(Om_f, init=prev, iters=3000, relax=0.08, tol=1e-9)
        if not r['ok']:
            print(f"   Ω_f = {Om_f}: not converged (change {r['err']:.1e}, min D {r['Dmin']:.4f})", flush=True)
        prev = r
    b = r['b']
    Theta_c = float(np.sum(b * np.sin(L.jj * np.pi / 2) / L.jj))
    th0 = np.array([np.pi / 2 - 1e-4]); _, Om0, _ = L.layer_fields(b, th0)
    c0 = float(Om0[0] / np.cos(th0[0]) / 2)
    a0 = L.phase(b); C = 2 * a0.real
    k = np.interp(1 + 1e-4, L.y, r['D']) / 1e-2
    t = np.array([0.5, 0.2, 0.1]); x, Oml, Thl = L.layer_fields(b, np.arccos(1 - t))
    print(f"Ω_f = {Om_f:7.2f}: ok={r['ok']} ({r['it']} it); c0 = {c0:.5f}, Θ(1) = {Theta_c:.5f}, k = {k:.4f}; "
          f"Θ(t=.5,.2,.1) = {np.round(Thl, 4).tolist()}, Ω = {np.round(Oml, 4).tolist()}; "
          f"a₀ = {a0.real:.5f}{a0.imag:+.5f}i, C = {C:.5f}, π/C = {np.pi / C:.5f}", flush=True)
    np.savez(f"hl_limit_Omf{Om_f:.3f}.npz", b=b, y=L.y, Om_o=r['Om_o'], D=r['D'], U=r['U'])
