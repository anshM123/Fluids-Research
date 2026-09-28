"""WKB phase of rapidly oscillating perturbations in the quasi-stagnant region of the Hou–Luo profile.
Local dispersion relation (leading order, Mellin symbol M ≈ −1/k for |k| → ∞):
    i D k² + k − Ω = 0... in the form  i D k² + k − Ω/D·D = 0  →  k(η) = (−1 + sqrt(1 + 4iΩ(η)))/(2 i D(η)),
with D = 1+λ+U/ξ the characteristic speed and Ω the vorticity (both along the profile, ξ = e^η).
Predicted: m − 2 ∝ Re[K exp(i∫k dη)]; with 1/ε ≈ 2z the crossings are spaced by Δz = π/(2 Re Φ₀) (or π/(4 Re Φ₀)
for a round trip), Φ₀ = ε∫k dη."""
import numpy as np, sys, glob, re
from hl_solver import HL
def phase(lam, q, N=8192, eta_cut=None):
    S = HL(lam, N=N)
    r = S.march(S.full(q))
    D, Om, eps = r['D'], r['Om'], r['eps']
    k = (-1 + np.sqrt(1 + 4j * Om)) / (2j * D)
    sel = slice(S.i0, None) if eta_cut is None else (S.eta >= S.eta[S.i0]) & (S.eta <= eta_cut)
    Phi = np.trapezoid(k[sel], S.eta[sel])
    return eps, Phi, r
for f in sys.argv[1:]:
    lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
    q = np.load(f)
    eps, Phi, r = phase(lam, q)
    # also a cut at the front: first η where D > 0.5
    S = HL(lam)
    etac = S.eta[np.argmax((r['D'] > 0.5) & (S.eta > -5))]
    _, Phic, _ = phase(lam, q, eta_cut=etac)
    print(f"λ={lam:.6f} z={1/(lam-1):7.3f} ε={eps:.5f} m={r['m']:.8f}  εΦ(full)={eps*Phi.real:+.5f}{eps*Phi.imag:+.5f}i   "
          f"εΦ(cut η<{etac:.2f})={eps*Phic.real:+.5f}{eps*Phic.imag:+.5f}i   Δz_pred π/(2ReΦ₀)={np.pi/(2*eps*Phic.real*r['m']/2):.4f}")
