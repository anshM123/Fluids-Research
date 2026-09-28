"""Test: Hilbert transform of even functions in log variable eta=ln(xi) as a Fourier multiplier.
H f(x) = (1/pi) p.v. int f(y)/(x-y) dy.  For even f:  H f(e^eta) = (1/pi) p.v. int f(e^mu)/sinh(eta-mu) dmu.
With f(e^eta)=e^{c eta} Psi(eta):  H f = e^{c eta} (K_c * Psi),  K_c^(k) = -i tanh(pi (k - i c)/2)."""
import numpy as np
def Hlog(Psi, eta, c):
    N = len(eta); P = N*(eta[1]-eta[0])
    k = 2*np.pi*np.fft.fftfreq(N, d=eta[1]-eta[0])
    mult = -1j*np.tanh(np.pi*(k - 1j*c)/2)
    return np.real(np.fft.ifft(mult*np.fft.fft(Psi)))
for c in [0.3, 0.5, 0.7]:
    L1, L2, N = 40.0, 60.0, 4096
    eta = -L1 + (L1+L2)*np.arange(N)/N
    xi = np.exp(eta)
    f = xi**2/(1+xi**2)                     # even; H f = -x/(1+x^2)
    Psi = np.exp(-c*eta)*f
    Hf = np.exp(c*eta)*Hlog(Psi, eta, c)
    exact = -xi/(1+xi**2)
    m = (eta > -20) & (eta < 20)
    print(f"c={c}: max err on |eta|<20: {np.max(np.abs(Hf-exact)[m]):.2e}")
# growing test: f = |x|^b (0<b<1): H f = -tan(pi b/2) sgn(x)|x|^b ; use f=(x^2)^(b/2) * x^2/(1+x^2) ... check pure power via analytic part only
