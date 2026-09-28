"""Fixed-lambda Picard iteration for CCF profile (integral form, avoids index-1 discretization issue):
   d(ln Psi)/d eta = b(eta) = lam/(1+lam+G) - c,   G = HTheta/xi = e^{(c-1)eta} (K_c * Psi)."""
import numpy as np, sys
from scipy.integrate import cumulative_trapezoid
def setup(L1, L2, N, c):
    h = (L1+L2)/N; eta = -L1 + h*np.arange(N)
    k = 2*np.pi*np.fft.fftfreq(N, d=h)
    mK = -1j*np.tanh(np.pi*(k-1j*c)/2)
    return eta, h, mK
def G_of(Psi, eta, mK, c):
    return np.exp((c-1)*eta)*np.real(np.fft.ifft(mK*np.fft.fft(Psi)))
def picard(lam, L1=30, L2=150, N=4096, c=0.75, p0=2.0, iters=400, relax=0.3, verbose=False):
    eta, h, mK = setup(L1, L2, N, c)
    beta = lam/(1+lam)
    i0 = np.argmin(np.abs(eta))
    # initial guess Theta ~ xi^p0 / (1+xi^2)^((p0-beta)/2), amplitude so that h1 = lam/p0-1-lam
    Theta = np.exp(p0*eta - 0.5*(p0-beta)*np.logaddexp(0, 2*eta))
    Psi = np.exp(-c*eta)*Theta
    h1 = -(2/np.pi)*h*np.sum(Psi*np.exp((c-1)*eta))
    Psi *= (lam/p0-1-lam)/h1
    lnP = np.log(Psi)
    hist = []
    for it in range(iters):
        Psi = np.exp(lnP)
        G = G_of(Psi, eta, mK, c)
        den = 1+lam+G
        if np.min(den) <= 0:
            return None, dict(fail="sonic", it=it, mind=np.min(den))
        b = lam/den - c
        phi = cumulative_trapezoid(b, eta, initial=0.0)
        phi = phi - phi[i0] + lnP[i0]          # keep value at eta=0 (normalisation by scaling symmetry)
        diff = np.max(np.abs(phi-lnP)[(eta>-20)&(eta<60)])
        lnP = (1-relax)*lnP + relax*phi
        h1 = -(2/np.pi)*h*np.sum(np.exp(lnP)*np.exp((c-1)*eta))
        hist.append((diff, h1))
        if verbose and it % 20 == 0:
            print(f" it {it}: dlnPsi={diff:.2e} h1={h1:.10f} p={lam/(1+lam+h1):.8f}")
        if diff < 1e-12:
            break
    return lnP, dict(h1=h1, p=lam/(1+lam+h1), it=it, diff=diff, eta=eta)
if __name__ == "__main__":
    for lam in [float(x) for x in sys.argv[1].split(',')]:
        lnP, info = picard(lam, verbose=False)
        print(f"lam={lam}: ", {k: v for k, v in info.items() if k != 'eta'}, flush=True)
