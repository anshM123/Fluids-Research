import numpy as np, time
from bq_logpolar import BQLogPolar, _march
lam = 1.4
t0 = time.time()
B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
print(f"setup: Ns={B.Ns} Nb={B.Nb} a={B.a:.4f} ({time.time()-t0:.1f}s)")
# ---- Poisson test: Ψ = f(r) sin 2β, f = r²/(1+r²)^q
q = 0.3
s = B.s[:, None]; b = B.beta[None, :]
r = np.exp(s)
f = r**2 / (1 + r**2)**q
# f' and f'' analytically
fp = 2*r/(1+r**2)**q - 2*q*r**3/(1+r**2)**(q+1)
fpp = (2/(1+r**2)**q - 4*q*r**2/(1+r**2)**(q+1) - 6*q*r**2/(1+r**2)**(q+1) + 4*q*(q+1)*r**4/(1+r**2)**(q+2))
Om = -(fpp + fp/r - 4*f/r**2) * np.sin(2*b)
t0 = time.time()
P = B.poisson(Om)
Psi = np.exp(B.a*s) * P
sel = (B.s > -30) & (B.s < 30)
err = np.abs(Psi[sel] - (f*np.sin(2*b))[sel]) / np.abs(f[sel]).max(axis=None)
print(f"Poisson: max rel err over s∈(-30,30): {np.abs(Psi[sel] - (f*np.sin(2*b))[sel]).max() / (np.abs(f[sel])).max():.2e};"
      f" pointwise rel at s=0: {np.abs(Psi[np.argmin(abs(B.s))] - (f*np.sin(2*b))[np.argmin(abs(B.s))]).max():.2e}  ({time.time()-t0:.2f}s)")
# velocity of the manufactured Ψ: U_r/r = e^{-2s}Ψ_β = e^{-2s} f 2cos2β ; w = -e^{-2s}Ψ_s = -e^{-2s} r f' sin 2β
Ur, w, wt, Urh, wh, wth = B.velocity(P)
Ur_ex = np.exp(-2*s) * f * 2*np.cos(2*b); w_ex = -np.exp(-2*s) * r * fp * np.sin(2*b)
print(f"velocity: max err Ur {np.abs(Ur-Ur_ex)[sel].max():.2e} (scale {np.abs(Ur_ex[sel]).max():.2f}), w {np.abs(w-w_ex)[sel].max():.2e}")
# ---- March test with pure strain U = (−A y1, A y2): U_r/r = −A cos2β, w = A sin2β
A = 2.2
m = (lam - 1) / (1 + lam - A); C = m / (A - 1)
Ns, Nb1 = B.Ns, B.Nb + 1
Urs = -A*np.cos(2*B.beta)[None, :] * np.ones((Ns, 1)); ws = A*np.sin(2*B.beta)[None, :] * np.ones((Ns, 1))
wts = B._wtan(ws)
s0 = B.s[B.i0]
t0 = time.time()
Th, Omh = _march(B.s, B.hs, np.full(Nb1, np.exp(m*s0)), np.full(Nb1, C*np.exp((m-1)*s0)), Urs, ws, wts, Urs, ws, wts,
                 B.Db, B.cb, B.sb, lam, m, B.i0)
t1 = time.time()
Th, Omh = _march(B.s, B.hs, np.full(Nb1, np.exp(m*s0)), np.full(Nb1, C*np.exp((m-1)*s0)), Urs, ws, wts, Urs, ws, wts,
                 B.Db, B.cb, B.sb, lam, m, B.i0)
j = np.argmin(abs(B.s - 0.0))
print(f"march (pure strain, m={m:.4f}): Θ̂/e^(ms) at s=0 range [{(Th[j]/np.exp(m*B.s[j])).min():.10f}, {(Th[j]/np.exp(m*B.s[j])).max():.10f}], "
      f"Ω̂/(C e^((m-1)s)) [{(Omh[j]/(C*np.exp((m-1)*B.s[j]))).min():.10f}, {(Omh[j]/(C*np.exp((m-1)*B.s[j]))).max():.10f}]"
      f"  (first call {t1-t0:.1f}s incl. compile, second {time.time()-t1:.2f}s)")
