import numpy as np
from bq_logpolar import BQLogPolar
from bq_newton import residual
lam = 1.92
B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
Y = np.load("bq_X_lam1.9200.npy")
R, info = residual(B, Y)
ss = B.s[B.i0:, None]; b = B.beta[None, :]
for name, v in (("uniform strain (all s)", -0.5*np.sin(2*b)*np.ones_like(ss)),
                ("strain localized s<-10", -0.5*np.sin(2*b)*(ss < -10)),
                ("strain decaying like far field", -0.5*np.sin(2*b)/(1+np.exp(2*ss))**(1/(2*(1+lam))))):
    for eps in (1e-6, 1e-8):
        Rp, ip = residual(B, Y + eps*v)
        Jv = (Rp - R)/eps
        # project Jv onto the strain mode shape in the window s∈[-20,-10]
        sel = (ss[:, 0] < -10)
        proj = np.sum(Jv[sel]*np.sin(2*b)) / np.sum((np.sin(2*b)*np.ones((sel.sum(), 1)))**2)
        print(f"{name:32s} eps={eps:.0e}: dA/dε={(ip['A']-info['A'])/eps:+.5f}  |Jv|max={np.abs(Jv).max():.3e}  "
              f"strain-projection of Jv (s<-10) = {proj:+.5f}", flush=True)
