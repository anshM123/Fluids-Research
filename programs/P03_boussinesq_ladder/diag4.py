import numpy as np
from bq_solver import BQ
from bq_newton import residual
B = BQ(1.92)
Y = np.load('Y2_lam1.9200_Nb32_hs0.025.npy')
R, info = residual(B, Y)
s = B.s[B.i0:]
mode = -0.5*np.sin(2*B.beta)
v = np.tile(mode, (len(s), 1)) * (0.5*(1 - np.tanh((s + 8)/1.0)))[:, None]   # localised strain, max 0.5
out = {}
for eta in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11):
    Rp, ip = residual(B, Y + eta*v)
    d = (Rp - R)/eta
    out[eta] = d
    print(f"eta={eta:.0e}: dA={ip['A']-info['A']:+.4e}  |dR/eta|max={np.abs(d).max():.8e}  at s<-15 mean strain coef={((d[s<-15]@mode)/(mode@mode)).mean():+.8e}")
ref = out[1e-6]
for eta, d in out.items():
    print(f"eta={eta:.0e}: |d - d(1e-6)|max = {np.abs(d-ref).max():.3e}")
