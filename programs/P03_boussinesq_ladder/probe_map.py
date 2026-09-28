import numpy as np, sys
from bq_logpolar import BQLogPolar
lam = 1.92
B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
s = B.s[:, None]; b = B.beta[None, :]; r = np.exp(s)
for A0 in (1.5, 2.0, 2.46, 3.0, 4.0):
    for amp_scale in (1.0,):
        P = np.exp(-B.a*s) * (-(A0/2) * r**2 * np.sin(2*b) / (1 + r**2) ** (1/(2*(1+lam))))
        r1 = B.march(P, return_all=True)
        Pn = B.poisson(r1['Omega'])
        Ur_n = B.velocity(Pn)[0]
        An = B.strain(Ur_n)
        j = np.argmin(abs(B.s))
        print(f"A_in={r1['A']:.4f} m={r1['m']:.4f} min(V_r/r)={r1['vrmin']:.3f} → A_out={An:.6f};  Θ̂(s=0) range [{r1['Th'][j].min():.3e},{r1['Th'][j].max():.3e}]  Ω̂(s=0) [{r1['Om'][j].min():.3e},{r1['Om'][j].max():.3e}]", flush=True)
