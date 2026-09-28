import numpy as np
from bq_logpolar import BQLogPolar
from bq_newton import residual, initial_guess
lam = 1.92
B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
X = initial_guess(B, (3 + lam) / 2)
R, info = residual(B, X)
Xn = X - R
print("A, m, vrmin:", info['A'], info['m'], info['vrmin'])
for sv in (-110, -60, -30, -20, -10, -5, -2, 0, 2, 5, 10, 20, 40, 60, 80, 95):
    j = np.argmin(abs(B.s - sv))
    print(f"s={sv:5.0f}: X(β=π/4)={X[j, B.Nb//2]:+.4e}  Xnew={Xn[j, B.Nb//2]:+.4e}  |R|max_β={np.abs(R[j]).max():.3e}   "
          f"Θ̂(π/4)={info['Th'][j, B.Nb//2]:+.3e} Ω̂(π/4)={info['Om'][j, B.Nb//2]:+.3e}")
