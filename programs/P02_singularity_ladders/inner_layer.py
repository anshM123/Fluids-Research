"""Universal inner-layer problem of the sonic cusp (IDEA117).
In the layer, den = δ f(X), X = (η−η_s)/ℓ, ℓ = δ²/k² (k: cusp coefficient, den → k|η−η_s|^{1/2}), and
      f(X) = 1 + ½ [H F(X) − H F(0)],   F' = 1/f,  F odd, f even,
i.e.  f(X) − 1 = (1/π) p.v.∫_0^∞ F(y) X² / (y (X² − y²)) dy,   f → |X|^{1/2} as |X| → ∞ (automatic).
Log variable X = e^t:  f(e^t) − 1 = (1/π) p.v.∫ F(e^s) k(t−s) ds,  k(τ) = 1/(1 − e^{−2τ}).
Output: f''(0) (curvature at the minimum) ⇒ predicted layer width  w = ℓ/√f''(0) = δ² /(k² √f''(0)).
"""
import numpy as np
from scipy.optimize import newton_krylov

T0, T1, h = -18.0, 22.0, 0.01
t = np.arange(T0, T1 + h / 2, h); N = t.size
I = np.arange(N)
D = t[:, None] - t[None, :]
mask = ((I[:, None] - I[None, :]) % 2) == 1
with np.errstate(over="ignore", divide="ignore"):
    Kmat = np.where(mask, 1.0 / (-np.expm1(-2.0 * np.where(mask, D, 1.0))), 0.0)
Kmat *= 2 * h / np.pi
nodes = np.arange(-3, 5); rhs = np.array([1.0 / (m + 1) for m in range(8)])
w8 = np.linalg.solve(np.vander(nodes, 8, increasing=True).T, rhs)
wl = [np.linalg.solve(np.vander(np.arange(-s, 8 - s), 8, increasing=True).T, rhs) for s in range(8)]


def cumint(Y):
    mid = np.zeros(N - 1)
    for j, off in enumerate(range(-3, 5)):
        mid[3:N - 4] += w8[j] * Y[3 + off:N - 4 + off]
    for i in list(range(0, 3)) + list(range(N - 4, N - 1)):
        s = i if i < 3 else 7 - (N - 1 - i)
        mid[i] = wl[s] @ Y[i - s:i - s + 8]
    C = np.zeros(N); C[1:] = np.cumsum(mid * h)
    return C


X = np.exp(t)
# left tail (s < T0): F ≈ X (f≈1), kernel ≈ 1 → ∫_{-∞}^{T0} e^s ds = e^{T0}; right tail (s > T1): F ≈ 2X^{1/2},
# k ≈ −e^{2(t−s)} → −(4/3) e^{2t} e^{−3T1/2}
tailL = np.exp(T0) / np.pi
tailR = -(4.0 / 3.0) * np.exp(2 * t) * np.exp(-1.5 * T1) / np.pi


def Fof(g):                      # g = ln f
    return np.exp(T0) + cumint(np.exp(t - g))


def resid(g):
    F = Fof(g)
    rhsv = 1.0 + Kmat @ F + tailL + tailR
    return g - np.log(np.maximum(rhsv, 1e-300))


g0 = 0.5 * np.log1p(X)            # initial guess f = √(1+X)
g = newton_krylov(resid, g0, f_tol=1e-11, method="lgmres", verbose=False, maxiter=200)
f = np.exp(g)
r = resid(g)
print(f"N={N}, max|residual|={np.abs(r[(t > T0 + 2) & (t < T1 - 6)]).max():.2e}")
sel = (t > -9) & (t < -5)
c2 = np.polyfit(X[sel] ** 2, f[sel] - 1.0, 1)[0]
print(f"f''(0) = 2·c2 = {2*c2:.8f}")
for tt in (8, 10, 12):
    j = int(np.argmin(np.abs(t - tt)))
    print(f"   X=e^{tt}: f − X^(1/2) = {f[j]-np.sqrt(X[j]):.6f}")
np.save("inner_layer_f.npy", np.vstack([t, f]))
print(f"prediction: w/δ² = 1/(k² √f''(0)) = {1/np.sqrt(2*c2):.6f}/k²")
