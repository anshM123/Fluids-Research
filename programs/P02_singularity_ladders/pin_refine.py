"""Refine a p=2 crossing (smooth self-similar profile) of the layer-pinned CCF branch in the sonic depth δ,
then compute the instability spectrum of the refined profile (dense method-of-lines generator on the same
graded grid).  Usage:  pin_refine.py START.npy DELTA_LO DELTA_HI [pts] [hs] [tag]"""
import numpy as np, sys, time
import scipy.linalg as sla
from scipy.interpolate import CubicSpline
from mapped_continuation import layer_info
from pinned_delta import Fpin, bordered, solve, new_grid

start, dlo, dhi = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
pts = int(sys.argv[4]) if len(sys.argv) > 4 else 24
hs = float(sys.argv[5]) if len(sys.argv) > 5 else 0.03
tag = sys.argv[6] if len(sys.argv) > 6 else "l3"
d = np.load(start)
eta0, phi0, lam = d[0], d[1], float(d[2][0])
sp = CubicSpline(eta0, phi0)
# centre = minimum of den on the source grid (source is a pinned solution: centre is a grid point)
from ccf_sinhgrid import CCFSinh
# rebuild a grid resolving the thinner end of the bracket (w ≈ 9.05 δ² near the cusp limit)
w_target = 9.05 * dlo ** 2
# find source centre: evaluate den on a temporary grid through the source points
Mt, it0 = new_grid(-0.96, 9.05 * dhi ** 2, pts, hs)
Ft, dent = Fpin(Mt, sp(Mt.eta), lam, it0)
etac = Mt.eta[int(np.argmin(dent))]
M, ic = new_grid(etac, w_target, pts, hs)
phi = sp(M.eta)
F, den = Fpin(M, phi, lam, ic)
delta0 = den[ic]
phi, lam, ok, hist, _ = solve(M, phi, lam, ic, delta0)
print(f"grid: centre {M.eta[ic]:.6f} hc={M.hc:.2e} N={M.N}; start δ={delta0:.5e} ok={ok} lam={lam:.12f} p={M.p_of(phi, lam):.12f}", flush=True)


def at_delta(dt, phi, lam, delta):
    F, den = Fpin(M, phi, lam, ic)
    A = bordered(M, phi, lam, den, ic)
    t = sla.solve(A, np.concatenate([np.zeros(M.N), [1.0]]), check_finite=False)
    # march in ≤15 % sub-steps
    cur_phi, cur_lam, cur = phi, lam, delta
    while abs(dt - cur) > 1e-15:
        nxt = dt if abs(dt / cur - 1) < 0.15 else cur * (0.85 if dt < cur else 1.15)
        ph, lm = cur_phi + (nxt - cur) * t[:-1], cur_lam + (nxt - cur) * t[-1]
        ph, lm, ok, hist, _ = solve(M, ph, lm, ic, nxt)
        if not ok:
            raise RuntimeError(f"solve failed at δ={nxt:.4e}")
        cur_phi, cur_lam, cur = ph, lm, nxt
        F, den = Fpin(M, cur_phi, cur_lam, ic)
        A = bordered(M, cur_phi, cur_lam, den, ic)
        t = sla.solve(A, np.concatenate([np.zeros(M.N), [1.0]]), check_finite=False)
    return cur_phi, cur_lam


# secant / regula falsi on g(δ) = p(δ) − 2 in ln δ
pts_list = []
for dt in (dhi, dlo):
    phi, lam = at_delta(dt, phi, lam, delta0 if not pts_list else pts_list[-1][0])
    pts_list.append((dt, lam, M.p_of(phi, lam), phi.copy()))
    print(f"  δ={dt:.6e} lam={lam:.12f} p-2={pts_list[-1][2]-2:+.3e}", flush=True)
a, b = pts_list[0], pts_list[1]
if (a[2] - 2) * (b[2] - 2) > 0:
    print("no sign change in bracket", flush=True)
    sys.exit(0)
cur = b
for k in range(30):
    la, lb = np.log(a[0]), np.log(b[0])
    lx = lb - (b[2] - 2) * (lb - la) / ((b[2] - 2) - (a[2] - 2))
    dx = np.exp(lx)
    phi, lam = at_delta(dx, cur[3], cur[1], cur[0])
    x = (dx, lam, M.p_of(phi, lam), phi.copy())
    print(f"  it {k}: δ={dx:.10e} lam={lam:.13f} p-2={x[2]-2:+.3e}", flush=True)
    if abs(x[2] - 2) < 1e-12:
        cur = x
        break
    if (x[2] - 2) * (a[2] - 2) < 0:
        b = x
    else:
        a = x
    cur = x
dstar, lam3, p3, phi3 = cur
F, den = Fpin(M, phi3, lam3, ic)
e, md, w = layer_info(M, den)
print(f"SMOOTH PROFILE ({tag}): lam={lam3:.13f} delta={dstar:.6e} w={w:.3e} layer eta={e:.6f} p={p3:.13f} "
      f"N={M.N} hc={M.hc:.2e} pts={pts} hs={hs}", flush=True)
np.save(f"pin_{tag}_pts{pts}_hs{hs}.npy", np.vstack([M.eta, phi3, np.full(M.N, lam3)]))

# ---- instability spectrum: dense MOL generator (3rd-order upwind in s, alternating-point Hilbert) -------
if len(sys.argv) > 7 and sys.argv[7] == "nospec":
    sys.exit(0)
t0 = time.time()
Th = np.exp(phi3 + M.c * M.eta); dd = 1 + lam3 + M.G(phi3); Th_eta = lam3 * Th / dd
pre, post = np.exp(-M.c * M.eta), np.exp((M.c - 1) * M.eta)
v = dd / (M.gp * M.hs); N = M.N
A = lam3 * np.eye(N)
for i in range(2, N - 1):
    A[i, i + 1] -= v[i] * 2 / 6; A[i, i] -= v[i] * 3 / 6; A[i, i - 1] -= v[i] * (-6 / 6); A[i, i - 2] -= v[i] * 1 / 6
A[1, 1] -= v[1]; A[1, 0] += v[1]; A[N - 1, N - 1] -= v[N - 1]; A[N - 1, N - 2] += v[N - 1]
A -= (post * Th_eta)[:, None] * (M.Hm * pre[None, :])
A[0, :] = 0.0; A[0, 0] = -50.0
ev = sla.eigvals(A, check_finite=False)
ev = ev[np.argsort(-ev.real)]
print(f"spectrum ({time.time()-t0:.0f}s): top eigenvalues (Re desc): " +
      ", ".join(f"{z.real:.5f}{z.imag:+.5f}i" for z in ev[:14]), flush=True)
real = ev[(np.abs(ev.imag) < 1e-6) & (ev.real > 0)]
print("real positive eigenvalues:", ", ".join(f"{z.real:.6f}" for z in real[:10]), flush=True)
np.save(f"pin_{tag}_spectrum.npy", ev)
