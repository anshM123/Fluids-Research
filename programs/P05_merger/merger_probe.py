"""IDEA005 LEVEL-2 probe: viscous merger of two co-rotating Gaussian vortices.
Measures separation b(t) of the two vorticity maxima and core size; merger time vs (a0/b0, Re)."""
import numpy as np, sys, time
sys.path.insert(0, '/home/user/Fluids-Research/tools')
from ps2d import PS2D
from scipy.ndimage import maximum_filter

def run(Re, a0, N=256, b0=1.0, Tmax=400.0, cfl=0.4):
    nu = 1.0 / Re                     # Gamma = 1 per vortex
    s = PS2D(N, nu=nu, workers=1)
    c = np.pi
    w = np.zeros_like(s.X)
    for sgn in (+1, -1):
        xc, yc = c + sgn * b0 / 2, c
        r2 = (s.X - xc)**2 + (s.Y - yc)**2
        w += 1.0 / (np.pi * a0**2) * np.exp(-r2 / a0**2)
    wh = s.fft(w)
    t, out = 0.0, []
    dt = min(s.cfl_dt(wh, cfl), 0.05)
    Trot = 2 * np.pi**2 * b0**2 / 1.0
    nextout = 0.0
    while t < Tmax:
        if t >= nextout:
            wf = s.ifft(wh)
            mx = (wf == maximum_filter(wf, size=9, mode='wrap')) & (wf > 0.2 * wf.max())
            idx = np.argwhere(mx)
            vals = wf[mx]
            order = np.argsort(-vals)
            pts = idx[order][:2] * (2 * np.pi / N)
            if len(pts) >= 2:
                d = pts[0] - pts[1]; d = (d + np.pi) % (2 * np.pi) - np.pi
                b = float(np.hypot(*d))
            else:
                b = 0.0
            out.append((t, b, float(wf.max())))
            if b < 0.05 * b0 and t > 1.0:
                break
            nextout += Trot / 20
        wh = s.step(wh, dt)
        t += dt
    return np.array(out)

if __name__ == "__main__":
    Re = float(sys.argv[1]); N = int(sys.argv[2])
    res = {}
    for a0 in [float(x) for x in sys.argv[3].split(',')]:
        t0 = time.time()
        o = run(Re, a0, N=N)
        # merger time: first time b < 0.5 b0
        tm = o[np.argmax(o[:, 1] < 0.5), 0] if np.any(o[:, 1] < 0.5) else np.nan
        # core size at merger onset estimate via diffusion law
        res[a0] = (tm, o)
        am = np.sqrt(a0**2 + 4 / Re * tm) if np.isfinite(tm) else np.nan
        print(f"Re={Re:.0f} N={N} a0={a0:.3f}: t_merge(b<0.5)={tm:.2f}  a_diff(t_m)/b0={am:.4f}  wallclock={time.time()-t0:.0f}s", flush=True)
        np.save(f"merger_Re{Re:.0f}_a{a0:.3f}_N{N}.npy", o)
