"""Independent check: direct time integration of the CCF equation θ_t + (Hθ) θ_x = 0 on a periodic domain
(pseudo-spectral, FFT Hilbert transform, RK4, CFL-adaptive dt), from smooth even data with θ(0)=0.
Measures self-similar scaling as t→T:  L(t) (location of max θ_x)  ∝ (T−t)^{1+λ},  θ at fixed ξ ∝ (T−t)^λ.
Local exponent estimate:  1+λ = d ln L / d ln(1/max θ_x)  (since max θ_x ∝ (T−t)^{-1})."""
import numpy as np, sys, time
import scipy.fft as sfft


def run(N=2**16, T_end=10.0, cfl=0.25, amp=1.0, out="td"):
    x = -np.pi + 2 * np.pi * np.arange(N) / N
    k = sfft.fftfreq(N, 1.0 / N)
    sgn = -1j * np.sign(k)            # periodic Hilbert: H e^{ikx} = -i sgn(k) e^{ikx}
    ik = 1j * k
    kmax = N // 2
    filt = np.exp(-36 * (np.abs(k) / kmax) ** 36)   # Hou–Li spectral filter
    th = amp * (1 - np.cos(x))
    th_h = sfft.fft(th)

    def rhs(th_h):
        u = np.real(sfft.ifft(sgn * th_h))
        thx = np.real(sfft.ifft(ik * th_h))
        return -sfft.fft(u * thx), np.max(np.abs(u)), thx
    t = 0.0
    rec = []
    while t < T_end:
        r1, umax, thx = rhs(th_h)
        g = np.max(np.abs(thx))
        dt = cfl * (2 * np.pi / N) / max(umax, 1e-12)
        dt = min(dt, 0.05 / g)             # resolve the 1/(T-t) growth
        r2, _, _ = rhs(th_h + 0.5 * dt * r1)
        r3, _, _ = rhs(th_h + 0.5 * dt * r2)
        r4, _, _ = rhs(th_h + dt * r3)
        th_h = (th_h + dt / 6 * (r1 + 2 * r2 + 2 * r3 + r4)) * filt
        t += dt
        i = int(np.argmax(np.abs(thx)))
        rec.append((t, g, abs(x[i]), np.real(sfft.ifft(th_h))[N // 2]))
        # stop when the gradient scale reaches ~8 grid spacings
        if abs(x[i]) < 8 * 2 * np.pi / N:
            break
        # spectral tail check
        tail = np.abs(th_h[kmax - kmax // 8:kmax]).max() / np.abs(th_h).max()
        if tail > 1e-10:
            break
    rec = np.array(rec)
    np.save(f"{out}_N{N}.npy", rec)
    return rec


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2**16
    t0 = time.time()
    rec = run(N)
    t, g, L = rec[:, 0], rec[:, 1], rec[:, 2]
    print(f"N={N}: steps={len(rec)} t_end={t[-1]:.8f} max|θ_x|={g[-1]:.3e} L={L[-1]:.3e} ({time.time()-t0:.0f}s)")
    # local exponent 1+λ from d lnL / d ln(1/g) over the last decades
    m = g > 20
    lg, lL = np.log(1 / g[m]), np.log(L[m])
    for lo, hi in [(1e1, 1e2), (1e2, 1e3), (1e3, 1e4), (1e4, 1e5)]:
        s = (g[m] > lo) & (g[m] < hi)
        if s.sum() > 20:
            p = np.polyfit(lg[s], lL[s], 1)
            print(f"   max|θ_x| in [{lo:.0e},{hi:.0e}]: fitted 1+λ = {p[0]:.4f}  → λ = {p[0]-1:.4f}")
