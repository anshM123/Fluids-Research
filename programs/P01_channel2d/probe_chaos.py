"""Probe: do finite-amplitude disturbances in 2D plane Poiseuille flow produce sustained chaos? (IDEA001, LEVEL 2)"""
import numpy as np, sys, time, json
sys.path.insert(0, '/home/user/Fluids-Research/tools')
from channel2d import Channel2D
Re = float(sys.argv[1]); Lx = float(sys.argv[2]); amp = float(sys.argv[3]); seed = int(sys.argv[4]); T = float(sys.argv[5])
Nx = int(sys.argv[6]) if len(sys.argv) > 6 else int(32 * Lx / (2*np.pi))
M = int(sys.argv[7]) if len(sys.argv) > 7 else 64
ch = Channel2D(Nx=Nx, M=M, Lx=Lx, Re=Re)
dt = 0.01; ch.set_dt(dt)
a = ch.random_perturbation(amp, kmax=6.0, seed=seed)
out = []; t0 = time.time()
nsteps = int(T/dt)
for n in range(nsteps+1):
    if n % 100 == 0:
        e = ch.perturbation_energy(a); ws = ch.wall_shear(a)
        out.append((n*dt, e, float(ws[0]), float(ws[1])))
        if not np.isfinite(e): break
    a = ch.step(a)
np.save(f"chaos_Re{Re:.0f}_L{Lx:.2f}_A{amp:g}_s{seed}.npy", np.array(out))
print(f"Re={Re} Lx={Lx:.2f} amp={amp} seed={seed} Nx={Nx} M={M}: final E={out[-1][1]:.3e}, wall-clock {time.time()-t0:.0f}s")
