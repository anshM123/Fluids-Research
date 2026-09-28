"""Transfer a solution Y = X[i0:] (X = Ψ/r²) between log-polar grids (different Nb, hs, s range) and re-solve.
Usage: bq_regrid.py Y.npy LAM NB_OLD HS_OLD NB_NEW HS_NEW [S_MIN S_MAX S_START] -> prints A, m; saves Y_new"""
import numpy as np, sys, time
from numpy.polynomial import chebyshev as C
from scipy.interpolate import CubicSpline
from bq_logpolar import BQLogPolar
from bq_newton import newton, full


def transfer(Bo, Yo, Bn):
    Xo = full(Bo, Yo)                                   # whole old grid
    xo = np.cos(np.pi * np.arange(Bo.Nb + 1) / Bo.Nb)
    xn = np.cos(np.pi * np.arange(Bn.Nb + 1) / Bn.Nb)
    coef = C.chebfit(xo, Xo.T, Bo.Nb)                   # exact interpolation on the Lobatto nodes
    Xb = C.chebval(xn, coef)                            # (Ns_old, Nb_new+1)
    cs = CubicSpline(Bo.s, Xb, axis=0)
    sn = np.clip(Bn.s, Bo.s[0], Bo.s[-1])
    Xn = cs(sn)
    Xn[:, 0] = 0.0; Xn[:, -1] = 0.0
    return Xn[Bn.i0:]


if __name__ == "__main__":
    f, lam = sys.argv[1], float(sys.argv[2])
    nbo, hso, nbn, hsn = int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6])
    smin, smax, sst = (float(sys.argv[7]), float(sys.argv[8]), float(sys.argv[9])) if len(sys.argv) > 9 else (-120., 100., -20.)
    Bo = BQLogPolar(lam, s_min=-120, s_max=100, hs=hso, Nb=nbo)
    Bn = BQLogPolar(lam, s_min=smin, s_max=smax, hs=hsn, Nb=nbn, s_start=sst)
    Y = transfer(Bo, np.load(f), Bn)
    t0 = time.time()
    Y, info, ok = newton(Bn, Y, tol=1e-6, maxit=15, verbose=True, t0=t0)
    print(f"RESULT lam={lam} Nb={nbn} hs={hsn} s=[{smin},{smax}] s_start={sst}: ok={ok} A={info['A']:.10f} m={info['m']:.10f}",
          flush=True)
    np.save(f"Yr_lam{lam:.4f}_Nb{nbn}_hs{hsn}.npy", Y)
