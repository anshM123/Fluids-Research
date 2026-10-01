"""Linear stability of self-similar IPM profiles through the eigen-condition ν(μ) = 1 of the linearised
"march + Biot–Savart" map T_μ. Same method as bq_stability.py for Boussinesq, with the vorticity slaved to the
density gradient.

Perturbations (r, ω', ψ') e^{μτ} of a profile (R, Ω, Ψ), τ = −ln(1−t), satisfy
    (μ − λ) r + V·∇r + u'·∇R = 0,      ω' = −∂₁r,      u' = ∇^⊥ψ',  −Δψ' = ω',      V = (1+λ)y + U.
All characteristics of V leave the stagnation point. For a given velocity perturbation X' = ψ'/r², r is marched
outward from the origin with vanishing data, ω' follows algebraically, and Biot–Savart returns X'_new; μ is an
eigenvalue iff T_μ has the eigenvalue ν = 1.

Log-polar form (D = V_r/r, w = V_β/r, U'_r/r = X'_β, U'_β/r = −(2X' + X'_s)):
    D r_s = −(μ − λ) r − w r_β − (U'_r/r) R_s − (U'_β/r) R_β,     ω' = −e^{−s}(cos β r_s − sin β r_β),
with r_s in ω' taken from the transport equation itself. Trivial modes: μ = 1 (time translation,
r = λR − (1+λ)R_s) and μ = 0 (dilation; a Jordan block with ∂_λ of the profile)."""
import numpy as np
from numba import njit
from scipy.sparse.linalg import LinearOperator, eigs
from bq_solver import _grid_deriv, _offset_values, C1, C2
from bq_newton import full
from bq_stability import BQStab
from ipm_solver import IPM


@njit(cache=True)
def _lin_march_r(s, h, mu, lam, Db, i0, jsw, D0, w0, D1, w1, D2, w2, Dh, wh, Dg1, wg1,
                 F0, F1, F2, Fh, Fg1):
    """complex linear march of r (base speeds real, forcings complex) — Gauss–Legendre for s < s_sw, RK4 beyond"""
    Ns, n = D0.shape
    R = np.zeros((Ns, n), dtype=np.complex128)
    Dbc = Db.astype(np.complex128)
    I = np.eye(n).astype(np.complex128)
    a11, a12, a21, a22 = 0.25, 0.25 - np.sqrt(3.0) / 6, 0.25 + np.sqrt(3.0) / 6, 0.25
    M = np.zeros((2 * n, 2 * n), dtype=np.complex128); rhs = np.zeros(2 * n, dtype=np.complex128)
    for j in range(i0, Ns - 1):
        r = R[j].copy()
        if j < jsw:
            P1 = np.empty((n, n), dtype=np.complex128); P2 = np.empty((n, n), dtype=np.complex128)
            for i in range(n):
                for k in range(n):
                    P1[i, k] = -w1[j, i] * Db[i, k] / D1[j, i]
                    P2[i, k] = -w2[j, i] * Db[i, k] / D2[j, i]
                P1[i, i] += -(mu - lam) / D1[j, i]
                P2[i, i] += -(mu - lam) / D2[j, i]
            M[:n, :n] = I - h * a11 * P1; M[:n, n:] = -h * a12 * P1
            M[n:, :n] = -h * a21 * P2;    M[n:, n:] = I - h * a22 * P2
            rhs[:n] = P1 @ r + F1[j] / D1[j]; rhs[n:] = P2 @ r + F2[j] / D2[j]
            K = np.linalg.solve(M, rhs)
            R[j + 1] = r + 0.5 * h * (K[:n] + K[n:])
        else:
            k1 = (-(mu - lam) * r - w0[j] * (Dbc @ r) + F0[j]) / D0[j]
            r2 = r + 0.5 * h * k1
            k2 = (-(mu - lam) * r2 - wh[j] * (Dbc @ r2) + Fh[j]) / Dh[j]
            r3 = r + 0.5 * h * k2
            k3 = (-(mu - lam) * r3 - wh[j] * (Dbc @ r3) + Fh[j]) / Dh[j]
            r4 = r + h * k3
            k4 = (-(mu - lam) * r4 - wg1[j] * (Dbc @ r4) + Fg1[j]) / Dg1[j]
            R[j + 1] = r + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    return R


class IPMStab(BQStab):
    def __init__(self, lam, Y, Nb=32, hs=0.025, s_sw=12.0, s_start=-20):
        self.B = B = IPM(lam, Nb=Nb, hs=hs, s_sw=s_sw, s_start=s_start)
        self.lam = lam
        X = full(B, Y)
        r = B.march(X / B.ea2[:, None], return_all=True)
        self.m, self.A, self.vrmin = r['m'], r['A'], r['vrmin']
        cp = np.where(B.cb > 0, B.cb, 0.0)
        self.Th = r['Th'] * cp[None, :] ** self.m           # R (unhatted)
        self.Om = r['Omega']
        h = B.hs; Db = B.Db
        V = B.velocity(X / B.ea2[:, None])
        pad = lambda Z: np.vstack([Z, Z[-1:]])
        self.D0 = (1 + lam) + V[0]; self.w0 = V[1]
        self.Dh = pad((1 + lam) + V[3]); self.wh = pad(V[4])
        self.D1 = pad((1 + lam) + V[6]); self.w1 = pad(V[7])
        self.D2 = pad((1 + lam) + V[9]); self.w2 = pad(V[10])
        self.Dg1 = np.vstack([self.D0[1:], self.D0[-1:]]); self.wg1 = np.vstack([self.w0[1:], self.w0[-1:]])

        def grads(F):
            Fs = _grid_deriv(F, h); Fb = F @ Db.T
            out = {'0': (Fs, Fb)}
            for key, c in (('h', 0.5), ('1', C1), ('2', C2)):
                Vv, Dd = _offset_values(F, h, c)
                out[key] = (pad(Dd), pad(Vv @ Db.T))
            return out
        self.gT = grads(self.Th)
        self.i0, self.jsw = B.i0, B.jsw
        self.n_unk = (B.Ns - B.i0) * (B.Nb + 1)
        self.ems = np.exp(-B.s)[:, None]

    def _forcing(self, Xp):
        u = self._pert_vel(Xp)
        F = {}
        for key in ('0', 'h', '1', '2'):
            Ur, w = u[key]
            Ts, Tb = self.gT[key]
            F[key] = -(Ur * Ts + w * Tb)
        return F

    def Tc(self, Yp, mu):
        B = self.B
        Yp = np.asarray(Yp, dtype=complex).reshape(B.Ns - B.i0, B.Nb + 1)
        Fr = self._forcing(full(B, Yp.real)); Fi = self._forcing(full(B, Yp.imag))
        F = {k: Fr[k] + 1j * Fi[k] for k in Fr}
        sh = lambda Z: np.vstack([Z[1:], Z[-1:]])
        mu = complex(mu)
        r = _lin_march_r(B.s, B.hs, mu, self.lam, B.Db, self.i0, self.jsw, self.D0, self.w0, self.D1, self.w1,
                         self.D2, self.w2, self.Dh, self.wh, self.Dg1, self.wg1,
                         F['0'], F['1'], F['2'], F['h'], sh(F['0']))
        rb = r @ B.Db.T
        rs = (-(mu - self.lam) * r - self.w0 * rb + F['0']) / self.D0
        om = -self.ems * (B.cb[None, :] * rs - B.sb[None, :] * rb)
        om[:B.i0] = 0.0
        Pn = B.poisson(om.real) + 1j * B.poisson(om.imag)
        return (B.ea2[:, None] * Pn)[B.i0:].ravel()

    def T(self, Yp, mu):
        return self.Tc(Yp, mu).real

    def time_translation_mode(self):
        """X' of the μ = 1 mode: r = λR − (1+λ)R_s, ω' = −∂₁r"""
        B = self.B
        Rs, Rb = self.gT['0']
        r = self.lam * self.Th - (1 + self.lam) * Rs
        rb = r @ B.Db.T; rs = _grid_deriv(r, B.hs)
        om = -self.ems * (B.cb[None, :] * rs - B.sb[None, :] * rb)
        om[:B.i0] = 0.0
        return (B.ea2[:, None] * B.poisson(om))[B.i0:].ravel()


def spectrum_c(S, mu, k=10, tol=1e-8, v0=None):
    op = LinearOperator((S.n_unk, S.n_unk), matvec=lambda v: S.Tc(v, mu), dtype=complex)
    vals, vecs = eigs(op, k=k, which='LM', tol=tol, v0=v0, maxiter=3000)
    o = np.argsort(-np.abs(vals))
    return vals[o], vecs[:, o]


if __name__ == "__main__":
    import sys, time
    f, lam = sys.argv[1], float(sys.argv[2])
    t0 = time.time()
    S = IPMStab(lam, np.load(f))
    print(f"IPM λ={lam}: m={S.m:.8f} A={S.A:.8f} vrmin={S.vrmin:.4f}", flush=True)
    x1 = S.time_translation_mode()
    Tx = S.T(x1, 1.0)
    print(f"time-translation check: |T_1 x − x|/|x| = {np.linalg.norm(Tx - x1) / np.linalg.norm(x1):.2e} "
          f"({time.time() - t0:.0f}s)", flush=True)
    for mu in (1.0, 0.5, 0.2):
        v, _ = spectrum_c(S, mu, k=8)
        print(f"μ={mu}: top ν = {np.round(v, 5).tolist()} ({time.time() - t0:.0f}s)", flush=True)
