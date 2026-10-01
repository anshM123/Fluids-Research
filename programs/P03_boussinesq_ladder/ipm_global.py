"""INDEPENDENT solver for smooth self-similar IPM profiles (cross-check of the marching solver ipm_solver.py), built
on the 2D Boussinesq global solver bq_global.py and different from the marching solver in every numerical
ingredient:
  * unknowns are the unhatted fields R, Ω and X = Ψ/r² on the (s, β) grid, all solved simultaneously (no marching
    along characteristics, no stagnation-point factorization R = cos^m β R̂);
  * s-derivatives by 6th-order finite differences (upwind-biased for the transport), β by Chebyshev collocation;
  * Biot–Savart as the local elliptic equation (∂_s + 2)²X + X_ββ + Ω = 0 (sparse, no FFT);
  * smoothness imposed directly (R ≈ +y₁², hence Ω = −∂₁R ≈ −2y₁, m = 2; the sign of ipm_solver, compressive
    strain A > 0), with λ an unknown eigenvalue fixed by the
    strain condition A = −∂₁U₁(0) = 1 + λ/2 (m = λ/(1 + λ − A) = 2);
  * Newton's method with the exact sparse Jacobian and a sparse direct solve.

Equations (D = (1+λ) + X_β, w = −(2X + X_s); ∂₁ = e^{−s}(cos β ∂_s − sin β ∂_β)):
    −λR + D R_s + w R_β = 0                     (V·∇R = λR)
    Ω + ∂₁R = 0                                 (IPM: Ω = −∂₁R, no transport equation for Ω)
    X_ss + 4X_s + 4X + X_ββ + Ω = 0
Boundary conditions: s_min: R = e^{2s}cos²β (inflow), X_s = e^{s} cos β sin²β (no singular harmonics for
Ω = −2e^{s}cos β); s_max: X_s + X/(1+λ) = 0; β = 0, π/2: X = 0.  Extra equation: −X_β(s_min, 0) − (1 + λ/2) = 0."""
import numpy as np
import scipy.sparse as sp
from bq_global import BQGlobal


class IPMGlobal(BQGlobal):
    def local_data(self, lam):
        s0 = self.s[0]
        R0 = np.exp(2 * s0) * np.cos(self.beta) ** 2                          # R = +y₁² (sign of ipm_solver, A > 0)
        Xs0 = np.exp(s0) * np.cos(self.beta) * np.sin(self.beta) ** 2         # for Ω = −2y₁
        return R0, Xs0

    def residual(self, U):
        R, Om, X, lam = self.unpack(U)
        Rs, Rb = self.Ds @ R, self.Dbeta @ R
        Xs, Xb = self.Ds @ X, self.Dbeta @ X
        D = (1 + lam) + Xb
        w = -(2 * X + Xs)
        Ru = self.Du @ R
        Rt = -lam * R + D * Ru + w * Rb
        Ro = Om + self.eS * (self.CB * Rs - self.SB * Rb)
        Rx = self.Dss @ X + 4 * Xs + 4 * X + self.Dbeta2 @ X + Om
        R0, Xs0 = self.local_data(lam)
        Rt[self.inflow] = R[self.inflow] - R0
        Rx[self.xbc] = X[self.xbc]
        Rx[self.xreg] = (Xs - np.tile(Xs0, self.Ns))[self.xreg]
        Rx[self.xfar] = (Xs + X / (1 + lam))[self.xfar]
        Rl = -Xb[self.k0_smin] - (1 + lam / 2)
        return np.concatenate([Rt, Ro, Rx, [Rl]])

    def jacobian(self, U):
        N = self.N
        R, Om, X, lam = self.unpack(U)
        Rb = self.Dbeta @ R
        Xs, Xb = self.Ds @ X, self.Dbeta @ X
        D = (1 + lam) + Xb
        w = -(2 * X + Xs)
        dg = sp.diags
        Ru = self.Du @ R
        Jtt = -lam * sp.identity(N) + dg(D) @ self.Du + dg(w) @ self.Dbeta
        Jtx = dg(Ru) @ self.Dbeta + dg(Rb) @ (-2 * sp.identity(N) - self.Ds)
        Jtl = -R + Ru
        Joo = sp.identity(N)
        Jot = dg(self.eS * self.CB) @ self.Ds - dg(self.eS * self.SB) @ self.Dbeta
        Jxx = self.Dss + 4 * self.Ds + 4 * sp.identity(N) + self.Dbeta2
        Jxo = sp.identity(N)
        Z = sp.csr_matrix((N, N))
        Jtt, Jtx = Jtt.tolil(), Jtx.tolil()
        Jxx, Jxo = Jxx.tolil(), Jxo.tolil()
        Jtl = Jtl.copy(); Jxl = np.zeros(N)
        for i in np.where(self.inflow)[0]:
            Jtt.rows[i] = [i]; Jtt.data[i] = [1.0]; Jtx.rows[i] = []; Jtx.data[i] = []
        Jtl[self.inflow] = 0.0
        Dsl = self.Ds.tolil()
        for i in np.where(self.xbc)[0]:
            Jxx.rows[i] = [i]; Jxx.data[i] = [1.0]; Jxo.rows[i] = []; Jxo.data[i] = []
        for i in np.where(self.xfar)[0]:
            rr = list(Dsl.rows[i]); dd = list(Dsl.data[i])
            if i in rr:
                dd[rr.index(i)] += 1.0 / (1 + lam)
            else:
                rr.append(i); dd.append(1.0 / (1 + lam))
                o = np.argsort(rr); rr = [int(rr[q]) for q in o]; dd = [dd[q] for q in o]
            Jxx.rows[i] = rr; Jxx.data[i] = dd; Jxo.rows[i] = []; Jxo.data[i] = []
            Jxl[i] = -X[i] / (1 + lam) ** 2
        for i in np.where(self.xreg)[0]:
            Jxx.rows[i] = list(Dsl.rows[i]); Jxx.data[i] = list(Dsl.data[i]); Jxo.rows[i] = []; Jxo.data[i] = []
        J = sp.bmat([[Jtt, Z, Jtx, sp.csr_matrix(Jtl[:, None])],
                     [Jot, Joo, Z, None],
                     [Z, Jxo, Jxx, sp.csr_matrix(Jxl[:, None])],
                     [None, None, sp.csr_matrix(-self.Dbeta[self.k0_smin].toarray()), sp.csr_matrix([[-0.5]])]],
                    format='csc')
        return J

    def initial_from_march(self, f, lam, Nb_src=32, hs_src=0.0125, s_start=-20.0):
        """interpolate a converged marching-solver state (ipm_solver) onto this grid"""
        from ipm_solver import IPM
        from bq_newton import full
        from numpy.polynomial import chebyshev as Cb
        B = IPM(lam, Nb=Nb_src, hs=hs_src, s_sw=12.0, s_start=s_start)
        Y = np.load(f); Xm = full(B, Y)
        r = B.march(Xm / B.ea2[:, None], return_all=True)
        cp = np.where(B.cb > 0, B.cb, 0.0)
        Rm = r['Th'] * cp[None, :] ** r['m']; Omm = r['Omega']
        xo = np.cos(np.pi * np.arange(B.Nb + 1) / B.Nb); xn = np.cos(np.pi * np.arange(self.Nb + 1) / self.Nb)
        out = []
        for F in (Rm, Omm, Xm):
            Fb = Cb.chebval(xn, Cb.chebfit(xo, F.T, B.Nb))
            Fi = np.array([np.interp(self.s, B.s, Fb[:, k]) for k in range(self.Nb + 1)]).T
            out.append(Fi.ravel())
        return np.concatenate(out + [[lam]])


if __name__ == "__main__":
    # Jacobian check against finite differences on a small grid
    G = IPMGlobal(s_min=-4.0, s_max=2.0, hs=0.1, Nb=8)
    rng = np.random.default_rng(1)
    U = np.concatenate([rng.standard_normal(3 * G.N) * 0.1, [0.4]])
    J = G.jacobian(U).toarray(); R0 = G.residual(U)
    Jfd = np.zeros_like(J)
    for k in range(len(U)):
        e = np.zeros_like(U); e[k] = 1e-6
        Jfd[:, k] = (G.residual(U + e) - G.residual(U - e)) / 2e-6
    print("max |J − J_fd| =", np.abs(J - Jfd).max(), " (max |J| =", np.abs(J).max(), ")")
