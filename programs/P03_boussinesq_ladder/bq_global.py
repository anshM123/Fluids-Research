"""
INDEPENDENT solver for smooth self-similar 2D Boussinesq profiles (cross-check of the marching solver).

Different in every numerical ingredient from bq_solver.py:
  * unknowns are the UNHATTED fields Θ, Ω and X = Ψ/r² on the (s, β) grid, all solved simultaneously
    (no marching along characteristics);
  * s-derivatives by 6th-order central finite differences (one-sided near the ends), β by Chebyshev collocation;
  * Biot–Savart as the local elliptic equation (∂_s + 2)²X + X_ββ + Ω = 0 (sparse, no FFT, no exponential
    weights), X = 0 on the wall and the axis;
  * smoothness imposed directly (Θ ≈ −y₁², Ω ≈ −C y₁ with C = 4/(1+λ), i.e. m = 2), and λ an unknown
    eigenvalue fixed by the strain condition A = −∂₁U₁(0) = (3+λ)/2;
  * Newton's method with the exact sparse Jacobian and a sparse direct (SuperLU) solve.

Equations (D = (1+λ) + X_β, w = −(2X + X_s)):
    (1−λ)Θ + D Θ_s + w Θ_β = 0
    Ω + D Ω_s + w Ω_β − e^{−s}(cosβ Θ_s − sinβ Θ_β) = 0
    X_ss + 4X_s + 4X + X_ββ + Ω = 0
Boundary conditions:
  * s_min: exact smooth local data Θ = −e^{2s}cos²β, Ω = −C e^{s} cosβ (inflow);
    regularity X_s = (C/2) e^{s} cosβ sin²β (no singular harmonics);
  * s_max: transport outflow (no condition); X = 0 (decay);
  * β = 0, π/2: X = 0.
Extra unknown λ; extra equation: −X_β(s_min, β=0) − (3+λ)/2 = 0.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from bq_logpolar import cheb


def fd_matrix(N, h, order=6, der=1, shift=0):
    """sparse FD matrix for the der-th derivative: (2w+1)-point stencils centred at i − shift inside (shift = 1:
    upwind-biased for flow towards +s), shifted one-sided stencils near the ends
    (w = order/2 for der = 1, order/2 + 1 for der = 2)"""
    import math
    w = order // 2 + (1 if der == 2 else 0)
    rows, cols, vals = [], [], []
    for i in range(N):
        lo = min(max(i - w - shift, 0), N - 2 * w - 1)
        idx = np.arange(lo, lo + 2 * w + 1)
        x = (idx - i).astype(float)
        V = np.vander(x, len(x), increasing=True).T
        rhs = np.zeros(len(x)); rhs[der] = math.factorial(der)
        c = np.linalg.solve(V, rhs) / h ** der
        rows += [i] * len(idx); cols += list(idx); vals += list(c)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N))


class BQGlobal:
    def __init__(self, s_min=-8.0, s_max=40.0, hs=0.025, Nb=24, order=6):
        self.Ns = int(round((s_max - s_min) / hs)) + 1
        self.hs = hs
        self.s = s_min + hs * np.arange(self.Ns)
        D, x = cheb(Nb)
        self.Nb = Nb
        self.beta = (np.pi / 4) * (1 - x)
        self.Db = -(4.0 / np.pi) * D
        self.Dbb = self.Db @ self.Db
        self.cb = np.cos(self.beta); self.cb[-1] = 0.0; self.sb = np.sin(self.beta)
        n = Nb + 1
        self.n = n
        Is = sp.identity(self.Ns, format='csr'); Ib = sp.identity(n, format='csr')
        self.Ds = sp.kron(fd_matrix(self.Ns, hs, order, 1), Ib, format='csr')      # ∂_s on the flattened grid
        self.Du = sp.kron(fd_matrix(self.Ns, hs, order, 1, shift=1), Ib, format='csr')   # upwind-biased ∂_s
        self.Dss = sp.kron(fd_matrix(self.Ns, hs, order, 2), Ib, format='csr')
        self.Dbeta = sp.kron(Is, sp.csr_matrix(self.Db), format='csr')             # ∂_β
        self.Dbeta2 = sp.kron(Is, sp.csr_matrix(self.Dbb), format='csr')
        self.N = self.Ns * n
        S, Bt = np.meshgrid(self.s, self.beta, indexing='ij')
        self.S = S.ravel(); self.CB = np.cos(Bt).ravel(); self.SB = np.sin(Bt).ravel()
        self.CB[np.isclose(Bt.ravel(), np.pi / 2)] = 0.0
        self.eS = np.exp(-self.S)
        # index sets
        j = np.arange(self.N) // n; k = np.arange(self.N) % n
        self.inflow = j == 0
        self.xbc = (k == 0) | (k == n - 1)                                          # Dirichlet X = 0 (wall, axis)
        self.xfar = (j == self.Ns - 1) & ~((k == 0) | (k == n - 1))                 # far field: X_s + X/(1+λ) = 0
        self.xreg = (j == 0) & ~((k == 0) | (k == n - 1))                            # regularity at s_min
        self.k0_smin = 0                                                             # index of (s_min, β=0)

    # ---------------- residual and Jacobian ----------------
    def unpack(self, U):
        N = self.N
        return U[:N], U[N:2 * N], U[2 * N:3 * N], U[3 * N]

    def local_data(self, lam):
        C = 4.0 / (1.0 + lam)
        s0 = self.s[0]
        Th0 = -np.exp(2 * s0) * np.cos(self.beta) ** 2
        Om0 = -C * np.exp(s0) * np.cos(self.beta)
        Xs0 = 0.5 * C * np.exp(s0) * np.cos(self.beta) * np.sin(self.beta) ** 2
        return C, Th0, Om0, Xs0

    def residual(self, U):
        Th, Om, X, lam = self.unpack(U)
        Ths, Thb = self.Ds @ Th, self.Dbeta @ Th
        Oms, Omb = self.Ds @ Om, self.Dbeta @ Om
        Xs, Xb = self.Ds @ X, self.Dbeta @ X
        D = (1 + lam) + Xb
        w = -(2 * X + Xs)
        Thu, Omu = self.Du @ Th, self.Du @ Om
        Rt = (1 - lam) * Th + D * Thu + w * Thb
        Ro = Om + D * Omu + w * Omb - self.eS * (self.CB * Ths - self.SB * Thb)
        Rx = self.Dss @ X + 4 * Xs + 4 * X + self.Dbeta2 @ X + Om
        C, Th0, Om0, Xs0 = self.local_data(lam)
        Rt[self.inflow] = Th[self.inflow] - Th0
        Ro[self.inflow] = Om[self.inflow] - Om0
        Rx[self.xbc] = X[self.xbc]
        Rx[self.xreg] = (Xs - np.tile(Xs0, self.Ns))[self.xreg]
        Rx[self.xfar] = (Xs + X / (1 + lam))[self.xfar]
        Rl = -Xb[self.k0_smin] - (3 + lam) / 2
        return np.concatenate([Rt, Ro, Rx, [Rl]])

    def jacobian(self, U):
        N = self.N
        Th, Om, X, lam = self.unpack(U)
        Ths, Thb = self.Ds @ Th, self.Dbeta @ Th
        Oms, Omb = self.Ds @ Om, self.Dbeta @ Om
        Xs, Xb = self.Ds @ X, self.Dbeta @ X
        D = (1 + lam) + Xb
        w = -(2 * X + Xs)
        dg = sp.diags
        # Θ equation
        Thu, Omu = self.Du @ Th, self.Du @ Om
        Jtt = (1 - lam) * sp.identity(N) + dg(D) @ self.Du + dg(w) @ self.Dbeta
        Jtx = dg(Thu) @ self.Dbeta + dg(Thb) @ (-2 * sp.identity(N) - self.Ds)
        Jtl = -Th + Thu
        # Ω equation
        Joo = sp.identity(N) + dg(D) @ self.Du + dg(w) @ self.Dbeta
        Jot = -(dg(self.eS * self.CB) @ self.Ds - dg(self.eS * self.SB) @ self.Dbeta)
        Jox = dg(Omu) @ self.Dbeta + dg(Omb) @ (-2 * sp.identity(N) - self.Ds)
        Jol = Omu
        # X equation
        Jxx = self.Dss + 4 * self.Ds + 4 * sp.identity(N) + self.Dbeta2
        Jxo = sp.identity(N)
        Z = sp.csr_matrix((N, N))
        # boundary rows
        C, Th0, Om0, Xs0 = self.local_data(lam)
        dC = -4.0 / (1 + lam) ** 2
        Jtt, Jtx = Jtt.tolil(), Jtx.tolil(); Joo, Jot, Jox = Joo.tolil(), Jot.tolil(), Jox.tolil()
        Jxx, Jxo = Jxx.tolil(), Jxo.tolil()
        Jtl = Jtl.copy(); Jol = Jol.copy(); Jxl = np.zeros(N)
        idx = np.where(self.inflow)[0]
        for i in idx:
            Jtt.rows[i] = [i]; Jtt.data[i] = [1.0]; Jtx.rows[i] = []; Jtx.data[i] = []
            Joo.rows[i] = [i]; Joo.data[i] = [1.0]; Jot.rows[i] = []; Jot.data[i] = []; Jox.rows[i] = []; Jox.data[i] = []
        Jtl[idx] = 0.0
        Jol[idx] = -(-dC * np.exp(self.s[0]) * np.cos(self.beta))
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
            kk = i % self.n
            Jxl[i] = -0.5 * dC * np.exp(self.s[0]) * np.cos(self.beta[kk]) * np.sin(self.beta[kk]) ** 2
        J = sp.bmat([[Jtt, Z, Jtx, sp.csr_matrix(Jtl[:, None])],
                     [Jot, Joo, Jox, sp.csr_matrix(Jol[:, None])],
                     [Z, Jxo, Jxx, sp.csr_matrix(Jxl[:, None])],
                     [None, None, sp.csr_matrix(-self.Dbeta[self.k0_smin].toarray()), sp.csr_matrix([[-0.5]])]],
                    format='csc')
        return J

    def newton(self, U, tol=1e-11, maxit=20, verbose=True):
        for it in range(maxit):
            R = self.residual(U)
            nr = np.abs(R).max()
            if verbose:
                print(f"   global Newton it {it}: |R|={nr:.3e}  λ={U[-1]:.12f}", flush=True)
            if nr < tol:
                return U, True
            J = self.jacobian(U)
            dU = splu(J, permc_spec='COLAMD').solve(-R)
            U = U + dU
        R = self.residual(U)
        return U, np.abs(R).max() < tol

    def initial_from_march(self, f, lam, Nb_src=32, hs_src=0.025):
        """interpolate a converged marching-solver state (bq_solver) onto this grid"""
        from bq_solver import BQ
        from bq_newton import full
        from numpy.polynomial import chebyshev as Cb
        B = BQ(lam, Nb=Nb_src, hs=hs_src)
        Y = np.load(f); Xm = full(B, Y)
        r = B.march(Xm / B.ea2[:, None])
        cp = np.where(B.cb > 0, B.cb, 0.0)
        Thm = r['Th'] * cp[None, :] ** r['m']; Omm = r['Omega']
        xo = np.cos(np.pi * np.arange(B.Nb + 1) / B.Nb); xn = np.cos(np.pi * np.arange(self.Nb + 1) / self.Nb)
        out = []
        for F in (Thm, Omm, Xm):
            Fb = Cb.chebval(xn, Cb.chebfit(xo, F.T, B.Nb))            # (Ns_src, Nb+1)
            Fi = np.array([np.interp(self.s, B.s, Fb[:, k]) for k in range(self.Nb + 1)]).T
            out.append(Fi.ravel())
        return np.concatenate(out + [[lam]])
