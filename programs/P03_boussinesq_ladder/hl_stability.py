"""Linear stability of self-similar Hou–Luo profiles (1D boundary model), by the same eigen-condition as for 2D
(bq_stability.py): perturbations e^{μτ}, τ = −ln(1−t), obey in η = ln ξ
    (μ+1−λ) θ' + D θ'_η + q' Θ̄_η = 0,
    (μ+1) ω' + D ω'_η + q' Ω̄_η = e^{−η} θ'_η,        q' = u'/ξ = Mellin multiplier of ω' (as in hl_solver.HL.Q).
For a given q', (θ', ω') are marched from zero data at η_start (regular perturbations) with the 2-stage Gauss–Legendre
scheme of the base solver; the new q' defines the linear map T_μ. μ is an eigenvalue iff T_μ has an eigenvalue 1.
Complex μ is supported (complex march, real Mellin map applied to real and imaginary parts)."""
import numpy as np
from numba import njit
from scipy.sparse.linalg import LinearOperator, eigs
from hl_solver import HL


@njit(cache=True)
def _lin_gl2(h, i0, mu, lam, eta, D1, D2, Tt1, Tt2, Ot1, Ot2, q1, q2):
    """march θ', ω' (complex) on [i0, N): θ'_η = −a θ' + r, ω'_η = −b ω' + s with a = (μ+1−λ)/D, r = −q'Θ̄_η/D,
    b = (μ+1)/D, s = (e^{−η}θ'_η − q'Ω̄_η)/D, coefficients at the two Gauss points of each interval"""
    N = len(D1) + 1
    th = np.zeros(N, dtype=np.complex128); om = np.zeros(N, dtype=np.complex128)
    A11, A12, A21, A22 = 0.25, 0.25 - np.sqrt(3.0) / 6, 0.25 + np.sqrt(3.0) / 6, 0.25
    c1, c2 = 0.5 - np.sqrt(3.0) / 6, 0.5 + np.sqrt(3.0) / 6
    for j in range(i0, N - 1):
        a1 = (mu + 1.0 - lam) / D1[j]; a2 = (mu + 1.0 - lam) / D2[j]
        r1 = -q1[j] * Tt1[j] / D1[j]; r2 = -q2[j] * Tt2[j] / D2[j]
        y = th[j]
        m11 = 1 + h * a1 * A11; m12 = h * a1 * A12; m21 = h * a2 * A21; m22 = 1 + h * a2 * A22
        f1 = r1 - a1 * y; f2 = r2 - a2 * y
        det = m11 * m22 - m12 * m21
        K1 = (f1 * m22 - m12 * f2) / det; K2 = (m11 * f2 - m21 * f1) / det
        th[j + 1] = y + 0.5 * h * (K1 + K2)
        e1 = np.exp(-(eta[j] + c1 * h)); e2 = np.exp(-(eta[j] + c2 * h))
        b1 = (mu + 1.0) / D1[j]; b2 = (mu + 1.0) / D2[j]
        s1 = (e1 * K1 - q1[j] * Ot1[j]) / D1[j]; s2 = (e2 * K2 - q2[j] * Ot2[j]) / D2[j]
        z = om[j]
        n11 = 1 + h * b1 * A11; n12 = h * b1 * A12; n21 = h * b2 * A21; n22 = 1 + h * b2 * A22
        g1 = s1 - b1 * z; g2 = s2 - b2 * z
        det = n11 * n22 - n12 * n21
        L1 = (g1 * n22 - n12 * g2) / det; L2 = (n11 * g2 - n21 * g1) / det
        om[j + 1] = z + 0.5 * h * (L1 + L2)
    return th, om


class HLStab:
    def __init__(self, lam, q, N, L1=40.0, L2=160.0, c=-0.25, eta_start=-30.0, sign=1.0):
        self.S = S = HL(lam, L1=L1, L2=L2, N=N, c=c, eta_start=eta_start, sign=sign)
        self.lam = lam
        q = S.full(np.asarray(q, float))
        r = S.march(q)
        self.m, self.A, self.eps = r['m'], r['A'], r['eps']
        self.q = q
        h = S.h; eta = S.eta
        Th, Om = r['Th'], r['Om']
        lnT = np.log(Th / sign)
        self.D1 = (1 + lam) + S.st.offset(q, S.c1); self.D2 = (1 + lam) + S.st.offset(q, S.c2)
        T1 = sign * np.exp(S.st.offset(lnT, S.c1)); T2 = sign * np.exp(S.st.offset(lnT, S.c2))
        self.Tt1 = (lam - 1) * T1 / self.D1; self.Tt2 = (lam - 1) * T2 / self.D2          # Θ̄_η at Gauss points
        O1 = S.st.offset(Om, S.c1); O2 = S.st.offset(Om, S.c2)
        e1 = np.exp(-(eta[:-1] + S.c1 * h)); e2 = np.exp(-(eta[:-1] + S.c2 * h))
        self.Ot1 = (e1 * self.Tt1 - O1) / self.D1; self.Ot2 = (e2 * self.Tt2 - O2) / self.D2   # Ω̄_η
        self.i0 = S.i0; self.n = N - S.i0
        self.Th, self.Om = Th, Om

    def _mellin(self, om):
        S = self.S
        return S.Ec * np.real(np.fft.ifft(S.M * np.fft.fft(S.Emc * om)))

    def T(self, v, mu):
        S = self.S; i0 = self.i0
        qp = np.empty(S.N, dtype=complex); qp[i0:] = v; qp[:i0] = v[0]
        q1 = S.st.offset(qp.real, S.c1) + 1j * S.st.offset(qp.imag, S.c1)
        q2 = S.st.offset(qp.real, S.c2) + 1j * S.st.offset(qp.imag, S.c2)
        th, om = _lin_gl2(S.h, i0, complex(mu), self.lam, S.eta, self.D1, self.D2, self.Tt1, self.Tt2,
                          self.Ot1, self.Ot2, q1, q2)
        qn = self._mellin(om.real) + 1j * self._mellin(om.imag)
        return qn[i0:]

    def spectrum(self, mu, k=10, tol=1e-10, v0=None):
        dt = complex
        op = LinearOperator((self.n, self.n), matvec=lambda v: self.T(v, mu), dtype=dt)
        vals, vecs = eigs(op, k=k, which='LM', tol=tol, v0=v0, maxiter=5000)
        o = np.argsort(-np.abs(vals))
        return vals[o], vecs[:, o]

    def time_translation(self):
        """q' of the μ = 1 mode: q' = q + (1+λ) q_η"""
        q = self.q; S = self.S
        qe = np.gradient(q, S.h)
        return (q + (1 + self.lam) * qe)[self.i0:]
