"""INDEPENDENT stability computation on the global discretisation of bq_global.py (different from the T_μ march of
bq_stability.py in every ingredient). Perturbations e^{μτ} of a converged global profile U (λ fixed) satisfy the
generalised eigenproblem
        J v + μ M v = 0,    M = diag(1 on the Θ, Ω transport rows, 0 on inflow / Biot–Savart / boundary rows),
where J is the exact sparse Jacobian of the steady residual (λ column and λ equation removed). Eigenvalues near a
(complex) shift σ are found by Arnoldi on (−J − σM)⁻¹M, whose eigenvalues are 1/(μ − σ).
usage: python3 glob_stab.py U_FILE SMIN SMAX HS NB TAG [shifts as re:im,...]"""
import numpy as np, scipy.sparse as sp, sys, time
from scipy.sparse.linalg import splu, LinearOperator, eigs
from bq_global import BQGlobal
fU = sys.argv[1]; smin, smax, hs, Nb = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
tag = sys.argv[6]
KEIG = int(sys.argv[8]) if len(sys.argv) > 8 else 12
BC = sys.argv[9] if len(sys.argv) > 9 else 'dirichlet'
shifts = ([complex(*map(float, p.split(':'))) for p in sys.argv[7].split(',')] if len(sys.argv) > 7 else
          [complex(x, y) for y in (0.0, 0.4) for x in (0.1, 0.3, 0.5, 0.7, 0.9, 1.1)])
G = BQGlobal(s_min=smin, s_max=smax, hs=hs, Nb=Nb)
U = np.load(fU); lam = U[-1]; N = G.N
log = open(f"gstab_{tag}.log", "w")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# global-discretisation stability {tag} (inflow BC {BC}): λ={lam:.10f} |R|={np.abs(G.residual(U)).max():.1e} "
    f"s∈[{smin},{smax}] hs={hs} Nb={Nb}")
J = G.jacobian(U)[:3 * N, :3 * N].tocsr()
if BC == 'robin':
    # smooth-perturbation inflow condition θ'_s = 2θ' (θ' ∝ r² at the stagnation point) instead of θ' = 0: it removes
    # the near-smooth homogeneous solutions r^{2−μ/ε} that a Dirichlet condition at finite s_min only damps
    inf = G.inflow.astype(float)
    Zb = sp.csr_matrix((N, N))
    Rt = sp.hstack([sp.diags(inf) @ (G.Ds - 2 * sp.identity(N)), Zb, Zb])
    keep = np.concatenate([1 - inf, np.ones(2 * N)])
    J = (sp.diags(keep) @ J + sp.vstack([Rt, sp.csr_matrix((2 * N, 3 * N))])).tocsc()
else:
    J = J.tocsc()
mdiag = np.concatenate([np.where(G.inflow, 0.0, 1.0), np.where(G.inflow, 0.0, 1.0), np.zeros(N)])
M = sp.diags(mdiag).tocsc()
found = []
for sig in shifts:
    t = time.time()
    A = (-J - (sig if sig.imag != 0 else sig.real) * M).tocsc()
    lu = splu(A, permc_spec='COLAMD')
    dt = complex if sig.imag != 0 else float
    op = LinearOperator((3 * N, 3 * N), matvec=lambda x: lu.solve(M @ x), dtype=dt)
    th, V = eigs(op, k=KEIG, which='LM', tol=1e-10, maxiter=5000)
    mu = sig + 1 / th
    # residual check of each pair
    good = []
    for i in range(len(mu)):
        v = V[:, i]; r = np.linalg.norm(J @ v + mu[i] * (M @ v)) / max(np.linalg.norm(M @ v), 1e-300)
        good.append((mu[i], r))
    good.sort(key=lambda p: -p[0].real)
    out(f"σ={sig.real:.2f}{sig.imag:+.2f}i: " + ", ".join(f"{m.real:.5f}{m.imag:+.5f}i(r={r:.0e})" for m, r in good)
        + f"  ({time.time()-t:.0f}s)")
    found += [m for m, r in good if r < 1e-6]
# unique eigenvalues with Re μ > −0.05
uniq = []
for m in sorted(found, key=lambda c: -c.real):
    if m.real > -0.05 and all(abs(m - u) > 1e-4 for u in uniq):
        uniq.append(m)
out("RESULT " + tag + ": eigenvalues with Re μ > −0.05: " + ", ".join(f"{u.real:.5f}{u.imag:+.5f}i" for u in uniq))
