"""Figure: linear stability of the eight smooth profiles.
(a) number of real eigenvalues of T_μ above 1 minus the trivial mode, versus μ (march-based method, bq_stability.py);
(b) the unstable eigenvalues μ_k of each profile from the independent global-discretisation eigen-solver (glob_stab.py),
    with the brackets of method (a)."""
import numpy as np, glob, re, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
MU_ART = 0.045          # below this the march-based count is affected by the origin truncation (see README)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
a, b = ax
cols = plt.cm.viridis(np.linspace(0, 0.95, 8))
for n in range(8):
    f = f"stab_l{n}.npy"
    if not os.path.exists(f):
        continue
    R = np.load(f)
    mu, N = R[:, 0], R[:, 1] - (R[:, 0] < 1.0)     # remove the trivial time-translation crossing at μ = 1
    sel = mu >= MU_ART
    a.step(mu[sel], N[sel] + 0.03 * n, where='post', color=cols[n], label=f"λ_{n}")
    # brackets for panel (b)
    for i in range(len(mu) - 1):
        if R[i + 1, 1] > R[i, 1] and R[i + 1, 0] >= MU_ART and R[i, 0] < 0.99:
            b.plot([n, n], [R[i + 1, 0], R[i, 0]], color='0.7', lw=6, solid_capstyle='butt', zorder=1)
a.axvspan(0, MU_ART, color='0.9'); a.set_xlim(0, 1.05)
a.set_xlabel('μ (growth rate of the perturbation, e^{μτ})'); a.set_ylabel('number of eigenvalues > μ (trivial mode removed)')
a.set_title('(a) march-based eigen-condition ν(μ) = 1: count of unstable eigenvalues', fontsize=9)
a.legend(fontsize=7, ncol=2)
for n in range(8):
    for f in (f"gstab_g{n}s20.log", f"gstab_g{n}.log"):
        if not os.path.exists(f):
            continue
        txt = open(f).read(); m = re.search(r"RESULT.*?: (.*)", txt)
        if not m:
            continue
        ev = [complex(t.replace('i', 'j')) for t in m.group(1).split(', ') if t.strip()]
        ev = [e for e in ev if abs(e.imag) < 1e-6 and 0.06 < e.real < 0.99]
        b.plot([n] * len(ev), [e.real for e in ev], 'o', color='tab:red', ms=5, zorder=3)
        break
b.plot(range(8), [1.0] * 8, 'k_', ms=12, label='trivial (time translation), μ = 1')
b.set_xlabel('profile n (continuation order)'); b.set_ylabel('unstable eigenvalue μ_k')
b.set_title('(b) unstable eigenvalues: global eigen-solver (red) and brackets of (a) (grey)', fontsize=9)
b.set_ylim(0, 1.1); b.legend(fontsize=7, loc='upper left')
plt.tight_layout(); plt.savefig('fig_bq_stability.png', dpi=150)
