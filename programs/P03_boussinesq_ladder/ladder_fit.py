"""The instability ladder: μ_{n,k}/(λ_n − 1) for fixed k (k = 0: lowest unstable eigenvalue) versus δ_n = λ_n − 1,
extrapolated to δ → 0, for 2D Boussinesq (method-2 eigenvalues; λ₇'s lowest from method 1) and Hou–Luo (real
eigenvalues from the parity-flip refinement)."""
import numpy as np
bq_lam = [1.9205593, 1.3990961, 1.2523487, 1.1842533, 1.1449857, 1.1194738, 1.1015817, 1.0883384]
bq = {1: [0.37379], 2: [0.55419, 0.22048], 3: [0.63437, 0.37541, 0.15492], 4: [0.68004, 0.45973, 0.28629, 0.11908],
      5: [0.70985, 0.51143, 0.36869, 0.23121, 0.09654], 6: [0.73139, 0.54481, 0.42577, 0.30749, 0.19415, 0.08154],
      7: [0.74632, 0.56676, 0.46713, 0.36192, 0.26371, 0.16678, 0.069854]}
hl_lam = [1.9987042, 1.44767467, 1.28676473, 1.21091642, 1.16668059, 1.13772677, 1.11731019, 1.10214553, 1.09044036]
hl = {1: [0.373847], 2: [0.573937, 0.226768], 3: [0.662131, 0.401562, 0.159959],
      4: [0.711427, 0.496887, 0.308705, 0.123332], 5: [0.742501, 0.553904, 0.402639, 0.250691, 0.100103],
      6: [0.763706, 0.587763, 0.468582, 0.337549, 0.210951, 0.084157],
      7: [0.779113, 0.602946, 0.522013, 0.39978, 0.291248, 0.182021, 0.072544],
      8: [0.790948, 0.443727, 0.352774, 0.255826, 0.160042, 0.063722]}      # + complex pair 0.587182 ± 0.022206i
def analyse(name, lam, ev, kmax=4):
    print(f"== {name}")
    for k in range(kmax):
        ns = [n for n in sorted(ev) if len(ev[n]) > k and n >= 2]
        d = np.array([lam[n] - 1 for n in ns]); y = np.array([sorted(ev[n])[k] / (lam[n] - 1) for n in ns])
        out = [f"k={k}: " + " ".join(f"{v:.4f}" for v in y)]
        for use, deg in ((3, 1), (4, 2)):
            if len(d) >= use:
                p = np.polyfit(d[-use:], y[-use:], deg); out.append(f"lim(deg{deg},last{use})={np.polyval(p, 0):.4f}")
        print("   " + "  ".join(out))
    for k in range(kmax - 1):
        ns = [n for n in sorted(ev) if len(ev[n]) > k + 1 and n >= 2]
        d = np.array([lam[n] - 1 for n in ns])
        y = np.array([(sorted(ev[n])[k + 1] - sorted(ev[n])[k]) / (lam[n] - 1) for n in ns])
        out = [f"Δ(k={k}→{k+1})/(λ−1): " + " ".join(f"{v:.4f}" for v in y)]
        for use, deg in ((3, 1), (4, 2)):
            if len(d) >= use:
                p = np.polyfit(d[-use:], y[-use:], deg); out.append(f"lim(deg{deg},last{use})={np.polyval(p, 0):.4f}")
        print("   " + "  ".join(out))
analyse("2D Boussinesq", bq_lam, bq)
analyse("Hou–Luo", hl_lam, hl)
