"""Asymptotics of the CCF self-similar branch at its terminal cusp: fit p(δ), λ(δ) to
      y(δ) = y* + δ [A ln δ + B + C cos(ω ln δ) + D sin(ω ln δ)]
with ω = 2τ, τ tanh(πτ/2) = 1/2 (cusp-linearisation prediction), and with ω free.
Uses the post-regrid (best-resolved) values parsed from a pinned_delta log."""
import numpy as np, re, sys
from scipy.optimize import brentq, least_squares

tau = brentq(lambda t: t * np.tanh(np.pi * t / 2) - 0.5, 0.1, 2)
OMEGA = 2 * tau


def parse(log):
    rows, cur = [], None
    for line in open(log):
        m = re.match(r"δ=([-\d.e+]+) lam=([-\d.e+]+) p=([-\d.e+]+)", line)
        if m:
            if cur is not None:
                rows.append(cur)
            cur = [float(m.group(1)), float(m.group(2)), float(m.group(3))]
            continue
        m = re.search(r"regrid: .* lam=([-\d.e+]+) p=([-\d.e+]+)", line)
        if m and cur is not None:
            cur[1], cur[2] = float(m.group(1)), float(m.group(2))
    if cur is not None:
        rows.append(cur)
    return np.array(rows)


def model(par, d, omega=None):
    ys, A, B, C, D = par[:5]
    om = OMEGA if omega is None else par[5]
    L = np.log(d)
    return ys + d * (A * L + B + C * np.cos(om * L) + D * np.sin(om * L))


def fit(d, y, free_omega=False, use_log=True):
    p0 = [y[-1], 0.0, 0.0, 0.0, 0.0] + ([OMEGA] if free_omega else [])

    def res(par):
        q = list(par)
        if not use_log:
            q[1] = 0.0
        return (model(q, d, omega=None if not free_omega else True) - y) / (1e-12 + 1e-3 * d)
    r = least_squares(res, p0, method="lm", max_nfev=20000)
    return r.x, np.sqrt(np.mean(r.fun ** 2))


if __name__ == "__main__":
    logs = sys.argv[1:] if len(sys.argv) > 1 else ["pin_C.log"]
    for log in logs:
        R = parse(log)
        R = R[np.argsort(-R[:, 0])]
        print(f"== {log}: {len(R)} points, δ ∈ [{R[:, 0].min():.2e}, {R[:, 0].max():.2e}]")
        for r in R:
            print(f"   δ={r[0]:.4e}  λ={r[1]:.12f}  p={r[2]:.12f}")
        for dcut in (1e-2, 3e-3):
            sel = R[:, 0] <= dcut
            if sel.sum() < 7:
                continue
            d = R[sel, 0]
            for name, col in (("p", 2), ("λ", 1)):
                y = R[sel, col]
                for fo in (False, True):
                    for ul in (True, False):
                        par, rms = fit(d, y, free_omega=fo, use_log=ul)
                        print(f"   δ≤{dcut:g} {name}: {'free ω' if fo else 'ω=2τ '} {'with' if ul else 'no  '} ln-term: "
                              f"{name}*={par[0]:.10f} A={par[1] if ul else 0:.4f} B={par[2]:.4f} C={par[3]:.4f} D={par[4]:.4f}"
                              + (f" ω={par[5]:.4f} (pred {OMEGA:.4f})" if fo else "") + f"  rms={rms:.2e}")
