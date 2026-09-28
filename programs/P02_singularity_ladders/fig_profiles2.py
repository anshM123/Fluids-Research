"""Fig. 2: the three smooth CCF self-similar profiles and a near-terminal (cusp) profile.
Θ normalised by its far-field amplitude C (Θ ~ Cξ^β); local slope d lnΘ/d lnξ = λ/den shows the sonic layer."""
import numpy as np, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

c = 0.7
prof = []
for f, lab in (("ladder_F16k_up_cross2.npy", "λ₀ (stable)"), ("ladder_F16k_up_cross1.npy", "λ₁ (1 unstable)")):
    d = np.load(f)
    eta = -30.0 + (150.0 / 16384) * np.arange(16384)
    prof.append((lab, float(d[0]), eta, d[1:]))
d = np.load("mapped_l2_hs0.02_eps0.025.npy")
prof.append(("λ₂ (2 unstable)", 0.471324227767, d[0], d[1]))
tail = sys.argv[1] if len(sys.argv) > 1 else "pin_C_last.npy"
d = np.load(tail)
prof.append((sys.argv[2] if len(sys.argv) > 2 else "near-terminal", float(d[2][0]), d[0], d[1]))

fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
for lab, lam, eta, phi in prof:
    beta = lam / (1 + lam)
    lnT = phi + c * eta
    iR = int(np.argmin(np.abs(eta - 50.0)))
    lnC = lnT[iR] - beta * eta[iR]
    m = (eta > -5) & (eta < 5)
    xi = np.exp(eta[m])
    ax[0].loglog(xi, np.exp(lnT[m] - lnC), label=f"{lab}: λ = {lam:.6f}")
    slope = np.gradient(lnT, eta)
    ax[1].semilogx(xi, slope[m], label=lab)
    ax[2].semilogx(xi, lam / slope[m], label=lab)
ax[0].set_xlabel("ξ"); ax[0].set_ylabel("Θ / C"); ax[0].legend(fontsize=8)
ax[0].set_title("profiles (far-field amplitude normalised)")
ax[1].set_xlabel("ξ"); ax[1].set_ylabel("d ln Θ / d ln ξ = λ / den"); ax[1].set_ylim(0, 12)
ax[1].set_title("local slope: sonic layer grows into a cusp")
ax[2].set_xlabel("ξ"); ax[2].set_ylabel("den = 1 + λ + HΘ/ξ"); ax[2].set_ylim(0, 1.5)
ax[2].set_title("sonic factor")
plt.tight_layout()
plt.savefig("fig_profiles2.png", dpi=140)
print("saved fig_profiles2.png")
