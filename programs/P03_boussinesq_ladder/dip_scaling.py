"""Scaling of the dip in front of the front as ε → 0: location s_d, depth D̂_min, half-width w (where D̂ < (1+D̂_min)/2)
and the front location s_f (first D̂ = 3 after the dip). Tests whether the dip is an inner layer of width ∝ ε."""
import numpy as np, warnings
warnings.filterwarnings('ignore')
from bq_wkb2d import boundary_params
files = ["Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.399096.npy", "Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.252349.npy",
         "Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.184253.npy"] + \
        [f"Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam{l}.npy" for l in ("1.144986", "1.119474", "1.101582", "1.088338")]
rows = []
for f in files:
    lam, m, eps, s, D, mu, G, Thy, Omb = boundary_params(f)
    Dh = D / eps
    sel = (s > -3) & (s < 1)
    ss, dd = s[sel], Dh[sel]
    kd = np.argmin(dd); sd, dmin = ss[kd], dd[kd]
    half = 0.5 * (1 + dmin)
    left = kd - np.argmax(dd[kd::-1] > half); right = kd + np.argmax(dd[kd:] > half)
    w = ss[right] - ss[left]
    kf = kd + np.argmax(dd[kd:] > 3); sf = ss[kf]
    # integrated 'dip deficit' ∫(1 − D̂)_+ ds over the dip
    deficit = np.trapezoid(np.clip(1 - dd[left:right + 1], 0, None), ss[left:right + 1])
    rows.append((eps, sd, dmin, w, sf, deficit))
    print(f"ε={eps:.5f}: dip at s={sd:.4f} (x={np.exp(sd):.4f}), D̂_min={dmin:.4f}, half-width={w:.4f}, front(D̂=3) s={sf:.4f}, gap s_f−s_d={sf-sd:.4f}, deficit={deficit:.4f}", flush=True)
R = np.array(rows)
for j, name in ((2, 'D̂_min'), (3, 'half-width'), (5, 'deficit')):
    p = np.polyfit(np.log(R[3:, 0]), np.log(R[3:, j]), 1)
    print(f"{name} ∝ ε^{p[0]:.3f} (fit on the four deepest)")
p = np.polyfit(np.log(R[3:, 0]), np.log(R[3:, 4] - R[3:, 1]), 1); print(f"gap s_f − s_d ∝ ε^{p[0]:.3f}")
