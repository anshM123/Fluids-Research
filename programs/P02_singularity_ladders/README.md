# P02 — Self-similar blow-up ladder of the CCF equation: reproduction guide

Everything below runs on a 4-core CPU with 15 GB RAM (NumPy/SciPy; numba only for `stability_mapped.py`).
Set `OMP_NUM_THREADS` to 1–2 per job. Times are wall-clock on that machine.

## Core solvers
| File | Purpose |
|---|---|
| `ccf_nk.py` | uniform FFT log-grid Newton–Krylov (fixed λ, p-constrained, arclength) |
| `ccf_mapped.py` | tanh-mapped adaptive grid, alternating-point Hilbert rule, dense/Krylov |
| `ccf_sinhgrid.py` | geometric grid for thin sonic layers (accurate local coordinates) |
| `pinned_delta.py` | layer-pinned continuation in the sonic depth δ (main tool near the cusp) |
| `pinned_dense.py` | same, with fixed δ ratio and re-gridding every step (dense sampling) |
| `fccf_nk.py` | fractional family θ_t + (HΛ^sθ)θ_x = 0 (exact Mellin multiplier) |

## Key numbers and how to regenerate them
| Result | Command | Time |
|---|---|---|
| λ₀ = 1.180777662899, λ₂ = 0.471324227767 (12 digits) | `python3 mapped_l2b.py` (λ₂); `ladder_scan.py` (λ₀, λ₁) | 20 s – few min |
| instability spectra λ₀–λ₂ | `python3 dense_spectrum.py l0|l1|l2`, `stab_scan3.py`, `linear_evolution.py` | ~1 min each |
| branch map p(λ), λ ∈ [0.4625, 2] | `python3 branch_map.py` | ~30 min |
| large-λ end (p ↓ 1, no crossing) | `python3 large_lambda_arc.py 0.9 500 65536` | ~15 min |
| arc λ₂ → cusp (pinned) | `python3 pinned_delta.py pin_start.npy C 24 0.03` | ~15 min (2 threads) |
| dense sampling of the arc | `python3 pinned_dense.py pin_start.npy E 0.85 8e-6 24 0.03` | ~40 min |
| gap λ ∈ [0.4564, 0.4644] | `python3 pinned_dense.py pin_start.npy G 1.08 0.0215 24 0.03` | ~7 min |
| grid-refinement at fixed δ | `python3 verify_points.py pdense_states/E_005.npy ...` | ~10 min / state |
| asymptotic fits (p*, λ*, ω) | `python3 cusp_fit.py pin_C.log` | seconds |
| cusp amplitude law B² = 2λΘ_s | `python3 cusp_check.py pin_C_final_state.npy` | ~1 min |
| universal inner layer (f''(0) = 0.08900) | `python3 inner_layer.py` | ~2 min |
| universal exponents (τ = 0.6494243) and fractional τ(s) | `python3 cusp_exponents_s.py` | seconds |
| operator accuracy on extreme grids | `python3 hilbert_extreme_test.py` | ~3 min |
| Fig. 1 / Fig. 2 | `python3 fig_ladder.py`; `python3 fig_profiles2.py pin_C_final_state.npy "near-terminal (δ = 1.3e−6)"` | seconds |

Reference data for fits are the **post-regrid** values in `pin_C.log` (fully resolved) and all values in
`pdense_E.npy`/`pdense_G.npy` (re-gridded at every step). `pin_*_oldgrid.log` and `pin_A_final.log` were produced
before the coordinate-precision fix (R019) and are kept only for the record (valid for δ ≳ 3e−4).
