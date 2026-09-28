# RESEARCH LOG (R-entries)

Format: project · hypothesis · computation · result · interpretation · cost · status · next step.
Costs are wall-clock on the 4-core container (CPU-h = cores × hours).

---
### R001 — Infrastructure / compute reconnaissance
- **Project:** all · **Hypothesis:** n/a
- **Computation:** hardware inventory; FFT micro-benchmarks (scipy.fft, 4 workers).
- **Result:** 4 × Xeon 2.8 GHz, 15 GB RAM, no GPU. 2D FFT pair: 256² 1.7 ms, 512² 4.3 ms, 1024² 16 ms, 2048² 58 ms. 3D: 64³ 2.8 ms, 128³ 20 ms, 256³ 306 ms.
- **Interpretation:** 1D/2D programs unconstrained; 3D limited to ≤128³ for long runs. Parallelism = 4 independent jobs.
- **Cost:** < 0.01 CPU-h · **Status:** DONE

### R002 — Shared solver validation
- **Project:** tools · **Hypothesis:** solvers reproduce exact/benchmark results.
- **Computation:** `tools/ps2d.py` (IF-RK4 pseudo-spectral 2D NS): Taylor–Green decay, Kolmogorov laminar state. `tools/channel2d.py` (Fourier × Legendre–Galerkin Shen basis, RK3/CN): Orr–Sommerfeld eigenvalue at Re=10⁴, α=1 (Orszag 1971) from the discrete operator, and growth rate of a 10⁻¹² perturbation in the nonlinear time-stepper.
- **Result:** TG error 1e−14; Kolmogorov laminar error 1.5e−9; OS eigenvalue c = 0.2375264888 + 0.0037396706 i (reference 0.23752649 + 0.00373967 i); time-stepper growth rate 0.003738–0.003740 (exact 0.00373967).
- **Interpretation:** both solvers verified to ≥ 8 digits on linear problems.
- **Cost:** 0.02 CPU-h · **Status:** VERIFIED
