# Phase-based resolution of ill-conditioned singularity ladders in nonlinear fluid models

**Ansh Mishra** (independent researcher) and **Aryan Senthilkumar** (Georgia Institute of Technology).

Submission package for *Nature Computational Science* (Article): manuscript, Supplementary Information, code, data
and submission documents. **Start with [`submission/SUBMISSION_CHECKLIST.md`](submission/SUBMISSION_CHECKLIST.md).**

## Summary

Unstable self-similar blow-up profiles of incompressible flow come in hierarchies. Each profile is a zero of a
smoothness defect that shrinks about tenfold per profile, which creates an exponential precision wall for direct
computation. In incompressible porous-media flow (IPM) we show five things:
- **The error floor.** An audit of every error source finds an erratic floor of 10⁻⁹ in the defect. Two
  independent solvers agree to that floor.
- **A pre-registered test.** A holdout test of the seventh profile could not be decided by the defect.
- **The resonance phase.** A quadrature over O(1)-accurate profile data moves 10 to over 1,000 times less under the
  same numerical changes. It locates the seventh profile more than ten times more precisely, and its imaginary part
  predicts the wall.
- **A tracking failure.** A failure that passed residual checks was exposed by a geometric invariant and repaired
  with an exact scaling.
- **A benchmark.** Three deeper profiles are registered for future tests.

## Contents

| Path | Content |
|---|---|
| `manuscript/` | `main.tex` and the compiled `main.pdf`: the main text, Methods, six figures and Extended Data Tables 1–2 |
| `supplementary/` | `supplementary.tex` and `.pdf`: Supplementary Notes 1–7, Figs 1–3 and Tables 1–12. Also `tables/`, generated from the data, and Supplementary Data 1 (the pre-registration record) and 2 (the GitHub push log). |
| `figures/` | Figs 1–6 and Supplementary Figs 1–3, as vector PDF and 600-dpi PNG, 183 mm wide |
| `code/` | solver library, production runs, analysis, figure script and tests; see [`code/README.md`](code/README.md) |
| `data/` | ten converged profiles (31 MB) and every logged output (0.7 MB); see [`data/README.md`](data/README.md) |
| `submission/` | cover letter (`.tex` and `.pdf`), checklist, referee template, and notes for the Reporting Summary and the Code and Software Submission Checklist |

## Quick start

```bash
python3 -m pip install -r requirements.txt
python3 -m pytest code/tests -q                    # 5 smoke tests, 15 s
python3 code/figures/make_figures.py               # all figures from data/outputs, 1 min
python3 code/analysis/make_si_tables.py            # all Supplementary Tables from data/outputs
cd manuscript && latexmk -pdf main.tex             # the manuscript (TeX Live 2023 or later)
```

Every number in the figures and Supplementary Tables is read from `data/outputs` by these two scripts; nothing is
transcribed by hand. `code/README.md` describes, with expected values and run times:
- the analyses that recompute the seventh-profile determination and the registered predictions from the logs,
  byte for byte;
- the recomputation of the phase, the independent global solver and profile location from the stored profiles;
- the full computational campaign (about 13 core-hours on a four-core machine, no GPU).

## Pre-registration

All predictions were committed to this public repository before the computations they concern, in
[`programs/P03_boussinesq_ladder/PREDICTIONS_IPM.md`](../programs/P03_boussinesq_ladder/PREDICTIONS_IPM.md). The
record is append-only: outcomes, deviations and corrections are added, never edited. The registration commits are:

| Commit | Content |
|---|---|
| `4dc82af` | stage 1 |
| `c922751`, `273ae28` | stage 3 |
| `f0b3f06` | the seventh-profile decision rule |
| `251bf1e` | stage 4, the open benchmark for λ₈–λ₁₀ |

Supplementary Data 2 lists the GitHub server time by which each commit was public.

## Licence

- Code: MIT ([`LICENSE`](LICENSE)).
- Data in `data/`: CC BY 4.0 ([`data/LICENSE-DATA.md`](data/LICENSE-DATA.md)).
- The manuscript, Supplementary Information and figures are under journal submission. Their licence follows the
  publication agreement.
