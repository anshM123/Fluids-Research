# Fluids-Research

Computational research on singularity formation in models of incompressible flow.

## Paper

**A resonance phase locates unstable self-similar singularities beyond an exponential precision wall.** This Article
was prepared for *Nature Computational Science*. Everything needed to read, check and submit it is in
[`ncs-submission/`](ncs-submission/):
- the manuscript and Supplementary Information;
- the figures;
- the code, with tests;
- the converged profiles and every logged output;
- the submission documents.

Start with [`ncs-submission/README.md`](ncs-submission/README.md).

## Repository layout

| Path | Content |
|---|---|
| `ncs-submission/` | the self-contained submission package (see above) |
| `programs/P03_boussinesq_ladder/` | working directory of the paper: self-similar blow-up hierarchies of 2D Boussinesq, the Hou–Luo model and incompressible porous-media flow, with all solvers, logs and notes. It also holds the pre-registration record [`PREDICTIONS_IPM.md`](programs/P03_boussinesq_ladder/PREDICTIONS_IPM.md), which is append-only. |
| `programs/P02_singularity_ladders/` | the Córdoba–Córdoba–Fontelos ladder: a finite hierarchy ending in a log-periodic square-root cusp. Includes code, method notes and a working draft. |
| `programs/P01_channel2d/`, `programs/P05_merger/` | exploratory probes (channel-flow chaos; vortex merger) |
| `ledger/` | idea registry (`IDEAS.md`), research log (`RLOG.md`), portfolio (`PORTFOLIO.md`) and program dossier (`DOSSIER.md`) |
| `literature/` | literature audits |
| `tools/` | shared solvers: 2D pseudo-spectral Navier–Stokes, a 2D channel solver and tensor-train analysis |

## Pre-registration

Predictions for the porous-media hierarchy were committed to this repository before the computations they concern.
Entries in `programs/P03_boussinesq_ladder/PREDICTIONS_IPM.md` are never edited; outcomes, deviations and corrections
are appended. The registration commits are listed in `ncs-submission/README.md`. The GitHub server time by which
each commit was public is in `ncs-submission/supplementary/SupplementaryData2_github_push_log.csv`. The history of
this repository is therefore part of the record: do not rewrite it with a rebase or force-push.
