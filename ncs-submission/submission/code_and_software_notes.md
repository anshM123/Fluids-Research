# Code and Software Submission Checklist: answers, and the Code Ocean capsule

*Nature Computational Science* peer-reviews custom code that is central to a paper. When the manuscript is sent for
review, the editors ask for:
- the Code and Software Submission Checklist;
- access to the code for the referees, normally through a private Code Ocean capsule.

## Checklist answers

| Item | Answer |
|---|---|
| Code location | https://github.com/anshM123/Fluids-Research, directory `ncs-submission/code`; Zenodo [DOI]; Code Ocean capsule [link] |
| Licence | MIT (OSI-approved), `ncs-submission/LICENSE`; data under CC BY 4.0 |
| System requirements | Linux x86-64 (tested on Ubuntu 24.04), Python 3.11, NumPy 2.4.6, SciPy 1.17.1, Numba 0.67.0, Matplotlib 3.11.2 (`requirements.txt`). No GPU or non-standard hardware. |
| Installation guide | `python3 -m pip install -r requirements.txt` (about 1 min on a normal connection) |
| Demo | `python3 -m pytest code/tests -q` (15 s). Then `python3 code/figures/make_figures.py` and `python3 code/analysis/make_si_tables.py` regenerate all figures and Supplementary Tables from `data/outputs` (1 min). |
| Expected output | 5 tests pass; `figures/` and `supplementary/tables/` are rewritten, and the tables are identical to the committed ones |
| Instructions for use | `code/README.md`: smoke tests, figures, analyses from logs, recomputation from the stored profiles (the resonance phase, the independent global solver and profile location, with expected values and run times), and the full campaign |
| Reproduction of the main claims | `code/README.md`, sections 3–4: the seventh-profile analysis, the registered predictions and the endpoint analysis are reproduced from the logs byte for byte; the phase, the global solver and profile location are reproduced from the stored profiles to all printed digits |

## Code Ocean capsule (recipe)

1. Create a capsule from the GitHub repository, or upload `ncs-submission/` (code, data, figures, supplementary).
2. Environment: the Python 3.11 base image; install the packages in `requirements.txt` with pip.
3. Put the contents of `ncs-submission/` in `/code`. Put `data/` in `/data`, or leave it in place: the scripts find
   it relative to `code/`.
4. Make `/code/run` an executable shell script:

   ```bash
   #!/usr/bin/env bash
   set -e
   cd /code
   python3 -m pytest code/tests -q
   NCS_FIGS=/results/figures   python3 code/figures/make_figures.py
   NCS_TABLES=/results/tables  python3 code/analysis/make_si_tables.py
   cd data/outputs
   for s in ipm_lambda7_final ipm_stage4_predict ipm_endpoint_geometry; do
       python3 ../../code/analysis/$s.py > /results/$s.out
   done
   cd /results
   python3 /code/code/analysis/ipm_phase3_table.py \
       /code/data/profiles/ipm_profile_n3_lam0.2415663353_Nb32_hs0.0125.npy > /results/phase_lambda3.out
   ```

5. A reproducible run takes about 3 minutes. The capsule stays private, with anonymous access for referees, until
   publication.
