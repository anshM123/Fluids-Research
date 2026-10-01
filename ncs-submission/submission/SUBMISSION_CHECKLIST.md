# Submission checklist: *Nature Computational Science*, Article

Everything in the package is filled in and the PDFs are built. What remains is the checks in section 3 and the upload.

## 1. What to upload

| Upload as | File |
|---|---|
| Manuscript (PDF; includes Methods, the 6 figures and Extended Data Tables 1–2) | `manuscript/main.pdf` |
| Cover letter | `submission/cover_letter.pdf` |
| Supplementary Information | `supplementary/supplementary.pdf` |
| Supplementary Data 1 | `supplementary/SupplementaryData1_preregistration_record.md`, a copy of the pre-registration record (plain text) |
| Supplementary Data 2 | `supplementary/SupplementaryData2_github_push_log.csv`, the GitHub server times of every commit |

The journal accepts initial submissions as PDF, Word or LaTeX. The LaTeX source is only needed at acceptance.

## 2. What is filled in

| Item | Content |
|---|---|
| Authors | Ansh Mishra¹ and Aryan Senthilkumar² |
| Affiliations | ¹ Independent researcher; ² Georgia Institute of Technology, Atlanta, GA, USA |
| Corresponding author | Ansh Mishra |
| Author contributions | A.M. carried out the mathematics and physics. A.S. carried out the aerodynamics and fluid mechanics. |
| Use of AI tools (Methods) | "AI tools were used continually throughout the development of this project." |
| Acknowledgements | none (the section is omitted) |
| Competing interests | none declared |
| Data and code availability | the public GitHub repository; no archive DOI |
| Licences | code MIT, data CC BY 4.0; copyright held by both authors |
| Cover letter | from and signed by both authors; dated automatically when compiled |

If you change anything, rebuild the PDFs. This needs TeX Live 2023 or later with `latexmk`:

```bash
cd manuscript    && latexmk -pdf main.tex          && cd ..
cd supplementary && latexmk -pdf supplementary.tex && cd ..
cd submission    && latexmk -pdf cover_letter.tex  && cd ..
```

## 3. Checks before submission

- [ ] **References.** Confirm authors, titles and identifiers of the two most recent preprints: Wang et al.,
  arXiv:2509.14185 (ref. 6), and Wang, Léger, Lai & Buckmaster, arXiv:2511.22819 (ref. 8). Replace any preprint
  that has since been published with its journal version: Chen & Hou I and II (refs 2, 3) and Buckmaster et al.
  (ref. 10).
- [ ] **Overlap with ref. 8.** If Wang et al. (arXiv:2511.22819) report IPM profiles, compare their λ values with
  the λ₁–λ₆ of this paper (`data/README.md`) and add one sentence to the main text on agreement or difference.
  The comparison was not possible from the environment in which the package was prepared.
- [ ] **Public record.** The paper states that predictions were committed to a public repository before the
  computations. GitHub's public event log covers the repository from its creation on 28 September 2026, which
  indicates it was public throughout. Confirm that it was never private. If it was, change "public repository" to
  "version-controlled repository, now public" in the main text, the Methods and Supplementary Note 3.
- [ ] **Push log.** `SupplementaryData2_github_push_log.csv` was retrieved from GitHub's public events interface,
  which keeps 90 days. Keep this file. It cannot be regenerated after late December 2026.
- [ ] **Read-through.** Read the full manuscript once against `manuscript/main.pdf`.

## 4. In the submission system

- **Article type:** Article.
- **Title:** *A resonance phase locates unstable self-similar singularities beyond an exponential precision wall*
  (12 words).
- **Abstract:** paste from `manuscript/main.tex` (148 words).
- **Authors:** as in section 2. Enter the corresponding author's e-mail address in the form. ORCIDs are optional.
- **Preprint:** the system offers to post the submission as a preprint on Research Square, which gives it a DOI.
  Opt in if you want it. The cover letter says the manuscript has not been posted as a preprint, which is true at
  submission.
- **Referees (optional):** suggested and excluded referees; see `suggested_reviewers.md`.
- **Competing interests:** none.
- **Funding:** answer in the form. The manuscript has no acknowledgements section.

## 5. If the paper goes out to review

- **Code review.** *Nature Computational Science* peer-reviews custom code that is central to a paper. The editors
  will ask for the Code and Software Submission Checklist and a private Code Ocean capsule for the referees;
  `code_and_software_notes.md` has the answers and the capsule recipe.
- **Default branch (optional).** The repository's default and only branch is
  `claude/fluid-dynamics-research-qo72sn`. To present a cleaner URL, rename it to `main` (GitHub: Settings →
  Branches → rename). GitHub redirects the old name, and the commit hashes cited in the paper are unchanged. Never
  rebase or force-push: the history is part of the pre-registration record.

## 6. Compliance with the Article format

The limits are from the journal's guide to authors (Content types; Formatting your initial submission). Check them
against the current guide at https://www.nature.com/natcomputsci/submission-guidelines, since they change.

| Item | Limit | This manuscript |
|---|---|---|
| Main text (excluding abstract, Methods, references, legends) | about 3,500 words | 3,290 |
| Abstract | 150 words, no references | 148 |
| Display items (figures and tables) | 6 | 6 figures |
| Extended Data items | 10 | 2 tables |
| References | about 50 | 37 |
| Methods | no strict limit | about 900 words |
| Figures | sized for 89 mm or 183 mm, sans-serif labels, bold lowercase panel letters | 183 mm, vector PDF and 600-dpi PNG |
| Required statements | Data availability, Code availability, Competing interests, Author contributions | present |

## 7. Later stages (for reference)

- **On revision:**
  - the Nature Portfolio Reporting Summary and editorial policy checklist (`reporting_summary_notes.md` has the
    answers for this computational study);
  - Source Data for the figures (the files in `data/outputs` are the source data);
  - figures as separate files (`figures/`).
- **On acceptance:**
  - the LaTeX source;
  - final figure files;
  - the publication of the Code Ocean capsule.
