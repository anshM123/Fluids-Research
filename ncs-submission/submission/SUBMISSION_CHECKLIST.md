# Submission checklist: *Nature Computational Science*, Article

This is everything between this package and a submitted manuscript, in order. Items marked **[authors]** are
decisions or facts that only the authors can supply. Everything else is already done and is listed so that it can be
checked.

## 1. What to upload at initial submission

| Upload as | File | Status |
|---|---|---|
| Manuscript (PDF; includes Methods, the 6 figures and Extended Data Tables 1–2) | `manuscript/main.pdf` (source `manuscript/main.tex`, figures in `figures/`) | ready, once the placeholders in section 2 are filled |
| Cover letter | `submission/cover_letter.pdf` (source `.tex`) | ready, once signed |
| Supplementary Information | `supplementary/supplementary.pdf` | ready |
| Supplementary Data 1 | `supplementary/SupplementaryData1_preregistration_record.md`, a copy of the pre-registration record (plain text) | ready |
| Supplementary Data 2 | `supplementary/SupplementaryData2_github_push_log.csv`, the GitHub server times of every commit | ready |

The journal accepts initial submissions in any reasonable format (PDF, Word or LaTeX), provided the manuscript
contains the names and affiliations of all authors. The PDF is the simplest option. LaTeX source can be uploaded
later, at acceptance.

## 2. Placeholders to fill **[authors]**

| Where | Placeholder |
|---|---|
| `manuscript/main.tex` lines 38–42 | `[Author names]`, `[Affiliations]`, `[corresponding author, e-mail]` |
| `manuscript/main.tex`, Methods, "Use of AI tools" (red) | the drafted disclosure: `[system, version and provider]`, `[dates]`, `[describe the checks made by the authors]`; then delete the red note |
| `manuscript/main.tex`, Data and Code availability | `[DOI]` (twice). It can stay as is at submission; see section 5 |
| `manuscript/main.tex`, end matter | `[Funding and acknowledgements.]`, `[Contributions, in the CRediT taxonomy.]`, `[corresponding author]` |
| `supplementary/supplementary.tex` line 31 | `[Author names]` |
| `submission/cover_letter.tex` | sender block, `[Date]`, preprint sentence (keep one alternative), signature |
| `LICENSE`, `data/LICENSE-DATA.md` | `[Author names]` (copyright holders); `[DOI]` |

The repository description on GitHub names Ansh Mishra and Aryan Senthilkumar. Check the author list, its order
and the corresponding author before filling these in.

After editing, rebuild the PDFs:

```bash
cd manuscript    && latexmk -pdf main.tex         && cd ..
cd supplementary && latexmk -pdf supplementary.tex && cd ..
cd submission    && latexmk -pdf cover_letter.tex  && cd ..
```

This needs TeX Live 2023 or later with `latexmk`; it was built with pdfTeX 1.40.25. The figures do not need to be
regenerated.

## 3. Disclosure of AI use **[authors]**

Springer Nature's policy has three requirements:
- large language models do not qualify as authors;
- their use must be documented in the Methods;
- generative-AI images are not accepted.

All figures in this package are plots of computed data, drawn by `code/figures/make_figures.py`; none is
AI-generated.

The Methods paragraph "Use of AI tools" is a draft. It has to state accurately which systems were used, for what
(code, computations, analysis, figures, text) and when, and how the authors verified the results. The authors are
accountable for every statement in the paper. Read the whole manuscript and Supplementary Information with that in
mind before submitting.

## 4. Checks before submission **[authors]**

- [ ] **References.** Confirm authors, titles and identifiers of the two most recent preprints: Wang et al.,
  arXiv:2509.14185 (ref. 6), and Wang, Léger, Lai & Buckmaster, arXiv:2511.22819 (ref. 8). Replace any preprint
  that has since been published with its journal version: Chen & Hou I and II (refs 2, 3) and Buckmaster et al.
  (ref. 10).
- [ ] **Overlap with ref. 8.** If Wang et al. (arXiv:2511.22819) report IPM profiles, compare their λ values with
  the λ₁–λ₆ of this paper (`data/README.md`) and add one sentence to the main text on agreement or difference.
  The comparison was not possible from the environment in which the package was prepared.
- [ ] **Public record.** The paper states that predictions were committed to a public repository before the
  computations. The repository is public now, and GitHub's public event log covers it from its creation on
  28 September 2026, which indicates it was public throughout. Confirm that it was never private. If it was, change
  "public repository" to "version-controlled repository, now public" in the main text, the Methods and Supplementary
  Note 3.
- [ ] **Push log.** `SupplementaryData2_github_push_log.csv` was retrieved from GitHub's public events interface,
  which keeps 90 days. Keep this file. It cannot be regenerated after late December 2026.
- [ ] **Read-through.** Read the full manuscript once against `manuscript/main.pdf`, for anything you cannot vouch
  for.

## 5. Archiving

- [ ] **Zenodo.** In Zenodo, enable the GitHub integration for `anshM123/Fluids-Research`, then create a GitHub
  release (for example `v1.0-ncs-submission`). Zenodo archives the release and assigns a DOI. Put the DOI into the
  Data and Code availability statements and the licence files. A DOI is not required at initial submission; "will be
  archived at Zenodo" is acceptable until acceptance.
- [ ] **Code Ocean.** *Nature Computational Science* peer-reviews custom code that is central to a paper. When the
  manuscript goes to review, the editors will ask for a Code Ocean capsule and the Code and Software Submission
  Checklist (`submission/code_and_software_notes.md` has the answers). The capsule is private to reviewers until
  publication.
- [ ] **Default branch (optional).** The repository's default and only branch is
  `claude/fluid-dynamics-research-qo72sn`. To present a cleaner URL, rename it to `main` (GitHub: Settings →
  Branches → rename). GitHub redirects the old name, and the commit hashes cited in the paper are unchanged.

## 6. In the submission system

- Article type: **Article**.
- Title: *A resonance phase locates unstable self-similar singularities beyond an exponential precision wall*
  (12 words).
- Abstract: paste from `manuscript/main.tex` (148 words).
- Authors, affiliations, ORCIDs and the corresponding author **[authors]**: see `author_information.md`.
- Suggested and excluded referees (optional) **[authors]**: see `suggested_reviewers.md`.
- Competing interests: none declared (confirm) **[authors]**.
- Preprint: posting on arXiv is compatible with the journal's policy. If you post one, give its identifier in the
  cover letter **[authors]**.

## 7. Compliance with the Article format

The limits are from the journal's guide to authors (Content types; Formatting your initial submission). Check them
against the current guide at https://www.nature.com/natcomputsci/submission-guidelines, since they change.

| Item | Limit | This manuscript |
|---|---|---|
| Main text (excluding abstract, Methods, references, legends) | about 3,500 words | 3,290 |
| Abstract | 150 words, no references | 148 |
| Display items (figures and tables) | 6 | 6 figures |
| Extended Data items | 10 | 2 tables |
| References | about 50 | 37 |
| Methods | no strict limit | about 1,000 words |
| Figures | sized for 89 mm or 183 mm, sans-serif labels, bold lowercase panel letters | 183 mm, vector PDF and 600-dpi PNG |
| Required statements | Data availability, Code availability, Competing interests, Author contributions | present (contributions to be filled in) |

## 8. Later stages (for reference)

- **On revision:**
  - the Nature Portfolio Reporting Summary and editorial policy checklist (`reporting_summary_notes.md` has the
    answers for this computational study);
  - Source Data for the figures (the files in `data/outputs` are the source data);
  - figures as separate files (`figures/`).
- **On acceptance:**
  - the LaTeX source;
  - final figure files;
  - the Zenodo DOI;
  - the publication of the Code Ocean capsule.
