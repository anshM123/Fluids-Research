# Nature Portfolio Reporting Summary: answers for this study

The journal asks for the Reporting Summary, usually at the first revision. This study is purely computational, so
most sections do not apply. These are the answers to copy into the form.

## Statistics

The paper has no experiments, samples or statistical hypothesis tests. Two items do apply:
- **Uncertainties.** Each quoted error is a stated numerical error (a grid-shift scatter, a fit error, or a
  systematic shift measured by varying one numerical parameter), not a statistical confidence interval. Their
  definitions are in the Methods ("Error audit and grid shifting") and in Supplementary Note 2.
- **Pre-registered decision rule.** It compares the measured profile position with the registered predictions,
  using the thresholds stated in Supplementary Note 3.

## Software and code

- **Data collection:** no software was used for data collection. All data are the outputs of the authors'
  numerical code.
- **Data analysis:** custom code in Python 3.11 with NumPy 2.4.6, SciPy 1.17.1, Numba 0.67.0 and Matplotlib 3.11.2.
  It is available at https://github.com/anshM123/Fluids-Research (directory `ncs-submission/code`) under the MIT
  licence and will be archived at Zenodo [DOI].
- **Custom code central to the claims:** yes. It is provided for peer review; see `code_and_software_notes.md`.

## Data

- **Data availability statement:** as in the manuscript.
  - The converged profiles and all logged outputs are at https://github.com/anshM123/Fluids-Research (directory
    `ncs-submission/data`) and will be archived at Zenodo [DOI].
  - The pre-registration record and the GitHub push log are Supplementary Data 1 and 2.
  - The source data for every figure are in `data/outputs`.
- **Restrictions:** none.
- **Human participants, human data, sex and gender, race and ethnicity:** not applicable.

## Field-specific reporting

Choose "Physical sciences" (or the closest available option).
- Life-science and behavioural-science items: not applicable.
- Ecological, evolutionary and environmental items: not applicable.

## Research involving humans, animals, cells, clinical data, dual-use research

Not applicable.
