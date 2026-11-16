# Rheumatic pain care cascade across seven ageing cohorts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24211369.svg)](https://doi.org/10.5281/zenodo.24211369)

Analysis code for a cross-national study of the care cascade for rheumatic pain — from
reporting pain, to carrying a diagnosis of arthritis or rheumatism, to receiving treatment —
in seven harmonised ageing cohorts: CHARLS (China), ELSA (UK), HRS (USA), KLoSA (Korea),
LASI (India), MHAS (Mexico) and SHARE (Europe, 19 countries).

The study treats measurement as the primary design constraint. Because the pain items are not
interchangeable, the code builds two panels that are analysed separately and never pooled:

| Panel | Cohorts | Pain item |
|---|---|---|
| A — frequent pain | HRS, MHAS, ELSA, SHARE waves 5–8 | "often troubled by pain" |
| B — joint pain | CHARLS, KLoSA, LASI, SHARE waves 1/2/4 | site list, or a single joint-pain item |

## What this repository contains

- `analysis/` — 32 Python scripts: variable extraction, cascade construction, weighted and
  age-standardised rates, cohort-specific logistic models with DerSimonian–Laird pooling,
  sensitivity analyses and reconciliation checks.
- `figures/` — three R scripts (ggplot2 + patchwork) that reproduce the figures from the
  CSV files written by the analysis scripts.
- `docs/analysis_notes.md` — definitions, denominators, weighting and standardisation, and a
  script-by-script map of what each file produces.

## What this repository deliberately does not contain

- **Any respondent-level data.** All seven cohorts are distributed by their own providers under
  registration or data-use agreements; none may be redistributed. `data/` is empty apart from a
  README describing the expected layout.
- **Manuscript, table and figure assembly code**, and the scripts that fetch and renumber
  references. Those files contain unpublished text and are not part of the analysis.
- **Derived result tables.** `outputs/` is created at run time and is git-ignored; every table and
  figure is regenerated from the code plus the harmonised files.

## Data

Download the harmonised files from the Gateway to Global Aging Data (https://g2aging.org/) after
registering with each cohort, then lay them out as described in `data/README.md`. Three
environment variables override the defaults:

| Variable | Default | Meaning |
|---|---|---|
| `AGING_DB_ROOT` | `<repo>/data` | Root of the harmonised files |
| `AGING_DB_RAW_ROOT` | `<repo>/data_raw` | CHARLS wave 6 raw module (validation only) |
| `OOPE02_OUT` | `<repo>/outputs` | Where CSV/TXT results are written |

## Reproducing the analysis

Python 3.11 with `pandas`, `numpy`, `pyreadstat`, `statsmodels`, `scipy`.
R 4.6 with `ggplot2`, `patchwork`, `dplyr`, `readr`, `svglite`, `ragg`, `scales`.

Run from the repository root, in this order:

```bash
python analysis/16_build_weights.py        # sampling weights per cohort-wave
python analysis/17_cascade_weighted.py     # cascade steps, weighted + age-standardised
python analysis/33_std_ci.py               # cluster-bootstrap confidence intervals
python analysis/42_fig1_funnel.py          # same-denominator cascade for Figure 1
python analysis/45_wave_weighted.py        # wave-level weighted rates (Table S3)
python analysis/19_determinants.py         # determinants of diagnosis and treatment
python analysis/31_outcome_gradient_v2.py  # functional limitation across cascade states
python analysis/24_incident.py             # incident diagnosis and follow-up function
python analysis/44_pain_bounds_weighted.py # selection bounds for pain prevalence
Rscript figures/fig1_cascade.R             # needs 42_
Rscript figures/fig2_coverage_forest.R     # needs 33_
Rscript figures/fig3_outcomes.R            # needs 31_
```

All scripts write to `outputs/`; the R scripts write PDF/SVG/TIFF files plus their source-data
CSV to `outputs/figures/`.

The remaining scripts are upstream exploration and data construction (`05`, `07`, `08`, `10`–`15`), reconciliation checks (`39`, `40`), or
alternative specifications reported as sensitivity analyses (`21`, `26`–`28`, `34`, `35`, `37`, `38`).
`docs/analysis_notes.md` maps every script to its output.

## Citation

If you use this code, cite the archived release:

> Wang K. *Analysis code for: The rheumatic pain care cascade across seven ageing cohorts*
> (v1.1.1). Zenodo. https://doi.org/10.5281/zenodo.24211369

Concept DOI (always the latest version): https://doi.org/10.5281/zenodo.24211369

## Versioning and the Zenodo record

The canonical archive is the concept DOI `10.5281/zenodo.24211369`, which always resolves to the
latest version (currently v1.1.1). This
repository also has the Zenodo–GitHub integration enabled, but the archived record was created
directly through the Zenodo API and cannot be adopted by that integration.

**Do not publish a new GitHub release in order to archive a new version** — Zenodo would create a
second, independent record with its own concept DOI, and the same code would end up with two DOIs.
To release a new version, open the record on Zenodo and use *New version*; the concept DOI then
keeps resolving to the latest version.

## License

MIT — see `LICENSE`.
