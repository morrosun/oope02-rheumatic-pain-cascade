# Analysis notes

## 1. Cascade definitions

**Pain.** Panel A uses the harmonised item asking whether the respondent is often troubled by
pain (`painfr`). Panel B uses any musculoskeletal site from a 15-site list (CHARLS, `da042s*`)
or a 13-site list (KLoSA, `pain_*`), after excluding headache, chest pain and stomach pain; and a
single joint-pain or stiffness item in LASI (`r1jointp`) and SHARE (`pain_s`).

**Diagnosis.** Self-reported physician diagnosis of arthritis or rheumatism (`arthre`;
`r1arthre` in LASI). The wording differs across cohorts — KLoSA asks about "rheumatism" and LASI
asks whether the respondent "has ever had" the condition — and this is carried into the
limitations rather than corrected.

**Treatment.** Current use of any medication for arthritis or rheumatism (`rxarthr`;
`rxarthr_c` in CHARLS; `r1rxarthr` in LASI). Two exceptions are structural, not missing at random:

- **ELSA** fielded no arthritis-medication item. ELSA therefore contributes to the diagnostic
  step only, and every statement about treatment is restricted to cohorts with a treatment item.
- **SHARE** asks only about analgesic use at least weekly (`rxpain`). SHARE treatment estimates
  are reported separately and never pooled with the others.

**Functional outcome.** Any limitation in activities of daily living: at least one difficulty
across the six harmonised ADL items (`adlab_c`). A 13-item ADL/IADL version is used only as a
sensitivity analysis (`26_w6_validation.py`), because the baseline and follow-up instruments must
measure the same construct for baseline adjustment to work.

## 2. Denominators

| Quantity | Denominator |
|---|---|
| Pain prevalence | All respondents aged 50–105 with a non-missing pain item and sampling weight |
| Diagnosis coverage | People with pain |
| Treatment coverage | People with pain **and** a diagnosis |
| Treated per 100 with pain | People with pain (the product of the two steps, computed from unrounded estimates) |

Respondents aged outside 50–105, without a sampling weight, or with a missing pain item are
excluded; the per-cohort attrition is audited in `40_missing_structure.py` and `44_pain_bounds_weighted.py`.

## 3. Weighting and standardisation

Cascade rates use each cohort's respondent weight (`wtresp` in HRS, `INDV_weight_post2` in
CHARLS, and the equivalent in the other cohorts; assembled by `16_build_weights.py`) and are
directly standardised to the WHO 2000–2025 standard population, renormalised over ages 50 and
above (eight five-year groups from 50–54 to 85+).

Confidence intervals come from a cluster bootstrap with 1,000 replicates: communities
(`communityID`, 450 units) in CHARLS, households elsewhere. Households are coarser than primary
sampling units, so intervals outside CHARLS may be understated and are not treated as exact.

**The regression models are not weighted.** They adjust for age, sex, education, wealth quintile,
residence, comorbidity count and smoking, and their estimates are conditional associations rather
than population-average effects. This is stated in the manuscript's Methods.

## 4. Pooling

Cohort-specific logistic models with individual-level cluster-robust standard errors are pooled on
the log-odds scale with DerSimonian–Laird random effects. Heterogeneity is reported as I², and
pooled estimates with I² above 85% are not interpreted as common effects — for example, the male
association with diagnosis is reported as a direction with cohort-specific magnitude, whereas the
male association with treatment (I² = 0) is reported as a common effect.

## 5. Script map

### Upstream exploration (not needed to reproduce the published estimates)

| Script | Purpose |
|---|---|
| `05_feasibility.py` | Which cohorts carry pain, diagnosis and treatment items |
| `07_wave_coverage.py` | Wave-level coverage of the cascade variables |
| `08_cascade_counts.py` | Unweighted descriptive cascade counts |
| `10`–`14` | CHARLS wave 6 structure, weights, ID format, skip logic, derived cascade |
| `15_weight_scan.py` | Candidate weight variables per cohort |

### Main analysis

| Script | Output |
|---|---|
| `16_build_weights.py` | `16_weights.csv` — (cohort, ID, wave, weight) |
| `17_cascade_weighted.py` | `17_cascade_tierA/B.csv`, `17_four_state.csv`, `17_missing.csv` |
| `33_std_ci.py` | `33_std_ci.csv` — weighted, age-standardised rates with bootstrap CIs |
| `19_determinants.py` | `19_determinants.csv`, `19_pooled.csv` |
| `31_outcome_gradient_v2.py` | `31_outcome_v2.csv`, `31_pooled_v2.csv` (main outcome model) |
| `24_incident.py` | `24_incident.csv`, `24_incident_pooled.csv`, `24_incident_desc.csv` |
| `42_fig1_funnel.py` | `42_fig1_funnel.csv` — same-denominator cascade for Figure 1 |
| `45_wave_weighted.py` | `45_wave_weighted.csv` — wave-level weighted rates (Table S3) |
| `44_pain_bounds_weighted.py` | `44_pain_bounds_weighted.csv` — selection bounds |

### Tables and figures

`18`, `20`, `23`, `25`, `28`, `32`, `35` assemble human-readable tables from the CSVs above.
`figures/fig1_cascade.R`, `fig2_coverage_forest.R` and `fig3_outcomes.R` read
`outputs/*.csv` and write PDF/SVG/TIFF.

### Sensitivity, audit and reconciliation

| Script | Purpose |
|---|---|
| `21_severity_sensitivity.py` | Adds pain severity to the diagnosis model |
| `26_w6_validation.py` | CHARLS wave 6 direct item vs 0→1 transition |
| `27_charls_w6_row.py` | CHARLS wave 6 as a comparator row |
| `34_cluster_se.py` | Alternative clustering |
| `37_four_state_v2.py` | Four-state description with "no pain, previously diagnosed" separated |
| `38_wq_linearity.py` | Wealth quintile as categorical rather than linear |
| `39_reconcile_share.py` | Reconciles the SHARE sample counts (484 weight-missing + 554 pain-missing = 1,038) |
| `40_missing_structure.py` | Missingness structure by cohort and country |

## 6. Conventions worth knowing

- **Two panels are never pooled.** SHARE appears in both, with different items and non-overlapping
  waves; its two estimates are an internal consistency check, not a replication.
- **`17_four_state.csv` uses three states by default** because ELSA and SHARE cannot be split by
  treatment; the four-state model is fitted only within cohorts that carry a treatment item.
- **The reference category excludes people with a previous diagnosis but no current pain.** They
  are a separate state, because in HRS they make up 45% of the otherwise "no pain" group and
  differ markedly in function.
- **Figure 1 rows are restricted to waves with a treatment item**, which is why its HRS diagnosis
  share (72.4%) is marginally below the all-wave value in the tables (72.8%) and its implied MHAS
  treatment share (41.9%) is below the all-wave value (44.0%).
