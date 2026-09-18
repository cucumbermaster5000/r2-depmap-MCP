# Version 3: clinical analyses

From the project root in RStudio:

```r
.libPaths(c(".R-library", .libPaths()))
shiny::runApp("shiny", launch.browser=TRUE)
```

Stop any old app before restarting. Enter HLX, MYC or SLC2A1 and click Run
analysis. Overview, Metastasis and Survival now include clinical results;
MB Subgroups and MB Subtypes preserve their Version 1/2 results and exports.
Clinical results and downloads use the same stored per-session objects.

## Inspection and reuse

Before editing, the full existing suite passed. `get_r2_expression()` supplies
patient expression, original metadata and provenance. `analyze_r2_subgroups()`
returns the V1 omnibus result but did not select associated subgroups.
`analyze_r2_subtypes()` supplies V2 and remains unchanged. Existing
`profile_group_comparison()` supplies BH-corrected rank comparisons;
`get_cavalli_survival()` and `verify_r2_survival_payload()` retrieve and verify
R2 patient survival data. `profile_survival()` supplies KM, log-rank, curve and
risk-table calculations. Its default behavior remains unchanged; V3 adds an
optional Cox switch and access to its underlying fitted objects.

The older general `clinical_associations()` service uses continuous standardized
expression and imported cohorts. It does not implement this request's mean split,
so V3 does not substitute that scientifically different analysis.

## Verified clinical metadata

Source: Cavalli GSE85217, R2 table `ps_avgpres_gse85217geo763_hugene11t`.
Original annotations and tracks are saved with each retrieval. For HLX the
source is `outputs/HLX/20260918T075622-76a46b3f52b4/source/`.

| Field | Original source type and values | Counts |
|---|---|---|
| `met_status_(1_met__0_m0)` | character: `0`, `1`, `na` | 397 / 176 / 190 |
| `dead` | character: `0`, `1`, `na` | 467 / 165 / 131 |
| `os_(years)` | numeric text, including `na`; observed range 0–25 | 625 recorded / 138 missing |

The exact metastasis field in both annotation and track documentation explicitly
defines 1 as metastatic and 0 as M0. V3 retains the raw code and derives
`Metastatic` and `M0 (non-metastatic)` only after validation. There is no additional
M1/M2/M3 stage field in the retrieved metadata. Usable metastasis N=573.

Survival is **overall survival**. R2 Kaplan data use follow-up in months;
V3 divides by 12 to express years. These values have more precision than the
rounded `os_(years)` annotation, which is retained. No separate validated date of
last contact or additional follow-up field is available.

Event meaning is not guessed from a numeric code or the name `dead`:

1. Read the original saved R2 Kaplan payload and the actual R2 `kaplan.js` renderer.
2. Check the payload declares overall survival and a month-based time axis.
3. Confirm the renderer reduces survival using
   `proportion *= (nrRemaining - status) / nrRemaining;`.
   Thus status 1 is an event and status 0 does not reduce survival (censoring).
4. Match unique sample IDs and verify every status equals that patient's `dead`
   value; verify the Kaplan sample set equals cases with both annotations.
5. Check months/12 against rounded annotation years (existing verifier tolerance
   0.11 years) and reject invalid times/statuses or mismatches.

V3 revalidates these saved source files on each analysis rather than trusting a
cached boolean. The verification object includes source paths, renderer URL and
MD5 hashes. HLX evidence is in the retrieval's `survival_source/r2-kaplan.html`
and `survival_source/kaplan.js`; MYC and SLC2A1 have their own verified sources.
Any failure blocks survival interpretation, while expression results remain.
The legacy R2 retrieval requests a median-rendered page solely to obtain patient
follow-up/status data; its split and statistical results are ignored. All V3
survival grouping and statistics are computed locally using the mean only.

There are **612 usable survival cases, 159 deaths and 453 censored**, with 151
patients excluded for incomplete survival information. Missing event and time
counts overlap and therefore must not be added. All three acceptance genes have
763 finite expression values. Excluded patients remain in patient CSVs with reasons.

## Predefined clinical cohort selection

The full MB analysis always runs when verified data permit. A subgroup qualifies
only when all these conditions hold:

- V1 overall subgroup p < 0.05.
- All three comparisons against the other broad groups have BH p < 0.05.
- Every comparison has absolute rank-biserial effect >= 0.30 and at least ten
  patients in both groups.
- Rank direction and median difference agree and consistently indicate higher
  or consistently lower expression against every other group.

BH covers the six broad subgroup pairs. These conservative thresholds are an
explicit exploratory operational definition, not a validated clinical classifier.
Subtypes are not used to select cohorts. Selected subgroups are tested separately,
never pooled. No clinical outcome is used for selection. If none qualify, only
the full cohort is tested. Criteria and the six comparisons are visible in the app.

Acceptance selection: **HLX: group3; MYC: none; SLC2A1: group3**.
Synthetic tests also verify that WNT can be selected: no Group 3 assumption.

## Statistical rules

Metastasis: compare expression between verified metastatic and M0 groups.
Inspect counts, finite values, mean, median, SD, variance and IQR. Require >=5 per
group and variable expression for inference; otherwise retain descriptive output.
The predefined distribution question uses two-sided Wilcoxon rank-sum with ties
and continuity correction, reusing the established pairwise implementation.
Report U, raw p, rank-biserial effect and median difference oriented as metastatic
minus M0. This is not necessarily a median-only test when distribution shapes
differ. No high/low expression cutoff is imposed on metastasis in V3.

Survival cutoff: arithmetic **mean of all finite expression values within the
analyzed cohort**, before exclusion for missing clinical data. High > mean;
Low <= mean, including exact ties. Each subgroup gets its own mean. The reference
N and full-precision cutoff are stored and displayed. No optimization or median
fallback is permitted.

KM/log-rank: require >=20 usable patients, >=10 in each expression group, >=5
deaths overall and variable follow-up. Use `survival::survfit` with pointwise
log-log 95% CIs and `survival::survdiff` (rho=0). Show censor marks and a separate
number-at-risk table. Cox additionally requires >=10 deaths overall and >=3 per
expression group. Fit only expression group, with Low as reference; report HR
High/Low, 95% CI and Wald p. No multivariable model is fitted.

Unstable or non-finite Cox results are withheld. `survival::cox.zph` checks the
proportional-hazards assumption using scaled Schoenfeld residuals; p < 0.05
produces an explicit warning. A nonsignificant diagnostic does not prove PH.
Insufficient clinical data produce a message, with no forced fit.

Raw clinical p-values remain visible. Additional BH values cover requested cohorts
separately for metastasis, log-rank and Cox within each queried gene, including
unavailable fits in family size. They do not correct outcome-independent but
data-adaptive subgroup selection or multiple gene queries.

References: [R log-rank documentation](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/survdiff.html)
and [R proportional-hazards diagnostics](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/cox.zph.html).

## Reusable API and artifacts

New Shiny-independent functions in `R/clinical_v3.R`:

- `select_r2_clinical_subgroups(retrieval, subgroup_analysis)`
- `verify_r2_clinical_metadata(retrieval, survival_source)`
- `analyze_r2_metastasis(retrieval, verification, cohort="all")`
- `analyze_r2_survival(retrieval, verification, cohort="all", cutoff="mean")`
- `analyze_r2_clinical_v3(cfg, retrieval, subgroup_analysis)`
- `write_r2_clinical_download(analysis, kind, file)`

Helpers `clinical_v3_data` and `clinical_v3_counts` preserve metadata and counts.
Objects contain raw patient rows, eligibility/exclusion reasons, verification,
counts, estimates, warnings, plot and (for survival) KM/log-rank/Cox/PH objects and
risk table. The shared core loader exposes these to direct R and future MCP use.

```r
source("R/load_core.R")
core <- load_r2_core(".")
cfg <- core$read_server_config("config/defaults.json")
r <- core$get_r2_expression(cfg, "HLX")
clinical <- core$analyze_r2_clinical_v3(cfg, r, core$analyze_r2_subgroups(r))
clinical$survival$all$statistics
```

Downloads for each endpoint: `GENE_metastasis_patient_data.csv`,
`GENE_metastasis_statistics.csv`, `GENE_metastasis_plot.pdf/png`; equivalent
`GENE_survival_*` files. Subgroups include their source label, e.g.
`HLX_group3_survival_statistics.csv`. PNG is 300 dpi. Statistics exports include
exclusion counts, warnings and coding evidence; patient exports retain all cohort
rows and original metadata. Unavailable figures are not offered. V1/V2 downloads
remain present. Exports serialize stored results and never rerun models.

## Files

Created: `R/clinical_v3.R`, `tests/testthat/test-clinical-v3.R`,
`scripts/check-shiny-v3.R`, `scripts/check-shiny-v3-browser.cjs`, `docs/SHINY_V3.md`.
Modified: `R/load_core.R`, `R/profile_statistics.R`, `shiny/interface.R`,
`shiny/README.md`.

## Validation and limitations

`scripts/check-shiny-v3.R` tests real HLX/MYC/SLC2A1 retrievals, source verification,
selection, mean cutoffs, counts, clinical exports and Shiny renderers. It checks
V1 H/p and V2 H/p and every pairwise row against saved acceptance results.
`tests/run-tests.R` includes missing/unverified metadata, exclusion accounting,
direct survival-package agreement, exact ties, no optimal cutoff, low-event
guards and subgroup-selection tests. The legacy opt-in live-R2 test is skipped;
the real-data acceptance script covers the clinical source separately.
Three existing isolated-environment `stats` serialization warnings persist.

Artifacts: `outputs/shiny-v3-acceptance/` (structured RDS, CSVs, figures and
acceptance table). The Playwright browser script targets port 3879 by default,
checks Overview, all four analysis tabs, cohort eligibility and old/new downloads.
Browser acceptance completed for HLX, MYC and SLC2A1: **70 downloads passed**
(30 existing V1/V2 and 40 clinical), with zero JavaScript errors. It explicitly
checked that MYC has no Group 3 clinical section. The full automated suite passed.

| Gene | Full-cohort metastasis p | Full-cohort log-rank p | HR High/Low |
|---|---:|---:|---:|
| HLX | 0.470101 | 0.0149048 | 1.51970 |
| MYC | 0.350306 | 0.00441659 | 1.57234 |
| SLC2A1 | 0.322947 | 0.802674 | 1.04040 |

These are unadjusted p-values, not causal or clinically validated conclusions.
For example HLX qualifies for Group 3 expression enrichment but its Group 3
log-rank p is 0.947525; expression selection is not clinical evidence.
The pooled cohort can reflect subgroup composition. Age, treatment and batch are
not adjusted. Informative missingness/censoring, independent-sample assumptions,
data-adaptive selection, approximate rank p-values and multiple queried genes
limit interpretation. No age/treatment multivariable model or optimal cutoff is
used. The unofficial R2 interface can change or be unavailable; unresolved coding
blocks interpretation rather than silently recoding. The Shiny process remains
synchronous, so uncached retrieval may temporarily block other sessions.

Enrichr, Functional Biology, TF targets and DepMap remain unimplemented in Shiny.
