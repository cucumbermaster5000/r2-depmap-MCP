# Integrated gene report

Ask Codex: **Use build_gene_page for HLX with median as the default cutoff.**
Restart Codex after updating the server so it reloads the seven tools. No global
configuration changes are needed if `r2-mb-v1` is already registered.

Or run from the project PowerShell terminal:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/gene-page.R HLX
```

Open the returned HTML path in a browser. It is self-contained. The cohort and
cutoff selectors switch between R-computed analyses; they do not recalculate
statistics in JavaScript. The report defaults to median, with mean and optimal
also included. `build_gene_page` accepts `cutoff="mean"` or `"optimal"` to change
the initial view. Publication parameters and optional explicit R2 reporters are
also supported. The first use downloads about 429 MB of DepMap data once.

## Evidence included

- Cavalli GSE85217: original patient expression, four broad subgroups, all 12
  source molecular subtypes, Kruskal-Wallis and pairwise Wilcoxon tests with BH
  correction. Pairwise median differences and rank-biserial effects are retained.
- Pfister: exactly `ps_avgpres_pfisterb272_informp3`, **Mixed Pediatric Pan Cancer -
  Pfister - 272 - fpkm - informp3**. All 272 source samples are retained. The CNS
  comparison uses 246 primary samples, excluding non-CNS diagnoses and non-primary
  samples. The original cancer labels and exclusions remain inspectable. Transform:
  local `log2(1 + FPKM)`; no pooling with Cavalli microarray values.
- Overall survival: Kaplan-Meier curves, pointwise confidence bounds where
  estimable, censor marks, risk tables, log-rank tests, Cox high-versus-low hazard
  ratios and PH diagnostics. Runs in the complete MB cohort and each original
  broad subgroup; insufficient samples/events produce explicit unavailable panels.
- Metastasis: M0/metastatic proportions and sample counts within low/high
  expression groups. Mean/median comparisons use Fisher exact tests, odds ratios
  and confidence intervals. Raw 0/1 codes remain intact.
- Enrichr: pinned GO Biological Process 2026, Reactome Pathways 2024, KEGG 2021
  Human and WikiPathways 2024 Human memberships. These are annotations, **not
  evidence of pathway deregulation**. No private transcriptomics data are used.
- TF evidence: Lambert human TF list v1.01; exact human ChEA 2022 target experiments
  when present. ChIP-associated targets are not called perturbation-confirmed or
  MB-specific targets. Target enrichment is separate for each experiment/library,
  with an explicit human ChEA annotation background and BH across all eligible
  terms. Absence of target evidence is reported, not replaced by invented targets.
- DepMap: official **24Q4 Public v1**, DOI `10.25452/figshare.plus.27993248.v1`.
  This is an archived release, **not the latest**. The live portal currently
  returns human verification to automated requests; the official Figshare archive
  supports reproducible downloads. Model and Chronos files are checked against
  archive MD5 hashes. All screened models provide context; ten models explicitly
  annotated as Medulloblastoma are highlighted. Missing gene effects stay missing;
  no subgroup is guessed from a cell-line name. Model-level expression/dependency
  correlation is not yet included.
- PubMed: recent symbol/alias matches, title, journal, citation/electronic dates,
  DOI and one labelled abstract excerpt, with query and retrieval provenance.

## Clinical definitions and limitations

Survival mapping is verified against R2's own patient-level Kaplan payload and
`kaplan.js`: status 1 decreases survival and agrees with `dead=1`; status 0 is
censored. Patient IDs are joined by unique case-normalised GEO accession. R2's
Kaplan follow-up is in months and retains precision lost in the rounded
`os_(years)` annotation. Models use these months divided by 12, checked against
the original years. For HLX there are 612 complete cases and 159 deaths.

Metastasis mapping is documented in R2's exact field
`met_status_(1_met__0_m0)`: 0 is M0; 1 is metastatic.

Mean/median thresholds use **log2 expression across all expression-complete
patients in the selected cohort**, before endpoint-specific exclusions. This
makes the same cutoff available for both clinical endpoints. High is strictly
greater than the threshold, Low is less than or equal; ties are never split.
Subgroup thresholds are recalculated within that subgroup. Results can differ
from R2's default survival-complete or raw-expression cutoffs.

Optimal means the smallest log-rank p-value among distinct expression cutoffs in
survival-complete patients, retaining at least ten patients and 10% per group.
The full search is saved. Bonferroni corrects across all eligible tested cutoffs.
Selected hazard ratios/CIs remain descriptive and potentially optimistic.
Metastasis at a survival-selected cutoff is **descriptive only**: its inferential
p-value and CI are omitted because the two endpoints may be correlated.

Clinical family correction uses BH across five cohorts times three modes,
separately per endpoint (including unavailable tests in family size). It does not
adjust across future gene queries. Results are exploratory and unadjusted for
age, subgroup, treatment and other confounders. Pooled associations can reflect
subgroup composition. No causal or therapeutic conclusions are automated.

## Reproducibility and testing

Each run creates `outputs/GENE/report-<UTC>-<id>/` with HTML, `report_data.rds`,
CSV tables, provenance, source-function snapshot, package versions and
`reproduce_report.R`. From that folder in RStudio, source the reproduction script
to rebuild the page offline. The saved object includes inputs and statistical
results; functions such as `profile_clinical` can rerun calculations explicitly.
R2 source snapshots and cached Enrichr/DepMap files are identified in provenance.

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla tests/run-tests.R
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/check-mcp-v1.R
```

Despite its historical filename, the second script tests the **current** stdio
server, checks seven tools, retrieves real cached HLX patient data and builds the
integrated report. It asserts successful main source modules. Source failures
are explicit in the report rather than silently replaced by empty results.

Primary source references: [R2](https://hgserver1.amc.nl/cgi-bin/r2/main.cgi),
[Enrichr](https://maayanlab.cloud/Enrichr/),
[Human TFs](https://humantfs.ccbr.utoronto.ca/download.php),
[DepMap archive](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248/1),
[NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/).
