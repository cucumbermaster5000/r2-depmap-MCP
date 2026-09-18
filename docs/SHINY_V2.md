# Molecular subtype expression, Version 2

Run from the project root in RStudio:

```r
.libPaths(c(".R-library", .libPaths()))
shiny::runApp("shiny")
```

Or from PowerShell, then open http://127.0.0.1:3876:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/serve-shiny.R
```

Enter a symbol and click **Run analysis**. **MB Subgroups** retains Version 1.
**MB Subtypes** displays the stored subtype results and **Downloads** adds six
subtype exports: patient CSV, descriptives CSV, pairwise CSV, overall test CSV,
13 x 7.5 inch PDF and 3900 x 2250 PNG (300 dpi). Downloads never rerun analysis.

## Metadata inspection

The retrieved R2 Cavalli dataset `ps_avgpres_gse85217geo763_hugene11t`, GSE85217,
links to PMID 28609654. Its only molecular subtype field is `subtype`;
`subgroup` contains broad group membership. Source annotation cross-tabulation
confirms the following membership without inference from expression or strings.

| subgroup | subtype | annotated N |
|---|---|---:|
| wnt | wnt_alpha | 49 |
| wnt | wnt_beta | 21 |
| shh | shh_alpha | 65 |
| shh | shh_beta | 35 |
| shh | shh_gamma | 47 |
| shh | shh_delta | 76 |
| group3 | group3_alpha | 67 |
| group3 | group3_beta | 37 |
| group3 | group3_gamma | 40 |
| group4 | group4_alpha | 98 |
| group4 | group4_beta | 109 |
| group4 | group4_gamma | 119 |

Total: 763 annotated, zero missing subtype annotations. These are the original
Cavalli 2017 labels, not a conversion to another molecular classification.
No second classification is available, so no selector is displayed. The function
accepts an explicit `subtype_field` for independently analyzing another source
field; it never merges systems or assigns an unidentified field a publication.

## Reusable core

```r
source("R/load_core.R")
core <- load_r2_core(".")
cfg <- core$read_server_config("config/defaults.json")
retrieval <- core$get_r2_expression(cfg, "HLX")
subtypes <- core$analyze_r2_subtypes(retrieval)
subtypes$overall
subtypes$pairwise
subtypes$associations
print(subtypes$plot)
```

`R/subtype_analysis.R` has no Shiny or network dependency. `load_core.R` makes
the same function available to R, Shiny and future MCP integration. It returns
gene, status, all patient data, metadata audit, missing/exclusion counts,
descriptives, diagnostics, overall and pairwise tests, directional association
rows, ggplot and warnings. Original metadata joins by exact sample ID. All raw
rows remain, including exclusions with reasons. The existing retrieval remains
unmodified. A subtype error is isolated from the successful subgroup result.

## Statistical choices

The predefined scientific target is an unadjusted rank-based comparison of
independent tumour expression distributions. Before inference the function
inspects group counts, sample sizes, finite values, variance, IQR, range and
ties. It does not choose tests by repeatedly testing normality. With fewer than
two represented subtypes, any represented group below five finite observations,
or constant expression, it provides descriptive results and a warning and
withholds inferential testing rather than silently dropping the small group.

For these 763-patient analyses every subtype has at least 21 observations.
The overall test is tie-corrected Kruskal-Wallis, with 11 degrees of freedom.
Its non-negative effect estimate is `(H - k + 1)/(n - k)` (epsilon-squared).
There is one omnibus test per queried gene; its adjusted p-value is not applicable.

All 66 planned pairs are reported regardless of the omnibus p-value. Each uses
a two-sided Wilcoxon rank-sum test: exact when both groups have fewer than 50
samples and there are no ties, otherwise a tie-adjusted normal approximation
with continuity correction. An identical constant pair is explicitly assigned
p=1 and zero effect. BH adjusts the complete 66-pair family once per gene.
The table includes Mann-Whitney U, sample sizes, raw/adjusted p, median difference
and rank-biserial effect `2U/(nA*nB)-1`. Positive effect indicates higher ranks in
A relative to B. Association rows expose this same evidence without triggering
other analyses or declaring specificity from a p-value threshold.

Methods: [R Kruskal-Wallis documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/kruskal.test.html)
and [R Wilcoxon documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html).

Ranks do not remove confounding or unequal distribution-shape effects. These
tests cannot establish tumour subtype specificity, causality or treatment value.
Age, batch and other covariates are not adjusted; samples are treated as
independent according to the source primary-tumour cohort. Multiple queried genes
need additional multiplicity control. Missingness may be informative. No
survival, metastasis, enrichment, TF-target or dependency analysis is called.

## Verification commands and files

Acceptance on 2026-09-18 passed for all three genes: 763 patients, 12 subtypes,
66 pairs each; all Shiny renderers and ten export formats per gene passed.
Headless Edge confirmed both plots and all 30 actual browser downloads with
zero JavaScript errors. HLX and SLC2A1 Version 1 H and p-values matched the
recorded baseline. The full automated suite passed (one opt-in legacy live-R2
test skipped; three existing isolated-environment `stats` serialization warnings).

| Gene | Subtype H (df=11) | Raw overall p | Epsilon-squared |
|---|---:|---:|---:|
| HLX | 240.891657 | 2.249787e-45 | 0.306114 |
| MYC | 375.571878 | 9.304785e-74 | 0.485449 |
| SLC2A1 | 159.509704 | 1.685952e-28 | 0.197749 |

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla tests/run-tests.R
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/check-shiny-v2.R
```

Real-data acceptance outputs are in `outputs/shiny-v2-acceptance/`, with each
gene's structured RDS, all ten exports, and `acceptance.csv`. The script drives
the actual Shiny server with `testServer` and checks rendered tables/plot, counts,
BH values, exported precision and Version 1 HLX/SLC2A1 baseline statistics.
`scripts/check-shiny-v2-browser.cjs` additionally exercises a running app in
headless Edge using Playwright (`NODE_PATH` must locate an installed Playwright).
It checks HLX, MYC and SLC2A1 tabs, plot rendering, all 30 browser downloads and
JavaScript errors. Default test URL is http://127.0.0.1:3878; override with
`SHINY_TEST_URL`. Browser artifacts are under the acceptance folder's `browser/`.

Changed application files: `R/subtype_analysis.R`, `R/load_core.R`,
`shiny/interface.R`. Tests: `tests/testthat/test-subtypes.R`,
`tests/testthat/test-shiny.R`. Added scripts: `scripts/check-shiny-v2.R`,
`scripts/check-shiny-v2-browser.cjs`, `scripts/serve-shiny.R`.
Documentation: this file and `shiny/README.md`.
