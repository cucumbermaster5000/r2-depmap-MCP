# Medulloblastoma Gene Explorer — Version 3

Version 3 activates **Metastasis**, **Survival**, and clinical summaries in
**Overview**, with verified R2 coding, a predefined mean-expression survival
cutoff, statistically selected subgroup cohorts, and clinical downloads.
See [Version 3 methods, coding evidence, thresholds and validation](../docs/SHINY_V3.md).

Version 2 adds all 12 original Cavalli molecular subtypes in **MB Subtypes**,
with descriptive statistics, an overall comparison, all 66 pairwise comparisons,
BH correction, effect sizes and subtype downloads. Broad subgroup results and
downloads retain Version 1 behavior. See [Version 2 methods, metadata and validation](../docs/SHINY_V2.md).

Open `r2-depmap.Rproj` in RStudio, then run:

```r
shiny::runApp("shiny")
```

Alternatively open `shiny/app.R` in RStudio and click **Run App**.
The project `.Rprofile` adds `.R-library` to the library search path. If RStudio
was already open before installation, first run:

```r
.libPaths(c(".R-library", .libPaths()))
```

For a new machine, install dependencies from the project root:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/install-dependencies.R --shiny
```

To start from a terminal on this machine:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla -e '.libPaths(c(".R-library", .libPaths())); shiny::runApp("shiny", host="127.0.0.1", port=3876)'
```

Enter **HLX** or **SLC2A1** and click **Run analysis**. Results do not update
while typing. Use official symbols (SLC2A1 for GLUT1). Overview displays the
dataset, retrieved/analyzed sample counts, metadata fields and factual omnibus
test summary. MB Subgroups shows the unchanged core ggplot, statistics and group
summaries. Downloads provides patient CSV, raw statistics CSV, PDF and 300 dpi
PNG of the same plot. At narrow widths the figure can be scrolled horizontally.

Functional Biology, TF Targets and DepMap remain
**Not implemented yet**. No corresponding service is called by this application.

## Architecture

`app.R` loads the shared core and application configuration. `interface.R`
contains UI, reactive orchestration, friendly error messages and download
serialization. The only analysis calls are:

1. `core$get_r2_expression(cfg, gene)`
2. `core$analyze_r2_subgroups(retrieval)`
3. `core$analyze_r2_subtypes(retrieval)`
4. `core$analyze_r2_clinical_v3(cfg, retrieval, subgroup_analysis)`

These are the same functions used by the MCP subgroup workflow. No MCP server,
AI model or coding tool is needed to run Shiny. Statistical methods, plot objects
and raw source labels are unchanged. Adjusted p-values are explicitly marked
not applicable for the existing single omnibus test; no correction is invented.

Results are held per Shiny session. Retrieval already uses the existing
dataset/gene/reporter/transformation cache. There is no new global results cache.
The bslib task button disables itself and shows progress while processing;
`withProgress()` reports retrieval/analysis stages. Technical errors go to the R
console/server log. Failed runs clear the previous result and downloads.

This prototype uses synchronous R execution: a slow external request blocks
the R process until it completes. Multi-user background execution and cache
invalidation controls can be added later without changing analysis functions.
Keep the server local for this prototype; deployment is outside Version 1.

## Verification

Run `Rscript --vanilla tests/run-tests.R`. `test-shiny.R` checks button-only
execution, error recovery, source/analysis delegation and four export formats.
The browser acceptance check uses the real cached R2 HLX and SLC2A1 datasets,
763 patients each, and confirms plot/table display and download controls.
See `docs/SHINY_V1_ACCEPTANCE.md` for recorded results and limits.

UI implementation follows Posit's [task button documentation](https://pkgs.rstudio.com/bslib/reference/input_task_button.html).
