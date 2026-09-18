# Integrated medulloblastoma gene-report MCP

This repository publishes the R-based Model Context Protocol (MCP) server. It
exposes medulloblastoma R2, DepMap, PubMed, enrichment and clinical analysis
tools over MCP's newline-delimited JSON-RPC stdio transport.

The server entry point is `R/server.R`; it loads the shared analysis core and
registers the tool implementations in `R/v1_tools.R`. The web app is not part
of this published repository.

The current report includes Cavalli subtypes, Pfister CNS expression, survival and metastasis with selectable cutoffs, Enrichr, archived public DepMap dependencies and PubMed. **Start with [the current report guide](docs/GENE_REPORT.md).** The initial Version 1 audit below is retained as historical context; its deferred scope has now been extended.

# Medulloblastoma R2 MCP — Version 1

Version 1 retrieves **one gene's patient-level expression and original metadata**
from the Cavalli R2 cohort, saves the source responses, and makes one basic R
subgroup comparison and plot. The requested extension adds recent PubMed
publications and a self-contained HTML page; other biological modules remain deferred.

## Open a gene page in your browser

Run from this project folder:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/gene-page.R HLX
```

Open the returned `html_path` in your browser. It includes the R subgroup plot,
statistics, searchable patient data and ten recent PubMed publications with title,
journal, publication/electronic dates, DOI link and one sourced abstract sentence.
Everything needed to read the report is embedded; no web server, PDF or JPG is needed.
In RStudio you can instead run `source("scripts/gene-page.R")`.

Ask the MCP client: **Build a browser gene page for HLX.** This calls
`build_gene_page` with `{"gene":"HLX"}`. For literature alone use
`get_gene_publications`. Both support explicit aliases and an optional disease
restriction (`publication_disease` for the page, `disease` for literature).
The script accepts `--medulloblastoma` and `--refresh-publications`.

Literature is searched across all contexts by default, using the symbol and
verified human NCBI Gene aliases. Results follow PubMed's descending publication
date order. A bounded set of recent candidates is screened for a standalone
symbol/alias in the title or abstract, excluding drug codes such as HLX-02.
The report records the query, screening counts, exclusions and retrieval date.
This is keyword evidence, not manual confirmation of each paper's gene identity.
Descriptions are labelled verbatim abstract excerpts, with title-only fallback
when an abstract is absent. Missing DOI/date fields stay missing. Citation and
electronic dates are shown separately, preserving partial dates.

PubMed results are cached for 24 hours; explicit refresh bypasses the cache.
Raw NCBI responses and structured results are saved under `publications/` beside
the HTML. Retrieval requires internet access to `eutils.ncbi.nlm.nih.gov`;
failures appear explicitly in the page. No patient data is sent to NCBI.
The connector uses documented [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/).

Read the [source and project audit](docs/VERSION1_AUDIT.md) first. It addresses
metadata, clinical coding, identifiers, deferred modules and the later architecture.
Broader R modules written during the earlier request remain on disk for review,
but **the current MCP server does not load or expose them**.

## Verified result

Live HLX retrieval on 18 September 2026 returned 763 expression values and 763
matching metadata records. Reporter: `7909890`; R2 transformation: `transform_log2`.
Dataset: `ps_avgpres_gse85217geo763_hugene11t`; GEO accession: `GSE85217`.

- Complete subgroup and subtype labels for all 763 samples.
- Metastasis: 397 `0`, 176 `1`, 190 `na`. The source field explicitly declares
  `met_status_(1_met__0_m0)`: 0 = M0, 1 = metastatic. Raw codes remain intact.
- Survival time: 625 available; `dead` status: 632 available; both: 612.
  Event semantics are not independently codebook-verified, so no event conversion
  or survival analysis is performed.
- HLX subgroup test: Kruskal–Wallis H(3) = 193.0617, p = 1.33111e-41,
  763 samples, no exclusions. This is an overall distribution comparison,
  not a claim of subgroup specificity or therapeutic value.

## 1. Open in RStudio

Open [r2-depmap.Rproj](r2-depmap.Rproj), then run:

```r
source("scripts/version1.R")
head(retrieval$data)
View(retrieval$metadata)
retrieval$audit
validation$statistics
validation$summaries
print(validation$plot)
```

The default gene is HLX. To request another single gene after sourcing the script:

```r
retrieval <- get_r2_expression(cfg, "SLC16A1")
validation <- plot_r2_subgroup_test(retrieval)
print(validation$plot)
```

Dependencies are installed for R 4.6.1 on this machine. For a new installation:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/install-dependencies.R --tests
```

Version 1 uses jsonlite, curl, xml2, rvest, data.table and ggplot2. No limma/edgeR
installation is required; the earlier optional `--analysis` installer is retained
for later matrix-analysis work.

## 2. Connect Codex

| Enabled MCP tool | Purpose |
|---|---|
| `describe_r2_dataset` | Show source and Version 1 scope |
| `get_r2_expression` | Retrieve expression, original metadata, mappings and file paths |
| `plot_r2_subgroup_test` | Compute the R plot and save an omnibus statistical test |
| `get_gene_publications` | Retrieve recent PubMed articles and provenance |
| `build_gene_page` | Build the local HTML report with embedded R plot and publications |

Register the local stdio server from PowerShell:

```powershell
codex mcp add r2-mb-v1 -- 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla 'C:\Users\dng\Documents\r2DepmapMCPtest\r2-depmap\R\server.R'
```

Alternatively, copy [config/codex-mcp.example.toml](config/codex-mcp.example.toml)
into your Codex MCP configuration; it also sets a 300-second tool timeout.
Use one registration method, restart/reload the client and check that `r2-mb-v1`
is listed. These instructions follow the
[official Codex MCP documentation](https://learn.chatgpt.com/docs/extend/mcp?surface=cli).
Your global Codex configuration has not been modified by this task.

Then ask:

> Get R2 expression for HLX.

Expected tool call: `get_r2_expression` with `{"gene":"HLX"}`. The default response
previews 20 rows; the full CSV/RDS contains every sample. Use `offset` and `limit`
to page through rows, or `refresh=true` to create a new source snapshot.

Next ask:

> Make the Version 1 HLX subgroup test plot, then stop.

The server can also start directly:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/run-server.ps1
```

It uses newline-delimited MCP JSON-RPC on stdin/stdout; diagnostics go to stderr.
The public guest session is negotiated automatically. First retrieval needs
network access to `hgserver1.amc.nl`; no credentials were needed for this cohort.

## 3. Inspect and reproduce outputs

Each retrieval creates `outputs/GENE/<UTC-run-id>/`:

```text
source/                       Exact R2 HTML and JSON responses
GENE_R2_patient_data.csv       Expression and explicitly derived fields
GENE_R2_metadata.csv           Complete original metadata, labels and IDs
metadata_audit.json            Fields, values, counts, mappings and limitations
provenance.json                Requests, transformations, timestamps and hashes
retrieval.rds                 Full structured R object
GENE_gene_page.html            Browser page with embedded vector plot
publications/                 PubMed results and raw source responses
GENE_subgroup_test.csv         Test statistic, effect estimate and p-value
GENE_subgroup_summaries.csv    n, mean, median, SD, IQR and range
subgroup_test_exclusions.csv   Explicit exclusions (empty for verified HLX run)
subgroup_test_method.txt       Assumptions and limits
analysis_functions.R          Snapshot of exact R test/plot functions
reproduce_plot.R               Offline regeneration from saved data
session-info.txt               Retrieval R/package versions
analysis_session-info.txt      Analysis R/package versions
```

For an offline rerun, set the RStudio working directory to the run folder and
source `reproduce_plot.R`. Results go into `rerun/` and the plot remains accessible
as an R object. PDF/PNG export is optional via `export_images = TRUE` in
`plot_r2_subgroup_test`; it is off by default. Existing earlier exports are retained.

Cache: `.cache/r2-depmap/`, overridable with `R2_DEPMAP_CACHE_DIR`.
Original R2 `samplenames` are the exact join key: they are methylation GEO IDs in
this cohort. `geo-gse85217` holds the separate expression accession; both are retained.

## 4. Validate, then stop

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla tests/run-tests.R
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/check-mcp-v1.R
```

The second check uses actual cached HLX data through the real stdio server;
run `scripts/version1.R HLX` first on a fresh checkout. Requests/responses/stderr
are saved under `outputs/audit/`. Tests also retain earlier-module coverage;
that does not enable those modules in V1. The old live differential-expression
test is opt-in and was not run; the new patient connector was verified separately live.

R2's unofficial interface can change. Checks cover dataset, gene, reporter,
transformation, sample count and exact sample-ID matches. Unexpected responses
fail clearly. Multiple reporters require explicit selection. Kruskal–Wallis
compares distributions; unequal shapes/variances prevent a pure median-test
interpretation. This basic unpaired check does not adjust for confounding.

**Historical Version 1 scope (superseded by the current report guide).** Survival, metastasis testing, subtype statistics,
Enrichr, TF targets, DepMap, subgroup selection, the complete biological report and multi-gene
analysis are deferred. Earlier arbitrary ranking scores are excluded from the
requested evidence-separated scientific design.
