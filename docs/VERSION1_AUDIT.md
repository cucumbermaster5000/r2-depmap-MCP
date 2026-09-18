# Source audit and Version 1 review

Audit date: 18 September 2026. Cohort counts were computed from retrieved R2 data.
Preserved audit: `outputs/audit/`. Validated gene run:
`outputs/HLX/20260918T075622-76a46b3f52b4/`.

## Existing implementation and R2 connection

The original repository contained an R stdio MCP server, five tools, an unofficial
two-group R2 differential-expression connector, local DepMap CSV readers, a
weighted expression/dependency ranking function, a cache, and fixture tests.
It did not retrieve patient-level single-gene expression or complete metadata.

| Files | Role and current status |
|---|---|
| `R/common.R`, `R/cache.R` | Configuration, JSON, local cache; used by V1 |
| `R/r2_client.R` | HTTP helpers and historical differential-expression workflow; V1 uses HTTP helpers only |
| `R/server.R` | MCP transport; now loads only V1 modules |
| `R/r2_expression.R` | New patient retrieval, metadata audit, exact joins, basic test and plot |
| `R/v1_tools.R` | Three V1 MCP tools with argument validation |
| `R/depmap.R`, `R/ranking.R`, `R/tools.R` | Original dependency/ranking tools; not exposed by V1 |
| `R/study.R`, `R/clinical.R`, `R/enrichment.R`, `R/integration.R`, `R/study_tools.R` | Draft extensions from the earlier request; retained for review, not loaded by V1 |
| `scripts/version1.R` | R/RStudio entry point for one gene and one plot/test |
| `scripts/inspect-r2.R` | Read-only audit of observed public interface requests |
| `tests/testthat/test-v1.R`, `test-mcp.R` | Parser, mapping, statistical export and protocol checks |

The draft clinical module uses continuous expression and assumed import
conventions; it is not the approved final clinical workflow for this cohort.
Future survival work must implement the requested predefined mean cutoff after
endpoint verification. Existing priority scores are also excluded from the new
evidence-separated design.

The historical connector submitted web forms for differential expression. Live
inspection found an initial guest-access page it did not handle. V1 follows the
site's explicit guest action, observed gene-reporter search and metadata JSON
requests, then the patient table embedded in the gene-view page. DataGrabber is
unavailable to a guest; the gene-view patient data are publicly accessible.

The audit preserves public JavaScript documenting the request names. R parses
the identified JSON argument without executing JavaScript. This is an unofficial
read-only interface, not a stable documented external API.

## Dataset, expression and identifiers

- R2 table: `ps_avgpres_gse85217geo763_hugene11t`.
- Label: `Tumor Medulloblastoma - Cavalli - 763 - rma_sketch - hugene11t`.
- Accession: `GSE85217`; platform: `hugene11t`; species: human.
- Source description: 763 fresh-frozen primary medulloblastoma samples.
- Dataset publication: PubMed `28609654`, Cavalli et al. (2017).
- HLX resolves to reporter `7909890`. Every requested gene is resolved anew;
  ambiguous reporters require explicit selection.
- All 763 HLX expression values are finite. Requested and returned transformation:
  `transform_log2`. These values are not RNA-seq counts or TPM.

R2 explicitly says sample names were renamed to the corresponding methylation GEO
IDs. Sample `GSM2260745`, for example, has separate expression GEO annotation
`gsm2261538` and lab ID `mb_subtypestudy_55001`. V1 joins expression `samplenames`
to annotation `samplenames` by exact equality, retaining every other ID unchanged.

[R2 cohort](https://hgserver1.amc.nl/cgi-bin/r2/main.cgi?table=ps_avgpres_gse85217geo763_hugene11t),
[source publication](https://pubmed.ncbi.nlm.nih.gov/28609654/).

## Metadata actually available

The patient endpoint returns these 16 fields, encoded as strings:

`samplenames`, `agegroup`, `dead`, `histology`, `met_status_(1_met__0_m0)`,
`subgroup`, `subtype`, `_data_set`, `age`, `gender`, `geo-gse85217`, `lab-id`,
`os_(years)`, `tissue`, `type`, `link`.

The track catalogue also mentions an internal `id` (0–762), but it is not in the
returned patient table and is not manufactured. Full per-field values/counts
are saved in `metadata_audit.json`. Raw `na` tokens remain intact; derived fields
explicitly treat `na`, `NA`, empty strings and JSON null as missing.

| Evidence | Available | Missing |
|---|---:|---:|
| HLX expression | 763 | 0 |
| Broad subgroup | 763 | 0 |
| Molecular subtype | 763 | 0 |
| Metastasis 0/1 annotation | 573 | 190 |
| Numeric `os_(years)` | 625 | 138 |
| `dead` 0/1 annotation | 632 | 131 |
| Both time and recorded status | 612 | 151 |
| Expression, subgroup, subtype, metastasis, time and status together | 538 | 225 |

Availability does not by itself establish clinically valid endpoint definitions.

## Broad subgroups and all observed subtypes

Broad labels: `wnt` (70), `shh` (223), `group3` (144), `group4` (326).
One molecular classification field, `subtype`, was exposed in this public metadata:
12 labels, no missing entries. No alternative system was returned; this does not
prove that none exists elsewhere.

| Broad subgroup | Exact subtype label | n |
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

Relationships come from cross-tabulating the source columns, not parsing label
prefixes or applying a classifier. The cohort is linked to Cavalli 2017 by R2.
R2 supplies no separate annotation-version identifier or subtype-specific citation
field. No relabelling to Greek characters, alternative classification mapping or
subtype statistical analysis was performed.

## Metastasis coding verification

Exact field: **`met_status_(1_met__0_m0)`**. Both the annotation response and track
catalogue carry this source label. It explicitly states `1_met` and `0_m0`;
the mapping is not inferred merely from seeing binary values.

| Original value | n | Separate derived label |
|---|---:|---|
| `0` | 397 | M0 (no metastasis) |
| `1` | 176 | Metastatic |
| `na` | 190 | Missing, no label assigned |

Evidence: `source/annotations.json` and `source/tracks.json` in the gene run.
The initial audit copies are named `annotations.html` and `track-info.html`,
but contain JSON. The interface does not document assessment time or detailed
M1/M2/M3 stage; these are not inferred. If the exact source declaration is absent
or unexpected codes appear, V1 does not silently apply the mapping. No metastasis
association analysis is performed.

## Survival variables and remaining uncertainty

`os_(years)` explicitly specifies years. Observed numeric values range from 0 to
25, with 138 `na` entries. `dead` contains `0` (467), `1` (165), `na` (131).

The field name suggests vital status, but the retrieved source does not include
an independent event-code dictionary. V1 therefore preserves
`survival_status_raw`, leaves `survival_event` unassigned, and does not run a
survival model. Future work needs authoritative event coding, endpoint and time
origin before mean-expression splitting, Kaplan–Meier/log-rank or Cox analysis.

## Version 1 data, statistics and output

`retrieval$data` holds queried gene, reporter, exact sample ID, numeric expression,
original expression string, original subgroup/subtype labels, metastasis raw code,
separate verified label, raw survival time, numeric years, raw status and an
unassigned event field. `retrieval$metadata` preserves the complete original
table. `retrieval$audit` and `$provenance` explain the source and derivations.

HLX: **Kruskal–Wallis H(3)=193.0617, p=1.33111e-41**, epsilon-squared 0.2504106,
763 samples, zero exclusions. Mean/median expression on the returned log2 scale:

| Subgroup | n | Mean | Median |
|---|---:|---:|---:|
| group3 | 144 | 7.2851 | 6.3659 |
| group4 | 326 | 5.5483 | 5.5251 |
| shh | 223 | 5.7407 | 5.7060 |
| wnt | 70 | 5.6992 | 5.5774 |

This predefined non-parametric test compares independent-sample distributions.
The unequal shapes/spreads shown in the plot prevent a pure median interpretation.
There are no pairwise tests, subgroup-specific claims, confounder adjustments,
therapeutic recommendations or clinical/subtype analyses at this stage.

Source responses, patient CSVs, original metadata, PDF/PNG, statistics, exclusions,
exact R analysis functions and session information are saved locally. A separate
stdio MCP process retrieved the real cached 763-sample result and generated the
plot. The README provides the remaining client-registration step; global Codex
configuration has not been changed.

## Later architecture and resource recommendations — review only

**Enrichr access:** discover and pin exact library names through the live catalogue.
Use downloaded GMT membership for a single-gene annotation, with no enrichment
p-value. Use the official `addList`/`enrich` API for defensible gene sets, retaining
complete service results; a custom background uses Speedrichr. Save query genes,
background, raw response, library name and retrieval date. Local over-representation
using Enrichr annotations must be distinguished from Enrichr combined-score output.
[API](https://maayanlab.cloud/Enrichr/help#api),
[R client and background support](https://github.com/wjawaid/enrichR).

**Defaults:** GO Biological Process, Molecular Function, Cellular Component and
Reactome give complementary functional and pathway views. Human KEGG and
WikiPathways can be optional additions to limit redundancy. Resolve exact
release-suffixed names from the current catalogue and freeze them per run.
Add ChEA/ENCODE only for the TF question. Do not hard-code an EMT analysis.
[Enrichr catalogue](https://maayanlab.cloud/Enrichr/datasetStatistics).

**TF identity:** use a versioned snapshot of the curated Human Transcription Factors
catalogue associated with Lambert et al., retaining identifier mapping and version.
Unmapped or unavailable coverage means unresolved, not automatically “not a TF.”
[Catalogue](https://humantfs.ccbr.utoronto.ca/),
[publication](https://pubmed.ncbi.nlm.nih.gov/29425488/).

**TF targets:** ChEA/ENCODE/ReMap provide binding-related evidence; perturbation
signatures supply expression-response evidence. Keep these distinct from
coexpression and motif predictions. Preserve study, species and cell context.
Binding near a locus does not alone prove direct functional regulation, and a
target set from another cell type is not automatically an MB program. ChEA3
documents these different evidence classes; avoid conflating them through an
integrated rank. [ChEA3 evidence description](https://maayanlab.cloud/chea3/).

**DepMap:** use local files from one pinned official release: `Model.csv`, gene
effect and optional expression. Join by ModelID and documented gene identifiers,
retain individual cell-line measurements and missing values, and record hashes.
Only report subgroup/subtype assignments with a defensible source. Never infer
classification from a familiar cell-line name. Dependency alone does not establish
a therapeutic window. [Official data portal](https://depmap.org/portal/data_page/?tab=allData).

**Structure:** keep the MCP transport thin and R functions independently callable.
Validate retrieval first, then broad/subtype comparisons, explicit subgroup-selection
rules, verified clinical modules, functional annotation, TF evidence, DepMap and
finally report orchestration. Preserve original sources separately from derived
tables and results. Revisit existing-DE imports later without forcing them into
the single-gene V1 chain. Do not implement multi-gene comparison yet.

**Reporting:** recommend Quarto for the later report: R execution via knitr,
RStudio support, and HTML/PDF from the same source with code inspectable. Start
with HTML; PDF adds document-system dependencies. No complete report template is
implemented in V1. [Quarto/RStudio guide](https://quarto.org/docs/get-started/hello/rstudio.html).

**Scientific risks:** wrong identifier joins; ambiguous reporters; double log
transformation; unverified clinical coding/time origin; invented subtype mappings;
missing values treated as a clinical class; arbitrary subgroup/cutoff selection;
small groups and non-independent samples; batch/age confounding; multiple testing;
single-gene enrichment misuse; incorrect backgrounds; conflating binding,
perturbation and prediction; and treating cell-line dependency as therapeutic
evidence. Address each before its module is enabled. No overall target score.

## Stop point

The smallest working V1 is now the three MCP tools, real patient-level retrieval,
source/metadata audit and one R plot/test. Review those outputs and test “Get R2
expression for HLX” in the connected Codex client before any later stage proceeds.
