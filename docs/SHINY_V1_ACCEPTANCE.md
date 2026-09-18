# Shiny Version 1 acceptance - 18 September 2026

Implemented `shiny/app.R` and `shiny/interface.R` using the existing shared
`get_r2_expression()` and `analyze_r2_subgroups()` functions. Scientific methods
and core plot construction were not changed. No MCP process or AI is involved.

## Real-data checks

| Gene | Retrieved | Analyzed | Kruskal-Wallis H | Raw p |
| --- | ---: | ---: | ---: | ---: |
| HLX | 763 | 763 | 193.06166 | 1.331110e-41 |
| SLC2A1 | 763 | 763 | 96.67672 | 8.052715e-21 |

Both used the existing validated R2 cache, not a new network retrieval. Both
were entered and run in the local browser; statistics, group summaries and
patient-level ggplots rendered. Changing HLX to SLC2A1 without clicking left
the HLX result intact. Clicking Run analysis showed the busy/disabled button
and progress notification, then displayed SLC2A1 results.

Browser downloads were exercised for patient CSV (HLX) and statistics CSV,
PDF and PNG (SLC2A1). The acceptance script independently exported all four
formats for both genes, checked nonempty files, verified patient counts, and
compared exported p-values with the core result. Artifacts are in
`outputs/shiny-v1-acceptance/`; rerun with:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --vanilla scripts/check-shiny-v1.R
```

## Automated checks

The full `tests/run-tests.R` suite passed, including the new Shiny server tests:

- No analysis at startup or on text edits; analysis on button click.
- Core analysis called with the retrieved patient object.
- CSV/PDF/PNG outputs use stored results.
- Friendly invalid-gene, retrieval-failure and missing-subgroup messages.
- Failed runs clear the old result; subsequent valid runs recover.

The optional live-network test was skipped. Three existing serialization
warnings in core-report tests and Windows locale startup warnings remain;
they did not fail tests. Browser tests do not assert live R2 availability.

## Issues resolved during implementation

- Shiny was absent: installed it and its dependencies into `.R-library`.
- `R/r2_client.R` and `R/depmap.R` had been emptied independently of this work.
  With user approval, restored their functions from the existing SLC2A1 report
  snapshot. This resolved six failing pre-existing tests. Newly created empty
  module files were preserved.
- Corrected Windows source-text encoding and provided horizontal scrolling
  for the unchanged figure at narrow browser widths.

## Scope limit

Overview, MB Subgroups and Downloads are working. MB Subtypes, Survival,
Metastasis, Functional Biology, TF Targets and DepMap show "Not implemented yet".
No calculations for those future modules are triggered. The prototype remains
synchronous; a slow R2 request occupies its R process. No public deployment
or multi-user scaling was attempted.
