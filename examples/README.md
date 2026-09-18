All numbers in `medulloblastoma-de.csv` and all sets in `demo-pathways.gmt`
are **synthetic software examples**, not biological results or annotations.
The DepMap files in `tests/fixtures` are synthetic too.

Run `Rscript --vanilla scripts/demo.R` to import three contrasts, rank candidates
using fixture dependencies and generate a local pathway overview without network access.
Outputs appear in `output/demo/`. Replace the input paths with real data for research.

For DESeq2 output, use `gene_column="gene"`,
`effect_column="log2FoldChange"`, `adjusted_p_column="padj"`.
For limma output, use `effect_column="logFC"`,
`adjusted_p_column="adj.P.Val"`. Add a `comparator` column to each complete
contrast table and concatenate them. Positive effects must consistently mean
target minus comparator. The server does not flip effects automatically.

For clinical analysis, an additional expression CSV/TSV has `gene,sample1,sample2,...`.
Metadata contains one row per independent sample:

```csv
sample_id,subgroup,survival_time,survival_event,metastasis,age,batch
sample1,group3,36,0,0,8,A
sample2,group4,12,1,1,10,B
```

This two-row snippet only illustrates the format, not a sufficient analysis sample.
Use consistent survival time units and document the endpoint (e.g. overall survival).
`survival_event`: 0=censored, 1=event; `metastasis`: 0=absent, 1=present.
Map stages M0/M1/M2/M3 explicitly before import and retain unknowns as missing.
The server never guesses whether a stage, relapse or death represents metastasis.
