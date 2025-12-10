# FlowSets

## Analysis and Visualization of Expression Patterns with Fuzzy Sets as FlowSets

---
![GitHub top language](https://img.shields.io/github/languages/top/mjoppich/FlowSets)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/mjoppich/FlowSets)

*FlowSets won the best poster award at ISMB/ECCB 2023 in the BioVis-Track!*

---

### Contact
- [Felix Offensperger](https://github.com/offenspergerfelix)
- [Markus Joppich](https://ibio.dev/) [![GitHub followers](https://img.shields.io/github/followers/mjoppich?style=social)](https://github.com/mjoppich) [![Twitter Follow](https://img.shields.io/twitter/follow/mjoppich?style=social&logo=twitter)](https://twitter.com/intent/follow?screen_name=mjoppich)

---

## Overview

**FlowSets** is a Python package for visualizing and analyzing gene expression patterns using fuzzy set theory. It enables the identification and visualization of gene expression flows across experimental conditions or clusters, and supports pathway enrichment analysis for genes according to a membership following specific expression patterns.

---

## Install

You can install **FlowSets** using pip:

```bash
pip install flowsets
```


## Installation Help

If you encounter issues installing polars, you may need to use the long-term support (LTS) CPU-only build:

```bash
pip uninstall polars
pip install polars-lts-cpu
```
For testing FlowSets, we recommend using a dedicated Conda environment:

```bash
conda create -n FlowSets_env python=3.12 
conda activate FlowSets_env
pip install flowsets         
```


## Quick Start Example

```python
from flowsets import *

# Read in data as polars dataframe
data = pl.read_csv(
    'small_example/deseq2_results_25deg_all_comparisons_cleaned.csv',
    null_values=['NA'],
    schema={
        "baseMean": pl.Float32,
        "log2FoldChange": pl.Float32,
        "lfcSE": pl.Float32,
        "stat": pl.Float32,
        "pvalue": pl.Float32,
        "padj": pl.Float32,
        "comparison": pl.Utf8,
        "gene_id": pl.Utf8
    }
)

# Fuzzify the log2FoldChange values for each gene and comparison
# Here all states are fuzzified with the same 
explDFWide, mfFuzzy = LegacyFuzzifier.fuzzify(
    data, #df
    stepsize=0.01,
    symbol_column="gene_id", # column name refering to feature
    meancolName="log2FoldChange", # column name refering to signal
    clusterColName="comparison", # column name refering to state
    mfLevels = ["strong_down","down","neutral","up", "strong_up"], # linguistic variables which should be created
    centers=[-2, -1, 0, 1, 2], # centers for the fuzzy sets
    sdcolName=None, exprcolName=None, # these parameters are not in use, they are meant for single cell
)
```

<img src="https://github.com/mjoppich/FlowSets/blob/main/small_example/plots/fuzzy_concept.png" width="500px" />

```python
# Create a FlowAnalysis (FlowSets) object for the fuzzified data
# The series is defined by tuples with the name in dataframe (clusterColName) and displayed name in FlowSets
def_series = (
    ("HSF1.KD vs Wildtype",'KO1 vs WT'), 
    ("Double.KDKO vs Wildtype",'KO1+2 vs WT'),
    ("MSN24.KO vs Wildtype",'KO2 vs WT')
)
fa = FlowAnalysis(explDFWide, "gene_id", def_series, mfFuzzy)

# Plot the flow memberships for all genes
fa.plot_flows(figsize=(15, 10), outfile="./small_example/plots/complete_flow.png")
```

<img src="https://github.com/mjoppich/FlowSets/blob/main/small_example/plots/complete_flow.png" width="500px" />

### Visualize only Specific Gene Sets

```python
solis_genes = ["YAL005C", "YBR101C", "YDR171W", "YDR214W", "YDR258C", "YFL016C", "YGR142W", "YLL024C", "YLL026W", "YLR216C", "YMR186W", "YNL007C", "YNL064C", "YNL281W", "YOR027W", "YOR298C-A", "YPL240C", "YPR158W"]

fa.plot_flows(genes=solis_genes, title="Solis et al. 2016 - KO1 dependent genes", figsize=(10, 8), outfile="./small_example/plots/geneset_flow.png")
```

<img src="https://github.com/mjoppich/FlowSets/blob/main/small_example/plots/geneset_flow.png" width="500px" />

### Pattern Search and Pathway Analysis

```python
# Find genes with specific flow patterns and perform pathway analysis
relFlow = fa.flow_finder(
    ["?","?"], 
    minLevels=[None,None,"down"], 
    maxLevels=["down","down","up"], 
    verbose=False
    )

fa.plot_flow_memberships(
    use_edges=relFlow, 
    color_genes=solis_genes, 
    outfile="./small_example/plots/pattern_memberships.png"
    )

pw_file = "small_example/goslim.gmt"

pwScores = fa.analyse_pathways(
    use_edges=relFlow, 
    genesets_file=pw_file, 
    additional_genesets=[("solis annotated genes", solis_genes)]
    )

pwScores_signif = pwScores.sort_values("pw_coverage_pval", ascending=True).head(20)
display(pwScores_signif)

# Show as ORA plot
fa.plotORAresult(pwScores_signif, "GOslim", numResults=10, figsize=(6,6), outfile="./small_example/plots/goslim_pathway_analysis.png")

```

<img src="https://github.com/mjoppich/FlowSets/blob/main/small_example/plots/pattern_memberships.png" width="250px" />
<img src="https://github.com/mjoppich/FlowSets/blob/main/small_example/plots/goslim_pathway_analysis.png" width="500px" />

---

## Paper Examples

- [Figure 2: Expression Analysis](https://github.com/mjoppich/FlowSets/blob/main/examples/sc_expression_data.asympt.ipynb)
- [Figure 3: Double-Differential Analysis](https://github.com/mjoppich/FlowSets/blob/main/examples/sc_ddiff.ipynb)

## other Examples
- [Basic functionalities](https://github.com/mjoppich/FlowSets/blob/main/small_example/Flowsets_example.ipynb)
- [Downstream analyses](https://github.com/mjoppich/FlowSets/blob/main/small_example/Flowsets_ORA.ipynb)
- [Custom import](https://github.com/mjoppich/FlowSets/blob/main/small_example/Flowsets_different_inputs.ipynb)

---

## Method Summary

- (Differential) Expression data are read in for each gene and each cluster (or state).
- Values are fuzzified by user-defined membership classes, min-max scaling, or quantiles.
- Relevant flows are defined using a simple grammar with `flow_finder`, specifying desired differences between levels.
- For each flow or group of flows, gene set enrichment analysis is performed. Gene sets are binned by size, and for each bin, flow memberships are calculated. A z-score is computed for each gene set (relative to others in the bin), which is transformed into a p-value for all positive-z-score (overrepresented) gene sets.

A more detailed description is available in the [working copy of our manuscript article](https://github.com/mjoppich/FlowSets/blob/main/examples/WorkingVersionFlowsets.pdf).

<!---
---



## Extracting Data from Seurat Objects

We provide a collection of useful R functions for extracting expression and differential data from Seurat objects. To access these methods, run:

```r
source("https://raw.githubusercontent.com/mjoppich/FlowSets/main/seurat_util_functions.R")
```


### Gene Expression Data

Calculate gene expression data and group by **TimePoint**:

```r
df.all = getExtendedExpressionData(obj, assay="RNA", group.by="TimePoint")
write.table(df.all, "expression_all.tsv", quote=F, sep="\t", row.names=F)
```

The data frame can then be read in Python and used for analysis. [Example analysis here.](https://github.com/mjoppich/FlowSets/blob/main/examples/sc_expression_data.all.ipynb)

### Differential Data

Prepare differential expression data for each **TimePoint**:

```r
celltype = "Monocytes-Immune-system"
print(celltype)

cells.celltype = cellIDForClusters(obj, "cellnamesread", c(celltype))
ignoreAnalysis = FALSE
tpDeList = list()
for (timep in c("1", "2", "3")) {
    print(timep)
    cells.timepoint = cellIDForClusters(obj, "TimePoint", c(timep))
    cells.comp.sympt = intersect(cells.sympt, intersect(cells.timepoint, cells.celltype))
    cells.comp.asympt = intersect(cells.asympt, intersect(cells.timepoint, cells.celltype))
    print(paste(length(cells.comp.sympt), length(cells.comp.asympt)))
    if ((length(cells.comp.sympt) < 3) || (length(cells.comp.asympt) < 3)) {
        ignoreAnalysis = TRUE
        next()
    }
    deResult = compareClusters(scdata=obj,
        cellsID1=cells.comp.sympt,
        cellsID2=cells.comp.asympt,
        prefix= paste("cluster", celltype, timep, sep="_"),
        suffix1="cells_sympt",
        suffix2="cells_asympt",
        test="t", fcCutoff=0.25, assay="RNA", outfolder=paste("./de_comparison_", celltype, sep=""))
    tpDeList[[timep]] = deResult
}
if (ignoreAnalysis) {
    print(paste("Skipping", celltype))
    next()
}
makeCombinedDF(tpDeList, paste("./de_comparison_", celltype, sep=""))
```
-->


---

## License

This project is licensed under the MIT License.

---

## Citation

If you use FlowSets in your research, please cite our manuscript (see [WorkingVersionFlowsets.pdf](https://github.com/mjoppich/FlowSets/blob/main/examples/WorkingVersionFlowsets.pdf)).
