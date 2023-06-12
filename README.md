# FlowSets

## FlowSets: an integrative approach to analyse discrete-series data using fuzzy concepts

---
![GitHub top language](https://img.shields.io/github/languages/top/mjoppich/FlowSets)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/mjoppich/FlowSets)

#### Contact

[Felix Offensperger](https://github.com/offenspergerfelix)

[Markus Joppich](https://ibio.dev/) [![GitHub followers](https://img.shields.io/github/followers/mjoppich?style=social)](https://github.com/mjoppich) [![Twitter Follow](https://img.shields.io/twitter/follow/mjoppich?style=social&logo=twitter)](https://twitter.com/intent/follow?screen_name=mjoppich)

*FlowSets will be presented at ISMB/ECCB 2023 in the BioVis-Track!*

---


<img src="https://github.com/mjoppich/FlowSets/blob/main/examples/ddiff_flows.png" width="500px" />

Check-out the expression example [here](https://github.com/mjoppich/FlowSets/blob/main/examples/sc_expression_data.all.ipynb) or the double-differential analysis [here](https://github.com/mjoppich/FlowSets/blob/main/examples/sc_ddiff.ipynb).

## Examples


## Extracting data from Seurat objects

We provide a collection of useful R-functions for extracting expression and differential data from Seurat objects. To access these methods you need to

    source("https://raw.githubusercontent.com/mjoppich/FlowSets/main/seurat_util_functions.R")



### Gene expression data

We first need to calculate gene expression data, and group them by **TimePoint**

    df.all = getExtendedExpressionData(obj, assay="RNA", group.by="TimePoint")

It is then possible to write the data frame to disk

    write.table(df.all, "expression_all.tsv", quote=F, sep="\t", row.names=F)

The data frame can then be read in python and used for analysis. [The example analysis is available here.](https://github.com/mjoppich/FlowSets/blob/main/examples/sc_expression_data.all.ipynb)



### Differential data

Similar to the expression data case, we first need to prepare the differential expression data for each **TimePoint**

    
    celltype = "Monocytes-Immune-system"
    print(celltype)
    
    cells.celltype = cellIDForClusters(obj, "cellnamesread", c(celltype))
    ignoreAnalysis = FALSE
    tpDeList = list()
    for (timep in c("1", "2", "3"))
    {
        
        print(timep)
        cells.timepoint = cellIDForClusters(obj, "TimePoint", c(timep))
        
        cells.comp.sympt = intersect(cells.sympt, intersect(cells.timepoint, cells.celltype))
        cells.comp.asympt = intersect(cells.asympt, intersect(cells.timepoint, cells.celltype))
        
        print(paste(length(cells.comp.sympt), length(cells.comp.asympt)))
        
        if ((length(cells.comp.sympt) < 3) || (length(cells.comp.asympt) < 3))
        {
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
    
    if (ignoreAnalysis)
    {
        print(paste("Skipping", celltype))
        next()
    }
    
    makeCombinedDF(tpDeList, paste("./de_comparison_", celltype, sep=""))


The combined dataframe is then ready for usage in the FlowSets framework. [The example analysis is available here.](https://github.com/mjoppich/FlowSets/blob/main/examples/sc_ddiff.ipynb)

### Brief Method description

(Differential) Expression data are read in for each gene and each cluster (or: state). The values are fuzzified either by user-defined membership classes, or equally distributed over the measurement range (min-max), or according to predefined quantiles.

Relevant flows can be defined using a simple grammar with the flow_finder function, where the desired difference between two levels can be specified.

For each flow, or a group of flows, gene set enrichment analysis can be performed. Here, the gene sets are binned according to their size. E.g. all gene sets with at least 2 and at most 5 genes are put together into one bin. For each bin, all flow memberships are calculated. For each membership a z-score is calculated (how different is a geneset from all other gene sets of that bin), which is transformed into a p-value for all positive-z-score (=more than expected) gene sets.

