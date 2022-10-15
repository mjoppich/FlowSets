# FlowSet

### FlowSet: an integrative approach to analyse discrete-series data using fuzzy concepts

---

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

Similar to the expression data case, we first need to prepare the differential expression data for each **TimePoint**:

    
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


The combined dataframe is then ready for usage in the FlowSet framework. [The example analysis is available here.](https://github.com/mjoppich/FlowSets/blob/main/examples/sc_ddiff.ipynb)



## Contact

[Felix Offensperger](https://github.com/offenspergerfelix)

[Markus Joppich](https://ibio.dev/) [![GitHub followers](https://img.shields.io/github/followers/mjoppich?style=social)](https://github.com/mjoppich) [![Twitter Follow](https://img.shields.io/twitter/follow/mjoppich?style=social&logo=twitter)](https://twitter.com/intent/follow?screen_name=mjoppich)
