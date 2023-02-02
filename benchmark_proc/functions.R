

library(Tempora)
library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)
library(dplyr)
library(tidyverse)
library(pbapply)
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(writexl)
library(EnhancedVolcano)
library(ggsunburst)
library(ggrepel)
library(RColorBrewer)
library(svglite)
library(ggpubr)
library(gtools)
library("purrr")
library(rstatix)
library(ComplexHeatmap)
library(circlize)
library(pbapply)
library(future)


#install.packages(c("Seurat", "devtools", "stringr", "ggplot2", "pbapply", "tidyverse", "BiocManager", "svglite", "writexl", "data.table", "ggrepel", "pbapply", "cowplot", "rstatix", "ggpubr"))
# BiocManager::install(c("EnhancedVolcano", "enrichplot", "clusterProfiler", "GSEABase", "GSVA", "ComplexHeatmap", "circlize", "biomaRt", "limma"))
# devtools::install_github("BaderLab/Tempora")
# devtools::install_github("didacs/ggsunburst")
# remotes::install_version("randomForest"), version="4.7-1")
# devtools::install_github("saeyslab/nichenetr")
#
##
### Global Helpers
##
#

patternList.human = list()
patternList.human[["MT"]] = "^MT-"
patternList.human[["RPL"]] = "^RPL"
patternList.human[["RPS"]] = "^RPS"


patternList.mouse = list()
patternList.mouse[["MT"]] = "^mt-"
patternList.mouse[["RPL"]] = "^Rpl"
patternList.mouse[["RPS"]] = "^Rps"



#
##
### Seurat Read-In
##
#


readMtxFiles = function(files, samplenamepos=3)
{

  allfiles.raw = list()
  allABs.raw = list()

  for (file in files)
  {
    samplename = str_split(dirname(file), "/")[[1]][samplenamepos]
    foldername = dirname(file)
    
    print(paste(samplename, foldername))
    
    h5file = Read10X(foldername,unique.features = TRUE)

    if (is.null(names(h5file)))
    {
      print(paste("WITHOUT AB", samplename))
      allfiles.raw[[samplename]] = h5file
    } else {
      print(paste("WITH AB", samplename))
      allfiles.raw[[samplename]] = h5file$`Gene Expression`
      allABs.raw[[samplename]] = h5file$`Antibody Capture`
    }

    print(paste(samplename,nrow(allfiles.raw[[samplename]]), ncol(allfiles.raw[[samplename]]), "genes x cells"))
  }

  return(list(gex=allfiles.raw, ab=allABs.raw))
}

toObjList = function(inputMatrices, patternlist, variable.features=3000)
{

objlist = list()

for (x in names(inputMatrices$gex))
{

    matrix = inputMatrices$gex[[x]]
    
    filteredObj = makeSeuratObj(matrix, x, patternlist)   
    
    filteredObj <- NormalizeData(filteredObj, verbose = FALSE)
    filteredObj <- FindVariableFeatures(filteredObj, nfeatures=variable.features, verbose = FALSE)
    
    objlist[[x]] = filteredObj

    print(x)
    print(filteredObj)
        
}

return(objlist)

}


makeSeuratObj = function(matrix, proj, pl)
{
    obj = CreateSeuratObject(matrix, project=proj)
    print("Renaming Cells")
    obj <- RenameCells(obj, add.cell.id=proj)
    
    print(paste("Seurat obj project", obj@project.name))
    
    mtPattern = pl[["MT"]]
    rplPattern = pl[["RPL"]]
    rpsPattern = pl[["RPS"]]
    rpPattern = paste(c(pl[["RPL"]], pl[["RPS"]]), sep="", collapse="|")

    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=mtPattern)]
    print(paste("Got a total of mt-Genes:", length(selGenes), paste(head(selGenes), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rplPattern)]
    print(paste("Got a total of Rpl-Genes:", length(selGenes), paste(head(selGenes), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rpsPattern)]
    print(paste("Got a total of Rps-Genes:", length(selGenes), paste(head(selGenes), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rpPattern)]
    print(paste("Got a total of Rp-Genes:", length(selGenes), paste(head(selGenes), collapse =", ")))
    
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mtPattern)
    obj[["percent.rpl"]] <- PercentageFeatureSet(obj, pattern = rplPattern)
    obj[["percent.rps"]] <- PercentageFeatureSet(obj, pattern = rpsPattern)
    obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = rpPattern)
    
    return(obj)
}

scatterAndFilter = function(objlist, nfeature_rna.lower=100, nfeature_rna.upper=6000, ncount_rna.lower=500, percent_mt.upper=7)
{

  for (name in names(objlist))
  {
    print(name)

    plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(plot1 + plot2, paste(name, "scatter_ncount_mt", sep="_"), fig.width=10, fig.height=6)

    plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
    plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(plot1 + plot2, paste(name, "scatter_ncount_rp", sep="_"), fig.width=10, fig.height=6)
  }


  objlist.new <- lapply(X = objlist, FUN = function(obj) {
    # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
    
    obj <- subset(obj, subset = nFeature_RNA > nfeature_rna.lower & nFeature_RNA < nfeature_rna.upper & nCount_RNA > ncount_rna.lower)
    obj <- subset(obj, subset = percent.mt < percent_mt.upper)
    print(obj)
    
    return(obj)
  })


  for (name in names(objlist.new))
  {
    p=VlnPlot(objlist.new[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
    save_plot(p, paste(name, "filtered_violins_qc", sep="_"), fig.width=10, fig.height=6)
    
    plot1 <- FeatureScatter(objlist.new[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(objlist.new[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(plot1 + plot2, paste(name, "filtered_scatter_ncount_mt", sep="_"), fig.width=10, fig.height=6)
    
    plot1 <- FeatureScatter(objlist.new[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
    plot2 <- FeatureScatter(objlist.new[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(plot1 + plot2, paste(name, "filtered_scatter_ncount_rp", sep="_"), fig.width=10, fig.height=6)
  }

  return(objlist.new)

}

#
##
### Hashtag Oligos
##
#

processHTO = function(inputMatrices, objlist, relevantHTOs)
{

htoObjList = list()
for (name in intersect(names(inputMatrices$ab), names(relevantHTOs)))
{

    relHTOs = relevantHTOs[[name]]
    print("Relevant HTOs:")
    print(relHTOs)

    htoMatrix = inputMatrices$ab[[name]]

    cellsGEX = length(colnames(objlist[[name]]))
    cellsAB = length(colnames(htoMatrix))

    gexnames = substring(colnames(objlist[[name]]), str_length(name)+2)
    print(head(gexnames))
    cellsJoint = intersect(gexnames, colnames(htoMatrix))
    cellsABsub = htoMatrix[, cellsJoint]
    cellsABsub = cellsABsub[relHTOs,]

    colnames(cellsABsub) = paste(objlist[[name]]@project.name, colnames(cellsABsub), sep="_")
    print(head(colnames(cellsABsub)))
    print(paste(cellsGEX, cellsAB, length(cellsJoint), ncol(cellsABsub)))

    # Normalize RNA data with log normalization
    xobj <- NormalizeData(objlist[[name]])
    # Find and scale variable features
    xobj <- FindVariableFeatures(xobj, selection.method = "mean.var.plot")
    xobj <- ScaleData(xobj, features = VariableFeatures(xobj))

    xobj[["HTO"]] <- CreateAssayObject(counts = cellsABsub)
    xobj <- NormalizeData(xobj, assay = "HTO", normalization.method = "CLR")
    xobj <- HTODemux(xobj, assay = "HTO", positive.quantile = 0.99)

    print(table(xobj$HTO_classification.global))
    print(table(xobj$HTO_classification))

    htoObjList[[name]] = xobj

}

return(htoObjList)

}

qcFilterHTO = function(htoObjList)
{

  for (name in names(htoObjList))
  {
    p=VlnPlot(htoObjList[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "HTO_classification.global")
    save_plot(p, paste(name, "hto_violins_qc", sep="_"), fig.width=10, fig.height=6)


    allHTOFeatures = rownames(htoObjList[[name]][["HTO"]])
    figHeight = round(length(allHTOFeatures)*0.5 * 6)
    r=RidgePlot(htoObjList[[name]], assay = "HTO", features = allHTOFeatures, ncol = 2)
    save_plot(r, paste(name, "hto_ridge_qc", sep="_"), fig.width=10, fig.height=figHeight)

    allHTOFeatures = rownames(htoObjList[[name]][["HTO"]])
    figHeight = round(length(unique(htoObjList[[name]]$HTO_classification))*0.5 * 4)
    r=RidgePlot(htoObjList[[name]], assay = "HTO", features = allHTOFeatures, group.by="HTO_classification", ncol = 2)
    save_plot(r, paste(name, "hto_ridge_detail_qc", sep="_"), fig.width=10, fig.height=figHeight)

  
    print(name)
    print(table(htoObjList[[name]]$HTO_classification))

  }

  htoFiltered <- lapply(X = htoObjList, FUN = function(obj) {
    
    print(obj)
    print(table(obj$HTO_classification))
    obj = subset(obj, subset=HTO_classification.global == "Singlet")
    print(obj)
    print(table(obj$HTO_classification)) 
    return(obj)
  })

  return(htoFiltered)

}

mergeRNAHTO = function(objlist, htolist)
{

  returnList = list()

  for (name in names(objlist))
  {
    if (name %in% names(htolist))
    {
      returnList[[name]] = htolist[[name]]
    } else {
      returnList[[name]] = objlist[[name]]
    }
  }
  return(returnList)
}

splitObjListByGroup = function(objlist, group.by)
{
  finalList = list()
  for (objname in names(objlist))
  {

    if (group.by %in% colnames(objlist[[objname]]@meta.data))
    {
        for (abTag in unique(objlist[[objname]][[group.by]][[group.by]]))
        {
            sampleName = paste(objname, abTag, sep="_")
            print(sampleName)

            selDF = objlist[[objname]][[group.by]]
            selCells = rownames(selDF)[selDF[[group.by]] == abTag]

            finalList[[sampleName]] = subset(objlist[[objname]], cells = selCells)
            finalList[[sampleName]]$orig_project = objname
            finalList[[sampleName]]$library = paste(objname, abTag, sep="_")
            
            print(sampleName)
            print(finalList[[sampleName]])
        }

    } else {
        finalList[[objname]] = objlist[[objname]]

        print(objname)
        print(finalList[[objname]])
    }

  }

  return(finalList)
}



#' Takes a list of antibody capture feature-cell-matrices, the list of Seurat objects and a list of relevant hashtag oligo IDs for each sample
#'
#'
#' @param imats list of antibody capture feature-cell-matrices
#' @param objlist list of Seurat objects
#' @param HTOs hashtag oligo IDs to ignore for processing CITE
#' @param assayName name of the newly created CITE-seq assay
#'
#' @return list of Seurat objects with HTO assay
#'
#'
#' @export
processCITE = function(objlist, imats, assayName="ADT", HTOs=NULL, run.parallel=FALSE, relevantCITEs=NULL, rownametransform=NULL)
{

    retlist = list()
    for (name in names(objlist))
    {
      print(name)
        obj.in = objlist[[name]]

        if (! (name %in% names(imats$ab)))
        {
          print(paste("Not processing", name))
          retlist[[name]] = obj.in
          next
        }
        ab.raw = imats$ab[[name]]

        #print(head(ab.raw))
        colnames(ab.raw) = paste(name, colnames(ab.raw), sep="_")

        #remove HTOs from ab.raw

        if (!is.null(HTOs))
        {
            ab.raw = ab.raw[ !(rownames(ab.raw) %in% HTOs), ]
            print(rownames(ab.raw))
        }

        if (!is.null(relevantCITEs))
        {
          if (!name %in% names(relevantCITEs))
          {
            print(paste("Skipping", name))
            retlist[[name]] = obj.in
            next
          }

          ab.raw = ab.raw[relevantCITEs[[name]],]
        }

        if (!is.null(rownametransform))
        {

          commonRows = intersect( names(rownametransform), rownames(ab.raw) )
          objRownameTransform = rownametransform[commonRows]

          ab.raw = ab.raw[names(objRownameTransform), ]
          rownames(ab.raw) = as.character(objRownameTransform)

          print("RownameTransform")
          print(rownames(ab.raw))
        }

        #print(head(ab.raw))

        ab.raw = ab.raw[, colnames(obj.in)]
        
        adt_assay <- Seurat::CreateAssayObject(counts = ab.raw)
        obj.in[[assayName]] = adt_assay


        DefaultAssay(obj.in) <- assayName
        if (!run.parallel)
        {
        t = plan()
        plan("sequential")
        }
        obj.in <- Seurat::NormalizeData(obj.in, normalization.method = "CLR", margin = 2)
        if (!run.parallel)
        {
        plan(t)
        }

        Seurat::DefaultAssay(obj.in) <- "RNA"
        retlist[[name]] = obj.in
    }

    return(retlist)
}



#
##
###
#### Seurat object integration
###
##
#



prepareIntegration = function(finalList, cc.use.genes, nfeatures.variable = 3000, nfeatures.scale=3000, scale.regress=c('percent.rp', 'percent.mt', "nCount_RNA","S.Score", "G2M.Score"), normalize=TRUE, findvariable=TRUE, run.parallel=TRUE, npcs=50)
{
 

    print("cells per experiment")
    print(mapply(sum, lapply(finalList, function(x) {dim(x)[2]})))
    print("total cells")
    print(sum(mapply(sum, lapply(finalList, function(x) {dim(x)[2]}))))


    objlist = list()
    for (objname in names(finalList))
    {
        x = finalList[[objname]]

        if (! "orig_project" %in% colnames(x@meta.data))
        {
          x$orig_project = objname
        }

        Project(x) = objname
        print(paste("Seurat obj project", x@project.name))

        DefaultAssay(x) = "RNA"

        if (normalize)
        {
          x <- NormalizeData(x, verbose = FALSE)
        }
        
        if (findvariable)
        {
          x <- FindVariableFeatures(x, nfeatures=nfeatures.variable, verbose = FALSE)
        }
      
        x$library = objname

        objlist[[objname]] = x
    }

    print("SelectIntegrationFeatures")
    features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = nfeatures.scale)

    objlist <- lapply(X = objlist, FUN = function(x) {

        print(paste("Seurat obj project", x@project.name))
        print(x)

        s.genes <- cc.use.genes$s.genes
        g2m.genes <- cc.use.genes$g2m.genes

        s.genes = intersect(s.genes, rownames(x))
        g2m.genes = intersect(g2m.genes, rownames(x))

        print(paste("CellCycle", length(s.genes), length(g2m.genes)))

        x <- CellCycleScoring(
        x,
        g2m.features = g2m.genes,
        s.features = s.genes)

        if (!run.parallel)
        {
        t = plan()
        plan("sequential")
        }

        x <- ScaleData(x, features = features, verbose = TRUE, assay="RNA", vars.to.regress = scale.regress)
        x <- RunPCA(x, npcs=npcs, verbose = FALSE, reduction.name="pca", assay="RNA")

        if (!run.parallel)
        {
        plan(t)
        }

        x$project = x@project.name

        return(x)
    })

    return(list("data"=objlist, "features"=features))
}

performIntegration = function(objlist, intname, features.integration = 3000, 
 gex.method.normalization="LogNormalize", gex.assay="RNA", gex.runpca=TRUE,
gex.dims=30, gex.method.integration="rpca", gex.k.filter=200, gex.k.anchor = 5,gex.k.weight=100, 
add.do=FALSE,add.assay="ABA",
add.dims=10, add.method.integration="rpca", add.k.filter=200, add.k.anchor = 5, add.features.integration=100, add.k.weight=5,
run.parallel=TRUE, fig.width=8, fig.height=6, add.cluster.resolution=1) 
{

    dir.create(intname, recursive = TRUE)

    if (add.do)
    {
      #
      # integrate based on ADDitional assay
      #
      objSamples = objlist

      objSamples = lapply(objSamples, function(x) {
        print(paste("Object", x@project.name))
        DefaultAssay(x) <- add.assay
        VariableFeatures(x) <- rownames(x) # all HTOs

        if (dim(x@assays[[add.assay]]@scale.data)[1] == 0)
        {
          if (!run.parallel)
          {
            t = plan()
            plan("sequential")
          }
          x = ScaleData(x, assay=add.assay)
          if (!run.parallel)
          {
            plan(t)
          }

        }
        
        if (("pca" %in% names(x@reductions)) && (x@reductions$pca@assay.used == add.assay))
        {
          print("PCA ALREADY THERE")
        } else {
          x = RunPCA(x,features = rownames(x),verbose = FALSE, reduction.name="pca", approx=FALSE, npcs=add.dims, assay=add.assay)
        }

        return(x)
      })

      print(objSamples)

      features_add <- rownames(objSamples[[1]][[add.assay]])

      if (!run.parallel)
      {
        t = plan()
        plan("sequential")
      }
      assayData = rep(add.assay, length(objSamples))
      objlist.anchors.add <- FindIntegrationAnchors(object.list = objSamples, assay=assayData, normalization.method = "LogNormalize",
                                                    anchor.features = add.features.integration, dims = 1:add.dims, reduction = add.method.integration, k.filter = add.k.filter, k.anchor=add.k.anchor)

      add.list.integrated <- IntegrateData(new.assay.name = "integrated_add", anchorset = objlist.anchors.add, normalization.method = "LogNormalize", dims=1:(add.dims-1), k.weight=add.k.weight)

      if (!run.parallel)
      {
        plan(t)
      }


      print("ADD integration done")
    }


    #
    # integrate based on RNA/GEX assay
    #
    objSamples = objlist
    print(objSamples)

    print("GEX integration features")
    print(objSamples)
    print(paste("Current integration mode:", gex.method.integration))

    objSamples = lapply(objSamples, function(x) {
        print(paste("Object", x@project.name))
        DefaultAssay(x) <- gex.assay
        return(x)
    })


    if (gex.method.integration!="SCT")
    {


      if (gex.method.normalization == "SCT")
      {
        print("SCTransform")
        objSamples <- lapply(X = objSamples, FUN = SCTransform, method = "glmGamPoi")

      }

      print("SelectIntegrationFeatures")
      if (is.numeric(features.integration))
      {
        print("RunPCA")
        objSamples <- lapply(X = objSamples, FUN = RunPCA, npcs=min(c(50, gex.dims)), verbose = FALSE, reduction.name="pca", assay=gex.assay)
        print("Select Integration Features")
        features_gex <- SelectIntegrationFeatures(object.list = objSamples, nfeatures = features.integration, assay=rep(gex.assay, length(objSamples)))
      } else {
        features_gex = features.integration
        print("RunPCA on given features")

        objSamples <- lapply(X = objSamples, FUN = function(x)
        {
          if (gex.runpca)
          {
            x = RunPCA(x,npcs=min(c(50, gex.dims)),verbose = FALSE, reduction.name="pca", features=features_gex, assay=gex.assay)
          }

          return(x)
      })

      }


      if (gex.method.normalization == "SCT")
      {
        print("PrepSCTIntegration")
        objSamples <- PrepSCTIntegration(object.list = objSamples, anchor.features = features_gex)
        print("Calculating PCAs on SCT")
        objSamples <- lapply(X = objSamples, FUN = RunPCA, features = features_gex)
      }

      print("FindIntegrationAnchors")

      if (!run.parallel)
      {
        t = plan()
        plan("sequential")
      }
      print(plan())

      objlist.anchors <- FindIntegrationAnchors(object.list = objSamples,  reduction = gex.method.integration, dims = 1:gex.dims, anchor.features = features_gex, normalization.method=gex.method.normalization, k.anchor=gex.k.anchor, k.filter = gex.k.filter)
      obj.list.integrated <- IntegrateData(new.assay.name = "integrated_gex", anchorset = objlist.anchors, dims = 1:gex.dims, verbose=T, normalization.method = gex.method.normalization, k.weight=gex.k.weight)
      
      if (!run.parallel)
      {
        plan(t)
      }


      print("IntegrateData")


    } else {

      objSamples = lapply(objSamples, function(x) {
          DefaultAssay(x) <- gex.assay

          x <- RunPCA(x, npcs=max(c(50, gex.dims)), verbose = FALSE, reduction.name="pca",  assay=gex.assay)
          suppressWarnings(x <- SCTransform(x,vars.to.regress = c('percent.mt', 'percent.rp'), verbose = T))

          return(x)
      })

      if (!run.parallel)
      {
        t = plan()
        plan("sequential")
      }

      features_gex <- SelectIntegrationFeatures(object.list = objSamples, nfeatures = features.integration)#, assay=rep("RNA", length(objSamples)))
      objSamples <- PrepSCTIntegration(object.list = objSamples, anchor.features = features_gex)

      objlist.anchors <- FindIntegrationAnchors(object.list = objSamples, normalization.method = "SCT", anchor.features = features_gex, k.anchor=gex.k.anchor,k.filter = add.k.filter)
      obj.list.integrated <- IntegrateData(anchorset = objlist.anchors, normalization.method = "SCT", new.assay.name = "integrated_gex",verbose=T, k.weight=gex.k.weight)

      if (!run.parallel)
      {
        plan(t)
      }

    }
    print("GEX integration done")

    #
    # integrated GEX viz
    #
    if (gex.method.normalization != "SCT")
    {
      if (!run.parallel)
      {
        t = plan()
        plan("sequential")
      }
      obj.list.integrated = ScaleData(obj.list.integrated, assay="integrated_gex")
    
      if (!run.parallel)
      {
        plan(t)
      }
    }
    
    obj.list.integrated <- RunPCA(obj.list.integrated, npcs = gex.dims, reduction.name="igpca", assay="integrated_gex")
    obj.list.integrated <- RunUMAP(obj.list.integrated, reduction = "igpca", dims = 1:gex.dims, reduction.name="ig.umap", reduction.key = "UMAPig_",)
    p=DimPlot(obj.list.integrated, group.by="orig_project", reduction="ig.umap", shuffle = TRUE, seed = 1)
    save_plot(p, paste(intname, "ig_dimplot", sep="/"), fig.width, fig.height)

    obj.list.gex_add = NULL

    if (add.do)
    {
      #
      # integrated ADT viz
      #
      obj.list.integrated[["integrated_add"]] = add.list.integrated[["integrated_add"]]

      if (!run.parallel)
      {
        t = plan()
        plan("sequential")
      }

      obj.list.integrated = ScaleData(obj.list.integrated, assay="integrated_add")

      if (!run.parallel)
      {
        plan(t)
      }

      obj.list.integrated <- RunPCA(obj.list.integrated, features = rownames(add.list.integrated[[add.assay]]), verbose = FALSE, approx=FALSE, npcs=add.dims, reduction.name="iapca", assay="integrated_add")
      obj.list.integrated <- RunUMAP(obj.list.integrated, reduction = "iapca", dims = 1:add.dims, reduction.name="ia.umap", reduction.key = "UMAPia_",)

      p=DimPlot(obj.list.integrated, group.by="orig_project", reduction="ia.umap", shuffle = TRUE, seed = 1)
      save_plot(p, paste(intname, "ia_dimplot", sep="/"), fig.width, fig.height)


      #
      # multi modal neighbors
      #

      obj.list.gex_add <- FindMultiModalNeighbors(
        obj.list.integrated, reduction.list = list("igpca", "iapca"), 
        dims.list = list(1:gex.dims, 1:add.dims), prune.SNN=1/20
      )
      #
      # multi modal viz
      #
      obj.list.gex_add <- RunUMAP(obj.list.gex_add, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
      obj.list.gex_add <- FindClusters(obj.list.gex_add, graph.name = "wsnn", algorithm = 3, resolution = add.cluster.resolution, verbose = FALSE)

      p <- DimPlot(obj.list.gex_add, reduction = 'wnn.umap', label = TRUE, repel = FALSE, label.size = 2.5)
      save_plot(p, paste(intname, "wnn_cluster_dimplot", sep="/"), fig.width, fig.height)

      p <- DimPlot(obj.list.gex_add, group.by="orig_project", reduction = 'wnn.umap', label = TRUE, repel = FALSE, label.size = 2.5)
      save_plot(p, paste(intname, "wnn_project_dimplot", sep="/"), fig.width, fig.height)


      p <- DimPlot(obj.list.gex_add, reduction = 'ig.umap', label = TRUE, repel = FALSE, label.size = 2.5)
      save_plot(p, paste(intname, "wnn_cluster_ig_dimplot", sep="/"), fig.width, fig.height)


      p <- DimPlot(obj.list.gex_add, reduction = 'ia.umap', label = TRUE, repel = FALSE, label.size = 2.5)
      save_plot(p, paste(intname, "wnn_cluster_ia_dimplot", sep="/"), fig.width, fig.height)


      obj.list.integrated$wnn_clusters = Idents(obj.list.gex_add)
    }

    return(list("integrated"=obj.list.integrated, "multimodal"=obj.list.gex_add))

}

#' Preprocesses an integrated Seurat object for downstream analysis
#'
#' @param obj.in Seurat object
#' @param useAssay which assay to use
#' @param inname output folder of all plots
#' @param do.scale whether to scale the object
#' @param num.pcs number of PCs to calculate and use for UMAP/neighbors
#' @param resolution resolution for the clustering
#' @param plot.reduction which reduction to use for plotting
#' @param dim.reduction name of the reduction in which the PCA will be stored, and which is used for UMAP and Neighbors
#' @param with.hto whether also HTO plots should be prepared
#' @param run.parallel whether the ScaleData function should run in parallel or sequential
#' @param run.umap_neighbors whether to run RunUMAP, FindNeighbors and FindClusters
#'
#' @return preprocessed Seurat object
#' @export
#'
preprocessIntegrated = function(obj.in, useAssay, inname, do.scale=T, num.pcs=50, resolution=0.5, plot.reduction="umap", dim.reduction="umap", with.hto=TRUE, run.parallel=TRUE, run.umap_neighbors=TRUE, clusters.graph.name=NULL)
{
  
  if (!dir.exists(inname))
  {
    print(paste("Creating DIR", inname))
    dir.create(inname, recursive = TRUE)
  }


  Seurat::DefaultAssay(obj.in) <- useAssay

  # Run the standard workflow for visualization and clustering
  if (do.scale)
  {
    # Scale is not required for SCT!
    print("Scale Data")

    if (!run.parallel)
    {
      t = future::plan()
      future::plan("sequential")
    }
    obj.in <- Seurat::ScaleData(obj.in, verbose = FALSE)

    if (!run.parallel)
    {
      future::plan(t)
    }
  }
    
  if ((is.null(dim.reduction)) || (!dim.reduction %in% names(obj.in@reductions)))
  {
    print("RunPCA Data")
    dim.reduction = "pca"

    obj.in <- Seurat::RunPCA(obj.in, npcs = max(c(num.pcs, 50)), verbose = FALSE, reduction.name=dim.reduction)
    
  }

  print(paste("dim.reduction", dim.reduction))

  if ("pca" %in% dim.reduction)
  {
    p=Seurat::ElbowPlot(obj.in, ndims=30, reduction = dim.reduction)
    save_plot(p, paste(inname, "elbowplot", sep="/"), 12, 6)
  }
  

  if (run.umap_neighbors == TRUE)
  {
    print("RunUMAP Data")
    obj.in <- Seurat::RunUMAP(obj.in, reduction = dim.reduction, dims = 1:num.pcs)
    print("FindNeighbors Data")
    obj.in <- Seurat::FindNeighbors(obj.in, reduction = dim.reduction, dims = 1:num.pcs)
  }
  print("FindClusters Data")
  obj.in <- Seurat::FindClusters(obj.in, resolution = resolution, graph.name=clusters.graph.name)

  obj.in$idents = Seurat::Idents(obj.in)

  p=Seurat::DimPlot(obj.in, pt.size = 0.001, label=T, reduction = plot.reduction)
  save_plot(p, paste(inname, "dimplot_umap", sep="/"), fig.width=12, fig.height=8)

  numProjects = length(unique(obj.in$orig_project))
  numRows = ceiling(numProjects/2)

  p=Seurat::DimPlot(obj.in, pt.size = 0.001, label=T, split.by="orig_project", reduction = plot.reduction, ncol=2)
  save_plot(p, paste(inname, "dimplot_umap_project", sep="/"), fig.width=24, fig.height=8*numRows)
  

  if (with.hto)
  {
    tryCatch({

      obj.in$libraryHTO = paste(obj.in$library, obj.in$HTO_classification, sep="_")

      numLibHTO = length(unique(obj.in$libraryHTO))
      numRows = ceiling(numLibHTO/3)

      p=Seurat::DimPlot(obj.in, split.by="libraryHTO", pt.size = 0.001, label=T, reduction = plot.reduction, ncol=3)
      save_plot(p, paste(inname, "dimplot_umap_libraryHTO", sep="/"), fig.width=16, fig.height=min(3*numRows, 60))
    })
  }


  return(obj.in)

}


makeQCPlots = function(inobj, outfolder, reduction="umap")
{
  
  p=FeaturePlot(inobj, "nCount_RNA", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_ncount_rna", sep="/"), fig.width=8, fig.height=6)

  p=FeaturePlot(inobj, "nFeature_RNA", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_nfeature_rna", sep="/"), fig.width=8, fig.height=6)

  inobj$log_ncount = log(inobj$nCount_RNA)
  inobj$log_nfeature = log(inobj$nFeature_RNA)

  p=FeaturePlot(inobj, "log_ncount", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_logncount_rna", sep="/"), fig.width=8, fig.height=6)

  p=FeaturePlot(inobj, "log_nfeature", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_lognfeature_rna", sep="/"), fig.width=8, fig.height=6)

  p=FeaturePlot(inobj, "percent.mt", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_percent_mt", sep="/"), fig.width=8, fig.height=6)

  p=FeaturePlot(inobj, "percent.rp", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_percent_rp", sep="/"), fig.width=8, fig.height=6)

  p=FeaturePlot(inobj, "G2M.Score", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_g2mscore", sep="/"), fig.width=8, fig.height=6)

  p=FeaturePlot(inobj, "S.Score", reduction=reduction)
  save_plot(p, paste(outfolder, "fplot_sscore", sep="/"), fig.width=8, fig.height=6)


  p=VlnPlot(inobj, "nCount_RNA", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_ncount_rna", sep="/"), fig.width=12, fig.height=4)

  p=VlnPlot(inobj, "nFeature_RNA", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_nfeature_rna", sep="/"), fig.width=12, fig.height=4)


  p=VlnPlot(inobj, "log_ncount", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_logncount_rna", sep="/"), fig.width=12, fig.height=4)

  p=VlnPlot(inobj, "log_nfeature", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_lognfeature_rna", sep="/"), fig.width=12, fig.height=4)

  p=VlnPlot(inobj, "percent.mt", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_percent_mt", sep="/"), fig.width=12, fig.height=4)

  p=VlnPlot(inobj, "percent.rp", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_percent_rp", sep="/"), fig.width=12, fig.height=4)

  p=VlnPlot(inobj, "G2M.Score", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_g2mscore", sep="/"), fig.width=12, fig.height=4)

  p=VlnPlot(inobj, "S.Score", group.by="idents")
  save_plot(p, paste(outfolder, "vplot_sscore", sep="/"), fig.width=12, fig.height=4)


  if ((reduction == "umap") && ("wnn_clusters" %in% names(inobj@meta.data)))
  {
    p=VlnPlot(inobj, "log_ncount", group.by="wnn_clusters")
    save_plot(p, paste(outfolder, "vplot_wnn_clusters_logncount_rna", sep="/"), fig.width=12, fig.height=4)

    p=VlnPlot(inobj, "log_nfeature", group.by="wnn_clusters")
    save_plot(p, paste(outfolder, "vplot_wnn_clusters_lognfeature_rna", sep="/"), fig.width=12, fig.height=4)

    p=VlnPlot(inobj, "percent.mt", group.by="wnn_clusters")
    save_plot(p, paste(outfolder, "vplot_wnn_clusters_percent_mt", sep="/"), fig.width=12, fig.height=4)
  }


  return(inobj)
}





splitObjListByLibrary = function(objlist)
{
    finalList = list()
    for (objname in names(objlist))
    { 

        xobj = objlist[[objname]]
        xobj$orig_project = objname
        finalList[[objname]] = xobj

    }

    return(finalList)
}



read_gmt <- function(filepath){
  reactome_raw <- readLines(filepath)
  Reactome <- list()
  print("STEP 1 : Read gmt file")
  pb <- txtProgressBar(min = 0, max = length(reactome_raw),style = 3,width = 50, char = "=") 
  for (i in (1:length(reactome_raw))) {
    split_string <- strsplit(reactome_raw[i], split = "\t")
    #print(split_string)
    #print(paste0("Identifier: "),as.character(split_string[[1]][1]))
    set_name<- split_string[[1]][1]
    set_id <- split_string[[1]][2]
    gene_list <- c()
    for (j in (3:length(split_string[[1]]))) {
      #print(split_string[[1]][j])
      gene_list <- c(gene_list,split_string[[1]][j])
    }
    Reactome[[set_id]] <- list(set_name,gene_list)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(Reactome)
}

toPathwayExpression = function(obj.in, genesets, assay="RNA", mean_method="rms", filter_ngenes=3)
{
  print("toPathwayExpression")
  scaled_expression <- obj.in@assays[[assay]]@data
  scaled_expression <- as.data.frame(scaled_expression)
  scaled_expression$gene_names <- rownames(scaled_expression)
  scaled_expression <- as.data.table(scaled_expression)
  res = getExpressionsPerGeneSets(genesets, scaled_expression, mean_method, filter_ngenes)
  mat=getGeneSetMatrixFromList(res)
  return(mat)
}

#input: Reactome = list_object obtained from read_gmt
#       scaled_expression = data.table of expressionvalues of seurat obj (rownames as a column named gene_names), f.e. 
#                           scaled_expression <- obj@assays[[assay]]@data
#                           scaled_expression <- as.data.frame(scaled_expression)
#                           scaled_expression$gene_names <- rownames(scaled_expression)
#                           scaled_expression <- as.data.table(scaled_expression)
#       mean_method = can be "harmonic","geometric" or "rms"
#       filter_ngenes = filter out pathways with only n genes from our data
#output:list object consisting of : n-th element -> (pathway name, gene names, data.table of mean expressionvalues)
getExpressionsPerGeneSets <- function(Reactome, scaled_expression, mean_method, filter_ngenes){
  ReactomeExpression <- list()
  counter <- 1
  gen_counter <- 1
  print("STEP 2 : Getting expression & mean expression values for every gene set")
  pb <- txtProgressBar(min = 0, max = length(Reactome),style = 3,width = 50, char = "=")   
  for (i in Reactome) {
    set_name <- i[[1]][1]
    gene_set <- i[[2]]
    expression_set <- scaled_expression[scaled_expression$gene_names %in% gene_set,]
    if(nrow(expression_set) >= filter_ngenes){ #groesser gleich ?
      es_without_genes  <- expression_set[, -"gene_names"]
      mean_expression <- NULL
      if(mean_method == "harmonic"){
        mean_expression <- apply(es_without_genes, 2,harmonic.mean)#, na.rm = TRUE#
      }else if(mean_method == "geometric"){
        mean_expression <- apply(es_without_genes, 2,geometric.mean)#, na.rm = TRUE#
      }else if(mean_method == "rms"){
        mean_expression <- apply(es_without_genes, 2, rms)#, na.rm = TRUE#
      }

      mean_expression <- as.data.table(t(mean_expression))
      ReactomeExpression[[counter]] <- list(set_name,expression_set$gene_names, mean_expression)
      counter <- counter + 1

    }
    setTxtProgressBar(pb, gen_counter)
    gen_counter <- gen_counter + 1
  }
  close(pb)
  return(ReactomeExpression)
}

#input: list object of getExpressionsPerGeneSets
#output: converts list into matrix: pathways & cells
getGeneSetMatrixFromList <- function(ReactomeExpression){
  MeltExp <- list()
  for (i in (1:length(ReactomeExpression))) {
    id <- ReactomeExpression[[i]][[1]]
    mean_exp <- ReactomeExpression[[i]][[3]]
    mean_exp$gene_set <- id
    MeltExp[[i]] <- mean_exp
  }
  MeanExpressions_dt <- rbindlist(MeltExp)
  MeanExpressions_df <- as.data.frame(MeanExpressions_dt)
  rownames(MeanExpressions_df) <- MeanExpressions_df[,"gene_set"]
  MeanExpressions_df[,"gene_set"] = NULL
  return(MeanExpressions_df)
}


#root mean sqared function
rms <- function(x){
  res <- sqrt(mean(x^2))
  return(res)
}



#
##
###
#### Plotting tools
###
##
#


combine_plot_grid = function(...)
{
  inplot = list(...)
  
  dataList = list()
  for (i in 1:length(inplot))
  {
    dataList[[length(dataList)+1]] = inplot[[i]]$data
  }
  
  p=cowplot::plot_grid(...)
  
  p$data = dataList
  
  return(p)
}

combine_plot_grid_list = function(plotlist, ...)
{
 
  dataList = list()
  for (i in 1:length(plotlist))
  {
    dataList[[i]] = plotlist[[i]]$data
  }
  
  p=cowplot::plot_grid(plotlist=plotlist, ...)
  
  p$data = dataList
  
  return(p)
}



save_plot = function(plotobj,outname, fig.width, fig.height, save.data=TRUE)
{
  print(paste(outname, fig.width, fig.height))
  
  fname=paste(outname, "png", sep=".")
  print(paste("Saving to file", fname))
  png(filename=fname, width = fig.width, height = fig.height, units = 'in', res = 300)#width = fig.width*100, height=fig.height*100)
  plot(plotobj)
  dev.off()
  
  fname=paste(outname, "pdf", sep=".")
  print(paste("Saving to file", fname))
  pdf(file=fname, width = fig.width, height=fig.height)
  plot(plotobj)
  dev.off()
  

  fname=paste(outname, "svg", sep=".")
  print(paste("Saving to file", fname))
  svglite::svglite(file = fname, width = fig.width, height = fig.height)
  plot(plotobj)
  dev.off()
  

  if (save.data)
  {
    if (class(plotobj$data) %in% c("list"))
    {
      print("list case")
      for (i in 1:length(plotobj$data))
      {
        fname = paste(outname,i, "data", sep=".")
        print(paste("Saving to file", fname))
        
        if (class(plotobj$data[[i]]) %in% c("list"))
        {
            print("multi list case")
            for (j in 1:length(plotobj$data[[i]]))
            {
                fname = paste(outname,i, j, "data", sep=".")
                print(paste("Saving to file", fname, class(plotobj$data[[i]][[j]])))

                if (class(plotobj$data[[i]][[j]]) %in% c("list", "waiver"))
                {
                  next()
                }
                write.table(plotobj$data[[i]][[j]], fname, row.names = TRUE, sep="\t")    

            }
        } else {
            
            tryCatch(write.table(plotobj$data[[i]], fname, row.names = TRUE, sep="\t"), error = function(e) NULL)
        }
        
      }
    } else {
      
        fname = paste(outname,"data", sep=".")
        print(paste("Saving to file", fname))

        write.table(plotobj$data, paste(outname, "data", sep="."), row.names = TRUE, sep="\t")
    }
  }
  
  return(plotobj)
}







#
##
###
#### Heatmap
###
##
#


makeComplexExprHeatmap = function( obj.in, plot_gois, group.by="idents", use.scaled.data=TRUE, include_all_clusters=FALSE, title=NULL, scale.by="GLOBAL", scale.limits=NULL, log.expression=TRUE)
{

combine_gois=TRUE
stopifnot(is.null(scale.by) || scale.by %in% c("GLOBAL", "ALL"))

if (!is.null(scale.by))
{
  cat("Scaling data", scale.by)
}


if (is.null(scale.by))
{
    print("Fetching average expression")
    avgexp = AverageExpression(obj.in, assays=c("RNA"), group.by=group.by, slot="data")$RNA

} else {

  if (scale.by %in% c("ALL"))
  {
    print("Fetching average expression")
    avgexp = AverageExpression(obj.in, assays=c("RNA"), group.by=group.by, slot="data")$RNA
  } else {
    print("Fetching global scaled average expression")
    avgexp = AverageExpression(obj.in, assays=c("RNA"), group.by=group.by, slot="scale.data")$RNA
  }
}

if (combine_gois)
{
  new_clusters = c()
  new_genes = c()
  for (gname in names(plot_gois))
  {

    clusters = plot_gois[[gname]]$clusters
    genes = plot_gois[[gname]]$genes

    new_clusters = c(new_clusters, clusters)
    new_genes = c(new_genes, genes)
  }

  plot_gois = list("combined_lists"=list(clusters=new_clusters, genes=new_genes))
}

goi2cluster_genes = list()
allClusters = c()
allGenes = c()

for (gname in names(plot_gois))
{

    clusters = plot_gois[[gname]]$clusters
    genes = plot_gois[[gname]]$genes

    print(gname)
    print(clusters)
    print(genes)

    print("Clusters Missing")
    print(setdiff(clusters, colnames(avgexp)))

    print("Genes Missing")
    print(setdiff(genes, rownames(avgexp)))
   
    clusters = as.character(clusters)

    if (include_all_clusters)
    {
      missingClusters = setdiff( colnames(avgexp), clusters )
      print("Adding missing clusters")
      print(missingClusters)
      clusters = c(clusters, missingClusters)
    }

    allClusters = c(allClusters, clusters)
    allGenes = c(allGenes, genes)

    goi2cluster_genes[[gname]] = list("clusters"=clusters, "genes"=genes)
}
allClusters = unique(allClusters)
allGenes = unique(allGenes)


allScaleMatrices = list()
if (!is.null(scale.by) && (scale.by=="ALL"))
{

  if (length(allClusters) == 1)
  {
    print("1 cluster case")
    mat = avgexp[genes, c(allClusters[1], allClusters[1])]
    mat[1:length(genes), 1] = rep(NA, length(genes))

    origColnames = colnames(mat)
    origColnames[1] = "NA"
    colnames(mat) = origColnames

  } else {
    mat = avgexp[allGenes,allClusters]
  }

  if (log.expression)
  {
    mat = log1p(mat)
  }
  

  matsd = sd(as.vector(mat))
  matmean = mean(mat)
  mat = (mat-matmean) / matsd

  for (gname in names(goi2cluster_genes))
  {

      clusters = goi2cluster_genes[[gname]]$clusters
      genes = goi2cluster_genes[[gname]]$genes

      allScaleMatrices[[gname]] = mat[genes, clusters]

  }
  
}



  clusters = goi2cluster_genes[["combined_lists"]]$clusters
  genes = goi2cluster_genes[["combined_lists"]]$genes

  if (is.null(scale.by) || scale.by=="GLOBAL")
  {
    if (length(clusters) == 1)
    {
        print("1 cluster case")
        mat = avgexp[genes, c(clusters[1], clusters[1])]
        mat[1:length(genes), 1] = rep(NA, length(genes))
        print(mat)

        origColnames = colnames(mat)
        origColnames[1] = "NA"
        colnames(mat) = origColnames

    } else {
        mat = avgexp[genes,clusters]
    }

    if (is.null(scale.by))
    {
      if (log.expression)
      {
        mat = log1p(mat)
      }
    }
    
  
  } else {
    # ALL case!

    mat = allScaleMatrices[[gname]]
  }



  valueTitle = "Average Expression"

  if (!is.null(scale.by))
  {
      valueTitle = paste("Average Scaled Expression", " (",scale.by, ")", sep="")
  }

  if (is.null(title))
  {
    title = paste("Heatmap", gname)
  }

  if (!is.null(scale.limits))
  {

    mat[mat > scale.limits[2]] = scale.limits[2]
    mat[mat < scale.limits[1]] = scale.limits[1]

  }

  p = Heatmap(mat, name=valueTitle, rect_gp = gpar(col = "white", lwd = 2), column_title=title, cluster_rows = FALSE, row_order = rownames(mat), column_order = colnames(mat))

  return(p)

}

library(circlize)


makeComplexExprHeatmapSplit = function( obj.in, plot_gois, split.by="condition", group.by="idents", use.scaled.data=TRUE, include_all_clusters=FALSE, title=NULL, scale.by="GLOBAL", scale.limits=c(-2, 0, 2), log.expression=TRUE)
{

  combine_gois=TRUE
  stopifnot(is.null(scale.by) || scale.by %in% c("GLOBAL", "GROUP", "ALL"))

  if (!is.null(scale.by))
  {
    cat("Scaling data", scale.by)
  }

  if (combine_gois)
  {

    orig_gois = plot_gois

    new_clusters = c()
    new_genes = c()
    for (gname in names(plot_gois))
    {

      clusters = plot_gois[[gname]]$clusters
      genes = plot_gois[[gname]]$genes

      new_clusters = c(new_clusters, clusters)
      new_genes = c(new_genes, genes)
    }

    plot_gois = list()
    plot_gois[["combined_lists"]] = list(clusters=unique(new_clusters), genes=new_genes)

  }

  clusters = plot_gois[["combined_lists"]]$clusters
  genes = plot_gois[["combined_lists"]]$genes

  print(clusters)
  print(genes)

  valueTitle = "Average Expression"

  if (!is.null(scale.by))
  {
    valueTitle = paste("Average Scaled\nExpression (", scale.by, ")", sep="")
  }






  processedMats = list()
  splitByValues = unique(obj.in@meta.data[[split.by]])
  plot = NULL
  for (splitName in splitByValues)
  {
    cells.sel = cellIDForClusters(obj.in, split.by, c(splitName))
    #
    ## Fetching AVG Expression per subset
    #
    if (is.null(scale.by) || scale.by %in% c("ALL", "GROUP"))
    {
      print("Fetching average expression")
      avgexp = AverageExpression(subset(obj.in, cells=cells.sel), assays=c("RNA"), group.by=group.by, slot="data")$RNA

    } else {

      print("Fetching global scaled average expression")
      avgexp = AverageExpression(subset(obj.in, cells=cells.sel), assays=c("RNA"), group.by=group.by, slot="scale.data")$RNA
    }  

    #
    ## Checking all clusters and genes there!
    #
    print("Clusters Missing")
    print(setdiff(clusters, colnames(avgexp)))

    print("Genes Missing")
    print(setdiff(genes, rownames(avgexp)))

    p.clusters = as.character(clusters)

    if (include_all_clusters)
    {
      missingClusters = setdiff( colnames(avgexp), p.clusters )
      print("Adding missing clusters")
      print(missingClusters)

      if (length(missingClusters) > 0)
      {
        p.clusters = c(p.clusters, missingClusters)
      }
    }

    #
    ## Subsetting mat as required
    #
    if (length(p.clusters) == 1)
    {
        print("1 cluster case")
        mat = avgexp[genes, c(p.clusters[1], p.clusters[1])]
        mat[1:length(genes), 1] = rep(NA, length(genes))
        print(mat)

        origColnames = colnames(mat)
        origColnames[1] = "NA"
        colnames(mat) = origColnames

    } else {

      print("Removing clusters because they're not represented in subset")
      print(setdiff(p.clusters, colnames(avgexp)))

      p.clusters = intersect(p.clusters, colnames(avgexp))
      mat = avgexp[genes,p.clusters]
    }

    if (is.null(scale.by) || scale.by%in%c("GROUP", "ALL"))
    {
      if (log.expression)
      {
        mat = log1p(mat)
      }
    }

    if (!is.null(scale.by) && scale.by == "GROUP")
    {
      matsd = sd(as.vector(mat))
      matmean = mean(mat)
      mat = (mat-matmean) / matsd
    }

    processedMats[[splitName]] = mat
  }

  if (!is.null(scale.by) && scale.by=="ALL")
  {


    valueVec = c()
    for (splitName in names(processedMats))
    {
      valueVec = c(valueVec, as.vector(processedMats[[splitName]]))
    }

    allsd = sd(valueVec)
    allmean = mean(valueVec)

    for (splitName in names(processedMats))
    {
      processedMats[[splitName]] = (processedMats[[splitName]]-allmean)/allsd

      print(processedMats[[splitName]])
    }

  }
  
  
  for (splitName in names(processedMats))
  {

    showlegend = splitName == splitByValues[1]

    if (is.null(title))
    {
      plottitle = paste("Heatmap", gname, splitName)
    } else {
      plottitle = paste(title, " (", splitName, ")", sep="")
    }

    mat = processedMats[[splitName]]

    if (!is.null(scale.limits))
    {

      mat[mat > scale.limits[3]] = scale.limits[3]
      mat[mat < scale.limits[1]] = scale.limits[1]

    }
  
    col_fun = colorRamp2(scale.limits, c("blue", "white", "red"))
    print(col_fun)

    p = Heatmap(mat, col=col_fun,show_heatmap_legend=showlegend, name=valueTitle, rect_gp = gpar(col = "white", lwd = 2), column_title=plottitle, cluster_rows = FALSE, row_order = rownames(mat), column_order = colnames(mat))
    
    if (showlegend)
    {
        plot = p
    } else {
        plot = plot + p
    }

  }


  return(plot)


}


#
##
###
#### Various tools to be sorted
###
##
#



write_cell_barcodes = function(scobj, outfolder, valid_orig_idents=NULL)
{
    dir.create(outfolder)

    orig_idents = unique(scobj$orig.ident)

    if (!is.null(valid_orig_idents))
    {
      orig_idents = intersect(valid_orig_idents, orig_idents)
      print("Performing writeout on orig_idents")
      print(orig_idents)
    }

    for (oident in orig_idents)
    {
      print(oident)
      scobjs = subset(scobj, orig.ident == oident)

      allCellNames = colnames(scobjs)

      cellBarcodes = unlist(lapply(strsplit(unlist(lapply(strsplit(colnames(scobjs), "_(?!.*_)", perl = TRUE), function(x){x[2]})), "-"), function(x){x[1]}))

      outfilename = paste(outfolder, paste(oident, "barcodes.tsv", sep="_"), sep="/")
      write.table(cellBarcodes, outfilename, sep="\t", row.names=FALSE, col.names=FALSE, quote=F)
    }

}


toCellPrefix = function(x) {
    datasetprefix = x
    datasetprefix = str_replace_all(datasetprefix, "[./]", "")
    datasetprefix = str_split(datasetprefix, "_")[[1]][1]
    
    return(paste(datasetprefix, sep=""))
}

makesum = function(a, suffix)
{
  out = {}
  out["sum"] = sum(a)

  return(out)
}


cellIDForClusters = function(obj.in, targetVar, clusters)
{

  targetVarDF = as.data.frame(obj.in[[targetVar]])
  #print(paste("orig:", length(rownames(targetVarDF))))
  cellNames = rownames(targetVarDF)[targetVarDF[[targetVar]] %in% clusters]

  #print(length(cellNames))
  return(cellNames)

}


HighlightedDimPlot = function(obj.in, highlightedClusters, group.by="idents", reduction=NULL)
{

p = DimPlot(obj.in, group.by=group.by, reduction=reduction)
g <- ggplot_build(p)

colorGroupDF = unique(g$data[[1]][,c("colour", "group")])
colorGroupDF$group = as.numeric(colorGroupDF$group)-1
print(colorGroupDF)

colorGroupDF[!(colorGroupDF$group %in% highlightedClusters),]$colour = "#808080"

colorList = as.list(colorGroupDF$colour)
names(colorList) = colorGroupDF$group

print(colorList)

p = DimPlot(obj.in, group.by=group.by, reduction=reduction, cols=colorList)
return(p)

}





#
##
### GENE CONVERSION
##
#


convertHumanGeneList <- function(x){
require("biomaRt")
httr::set_config(httr::config(ssl_verifypeer = FALSE))

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}
convertMouseGeneList <- function(x){
require("biomaRt")
httr::set_config(httr::config(ssl_verifypeer = FALSE))

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
mousex <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(mousex))
return(mousex)
}


#
##
### CELL COUNTING
##
#

getCellCountDF = function(scobj, prefix="", group_by="orig.ident", select_by="idents", split_by=NULL, relative=F, outname=NULL,show.percent=F)
{
  
  allClusters = as.character(sort(unique(scobj[[select_by]][,])))
  
  allSplits = NULL
  if (!is.null(split_by))
  {
    allSplits = as.character(sort(unique(scobj[[split_by]][,])))
  }
  
  cellCounts = list()
  
  print("All Clusters")
  print(allClusters)
  print("All Splits")
  print(allSplits)


  for (clusterID in allClusters)
  {
    print(clusterID)
    
    if (!is.null(split_by))
    {
      clusterList = list()
      clusterList[["cluster"]] = clusterID
      cs = scobj[,scobj[[select_by]] == clusterID]
      clusterList[["all"]] = nrow(cs[[group_by]])
      for (split_cat in allSplits)
      {
          print(paste(clusterID, split_cat))
        
          cs = scobj[,scobj[[select_by]] == clusterID & scobj[[split_by]] == split_cat]
          allElems = table(cs[[group_by]])
          cs = NULL

        
          for (grp in names(allElems))
          {
            if (!is.null(prefix))
            {
              ngrp = paste(prefix, paste(split_cat, grp, sep="."), sep=".")
  
            } else {
              ngrp = paste(split_cat, grp, sep=".")
            }
            clusterList[[ngrp]] = allElems[[grp]]
          }
      }
      
    } else {
      
        cs = scobj[,scobj[[select_by]] == clusterID]
      
        allElems = table(cs[[group_by]])
        relevantColumns = as.character(unique(scobj[[group_by]][[group_by]]))

        #print(allElems)
        #print(relevantColumns)

        allElems = allElems[intersect(relevantColumns, names(allElems))]

        clusterList = list()
        clusterList[["cluster"]] = clusterID
        clusterList[["all"]] = nrow(cs[[group_by]])
        cs = NULL
        
        for (grp in names(allElems))
        {
          if (!is.null(prefix))
          {
            ngrp = paste(prefix, grp, sep=".")
          } else {
            ngrp = grp
          }
          clusterList[[ngrp]] = allElems[[grp]]
        }
    }
    
    
    
    cellCounts[[clusterID]] = clusterList
    
  }
  df_bysamplerep = cellCounts %>% map(as.data.frame) %>% bind_rows()
  df_bysamplerep[is.na(df_bysamplerep)] <- 0
 
  rownames(df_bysamplerep) = df_bysamplerep$cluster
  df_bysamplerep$cluster = NULL
  
  
  if (relative)
  {
    df_bysamplerep = sweep(df_bysamplerep,2,colSums(df_bysamplerep),"/")
    
    if (show.percent)
    {
      df_bysamplerep = df_bysamplerep*100;
    }
  }
  
  totals=t(colSums(df_bysamplerep))
  totals.df = data.frame(totals)
  rownames(totals.df) = "Total"
  df_bysamplerep=rbind(df_bysamplerep, totals.df)
  
  df_bysamplerep = cbind("cluster"=rownames(df_bysamplerep), df_bysamplerep)
  rownames(df_bysamplerep) = NULL
  
  if (!is.null(outname))
  {
    write.table(df_bysamplerep, file=outname, row.names = F,  quote=FALSE, sep='\t')
    write_xlsx( df_bysamplerep, path = paste(outname, ".xlsx", sep="") )

  }
  
  return(df_bysamplerep)
  #  
}



makeBarCountPlot = function( seuratObj, outpath, select_by, group_by, size.text = 10, repel.direction="x", max.overlaps=10, fig.height=16)
{

outcheckpath = dirname(paste(outpath, ".tsv", sep=""))
print(outcheckpath)
print(dir.exists(outcheckpath))
if (!dir.exists(outcheckpath))
{
  dir.create(outcheckpath, recursive = TRUE)
}

countByManualCellname = getCellCountDF(seuratObj, prefix="", select_by = select_by, group_by=group_by, relative=TRUE, show.percent=T, outname=paste(outpath, "tsv", sep="."))

cbc_cname = melt(countByManualCellname)
cbc_cname$variable_factors = factor(cbc_cname$variable)#, levels=c("all", ".FIRE_NI_CTRL",".FIRE_NI_KO"))
cbc_cname = cbc_cname[!cbc_cname$cluster=="Total",]

p=DimPlot(seuratObj, group.by = select_by)
g <- ggplot_build(p)
#unique(g$data[[1]]["colour"])
#cellname2color = data.frame(group=g$data[[1]]$group, color=g$data[[1]]$colour)# , cellnames=seuratObj[[group_by]][g$data[[1]]$group])

cbc_cname$cluster = factor(cbc_cname$cluster, levels=levels(p$data[[select_by]]))
cbc_cname$total = 100
cbc_cname$label = paste0(cbc_cname$cluster, " (", round(cbc_cname$value, 2), "% )")

numColumns = length(unique(cbc_cname$variable))

print(cbc_cname)

plot_distribution_celltype = cbc_cname[cbc_cname$variable != "all",] %>%
ggplot(aes(variable_factors, value, fill=cluster)) + 
geom_col(width=0.75)+#scale_fill_manual()+
geom_text_repel(aes(label = paste(round(value, digits = 1))), max.overlaps=max.overlaps, size=size.text, position = position_stack(vjust = .5), direction =repel.direction, segment.size  = 1, segment.color = "black", min.segment.length=0, arrow = arrow(length = unit(0.015, "cm"), angle = 0, type = "closed", ends = "first"), force=2)+
  ylab("Percentage")+xlab("Condition")+
  theme(
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5, size=18),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())


use_width = numColumns * 3

save_plot(plot_distribution_celltype, outpath, fig.height=fig.height, fig.width=use_width)


}

#
##
### DIFFERENTIAL EXPRESSION
##
#


makesummary_getExtExpr = function(a, suffix)
{
out = {}
out["count_expr"] = length(a)

if (length(a) == 0)
{
    f = c(0,0,0,0,0)
    meanA = 0
    stdA = 0
} else {
    f = fivenum(a)
    meanA = mean(a)
    stdA = sd(a)
}

out["min"] = f[1]
out["lower_hinge"] = f[2]
out["median"] = f[3]
out["upper_hinge"] = f[4]
out["max"] = f[5]
out["mean"] = meanA
out["sd"] = stdA

names(out) = paste(names(out), suffix, sep=".")

return(out)
}

getExprData_getExtExpr = function(markerObj, markerCells, sampleSuffix, slot="data", assay="RNA")
{
expTable = GetAssayData(object = subset(x=markerObj, cells=markerCells), slot = slot, assay=assay)
allgenes = rownames(expTable)
cellnames = colnames(expTable)

expt.r = as(expTable, "TsparseMatrix") #TsparseMatrix
expt.df = data.frame(r = expt.r@i + 1, c = expt.r@j + 1, x = expt.r@x)

DT <- data.table(expt.df)
res = DT[, as.list(makesummary_getExtExpr(x, sampleSuffix)), by = r]
anumCol = paste("count_all", sampleSuffix, sep=".")
res[[anumCol]] = length(cellnames)
res$gene = allgenes[res$r]

res = res[,r:=NULL]

return(res)
}

getExtendedExpressionData = function ( scdata, assay="RNA", group.by=NULL)
{
outDF = NULL
DefaultAssay(object=scdata) = assay
print(group.by)
if (is.null(group.by))
{
  clusterIDs = as.character(sort(unique(Idents(scdata))))
} else {
  clusterIDs = as.character(sort(unique(scdata[[group.by]][,])))
}
scCells = Idents(scdata)
scCells = names(scCells)
scCells = unlist(as.character(scCells))

for (clusterID in clusterIDs){
    
    print(clusterID)
    
    cellIdents = Idents(scdata)
    
    if (is.null(group.by))
    {
      cellIdents.c = names(cellIdents[cellIdents == clusterID])
    } else {
      cellIdents.c = colnames(scdata[,scdata[[group.by]] == clusterID])
    }
    
    cellIdents.c = unlist(lapply(cellIdents.c, as.character))  
    
    if (length(cellIdents.c) < 3)
    {
      print(paste("Skipping cluster", clusterID, "due to < 3 cells!"))  
      next
    }
    
    expvals = getExprData_getExtExpr(scdata, cellIdents.c, "cluster", assay=assay)
    expvals$not_expr.cluster = 1- (expvals$count_expr.cluster/expvals$count_all.cluster)
    expvals$expr.cluster = (expvals$count_expr.cluster/expvals$count_all.cluster)

    origColNames = colnames(expvals)
    expvals = as.data.frame(cbind(paste("cluster", clusterID, sep="."), expvals))
    colnames(expvals) = c(c("cluster"), origColNames)
    
    if (!is.data.frame(outDF) || nrow(outDF)==0)
    {
    outDF = expvals
    } else {
    outDF = as.data.frame(rbind(outDF, expvals))
    }
    
}
return(outDF)
}






makesummary_getPop = function(a, suffix)
{
out = {}
out["num"] = length(a)

if (length(a) == 0)
{
    f = c(0,0,0,0,0)
    meanA = 0
} else {
    f = fivenum(a)
    meanA = mean(a)
}

out["min"] = f[1]
out["lower_hinge"] = f[2]
out["median"] = f[3]
out["upper_hinge"] = f[4]
out["max"] = f[5]
out["mean"] = meanA

names(out) = paste(names(out), suffix, sep=".")

return(out)
}

getExprData_getPop = function(markerObj, markerCells, sampleSuffix, slot="data", assay="RNA")
{
expTable = GetAssayData(object = subset(x=markerObj, cells=markerCells), slot = slot, assay=assay)
allgenes = rownames(expTable)
cellnames = colnames(expTable)

expt.r = as(expTable, "TsparseMatrix") # was dgTMatrix
expt.df = data.frame(r = expt.r@i + 1, c = expt.r@j + 1, x = expt.r@x)

DT <- data.table(expt.df)
res = DT[, as.list(makesummary_getPop(x, sampleSuffix)), by = r]
anumCol = paste("anum", sampleSuffix, sep=".")
res[[anumCol]] = length(cellnames)
res$gene = allgenes[res$r]

res = res[,r:=NULL]

return(res)
}

getDEXpressionDF = function ( scdata, markers, assay="SCT", group.by=NULL)
{
outDF = NULL
DefaultAssay(object=scdata) = assay
print(group.by)
if (is.null(group.by))
{
  clusterIDs = as.character(sort(unique(Idents(scdata))))
} else {
  clusterIDs = as.character(sort(unique(scdata[[group.by]][,])))
}
scCells = Idents(scdata)
scCells = names(scCells)
scCells = unlist(as.character(scCells))
for (clusterID in clusterIDs){
    
    print(clusterID)
    
    cellIdents = Idents(scdata)
    
    if (is.null(group.by))
    {
      cellIdents.c = names(cellIdents[cellIdents == clusterID])
    } else {
      cellIdents.c = colnames(scdata[,scdata[[group.by]] == clusterID])
    }
    
    cellIdents.c = unlist(lapply(cellIdents.c, as.character))  
    
    cellIdents.bg = setdiff(unlist(lapply(names(cellIdents), as.character)), cellIdents.c)

    if (length(cellIdents.c) < 3)
    {
      print(paste("Skipping cluster", clusterID, "due to < 3 cells!"))  
      next
    }
    
    expvals = getExprData_getPop(scdata, cellIdents.c, "cluster", assay=assay)
    expvals.bg = getExprData_getPop(scdata, cellIdents.bg, "bg", assay=assay)
    modmarkers = markers[[clusterID]]
    modmarkers$gene = rownames(modmarkers)
    
    markerdf = as.data.frame(modmarkers)
    
    if ((nrow(markerdf) > 0) && (nrow(expvals) > 0))
    {
      expvals = merge(markerdf, expvals, all.x=T, by.x="gene", by.y = "gene")  
    }
    
    if ((nrow(expvals) > 0) && (nrow(expvals.bg) > 0))
    {
      expvals = merge(expvals, expvals.bg, all.x=T, by.x="gene", by.y = "gene")  
    }
    
    expvals = as.data.frame(cbind(clusterID, expvals))
    
    if (!is.data.frame(outDF) || nrow(outDF)==0)
    {
    outDF = expvals
    } else {
    outDF = as.data.frame(rbind(outDF, expvals))
    }
    
}
return(outDF)
}



makeDEResults = function(inobj, group.by=NULL, assay="SCT", test="wilcox")
{
  if (is.null(group.by))
  {
    clusterIDs = as.character(sort(unique(Idents(inobj))))
  } else {
    clusterIDs = as.character(sort(unique(inobj[[group.by]][,])))
  }
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
  
      cellIdents = Idents(inobj)
      
      if (is.null(group.by))
      {
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
  
      } else {
        cellIdents.c = colnames(inobj[,inobj[[group.by]] == clusterID])
      }
      
      cellIdents.c = unlist(lapply(cellIdents.c, as.character))
  
      print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
      print(length(cellIdents.c) < 3)
  
      if (length(cellIdents.c) < 3)
      {
        print(paste("Skipping cluster", clusterID, "due to < 3 cells!"))  
        next
      }
      deMarkers = FindMarkers(inobj, assay=assay, ident.1 = cellIdents.c, test.use=test) 
      retList[[clusterID]] = deMarkers
  
  }

  retList = getDEXpressionDF(inobj, retList, assay=assay, group.by=group.by)
  
  return(retList)
}



get_population_expression_data = function(scobj, group, outname, assay="RNA", slot="counts", addScores = NULL)
{

if (!is.null(outname))
{
if (!dir.exists(dirname(paste(outname, ".tsv", sep=""))))
{
  dir.create(outname, recursive = TRUE)
}
}

exprData = list()

for (cellPop in unique(as.character(unlist(scobj[[group]]))))
{
    varPop = str_to_lower( str_replace_all(
                            str_replace_all(#
                            str_replace_all( cellPop, "\\(|\\)| |,", "_"),
                            "__", "_"),
                            "_$", "")
                        )
    print(paste(cellPop, varPop))
    
    allPopCells = scobj[[group]]
    allPopCells$cnames = rownames(allPopCells)
    cellPopCells = allPopCells[allPopCells[[group]] == cellPop, ]$cnames
    print(paste("Number of cells: ", length(cellPopCells)))

    exprData[[varPop]] = getExprData_getPop(markerObj=scobj, markerCells=cellPopCells, sampleSuffix=varPop, slot=slot, assay=assay)
}


meanExprData = list()

for (name in names(exprData))
{
    
    exprDf = as.data.frame(exprData[[name]])
    subdf = exprDf[ ,c("gene", paste("mean", name, sep=".")) ]

    meanExprData[[name]] = subdf
}

cellnames_manualExprDF = Reduce(function(x,y) merge(x = x, y = y, by = "gene", all.x=T, all.y=T), meanExprData)
cellnames_manualExprDF[is.na(cellnames_manualExprDF)] = 0

if (!is.null(addScores))
{

  for (moduleName in addScores)
  {
    fModuleName = paste(moduleName, 1, sep="")

    if (moduleName %in% colnames(scobj@meta.data))
    {
      fModuleName = moduleName
      print(paste("Using fModuleName", fModuleName))
    } else {
      print(paste("Using fModuleName", fModuleName))
    }

    print(moduleName)
    print(fModuleName)

    module_values = list()
    module_values[["gene"]] = moduleName
    for (cellPop in unique(as.character(unlist(scobj[[group]]))))
    {
        varPop = str_to_lower( str_replace_all(
                                str_replace_all(#
                                str_replace_all( cellPop, "\\(|\\)| |,", "_"),
                                "__", "_"),
                                "_$", "")
                            )

        allPopCells = scobj[[group]]
        allPopCells$cnames = rownames(allPopCells)
        cellPopCells = allPopCells[allPopCells[[group]] == cellPop, ]$cnames
        print(paste("Number of cells: ", length(cellPopCells)))

        moduleExpression = scobj@meta.data[fModuleName]
        print(head(moduleExpression))
        module_expression = moduleExpression[rownames(moduleExpression) %in% cellPopCells,]
        module_expression = mean(as.numeric(unlist(module_expression)))
        moduleScoreName = paste("mean.", varPop, sep="")
        print(moduleScoreName)
        print(module_expression)

        module_values[[moduleScoreName]] = module_expression
    }

    print(module_values)
    module_values = module_values[ colnames(cellnames_manualExprDF) ]
    print(module_values)
    cellnames_manualExprDF = rbind(cellnames_manualExprDF, module_values)   
  }
  
}

if (!is.null(outname))
{
  write.table(cellnames_manualExprDF, file = paste(outname, ".tsv", sep=""), quote=FALSE, sep = "\t", row.names = F)
  write_xlsx( cellnames_manualExprDF, path = paste(outname, ".xlsx", sep="") )
}

return(cellnames_manualExprDF)


}



compareClusters = function(scdata, cellsID1, cellsID2, suffix1, suffix2, prefix="cluster", test="MAST", assay="RNA", outfolder="./", fcCutoff=0.25, all=FALSE, heatmap.plot=FALSE, heatmap.pcutoff=0.05, heatmap.addgenes=NULL)
{
    logfc.threshold = fcCutoff
    
    if (all==TRUE)
    {
    logfc.threshold = 0.01  
    }

    if (!dir.exists(outfolder)){
        dir.create(outfolder, recursive = TRUE)
        print(outfolder)
        print("Dir created!")

    } else {
        print(outfolder)
        print("Dir already exists!")
    }
    
    markers = FindMarkers(scdata, assay=assay, ident.1 = cellsID1, ident.2 = cellsID2, test.use=test, logfc.threshold=logfc.threshold)
    
    outvalues1 = getExprData_getPop(scdata, cellsID1, suffix1, assay=assay)
    outvalues2 = getExprData_getPop(scdata, cellsID2, suffix2, assay=assay) 
    
    
    markers$gene = rownames(markers)
    joinedData = merge(markers, outvalues1, by="gene", all=T)
    joinedData = merge(joinedData, outvalues2, by="gene", all=T)  
    
    joinedData = joinedData[!is.na(joinedData$p_val),]
    
    suffix1=str_replace_all(str_replace_all(suffix1, "\\/", "_"), " ", "_")
    suffix2=str_replace_all(str_replace_all(suffix2, "\\/", "_"), " ", "_")

    
    outfile = paste(outfolder, "/", prefix, ".", suffix1, "_", suffix2, ".tsv", sep="")
    
    message(outfile)
    write.table(joinedData, file=outfile, row.names = F,  quote=FALSE, sep='\t')
    
    outfile = paste(outfolder, "/", prefix, ".", suffix1, "_", suffix2, ".xlsx", sep="")
    
    message(outfile)
    write_xlsx(joinedData, path=outfile)

    if (heatmap.plot)
    {
      genes.interest = (joinedData[joinedData$p_val_adj < heatmap.pcutoff,] %>% arrange(p_val_adj) %>% head(40) %>% top_n(40, avg_log2FC))$gene

      if (!is.null(heatmap.addgenes))
      {
        genes.interest = unique(c(genes.interest, heatmap.addgenes))
      }

      if (length(genes.interest) > 0)
      {
        scaleColors = c("#8b0000", "grey", "#008b2b")

        obj.rel = subset(scdata, cells=c(cellsID1, cellsID2))
        cellAnnot = colnames(obj.rel)
        names(cellAnnot) = cellAnnot

        cellAnnot[cellsID1] = suffix1
        cellAnnot[cellsID2] = suffix2
        obj.rel$heatmap_annot = cellAnnot

        obj.rel = ScaleData(obj.rel, features=genes.interest)

        p=DoHeatmap(obj.rel, genes.interest, group.by="heatmap_annot")+ scale_fill_gradientn(colors = scaleColors)
        p=p + ggtitle("Heatmap of DE genes; Data scaled by shown genes.")
        save_plot(p, paste(outfolder, paste("hplot_", prefix, ".", suffix1, "_", suffix2, sep=""), sep="/"), fig.width=7, fig.height=0.3*length(genes.interest))
      }
    }
    
      
    return(joinedData)
}

compareCellsByCluster = function(inobj, cellsID1, cellsID2, suffix1, suffix2, group.by=NULL, assay="RNA", test="t", outfolder="./", fcCutoff=0.25, all=FALSE, heatmap.plot=FALSE, heatmap.pcutoff=0.05, heatmap.addgenes=NULL)
{
  if (is.null(group.by))
  {
    clusterIDs = as.character(sort(unique(Idents(inobj))))  
  } else {
    clusterIDs = as.character(sort(unique(inobj[[group.by]][,])))
  }
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
  
      cellIdents = Idents(inobj)
      
      if (is.null(group.by))
      {
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
  
      } else {
        cellIdents.c = colnames(inobj[,inobj[[group.by]] == clusterID])
      }
      
      cellIdents.c = unlist(lapply(cellIdents.c, as.character))
  
      print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
      
      cells1 = intersect(cellsID1, cellIdents.c)
      cells2 = intersect(cellsID2, cellIdents.c)
      
      print(c(length(cells1), length(cells2)))
      
      if (length(cells1) < 3)
      {
        print("Cells1 too few")
        next
      }
      
      if (length(cells2) < 3)
      {
        print("Cells2 too few")
        next
      }

      clusterID_file=str_replace_all(str_replace_all(clusterID, "\\/", "_"), " ", "_")

  
      deMarkers = compareClusters(scdata=inobj,
                                  cellsID1=cells1,
                                  cellsID2=cells2,
                                  prefix= paste("cluster", clusterID_file, sep="_"),
                                  suffix1=suffix1,#paste(suffix1, clusterID, sep="_"),
                                  suffix2=suffix2, #paste(suffix2, clusterID, sep="_"),
                                  test=test, fcCutoff=fcCutoff, assay=assay, outfolder=outfolder,
                                  heatmap.plot=heatmap.plot, heatmap.pcutoff=heatmap.pcutoff, heatmap.addgenes=heatmap.addgenes)
  
  
      retList[[clusterID]] = deMarkers
  
  }
  
  return(retList)
}













makeVolcanos = function(loMG, titlePrefix, outname, restrict_labels=NULL, turnExpression=F, colors = list(neg=list(sig="#448CCA", nosig="#B4D1E9"), pos=list(sig="#F47B78", nosig="#FACAC9")), FCcutoff=0.5, pCutoff = 0.05, highlightGene=NULL)
{

    outfolder = dirname(outname)[1]

    if (!dir.exists(outfolder)){
        dir.create(outfolder, recursive = TRUE)
        print(outfolder)
        print("Dir created!")

    } else {
        print(outfolder)
        print("Dir already exists!")
    }
  
  for (cName in names(loMG))
  {
    print(cName)
    
    title=paste(titlePrefix, "(", cName, ")", sep=" ")
    subtitle=""
    
    
    indf = loMG[[cName]]

    print(dim(indf))

    if (dim(indf)[1] == 0)
    {
      print("Skipping sample for 0 rows:")
      print(cName)
      print(dim(indf))
      next()
    }
    
    cName = str_replace_all(str_replace_all(cName, "\\/", "_"), " ", "_")
    popName = str_to_lower( str_replace_all(str_replace_all(str_replace_all( cName, "\\(|\\)| ", "_"), "__", "_"), "_$", "") )

    plotlabels = NULL
    
    if (!is.null(restrict_labels))
    {
      if (cName %in% names(restrict_labels))
      {
          cInfoElem = restrict_labels[[cName]]
        
          popName = cInfoElem$fname
          plotlabels = cInfoElem$genes
      } else {
        selDF = indf[indf$p_val_adj < pCutoff,] %>% arrange(p_val_adj) %>% head(40)
        plotlabels = selDF$gene
      }
    }    

    if (!is.null(highlightGene))
    {

      possibleGenes = intersect(indf$gene, highlightGene)
      print("possibleGenes")
      print(possibleGenes)

      if (length(possibleGenes) > 0)
      {
        subtitle = "including highlight gene"
        
        plotlabels = c(plotlabels, possibleGenes)
      }      
    }


    if (turnExpression)
    {
      indf$avg_log2FC = -indf$avg_log2FC
    }
    

    
    keyvals <- ifelse(indf$avg_log2FC > 0,
                      
                      ifelse(indf$p_val_adj < pCutoff & indf$avg_log2FC > FCcutoff, colors$pos$sig, colors$pos$nosig)
                      ,
                      ifelse(indf$avg_log2FC <= 0,
                             ifelse(indf$p_val_adj < pCutoff & indf$avg_log2FC < -FCcutoff, colors$neg$sig, colors$neg$nosig)
                             ,
                             'black'))
    


    keyvals[is.na(keyvals)] <- 'black'
    names(keyvals)[keyvals == '#F47B78'] <- 'UpReg sig'
    names(keyvals)[keyvals == '#FACAC9'] <- 'UpReg non-sig'
    names(keyvals)[keyvals == '#448CCA'] <- 'DownReg sig'
    names(keyvals)[keyvals == '#B4D1E9'] <- 'DownReg non-sig'
    
    txtLabelSize = 5
    axisLabelSize=12
    legendLabSize=12
    legendIconSize=5.0
  
    filename=paste(outname, popName, "png", sep=".")
    print(filename)
    png(filename=filename,width = 1200, height = 700)
    p=EnhancedVolcano(indf[,c('avg_log2FC', 'p_val_adj')],
      lab = indf$gene,
      selectLab=plotlabels,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      colCustom = keyvals,
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2.0,
      legendPosition = 'right',
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      title = title,
      subtitle = subtitle,
      axisLabSize=axisLabelSize,
      labSize = txtLabelSize,
      legendLabSize = legendLabSize,
      legendIconSize = legendIconSize,
      ylab = bquote(~-Log[10]~italic(adj~P-Value)),
      gridlines.major = FALSE,
      gridlines.minor = FALSE,
      max.overlaps=100
     )
    plot(p)
    dev.off()
    
    
        filename=paste(outname, popName, "pdf", sep=".")
    print(filename)
        pdf(filename,width = 12, height = 7)
    p=EnhancedVolcano(indf,
      lab = indf$gene,
      selectLab=plotlabels,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      colCustom = keyvals,
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2.0,
      legendPosition = 'right',
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      title = title,
      subtitle = subtitle,
      axisLabSize=axisLabelSize,
      labSize = txtLabelSize,
      legendLabSize = legendLabSize,
      legendIconSize = legendIconSize,
      ylab = bquote(~-Log[10]~italic(adj~P-Value)),
      gridlines.major = FALSE,
      gridlines.minor = FALSE,
      max.overlaps=100
     )
    plot(p)
    dev.off()
    
    filename=paste(outname, popName, "svg", sep=".")
    print(filename)
    svglite(file = filename, width = 12, height = 7)
    p=EnhancedVolcano(indf,
      lab = indf$gene,
      selectLab=plotlabels,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      colCustom = keyvals,
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2.0,
      legendPosition = 'right',
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      title = title,
      subtitle = subtitle,
      axisLabSize=axisLabelSize,
      labSize = txtLabelSize,
      legendLabSize = legendLabSize,
      legendIconSize = legendIconSize,
      ylab = bquote(~-Log[10]~italic(adj~P-Value)),
      gridlines.major = FALSE,
      gridlines.minor = FALSE,
      max.overlaps=100
     )
    plot(p)
    dev.off()
      
  }
  
}

log_both <- function(x){ifelse(x == 0, 0, log(abs(x)) * sign(x))}

log_both_trans <- 
  function(){
    trans_new(name = 'log_both', 
              transform = log_both,
              inverse = log_both) #not clear what `inverse` does
}


cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

make_descr_label = function(plot, descr)
{
  descrLabel <- ggdraw() + draw_label(descr, fontface='bold', angle = 0)
  
  pe = combine_plot_grid_list(plotlist=list("a"=descrLabel, "b"=plot), ncol=1, nrow=2, labels=NULL,rel_heights = c(0.1, 1), align = "h", axis = "l")
  
  return(pe)
}










enhancedDotPlot = function(scobj, plotElems, featureGenes = c(""), group.by="cellnames_manual", col.min = -3, col.max = 3, cols = c("grey", "blue"), title="", scaled=T, scale.by="GROUP", rotate.x=F, abundance.perelem=FALSE, assay="RNA")
{

  stopifnot(is.null(scale.by) ||scale.by %in% c("GROUP", "FEATURE", "ALL", "GLOBAL"))

  use.slot="data"

  if (!is.null(scale.by) && scale.by == "GLOBAL")
  {
    use.slot="scale.data"
  }

  plotData = list()
  allFoundFeatures = c() 
  ctFractions = list()
  
  scFullTable = table(scobj[[group.by]])
  scFullDf = as.data.frame(scFullTable)

  fillLimitsMax = 0
  elemCount = 0
  for (plotName in names(plotElems))
  {
    plotDef = plotElems[[plotName]]
    plotCells = plotDef$cells
    plotDescr = plotDef$label
    
    scobj_subset = subset(scobj, cells=plotCells)

    avgexpMat = AverageExpression(scobj_subset, features=featureGenes, group.by=group.by, assay=assay, slot=use.slot)$RNA
    
    if (use.slot=="data")
    {
      avgexpMat = log1p(avgexpMat)
    }
    

    avgExpr = as.data.frame(data.table(features.plot = rownames(avgexpMat), id = colnames(avgexpMat), avg.exp = c(as.matrix(avgexpMat))))

    p=DotPlot(scobj_subset, features=featureGenes, group.by=group.by, assay=assay)

    avgExpr = merge(avgExpr, p$data[, c("id", "features.plot", "pct.exp")], by=c("id", "features.plot"))
    
    scTable = table(scobj_subset[[group.by]])
    scDf = as.data.frame(scTable)
    scobj_subset = NULL
    
    if (abundance.perelem)
    {
      scDf$perc = scDf$Freq / sum(scDf$Freq) # per elem abundance!
    } else {
      scDf$perc = scDf$Freq / sum(scFullDf$Freq) # total abundance!
    }

    fillLimitsMax = max(c(fillLimitsMax,max(scDf$perc)))

    ctFractions[[plotName]] = scDf   
    plotData[[plotName]] = avgExpr
  }

  idLevels = levels(scobj@meta.data[, c(group.by)])
  featureLevels = featureGenes

  
  # initialize Values
  for (plotName in names(plotData))
  {
    #plotData[[plotName]]$avg.exp.scaled2 = 0
    plotData[[plotName]]$id = as.character(plotData[[plotName]]$id)
    plotData[[plotName]]$features.plot = as.character(plotData[[plotName]]$features.plot)
  }

  #print(plotData[["Young"]])
  
  allIDs = c()
  allFeatures = c()
  for (plotName in names(plotData))
  {
      allIDs = c(allIDs, plotData[[plotName]]$id)
      allFeatures = c(allFeatures, plotData[[plotName]]$features.plot)
  }
  allIDs = unique(allIDs)
  allFeatures = unique(allFeatures)

  if ((is.null(scale.by) || scale.by == "GLOBAL"))
  {

    for (plotName in names(plotData))
    {

      plotDF = plotData[[plotName]]

      #plotDF$id = factor(plotDF$id, levels = mixedsort(as.character(unique(plotDF$id))))
      #plotDF$idn= as.numeric(plotDF$id)
      #plotDF$features.plot = factor(plotDF$features.plot, levels=unique(plotDF$features.plot))
      
      #just copy old values
      plotDF[, "avg.exp.scaled2"] = plotDF$avg.exp

      plotData[[plotName]] = plotDF
    }

  } else if (scale.by == "GROUP")
  {

    for (plotName in names(plotData))
    {

      plotDF = plotData[[plotName]]

      #just copy old values
      plotDF[, "avg.exp.scaled2"] = scale(plotDF$avg.exp)

      plotData[[plotName]] = plotDF
    }


  } else if (scale.by == "FEATURE")
  {
    print("SCALING BY FEATURE")
    # calculate avg.exp.scaled2 for each feature
    for (featureName in allFoundFeatures)
    {
      print(featureName)
      allUnscaledValues = NULL
      for (plotName in names(plotData))
      {
        pData = plotData[[plotName]]
        missingIDs = setdiff(allIDs, unique(pData$id))
                
        for (mid in missingIDs)
        {
          for (feature in allFeatures)
          {
            dfrow = data.frame(avg.exp=0.0,pct.exp=0.0,features.plot=feature, id=mid, avg.exp.scaled=0.0, avg.exp.scaled2=0.0)
            pData = rbind.data.frame(pData, dfrow, stringsAsFactors = F)
          }
        }
        
        #pData$id = factor(pData$id, levels = allIDs)
        #pData$features.plot = factor(pData$features.plot, levels=allFeatures)
        
        plotData[[plotName]] = pData
        
        
        if (is.null(allUnscaledValues))
        {
          allUnscaledValues = data.frame(plotName=pData[ pData$features.plot==featureName, ]$avg.exp)
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
          allUnscaledValues[["plotName"]] = NULL
          
        } else {
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
        }
        
      }

      allUnscaledValues$rnames = as.numeric(rownames(allUnscaledValues))
      allUnscaledLong = allUnscaledValues %>% gather(Type, Value, names(plotData))
      allUnscaledLong$Value = scale(allUnscaledLong$Value)
      allScaledValues = allUnscaledLong %>% spread(Type, Value) %>% arrange( order(rnames))
      
      for (plotName in names(plotData))
      {

        pData = plotData[[plotName]]
               
        origData = pData[pData$features.plot==featureName, ]
        pData[pData$features.plot==featureName, "avg.exp.scaled2"] = pData[pData$features.plot==featureName, "avg.exp"]
        #pData$idn=as.numeric(pData$id)

        plotData[[plotName]] = pData
      }


      
    }

  } else if (scale.by == "ALL") {

    combinedDataDF = data.frame()

    for (plotName in names(plotData))
    {
      pData = plotData[[plotName]]
      pData["plotpart"] = plotName

      combinedDataDF = rbind(combinedDataDF, pData)
    }

    combinedDataDF[, c("avg.exp.scaled2")] = scale(combinedDataDF$avg.exp)


    # reorder
    originalGroups = scobj@meta.data[,group.by]
    if (is.factor(originalGroups))
    {
      originalSorting = levels(originalGroups)
      combinedDataDF$id = factor(combinedDataDF$id, levels = originalSorting)
    } else {
      combinedDataDF$id = factor(combinedDataDF$id, levels = mixedsort(as.character(unique(combinedDataDF$id))))
    }


    for (plotName in names(plotData))
    {
      subDF = combinedDataDF[combinedDataDF$plotpart == plotName, c("features.plot", "id", "avg.exp.scaled2", "plotpart")]

      pData = plotData[[plotName]]

      pData = merge(pData, subDF, by.x=c("features.plot", "id"), by.y=c("features.plot", "id"))

      #pData$avg.exp.scaled2 = subDF[pData$features.plot,]$avg.exp.scaled2
      #pData$id = plotData$id
      #pData$idn= as.numeric(pData$id)
      #pData$features.plot = factor(pData$features.plot, levels=unique(pData$features.plot))
      
      plotData[[plotName]] = pData

    }
      
  } else {
    stopifnot(FALSE)
  }


  for (plotName in names(plotData))
  {
    pData = plotData[[plotName]]

    pData$id = factor(pData$id, levels=idLevels)
    pData$idn = as.numeric(pData$id)
    pData$features.plot = factor(pData$features.plot, levels=featureLevels)

    plotData[[plotName]] = pData
  }


  # prepare ID DF
  idDF = data.frame()
  for (plotName in names(plotData))
  {
    pData = plotData[[plotName]]
    idDF = rbind(idDF, pData[, c("id", "idn")])
  }

  rownames(idDF) = NULL
  idDF = unique(idDF)

  if (!is.null(scale.by))
  {
    title.expression = paste("Avg. Expression\n(scaled by ", scale.by, ")", sep="")
  } else {
    title.expression = "Avg. Expression"
  }

  if (abundance.perelem)
  {
    title.cellabundance = "Cell Abundance\n(per condition)"
  } else {
    title.cellabundance = "Cell Abundance\n(all object cells)"
  }

  plotList = list()
  
  for (plotName in names(plotData))
  {
    pData = plotData[[plotName]]

    print(pData)

    pData2 = merge(x=pData,y=ctFractions[[plotName]],by.x="id", by.y=group.by,all.x=TRUE)
    pData <-pData2 %>% mutate(featuren=as.numeric(features.plot), percn=100*perc)  

    print(pData)
    
    pData[pData$avg.exp.scaled2>col.max, 'avg.exp.scaled2'] = col.max
    pData[pData$avg.exp.scaled2<col.min, 'avg.exp.scaled2'] = col.min

    fillLimits = c(0, ceiling(fillLimitsMax*10)*10)

    print(pData)
    
    plotElem <- ggplot(pData) +
      scale_x_continuous(breaks=pData$featuren, labels=pData$features.plot) +
      scale_y_continuous(breaks=idDF$idn, labels=idDF$id, limits = c(min(idDF$idn)-0.6, max(idDF$idn)+0.6))+
      geom_rect(aes(xmin=featuren-.5, xmax=featuren+.5, ymin = idn-0.5, ymax = idn+0.5, fill=percn), alpha = 0.4, linetype="blank") +
      scale_fill_distiller(palette='Spectral', limits = fillLimits)+
      scale_size_continuous(range = c(0, 10))+
      geom_point(aes(x=featuren, y=idn, colour = avg.exp.scaled2, size = pct.exp)) +
      scale_color_gradient2(limits=c(col.min, col.max), low = cols[1], mid=cols[2], high = cols[3])+
      guides(color=guide_colourbar(title=title.expression, order = 1),
      size=guide_legend(title="Percent Expressing", order = 2),
      fill=guide_colourbar(title=title.cellabundance, order = 3))
    
    plotElem = plotElem + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line())

    if (rotate.x)
    {
      plotElem = plotElem + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }

    plotList[[plotName]] = plotElem
  }
    
  print("Final plot checks")

  plot_list=list()
  for (plotName in names(plotList))
  {
    
    plotElem = plotList[[plotName]]      
    
    print("descr label")
    pe = make_descr_label(plotElem, plotName) #plotElems[[plotName]]$label
  
    plot_list[[plotName]] = pe
 
  }

  print("Preparing Legend")
  
  # extract a legend that is laid out horizontally #
  legend_b <- get_legend(
    plotElem + 
      guides(color = guide_colorbar(title = title.expression, direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      theme(legend.position = "bottom", legend.box = "vertical")
  )
  
  title <- ggdraw() + draw_label(title, fontface='bold')
  print("Preparing plot")

  dataList = list()
  for (i in 1:length(plot_list))
  {
    dataList[[i]] = plot_list[[i]]$data
  }
  ap=combine_plot_grid_list(
    plotlist = plot_list,
    labels = NULL,
    nrow=1,
    align = "v", axis="bt"
  )

  finalPlotList = list(title, ap, legend_b)

  fp = combine_plot_grid_list(plotlist = finalPlotList, ncol = 1, rel_heights = c(.05, 1, 0.2) )
  fp$data = dataList

  return(fp)
}




















#' Creates an enhanced DotPlot
#'
#' @param scobj Seurat object for plotting
#' @param plotElems list of list(cells,label) where each entry defines one condition/split of the plot
#' @param featureGenes genes to plot
#' @param group.by name of the meta.data column used for grouping cells
#' @param scale.by how to show/scale the expression values (GROUP; FEATURE; ALL; GLOBAL)
#' @param col.min lower bound of the scale
#' @param col.max upper bound of the scale
#' @param cols color for the expression values
#' @param title title of the plot
#' @param rotate.x whether to rotate x-axis labels
#' @param abundance.perelem whether the cell group abundance is to be calculated on the global Seurat object, or per group
#' @param assay which assay to use for retrieving the expression values
#'
#' @return ggplot2 object
#'
#'
#' @export
enhancedDotPlot = function(scobj, plotElems, featureGenes = c(""), group.by="cellnames_manual", col.min = -3, col.max = 3, cols = c("grey", "blue"), title="", scale.by="GROUP", rotate.x=F, abundance.perelem=FALSE, assay="RNA")
{

  stopifnot(is.null(scale.by) ||scale.by %in% c("GROUP", "FEATURE", "ALL", "GLOBAL"))

  use.slot="data"

  if (!is.null(scale.by) && scale.by == "GLOBAL")
  {
    use.slot="scale.data"
  }

  plotData = list()
  allFoundFeatures = c() 
  ctFractions = list()
  
  scFullTable = table(scobj[[group.by]])
  scFullDf = as.data.frame(scFullTable)

  fillLimitsMax = 0
  elemCount = 0
  for (plotName in names(plotElems))
  {
    plotDef = plotElems[[plotName]]
    plotCells = plotDef$cells
    plotDescr = plotDef$label
    
    scobj_subset = subset(scobj, cells=plotCells)

    avgexpMat = Seurat::AverageExpression(scobj_subset, features=featureGenes, group.by=group.by, assay=assay, slot=use.slot)$RNA
    
    if (use.slot=="data")
    {
      avgexpMat = log1p(avgexpMat)
    }
    

    avgExpr = as.data.frame(data.table::data.table(features.plot = rownames(avgexpMat), id = colnames(avgexpMat), avg.exp = c(as.matrix(avgexpMat))))

    p=Seurat::DotPlot(scobj_subset, features=featureGenes, group.by=group.by, assay=assay)

    avgExpr = merge(avgExpr, p$data[, c("id", "features.plot", "pct.exp")], by=c("id", "features.plot"))
    
    scTable = table(scobj_subset[[group.by]])
    scDf = as.data.frame(scTable)
    scobj_subset = NULL
    
    if (abundance.perelem)
    {
      scDf$perc = scDf$Freq / sum(scDf$Freq) # per elem abundance!
    } else {
      scDf$perc = scDf$Freq / sum(scFullDf$Freq) # total abundance!
    }

    fillLimitsMax = max(c(fillLimitsMax,max(scDf$perc)))

    ctFractions[[plotName]] = scDf   
    plotData[[plotName]] = avgExpr
  }

  idLevels = levels(scobj@meta.data[, c(group.by)])
  featureLevels = featureGenes

  
  # initialize Values
  for (plotName in names(plotData))
  {
    #plotData[[plotName]]$avg.exp.scaled2 = 0
    plotData[[plotName]]$id = as.character(plotData[[plotName]]$id)
    plotData[[plotName]]$features.plot = as.character(plotData[[plotName]]$features.plot)
  }

  #print(plotData[["Young"]])
  
  allIDs = c()
  allFeatures = c()
  for (plotName in names(plotData))
  {
      allIDs = c(allIDs, plotData[[plotName]]$id)
      allFeatures = c(allFeatures, plotData[[plotName]]$features.plot)
  }
  allIDs = unique(allIDs)
  allFeatures = unique(allFeatures)

  if ((is.null(scale.by) || scale.by == "GLOBAL"))
  {
    print("SCALING BY GLOBAL")

    for (plotName in names(plotData))
    {

      plotDF = plotData[[plotName]]

      #plotDF$id = factor(plotDF$id, levels = mixedsort(as.character(unique(plotDF$id))))
      #plotDF$idn= as.numeric(plotDF$id)
      #plotDF$features.plot = factor(plotDF$features.plot, levels=unique(plotDF$features.plot))
      
      #just copy old values
      plotDF[, "avg.exp.scaled2"] = plotDF$avg.exp

      plotData[[plotName]] = plotDF
    }

  } else if (scale.by == "GROUP")
  {
    print("SCALING BY GROUP")

    for (plotName in names(plotData))
    {

      plotDF = plotData[[plotName]]

      #just copy old values
      plotDF[, "avg.exp.scaled2"] = scale(plotDF$avg.exp)

      plotData[[plotName]] = plotDF
    }


  } else if (scale.by == "FEATURE")
  {
    print("SCALING BY FEATURE")
    # calculate avg.exp.scaled2 for each feature
    for (featureName in allFoundFeatures)
    {
      print(featureName)
      allUnscaledValues = NULL
      for (plotName in names(plotData))
      {
        pData = plotData[[plotName]]
        missingIDs = setdiff(allIDs, unique(pData$id))
                
        for (mid in missingIDs)
        {
          for (feature in allFeatures)
          {
            dfrow = data.frame(avg.exp=0.0,pct.exp=0.0,features.plot=feature, id=mid, avg.exp.scaled=0.0, avg.exp.scaled2=0.0)
            pData = rbind.data.frame(pData, dfrow, stringsAsFactors = F)
          }
        }
        
        #pData$id = factor(pData$id, levels = allIDs)
        #pData$features.plot = factor(pData$features.plot, levels=allFeatures)
        
        plotData[[plotName]] = pData
        
        
        if (is.null(allUnscaledValues))
        {
          allUnscaledValues = data.frame(plotName=pData[ pData$features.plot==featureName, ]$avg.exp)
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
          allUnscaledValues[["plotName"]] = NULL
          
        } else {
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
        }
        
      }

      allUnscaledValues$rnames = as.numeric(rownames(allUnscaledValues))
      allUnscaledLong = allUnscaledValues %>% gather(Type, Value, names(plotData))
      allUnscaledLong$Value = scale(allUnscaledLong$Value)
      allScaledValues = allUnscaledLong %>% spread(Type, Value) %>% arrange( order(rnames))
      
      for (plotName in names(plotData))
      {

        pData = plotData[[plotName]]
               
        origData = pData[pData$features.plot==featureName, ]
        pData[pData$features.plot==featureName, "avg.exp.scaled2"] = pData[pData$features.plot==featureName, "avg.exp"]
        #pData$idn=as.numeric(pData$id)

        plotData[[plotName]] = pData
      }


      
    }

  } else if (scale.by == "ALL") {

    print("SCALING BY ALL")


    combinedDataDF = data.frame()

    for (plotName in names(plotData))
    {
      pData = plotData[[plotName]]
      pData["plotpart"] = plotName

      combinedDataDF = rbind(combinedDataDF, pData)
    }

    combinedDataDF[, c("avg.exp.scaled2")] = scale(combinedDataDF$avg.exp)


    # reorder
    originalGroups = scobj@meta.data[,group.by]
    if (is.factor(originalGroups))
    {
      originalSorting = levels(originalGroups)
      combinedDataDF$id = factor(combinedDataDF$id, levels = originalSorting)
    } else {
      combinedDataDF$id = factor(combinedDataDF$id, levels = gtools::mixedsort(as.character(unique(combinedDataDF$id))))
    }


    for (plotName in names(plotData))
    {
      subDF = combinedDataDF[combinedDataDF$plotpart == plotName, c("features.plot", "id", "avg.exp.scaled2", "plotpart")]

      pData = plotData[[plotName]]

      pData = merge(pData, subDF, by.x=c("features.plot", "id"), by.y=c("features.plot", "id"))

      #pData$avg.exp.scaled2 = subDF[pData$features.plot,]$avg.exp.scaled2
      #pData$id = plotData$id
      #pData$idn= as.numeric(pData$id)
      #pData$features.plot = factor(pData$features.plot, levels=unique(pData$features.plot))
      
      plotData[[plotName]] = pData

    }
      
  } else {
    stopifnot(FALSE)
  }


  for (plotName in names(plotData))
  {
    pData = plotData[[plotName]]

    pData$id = factor(pData$id, levels=idLevels)
    pData$idn = as.numeric(pData$id)
    pData$features.plot = factor(pData$features.plot, levels=featureLevels)

    plotData[[plotName]] = pData
  }


  # prepare ID DF
  idDF = data.frame()
  for (plotName in names(plotData))
  {
    pData = plotData[[plotName]]
    idDF = rbind(idDF, pData[, c("id", "idn")])
  }

  rownames(idDF) = NULL
  idDF = unique(idDF)

  if (!is.null(scale.by))
  {
    title.expression = paste("Avg. Expression\n(scaled by ", scale.by, ")", sep="")
  } else {
    title.expression = "Avg. Expression"
  }

  if (abundance.perelem)
  {
    title.cellabundance = "Cell Abundance\n(per condition)"
  } else {
    title.cellabundance = "Cell Abundance\n(all object cells)"
  }

  plotList = list()
  
  for (plotName in names(plotData))
  {
    pData = plotData[[plotName]]

    print(head(pData))
    print(head(ctFractions[[plotName]]))

    pData2 = merge(x=pData,y=ctFractions[[plotName]],by.x="id", by.y="Var1",all.x=TRUE)
    pData <-pData2 %>% dplyr::mutate(featuren=as.numeric(features.plot), percn=100*perc)  
   
    pData[pData$avg.exp.scaled2>col.max, 'avg.exp.scaled2'] = col.max
    pData[pData$avg.exp.scaled2<col.min, 'avg.exp.scaled2'] = col.min

    print(pData)

    minFeatureN = min(pData$featuren)
    maxFeatureN = max(pData$featuren)

    fillLimits = c(0, ceiling(fillLimitsMax*10)*10)
        
    plotElem <- ggplot2::ggplot(pData) +
      ggplot2::scale_x_continuous(breaks=pData$featuren, labels=pData$features.plot) +
      ggplot2::scale_y_continuous(breaks=idDF$idn, labels=idDF$id, limits = c(min(idDF$idn)-0.6, max(idDF$idn)+0.6))+
      ggplot2::geom_rect(ggplot2::aes(xmin=minFeatureN-.5, xmax=maxFeatureN+.5, ymin = idn-0.5, ymax = idn+0.5, fill=percn), alpha = 0.4, linetype="blank") +
      ggplot2::scale_fill_distiller(palette='Spectral', limits = fillLimits)+
      ggplot2::scale_size_continuous(range = c(0, 10))+
      ggplot2::geom_point(ggplot2::aes(x=featuren, y=idn, colour = avg.exp.scaled2, size = pct.exp)) +
      ggplot2::scale_color_gradient2(limits=c(col.min, col.max), low = cols[1], mid=cols[2], high = cols[3])+
      ggplot2::guides(color=ggplot2::guide_colourbar(title=title.expression, order = 1),
      size=ggplot2::guide_legend(title="Percent Expressing", order = 2),
      fill=ggplot2::guide_colourbar(title=title.cellabundance, order = 3))
    
    plotElem = plotElem + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), legend.position = "none", panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line())

    if (rotate.x)
    {
      plotElem = plotElem + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    }

    plotList[[plotName]] = plotElem
  }
    
  print("Final plot checks")

  plot_list=list()
  for (plotName in names(plotList))
  {
    
    plotElem = plotList[[plotName]]      
    
    print("descr label")
    pe = make_descr_label(plotElem, plotName) #plotElems[[plotName]]$label
  
    plot_list[[plotName]] = pe
 
  }

  print("Preparing Legend")
  
  # extract a legend that is laid out horizontally #
  legend_b <- ggpubr::get_legend(
    plotElem + 
      ggplot2::guides(color = ggplot2::guide_colorbar(title = title.expression, direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  )
  
  title <- cowplot::ggdraw() + cowplot::draw_label(title, fontface='bold')
  print("Preparing plot")

  dataList = list()
  for (i in 1:length(plot_list))
  {
    dataList[[i]] = plot_list[[i]]$data
  }
  ap=combine_plot_grid_list(
    plotlist = plot_list,
    labels = NULL,
    nrow=1,
    align = "v", axis="bt"
  )

  finalPlotList = list(title, ap, legend_b)

  fp = combine_plot_grid_list(plotlist = finalPlotList, ncol = 1, rel_heights = c(.05, 1, 0.2) )
  fp$data = dataList

  return(fp)
}





makeSideBySideDotPlot = function(scobj, plotElems, featureGenes = c(""), group.by="cellnames_manual", col.min = -3, col.max = 3, cols = c("grey", "blue"), title="", scaled=T, scale.by="GROUP", rotate.x=F, abundance.perelem=FALSE)
{

  print("SBS DP")
  print(scale.by)

  stopifnot(is.null(scale.by) ||scale.by %in% c("GROUP", "FEATURE", "ALL", "GLOBAL"))

  use.assay="data"

  if (!is.null(scale.by) && scale.by == "GLOBAL")
  {
    use.assay="scale.data"
  }

  plot_list = list()
  plot_orig_list = list()
  allFoundFeatures = c()
  
  ctFractions = list()
  
  

  scFullTable = table(scobj[[group.by]])
  scFullDf = as.data.frame(scFullTable)

  fillLimitsMax = 0
  elemCount = 0
  for (plotName in names(plotElems))
  {
    plotData = plotElems[[plotName]]
    plotCells = plotData$cells
    plotDescr = plotData$label
    
    scobj_subset = subset(scobj, cells=plotCells)
    
    plotElem_orig = DotPlot(scobj_subset, features=featureGenes, group.by = group.by, col.min = col.min, col.max=col.max, cols = cols, scale=FALSE) # scale is always false! scale.data => already scaled

    scTable = table(scobj_subset[[group.by]])
    scDf = as.data.frame(scTable)

    
    if (abundance.perelem)
    {
      scDf$perc = scDf$Freq / sum(scDf$Freq) # per elem abundance!
    } else {
      scDf$perc = scDf$Freq / sum(scFullDf$Freq) # total abundance!
    }

    fillLimitsMax = max(c(fillLimitsMax,max(scDf$perc)))

    ctFractions[[plotName]] = scDf   
    scobj_subset = NULL
    
    allFoundFeatures = unique(c(allFoundFeatures, as.character(plotElem_orig$data$features.plot)))
    
    plotElem_orig = plotElem_orig + scale_color_gradient(limits=c(col.min, col.max), low = cols[1], high = cols[2])
    plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")
    
    plot_orig_list[[plotName]] = plotElem_orig
  }  
  
  # initialize Values
  for (plotName in names(plot_orig_list))
  {
      plot_orig_list[[plotName]]$data$avg.exp.scaled2 = 0
      plot_orig_list[[plotName]]$data$id = as.character(plot_orig_list[[plotName]]$data$id)
      plot_orig_list[[plotName]]$data$features.plot = as.character(plot_orig_list[[plotName]]$data$features.plot)
  }
  
  allIDs = c()
  allFeatures = c()
  for (plotName in names(plot_orig_list))
  {
      allIDs = c(allIDs, plot_orig_list[[plotName]]$data$id)
      allFeatures = c(allFeatures, plot_orig_list[[plotName]]$data$features.plot)
  }
  allIDs = unique(allIDs)
  allFeatures = unique(allFeatures)

  if ((is.null(scale.by) || scale.by == "GLOBAL"))
  {

    for (plotName in names(plot_orig_list))
    {

      plotElem_orig = plot_orig_list[[plotName]]
                
      # https://github.com/satijalab/seurat/issues/2798
      plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression

      plotElem_orig$data$id = factor(plotElem_orig$data$id, levels = mixedsort(as.character(unique(plotElem_orig$data$id))))
      plotElem_orig$data$idn= as.numeric(plotElem_orig$data$id)
      plotElem_orig$data$features.plot = factor(plotElem_orig$data$features.plot, levels=unique(plotElem_orig$data$features.plot))
      
      #just copy old values
      plotElem_orig$data[, "avg.exp.scaled2"] = plotElem_orig$data$avg.exp

      plot_orig_list[[plotName]] = plotElem_orig
    }

  } else if (scale.by == "FEATURE")
  {
    print("SCALING BY FEATURE")
    # calculate avg.exp.scaled2 for each feature
    for (featureName in allFoundFeatures)
    {
      print(featureName)
      allUnscaledValues = NULL
      for (plotName in names(plot_orig_list))
      {
        pData = plot_orig_list[[plotName]]$data
        missingCelltypes = setdiff(allIDs, unique(pData$id))
                
        for (celltype in missingCelltypes)
        {
          for (feature in allFeatures)
          {
            dfrow = data.frame(avg.exp=0.0,pct.exp=0.0,features.plot=feature, id=celltype, avg.exp.scaled=0.0, avg.exp.scaled2=0.0)
            pData = rbind.data.frame(pData, dfrow, stringsAsFactors = F)
            
          }
        }
        
        pData$id = factor(pData$id, levels = allIDs)
        pData$features.plot = factor(pData$features.plot, levels=allFeatures)
        

        plot_orig_list[[plotName]]$data = pData
        
        
        if (is.null(allUnscaledValues))
        {
          allUnscaledValues = data.frame(plotName=pData[ pData$features.plot==featureName, ]$avg.exp)
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
          allUnscaledValues[["plotName"]] = NULL
          
        } else {
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
        }
        
      }

      allUnscaledValues$rnames = as.numeric(rownames(allUnscaledValues))
      allUnscaledLong = allUnscaledValues %>% gather(Type, Value, names(plot_orig_list))
      allUnscaledLong$Value = scale(allUnscaledLong$Value)
      allScaledValues = allUnscaledLong %>% spread(Type, Value) %>% arrange( order(rnames))
      
      for (plotName in names(plot_orig_list))
      {

        plotElem_orig = plot_orig_list[[plotName]]
        pData = plotElem_orig$data
        
        # https://github.com/satijalab/seurat/issues/2798
        plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
        
        origData = plotElem_orig$data[plotElem_orig$data$features.plot==featureName, ]
        plotElem_orig$data[plotElem_orig$data$features.plot==featureName, "avg.exp.scaled2"] = plotElem_orig$data[plotElem_orig$data$features.plot==featureName, "avg.exp"]
        plotElem_orig$data$idn=as.numeric(plotElem_orig$data$id)

        plot_orig_list[[plotName]] = plotElem_orig
      }


      
    }

  } else if (scale.by == "ALL") {

    combinedDataDF = data.frame()

    for (plotName in names(plot_orig_list))
    {

      plotElem_orig = plot_orig_list[[plotName]]
                
      # https://github.com/satijalab/seurat/issues/2798
      plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression

      storeDF = plotElem_orig$data
      storeDF["plotpart"] = plotName

      combinedDataDF = rbind(combinedDataDF, storeDF)
    }


    combinedDataDF[, "avg.exp.scaled2"] = scale(combinedDataDF$avg.exp)

    # reorder
    originalGroups = scobj@meta.data[,group.by]
    if (is.factor(originalGroups))
    {
      originalSorting = levels(originalGroups)
      combinedDataDF$id = factor(combinedDataDF$id, levels = originalSorting)
    } else {
      combinedDataDF$id = factor(combinedDataDF$id, levels = mixedsort(as.character(unique(combinedDataDF$id))))
    }


    for (plotName in names(plot_orig_list))
    {

      plotData = combinedDataDF[combinedDataDF$plotpart == plotName,]
      plotElem_orig = plot_orig_list[[plotName]]
      plotElem_orig$data$avg.exp.scaled2 = plotData$avg.exp.scaled2

      plotElem_orig$data$id = plotData$id
      plotElem_orig$data$idn= as.numeric(plotData$id)
      plotElem_orig$data$features.plot = factor(plotElem_orig$data$features.plot, levels=unique(plotElem_orig$data$features.plot))
      
      plot_orig_list[[plotName]] = plotElem_orig

    }
      
  } else {
    stopifnot(FALSE)
  }



  # prepare ID DF
  idDF = data.frame()
  for (plotName in names(plot_orig_list))
  {
    plotElem_orig = plot_orig_list[[plotName]]
    idDF = rbind(idDF, plotElem_orig$data[, c("id", "idn")])
  }
  rownames(idDF) = NULL
  idDF = unique(idDF)

  if (!is.null(scale.by))
  {
    title.expression = paste("Avg. Expression\n(scaled by ", scale.by, ")", sep="")
  } else {
    title.expression = "Avg. Expression"
  }

  if (abundance.perelem)
  {
    title.cellabundance = "Cell Abundance (per condition)"
  } else {
    title.cellabundance = "Cell Abundance (all object cells)"
  }
  
  for (plotName in names(plot_orig_list))
  {
    plotElem_orig = plot_orig_list[[plotName]]
    pData = plotElem_orig$data

    pData2 = merge(x=pData,y=ctFractions[[plotName]],by.x="id", by.y=group.by,all.x=TRUE)
    plotElem_orig$data <-pData2 %>% mutate(featuren=as.numeric(features.plot), percn=100*perc)  
    
    plotElem_orig$data[plotElem_orig$data$avg.exp.scaled2>col.max, 'avg.exp.scaled2'] = col.max
    plotElem_orig$data[plotElem_orig$data$avg.exp.scaled2<col.min, 'avg.exp.scaled2'] = col.min

    pData = plotElem_orig$data

    fillLimits = c(0, ceiling(fillLimitsMax*10)*10)
    
    plotElem_orig <- ggplot(pData) +
      scale_x_continuous(breaks=plotElem_orig$data$featuren, labels=plotElem_orig$data$features.plot) +
      scale_y_continuous(breaks=idDF$idn, labels=idDF$id, limits = c(min(idDF$idn)-0.6, max(idDF$idn)+0.6))+
      geom_rect(aes(xmin=featuren-.5, xmax=featuren+.5, ymin = idn-0.5, ymax = idn+0.5, fill=percn), alpha = 0.4, linetype="blank") +
      scale_fill_distiller(palette='Spectral', limits = fillLimits)+
      scale_size_continuous(range = c(0, 10))+
      geom_point(aes(x=featuren, y=idn, colour = avg.exp.scaled2, size = pct.exp)) +
      scale_color_gradient2(limits=c(col.min, col.max), low = cols[1], mid=cols[2], high = cols[3])+
      guides(color=guide_colourbar(title=title.expression, order = 1),
      size=guide_legend(title="Percent Expressing", order = 2),
      fill=guide_colourbar(title=title.cellabundance, order = 3))
    
    plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line())

    if (rotate.x)
    {
      plotElem_orig = plotElem_orig + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }

    plot_orig_list[[plotName]] = plotElem_orig
  }
    
  print("Final plot checks")

  plot_list=list()
  for (plotName in names(plot_orig_list))
  {
    
    plotElem_orig = plot_orig_list[[plotName]]      
    plotElem = plotElem_orig
    
    print("descr label")
    pe = make_descr_label(plotElem, plotElems[[plotName]]$label)
  
    plot_list[[plotName]] = pe
 
  }

  print("Preparing Legend")
  
  # extract a legend that is laid out horizontally #
  legend_b <- get_legend(
    plotElem_orig + 
      guides(color = guide_colorbar(title = title.expression, direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      theme(legend.position = "bottom", legend.box = "vertical")
  )
  
  title <- ggdraw() + draw_label(title, fontface='bold')
  print("Preparing plot")

  dataList = list()
  for (i in 1:length(plot_list))
  {
    dataList[[i]] = plot_list[[i]]$data
  }
  ap=combine_plot_grid_list(
    plotlist = plot_list,
    labels = NULL,
    nrow=1,
    align = "v", axis="bt"
  )

  finalPlotList = list(title, ap, legend_b)

  fp = combine_plot_grid_list(plotlist = finalPlotList, ncol = 1, rel_heights = c(.05, 1, 0.2) )
  fp$data = dataList

  return(fp)
}















makeTopDownDotplot = function(scobj, plotElems, featureGenes = c(""), group.by="cellnames_manual", col.min = -3, col.max = 3, cols = c("grey", "blue"), title="", scaled=T, scale.by="GROUP", rotate.x=FALSE)
{

  stopifnot(scale.by %in% c("GROUP", "ALL"))

  plot_list = list()
  plot_orig_list = list()
  allFoundFeatures = c()
  
  ctFractions = data.frame()

  scFullTable = table(scobj[[group.by]])
  scFullDf = as.data.frame(scFullTable)

  fillLimitsMax = 0
  elemCount = 0
  for (plotName in names(plotElems))
  {
    print(plotName)
    plotData = plotElems[[plotName]]
    plotCells = plotData$cells
    plotDescr = plotData$label
    
    scobj_subset = subset(scobj, cells=plotCells)
    plotElem_orig = DotPlot(scobj_subset, features=featureGenes, group.by = group.by, col.min = col.min, col.max=col.max, cols = cols)
    scTable = table(scobj_subset[[group.by]])
    scDf = as.data.frame(scTable)
    
    scDf$perc = scDf$Freq / sum(scFullDf$Freq) # total abundance!
    scDf$plotpart = plotName

    fillLimitsMax = max(c(fillLimitsMax,max(scDf$perc)))
    ctFractions = rbind(ctFractions, scDf)
    
    scobj_subset = NULL
    
    allFoundFeatures = unique(c(allFoundFeatures, as.character(plotElem_orig$data$features.plot)))
    
    plotElem_orig = plotElem_orig + scale_color_gradient(limits=c(col.min, col.max), low = cols[1], high = cols[2])
    plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")
    
    plot_orig_list[[plotName]] = plotElem_orig
  }
  
  
    
    allIDs = c()
    allFeatures = c()
    for (plotName in names(plot_orig_list))
    {
        plot_orig_list[[plotName]]$data$id = as.character(plot_orig_list[[plotName]]$data$id)
        plot_orig_list[[plotName]]$data$features.plot = as.character(plot_orig_list[[plotName]]$data$features.plot)
    

        allIDs = c(allIDs, plot_orig_list[[plotName]]$data$id)
        allFeatures = c(allFeatures, plot_orig_list[[plotName]]$data$features.plot)
    }
    allIDs = unique(allIDs)
    allFeatures = unique(allFeatures)
    
    print(allIDs)
    print(allFeatures)
    
    if ((scale.by == "ALL") && (scaled == T)) {

        combinedDataDF = data.frame()

        for (plotName in names(plot_orig_list))
        {

          plotElem_orig = plot_orig_list[[plotName]]
                   
          # https://github.com/satijalab/seurat/issues/2798
          plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression



          storeDF = plotElem_orig$data
          storeDF["plotpart"] = plotName

          storeDF$id = factor(storeDF$id, levels = mixedsort(as.character(unique(storeDF$id))))
          storeDF$idn= as.numeric(storeDF$id)

          combinedDataDF = rbind(combinedDataDF, storeDF)
        }

        combinedDataDF[, "avg.exp.scaled2"] = scale(combinedDataDF$avg.exp)
        combinedDataDF$id = factor(combinedDataDF$id, levels = mixedsort(as.character(unique(combinedDataDF$id))))
        combinedDataDF$idn= as.numeric(combinedDataDF$id)
        combinedDataDF$features.plot = factor(combinedDataDF$features.plot, levels=unique(combinedDataDF$features.plot))

       
    } else if ((scale.by == "GROUP") || ((scale.by == "ALL") && (scaled == F))) {
        combinedDataDF = data.frame()

        for (plotName in names(plot_orig_list))
        {

          plotElem_orig = plot_orig_list[[plotName]]
                   
          # https://github.com/satijalab/seurat/issues/2798
          plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression

          plotElem_orig$data["plotpart"] = plotName

          plotElem_orig$data$id = factor(plotElem_orig$data$id, levels = mixedsort(as.character(unique(plotElem_orig$data$id))))
          plotElem_orig$data$idn= as.numeric(plotElem_orig$data$id)
          plotElem_orig$data$features.plot = factor(plotElem_orig$data$features.plot, levels=unique(plotElem_orig$data$features.plot))
          

          if (scaled)
          {
            plotElem_orig$data[, "avg.exp.scaled2"] = scale(plotElem_orig$data$avg.exp)
          } else {
            plotElem_orig$data[, "avg.exp.scaled2"] = plotElem_orig$data$avg.exp
          }

          combinedDataDF = rbind(combinedDataDF, plotElem_orig)
        }
    } else {
      stopifnot(FALSE)
    }

    print(combinedDataDF)




    pData2 = merge(x=combinedDataDF,y=ctFractions,by.x=c("id", "plotpart"), by.y=c("Var1", "plotpart"),all.x=TRUE)
    combinedDataDF <-pData2 %>% mutate(featuren=as.numeric(features.plot), percn=100*perc)
    
    combinedDataDF[combinedDataDF$avg.exp.scaled2>col.max, 'avg.exp.scaled2'] = col.max
    combinedDataDF[combinedDataDF$avg.exp.scaled2<col.min, 'avg.exp.scaled2'] = col.min


    if (scaled==T)
    {
      exprTitle = paste("Avg. Expression (scaled by ", scale.by, ")", sep="")
    } else {
      exprTitle = "Avg. Expression"
    }

    combinedDataDF$id = paste(combinedDataDF$id, combinedDataDF$plotpart, sep=" ")

    combinedDataDF$id = factor(combinedDataDF$id, levels = mixedsort(as.character(unique(combinedDataDF$id))))
    combinedDataDF$idn= as.numeric(combinedDataDF$id)

    idDF = combinedDataDF[, c("id", "idn")]
    print("idDF")
    rownames(idDF) = NULL
    idDF = unique(idDF)

    print(plotName)
    print(idDF)
    print(fillLimitsMax)
    print(combinedDataDF)







    fillLimits = c(0, ceiling(fillLimitsMax*10)*10)
        
  mainPlot <- ggplot(combinedDataDF) +
    scale_x_continuous(breaks=combinedDataDF$featuren, labels=combinedDataDF$features.plot) +
    scale_y_continuous(breaks=idDF$idn, labels=idDF$id, limits = c(min(idDF$idn)-0.6, max(idDF$idn)+0.6))+
    geom_rect(aes(xmin=featuren-.5, xmax=featuren+.5, ymin = idn-0.5, ymax = idn+0.5, fill=percn), alpha = 0.4, linetype="blank") +
    scale_fill_distiller(palette='Spectral', limits = fillLimits)+
    scale_size_continuous(range = c(0, 10))+
    geom_point(aes(x=featuren, y=idn, colour = avg.exp.scaled2, size = pct.exp)) +
    scale_color_gradient2(limits=c(col.min, col.max), low = cols[1], mid=cols[2], high = cols[3])+
    guides(color=guide_colourbar(title=exprTitle, order = 1),
    size=guide_legend(title="Percent Expressing", order = 2),
    fill=guide_colourbar(title="Cell Abundance (selected cells)", order = 3))
  
  mainPlot = mainPlot + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line())
  

  if (rotate.x)
  {
    mainPlot = mainPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }

  legendDescr = 'Average Expression'
  if (scaled)
  {
    legendDescr = 'Average Scaled Expression'  
  }
  
  # extract a legend that is laid out horizontally #
  legend_b <- get_legend(
    mainPlot + 
      guides(color = guide_colorbar(title = legendDescr, direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      theme(legend.position = "bottom")
  )
  
  title <- ggdraw() + draw_label(title, fontface='bold')

  print("Combining plots")
  fp = combine_plot_grid_list(plotlist=list(title, mainPlot, legend_b), ncol = 1, rel_heights = c(.05, 1, 0.1) )


  return(fp)
}







splitFeaturePlot = function(obj, feature, split.by, title, filename,limits=c(-1,1), reduction="umap", low="lightgrey", high="blue", mid="white", mirrorLimits=TRUE)
{
  abLimit = max(abs(limits))

  if (mirrorLimits)
  {
    limits = c(-abLimit, abLimit)
  }

  print(paste("limits", limits))
  
  pds = FeaturePlot(obj, features = feature, reduction = reduction, split.by=split.by, combine=F,min.cutoff=limits[1], max.cutoff=limits[2],order=T)
  pds[[1]] = pds[[1]] + ggtitle(NULL)

  print(paste(length(pds)))
  pdsRange = c(1:(length(pds)))

  for (i in pdsRange)
  {
    print(i)

    if (is.null(mid))
    {
      cGradient = scale_color_gradient(limits = limits, low = low, high = high)
    } else {
      cGradient = scale_color_gradient2(limits = limits, low = low, high = high, mid=mid)
    }


    pds[[i]] = pds[[i]] + cGradient + theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0), legend.position = c(0.05,0.1), legend.key.height = unit(0.25, 'cm'), legend.text = element_text(size=8  ))+labs(x="", y="")  
  }

  #pds[[length(pds)]] = (pds[[length(pds)]] + scale_color_gradient2(limits = limits, low = low, high = high, mid=mid )) + theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0), legend.position = c(0.05,0.15))+labs(x="", y="")


  prow =combine_plot_grid_list(plotlist=pds, label_x = "a", ncol=length(pds), align="hv")
  # now add the title
  title <- ggdraw() + 
    draw_label(
      title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )

  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
  fplot = plot_grid(
    title, prow,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )

  prow = combine_plot_grid_list(plotlist=pds, label_x = "a", ncol=length(pds), align="hv")
  fplot= combine_plot_grid_list(plotlist=list(title, prow), ncol=1, rel_heights=c(0.1, 1))

  if (!is.null(filename))
  {
    save_plot(fplot, filename, fig.width=length(pds) *6 + 2, fig.height=6)

  }
  return(fplot)
}





VlnBoxPlot = function( obj.sc, gene, group.by="idents", split.by=NULL, pt.size=0, assay="RNA", min.threshold=3, col=NULL, per.sample=FALSE)
{

    remainingCells = NULL

    if (is.null(split.by))
    {

        cellGroupDF = obj.sc[[ group.by ]]

        countDF = table(cellGroupDF)
        countDF = countDF[countDF >= min.threshold]

        remainingClusters = rownames(countDF)
        remainingCells = rownames(cellGroupDF)[ cellGroupDF[[group.by]] %in% remainingClusters ]

    } else {

        cellGroupDF = obj.sc[[ c(group.by, split.by) ]]

        countDF = as.data.frame(table(cellGroupDF))
        countDF = countDF[countDF$Freq >= min.threshold,]

        if (per.sample)
        {

          snames = unique(obj.sc[[ split.by ]][[split.by]])
          remainingCells = c()

          for (sname in snames)
          {
            snameDF = countDF[ countDF[split.by] == sname,]
            snameDF = snameDF[snameDF$Freq >= min.threshold,]

            sCells = rownames(obj.sc[[split.by]])[ obj.sc[[split.by]][[split.by]] == sname ]
            gCells = rownames(obj.sc[[group.by]])[ obj.sc[[group.by]][[group.by]] %in% snameDF[[group.by]]]

            remainingCells = c(remainingCells, intersect(sCells, gCells))
          }

        } else {
          remainingClusters = countDF[duplicated(countDF[group.by]),c(group.by)]

          cellGroupDF = obj.sc[[ group.by ]]
          remainingCells = rownames(cellGroupDF)[ cellGroupDF[[group.by]] %in% remainingClusters ]
        }



    }
    

    print(head(remainingCells))
    obj.sc.subset = subset(obj.sc, cells=remainingCells)
    print(obj.sc.subset)


    if (is.null(split.by))
    {

      ucols = NULL
      if (!is.null(col))
      {
        ncolors = length(unique(obj.sc[[group.by]][[group.by]]))
        ucols = rep(col, ncolors)
      }
      

    plots=VlnPlot(obj.sc.subset, gene, group.by=group.by, split.by=split.by, pt.size=pt.size, assay=assay, col=ucols)
    plots = plots + geom_boxplot(color="grey", alpha=0.4) + stat_summary(fun=mean, geom="point", color="black", size=4)

    } else {

    plots=VlnPlot(obj.sc.subset, gene, group.by=group.by, split.by=split.by, pt.size=pt.size, assay=assay, col=col)
    plots = plots + geom_boxplot(color="grey", alpha=0.4, position =position_dodge(width = 0.9)) + stat_summary(fun=mean, geom="point", aes(group=split), position=position_dodge(.9), color="black", size=4)

    }
    

    return(plots)
}


SplitVlnBoxPlot = function( obj.sc, gene, group.by="idents", split.by=NULL, pt.size=0, assay="RNA", min.threshold=3, col=NULL, per.sample=TRUE)
{

    splitValues = unique(obj.sc[[split.by]][[split.by]])

    vplots = list()

    for (sname in splitValues)
    {
      print(sname)
      subsetcellnames = rownames(obj.sc[[split.by]])[ obj.sc[[split.by]][split.by] == sname ]

      if (is.null(col))
      {
        scol = NULL
      } else {
        scol = col[[sname]]
      }

      p=VlnBoxPlot( subset(obj.sc, cells=subsetcellnames), gene, group.by=group.by, split.by=NULL, pt.size=pt.size, assay=assay, min.threshold=min.threshold, col=scol, per.sample=per.sample)
      p = p + ggtitle(paste(gene, "-", sname))
      vplots[[sname]] = p
    }

    ps = combine_plot_grid_list(plotlist=vplots, nrow=1)

    return(ps)
}




plotPValueViolentBoxPlot = function( obj.sc, feature, group.by, split.by=NULL, split.values=NULL, dsrCols=NULL, onelineLabel=FALSE, dot.size=0, min.threshold=3, yStepIncrease=0.25, override=FALSE)
{
  # dsrCols = list("Thrombus"="grey", "Blood="red")
  


  remainingCells = NULL

  if (is.null(split.by))
  {

      cellGroupDF = obj.sc[[ group.by ]]

      countDF = table(cellGroupDF)
      countDF = countDF[countDF >= min.threshold]

      remainingClusters = rownames(countDF)
      remainingCells = rownames(cellGroupDF)[ cellGroupDF[[group.by]] %in% remainingClusters ]

  } else {

    cellGroupDF = obj.sc[[ c(group.by, split.by) ]]

    countDF = as.data.frame(table(cellGroupDF))

    groupValues = unique(cellGroupDF[, group.by])
    splitValues = unique(cellGroupDF[, split.by])
    remainingCells = c() 
    for (gv in groupValues)
    {
      
      keepGroup = TRUE
      for (sv in splitValues)
      {
        counts = (countDF %>% dplyr::filter(.data[[group.by]] == gv) %>% dplyr::filter(.data[[split.by]] == sv))
        
        countsValue = 0
        if (dim(counts)[1] > 0)
        {
          countsValue = counts$Freq
        }

        if ((override==FALSE) && (countsValue < min.threshold))
        {
          keepGroup = FALSE
          print(paste(gv, sv, countsValue))
        }
      }

      if (keepGroup)
      {
        sCells = cellIDForClusters(obj.sc, split.by, splitValues)
        gCells = cellIDForClusters(obj.sc, group.by, gv)
        remainingCells = c(remainingCells, intersect(sCells, gCells))
      } else {
        print(paste("Removing group", gv))
      }
    }

  }


  print(paste("Existing Cells", length(colnames(obj.sc))))
  print(paste("New Cells", length(remainingCells)))

  obj.sc = subset(obj.sc, cells=remainingCells)

  if ((is.null(split.values) && !is.null(split.by)))
  {
    split.values = as.character(unique(obj.sc@meta.data[[split.by]]))
    print(split.values)
  }

  if (feature %in% colnames(obj.sc@meta.data))
  {
    dataDF = obj.sc@meta.data[,c(feature, group.by, split.by)]
  } else {
    dataDF = obj.sc@meta.data[,c(group.by, split.by)]
    dataDF[[feature]] = obj.sc@assays$RNA@data[feature, rownames(dataDF)]
  }
  
  if (!is.null(split.by))
  {
    dataDF = dataDF[dataDF[[split.by]] %in% split.values, ]
  }
  

  if (!is.null(levels(dataDF[[group.by]])))
  {
    dataDF[[group.by]]= factor(dataDF[[group.by]], levels=intersect( levels(dataDF[[group.by]]), as.character(dataDF[[group.by]]) ))
  } else {
    dataDF[[group.by]] = as.factor(dataDF[[group.by]])
  }

  if ((!is.null(split.by)) && (is.null(levels(dataDF[[split.by]]))))
  {
    dataDF[[split.by]] = as.factor(dataDF[[split.by]])
  }


  keepAC = c()
  for (ac in unique(dataDF[[group.by]]))
  {
    subdf = dataDF %>% filter(.data[[group.by]] == ac)

    if (!is.null(split.by))
    {
      numDiffGroups = length(unique(subdf[[split.by]]))
    } else {
      numDiffGroups = 2
    }
    

    print(paste(ac, numDiffGroups))

    if (numDiffGroups > 1)
    {
      keepAC = c(keepAC, ac)
    }
  }

  stat.test = dataDF %>% filter((!!as.symbol(group.by)) %in% keepAC) %>% group_by_at(group.by, .add=T)

  if (!is.null(split.by))
  {
   stat.test = stat.test %>%
  pairwise_t_test(
      as.formula(paste(feature, " ~ ", split.by, sep="")), paired = FALSE, 
      p.adjust.method = "BH"
      )
  } else {

        print(head(as.data.frame(stat.test)))


    stat.test = as.data.frame(stat.test) %>%
  pairwise_t_test(
      as.formula(paste(feature, " ~ ", group.by, sep="")), paired = FALSE, 
      p.adjust.method = "BH"
      )

    print(stat.test)
  }
  


  if (onelineLabel)
  {
    stat.test$label = paste("n1=",stat.test$n1," ", "n2=",stat.test$n2,", p < ", formatC(stat.test$p.adj, format = "e", digits = 2), sep="")
  } else {
    stat.test$label = paste("n1=",stat.test$n1,",\n", "n2=",stat.test$n2,"\n","p < ", formatC(stat.test$p.adj, format = "e", digits = 2), sep="")
  }

  maxValue = max(dataDF[,c(feature)])
  minValue = min(dataDF[,c(feature)])


  if (!is.null(split.by))
  {

    groupValues = unique(dataDF[, group.by])
    splitValues = unique(dataDF[, split.by])

    remainingCells = c() 
    for (gv in groupValues)
    {
      
      keepGroup = TRUE
      for (sv in splitValues)
      {
        subsetExprDF = (dataDF %>% dplyr::filter(.data[[group.by]] == gv) %>% dplyr::filter(.data[[split.by]] == sv))

        featureSD = sd(subsetExprDF[ , c(feature)])
        print(paste(sv, featureSD))

        if (!is.na(featureSD) && (featureSD == 0))
        {
          dataDF[rownames(subsetExprDF)[1], c(feature)] = dataDF[rownames(subsetExprDF)[1], c(feature)] + 0.000001
        }
        
      }
    }
  }

  print("Creating Plot")

  if (is.null(split.by))
  {
    split.by=group.by
  }

  bxp <- ggplot(dataDF, aes_string(x=group.by, y=feature)) +
        geom_violin(aes_string(fill=split.by), position =position_dodge(width = 0.9), trim=TRUE)

  if (dot.size > 0)
  {
    bxp = bxp + geom_jitter(aes_string(fill=split.by), position =position_jitterdodge(dodge.width=0.9), size=dot.size)
  }

  bxp = bxp +
        geom_boxplot(aes_string(fill=split.by),color="grey", alpha=0.4, position =position_dodge(width = 0.9)) +
        theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1.0),
              panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),panel.background = element_blank())+labs(fill='Group')
        
  if (!is.null(dsrCols))
  {
    bxp = bxp + scale_fill_manual(values = dsrCols[names(dsrCols) %in% split.values])
  }

  stat.test = stat.test[stat.test$p.adj < 0.05,]

  stat.test <- stat.test %>% rstatix::add_xy_position(x = group.by, dodge = 0.8, step.increase=yStepIncrease)#, y.trans=function(x){x*1.1})
  bxp = bxp + stat_pvalue_manual(stat.test,  label = "label", tip.length = 0)+ylim(c(floor(minValue*1.1), ceiling(maxValue*1.2)))
  return(bxp)
}



makeUMAPPlot = function(obj.in, dim1, dim2, reduction="umap", xmin = -11,xmax = 15,ymin = -15,ymax = 10.5, group.by="cellnames", downsample=FALSE)
{
  DefaultAssay(obj.in) = "RNA"

  targetElemsDim1 = sort(unique(obj.in[[dim1]][[dim1]]))
  targetElemsDim2 = sort(unique(obj.in[[dim2]][[dim2]]))

  print(targetElemsDim1)
  print(targetElemsDim2)

  allplots = list()
  lastrealplot = NULL

  minCells = length(colnames(obj.in))

  if (downsample)
  {
    for (di1 in targetElemsDim1)
    {
      for (di2 in targetElemsDim2)
      {
        dimcells = intersect(
          cellIDForClusters(obj.in, dim1, di1),
          cellIDForClusters(obj.in, dim2, di2)
        )

        if (length(dimcells) > 10)
        {
          minCells = min(c(minCells, length(dimcells)))
        }
      }
    }

    print(paste("Downsampling", minCells))
  }

  for (di1 in targetElemsDim1)
  {

    for (di2 in targetElemsDim2)
    {

      pname = paste(di1, di2, sep="_")
      print(pname)

      dimcells = intersect(
        cellIDForClusters(obj.in, dim1, di1),
        cellIDForClusters(obj.in, dim2, di2)
      )

      if ((downsample) && (length(dimcells) > minCells))
      {
        dimcells = sample(dimcells, minCells)
      }

      print(length(dimcells))

      if (length(dimcells) == 0)
      {
        allplots[[pname]] = ggplot() + theme_void()
      } else {
        allplots[[pname]] = DimPlot(subset(obj.in, cells=dimcells), group.by=group.by)  + xlim(xmin, xmax) + ylim(ymin, ymax)+theme(legend.position = "none")+ggtitle(NULL)
        lastrealplot = pname
      }
      

    }

  }

  print("Finishing Plot")

  legend_b <- get_legend(
    allplots[[lastrealplot]] + 
      guides(color = guide_legend(nrow = 2, override.aes = list(size=10)))+
      theme(legend.position = "bottom")
  )

  print("Combining Plots")

  ap = combine_plot_grid_list(
    plotlist = allplots, align="hv",
    labels = names(allplots), ncol=length(targetElemsDim2)
  )


  titletext = paste("DimPlot", dim1, dim2, sep=" ")
  if (downsample)
  {
    titletext = paste(titletext, "(downsampled to", minCells, "cells)", sep=" ")
  }
  title <- ggdraw() + draw_label(titletext, fontface='bold')
  print("Combining Legend")
  fp = combine_plot_grid_list(plotlist=list(title, ap, legend_b), ncol = 1, rel_heights = c(.1, 1, .1))


  return(fp)
}



compareCellClusterAssignments = function (inobj1, inobj2, outname="cluster_comparison", inName1="pat", inName2="lib", return.matrix=FALSE)
{


in1CellClusters = inobj1@meta.data[,c(inName1)]
names(in1CellClusters) = rownames(inobj1@meta.data)

in2CellClusters = inobj2@meta.data[,c(inName2)]
names(in2CellClusters) = rownames(inobj2@meta.data)

in1Clusters = length(unique(in1CellClusters))
in2Clusters = length(unique(in2CellClusters))

clusterComparison = matrix(0, nrow=in1Clusters, ncol=in2Clusters)
rownames(clusterComparison) = unique(in1CellClusters)
colnames(clusterComparison) = unique(in2CellClusters)

for (cellname in intersect(names(in1CellClusters), names(in2CellClusters)))
{

  in1Cluster = in1CellClusters[cellname]
  in2Cluster = in2CellClusters[cellname]
  
  clusterComparison[in1Cluster, in2Cluster] = clusterComparison[in1Cluster, in2Cluster] + 1
}


melted_comparison <- melt(clusterComparison)
colnames(melted_comparison) = c(inName1, inName2, "value")

melted_comparison[, c(inName1)] = as.character(melted_comparison[, c(inName1)])
melted_comparison[, c(inName2)] = as.character(melted_comparison[, c(inName2)])

#melted_comparison[[inName1]] = melted_comparison[[inName1]]-1
#melted_comparison[[inName2]] = melted_comparison[[inName2]]-1

p=ggplot(melted_comparison, aes_string("x" = (inName1), "y" = (inName2), "col" = "value", "fill" = "value", "label" = "value")) +
  geom_tile() +
  geom_text(col = "black") +
  theme_minimal() +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red") +
  scale_color_gradient2(low = "white", mid = "yellow", high = "red")
save_plot(p, outname, 8,8)

if (return.matrix)
{
  return(clusterComparison)
}
}

#
##
### Annotation
##
#

annotateByCellnamePattern = function(obj.in, group.name, annotation, order=NULL, use.base=NULL)
{

cellList = colnames(obj.in)
featVec <- vector(mode="character", length=length(cellList))

if (is.null(use.base) || (use.base=="cells"))
{
  cellList = colnames(obj.in)
} else {
  cellList = obj.in@meta.data[, use.base]
}

names(featVec) = cellList

for (elem in annotation)
{
  annotName = elem$name
  annotSelector = elem$selector

  print(paste(annotName, is.vector(annotSelector) & length(annotSelector) == 1))

  if (is.vector(annotSelector) & length(annotSelector) == 1) { #string
    featVec[grep(x = cellList, pattern = annotSelector)] = annotName
  } else if (is.vector(annotSelector) & length(annotSelector) > 1) { #vector
    featVec[ annotSelector ] = annotName
  } else {
    print(paste("Ignored selector", annotName, annotSelector))
  }

}
obj.in = AddMetaData(obj.in, featVec, col.name=group.name)

if (!is.null(order))
{

obj.in@meta.data[,c(group.name)] = factor(obj.in@meta.data[,c(group.name)], levels=order)

}

print(unique(obj.in@meta.data[,c(group.name)]))

return(obj.in)
}



getClusterColors = function(seuratObj, group_by)
{

  # create a dimplot and extract colorscheme from there!
  p=DimPlot(seuratObj, group.by = group_by)
  g <- ggplot_build(p)

  g2c = unique(g$data[[1]][, c("group", "colour")])
  g2c = g2c[order(g2c$group),]

  color2label = data.frame(colours = g2c$colour, group = g2c$group,
              label = g$plot$scales$scales[[1]]$get_labels()) 


  colorlabel = color2label$colour
  names(colorlabel) = color2label$label
  
  return(colorlabel)

}





makePieTrieCounts = function( seuratObj, outpath, select_by, group_by, size.text = 10, repel.direction="x", max.overlaps=10, fig.height=16, use.palette=NULL)
{

  # pallette should be something rcolorbrewer understands, e.g. set2

outcheckpath = dirname(paste(outpath, ".tsv", sep=""))
print(outcheckpath)
print(dir.exists(outcheckpath))
if (!dir.exists(outcheckpath))
{
  dir.create(outcheckpath, recursive = TRUE)
}

countByManualCellname = getCellCountDF(seuratObj, prefix="", select_by = select_by, group_by=group_by, relative=TRUE, show.percent=T, outname=paste(outpath, "tsv", sep="."))

cbc_cname = melt(countByManualCellname)

cbc_cname$variable_factors = factor(cbc_cname$variable)
cbc_cname = cbc_cname[!cbc_cname$cluster=="Total",]
cbc_cname$cluster = factor(cbc_cname$cluster, levels=levels(seuratObj@meta.data[,select_by]))
cbc_cname$total = 100
cbc_cname$label = paste0(cbc_cname$cluster, " (", round(cbc_cname$value, 2), "% )")


if (!is.null(use.palette))
{

  if (length(use.palette)==1)
  {
    colourCount = length(unique(cbc_cname$cluster))
    colorPalette = colorRampPalette(brewer.pal(colourCount, use.pallete))(colourCount)
  } else {
    colorPalette = use.palette[levels(cbc_cname$cluster)]
  }

} else {

  colorlabel = getClusterColors(seuratObj, select_by)
  colorPalette = colorlabel[levels(cbc_cname$cluster)]
}

print(colorPalette)



numColumns = length(unique(cbc_cname$variable))

condplots = list()
plotvars = unique(cbc_cname$variable)



  for (plotcond in plotvars)
  {

    print(plotcond)
    subdf = cbc_cname[cbc_cname$variable == plotcond,]
    #print(subdf)

    data = subdf %>% mutate(ypos = cumsum(value)- 0.5*value )

    p = ggplot(data, aes(x="", y=value, fill=cluster)) +
        geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) +
        theme_void() + theme(legend.position="none") +
        scale_fill_manual(values=colorPalette) + ggtitle(plotcond)
        #        geom_text(aes(y = ypos, label = cluster), color = "white", size=6) +
        #        theme(legend.position="none") +

    if (!(plotcond == plotvars[length(plotvars)]))
    {
        p = p + theme(legend.position="none")
    }


    condplots[[plotcond]] = p
    
  }

  legend_b <- get_legend(
    condplots[[plotcond]] + 
      guides(color = guide_colorbar(title = "Group", direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      theme(legend.position = "bottom")
  )

  topPlot = combine_plot_grid_list(condplots, ncol=length(condplots))
  
  topLegend = combine_plot_grid_list(list(topPlot, legend_b), ncol=1)

  numcols = length(condplots)
  use_width = numcols * 6
  save_plot(topLegend, outpath, fig.height=12, fig.width=use_width)


}



#
##
### NICHENET
##
#

library("nichenetr")

runNicheNet = function(obj.in, outpath="./nichenet/", group.by="idents", use.assay="RNA", use.receivers=NULL, use.senders=NULL, use.condition="condition", use.condition.coi=c("sepsis"), use.condition.bg=c("control"), number.ligands = 30, nichenet.folder="../nichenet/", fig.width=10, fig.height=10, add.ligands=NULL, add.receptors=NULL, add.targets=NULL)
{

#obj.in = obj.integrated.patient_thr
#outpath="./nichenet/"
#group.by="idents"
#use.assay="RNA"
#fig.width=10
#fig.height=10
#use.condition.coi=c("sepsis")
#use.condition.bg=c("control")
#use.condition="condition"
#use.receivers=c(1)
#use.senders=c(0,2,3,5,6,8,9)
#nichenet.folder="../../nichenet/"
#number.ligands = 50

#relgenes = c("CD177", "PECAM1", "CD38")
#add.ligands=relgenes
#add.receptors=relgenes
#add.targets=relgenes


  allGroupElems = unique(obj.in@meta.data[, c(group.by)])
  if (is.null(use.receivers))
  {
    use.receivers = allGroupElems
  }

  if (is.null(use.senders))
  {
    use.senders = setdiff(allGroupElems, use.receivers)
  }

  for (recv in use.receivers)
  {
    stopifnot( recv %in%  allGroupElems)
  }

  for (sndr in use.senders)
  {
    stopifnot( sndr %in%  allGroupElems)
  }


  ligand_target_matrix = readRDS(paste(nichenet.folder, "ligand_target_matrix.rds", sep="/"))
  lr_network = readRDS(paste(nichenet.folder, "lr_network.rds", sep="/"))
  sig_network = readRDS(paste(nichenet.folder, "signaling_network.rds", sep="/"))
  gr_network = readRDS(paste(nichenet.folder, "gr_network.rds", sep="/"))
  ligand_tf_matrix = readRDS(paste(nichenet.folder, "ligand_tf_matrix.rds", sep="/"))

  weighted_networks = readRDS(paste(nichenet.folder, "weighted_networks.rds", sep="/"))
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))


  ligand_target_matrix = t(ligand_target_matrix)

  receiver = use.receivers
  sender_celltypes = use.senders

  DefaultAssay(obj.in) = use.assay

  expressed_genes_receiver = get_expressed_genes(receiver, obj.in, pct = 0.10, assay_oi=use.assay)
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

  list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, obj.in, 0.10, assay_oi=use.assay) # lapply to get the expressed genes of every sender cell type separately here
  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

  outfolder = paste(outpath, paste(receiver, collapse="_"), sep="/")

  if (!dir.exists(outfolder))
  {
    print(paste("Creating folder", outfolder))
    dir.create(outfolder, recursive = TRUE)
  }

  if (!is.null(add.receptors))
  {
    expressed_genes_receiver = c(expressed_genes_receiver, add.receptors)
  }

  if (!is.null(add.ligands))
  {
    expressed_genes_sender = c(expressed_genes_sender, add.ligands)
  }

  #
  ## 2
  #

  cells.recv = cellIDForClusters(obj.in, group.by, c(receiver))

  seurat_obj_receiver= subset(obj.in, cells = cells.recv)
  #seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["aggregate"]])

  oiCond = use.condition.coi
  refCond = use.condition.bg

  cells.oi = intersect(cellIDForClusters(obj.in, use.condition, oiCond), colnames(seurat_obj_receiver))
  cells.ref = intersect(cellIDForClusters(obj.in, use.condition, refCond), colnames(seurat_obj_receiver))

  conditionNames = paste(paste(oiCond, collapse="_"), paste(refCond, collapse="_"), sep="__")
    
  DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = cells.oi, ident.2 = cells.ref, min.pct = 0.10, logfc.threshold=0.1) %>% rownames_to_column("gene")

  print("DE Results Receiver")
  print(dim(DE_table_receiver))


  # these are the genes in the “receiver/target” cell population that are potentially affected by ligands expressed by interacting cells
  geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

  if (!is.null(add.targets))
  {
    geneset_oi = c(geneset_oi, add.targets)
  }

  #
  ## 3
  #

  print("lr_network from missing add.ligands")
  print(setdiff(add.ligands, unique(lr_network$from)))
  print("lr_network from missing add.receptors")
  print(setdiff(add.receptors, unique(lr_network$from)))
  print("lr_network from missing add.targets")
  print(setdiff(add.targets, unique(lr_network$from)))

  print("lr_network to missing add.ligands")
  print(setdiff(add.ligands, unique(lr_network$to)))
  print("lr_network to missing add.receptors")
  print(setdiff(add.receptors, unique(lr_network$to)))
  print("lr_network to missing add.targets")
  print(setdiff(add.targets, unique(lr_network$to)))

  if (!is.null(add.ligands))
  {
    expressed_genes_sender = c(expressed_genes_sender, add.ligands)
  }

  if (!is.null(add.receptors))
  {
    expressed_genes_receiver = c(expressed_genes_receiver, add.receptors)
  }

  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()

  #Define a set of potential ligands: these are ligands that are expressed by the “sender/niche” cell population and bind a (putative) receptor expressed by the “receiver/target” population
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)

  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

  print("Potential Receptors")
  print(expressed_receptors)

  print("Potential Ligands")
  print(potential_ligands)

  #
  ## 4 Perform NicheNet ligand activity analysis: rank the potential ligands based on the presence of their target genes in the gene set of interest
  #

  print("ligand_target_matrix rows missing add.ligands")
  print(setdiff(add.ligands, rownames(ligand_target_matrix)))
  print("ligand_target_matrix rows missing add.receptors")
  print(setdiff(add.receptors, rownames(ligand_target_matrix)))
  print("ligand_target_matrix rows missing add.targets")
  print(setdiff(add.targets, rownames(ligand_target_matrix)))

  print("ligand_target_matrix cols missing add.ligands")
  print(setdiff(add.ligands, colnames(ligand_target_matrix)))
  print("ligand_target_matrix cols missing add.receptors")
  print(setdiff(add.receptors, colnames(ligand_target_matrix)))
  print("ligand_target_matrix cols missing add.targets")
  print(setdiff(add.targets, colnames(ligand_target_matrix)))

  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands) 

  ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
  best_upstream_ligands = ligand_activities %>% top_n(number.ligands, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()


  if (!is.null(add.ligands))
  {
    best_upstream_ligands = unique(c(best_upstream_ligands, add.ligands))
  }

  print("Best Upstream Ligands")
  print(best_upstream_ligands)

  p=DotPlot(obj.in, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
  save_plot(p, paste(outfolder, paste("dotplot", conditionNames, sep="_"), sep="/"), fig.width=fig.width, fig.height=fig.height)



  #
  ## 5 Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
  #

  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows() %>% drop_na()
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.0)

  order_ligands_raw = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_ligands = order_ligands_raw %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()



  plt <- ggplot(active_ligand_target_links_df,aes( target, ligand,fill=weight))
  plt <- plt + geom_tile(colour = "white",  size=1.5) + scale_fill_gradient2(low = "white", high = "red") + scale_x_discrete(
    expand = expansion(mult = c(0,0)), guide = guide_axis(angle = 45),
    position = "top"
  )+ theme_minimal()
  save_plot(plt, paste(outfolder, paste("ligand_target2", conditionNames, sep="_"), sep="/"), fig.width=fig.width, fig.height=fig.height)



  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

  p_ligand_target = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  save_plot(p_ligand_target, paste(outfolder, paste("ligand_target", conditionNames, sep="_"), sep="/"), fig.width=fig.width, fig.height=fig.height)

  #
  ## Prepare the ligand activity matrix
  #

  ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

  print(order_ligands_raw)

  vis_ligand_pearson = ligand_pearson_matrix[intersect(order_ligands_raw, rownames(ligand_pearson_matrix)), ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
  save_plot(p_ligand_pearson, paste(outfolder, paste("ligand_pearson", conditionNames, sep="_"), sep="/"), fig.width=fig.width, fig.height=fig.height)




  #
  ## Receptors of top-ranked ligands
  #

  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

  print("weighted_networks_lr$from missing add.ligands")
  print(setdiff(add.ligands, weighted_networks_lr$from))
  print("weighted_networks_lr$from missing add.receptors")
  print(setdiff(add.receptors, weighted_networks_lr$from))
  print("weighted_networks_lr$from missing add.targets")
  print(setdiff(add.targets, weighted_networks_lr$from))

  print("weighted_networks_lr$to missing add.ligands")
  print(setdiff(add.ligands, weighted_networks_lr$to))
  print("weighted_networks_lr$to missing add.receptors")
  print(setdiff(add.receptors, weighted_networks_lr$to))
  print("weighted_networks_lr$to missing add.targets")
  print(setdiff(add.targets, weighted_networks_lr$to))

  if (!is.null(add.ligands))
  {
    best_upstream_ligands = unique(c(best_upstream_ligands, add.ligands))
  }
  if (!is.null(add.receptors))
  {
    best_upstream_receptors = unique(c(best_upstream_receptors, add.receptors))
  }


  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
      
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
  save_plot(p_ligand_receptor_network, paste(outfolder, paste("ligand_target_network", conditionNames, sep="_"), sep="/"), fig.width=fig.width, fig.height=fig.height)

  #
  ## Bona Fide
  #

  lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

  dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]

  dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

  vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

  p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
  save_plot(p_ligand_receptor_network_strict, paste(outfolder, paste("ligand_target_network_strict", conditionNames, sep="_"), sep="/"), fig.width=fig.width, fig.height=fig.height)

  #
  ## 6
  #


  DE_table_all = Idents(obj.in) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = obj.in, condition_colname = "condition", condition_oi = oiCond, condition_reference = refCond, expression_pct = 0.10, celltype_col = "idents") %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
  DE_table_all[is.na(DE_table_all)] = 0

  # Combine ligand activities with DE information
  ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
  ligand_activities_de[is.na(ligand_activities_de)] = 0

  # make LFC heatmap
  lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
  rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

  order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
  vis_ligand_lfc = lfc_matrix[order_ligands,]

  colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

  p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
  save_plot(p_ligand_lfc, paste(outfolder, paste("ligand_lfc", conditionNames, sep="_"), sep="/"), fig.width=fig.width, fig.height=fig.height)


  #
  ##
  ### Ligand signaling path
  ##
  #

  save_png <- function(plot, path){
    DiagrammeRsvg::export_svg(plot) %>%
      charToRaw() %>%
      rsvg::rsvg() %>%
      png::writePNG(path)
  }

  #ligands_all = "TGFB3" # this can be a list of multiple ligands if required
  #targets_all = c("TGFBI","LAMC2","TNC")

  #active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)

  # For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
  #active_signaling_network_min_max = active_signaling_network
  #active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
  #active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

  #graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")
  #save_png(graph_min_max, paste(outfolder, paste("graph", "sepsis_ctrl", sep="_"), sep="/"))

  return(list( "lr_network"=vis_ligand_receptor_network, "lr_network_strict"=vis_ligand_receptor_network_strict, "ligand_target"=vis_ligand_target, "active_ligand_target_links"=active_ligand_target_links, "active_ligand_target_links_df"=active_ligand_target_links_df))
}





scatterAnalysisPlot = function(mat, gene1, gene2, logg1=TRUE, logg2=TRUE, formula=y ~ x)
{

  smat = mat[mat$gene %in% c(gene1, gene2),]

  rownames(smat) = smat$gene
  smat$gene = NULL

  if (logg1)
  {
    smat[gene1, ] = log1p(smat[gene1, ])
  }

  if (logg2)
  {
    smat[gene2, ] = log1p(smat[gene2, ])
  }
  

  df = as.data.frame(t(smat))

  lowerlimX = plyr::round_any(min(df[[gene1]]), 0.1, f = floor)-1
  upperlimX = plyr::round_any(max(df[[gene1]]), 0.1, f = ceiling)+1

  lowerlimY = plyr::round_any(min(df[[gene2]]), 0.1, f = floor)-1
  upperlimY = plyr::round_any(max(df[[gene2]]), 0.1, f = ceiling)+1


  p = ggplot(df, aes_string(x=gene1, y=gene2)) +
  geom_point() + 
  geom_text(label=rownames(df))+stat_smooth(formula=formula, method = "lm")+
  stat_cor(r.accuracy = 0.01, label.y = 0.9*upperlimY, size = 4, label.sep='\n') +
  stat_regline_equation(aes(label = ..eq.label..), size = 4, label.y = 0.8*upperlimY)

  p = p + xlim(lowerlimX, upperlimX)
  p = p + ylim(lowerlimY, upperlimY)

  return(p)
}







