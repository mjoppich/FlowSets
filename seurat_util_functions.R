
library(ggplot2)
library(data.table)
library(writexl)
library(purrr)
library(dplyr)
library(cowplot)
library(stringr)

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
  
  expt.r = as(expTable, "dgTMatrix")
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

cellIDForClusters = function(obj.in, targetVar, clusters)
{
  
  targetVarDF = as.data.frame(obj.in[[targetVar]])
  #print(paste("orig:", length(rownames(targetVarDF))))
  cellNames = rownames(targetVarDF)[targetVarDF[[targetVar]] %in% clusters]
  
  #print(length(cellNames))
  return(cellNames)
  
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



#
#
#
#
# FELIX: this is stupid simple fuzzification of logFC values. 
#
#
#
#
#


df_merge_outer <- function(df1, df2){                                # Create own merging function
  merge(df1, df2, by = "gene", all=TRUE)
}

createBins = function(indf, col.fc, col.sig, bounds, sigThreshold=0.05)
{
  
  smalldf = indf[, c("gene", col.fc, col.sig)]
  colnames(smalldf) = c("gene", "fc", "sig")
  
  binFun = function(x)
  {
    
    fc = as.numeric(x[2])
    sig = as.numeric(x[3])
    
    if ((is.na(sig)) || is.nan(sig) ||(sig > sigThreshold) || is.na(fc) || is.nan(fc))
    {
      return(0)
    }
    
    curBin = 0
    
    if (fc > bounds[length(bounds)])
    {
      return(length(bounds))
    }
    
    for (i in 1:length(bounds))
    {
      if (fc < bounds[i])
      {
        return(curBin + i -1)
      }
    }
    
    return(10)
  }
  
  bins = apply(smalldf, 1, binFun)
  return(bins)
}





prepareFuzzyData = function(deObj, celltype, type="expr", use.quantiles=c(0.25, 0.75), column.name = "gene", column.pval="", column.data="avg_log2FC")
{
  #deObj is a list with "series data" as name
  
  if (! type %in% c("expr", "logFC"))
  {
    print("invalid type")
    return(NULL)
  }
  
  deDFs = deObj[[celltype]]  
  
  celltype.logfcs = c()
  short.results = list()
  
  
  for (condtp in names(deDFs))
  {
    colname = column.data
    celltype.logfcs = c(celltype.logfcs, as.vector(deDFs[[condtp]][[colname]]))
  }
  
  if (type == "expr")
  {
    # only expressed data
    celltype.logfcs = celltype.logfcs[celltype.logfcs>0]
  }
  
  quantileVec = as.numeric(quantile(abs(celltype.logfcs), probs = use.quantiles, na.rm=T))
  
  if (type == "logFC")
  {
    # up and down regulated!
    quantileVec = c(-rev(quantileVec), quantileVec)
  }
  
  print(paste("Quantile Vector"))
  print(quantileVec)
  
  condtp.results = list()
  for (tp in names(deDFs))
  {
    
    celltype.logfcs = c(celltype.logfcs, deDFs[[condtp]][[tp]][[column.data]])
    
    resDF = data.frame(gene=deDFs[[tp]][[column.name]])
    resDF[[paste("data", tp, sep=".")]] = as.vector(deDFs[[tp]][[column.data]])
    
    resDF[[paste("p_val_adj", tp, sep=".")]] = 0
    resDF[[paste("bins", tp, sep=".")]] = createBins(resDF, paste("data", tp, sep="."), paste("p_val_adj", tp, sep="."), quantileVec)
    
    condtp.results[[tp]] = resDF
  }
  
  celltype.results = Reduce(my_merge_outer, condtp.results) 
  
  
  
  for (binID in grep("^bins.", colnames(celltype.results), value = T))
  {
    celltype.results[[binID]][is.na(celltype.results[[binID]])] = 0
    
  }
  
  return(celltype.results)
  
}



compareClusters = function(scdata, cellsID1, cellsID2, suffix1, suffix2, prefix="cluster", test="MAST", assay="RNA", outfolder="./", fcCutoff=0.25, all=FALSE)
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
  
  outvalues1 = getExprData_getExtExpr(scdata, cellsID1, suffix1, assay=assay)
  outvalues2 = getExprData_getExtExpr(scdata, cellsID2, suffix2, assay=assay) 
  
  
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
  
  return(joinedData)
}

