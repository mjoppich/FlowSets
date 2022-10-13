
import os, sys

import numpy as np
import pandas as pd
import skfuzzy as fuzz
from skfuzzy import control as ctrl
from statsmodels.stats.multitest import multipletests

from natsort import natsorted
from scipy.stats import hypergeom
from scipy.stats import chi2_contingency

from collections import defaultdict, OrderedDict

from scipy.special import expit
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable, hsv


import progressbar
def makeProgressBar() -> progressbar.ProgressBar:
    return progressbar.ProgressBar(widgets=[
        progressbar.Bar(), ' ', progressbar.Percentage(), ' ', progressbar.AdaptiveETA()
        ])



class SankeyPlotter:

    @classmethod
    def generate_colormap(cls, N):
        arr = np.arange(N)/N
        N_up = int(np.ceil(N/7)*7)
        arr.resize(N_up)
        arr = arr.reshape(7,N_up//7).T.reshape(-1)
        ret = hsv(arr)
        n = ret[:,3].size
        a = n//2
        b = n-a
        for i in range(3):
            ret[0:n//2,i] *= np.arange(0.2,1,0.8/a)
        ret[n//2:,3] *= np.arange(1,0.1,-0.9/b)
    
        return ret
        
    @classmethod
    def sigmoid_curve(cls, p1, p2, resolution=0.1, smooth=0):
        x1, y1 = p1
        x2, y2 = p2
        
        xbound = 6 + smooth

        fxs = np.arange(-xbound,xbound+resolution, resolution)
        fys = expit(fxs)
        
        x_range = x2 - x1
        y_range = y2 - y1
        
        xs = x1 + x_range * ((fxs / (2*xbound)) + 0.5)
        ys = y1 + y_range * fys
        
        return xs, ys


    @classmethod 
    def sigmoid_arc(cls, p1, w1, p2, w2=None, resolution=0.1, smooth=0, ax=None):
        
        xs, ys1 = cls.sigmoid_curve(p1, p2, resolution, smooth)
        
        if w2 is None:
            w2 = w1
        
        p1b = p1[0], p1[1] - w1
        p2b = p2[0], p2[1] - w2

        xs, ys2 = cls.sigmoid_curve(p1b, p2b, resolution, smooth)
        
        return xs, ys1, ys2

    @classmethod
    def sankey(cls, flow_matrix=None, node_positions=None, link_alpha=0.5, colours=None, 
            colour_selection="source", resolution=0.1, smooth=0, **kwargs):
        #node_widths = [np.max([i, o]) for i, o in zip(in_totals, out_totals)]
        n = np.max(flow_matrix.shape)
        in_offsets = [0] * n
        out_offsets = [0] * n

        ax = kwargs.get("ax", plt.gca())
        
        for i, b1 in enumerate(node_positions):
            outputs = flow_matrix[i,:]
            for j, (w, b2) in enumerate(zip(outputs, node_positions)):
                if w:
                    p1 = b1[0], b1[1] - out_offsets[i]
                    p2 = b2[0], b2[1] - in_offsets[j]
                    xs, ys1, ys2 = cls.sigmoid_arc(p1, w, p2, resolution=resolution, smooth=smooth, ax=ax)
                    out_offsets[i] += w
                    in_offsets[j] += w
                
                    c = 'grey'

                    if type(colours) == str:
                        c = colours
                    elif type(colours) == list:
                        if colour_selection == "sink":
                            c = colours[j]
                        elif colour_selection == "source":
                            c = colours[i]
                    plt.fill_between(x=xs, y1=ys1, y2=ys2, alpha=link_alpha, color=c, axes=ax)

    @classmethod
    def _make_plot( cls, nodeWeigthSequence, series2name, levelOrder, seriesOrder, specialColors=None, transformCounts = lambda x: x):

        nodePositions = {}

        for nodeName in series2name:
            for nli, nodeLevel in enumerate(levelOrder):

                nodePositions[ (nodeName, nodeLevel) ] = (seriesOrder.index(nodeName), 2*nli)

        minNodeLevel = min([nodePositions[x][1] for x in nodePositions])
        maxNodeLevel = max([nodePositions[x][1] for x in nodePositions])

        #nodeWeigthSequence = [ (("WT", 2), ("KO", 0), 1), (("WT", 2), ("KO", -2), 1), (("WT", -1), ("KO", -2), 1) ]
        nodeOffsets = defaultdict(lambda: 0)

        colours = cls.generate_colormap(len(nodeWeigthSequence))

        maxFlowPerNode = defaultdict(lambda: 0)

        for si, fIDWeights in enumerate(nodeWeigthSequence):
            fid, nws = fIDWeights
            weight = transformCounts(nws[-1])
            nodes = nws[0:-1]

            for node in nodes:
                maxFlowPerNode[node] += weight

        maxFlowInAllNodes = max([x for x in maxFlowPerNode.values()])

        fig, ax = plt.subplots(figsize=(20,10))
        ax.axis('off')
        plt.title("")

        for si, fIDWeights in enumerate(nodeWeigthSequence):

            fid, nws = fIDWeights

            weight = transformCounts(nws[-1]) / maxFlowInAllNodes
            nodes = nws[0:-1]

            for i in range(1, len(nodes)):

                src = nodes[i-1]
                tgt = nodes[i]

                p1 = nodePositions[src]
                p2 = nodePositions[tgt]

                p1 = p1[0], p1[1] - nodeOffsets[src] + maxFlowPerNode[src]/maxFlowInAllNodes/2.0
                p2 = p2[0], p2[1] - nodeOffsets[tgt] + maxFlowPerNode[tgt]/maxFlowInAllNodes/2.0

                if tgt == ("KO", 0):
                    print(p2)

                xs, ys1, ys2 = cls.sigmoid_arc(p1, weight, p2, resolution=0.1, smooth=0, ax=ax)

                nodeOffsets[src] += weight

                if tgt[0] == seriesOrder[-1]:
                    nodeOffsets[tgt] += weight

                if not specialColors is None:
                    if fid in specialColors:
                        c = specialColors[fid]
                    else:
                        c = "grey"                    
                else:
                    c = colours[si % len(colours)]
                plt.fill_between(x=xs, y1=ys1, y2=ys2, alpha=0.5, color=c, axes=ax)



        props = dict(boxstyle='round', facecolor='grey', alpha=1.0, pad=1)

        for npi, nn in enumerate(nodePositions):

            nodeStr = "{lvl}".format(cond=series2name[nn[0]], lvl=nn[1])
            nodePosition = nodePositions[nn]

            ax.text(nodePosition[0], nodePosition[1], nodeStr, transform=ax.transData, fontsize=14,rotation=90,
                verticalalignment='center', ha='center', va='center', bbox=props)

        # place a text box in upper left in axes coords

        for si, series in enumerate(series2name):

            ax.text(si, minNodeLevel-1, series2name[series], transform=ax.transData, fontsize=14,rotation=0,
                verticalalignment='center', ha='center', va='center', bbox=props)


        plt.ylim((minNodeLevel-1.5, maxNodeLevel+0.5))

        plt.tight_layout()
        plt.show()

def to_fuzzy(value, fzy):
    return np.array(
        [fuzz.interp_membership(fzy.universe, fzy[x].mf, value) for x in fzy.terms]
    )


def distribution_to_fuzzy(meanValue, sdValue, exprCells, fzMFs, threshold=0.0):

    fuzzySet = [0] * len(fzMFs.terms)

    if sdValue is None:
        return np.array([x for x in zip(fzMFs.terms, to_fuzzy(meanValue, fzMFs))])

    if np.isnan(meanValue) or np.isnan(sdValue):
        fuzzySet[0] = 1
        return [x for x in zip(fzMFs.terms, fuzzySet)]

    normValues = np.random.normal(meanValue, sdValue, 100)
    
    for v in normValues:
        fuzzySet += to_fuzzy(v, fzMFs)

    fuzzySet = fuzzySet / len(normValues)   

    fuzzySetNoExpr = [0] * len(fzMFs.terms)
    fuzzySetNoExpr[1] = 1

    fuzzySet = (1-exprCells) * np.array(fuzzySetNoExpr) + exprCells * fuzzySet
    fuzzySet[fuzzySet < threshold] = 0
    fuzzySet = fuzzySet / np.sum(fuzzySet)
    outset = np.array([x for x in zip(fzMFs.terms, fuzzySet)])

    return outset

def toWideDF( indf ):
    dfPivot = pd.pivot(indf, index=["gene"], columns=["cluster"], values=["fuzzy_set"])
    dfWide = dfPivot.copy()
    dfWide.columns = dfWide.columns.droplevel(0)
    dfWide.reset_index(inplace=True)
    dfWide.reset_index(drop=True, inplace=True)
    return dfWide

def to_homogeneous(df, exprMFs, is_foldchange=False):

    if not is_foldchange:
        emptyElem = [x for x in zip(exprMFs.terms, to_fuzzy(0, exprMFs))]
    else:
        emptyElem = [x for x in zip(exprMFs.terms, to_fuzzy(0, exprMFs))]

    dfc = df.copy()
    dfc = dfc.applymap(lambda x: x if not x is np.nan else emptyElem)

    return dfc

def identify_threshold_level(dfWideNew, force_calculation = False):

    clusterColumns = [x for x in dfWideNew.columns if x.startswith("cluster.")]
    explDF = dfWideNew.copy()
    explDF.reset_index(drop=True, inplace=True)

    for thresholdIteration in range(0,11, 1):
        currentThreshold = thresholdIteration * 0.05

        explDF = dfWideNew.copy()
        explDF.reset_index(drop=True, inplace=True)

        print("Current Threshold", currentThreshold)
        finished = True

        for cidx, ccol in enumerate(clusterColumns):
            print(ccol)
            print(explDF.shape)

            initialCount = explDF.shape[0]

            explDF = explDF.explode(ccol)
            explDF.reset_index(drop=True, inplace=True)
            clDF = pd.DataFrame(explDF[ccol].to_list())

            binColName = "bin.{}".format(ccol)
            mfColName = "mf.{}".format(ccol)

            clDF.columns = [binColName, mfColName]

            clDF[binColName] = clDF[binColName].astype("category")
            clDF[mfColName] = clDF[mfColName].astype('float16')

            explDF = pd.concat([explDF, clDF], axis=1)
            explDF = explDF[ explDF[mfColName] > currentThreshold]

            del explDF[ccol]

            afterCount = explDF.shape[0]

            increase = afterCount/initialCount
            remaining = len(clusterColumns)-cidx-1
            prospIncrease = increase ** remaining

            print("increase", increase)
            print("remaining", remaining, prospIncrease, afterCount * prospIncrease)

            if not force_calculation and afterCount * prospIncrease > 20000000:
                print("Should abort this run", afterCount * prospIncrease)
                finished = False
                break
            
        if finished:    
            print("Finished")
            break

    return explDF, currentThreshold



class CustomFuzzyVar(ctrl.Antecedent):

    def __init__(self, universe, label):
        super().__init__(universe, label)


    def automf(self, number=5, names=None, centers=None, shape="tri"):

        assert(shape in ("tri", "sig"))
        assert(len(names) == number)

        if not centers is None:
            assert(len(centers) == len(names))

        limits = [self.universe.min(), self.universe.max()]
        universe_range = limits[1] - limits[0]
        widths = [universe_range / ((number - 1) / 2.)] * int(number)

        if centers is None:
            centers = np.linspace(limits[0], limits[1], number)
            widths = [universe_range / ((number - 1) / 2.)] * int(number)
        else:
            widths = []
            for i in range(0, len(centers)-1):
                leftWidth = None
                rightWidth = None

                if i > 0:
                    leftWidth = centers[i]-centers[i-1]
                
                rightWidth = centers[i+1]-centers[i]

                if not leftWidth is None:
                    avgWidth = (rightWidth+leftWidth)/2
                else:
                    avgWidth = rightWidth

                widths.append(2*avgWidth)
            

        abcs = [[c - w / 2, c, c + w / 2] for c, w in zip(centers, widths)]
        cws = [[c, w] for c, w in zip(centers, widths)]

        # Clear existing adjectives, if any
        self.terms = OrderedDict()

        # Repopulate
        for name, abc in zip(names, abcs):

            if shape == "tri":
                self[name] = fuzz.trimf(self.universe, abc)
            elif shape == "sig":
                self[name] = fuzz.sigmf(self.universe, cws[0], cws[1])


def to_arg_max(df, exprMFs):

    emptyElem = [0] * len(exprMFs.terms)

    dfc = df.copy()
    dfc = dfc.applymap(lambda x: x if not x is np.nan else emptyElem)

    dfIdx = dfc.applymap( lambda x: np.argmax(x))
    dfValue = dfc.applymap( lambda x: np.max(x))
    dfValue = dfValue.fillna(0)

    return dfIdx, dfValue













class FlowAnalysis:

    @classmethod
    def make_fuzzy_concepts(cls, exprData, mfLevels, centers, meancolName, mfLevelsMirrored):


        if "max.cluster" in exprData.columns:
            minValue = np.floor(exprData["max.cluster"].min())
            maxValue = np.ceil(exprData["max.cluster"].max())
        else:
            minValue = np.floor(exprData[meancolName].min())
            maxValue = np.ceil(exprData[meancolName].max())

        if mfLevelsMirrored:
            absValue = max(abs(minValue), abs(maxValue))
            minValue = -absValue
            maxValue = absValue

            assert(len(mfLevels) % 2 == 1)

        print("Creating Range", minValue, "->", maxValue)

        stepSize = min(0.1, (maxValue-minValue)/100)

        exprMFs = CustomFuzzyVar(np.arange(minValue, maxValue, stepSize), 'exprMFs')
        exprMFs.automf(len(mfLevels), names=mfLevels, centers=centers) 
        # You can see how these look with .view()
        exprMFs.view()

        return exprMFs

    @classmethod
    def exprDF2LongDF(cls, indf, mfLevels = ["NO", "LOW", "med", "HIGH"], mfLevelsMirrored=False, centers=None, meancolName="mean.cluster", sdcolName="sd.cluster", exprcolName="expr.cluster"):

        exprData = indf.copy()

        exprMFs = cls.make_fuzzy_concepts(exprData, mfLevels, centers, meancolName, mfLevelsMirrored)

        meanExprCol = list(exprData.columns).index(meancolName)
        exprCountCol = list(exprData.columns).index(exprcolName)

        if not sdcolName is None:
            sdExprCol = list(exprData.columns).index(sdcolName)

        df = exprData.copy()

        if not sdcolName is None:
            df["fuzzy_set"] = df.apply(lambda row : distribution_to_fuzzy(row[meanExprCol], row[sdExprCol], row[exprCountCol], exprMFs, threshold=0.0), axis=1)
        else:
            print("No SD col name given")
            df["fuzzy_set"] = df.apply(lambda row : distribution_to_fuzzy(row[meanExprCol], None, row[exprCountCol], exprMFs, threshold=0.0), axis=1)

        dfWide = to_homogeneous(toWideDF(df), exprMFs)
        _, threshold = identify_threshold_level(dfWide)

        print("Identified Threshold", threshold)

        dfNew = exprData.copy()
        if not sdcolName is None:
            dfNew["fuzzy_set"] = dfNew.apply(lambda row : distribution_to_fuzzy(row[meanExprCol], row[sdExprCol], row[exprCountCol], exprMFs, threshold=0.0), axis=1)
        else:
            dfNew["fuzzy_set"] = dfNew.apply(lambda row : distribution_to_fuzzy(row[meanExprCol], None, row[exprCountCol], exprMFs, threshold=0.0), axis=1)
        
        dfNewWide = to_homogeneous(toWideDF(dfNew), exprMFs)

        explDF, _ = identify_threshold_level(dfNewWide, force_calculation=True)

        return explDF, exprMFs

    @classmethod
    def toFlowsDF(cls, indf):

        explDF = indf.copy()
        binColumns = [x for x in explDF.columns if x.startswith("bin.")]
        print(binColumns)

        mfColumns = [x.replace("bin.", "mf.") for x in binColumns]
        print(mfColumns)

        explDF['group.flow'] = explDF.apply( lambda row: tuple( [row[x] for x in binColumns] ), axis=1)
        explDF['mf.flow'] = explDF.apply( lambda row: np.prod( [row[x] for x in mfColumns] ), axis=1)
        allgroups = list(set(explDF["group.flow"]))
        explDF['id.flow'] = explDF.apply( lambda row: allgroups.index(row["group.flow"]), axis=1)

        return explDF


    def __init__(self, flows, symbol_column, series2name, exprMF):
        
        self.flows = flows
        self.seriesOrder = [x[0] for x in series2name]
        self.series2name = {x[0]: x[1] for x in series2name}
        self.symbol_column = symbol_column
        self.exprMF = exprMF
        self.levelOrder =  [x for x in self.exprMF.terms]

        self.flows[self.symbol_column] = self.flows[self.symbol_column].str.upper()

        self.flowgroup_flow, self.flowgroup_route, self.flowgroup_genes = self.prepare_flows()

    def prepare_flows(self, flowDF=None):

        if flowDF is None:
            flowDF = self.flows

        flowgroup_flow = defaultdict(lambda: 0)
        flowgroup_route = {}
        flowgroup_genes = defaultdict(set)

        for ri, row in flowDF.iterrows():
            fgid = row["id.flow"]
            fgflow = row["mf.flow"]
            fgroute = row["group.flow"]
            fggene = row[self.symbol_column]

            flowgroup_route[fgid] = [(x,y) for x,y in zip(self.seriesOrder, fgroute)]
            flowgroup_flow[fgid] += fgflow
            flowgroup_genes[fgid].add(fggene)

        return flowgroup_flow, flowgroup_route, flowgroup_genes


    def plot_flows(self, use_flows = None):

        if use_flows is None:
            use_flows = [x for x in self.flowgroup_route]

        weightSequence = self._to_weight_sequence( self.flowgroup_route, self.flowgroup_flow)

        SankeyPlotter._make_plot(weightSequence, self.series2name, self.levelOrder, self.seriesOrder, transformCounts=lambda x: np.sqrt(x))

    def plot_genes(self, genes):

        if not isinstance(genes, (tuple, list)):
            genes = [genes]

        subsetFlows = self.flows.copy()
        subsetFlows = subsetFlows[subsetFlows[self.symbol_column].isin(genes)]

        print(subsetFlows.shape)

        flowgroup_flow, flowgroup_route, flowgroup_genes = self.prepare_flows(subsetFlows)

        weightSequence = self._to_weight_sequence( flowgroup_route, flowgroup_flow)

        SankeyPlotter._make_plot(weightSequence, self.series2name, self.levelOrder, self.seriesOrder, transformCounts=lambda x: np.sqrt(x))

    def _to_weight_sequence(self, flowgroup_route, flowgroup_flow):
        weightSequence = []
        for fgid in [x for x in flowgroup_route]:
            weightSequence.append( (fgid, tuple(flowgroup_route[fgid] + [flowgroup_flow.get(fgid, 0)])) )

        return weightSequence

    def highlight_genes(self, genes):

        if not isinstance(genes, (tuple, list)):
            genes = [genes]

        subsetFlows = self.flows.copy()
        subsetFlows.loc[subsetFlows[self.symbol_column].isin(genes),"id.flow"] = -1
        
        subsetFlowsBG = subsetFlows[subsetFlows["id.flow"] != -1].copy()
        subsetFlowsFG = subsetFlows[subsetFlows["id.flow"] == -1].copy()

        flowgroup_flow, flowgroup_route, flowgroup_genes = self.prepare_flows(subsetFlowsBG)
        weightSequenceBG = self._to_weight_sequence( flowgroup_route, flowgroup_flow)

        allgroups = list(set(subsetFlowsFG["group.flow"]))
        subsetFlowsFG['id.flow'] = subsetFlowsFG.apply( lambda row: (-1)*(allgroups.index(row["group.flow"])+1), axis=1)

        flowgroup_flow, flowgroup_route, flowgroup_genes = self.prepare_flows(subsetFlowsFG)
        weightSequenceFG = self._to_weight_sequence( flowgroup_route, flowgroup_flow)

        specialColors = {x: "red" for x in set(subsetFlowsFG['id.flow'])}

        SankeyPlotter._make_plot(weightSequenceBG+weightSequenceFG, self.series2name, self.levelOrder, self.seriesOrder, specialColors=specialColors, transformCounts=lambda x: np.sqrt(x))

    def flow_finder( self, sequence, verbose=True ):

        seriesOrder = list(self.exprMF.terms.keys())

        matchingFIDs = set()

        fidGroups = self.flows[["id.flow", "group.flow"]]
        fidGroups = fidGroups.drop_duplicates()

        for fid, fgroup in zip(fidGroups["id.flow"], fidGroups["group.flow"]):

            acceptFlow = True
            
            for ci, comp in zip(range(0, len(sequence)), sequence):

                startIdx = seriesOrder.index(fgroup[ci])
                endIdx = seriesOrder.index(fgroup[ci+1])

                if comp == "<":
                    if not startIdx < endIdx:
                        acceptFlow=False
                        break
                elif comp == "<=":
                    if not startIdx <= endIdx:
                        acceptFlow=False
                        break
                elif comp == ">=":
                    if not startIdx >= endIdx:
                        acceptFlow=False
                        break
                elif comp == ">":
                    if not startIdx > endIdx:
                        acceptFlow=False
                        break
                elif comp == "=":
                    if not startIdx == endIdx:
                        acceptFlow=False
                        break

                if not acceptFlow:
                    break

            if acceptFlow:
                if verbose:
                    print(fgroup)
                matchingFIDs.add(fid)

        return matchingFIDs

    def read_gmt_file(self, filepath):

        geneset2genes = {}

        with open(filepath) as fin:

            for line in fin:

                line = line.strip().split("\t")
                pwName = line[0]
                pwID = line[1]
                pwGenes = line[2:]
                pwGenes = [x.upper() for x in pwGenes]

                geneset2genes[pwID] = (pwName, pwGenes)

        return geneset2genes

    def analyse_pathways(self):


        rp = self.read_gmt_file("ReactomePathways.gmt")
        allDFs = []

        bar = makeProgressBar()

        for fgid in bar(self.flowgroup_genes):
            #print(fgid)

            fg_df = self.flows[self.flows["id.flow"] == fgid]

            df = self.analyse_genes_for_genesets(rp, fg_df)
            df["fgid"] = fgid

            #print(df[df["pwsize"] > 1].sort_values(["pval"], ascending=True).head(3))

            allDFs.append(df)

        allFGDFs = pd.concat(allDFs, axis=0)

        _ , elemAdjPvals, _, _ = multipletests(allFGDFs["pval"], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        allFGDFs["adj_pval"] = elemAdjPvals

        return allFGDFs

    def analyse_genes_for_genesets(self, pathways, flowDF, populationSize=None):
        """
        
        pathways: pathway object
        genes: genes of interest

        populationSize: size of universe. if None, all genes in pathways will be chosen
        
        """

        flowDF._is_copy = None
        flowDF.loc[:, self.symbol_column] = list(flowDF.loc[:,self.symbol_column].str.upper())

        allPathwayGenes = set()
        for x in pathways:
            pwname, pwGenes = pathways[x]
            for gene in pwGenes:
                allPathwayGenes.add(gene)

        allPathwayGenes = allPathwayGenes.intersection(self.flows[self.symbol_column])

        pwGeneFlow = self.flows[self.flows[self.symbol_column].isin(allPathwayGenes)]
        totalFlowSumOfPathwayGenes = pwGeneFlow["mf.flow"].sum()

        if populationSize is None:
            populationSize = len(set(self.flows[self.symbol_column]))

        flowGenes = set(flowDF[self.symbol_column])
        flowScore = flowDF["mf.flow"].sum()
        numSuccInPopulation = len(allPathwayGenes)
        
        outData = defaultdict(list)

        for pwID in pathways:

            pwName, pwGenes = pathways[pwID]
            inflow_inset = flowGenes.intersection(pwGenes)

            pathwayDF = flowDF[flowDF[self.symbol_column].isin(pwGenes)]
            flowInPathwayScore = pathwayDF["mf.flow"].sum()


            if abs(flowInPathwayScore) < 1 or abs(flowScore) < 1:
                chi2 = 0
                pval = 1.0
            else:

                table=np.array([[flowScore,flowInPathwayScore], #len(inflow_inset)
                                [populationSize-flowScore,len(pwGenes)-flowInPathwayScore]
                            ])
                chi2, pval, dof, expected=chi2_contingency(table, correction=True)

            # population: all genes
            # condition: genes
            # subset: pathway

            genes_coverage = flowInPathwayScore / flowScore if flowScore > 0 else 0
            pathway_coverage = flowInPathwayScore / len(pwGenes) if len(pwGenes) > 0 else 0

            outData["pwid"].append(pwID)
            outData["pwname"].append(pwName)
            outData["flow_pw_score"].append(flowInPathwayScore)
            outData["pwsize"].append(len(pwGenes))
            outData["flow_score"].append(flowScore)
            outData["flow_size"].append(len(flowGenes))

            outData["pw_gene_intersection"].append(len(inflow_inset))
            outData["pw_coverage"].append(pathway_coverage)
            outData["genes_coverage"].append(genes_coverage)
            outData["pval"].append(pval)
            outData["chi2"].append(chi2)
            outData["mean_coverage"].append(pathway_coverage*genes_coverage)

        return pd.DataFrame.from_dict(outData)