
import os, sys

import numpy as np
import pandas as pd
import skfuzzy as fuzz
from skfuzzy.control.fuzzyvariable import FuzzyVariable
from skfuzzy.control.visualization import FuzzyVariableVisualizer
from statsmodels.stats.multitest import multipletests

from natsort import natsorted
from scipy.stats import hypergeom, norm
from scipy.stats import chi2_contingency, fisher_exact

from collections import defaultdict, OrderedDict
import itertools

from scipy.special import expit
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable, hsv, get_cmap
import seaborn as sns
import matplotlib.patches as patches
import matplotlib.colors as mcolors

import colorsys

import polars as pl

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
        arr.resize(N_up,refcheck=False)
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
    def _make_plot( cls, nodeWeigthSequence, series2name, levelOrder, seriesOrder, specialColors=None, nodeColors=None, fsize=None, transformCounts = lambda x: x, cmap=None, norm=None, outfile=None, title=None, verbose=False, linewidth=0.01, seriesColorMap=None, independentEdges=False):

        
        levelHeight=2

        maxNumLevels = max([len(levelOrder[x]) for x in levelOrder])
        maxLevelPos = (maxNumLevels-1)*levelHeight

        nodePositions = {}
        for nodeName in series2name:
            
            numLevelsInState = len(levelOrder[nodeName])
            
            stateHeight = (numLevelsInState-1)*levelHeight
            levelOffset = 0.5*(maxLevelPos-stateHeight)                        
            
            for nli, nodeLevel in enumerate(levelOrder[nodeName]):
                nodePositions[ (nodeName, nodeLevel) ] = (seriesOrder.index(nodeName), levelOffset+ (levelHeight*nli))

        if verbose:
            for x in nodePositions:
                print(x, nodePositions[x])

        minXValue = min([nodePositions[x][0] for x in nodePositions])
        maxXValue = max([nodePositions[x][0] for x in nodePositions])
        minYValue = min([nodePositions[x][1] for x in nodePositions])
        maxYValue = max([nodePositions[x][1] for x in nodePositions])

        #nodeWeigthSequence = [ (("WT", 2), ("KO", 0), 1), (("WT", 2), ("KO", -2), 1), (("WT", -1), ("KO", -2), 1) ]
        nodeOffsets = defaultdict(lambda: 0)
        maxFlowPerNode = defaultdict(lambda: 0)

        for si, fIDWeights in enumerate(nodeWeigthSequence):
            fid, nws = fIDWeights
            weight = transformCounts(nws[-1])
            nodes = nws[0:-1]

            for ni, node in enumerate(nodes):
                
                if independentEdges:
                    
                    if ni % 2 == 0:
                        maxFlowPerNode["{}_in".format(node)] += weight
                    else:
                        maxFlowPerNode["{}_out".format(node)] += weight

                else:
                    maxFlowPerNode[node] += weight


        numNodes = len(series2name)
        if independentEdges: 
            numNodes=(numNodes*2)-2

        maxFlowInAllNodes = sum([x for x in maxFlowPerNode.values()])/numNodes
        if verbose:
            print("Max Flow Per Node")
            print( sum([x for x in maxFlowPerNode.values()]))
            print(maxFlowInAllNodes)

            print(numNodes)
            for node in maxFlowPerNode:
                print(node, maxFlowPerNode[node])


        #reorder nodeWeigthSequence such that paths are ordered
        nodeWeigthSequence = sorted(nodeWeigthSequence, key=lambda x: [nodePositions[ne][1] for ne in x[1][0:-1]], reverse=True)

        nodeColors={}
        if not seriesColorMap is None:
            
            for npi, nn in enumerate(nodePositions):
                nodePosition = nodePositions[nn]
                nodeColor = seriesColorMap[ series2name[nn[0]] ]( nodePosition[1]/(levelHeight*(maxNumLevels-1)) )[:3]
                nodeColors[nn] = nodeColor

            colours = []
            for si, fIDWeights in enumerate(nodeWeigthSequence):

                fid, nws = fIDWeights
                startNode = nws[0]
                
                nodeColor = nodeColors[startNode]
                if cls.get_color_is_bright(nodeColor):
                    patchColor = cls.scale_lightness(nodeColor, 0.75)
                else:
                    patchColor = cls.scale_lightness(nodeColor, 1.25)
                
                colours.append( patchColor )
                
        else:
            colours = cls.generate_colormap(len(nodeWeigthSequence))


        # close any maybe still open plot
        plt.close()
        
        if fsize is None:
            fsize = (4 * len(seriesOrder), 2*(len(levelOrder)+1))

        if verbose:
            print("Figure Size", fsize)
            
        #create new figuer
        fig, ax = plt.subplots(figsize=fsize)
        ax.axis('off')
        plt.title("")

        for si, fIDWeights in enumerate(nodeWeigthSequence):

            fid, nws = fIDWeights
            
            weight = transformCounts(nws[-1]) / maxFlowInAllNodes
            nodes = nws[0:-1]
            
            if weight == 0:
                if verbose:
                    print("Skipping flow", fid, "for zero weight")
                continue

            for i in range(1, len(nodes)):

                src = nodes[i-1]
                tgt = nodes[i]

                p1 = nodePositions[src]
                p2 = nodePositions[tgt]
                
                if independentEdges:
                    srcName = "{}_in".format(src)
                    tgtName = "{}_out".format(tgt)
                else:
                    srcName = src
                    tgtName = tgt
                
                p1 = p1[0], p1[1] - nodeOffsets[srcName] + (maxFlowPerNode[srcName]/maxFlowInAllNodes)/2.0
                p2 = p2[0], p2[1] - nodeOffsets[tgtName] + (maxFlowPerNode[tgtName]/maxFlowInAllNodes)/2.0

                #left/right displacement so paths appear to hit outer box area
                p1 = (p1[0]+0.1, p1[1])
                p2 = (p2[0]-0.1, p2[1])

                xs, ys1, ys2 = cls.sigmoid_arc(p1, weight, p2, resolution=0.1, smooth=0, ax=ax)

                nodeOffsets[srcName] += weight

                if (tgt[0] == seriesOrder[-1]) or independentEdges:
                    nodeOffsets[tgtName] += weight

                if not specialColors is None:
                    if fid in specialColors:
                        c = specialColors[fid]
                    else:
                        c = "grey"                    
                else:
                    c = colours[si % len(colours)]
                                        
                    
                plt.fill_between(x=xs, y1=ys1, y2=ys2, alpha=0.7, color=c, axes=ax, linewidth=linewidth)


        for npi, nn in enumerate(nodePositions):

            nodeStr = "{lvl}".format(cond=series2name[nn[0]], lvl=nn[1])
            nodePosition = nodePositions[nn]
                                    
            nodeColor = nodeColors.get(nn, (0.7,0.7,0.7))

            rect = patches.FancyBboxPatch( (nodePosition[0]-0.1, nodePosition[1]-0.75), width=0.2, height=1.5, facecolor=nodeColor,linewidth=0)
            rect.set_boxstyle("round", rounding_size=0.1, pad=0)
            ax.add_patch(rect)
            
            
            if cls.get_color_is_bright(nodeColor):
                patchColor = cls.scale_lightness(nodeColor, 0.75)
            else:
                patchColor = cls.scale_lightness(nodeColor, 1.25)
            
            rect = patches.Rectangle( (nodePosition[0]-0.1, nodePosition[1]-0.5), width=0.2, height=1, linewidth=0, facecolor=patchColor )
            ax.add_patch(rect)

            textColor = cls.get_text_color_based_on_background_color(nodeColor, "#FFFFFF", "#000000")
            

            t = ax.text(nodePosition[0], nodePosition[1], nodeStr, transform=ax.transData, fontsize=14,rotation=90,
                verticalalignment='center', ha='center', va='center', bbox=None, color=textColor)
            

        # place a text box in upper left in axes coords

        for si, series in enumerate(series2name):
            props = dict(boxstyle='round', facecolor="lightgrey", alpha=1.0, pad=0.5)

            bwidth = 0.25
            bheight= 0.4
            
            if not seriesColorMap is None:
                faceColor = seriesColorMap[ series2name[series] ]( 0.9 )[:3]
            else:
                faceColor = mcolors.to_rgb("lightgrey")
            
            faceColorRGB = mcolors.to_rgb(faceColor)
            textColor = cls.get_text_color_based_on_background_color(faceColorRGB, "#FFFFFF", "#000000")
            
            rect = patches.FancyBboxPatch( (si-0.5*bwidth, -2 - 0.5*bheight), width=bwidth, height=bheight, facecolor=faceColorRGB, linewidth=0)
            rect.set_boxstyle("round", rounding_size=0.05, pad=0)
            ax.add_patch(rect)

            ax.text(si, -2, series2name[series], transform=ax.transData, fontsize=14,rotation=0,
                verticalalignment='center', ha='center', va='center', bbox=None, color=textColor)


        #plt.ylim((minNodeLevel-1.5, maxNodeLevel+0.5))
        #plt.subplots_adjust(left=0.4, right=0.6, bottom=0.4, top=0.6)
        #plt.tight_layout()
        
        if not cmap is None:
            #cb_ax = fig.add_axes([0.27, 0.8, 0.5, 0.05])            
            cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='horizontal', shrink=0.25, fraction=0.005)
            cb.set_label("Membership")
        
        xpadding=0.25
        ypadding=0.5
        plt.xlim(minXValue-xpadding, maxXValue+xpadding)
        plt.ylim(-2-ypadding, maxYValue+1+ypadding)
                
        if not title is None:
            print("Adding title", title)
            plt.title(title)
        
        if not outfile is None:
            plt.savefig(outfile + ".png", bbox_inches='tight')
            plt.savefig(outfile + ".pdf", bbox_inches='tight')

        plt.show()
        plt.close()


    @classmethod
    def get_color_is_bright(cls, bgColor):
        r, g, b = bgColor
        uicolors = [1.0-r, 1.0-g, 1.0-b]
        adjusted = []
        for col in uicolors:
            col2 = col
            if col <= 0.03928:
                col2 = col/12.92
            col2 = pow((col2 + 0.055)/1.055,2.4)
            adjusted.append(col2)
        L = (0.2126 * adjusted[0] + 0.7152 * adjusted[1] + (0.072 * adjusted[2]))
        
        return L < 0.179

    @classmethod
    def get_text_color_based_on_background_color(cls, bgColor, lightColor, darkColor):

        if cls.get_color_is_bright(bgColor):
            return darkColor
        else:
            return lightColor
        
    @classmethod
    def scale_lightness(cls, rgb, scale_l):
        # convert rgb to hls
        h, l, s = colorsys.rgb_to_hls(*rgb)
        # manipulate h, l, s values and return as rgb
        return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

    @classmethod
    def createColorMap(cls, inColor, mode="scaling"):
        
        inColorRGB = mcolors.to_rgb(inColor)
        
        if mode == "scaling":
            

            
            brightInColor = cls.scale_lightness(inColorRGB, 2)
            colorList = [brightInColor, inColorRGB]
        
        elif mode == "diverging":
            
            def get_complementary(color):
                color = "#%02x%02x%02x" % color
                color = color[1:]
                color = int(color, 16)
                comp_color = 0xFFFFFF ^ color
                comp_color = "#%06X" % comp_color
                return mcolors.to_rgb(comp_color)
            
            inCColorRGB = get_complementary(tuple([int(x*255) for x in inColorRGB]))
            colorList = [cls.scale_lightness(inCColorRGB, 0.5), inCColorRGB,  cls.scale_lightness(inColorRGB, 1.5), inColorRGB]
            
        cmap = mcolors.LinearSegmentedColormap.from_list("custom", colorList)
        
        return cmap
        

def to_fuzzy(value, fzy):

    return np.array(
        [fuzz.interp_membership(fzy.universe, fzy[x].mf, value) for x in fzy.terms], dtype=np.float32
    )

def to_crisp(value, fzy):
    value=round(value,1)
    return np.array(
        [round(fuzz.interp_membership(fzy.universe, fzy[x].mf, value),0) for x in fzy.terms]
    )

def distribution_to_crisp(meanValue,fzMFs, threshold=0.0):
    return [x for x in to_crisp(meanValue, fzMFs)]



def distribution_to_fuzzy(meanValue, sdValue, exprCells, fzMFs, threshold=0.0):

    #fuzzySetNoExpr = [0.0] * len(fzMFs.terms)
    #fuzzySetNoExpr[0] = 1.0
    fuzzySetNoExpr = to_fuzzy(0, fzMFs)
    fuzzySet = [0.0] * len(fzMFs.terms)

    if sdValue is None:
        fuzzySet =  to_fuzzy(meanValue, fzMFs)
        #fuzzySet = (1-exprCells) * np.array(fuzzySetNoExpr) + exprCells * fuzzySet
        #fuzzySet[fuzzySet < threshold] = 0
        #return [list(x) for x in zip(fzMFs.terms, to_fuzzy(meanValue, fzMFs))]


    if np.isnan(meanValue) or ((not sdValue is None) and np.isnan(sdValue)):
        fuzzySet = fuzzySetNoExpr
        #return [list(x) for x in zip(fzMFs.terms, fuzzySet)]
        return [x for x in fuzzySet]


    if not sdValue is None:
        normValues = np.random.normal(meanValue, sdValue, 100)
        for v in normValues:
            fuzzySet += to_fuzzy(v, fzMFs)

        fuzzySet = fuzzySet / len(normValues)   

    fuzzySet = (1-exprCells) * np.array(fuzzySetNoExpr) + exprCells * fuzzySet
    fuzzySet[fuzzySet < threshold] = 0
    
    fuzzySet = np.array(fuzzySet, dtype=np.float32) # this should always match identify_threshold_level 
    fuzzySet = fuzzySet / np.sum(fuzzySet)
    
    outset = [x for x in fuzzySet]
    #return [list(x) for x in zip(fzMFs.terms, fuzzySet)]
    return outset

def toWideDF( df,symbol_column ):
    #dfPivot = pd.pivot(indf, index=["gene"], columns=["cluster"], values=["fuzzy_set"])
    #dfWide = dfPivot.copy()
    #dfWide.columns = dfWide.columns.droplevel(0)
    #dfWide.reset_index(inplace=True)
    #dfWide.reset_index(drop=True, inplace=True)
    
    dfPivot = df.pivot(values=["fuzzy.mfs"], index=symbol_column, columns="cluster")
    return dfPivot

def to_homogeneous(df, exprMFs, is_foldchange=False):
        
    #print(fuzzyBins)
    #print(fuzzyValues)

    for col in exprMFs:
        
        print(col)
                
        exprMF = exprMFs[col]
        
        if not is_foldchange:
            #fuzzyBins = [x for x in exprMF.terms]
            fuzzyValues = [float(x) for x in to_fuzzy(0, exprMF)]
        else:
            #fuzzyBins = [x for x in exprMF.terms]
            fuzzyValues = [float(x) for x in to_fuzzy(0, exprMF)]
                    
        df=df.with_column( 
            pl.when(pl.col(col).is_null())
            .then( fuzzyValues )
            .otherwise(pl.col(col)).alias(col)
        )
        
    print(df.head())
        
    return df

def identify_threshold_level(dfWideNew, force_calculation = False, clusterColumns=None):

    if clusterColumns is None:
        clusterColumns = [x for x in dfWideNew.columns if x.startswith("cluster.")]
    else:
        for x in clusterColumns:
            assert x in dfWideNew.columns
        
    explDF = dfWideNew.clone()
    #explDF.reset_index(drop=True, inplace=True)

    for thresholdIteration in range(0, 11, 1):
        currentThreshold = thresholdIteration * 0.05

        explDF = dfWideNew.clone()
        #explDF.reset_index(drop=True, inplace=True)

        print("Current Threshold", currentThreshold)
        finished = True

        for cidx, ccol in enumerate(clusterColumns):
            print(ccol)
            print(explDF.shape)

            initialCount = explDF.shape[0]

            explDF = explDF.explode(ccol)


            clDF = explDF.select(pl.col(ccol))

            binColName = "bin.{}".format(ccol)
            mfColName = "mf.{}".format(ccol)

            clDF.columns = [binColName, mfColName]

            clDF.with_column(pl.col(binColName).cast(pl.Categorical))
            clDF.with_column(pl.col(mfColName).cast(pl.Float32))

            explDF = pl.concat([explDF, clDF])
            explDF = explDF.filter(pl.col(mfColName) > currentThreshold)
            explDF.drop(ccol, inplace=True)

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
            
            print("Finished")
            break

    return explDF, currentThreshold



class CustomFuzzyVar(FuzzyVariable):

    def __init__(self, universe, label):
        super().__init__(universe, label)
        self.__name__ = "CustomFuzzyVar"

    def view(self, title=None, *args, **kwargs):
        """""" + FuzzyVariableVisualizer.view.__doc__
        fig, ax = FuzzyVariableVisualizer(self).view(*args, **kwargs)
        if not title is None:
            ax.set_title(title)
        fig.show()

    def __str__(self) -> str:
        return "CFuzzyVar [{}]".format("")


    def membership(self, value):        
        return [fuzz.interp_membership(self.universe, self.terms[x].mf, value) for x in self.terms]


    def automf(self, number=5, names=None, centers=None, shape="tri"):

        assert(shape in ("tri", "gauss","crisp"))
        assert(len(names) == number)

        if not centers is None:
            assert(len(centers) == len(names))

        limits = [self.universe.min(), self.universe.max()]
        universe_range = limits[1] - limits[0]
        widths = [universe_range / ((number - 1) / 2.)] * int(number)


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

        widths.append(2*(centers[-1]-centers[-2]))
            

        abcs = [[c - w / 2, c, c + w / 2] for c, w in zip(centers, widths)]
        cws = [[c, w/2] for c, w in zip(centers, widths)]  

        # Clear existing adjectives, if any
        self.terms = OrderedDict()
        unscaledValues = np.array([0.0]*len(self.universe))
   

        for name, abc, center_width in zip(names, abcs, cws):
            if shape == "tri":
                unscaledValues += fuzz.trimf(self.universe, abc)
            elif shape == "gauss":
                unscaledValues += fuzz.gaussmf(self.universe, center_width[0], center_width[1])
            elif shape == "crisp":
                pass
            else:
                raise ValueError("Shape not implemented" + str(shape))

        # Repopulate
        for e, (name, abc, center_width) in enumerate(zip(names, abcs, cws)):

            if shape == "tri":
                values = fuzz.trimf(self.universe, abc)/unscaledValues
            elif shape == "gauss":
                values = fuzz.gaussmf(self.universe, center_width[0], center_width[1])/unscaledValues
            elif shape == "crisp":  
                values= np.array( [ 1 if (x>=(center_width[0]-center_width[1]/2) and x<(center_width[0]+center_width[1]/2)) else float("nan") for x in self.universe ] )
            else:
                raise ValueError("Shape not implemented" + str(shape))
                
            half = int(len(values)/2)

            closest_to_center=self.universe.tolist()[min(range(len(self.universe.tolist())), key = lambda i: abs(self.universe.tolist()[i]-abc[1]))]
            index_center = self.universe.tolist().index(closest_to_center)

            if e == 0:
                values[:index_center][np.isnan(values[:index_center])] = 1    
            elif e == len(names)-1:
                values[index_center:][np.isnan(values[index_center:])] = 1
            
            values[np.isnan(values)] = 0
            
            self[name] =  values


        unscaledValues = np.array([0.0]*len(self.universe))

        for name in self.terms:

            unscaledValues += self.terms[name].mf

        assert(abs(sum(unscaledValues)-len(self.universe)) < 1)



def to_arg_max(df, exprMFs):

    emptyElem = [0] * len(exprMFs.terms)

    dfc = df.copy()
    dfc = dfc.applymap(lambda x: x if not x is np.nan else emptyElem)

    dfIdx = dfc.applymap( lambda x: np.argmax(x))
    dfValue = dfc.applymap( lambda x: np.max(x))
    dfValue = dfValue.fillna(0)

    return dfIdx, dfValue



def filter_weightSequence(weightSequence,cutoff=None):
    if not cutoff is None and cutoff > 0:
        return [x for x in weightSequence if x[1][-1] > cutoff]
    else:
        return weightSequence


def top_weightSequence(weightSequence,top=10):
    return [weightSequence[i] for i in range(top)]
def confusion_matrix(Trues,Results,All,outfile=None):
    Trues=set(Trues)
    Results=set(Results)
    All=set(All)

    tp=len(Trues.intersection(Results))
    fp=len(Results)-tp
    tn=len(All)-len(Trues.union(Results))
    fn=len(Trues)-tp

    confusion_matrix=[[tp,fp],[fn,tn]]

    if not outfile is None:
        L=pd.DataFrame(confusion_matrix)
        L.index=[1,0]
        L.columns=[1,0]
        pd.DataFrame(L).to_csv(outfile)
    
    return(confusion_matrix)

class FlowAnalysis:

    @classmethod
    def make_fuzzy_concepts(cls, exprData, mfLevels, centers, clusterColName, meancolName, mfLevelsMirrored, stepsize=None, shape="tri", series = None, perSeriesFuzzy=False, **kwargs):
        
        exprMFs = {}
        availableClusters =  set([x for x in exprData.get_column(clusterColName)])
        
        if "max.cluster" in exprData.columns:
            minValue = np.floor(exprData.select(pl.col("max.cluster")).min()[0])[0][0]
            maxValue = np.ceil(exprData.select(pl.col("max.cluster")).max()[0])[0][0]
        else:
            minValue = np.floor(exprData.select(pl.col(meancolName)).min()[0])[0][0]
            maxValue = np.ceil(exprData.select(pl.col(meancolName)).max()[0])[0][0]


        minValue=minValue-1
        maxValue=maxValue+1
        if not perSeriesFuzzy:
            
            print(kwargs)
            print(minValue, maxValue)
            
            if (not "centerMode" in kwargs and centers is None) or kwargs.get("centerMode", None) == "minmax":

                centers = np.linspace(minValue, maxValue, len(mfLevels))
                
                #widths = [universe_range / ((number - 1) / 2.)] * int(number)                                
            elif kwargs.get("centerMode", None) == "quantile_ends":
                
                quantileLimits = kwargs["centerQuantiles"]
                
                limits = np.quantile(exprData.select(pl.col(meancolName)), quantileLimits)
                
                limitRange = limits[1]-limits[0]
                limitStep = limitRange / len(mfLevels)
                centers = [limits[0] + x for x in range(0, limitRange, limitStep)]
                
                print("Limit Range", limitRange)
                print("Limit Step", limitStep)
            elif kwargs.get("centerMode", None) == "quantiles":
                
                quantileLimits = kwargs["centerQuantiles"]
                
                values = exprData.select(pl.col(meancolName))
                values = np.array(values)
                
                print(quantileLimits)
                print(values[:10])
                
                limits = np.quantile(values, quantileLimits)
                
                if not len(limits) == len(mfLevels):
                    raise ValueError("centerQuantiles length must be mfLevels length!")
                
                #minValue = limits[0]
                #maxValue = limits[-1]
                
                centers = [x for x in limits]
                                
            print("centers", centers)
            
            if mfLevelsMirrored:
                absValue = max(abs(minValue), abs(maxValue))
                minValue = -absValue
                maxValue = absValue
                assert(len(mfLevels) % 2 == 1)
            
            if stepsize is None:
                stepsize = min(0.1, (maxValue-minValue)/200)
                print("Fuzzy Step Size", stepsize)

            print("Creating Universe Range", minValue, "->", maxValue, "with step size", stepsize)


            exprMF = CustomFuzzyVar(np.arange(minValue, maxValue, stepsize), 'exprMFs')
            exprMF.automf(len(mfLevels), names=mfLevels, centers=centers, shape=shape) 
            
            for series in availableClusters:
                exprMFs[series] = exprMF
            
            exprMF.view("All MFs")
            
        else:
            
            # per series fuzzyfication
            for series in availableClusters:
                
                subsetExprData = exprData.filter(pl.col(clusterColName) == series)
                
                subsetMFs = cls.make_fuzzy_concepts(exprData=subsetExprData, mfLevels=mfLevels, centers=centers, clusterColName=clusterColName, meancolName=meancolName, mfLevelsMirrored=mfLevelsMirrored, stepsize=stepsize, shape=shape, series=series, perSeriesFuzzy=False, **kwargs)
                
                exprMFs[series] = subsetMFs
                exprMFs.view(series)
                plt.show()
                plt.close()

        return exprMFs


    @classmethod
    def fuzzify_exprvalues(cls, indf:pl.DataFrame, series = None, perSeriesFuzzy=False, mfLevels = ["NO", "LOW", "med", "HIGH"], mfLevelsMirrored=False, centers=None,symbol_column="gene", meancolName="mean.cluster", sdcolName="sd.cluster", exprcolName="expr.cluster", clusterColName="cluster", shape="tri", stepsize=None,combineOverState=False, **kwargs):

        exprData = indf.clone()

        exprMFs = cls.make_fuzzy_concepts(exprData, mfLevels, centers, clusterColName, meancolName, mfLevelsMirrored, stepsize=stepsize, shape=shape, series=series, perSeriesFuzzy=perSeriesFuzzy, **kwargs)
        meanExprCol = exprData.columns.index(meancolName)
        clusterCol = exprData.columns.index(clusterColName)

        if not exprcolName is None:
            exprCountCol = exprData.columns.index(exprcolName)
        else:

            exprData=exprData.with_column(pl.lit(1).alias("cell_expr"))
            exprCountCol = exprData.columns.index("cell_expr")
            exprcolName ="cell_expr"

    
        sdExprCol=-1
        if not sdcolName is None:
            sdExprCol = exprData.columns.index(sdcolName)
        else:
            sdExprCol=None

        print("Mean Expr", meancolName, "col", meanExprCol)
        print("Expr Count", exprcolName, "col", exprCountCol)
        print("SD", sdcolName, "col", sdExprCol)
        print("Cluster", clusterColName, "col", clusterCol)
        print("Combining over state: " ,combineOverState)


        df = exprData.clone()
                
        availableClusters =  set([x for x in exprData.get_column(clusterColName)])
        fuzzyOuts = []
        for seriesName in availableClusters:
            
            indf = df.filter(pl.col(clusterColName) == seriesName)

            if shape=="crisp":
                print("Crisp mode, different assignment")
      
                # weighted mean is calculated with all (not only expressed) cells
                indf = indf.with_columns(
                    (pl.col(meancolName)*pl.col(exprcolName)).alias(meancolName)
                )

                seriesOut = indf.select(
                    pl.struct([meancolName]).apply(lambda x:
                        distribution_to_crisp(x[meancolName], exprMFs[seriesName], threshold=0.0)
                        ).alias("fuzzy.mfs")
                )         
            else:
                if not sdcolName is None:          
                    seriesOut = indf.select(
                        pl.struct([meancolName, sdcolName, exprcolName]).apply(lambda x:
                            distribution_to_fuzzy(x[meancolName], x[sdcolName], x[exprcolName], exprMFs[seriesName], threshold=0.0)
                            ).alias("fuzzy.mfs")
                    )
                else:                    
                    #print(len(exprMFs[seriesName].terms))
                    
                    seriesOut = indf.select(
                        pl.struct([meancolName, exprcolName]).apply(lambda x:
                            distribution_to_fuzzy(x[meancolName], None, x[exprcolName], exprMFs[seriesName], threshold=0.0)
                            ).alias("fuzzy.mfs")
                    )    
            if combineOverState == True:
                
                #Here values are combined for each symbol_col + cluster entry
                seriesOut_wide=(
                    seriesOut
                    .with_row_count('id')
                    .explode("fuzzy.mfs")
                    .with_column(
                        ("FV_" + pl.arange(0, pl.count()).cast(pl.Utf8).str.zfill(2))
                        .over("id")
                        .alias("col_nm")
                    )
                    .pivot(
                        index=['id'],
                        values="fuzzy.mfs",
                        columns='col_nm',
                    )
                )

                new_indf = indf.hstack(seriesOut_wide)
                                

                FV_columns=list(filter(lambda x:'FV_' in x, new_indf.columns))
                new_indf=new_indf.groupby([symbol_column,clusterColName], maintain_order=True).agg([
                    pl.col(FV_columns).mean()
                ])


                new_indf=(new_indf
                .melt(
                    id_vars = [symbol_column,clusterColName], 
                    variable_name = 'FV',
                    value_name="fuzzy.mfs"
                    )
                )

                new_indf=new_indf.groupby([symbol_column,clusterColName], maintain_order=True).agg(pl.col("fuzzy.mfs"))

                seriesOut=new_indf.select(pl.col('fuzzy.mfs'))
                if sdcolName is None:
                    indf=indf.groupby([symbol_column,clusterColName], maintain_order=True).agg([
                        pl.col(meancolName),
                        pl.col(exprcolName)
                        ])
                else:
                    indf=indf.groupby([symbol_column,clusterColName], maintain_order=True).agg([
                        pl.col(meancolName),
                        pl.col(sdcolName),
                        pl.col(exprcolName)
                        ])      
            fuzzyOuts.append((indf, seriesOut))
            
        allExpr = pl.concat([x[0] for x in fuzzyOuts], how="vertical")
        allFuzzy = pl.concat([x[1] for x in fuzzyOuts], how="vertical")

        df = pl.concat([allExpr, allFuzzy], how="horizontal")           
        dfWide = to_homogeneous(toWideDF(df,symbol_column), exprMFs)
        
        return dfWide, exprMFs



    @classmethod
    def exprDF2LongDF(cls, indf:pl.DataFrame, seriesOrder = None, mfLevels = ["NO", "LOW", "med", "HIGH"], mfLevelsMirrored=False, centers=None, meancolName="mean.cluster", sdcolName="sd.cluster", exprcolName="expr.cluster", shape="tri", stepsize=None):

        return cls.fuzzify_exprvalues(indf, seriesOrder=seriesOrder, mfLevels=mfLevels, mfLevelsMirrored=mfLevelsMirrored, centers=centers, meancolName=meancolName, sdcolName=sdcolName, exprcolName=exprcolName, shape=shape, stepsize=stepsize)

        
    
    @classmethod
    def to_vwide(cls, indf, mfFuzzy):
        
        clusterCols = [x for x in indf.columns if x.startswith("cluster.")]
        
        for col in clusterCols:
            assert(col in mfFuzzy)
        
        outDF = indf.clone()
        
        for col in clusterCols:
            outDF = outDF.with_column(
                pl.struct([col]).apply(lambda x:
                                dict(zip( [x+".{}".format(col) for x in mfFuzzy[col].terms] , x[col] ))
                                ).alias("fuzzy.mfs")
            ).unnest("fuzzy.mfs")
            
            outDF = outDF.drop(col)
            
        return outDF

    @classmethod
    def toFlowsDF(cls, indf):

        explDF = indf.clone()
        binColumns = [x for x in explDF.columns if x.startswith("bin.")]
        binColumnsIndex = [explDF.columns.index(x) for x in binColumns]
        print(binColumns)

        mfColumns = [x.replace("bin.", "mf.") for x in binColumns]
        print(mfColumns)
        
        df=df.with_column(df.apply( lambda row: tuple( [row[x] for x in binColumnsIndex] )).to_series().alias("group.flow"))
        df=df.with_column(pl.struct(mfColumns).apply( lambda row: np.prod(row.values())).to_series().alias("mf.flow"))
        
        allgroups = list(set(df.select(pl.col("group.flow"))))
        
        flow2id = {}
        for idx, flow in enumerate(allgroups):
            flow2id[flow] = idx
        
        df=df.with_column(pl.struct(["group.flow"]).apply( lambda row: allgroups[row["group.flow"]] ).to_series().alias("id.flow"))

        return explDF


    def __init__(self, flows, symbol_column, series2name, exprMF):
        
        nExprMF = {}
        for x in exprMF:
            idx = x
            if idx.startswith("cluster."):
                idx = idx.replace("cluster.", "")
            nExprMF[idx] = exprMF[x]
            
        for x in series2name:
            for term in nExprMF[x[0]].terms:
                checkCol = "{}.cluster.{}".format(term, x[0])
                if not checkCol in flows.columns:
                    print("Missing column", checkCol)
                    raise ValueError("Missing column {}".format(checkCol))
        
        
        self.flows = flows
        self.seriesOrder = [x[0] for x in series2name]
        self.series2name = {x[0]: x[1] for x in series2name}
        self.symbol_column = symbol_column
        self.exprMF = nExprMF
        
        self.levelOrder = {}
        for state in self.exprMF:
            self.levelOrder[state] = [x for x in self.exprMF[state].terms]
            
        notEqualLevelOrder = False
        for state1 in self.exprMF:
            for state2 in self.exprMF:
                if len(set(self.levelOrder[state1]).difference(self.levelOrder[state2])):
                    notEqualLevelOrder = True        
        if notEqualLevelOrder:
            print("WARNING: Your level orders are not equal. Scales in some plots may not reflect this appropriately.")
        

    @property
    def flowid2flow(self):
        
        if not hasattr(self,"_flowid2flow"):
            print("Creating FlowIDs")
            self._flowid2flow = {}
            
            def createFlows(states):
                if len(states) == 0:
                    return [[]]
                
                rems = createFlows(states[1:])
                curstate = states[0]
                rets = []
                for level in self.levelOrder[curstate]:
                    for exflow in rems:
                        flow = [(curstate, level)] + exflow
                        rets.append(flow)                
                return rets               
            
            for flow in createFlows(self.seriesOrder):
                self._flowid2flow[len(self._flowid2flow)] = flow
                
                
        return self._flowid2flow

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


    def create_series_color_map(self, seriesColors, mode="scaling"):
        if seriesColors is None:
            seriesColors = {}
            for si, x in enumerate(self.series2name):
                seriesColors[self.series2name[x]] = plt.get_cmap("viridis")(si/len(self.series2name))
        seriesColorMap = {}
        for x in seriesColors:
            scmap = SankeyPlotter.createColorMap(seriesColors[x], mode)
            seriesColorMap[x] = scmap
            
        return seriesColorMap

    def plot_flows(self, use_flows = None, figsize=None, outfile=None, min_flow=None, min_gene_flow=None, transformCounts = lambda x: np.sqrt(x), verbose=False, seriesColors=None, colorMode="scaling"):

        if use_flows is None:
            use_flows = [x for x in self.flowid2flow]

        weightSequence = self._to_weight_sequence( flows=self.flows, use_flows=use_flows)

        weightSequence = filter_weightSequence(weightSequence, min_flow)          
                    
        if verbose:
            for x in weightSequence:
                print(x)

        
        seriesColorMap = self.create_series_color_map(seriesColors, colorMode)

        SankeyPlotter._make_plot(weightSequence, self.series2name, self.levelOrder, self.seriesOrder, transformCounts=transformCounts, fsize=figsize, outfile=outfile, verbose=verbose, seriesColorMap=seriesColorMap)


    def plot_coarse_flows(self, use_flows = None, genes=None,figsize=None, outfile=None, min_flow=None,  transformCounts = lambda x: x, verbose=False,specialColors=None,sns_palette="icefire", seriesColors=None, colorMode="scaling",title=None):
        
        
        if genes is None:
            flowDF=self.flows
        else:
            flowDF=self.flows.filter(pl.col(self.symbol_column).is_in(genes) )


        allSeries = [x for x in self.series2name]
        Cflowid2flow = {}
        for i in range(len(allSeries) - 1):
            
            srcSeries = allSeries[i]
            tgtSeries = allSeries[i+1]
        
            for comb in list(itertools.product(
                reversed(self.levelOrder[srcSeries]),
                reversed(self.levelOrder[tgtSeries])
                )):
                largeComp = [x for x in zip([allSeries[i],allSeries[i+1]], comb)]          
                Cflowid2flow[len(Cflowid2flow)] = largeComp
        
        if use_flows is None:
            use_flows = [x for x in Cflowid2flow]
        weightSequence = []

        for fgid in Cflowid2flow:

            if not fgid in use_flows:
                continue

            flowCols=["{}.cluster.{}".format(fclass, state) for state, fclass  in Cflowid2flow[fgid] ]
            flowScoreDF = flowDF.select(
                    pl.struct(flowCols).apply(lambda x: np.prod(list(x.values()))).alias("pwscore")
                    )
            flowScore = flowScoreDF.sum()[0,0]


            if not transformCounts is None:
                fgid = transformCounts(fgid)

            outlist = list(Cflowid2flow[fgid])
            outlist.append(flowScore )
            weightSequence.append(
                        (fgid, outlist)
                    )
        if specialColors is None:
            indices=[w[0] for w in weightSequence ]
            startingnodes=[w[1][0][1] for w in weightSequence ]
            #maxLevels = max([len(self.levelOrder[x]) for x in self.levelOrder])
            levelColors=dict.fromkeys(startingnodes)
            for i,k in enumerate(list(levelColors.keys())):
                levelColors[k] =i 
            maxLevels=len(levelColors.keys())

            colours=sns.color_palette(sns_palette,maxLevels)
            # color each state transition:             c = [colours[j] for j in range(maxLevels) for i in range(int(len(indices)/maxLevels))  ]

            # color end node transition:               c = [colours[j] for i in range(int(len(indices)/maxLevels)) for j in range(maxLevels) ]
            # or  [colours[k] for j in range(len(self.seriesOrder)-1)  for i in range(maxLevels) for k in range(maxLevels) ]

            # color every start node:            colours=sns.color_palette(sns_palette,max(int(len(indices)/maxLevels),maxLevels))     c = [colours[j] for j in range(int(len(indices)/maxLevels)) for i in range(maxLevels)  ]

            # color start node transition:

            #c = [colours[i] for j in range(len(self.seriesOrder)-1)  for i in range(maxLevels) for k in range(maxLevels) ]

            c=[colours[levelColors[s]] for s in startingnodes]
            specialColors=pd.DataFrame({
                'values':c},
                index=indices).to_dict()['values']
        weightSequence=filter_weightSequence(weightSequence,cutoff=min_flow)
        new_indices=[w[0] for w in weightSequence ]
        specialColors={k: v for k, v in specialColors.items() if k in new_indices}
        
        
        if not seriesColors is None:
            seriesColorMap = self.create_series_color_map(seriesColors, colorMode)
        else:
            seriesColorMap=None
        SankeyPlotter._make_plot(weightSequence, self.series2name, self.levelOrder, self.seriesOrder, specialColors=specialColors, seriesColorMap=seriesColorMap, independentEdges=True,title=title)



    def hist_level_membershipsum(self):

        flowDF=self.flows
        colmeans=flowDF.select(pl.col(pl.Float64)).transpose(include_header=True).with_column(
            pl.fold(0, lambda acc, s: acc + s, pl.all().exclude("column")).alias("horizontal_sum")
        )
                
        order1=[self.levelOrder[state].index(level) for level,state in [(x.split('.')[0],x.split('.')[2])  for x in colmeans.select(['column']).to_series()]]
        order2=[self.series2name[i] for i in [x.split('.')[2] for x in colmeans.select(['column']).to_series()]]

        colmeans=colmeans.with_column(pl.Series(name="order1", values=order1))
        colmeans=colmeans.with_column(pl.Series(name="order2", values=order2))
        
        colmeans=colmeans.sort("order2")
        
        fig, ax = plt.subplots()
        bar=sns.barplot(data=colmeans.to_pandas(), x="order1", y="horizontal_sum", hue="order2")
        
        xtickLabels = []
        for tickpos in ax.get_xticks():
            posLabels=[]
            for state in sorted(self.levelOrder):
                posLabels.append( "{state}: {label}".format(state=self.series2name[state], label=self.levelOrder[state][tickpos]) )
            xlabel = "\n".join(posLabels)
            xtickLabels.append(xlabel)
            
        ax.set_xticklabels(xtickLabels)
        
        ax.set_xlabel("Levels")
        ax.set_ylabel("Membership sum")
        ax.legend(title='State')
        fig.show()
        
    def plot_state_memberships(self,genes,name="", cluster_genes=False, outfile=None, limits=(0,1), annot=True, annot_fmt=".2f", prefix="Cluster", verbose=False, figsize=(6,6), font_scale=0.4):
        
        
        #firstLevelOrder = self.levelOrder[list(self.levelOrder.keys())[0]]
        
        filtered_flow=self.flows.filter(pl.col(self.symbol_column).is_in(genes) )

        pd_filtered_flow=pd.DataFrame(filtered_flow[:,1:], columns=filtered_flow[:,0].to_pandas().tolist(), index=filtered_flow[:,1:].columns)      
        
        # OR sort by states; switch hlines for correct white lines
        #pd_filtered_flow=pd_filtered_flow.sort_index(ascending=False)
        pd_filtered_flow["orderLevel"]=[self.levelOrder[state].index(level) for level,state in [(x.split('.')[0],x.split('.')[2]) for x in pd_filtered_flow.index]]
        pd_filtered_flow["orderState"]=[self.seriesOrder.index(i) for i in [x.split('.')[2] for x in pd_filtered_flow.index]]
        pd_filtered_flow= pd_filtered_flow.sort_values(["orderState", "orderLevel"])
        #pd_filtered_flow = pd_filtered_flow.drop(["order1","order2"], axis=1)
        
        col_order = ["orderLevel", "orderState"]+[x for x in genes if x in pd_filtered_flow.columns]
        pd_filtered_flow = pd_filtered_flow[col_order]
        
        if verbose:
            print(pd_filtered_flow)
        
        
        newIndex = []
        
        for ri, row in pd_filtered_flow.iterrows():
            
            rorder1 = int(row["orderLevel"])
            rorder2 = int(row["orderState"])
            
            state = self.series2name[self.seriesOrder[rorder2]]
            level = self.levelOrder[ self.seriesOrder[rorder2] ][rorder1]
            
            if prefix is None:
                newIndex.append("{} {}".format(state, level))
            else:
                newIndex.append("{} {} ({})".format(prefix, state, level))
        
        
        # number of elements until hline
        stateCounts = []
        states = list(pd_filtered_flow["orderState"])
        lastState = states[0]
        curCount = 1
        for x in states[1:]:
            if x == lastState:
                curCount += 1
            else:
                stateCounts.append(curCount)
                lastState = x
                curCount = 1
        
        stateCounts=np.cumsum(stateCounts)
        pd_filtered_flow = pd_filtered_flow.drop(["orderLevel","orderState"], axis=1)
        
        #for x in pd_filtered_flow.index:
        #    x = x.split(".")
        #    
        #    if prefix is None:
        #        newIndex.append("{} {} ({}) ".format(".".join(x[1:len(x)-1]).title(), x[-1], x[0]))
        #    else:
        #        newIndex.append("{} {} ({}) ".format(prefix, x[-1], x[0]))
            
        pd_filtered_flow.index = newIndex
                
        plt.figure(figsize=figsize)
        sns.set(font_scale=font_scale)
        
        pd_filtered_flow = pd_filtered_flow.iloc[::-1]
        
        g=sns.clustermap(pd_filtered_flow, row_cluster=False, col_cluster=cluster_genes, vmin=limits[0], vmax=limits[1], annot=annot, fmt=annot_fmt, cbar_pos=None, dendrogram_ratio=0.1)
        
        g.ax_col_dendrogram.set_visible(False)
        
        g.fig.suptitle(name) 
        g.ax_heatmap.yaxis.set_ticks_position("left")
        g.ax_heatmap.hlines(
            stateCounts,
            *g.ax_heatmap.get_xlim(),
            colors="white")
        #ax.set(title=name)
        #ax.hlines([[ int(pd_filtered_flow.shape[0]/len(self.levelOrder))  *x for x in range(len(self.levelOrder))]], *ax.get_xlim(),colors="white")

        if not outfile is None:
            plt.savefig(outfile + ".png", bbox_inches='tight')
            plt.savefig(outfile + ".pdf", bbox_inches='tight')


        plt.show()
        plt.close()
        
        return g

    def plot_genes(self, genes, figsize=None, outfile=None, min_flow=None, min_gene_flow=None, use_flows=None, transformCounts=lambda x: np.sqrt(x)):

        if not isinstance(genes, (tuple, list)):
            genes = [genes]

        #useFlows = self.flows.filter( self.flows["gene"].to_pandas().isin(genes).tolist())
        useFlows = self.flows.filter(pl.col(self.symbol_column).is_in(genes) )

        weightSequence = self._to_weight_sequence( flows=useFlows, use_flows=use_flows)
        weightSequence=filter_weightSequence(weightSequence,cutoff=min_flow)

        SankeyPlotter._make_plot(weightSequence, self.series2name, self.levelOrder, self.seriesOrder, transformCounts=transformCounts, fsize=figsize, outfile=outfile)


    def highlight_genes(self, genes, figsize=None, outfile=None, min_flow=None, min_gene_flow=None, transformCounts=lambda x: np.sqrt(x)):

        if not isinstance(genes, (tuple, list)):
            genes = [genes]
            
            
        bgData = self.flows.filter( ~pl.col(self.symbol_column).is_in(genes) )
        bgWeightSequence = self._to_weight_sequence( flows=bgData, use_flows=None, min_gene_flow=min_gene_flow)
        bgWeightSequence = filter_weightSequence(bgWeightSequence, min_flow)
        bgWeightSequence = filter_weightSequence(bgWeightSequence, min_flow)          

        survivedFlows = [x[0] for x in bgWeightSequence]
        
        fgData = self.flows.filter( pl.col(self.symbol_column).is_in(genes) )
        fgWeightSequence = self._to_weight_sequence( flows=fgData, use_flows=None, min_gene_flow=min_gene_flow, flowIDMod=lambda x: x*(-1))
        
        fgWeightSequence = [x for x in fgWeightSequence if x[0]*(-1) in survivedFlows]
        fgWeightSequence = filter_weightSequence(fgWeightSequence, min_flow)          

        specialColors = {x[0]: "red" for x in fgWeightSequence}
        
        SankeyPlotter._make_plot(bgWeightSequence+fgWeightSequence, self.series2name, self.levelOrder, self.seriesOrder, specialColors=specialColors, transformCounts=transformCounts, fsize=figsize, outfile=outfile)


    def visualize_genes(self, genes, figsize=None, min_flow=None, min_gene_flow=None, use_flows=None, title=None, outfile=None, score_modifier=lambda x: x):

        if not isinstance(genes, (tuple, list)):
            genes = [genes]
            
        if use_flows is None:
            use_flows = [x for x in self.flowid2flow]
            
            
        #bgData = self.flows.filter( ~pl.col("gene").is_in(genes) )
        bgWeightSequence = self._to_weight_sequence( flows=self.flows, use_flows=use_flows, min_gene_flow=min_gene_flow)
        
        bgWeightSequence = filter_weightSequence(bgWeightSequence, min_flow)          

        
        fgData = self.flows.filter( pl.col(self.symbol_column).is_in(genes) )
        fgWeightSequence = self._to_weight_sequence( flows=fgData, use_flows=use_flows, min_gene_flow=min_gene_flow, flowIDMod=None)#lambda x: x*(-1))
        
        for x in fgWeightSequence:
            x[1][-1] = score_modifier(x[1][-1])
            
        
        #print(fgWeightSequence)
        maxFlowValue = max([ x[1][-1] for x in fgWeightSequence])
        #print(maxFlowValue)

        cmap = get_cmap("cividis")
        norm = mpl.colors.Normalize(vmin=0, vmax=maxFlowValue)
        
        specialColors = {x[0]: cmap(x[1][-1]/maxFlowValue) for x in fgWeightSequence}
        SankeyPlotter._make_plot(bgWeightSequence, self.series2name, self.levelOrder, self.seriesOrder, specialColors=specialColors, transformCounts=lambda x: np.sqrt(x), fsize=figsize, cmap=cmap, norm=norm, outfile=outfile, title=title)
                       

    def plot_flow_memberships(self,use_flows, n_genes=30,color_genes=None, gene_include_patterns=[],gene_exclude_patterns=[], figsize=(2,5), outfile=None, plot_histogram=True,violin=False, labelsize=4,min_gene_flow=0.0001,parallel=True):
        flowDF=self.flows.clone()
        
        for gene_exclude_pattern in gene_exclude_patterns:
            flowDF = flowDF.filter(~pl.col(self.symbol_column).str.starts_with(gene_exclude_pattern))
        if(len(gene_include_patterns)>0):
            flowDF = flowDF.filter(pl.col(self.symbol_column).str.starts_with(pl.Series(gene_include_patterns)))
        
        flowScores = pl.DataFrame()

        
        def prepare_flow_calculation(fa, fgid):
            flowCols, flow = self._get_flow_columns(fgid)
            #filter(~pl.all(pl.col(flowCols) > min_gene_flow)).
            flowScore_perGene = flowDF.filter(pl.all(pl.col(flowCols) > min_gene_flow)).select([self.symbol_column,
                            pl.struct(flowCols).apply(lambda x: np.prod(list(x.values()))).alias("pwscore")
                        ]       
                        )
            flowScore_perGene.columns=[self.symbol_column,str(fgid) ]

        
            return(flowScore_perGene)
        
        bar = makeProgressBar()

        
        if parallel:
            
            print("Starting Event Loop")
            from joblib import Parallel, delayed, parallel_backend
            with parallel_backend('loky', n_jobs=20):
                results = Parallel()(delayed(prepare_flow_calculation)(self, fgid) for fgid in bar(use_flows))
            print("Event Loop Completed")
            
        else:
            results=[]
            for fgid in bar(use_flows):
                results.append(prepare_flow_calculation(self,fgid))


        for idx, res in enumerate(results):
                    
            if flowScores.shape == (0,0):
                flowScores=res
            else:
                flowScores=flowScores.join(res,on=self.symbol_column,how="outer")
        
        flowScores=flowScores.fill_null(0)  
        ## here also coloring bars content from flows would be possible
        flowScores_df= flowScores.select([self.symbol_column,
                        pl.struct(flowScores[:,1:].columns).apply(lambda x: np.sum(list(x.values()))).alias("pwscore")
                    ]       
        )

        flowScores_df=flowScores_df.sort(["pwscore",self.symbol_column],reverse=True)

        if plot_histogram:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={'height_ratios': [4, 1]}, sharex=True)
        else:
            fig, ax1 = plt.subplots(1, 1, figsize=figsize)
            ax2 = None

        n_genes = min(n_genes, flowScores_df.shape[0])

        flowScores_df_top=flowScores_df[:n_genes]
        colormap=['blue']*n_genes

         
        if color_genes:
            not_in_genes=flowScores_df_top.select(~pl.col(self.symbol_column).is_in(list(color_genes)))
            print("Found "+str(n_genes-sum(not_in_genes.to_series()))+" in "+str(n_genes))
            colormap=np.where(not_in_genes == True, 'red', colormap)[0]
            
        ax1.set_facecolor('#FFFFFF')

        
        # Style the grid.
        ax1.grid(which='major', color='#EBEBEB', linewidth=1.2)
        ax1.grid(which='minor', color='#EBEBEB', linewidth=0.6)
            
        #ax1.barh(range(n_genes),flowScores_df_top["pwscore"],color=colormap)

        ax1.hlines(y=range(n_genes), xmin = 0 , xmax = flowScores_df_top["pwscore"], color=colormap)
        ax1.plot(list(flowScores_df_top["pwscore"]), range(n_genes), "o")
        ax1.set_yticks(range(n_genes),flowScores_df_top[self.symbol_column])
        ax1.tick_params(axis="y",labelsize=labelsize)
        ax1.set_title("Top "+str(n_genes)+" memberships",fontsize = labelsize)
        ax1.invert_yaxis()

        # Only show minor gridlines once in between major gridlines.
        from matplotlib.ticker import AutoMinorLocator
        ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
        #ax2.yaxis.set_minor_locator(AutoMinorLocator(2))


        if not ax2 is None:
            flowScores_df_rounded=flowScores_df.with_column(
                pl.col("pwscore").round(1)
            )
            counted_scores=flowScores_df_rounded.select("pwscore").to_series().value_counts()

            if violin:
                sns.set_style('whitegrid')
                #sns.kdeplot(x=np.array(flowScores_df["pwscore"]),ax=ax2)
                sns.violinplot(x=np.array(flowScores_df["pwscore"]),ax=ax2,inner="stick")
                ax2.set_title("Flow membership distribution")
                ax2.set_xlim(0,1)
            else:
                bars=ax2.bar(counted_scores['pwscore'], counted_scores['counts'], width = 0.1)
                for rect in bars:
                    real_height=rect.get_height()
                    ax2.text(rect.get_x() + rect.get_width()/2., (real_height*.01).clip(min=2.5),
                            '%d' % int(real_height),c='w',
                            ha='center', va='bottom',rotation=90)
                ax2.set_title("Binned membership histogram",fontsize = labelsize)
                ax2.set_yscale('log')
                ax2.set_facecolor('white')
                ax2.grid(which='major', color='#EBEBEB', linewidth=1)
                ax2.grid(which='minor', color='#EBEBEB', linewidth=1)
                ax2.tick_params(axis="y",labelsize=labelsize)
                ax2.axis(ymin=1)

        fig.tight_layout()
        
        if not outfile is None:
            plt.savefig(outfile + ".png", bbox_inches='tight')
            plt.savefig(outfile + ".pdf", bbox_inches='tight')
        
        plt.show()
        
        return list(flowScores_df_top.select(pl.col(self.symbol_column)))[0].to_list(), flowScores_df, (ax1, ax2)

    def get_confusion_matrix(self,scores_df,Trues,num_true=None,outfile=None):
        if not num_true is None:
            scores_df=scores_df[:num_true]
        scores_df=scores_df.filter(pl.col("pwscore")>0.0)

        Results=list(scores_df.select(pl.col(self.symbol_column)))[0].to_list()
        All=list(self.flows.select(pl.col(self.symbol_column)))[0].to_list()
        cm=confusion_matrix(Trues,Results,All,outfile)


        print(cm)

    def plot_genes_membership(self,genes, n_genes=30,use_flows=None,plot_hist=False, outfile=None):
        genes_flow=self.flows.filter(pl.col(self.symbol_column).is_in(genes) )
        if use_flows is None:
            use_flows = [x for x in self.flowid2flow]

        weightSequence = self._to_weight_sequence( flows=genes_flow, use_flows=use_flows)

        fgid=[x[0] for x in weightSequence]
        membership=[x[1][-1] for x in weightSequence]

        flowScores_df=pl.DataFrame({'fgid':fgid, 'membership':membership})

        flowScores_df=flowScores_df.sort("membership",reverse=True)
        fig, (ax1, ax2) = plt.subplots(1, 2)
        
        n_genes = min(n_genes, flowScores_df.shape[0])
        
        ax2.barh(range(n_genes),flowScores_df["membership"][range(n_genes)])
        
        flowIDList = list(flowScores_df["fgid"][range(n_genes)])
               
        ax2.set_yticks(range(n_genes),[" -> ".join([str("+".join(x)) for x in self.flowid2flow[y]]) for y in flowIDList])
        
        
        ax2.tick_params(axis="y",labelsize=4)
        ax2.set_title("Top {} memberships".format(n_genes))
        ax2.invert_yaxis()

        flowScores_df=flowScores_df.with_column(
            pl.col("membership").round(1)
        )
        if plot_hist:
            counted_scores=flowScores_df.select("membership").to_series().value_counts()
            #print(counted_scores)
            ax1.bar(counted_scores['membership'], counted_scores['counts'], width = 10)
            ax1.set_title("Binned membership histogram")
        else:
            ax1.axis('off')

        
        if not outfile is None:
            plt.savefig(outfile + ".png", bbox_inches='tight')
            plt.savefig(outfile + ".pdf", bbox_inches='tight')
        
        plt.show()


    
    def _get_flow_columns(self, flowID):
        
        flow = self.flowid2flow[flowID]
        # flow = [('control', 'NO'), ('sepsis', 'LOW')]
        
        flowCols = ["{}.cluster.{}".format(y,x) for x,y in flow]
        #flowCols = ['NO.cluster.control', 'LOW.cluster.sepsis']
        
        return flowCols, flow

    def _to_weight_sequence(self, flows, use_flows=None, flowIDMod=None, min_gene_flow=None):
        """Generates weight-sequence for plotting flows

        Args:
            flows (pl.DataFrame): Flow-Dataframe used for plotting
            use_flows (list, optional): flows to consider. Defaults to None. If None, all Flows will be plotted
            flowIDMod (lambda, optional): if not None, applied to all flow ID results

        Returns:
            list: weight-sequence
        """
        

        if use_flows is None:
            use_flows = [x for x in self.flowid2flow]
        
        weightSequence = []
        for fgid in use_flows:            
            # (fgid, [(("WT", 2), ("KO", 0), 1), .... )]
            flowScore, flow = self._calculate_flow_score(flows, fgid, min_gene_flow=min_gene_flow)
            
            if not flowIDMod is None:
                fgid = flowIDMod(fgid)
            #print(fgid, flow, flowScore)
            
            outlist = list(flow)
            outlist.append(flowScore )
            weightSequence.append(
                (fgid, outlist)
            )
                                    
        return weightSequence

    def flow_finder( self, sequence, minLevels=None, maxLevels=None, verbose=True ):

        firstState = list(self.exprMF.keys())[0]
        seriesOrder = list(self.exprMF[firstState].terms.keys())

        if not minLevels is None:
            for x in minLevels:
                if not x is None:
                    assert(x in seriesOrder)
        if not maxLevels is None:
            for x in maxLevels:
                if not x is None:
                    assert(x in seriesOrder)


        matchingFIDs = set()

        for fid in self.flowid2flow:
            
            fgroup = self.flowid2flow[fid]

            acceptFlow = True
            
            for ci, comp in zip(range(0, len(sequence)), sequence):

                startIdx = seriesOrder.index(fgroup[ci][1])
                endIdx = seriesOrder.index(fgroup[ci+1][1])

                if not minLevels is None:
                    if (not minLevels[ci] is None and startIdx < seriesOrder.index(minLevels[ci])) or (not minLevels[ci+1] is None and endIdx < seriesOrder.index(minLevels[ci+1])):
                        acceptFlow=False
                        break

                if not maxLevels is None:
                    if (not maxLevels[ci] is None and startIdx > seriesOrder.index(maxLevels[ci])) or (not maxLevels[ci+1] is None and endIdx > seriesOrder.index(maxLevels[ci+1])):
                        acceptFlow=False
                        break

                if comp == "<":
                    if not startIdx < endIdx:
                        acceptFlow=False
                        break
                if comp == "<<":
                    if not startIdx+1 < endIdx:
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
                elif comp == ">>":
                    if not startIdx > endIdx+1:
                        acceptFlow=False
                        break
                elif comp == "=":
                    if not startIdx == endIdx:
                        acceptFlow=False
                        break
                elif comp == "~":
                    #circa
                    if not (startIdx == endIdx or startIdx == endIdx+1 or startIdx == endIdx -1) :
                        acceptFlow=False
                        break
                elif comp == "?":
                    #always accept
                    pass

                if not acceptFlow:
                    break

            if acceptFlow:
                if verbose:
                    print(fid, fgroup)
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
    
    def read_gaf_file(self, filepath):

        geneset2genes = defaultdict(set)
        geneset2name = {}

        with open(filepath) as fin:

            for line in fin:
                
                if line.startswith("!"):
                    continue

                line = line.strip().split("\t")
                pwName = line[4]
                pwID = line[4]
                pwGene = line[2].upper()
                
                geneset2name[pwID] = pwName
                geneset2genes[pwID].add(pwGene)
                
        output = {}
        for x in geneset2name:
            output[x] = (x, geneset2genes[x])

        return output

    def get_pathways(self, pathways_file):
        
        print("Loading pathways from", pathways_file)
        
        if pathways_file.endswith("gmt"):
            rp = self.read_gmt_file(pathways_file)
        elif pathways_file.endswith("gaf"):
            rp = self.read_gaf_file(pathways_file)
        else:
            raise ValueError("Invalid File Format")
        
        print("Identified", len(rp), "pathways")
        
        return rp


    def analyse_pathways(self, pathways_file="ReactomePathways.gmt", additional_pathways=None, use_flows=None, parallel=True):

        rp = self.get_pathways(pathways_file)
        
        if not additional_pathways is None:
            for pname, pgenes in additional_pathways:
                rp[pname] = (pname, pgenes)

        allDFs = []

        bar = makeProgressBar()

        if use_flows is None:
            use_flows = [x for x in self.flowid2flow]


        if parallel:
            
            
            def prepare_flow_analysis_df(fa, fgid):
                flowCols, _ = fa._get_flow_columns(fgid)
                flowDF = fa.flows.select(pl.col(flowCols + [fa.symbol_column]))
                        
                df = fa.analyse_genes_for_genesets(rp, flowDF, bgFlowDF=fa.flows, considerFlows=[fgid])
                df["fgid"] = fgid
                return df
            
            
            print("Starting Event Loop")
            from joblib import Parallel, delayed, parallel_backend
            
            with parallel_backend('loky', n_jobs=8):
                results = Parallel()(delayed(prepare_flow_analysis_df)(self, fgid) for fgid in use_flows)

                for idx, res in enumerate(results):
                    allDFs.append(res)
            
            print("Event Loop Completed")
            
            
        else:

            for fgid in bar(use_flows):
                
                flowCols, _ = self._get_flow_columns(fgid)
                flowDF = self.flows.select(pl.col(flowCols + [self.symbol_column]))
                        
                df = self.analyse_genes_for_genesets(rp, flowDF, bgFlowDF=self.flows, considerFlows=[fgid])
                df["fgid"] = fgid

                #print(df[df["pwsize"] > 1].sort_values(["pval"], ascending=True).head(3))

                allDFs.append(df)



        allFGDFs = pd.concat(allDFs, axis=0)
        allFGDFs = self._calculate_pvalues(allFGDFs)

        return allFGDFs

    def analyse_pathways_grouped(self, use_flows, pathways_file="ReactomePathways.gmt", additional_pathways=None, set_size_threshold=[ 2, 10, 50, 100]):

        rp = self.get_pathways(pathways_file)
        
        if not additional_pathways is None:
            for pname, pgenes in additional_pathways:
                rp[pname] = (pname, pgenes)


        flowCols = set()
        
        for fgid in use_flows:
            fc, _ = self._get_flow_columns(fgid)
            flowCols.update(fc)
            
        flowCols=list(flowCols)
            
        flowDF = self.flows.select(pl.col(flowCols + [self.symbol_column]))
    
        allFGDFs = self.analyse_genes_for_genesets(rp, flowDF, bgFlowDF=self.flows, considerFlows=use_flows)
        
        allFGDFs = self._calculate_pvalues(allFGDFs, set_size_threshold=set_size_threshold)
        
        return allFGDFs


    def _calculate_flow_score(self, flowDF, flowID, min_gene_flow=None):
        
        
        flowCols, flow = self._get_flow_columns(flowID)
        
        if flowDF.shape[0] == 0:
            return 0.0, flow

        for col in flowCols:
            if not col in flowDF.columns:
                return 0.0, flow
        
        flowScoreDF = flowDF.select(
            pl.struct(flowCols).apply(lambda x: np.prod(list(x.values()))).alias("pwscore")
        )
        
        if not min_gene_flow is None and min_gene_flow > 0.0:
            flowScoreDF = flowScoreDF.filter(pl.col("pwscore") > min_gene_flow)
        
        
        if flowScoreDF.shape[0] == 0:
            return 0, flow
        
        flowScore = flowScoreDF.sum()[0,0]
        
        return flowScore, flow

    def _calculate_pvalues(self, df, set_size_threshold=None):
        
        if set_size_threshold is None:
            set_size_threshold = [ 2, 10, 50, 100]
        
        inDF = df.copy()
        
        resultDFs = []
        lastThreshold=-1
        
        set_size_threshold = list(sorted(set(set_size_threshold)))
        
        if set_size_threshold[-1] < inDF.pwGenes.max():
            set_size_threshold.append( inDF.pwGenes.max()+1 )
        
        print("Calculating p-values for groups", set_size_threshold)
        
        for t in set_size_threshold:
            
            curDF = inDF[(inDF.pwGenes > lastThreshold) & (inDF.pwGenes <= t)].copy()
        
            pwC_mean = curDF[curDF.pw_coverage != 0].pw_coverage.mean()
            pwC_std = curDF[curDF.pw_coverage != 0].pw_coverage.std(ddof=0)
            
            curDF["pw_coverage_zscore"] = (curDF["pw_coverage"]-pwC_mean)/pwC_std
            curDF["pval"] = norm.sf(abs(curDF["pw_coverage_zscore"]))

            curDF.loc[curDF["pw_coverage_zscore"] < 0, "pval"] = 1.0
            lastThreshold = t
            
            resultDFs.append(curDF)

        outDF = pd.concat(resultDFs, axis=0)

        _ , elemAdjPvals, _, _ = multipletests(outDF["pval"], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        outDF["adj_pval"] = elemAdjPvals
        
        return outDF


    def analyse_genes_for_genesets(self, pathways, flowDF, bgFlowDF, considerFlows=None, populationSize=None):
        """
        
        pathways: pathway object
        genes: genes of interest

        populationSize: size of universe. if None, all genes in pathways will be chosen
        
        """
        
        if considerFlows is None:
            considerFlows = [x for x in self.flowid2flow]


        allPathwayGenes = set()
        for x in pathways:
            pwname, pwGenes = pathways[x]
            for gene in pwGenes:
                allPathwayGenes.add(gene)

        allPathwayGenes = list(allPathwayGenes.intersection(list(bgFlowDF.select(pl.col(self.symbol_column)).to_series())))
        
        #print("Pathway Genes", len(allPathwayGenes))
        #print("Measured Pathways Genes", len(allPathwayGenes))
        
               
        if populationSize is None:
            allPWGeneFlowDF = bgFlowDF.filter(pl.col(self.symbol_column).is_in(allPathwayGenes))

            totalFlowSumOfPathwayGenes = 0

            for flowID in self.flowid2flow:
                flowScore, _ = self._calculate_flow_score( allPWGeneFlowDF, flowID)
                totalFlowSumOfPathwayGenes += flowScore
            
            populationSize = totalFlowSumOfPathwayGenes

            #print("populationSize", populationSize)


        flowGenes = set(flowDF.select(self.symbol_column).to_series())       
        
        allFlowsScore = 0
        for flowID in considerFlows:
            flowScore, _ = self._calculate_flow_score( flowDF.filter(pl.col(self.symbol_column).is_in(allPathwayGenes)), flowID)
            #
            allFlowsScore += flowScore
            
        #print("AllFlowScore", allFlowsScore)
        
        outData = defaultdict(list)

        for pwID in pathways:

            pwName, pwGenes = pathways[pwID]
            pwGenes = list(pwGenes)
            inflow_inset = list(flowGenes.intersection(pwGenes))

            pathwayDF = flowDF.filter(pl.col(self.symbol_column).is_in(pwGenes))
            
            #print(pathwayDF.shape)
            #print(list(pathwayDF.select(self.symbol_column).to_series()))
        
            flowInPathwayScore = 0
            for flowID in considerFlows:
                flowScore, _ = self._calculate_flow_score( pathwayDF, flowID)
                flowInPathwayScore += flowScore


            if False:
                if abs(flowInPathwayScore) < 0 or abs(allFlowsScore) < 0:
                    chi2 = 0
                    pval = 1.0
                else:

                    #inflow inflow_inset
                    #not-inflow not-inflow_inset
                    
                    #TBL = rbind(c(x1,x2), c(n1-x1, n2-x2))
                    #x1 = 1680000;  n1 = 12000000
                    #x2 = 3;  n2 = 30
                    
                    # n1/n2 = populations (n1: global, n2: pathway)
                    # x1/x2 = ofinterest (x1: flow-score-all-pathways, flow-score-of-pathway)
                    

                    table=np.array([
                                    [allFlowsScore,flowInPathwayScore], #len(inflow_inset)
                                    [populationSize-allFlowsScore,len(pwGenes)-flowInPathwayScore] # 
                                ])
                    
                    if np.any(table < 0) or pwID == "hsa04740" or pwName == "hsa04740":
                        print(pwID)
                        print(table)
                        print(pwName)

                        print("pwFlow, x2", flowInPathwayScore)
                        print("pwGenes, n2", len(pwGenes))
                        print("allPwFlow, x1", allFlowsScore)
                        print("allPwGenes, n2", populationSize)
                        
                        #return None
                        
                    
                    chi2, pval, dof, expected=chi2_contingency(table, correction=True)
                    
                    if flowInPathwayScore <= expected[0][1]:
                        pval=1.0
                    
                    if pwName == "hsa04740":
                        print("expected", expected)

            # population: all genes
            # condition: genes
            # subset: pathway

            genes_coverage = flowInPathwayScore / allFlowsScore if allFlowsScore > 0 else 0
            pathway_coverage = flowInPathwayScore / len(pwGenes) if len(pwGenes) > 0 else 0

            outData["pwid"].append(pwID)
            outData["pwname"].append(pwName)
            outData["pwFlow"].append(flowInPathwayScore)
            outData["pwGenes"].append(len(pwGenes))
            outData["allPwFlow"].append(allFlowsScore)
            outData["allPwGenes"].append(populationSize)

            outData["pw_gene_intersection"].append(len(inflow_inset))
            outData["pw_coverage"].append(pathway_coverage)
            outData["genes_coverage"].append(genes_coverage)
            #outData["pval"].append(pval)
            #outData["chi2"].append(chi2)
            outData["mean_coverage"].append(pathway_coverage*genes_coverage)

        outdf = pd.DataFrame.from_dict(outData)       

        return outdf












    def custom_div_cmap(self, numcolors=11, name='custom_div_cmap',
                        mincol='blue', midcol='white', maxcol='red'):
        """ Create a custom diverging colormap with three colors
        
        Default is blue to white to red with 11 colors.  Colors can be specified
        in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
        """

        from matplotlib.colors import LinearSegmentedColormap 
        
        cmap = LinearSegmentedColormap.from_list(name=name, 
                                                colors = [x for x in [mincol, midcol, maxcol] if not x is None],
                                                N=numcolors)
        return cmap

    def plotORAresult( self, dfin, title, numResults=10, figsize=(10,10), outfile=None, sep=" ", entryformat="{x}{sep}({y}, pw_cov={z:.3f}/{s})"):
        #https://www.programmersought.com/article/8628812000/
        
        def makeTitle(colDescr, colID, colSize, setSize):
            out = []
            for x,y,z, s in zip(colDescr, colID, colSize, setSize):
                out.append(entryformat.format(x=x, sep=sep, y=y, z=z, s=s))

            return out


        df_raw = dfin.copy()

        # Prepare Data
        #determine plot type

        termIDColumn = "pwid"
        termNameColumn = "pwname"
        qvalueColumn = "adj_pval"
        df = df_raw.copy()

        print(df_raw.columns)
        print("GeneRatio" in df_raw.columns)
        print("BgRatio" in df_raw.columns)

        rvb = None

        color1 = "#883656"
        color2 = "#4d6841"
        color3 = "#afafaf"

        if "pw_coverage" in df_raw.columns and "genes_coverage" in df_raw.columns:
            #ORA
            rvb = self.custom_div_cmap(150, mincol=color2, maxcol=color3, midcol=None)
            colorValues = [rvb(x/df.pwGenes.max()) for x in df.pwGenes]

        else:
            raise ValueError()

        df['termtitle'] = makeTitle(df_raw[termNameColumn], df_raw[termIDColumn], df["pwFlow"],df["pwGenes"])


        #df.sort_values('adj_pval', inplace=True, ascending=True)
        df.reset_index()
        df = df[:numResults]
        colorValues = colorValues[:numResults]
        
        df = df.iloc[::-1]
        colorValues = colorValues[::-1]

        print(df.shape)

        if df.shape[0] == 0:
            return
        
        maxNLog = max(-np.log(df.adj_pval))
        maxLine = ((maxNLog// 10)+1)*10       
        
        # Draw plot
        fig, ax = plt.subplots(figsize=figsize, dpi= 80)
        ax.hlines(y=df.termtitle, xmin=0, xmax=maxLine, color='gray', alpha=0.7, linewidth=1, linestyles='dashdot')
        ax.vlines(x=-np.log(0.05), ymin=0, ymax=numResults, color='red', alpha=0.7, linewidth=1, linestyles='dashdot')
        
        sizeFactor = 10    
        scatter = ax.scatter(y=df.termtitle, x=-np.log(df.adj_pval), s=df.pwGenes*sizeFactor, c=colorValues, alpha=0.7, )

        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6, func=lambda x: x/sizeFactor)
        labels = [x for x in labels]

        # Title, Label, Ticks and Ylim
        ax.set_title(title, fontdict={'size':12})
        ax.set_xlabel('Neg. Log. Adj. p-Value', fontdict={'size':12})
        ax.set_yticks(df.termtitle)
        ax.set_yticklabels(df.termtitle, fontdict={'horizontalalignment': 'right'})
        
        plt.tick_params(axis='x', which="major", length=7, width=2)
        plt.tick_params(axis='x', which="minor", length=5, width=2)
        
        if (0-np.log10(df["adj_pval"].min())) > 50:
            ax.set_xscale('log')
            
            
        from matplotlib.ticker import ScalarFormatter
        ax.xaxis.set_major_formatter(ScalarFormatter())
        
        plt.grid(b=None)
        plt.tight_layout()
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=8)
        plt.setp(ax.get_xmajorticklabels(), fontsize=8)
        plt.setp(ax.get_xminorticklabels(), fontsize=8)

        if not outfile is None:
            plt.savefig(outfile + ".png", bbox_inches='tight')
            plt.savefig(outfile + ".pdf", bbox_inches='tight')
        
        plt.show()
