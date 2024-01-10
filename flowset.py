
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
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import AutoMinorLocator


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
    def _make_plot( cls, nodeWeigthSequence, series2name, levelOrder, seriesOrder, specialColors=None, nodeColors=None, fsize=None, transformCounts = lambda x: x, cmap=None, norm=None, outfile=None, title=None, verbose=False, linewidth=0.01, seriesColorMap=None, independentEdges=False, seriesFontsize=22, classFontsize=22):

        
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
            

            t = ax.text(nodePosition[0], nodePosition[1], nodeStr, transform=ax.transData, fontsize=classFontsize,rotation=90,
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
            
            rect = patches.FancyBboxPatch( (si-1.0*bwidth, -2 - 0.5*bheight), width=2*bwidth, height=bheight, facecolor=faceColorRGB, linewidth=0)
            rect.set_boxstyle("round", rounding_size=0.05, pad=0)
            ax.add_patch(rect)

            ax.text(si, -2, series2name[series], transform=ax.transData, fontsize=seriesFontsize,rotation=0,
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
            plt.title(title,fontsize = 30)
        
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
        
        #return L < 0.179
        return L < 0.3

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
            
            brightInColor = cls.scale_lightness(inColorRGB, 1.5)
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
            
        elif mode == "constant":
            
            return lambda x: inColorRGB
            
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

def toWideDF( df,symbol_column, cluster_column="cluster"):
    #dfPivot = pd.pivot(indf, index=["gene"], columns=["cluster"], values=["fuzzy_set"])
    #dfWide = dfPivot.copy()
    #dfWide.columns = dfWide.columns.droplevel(0)
    #dfWide.reset_index(inplace=True)
    #dfWide.reset_index(drop=True, inplace=True)
    
    dfPivot = df.pivot(values=["fuzzy.mfs"], index=symbol_column, columns=cluster_column)
    return dfPivot

def to_homogeneous(df:pl.DataFrame, exprMFs, is_foldchange=False):
        
    #print(fuzzyBins)
    #print(fuzzyValues)

    for col in exprMFs:
        
        print("to_homogeneous:", col)
        
        if not col in df.columns:
            df=df.with_columns( pl.lit(None).alias(col) )
                
        exprMF = exprMFs[col]
        
        if not is_foldchange:
            #fuzzyBins = [x for x in exprMF.terms]
            fuzzyValues = [float(x) for x in to_fuzzy(0, exprMF)]
        else:
            #fuzzyBins = [x for x in exprMF.terms]
            fuzzyValues = [float(x) for x in to_fuzzy(0, exprMF)]
                 
        df=df.with_columns( 
            pl.when(pl.col(col).is_null())
            .then( fuzzyValues )
            .otherwise(pl.col(col)).alias(col)
        )
                
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

            clDF.with_columns(pl.col(binColName).cast(pl.Categorical))
            clDF.with_columns(pl.col(mfColName).cast(pl.Float32))

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



import bisect

class Closest:

    """Assumes *no* redundant entries - all inputs must be unique"""

    def __init__(self, numlist=None, firstdistance=0):
        if numlist is None:
            numlist=[]
        self.numindexes = dict((val, n) for n, val in enumerate(numlist))
            #Only if unsorted
        self.nums = sorted(self.numindexes)
        self.firstdistance = firstdistance

    def append(self, num):
        if num in self.numindexes:
            raise ValueError("Cannot append '%s' it is already used" % str(num))
        self.numindexes[num] = len(self.nums)
        bisect.insort(self.nums, num)

    def rank(self, target):
        rank = bisect.bisect(self.nums, target)
        if rank == 0:
            pass
        elif len(self.nums) == rank:
            rank -= 1
        else:
            dist1 = target - self.nums[rank - 1]
            dist2 = self.nums[rank] - target
            if dist1 < dist2:
                rank -= 1
        return rank

    def closest(self, target):
        try:
            return self.numindexes[self.nums[self.rank(target)]]
        except IndexError:
            return 0

    def distance(self, target):
        rank = self.rank(target)
        try:
            dist = abs(self.nums[rank] - target)
        except IndexError:
            dist = self.firstdistance
        return dist

def to_fuzzy_fast(targets, fzy):
    a = fzy.universe
    cl = Closest(a)


    #for x in targets:
    #    print(x)
    #    print(cl.closest(x))


    indices=[cl.closest(x) for x in targets ]

    out=pl.DataFrame()
    out=out.with_columns(pl.Series(name='universe', values=fzy.universe) )
    for n,v in fzy.terms.items():

        out=out.with_columns(pl.Series(name=n, values=fzy[n].mf) )
    out=out.with_column( pl.all().round(3) )
    out=out.with_row_count()
    out = out.select(
        [
            pl.col("row_nr").cast(pl.Int64),
            pl.all().exclude("row_nr")
        ]
    )
    out=out.join(pl.DataFrame(indices),left_on='row_nr',right_on='column_0')
    out=out.drop(['row_nr','universe'])
    return out


from abc import ABC, abstractmethod

class AbstractFuzzifier(ABC):
    
    @abstractmethod
    def fuzzify(Self):
        pass
    
    @classmethod
    def exprDF2LongDF(cls, indf:pl.DataFrame, seriesOrder = None, mfLevels = ["NO", "LOW", "med", "HIGH"], mfLevelsMirrored=False, centers=None, meancolName="mean.cluster", sdcolName="sd.cluster", exprcolName="expr.cluster", shape="tri", stepsize=None):

        return cls.fuzzify_exprvalues(indf, seriesOrder=seriesOrder, mfLevels=mfLevels, mfLevelsMirrored=mfLevelsMirrored, centers=centers, meancolName=meancolName, sdcolName=sdcolName, exprcolName=exprcolName, shape=shape, stepsize=stepsize)

        
    
    @classmethod
    def to_vwide(cls, indf, mfFuzzy, meta_columns=["gene"]):
        
        clusterCols = [x for x in indf.columns if not x in meta_columns]
        
        for col in clusterCols:
            if not col in mfFuzzy:
                print("Expected col", col, "in fuzzy concepts.")
            assert(col in mfFuzzy)
        
        outDF = indf.clone()
        
        
        def listcol_to_cols(x):
            retDict = OrderedDict(zip( ["{}.{}".format(y, col) for y in mfFuzzy[col].terms] , x[col] ))
            return retDict
                                
        
        
        for col in clusterCols:
            
            structAlias = "{}.mfs".format(col)
            
            outDF = outDF.with_columns(
                pl.struct([col]).apply(listcol_to_cols).alias(structAlias)
            ).unnest(structAlias)
            
            outDF = outDF.drop(col)
            
        return outDF
    
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

        if not perSeriesFuzzy:
                    
            
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
    
    
class BaseFuzzifier(AbstractFuzzifier):
    
    def __init__(self) -> None:
        super().__init__()
        
        
    def fuzzify(self, indf:pl.DataFrame, exprMFs, mean_col, sd_col, expr_col, cluster_col, feature_col, num_features=None, missingLevelDefaults=None, sep="."):
        
        #mean_col = "mean.cluster"
        #sd_col = "sd.cluster"
        #expr_col = "expr.cluster"
        #cluster_col = "condition_block"
        #feature_col = "gene"
    
        all_features = set(indf[feature_col])
    
        exprdata = defaultdict(lambda: defaultdict(list))
    
        for row in indf.iter_rows(named=True):
    
            mean_val = row[mean_col]
            sd_val = row[sd_col]
            expr_val = row[expr_col]
    
            cluster_val = row[cluster_col]
            feature_val = row[feature_col]
    
            exprdata[cluster_val][feature_val].append( (mean_val, sd_val, expr_val) )

        fuzzydata = defaultdict(list)
        fuzzydfs = []
    
        for cluster in exprdata:
            
            if (not num_features is None) and (cluster in num_features):
                numSamples = num_features[cluster]
            else:
                #numSamples = len([x for x in exprdata[cluster][feature]])
                numSamples = 1
                
            print("Processing cluster", cluster, "with", numSamples, "samples")
            
            
            levelNames = [x for x in exprMFs[cluster].terms]
            fuzzyHeader = tuple([feature_col] + ["{}{}{}".format(x, sep, cluster) for x in levelNames])
            
            for feature in all_features:

                featureElems = [(0,0,0)]
    
                if feature in exprdata[cluster]:
                    featureElems = [x for x in exprdata[cluster][feature]]
                    while len(featureElems) < numSamples:
                        featureElems.append( (0,0,0) )
            
                res = None

                for elem in featureElems:
                
                    fuzzydat = distribution_to_fuzzy( elem[0], elem[1], elem[2], exprMFs[cluster], threshold=0.0)
                
                    if res is None:
                        res = fuzzydat
    
                    else:
                        res = [x+y for x, y in zip(res, fuzzydat)]
    
                usedSamples = len(featureElems)
                res = [x/usedSamples for x in res]
    
                resElem = tuple( [feature] + list(res) )
                fuzzydata[cluster].append(resElem)

            clusterdf = pd.DataFrame.from_records(fuzzydata[cluster], columns=fuzzyHeader)
            fuzzydfs.append(clusterdf)
            
        from functools import reduce
        fuzzyDF = reduce(lambda x,y: pd.merge(x, y, on=feature_col), fuzzydfs)
        fuzzyPlDF = pl.from_pandas(fuzzyDF)
        
        all_classes = []
        for state in exprMFs:
            for level in exprMFs[state].terms:
                all_classes.append("{}{}{}".format(level, sep, state))
                
        
        for x in all_classes:
            if not x in fuzzyPlDF.columns:
                defaultVal = missingLevelDefaults.get(x, 0.0)
                print("Need to add", x, "with default", defaultVal)
                fuzzyPlDF = fuzzyPlDF.with_columns(pl.lit( defaultVal ).alias(x))
                
        
        return fuzzyPlDF

class LegacyFuzzifier(AbstractFuzzifier):
    
    def __init__(self) -> None:
        super().__init__()
        
    @classmethod
    def fuzzify(self, indf:pl.DataFrame, series = None, perSeriesFuzzy=False, mfLevels = ["NO", "LOW", "med", "HIGH"], mfLevelsMirrored=False, centers=None,symbol_column="gene", meancolName="mean.cluster", sdcolName="sd.cluster", exprcolName="expr.cluster", clusterColName="cluster", shape="tri", stepsize=None,combineOverState=False, fuzzifiers=None, **kwargs):

        exprData = indf.clone()

        if fuzzifiers is None:
            exprMFs = self.make_fuzzy_concepts(exprData, mfLevels, centers, clusterColName, meancolName, mfLevelsMirrored, stepsize=stepsize, shape=shape, 
                                            series=series, perSeriesFuzzy=perSeriesFuzzy, **kwargs)
        else:
            exprMFs = fuzzifiers
        
        
        meanExprCol = exprData.columns.index(meancolName)
        clusterCol = exprData.columns.index(clusterColName)

        if not exprcolName is None:
            exprCountCol = exprData.columns.index(exprcolName)
        else:

            exprData=exprData.with_columns(pl.lit(1).alias("cell_expr"))
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
                    .with_columns(
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
                 
        dfWide = to_homogeneous(
                                    toWideDF(df,symbol_column, cluster_column=clusterColName)
                                    , exprMFs)
        
        
        explDFWide = self.to_vwide(dfWide, exprMFs, meta_columns=[ symbol_column ])

        return explDFWide, exprMFs



class FastFuzzifier(AbstractFuzzifier):
    
    def __init__(self) -> None:
        super().__init__()
        

    @classmethod
    def fuzzify(cls, indf:pl.DataFrame, series = None, perSeriesFuzzy=False, mfLevels = ["NO", "LOW", "med", "HIGH"], mfLevelsMirrored=False, centers=None,symbol_column="gene", meancolName="mean.cluster", sdcolName="sd.cluster", exprcolName="expr.cluster", clusterColName="cluster", shape="tri", stepsize=None,combineOverState=False, fuzzifiers=None, **kwargs):
        import datetime

        a0 = datetime.datetime.now()

        exprData = indf.clone()

        #exprMFs = cls.make_fuzzy_concepts(exprData, mfLevels, centers, clusterColName, meancolName, mfLevelsMirrored, stepsize=stepsize, shape=shape, series=series, perSeriesFuzzy=perSeriesFuzzy, **kwargs)
        
        if fuzzifiers is None:
            exprMFs = cls.make_fuzzy_concepts(exprData, mfLevels, centers, clusterColName, meancolName, mfLevelsMirrored, stepsize=stepsize, shape=shape, series=series, perSeriesFuzzy=perSeriesFuzzy, **kwargs)
        else:
            exprMFs = fuzzifiers
        
        meanExprCol = exprData.columns.index(meancolName)
        clusterCol = exprData.columns.index(clusterColName)

        if not exprcolName is None:
            exprCountCol = exprData.columns.index(exprcolName)
        else:

            exprData=exprData.with_columns(pl.lit(1).alias("cell_expr"))
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
                
        availableClusters =  list(exprData.select(pl.col(clusterColName)).unique().to_series())
        fuzzyOuts = []
        a1 = datetime.datetime.now()
        print("Time for setting up fuzzy concept:"+str(a1-a0))

        for seriesName in availableClusters:


            a = datetime.datetime.now()
            print("Fuzzifing: "+seriesName)

            indf = df.filter(pl.col(clusterColName) == seriesName)
            signal=list(indf.select(pl.col(meancolName)).to_series())
            signal=[0 if v is None else v for v in signal]


            
            if shape=="crisp":
                print("Crisp mode, different assignment")
      
                # weighted mean is calculated with all (not only expressed) cells
                indf = indf.with_columns(
                    (pl.col(meancolName)*pl.col(exprcolName)).alias(meancolName)
                )
                seriesOut=to_fuzzy_fast(signal, exprMFs[seriesName])
            else:
                if not sdcolName is None:    
                # TODO - implement cases with distribution
      
                    seriesOut = indf.select(
                        pl.struct([meancolName, sdcolName, exprcolName]).apply(lambda x:
                            distribution_to_fuzzy(x[meancolName], x[sdcolName], x[exprcolName], exprMFs[seriesName], threshold=0.0)
                            ).alias("fuzzy.mfs")
                    )
                else:                    
                    #print(len(exprMFs[seriesName].terms))
                    # TODO - implement cases with distribution

                    seriesOut=to_fuzzy_fast(signal, exprMFs[seriesName])

            b = datetime.datetime.now()
            print("Time for fuzzification:"+str(b-a))

            seriesOut=seriesOut.with_columns(pl.all().prefix("FV_"))
            FV_columns=list(filter(lambda x:'FV_' in x, seriesOut.columns))
            seriesOut=seriesOut.select(pl.col(FV_columns))
            new_indf = indf.hstack(seriesOut)
            if combineOverState == True:
                

                #print(new_indf)
         

                new_indf=new_indf.groupby([symbol_column,clusterColName], maintain_order=True).agg([
                    pl.col(FV_columns).mean()
                ])

                
                c = datetime.datetime.now()
                print("Time for combining:"+str(c-b))


            # Loose here everything besides fuzzy values
            new_indf=new_indf.select(pl.col([symbol_column,clusterColName]+FV_columns))
            # Loose here multiple given features
            new_indf=new_indf.unique(subset=[symbol_column,clusterColName], maintain_order=True)
            new_indf=(new_indf
                .melt(
                    id_vars = [symbol_column,clusterColName], 
                    variable_name = 'FV',
                    value_name="fuzzy.mfs"
                    )
                )

            new_indf=new_indf.select(pl.col([symbol_column,clusterColName,'fuzzy.mfs'])).groupby([symbol_column,clusterColName], maintain_order=True).agg(pl.col("fuzzy.mfs"))
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
                 
        dfWide = to_homogeneous(toWideDF(df,symbol_column, cluster_column=clusterColName), exprMFs)
        
        return dfWide, exprMFs



class FlowAnalysis:

    @classmethod
    def _filter_weightSequence(cls, weightSequence,cutoff=None):
        """Filters weight sequence for sequences with weight > cutoff.

        Args:
            weightSequence (list): weight sequence
            cutoff (float, optional): cutoff for filter. Defaults to None.

        Returns:
            list: filtered weightSequence
        """
        if not cutoff is None and cutoff > 0:
            return [x for x in weightSequence if x[1][-1] > cutoff]
        else:
            return weightSequence
    

    @classmethod
    def fuzzify_exprvalues(cls, *args, **kwargs):
        import sys
        
        raise NotImplementedError("This functionality has been moved into the LegacyFuzzifier. First create a LegacyFuzzifier lfz = LegacyFuzzifier(), then call lfz.fuzzify(...). This already included the calls to to_homogeneous and toWideDF.") 
        


    """



    @classmethod
    def toFlowsDF(cls, indf):

        explDF = indf.clone()
        binColumns = [x for x in explDF.columns if x.startswith("bin.")]
        binColumnsIndex = [explDF.columns.index(x) for x in binColumns]
        print(binColumns)

        mfColumns = [x.replace("bin.", "mf.") for x in binColumns]
        print(mfColumns)
        
        df=df.with_columns(df.apply( lambda row: tuple( [row[x] for x in binColumnsIndex] )).to_series().alias("group.flow"))
        df=df.with_columns(pl.struct(mfColumns).apply( lambda row: np.prod(row.values())).to_series().alias("mf.flow"))
        
        allgroups = list(set(df.select(pl.col("group.flow"))))
        
        flow2id = {}
        for idx, flow in enumerate(allgroups):
            flow2id[flow] = idx
        
        df=df.with_columns(pl.struct(["group.flow"]).apply( lambda row: allgroups[row["group.flow"]] ).to_series().alias("id.flow"))

        return explDF
    
    """


    def __init__(self, flows, symbol_column, series2name, exprMF, sep="."):
        
        nExprMF = {}
        for x in exprMF:
            idx = x
            if idx.startswith("cluster."):
                idx = idx.replace("cluster.", "")
            nExprMF[idx] = exprMF[x]
            
        for state in series2name:
            for level in nExprMF[state[0]].terms:
                checkCol = "{}{}{}".format(level, sep, state[0])
                if not checkCol in flows.columns:
                    print("Missing column", checkCol)
                    raise ValueError("Missing column {}".format(checkCol))
        
        self.sep = sep
        self.flows = flows
        self.seriesOrder = [x[0] for x in series2name]
        self.series2name = {x[0]: x[1] for x in series2name}
        self.symbol_column = symbol_column
        self.exprMF = nExprMF
        
        self.levelOrder = {}
        for state in self.exprMF:
            self.levelOrder[state] = [x for x in self.exprMF[state].terms]
            
        self._edgeid2flow = None            
            
        notEqualLevelOrder = False
        for state1 in self.exprMF:
            for state2 in self.exprMF:
                if len(set(self.levelOrder[state1]).difference(self.levelOrder[state2])):
                    notEqualLevelOrder = True        
        if notEqualLevelOrder:
            print("WARNING: Your level orders are not equal. Scales in some plots may not reflect this appropriately.")
        


    """
    # TODO - wozu?
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
    """

    def copy(self):
        """Creates a copy of the current FlowAnalysis object (deep: flows, shallow: series2name, exprMF)

        Returns:
            FlowAnalysis: Copy of the current FlowAnalysis object
        """
        s2nvec = [(x, self.series2name[x]) for x in self.series2name]
        fa = FlowAnalysis(flows=self.flows.clone(), symbol_column=self.symbol_column, series2name=s2nvec, exprMF=self.exprMF.copy(), sep=self.sep)
        
        return fa        


    def subset_states(self, keep_states):
        """ Subsets flowset object to conly contain state from the keep_states attribute.

        Args:
            keep_states (list): list of descriptors for states contained in this flowsets object

        Returns:
            FlowsetAnalysis: object only containing the states given by keep_states
        """
        
        fa = self.copy()
        fa.series2name = {x: fa.series2name[x] for x in fa.series2name if x in keep_states}
        fa.seriesOrder = [x for x in fa.series2name]
        
        remainingColumns = [fa.symbol_column] + [x for x in fa.flows.columns if x.split(fa.sep)[-1] in keep_states]
        fa.flows = fa.flows.select (pl.col (remainingColumns))
    
        fa._edgeid2flow = None
        
        return fa

    def filter_genes(self, symbol_prefixes=["AC0", "AL0", "AP0", "AC1", "AL1", "AP1"]):
        """Filters out genes starting with the given prefixes

        Args:
            symbol_prefixes (list, optional): _description_. Defaults to ["AC0", "AL0", "AP0", "AC1", "AL1", "AP1"].
        """
        
        print(self.flows.shape)
        for gp in symbol_prefixes:
            self.flows = self.flows.filter( ~pl.col(self.symbol_column).str.starts_with(gp) )
            
        print(self.flows.shape)

    #
    ##
    ### EDGE BASED
    ##
    #


    def plot_flows(self, use_edges:set = None, genes:list=None,figsize:tuple=None, outfile:str=None, min_flow:float=None,  transformCounts = lambda x: x, verbose=False,specialColors=None,sns_palette="icefire", seriesColors=None, colorMode="scaling",title:str=None,linewidth=0, seriesFontsize=10, classFontsize=12):
        """Plots the FlowSets system as sankey plot.

        Args:
            use_flows (set, optional): Restricts the system to only selected flows. Defaults to None.
            genes (list, optional): Only considers specified genes. Defaults to None.
            figsize (tuple, optional): Tuple (width, height) of figure. Defaults to None.
            outfile (str, optional): Path to file to save plot in. Defaults to None.
            min_flow (float, optional): Remove all edges with less membership. Defaults to None.
            transformCounts (_type_, optional): Function to transform weights with for plotting. Defaults to lambdax:x.
            verbose (bool, optional): Verbose output. Defaults to False.
            specialColors (_type_, optional): _description_. Defaults to None.
            sns_palette (str, optional): Color-palette to use for coloring states. Defaults to "icefire".
            seriesColors (dict, optional): Dictionary with seriesname to color. Defaults to None.
            colorMode (str, optional): Which coloring mode for states. Either scaling or diverging. Defaults to "scaling".
            title (str, optional): Title of the plot. Defaults to None.
        """
        
        
        if genes is None:
            flowDF=self.flows
        else:
            flowDF=self.flows.filter(pl.col(self.symbol_column).is_in(genes) )

        
        if use_edges is None:
            use_edges = [x for x in self.edgeid2flow]
        
            
        _, node_memberships, used_edges = self.calc_coarse_flow_memberships(use_edges=use_edges, genes=genes,backtracking=True)
        weightSequence = self._backtrack_coarse_flows(genes=genes, node_memberships=node_memberships, used_edges=used_edges)

            
        if specialColors is None:
            indices=[w[0] for w in weightSequence ]
            startingnodes=[w[1][0][1] for w in weightSequence ]
            #maxLevels = max([len(self.levelOrder[x]) for x in self.levelOrder])
            levelColors=dict.fromkeys(startingnodes)
            for i,k in enumerate(list(levelColors.keys())):
                levelColors[k] =i 
            maxLevels=len(levelColors.keys())#
            colours=sns.color_palette(sns_palette,maxLevels)
            c=[colours[levelColors[s]] for s in startingnodes]
            specialColors=pd.DataFrame({
                'values':c},
                index=indices).to_dict()['values']
            
            weightSequence = self._filter_weightSequence(weightSequence,cutoff=min_flow)
            
            new_indices=[w[0] for w in weightSequence ]
            specialColors={k: v for k, v in specialColors.items() if k in new_indices}
        else:
            weightSequence = self._filter_weightSequence(weightSequence,cutoff=min_flow)

                
        if not seriesColors is None:
            seriesColorMap = self._create_series_color_map(seriesColors, colorMode)
        else:
            seriesColorMap=None
            
        SankeyPlotter._make_plot(weightSequence, self.series2name, self.levelOrder, self.seriesOrder, specialColors=specialColors, seriesColorMap=seriesColorMap, independentEdges=True,outfile=outfile,title=title,fsize=figsize,linewidth=linewidth, seriesFontsize=seriesFontsize, classFontsize=classFontsize)



    def hist_level_membershipsum(self):
        """
            Plots the membership sum for each state in each level.
        """

        flowDF=self.flows
        colmeans=flowDF.select(pl.col(pl.Float64)).transpose(include_header=True).with_columns(
            pl.fold(0, lambda acc, s: acc + s, pl.all().exclude("column")).alias("horizontal_sum")
        )
                
        order1=[self.levelOrder[state].index(level) for level,state in [(x.split( self.sep )[0],x.split( self.sep )[-1])  for x in colmeans.select(['column']).to_series()]]
        order2=[self.series2name[i] for i in [x.split( self.sep )[-1] for x in colmeans.select(['column']).to_series()]]

        colmeans=colmeans.with_columns(pl.Series(name="order1", values=order1))
        colmeans=colmeans.with_columns(pl.Series(name="order2", values=order2))
        
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
        """Plots for the given genes, a matrix visualizing these genes' membership in each state and level.

        Args:
            genes (list): List of features to look at.
            name (str, optional): Title of the plot. Defaults to "".
            cluster_genes (boolean, optional): Whether the genes should be clustered by similarity in their expression pattern. Defaults to False.
            outfile (str, optional): Path to save plot. Defaults to None.
            limits (tuple, optional): Limits for the membership colormap. Defaults to (0,1).
            annot (bool, optional): Whether to plot memberships into the cells. Defaults to True.
            annot_fmt (str, optional): Str-format for memberships plotted into cells. Defaults to ".2f".
            prefix (str, optional): Prefix for the states. Defaults to "Cluster".
            verbose (bool, optional): Verbose output. Defaults to False.
            figsize (tuple, optional): Size of the created figure. Defaults to (6,6).
            font_scale (float, optional): Scale for the fonts shown in plot. Defaults to 0.4.

        Returns:
            seaborn clustermap: Clustermap object for further manipulation.
        """
        
        
        filtered_flow=self.flows.filter(pl.col(self.symbol_column).is_in(genes) )      
        #pd_filtered_flow=pd.DataFrame(filtered_flow[:,1:], columns=filtered_flow[:,1:].columns, index=filtered_flow[:,0].to_pandas().tolist())
        pd_filtered_flow=filtered_flow.to_pandas().set_index(self.symbol_column)

        
        pd_filtered_flow = pd_filtered_flow.transpose()     
        
        pd_filtered_flow["orderLevel"]=[self.levelOrder[state].index(level) for level,state in [(x.split(self.sep)[0],x.split(self.sep)[1]) for x in pd_filtered_flow.index]]
        pd_filtered_flow["orderState"]=[self.seriesOrder.index(i) for i in [x.split(self.sep)[1] for x in pd_filtered_flow.index]]
        pd_filtered_flow= pd_filtered_flow.sort_values(["orderState", "orderLevel"])
        
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
        for x in states:
            if x == lastState:
                curCount += 1
            else:
                stateCounts.append(curCount)
                lastState = x
                curCount = 1
        stateCounts.append(curCount)

        stateCounts=np.cumsum(stateCounts[::-1])
        pd_filtered_flow = pd_filtered_flow.drop(["orderLevel","orderState"], axis=1)
        
        #for x in pd_filtered_flow.index:
        #    x = x.split(".")
        #    
        #    if prefix is None:
        #        newIndex.append("{} {} ({}) ".format(".".join(x[1:len(x)-1]).title(), x[-1], x[0]))
        #    else:
        #        newIndex.append("{} {} ({}) ".format(prefix, x[-1], x[0]))
            
        pd_filtered_flow.index = newIndex
                
        sns.set(font_scale=font_scale)
        
        pd_filtered_flow = pd_filtered_flow.iloc[::-1]
        
        g=sns.clustermap(pd_filtered_flow, figsize=figsize, row_cluster=False, col_cluster=cluster_genes, vmin=limits[0], vmax=limits[1], annot=annot, fmt=annot_fmt, cbar_pos=None, dendrogram_ratio=0.1, yticklabels=True, xticklabels=True)
        
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





    def highlight_genes(self, genes, figsize=None, outfile=None, min_flow=None, min_gene_flow=None, transformCounts=lambda x: np.sqrt(x)):
        """Visualizes the selected genes in a different color within the system.

        Args:
            genes (list): Genes to visualize.
            figsize (tuple, optional): plot size. Defaults to None.
            outfile (str, optional): path to write out plot. Defaults to None.
            min_flow (float, optional): only plot edges with flow > min_flow. Defaults to None.
            transformCounts (function, optional): function with which the counts are transformed for plotting. Defaults to lambdax:np.sqrt(x).
        """

        if not isinstance(genes, (tuple, list, set)):
            genes = [genes]
            
            
                        
        # full system
        bgData = self.flows.filter( ~pl.col(self.symbol_column).is_in(genes) )    

        print("Background")       
        _, node_memberships, used_edges = self.calc_coarse_flow_memberships(flowDF=bgData, use_edges=None,genes=None,backtracking=True)
        bgWeightSequence = self._backtrack_coarse_flows(flowDF=bgData, genes=None, node_memberships=node_memberships, used_edges=used_edges)
        bgWeightSequence = self._filter_weightSequence(bgWeightSequence, min_flow)

        survivedFlows = [x[0] for x in bgWeightSequence]
        
        fgData = self.flows.filter( pl.col(self.symbol_column).is_in(genes) )
       
        print("Foreground")
        _, node_memberships, used_edges = self.calc_coarse_flow_memberships(flowDF=fgData, use_edges=None,genes=None,backtracking=True)
        fgWeightSequence = self._backtrack_coarse_flows(flowDF=fgData, genes=None, node_memberships=node_memberships, used_edges=used_edges, edgeIDMod=lambda x: x*(-1))
        fgWeightSequence = self._filter_weightSequence(fgWeightSequence, min_flow)


        specialColors = {x[0]: "red" for x in fgWeightSequence}
        
        SankeyPlotter._make_plot(bgWeightSequence+fgWeightSequence, self.series2name, self.levelOrder, self.seriesOrder, specialColors=specialColors, transformCounts=transformCounts, fsize=figsize, outfile=outfile, independentEdges=True)


    def visualize_genes(self, genes, figsize=None, min_flow=None, use_edges=None, title=None, outfile=None, score_modifier=lambda x: x, colormap="cividis", seriesColors=None, colorMode="scaling"):
        """Plots the flow system and colors edges by the genes' memberships.

        Args:
            genes (list): list of genes
            figsize (tuple, optional): Size of the plot. Defaults to None.
            min_flow (float, optional): only plot edges with flow > min_flow. Defaults to None.
            use_edges (list, optional): restrict edges of the system. Defaults to None.
            title (str, optional): title of the plot. Defaults to None.
            outfile (str, optional): path to save plot to. Defaults to None.
            score_modifier (function, optional): Function to modify score with. Defaults to lambdax:x.
        """

        if not isinstance(genes, (tuple, list, set)):
            genes = [genes]
            
        if use_edges is None:
            use_edges = [x for x in self.edgeid2flow]
            
            
        #bgData = self.flows.filter( ~pl.col("gene").is_in(genes) )
        _, node_memberships, used_edges = self.calc_coarse_flow_memberships(flowDF=self.flows, use_edges=use_edges,genes=None,backtracking=True)
        bgWeightSequence = self._backtrack_coarse_flows(flowDF=self.flows, genes=None, node_memberships=node_memberships, used_edges=used_edges)
        bgWeightSequence = self._filter_weightSequence(bgWeightSequence, min_flow)

        
        fgData = self.flows.filter( pl.col(self.symbol_column).is_in(genes) )

        _, node_memberships, used_edges = self.calc_coarse_flow_memberships(flowDF=fgData, use_edges=use_edges,genes=None,backtracking=True)
        fgWeightSequence = self._backtrack_coarse_flows(flowDF=fgData, genes=None, node_memberships=node_memberships, used_edges=used_edges)
        fgWeightSequence = self._filter_weightSequence(fgWeightSequence, min_flow)

        for x in fgWeightSequence:
            x[1][-1] = score_modifier(x[1][-1])
            
        
        #print(fgWeightSequence)
        maxFlowValue = max([ x[1][-1] for x in fgWeightSequence])
        #print(maxFlowValue)

        cmap = get_cmap(colormap)
        norm = mpl.colors.Normalize(vmin=0, vmax=maxFlowValue)
        
        seriesColorMap=None
        if not seriesColors is None:
            seriesColorMap = self._create_series_color_map(seriesColors, colorMode)
        
        specialColors = {x[0]: cmap(x[1][-1]/maxFlowValue) for x in fgWeightSequence}
        SankeyPlotter._make_plot(bgWeightSequence, self.series2name, self.levelOrder, self.seriesOrder,
                                 specialColors=specialColors, transformCounts=lambda x: np.sqrt(x),
                                 fsize=figsize, cmap=cmap, norm=norm, outfile=outfile,
                                 title=title, independentEdges=True, seriesColorMap=seriesColorMap)
                       

    def make_plot_flow_memberships(self, flowScores_df, n_genes=30, color_genes=None, figsize=(2,5), outfile=None, plot_histogram=True,violin=False, labelsize=4, countsize=8,draw_zscores=True):
               
        flowScores_df=flowScores_df.sort(["membership",self.symbol_column],descending=[True, False])
        
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
        ax1.grid(which='major', color='#EBEBEB', linewidth=1)
        
        
        ax1.hlines(y=range(n_genes), xmin = 0 , xmax = flowScores_df_top["membership"], color=colormap)
        ax1.plot(list(flowScores_df_top["membership"]), range(n_genes), "o")
        ax1.set_yticks(range(n_genes),flowScores_df_top[self.symbol_column])
        ax1.tick_params(axis="y",labelsize=labelsize)
        ax1.set_title("Top "+str(n_genes)+" memberships",fontsize = labelsize)
        ax1.invert_yaxis()

        # Only show minor gridlines once in between major gridlines.
        ax1.xaxis.set_minor_locator(AutoMinorLocator(2))


        if draw_zscores:
            mean= flowScores_df.filter(pl.col("membership")!= 0).select(pl.col("membership")).to_series().mean()
            std= flowScores_df.filter(pl.col("membership")!= 0).select(pl.col("membership")).to_series().std()
            
            if not mean is None:
                ax1.axvline(x=mean, color='r',linestyle="solid")
                ax1.axvline(x=mean+std, color='r',linestyle="dashed")
                ax1.axvline(x=mean+10*std, color='r',linestyle="dashdot")
                ax1.axvline(x=mean+50*std, color='r',linestyle="dotted")

        if not ax2 is None:
            flowScores_df_rounded=flowScores_df.with_columns(
                pl.col("membership").round(1)
            )
            counted_scores=flowScores_df_rounded.select("membership").to_series().value_counts()

            if violin:
                sns.set_style('whitegrid')
                #sns.kdeplot(x=np.array(flowScores_df["pwscore"]),ax=ax2)
                sns.violinplot(x=np.array(flowScores_df["membership"]),ax=ax2,inner="stick")
                ax2.set_title("Flow membership distribution")
                ax2.set_xlim(0,1)
            else:
                bars=ax2.bar(counted_scores['membership'], counted_scores['counts'], width = 0.1)
                for rect in bars:
                    real_height=rect.get_height()
                    ax2.text(rect.get_x() + rect.get_width()/2., (real_height*.01).clip(min=2.5),
                            '%d' % int(real_height),c='navy',
                            ha='center', va='bottom',rotation=90,fontsize=countsize,weight='bold')
                ax2.set_title("Binned membership histogram",fontsize = labelsize)
                ax2.set_yscale('log')
                ax2.set_facecolor('white')
                ax2.grid(which='major', color='#EBEBEB', linewidth=.1)
                ax2.grid(which='minor', color='#EBEBEB', linewidth=.1)
                ax2.tick_params(axis="y",labelsize=labelsize)
                ax2.axis(ymin=1)
        fig.tight_layout()
        plt.xlim([-0.05, 1.05])

        if not outfile is None:
            plt.savefig(outfile + ".png", bbox_inches='tight')
            plt.savefig(outfile + ".pdf", bbox_inches='tight')
        
        plt.show()
        
        if ax2 is None:
            return flowScores_df_top, ax1
        
        return flowScores_df_top, ax1, ax2
        
        #return list(flowScores_df_top.select(pl.col(self.symbol_column)))[0].to_list(), flowScores_df, (ax1, ax2)


    



    def plot_flow_memberships(self,use_edges, genes=None,n_genes=30,color_genes=None, figsize=(2,5), outfile=None, plot_histogram=True, violin=False, labelsize=4, countsize=8, gene_exclude_patterns=[]):

        flowDF=self.flows.clone()
        
        for gene_exclude_pattern in gene_exclude_patterns:
            flowDF = flowDF.filter(~pl.col(self.symbol_column).str.starts_with(gene_exclude_pattern))

        flowScores_df = self.calc_coarse_flow_memberships(flowDF=flowDF, use_edges=use_edges,genes=genes,backtracking=False)
        
        return self.make_plot_flow_memberships(flowScores_df, n_genes=n_genes,color_genes=color_genes, figsize=figsize, outfile=outfile, plot_histogram=plot_histogram,violin=violin, labelsize=labelsize, countsize=countsize)


    def calc_coarse_flow_memberships(self, flowDF=None, use_edges=None, genes=None, backtracking=False):
        """Calculates memberships of genes given flows as edgeIDs or genes.

        Args:
            use_edges (set, optional): Set of edgeIDs to consider. If None, all edges are considered. Defaults to None.
            genes (list, optional): List of genes to consider. If None, all genes are considered. Defaults to None.
            backtracking (bool, optional): If true, all data required for backtrack are returned. Defaults to False.

        Returns:
            polars.DataFrame: DataFrame which contains for each feature a membership
        """
        

        if flowDF is None:
            flowDF = self.flows
            
        if not genes is None:        
            flowDF = self.flows.filter( pl.col(self.symbol_column).is_in(genes) )
    

        allSeries = [x for x in self.series2name]

        node_memberships={}
        weightSequence_before = []
        
        for edgeID in self.edgeid2flow:
            
            largeComp = self.edgeid2flow[edgeID]
            
            src, tgt = largeComp
            
            srcSeries = src[0]
            tgtSeries = tgt[0]
            
            srcLevel = src[1]
            tgtLevel = tgt[1]
            
            i = allSeries.index(srcSeries)
                   
            srcSeries = allSeries[i]
            tgtSeries = allSeries[i+1]
            
            # TODO this could be done too often
            if not srcSeries in node_memberships:
                node_memberships[srcSeries]={}
                for l in self.levelOrder[srcSeries]:
                    colname = "{}{}{}".format(l, self.sep, srcSeries)
                    node_memberships[srcSeries][l]=flowDF.select(pl.col(colname)).to_series()


            # TODO this could be done too often
            if not tgtSeries in node_memberships: 
                node_memberships[tgtSeries]={}
                for l in self.levelOrder[tgtSeries]:
                    node_memberships[tgtSeries][l]=pl.Series([0]*flowDF.shape[0])


            if not use_edges is None:
                if not edgeID in use_edges:
                    continue

            flowCols=["{}{}{}".format(fclass, self.sep, state) for state, fclass  in largeComp ]
            
            
            temp_df=pl.DataFrame({
                                    "start":node_memberships[srcSeries][ srcLevel ],
                                    "target":node_memberships[tgtSeries][ tgtLevel ],
                                    "factor":flowDF.select(pl.col(flowCols[1])).to_series()
                                })

            temp_df=temp_df.with_columns([
                pl.map(["start", "target", "factor"], lambda s: s[0] * s[2] ).alias("flow"),
                pl.map(["start", "target", "factor"], lambda s: s[1] + s[0] * s[2] ).alias("pwscore")
            ])         
                          
            node_memberships[tgtSeries][ tgtLevel ]=temp_df['pwscore']


            outlist = list(largeComp)
            outlist.append(temp_df["flow"].sum() )
            weightSequence_before.append( (edgeID, outlist) )

        end_memberships=node_memberships[allSeries[-1] ]
    
        pattern_membership=sum(end_memberships.values())
        flowScores_df=pl.DataFrame({
            self.symbol_column : flowDF[self.symbol_column],
            'membership':pattern_membership
            })
        used_edges=[x[0] for x in weightSequence_before if x[1][2] >0]

        if backtracking:
            return flowScores_df, node_memberships, used_edges

        return flowScores_df
                        
    def _backtrack_coarse_flows(self, node_memberships, used_edges, flowDF=None, genes=None, edgeIDMod=None):
        
        allSeries = [x for x in self.series2name]
        
        if flowDF is None:
            flowDF = self.flows
            
        if not genes is None:        
            flowDF = self.flows.filter( pl.col(self.symbol_column).is_in(genes) )
    

        numGenes = flowDF.shape[0]
        end_memberships=node_memberships[allSeries[-1] ]
        
        
        backtracking_memberships={}
        weightSequence_after = []

        selected_edges=[self.edgeid2flow[x] for x in used_edges ]
        start_nodes=[x[0] for x  in selected_edges]
        target_nodes=[x[1] for x  in selected_edges]
        for i in reversed(range(len(allSeries) - 1)):

            srcSeries = allSeries[i]
            tgtSeries = allSeries[i+1]


            if not srcSeries in backtracking_memberships:
                backtracking_memberships[srcSeries]={}
                for l in self.levelOrder[srcSeries]:
                    backtracking_memberships[srcSeries][l]=pl.Series([0]* numGenes)


            if not tgtSeries in backtracking_memberships: 
                backtracking_memberships[tgtSeries]={}
                for l in self.levelOrder[tgtSeries]:
                    backtracking_memberships[tgtSeries][l]=end_memberships[l]

            for comb1 in reversed(self.levelOrder[tgtSeries]):

                edgeIDs= [list(used_edges)[x] for x in range(len(selected_edges))  if (tgtSeries,comb1) == target_nodes[x] ] 
                edges=[self.edgeid2flow[x] for x in edgeIDs]
                nodes=[x[0] for x  in edges]

                #print(edges)
                if len(edges)>0:
                    flowCols=["{}{}{}".format(fclass, self.sep, state) for state, fclass  in nodes ]


                    relative_flow=flowDF.select(pl.col(flowCols))

                    relative_flow=relative_flow.with_columns(
                        pl.fold(0, lambda acc, s: acc + s,pl.all()).alias("horizontal_sum")
                            )
                    
                    relative_flow=relative_flow.with_columns(
                        pl.all().exclude("horizontal_sum") / pl.col(("horizontal_sum"))
                    )

                    relative_flow=relative_flow.fill_nan(0)  


                    for edgeID in edgeIDs:
                        node=self.edgeid2flow[edgeID][0]
                        flowCol= "{}{}{}".format(node[1], self.sep, node[0]) #  node[1]+".cluster."+node[0]
                        
                        temp_df=pl.DataFrame({
                                              "start":backtracking_memberships[srcSeries][node[1]],
                                              "target":backtracking_memberships[tgtSeries][comb1],
                                              "factor":relative_flow.select(pl.col(flowCol)).to_series()
                                              })

                        #print(temp_df)
                        temp_df=temp_df.with_columns([
                            pl.map(["start", "target", "factor"], lambda s: s[1] * s[2] ).alias("flow"),
                            pl.map(["start", "target", "factor"], lambda s:s[0] + s[1] * s[2] ).alias("membership")
                        ])           
                        backtracking_memberships[srcSeries][node[1]]=temp_df['membership']

                        outlist = list(self.edgeid2flow[edgeID])
                        outlist.append(temp_df["flow"].sum() )
                        
                        if not edgeIDMod is None:
                            edgeID = edgeIDMod(edgeID)
                        
                        weightSequence_after.append(
                                        (edgeID, outlist)
                                    )
            
        return weightSequence_after               

            #SankeyPlotter._make_plot(weightSequence_after, self.series2name, self.levelOrder, self.seriesOrder,independentEdges=True)






    # TODO wozu
    def get_confusion_matrix(self,scores_df,Trues,num_true=None,outfile=None):
        scores_df=scores_df.sort(["membership"],reverse=True)

        if not num_true is None:
            scores_df=scores_df[:num_true]
        scores_df=scores_df.filter(pl.col("membership")>0.0)

        Results=list(scores_df.select(pl.col(self.symbol_column)))[0].to_list()
        All=list(self.flows.select(pl.col(self.symbol_column)))[0].to_list()
        cm=confusion_matrix(Trues,Results,All,outfile)

        print(cm)


    
    @property
    def edgeid2flow(self):
        
        if self._edgeid2flow is None:
            self._edgeid2flow = {}
            allSeries = [x for x in self.series2name]

            for i in range(len(allSeries) - 1):

                srcSeries = allSeries[i]
                tgtSeries = allSeries[i+1]

                for comb in list(itertools.product(
                    reversed(self.levelOrder[srcSeries]),
                    reversed(self.levelOrder[tgtSeries]) )):
                            
                            largeComp = [x for x in zip([allSeries[i],allSeries[i+1]], comb)]        
                            fgid=len(self._edgeid2flow)
                            self._edgeid2flow[fgid] = largeComp
                            
        return self._edgeid2flow
    

    def flow_finder( self, pattern:list, minLevels:list=None, maxLevels:list=None, verbose=True ):
        """Select flow according to a specific pattern with minimal and maximal levels per state.

        Args:
            pattern (list): pattern between states
            minLevels (list, optional): Minimal levels per state. None for no restriction. Defaults to None.
            maxLevels (list, optional): Maximal levels per state. None for no restriction. Defaults to None.
            verbose (bool, optional): Verbose output. Defaults to True.

        Returns:
            set: set of matching edge ids
        """
        

        allSeries = [x for x in self.series2name]
        matchingFIDs = set()

        if minLevels is None:
            minLevels=[None]*len(allSeries)
        if maxLevels is None:
            maxLevels=[None]*len(allSeries)
            
            
        for edgeID in self.edgeid2flow:
            
            src, tgt = self.edgeid2flow[edgeID]
            
            srcSeries = src[0]
            tgtSeries = tgt[0]
            
            srcLevel = src[1]
            tgtLevel = tgt[1]
            
            i = allSeries.index(srcSeries)
            comp = pattern[ i ]

            startIdx = self.levelOrder[srcSeries].index(srcLevel)
            endIdx = self.levelOrder[tgtSeries].index(tgtLevel)

            if not minLevels[i] is None:

                if startIdx<self.levelOrder[srcSeries].index(minLevels[i]):
                    continue
            if not minLevels[i+1] is None:
                if endIdx<self.levelOrder[tgtSeries].index(minLevels[i+1]):
                    continue
            if not maxLevels[i] is None:
                if startIdx>self.levelOrder[srcSeries].index(maxLevels[i]):
                    continue
            if not maxLevels[i+1] is None:
                if endIdx>self.levelOrder[tgtSeries].index(maxLevels[i+1]):
                    continue

            if comp == "<":
                if not startIdx < endIdx:
                    continue
            elif comp == "<<":
                if not startIdx+1 < endIdx:
                    continue
            elif comp == "<=":
                if not startIdx <= endIdx:
                    continue
            elif comp == ">=":
                if not startIdx >= endIdx:

                    continue
            elif comp == ">":
                if not startIdx > endIdx:
                    continue
            elif comp == ">>":
                if not startIdx > endIdx+1:
                    continue
            elif comp == "=":
                if not startIdx == endIdx:
                    continue
            elif comp == "x":
                #change
                if startIdx == endIdx:
                    continue
            elif comp == "~":
                #circa
                if not (startIdx == endIdx or startIdx == endIdx+1 or startIdx == endIdx -1) :
                    continue
            elif comp == "?":
                #always accept
                pass


            if verbose:
                print(edgeID, (src, tgt))
                
            matchingFIDs.add(edgeID)

        return matchingFIDs









    #
    ##
    ### PATH BASED
    ##
    #

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


    def path_finder( self, sequence, minLevels=None, maxLevels=None, verbose=True ):

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


    def _paths_to_weight_sequence(self, flows, use_flows=None, flowIDMod=None, min_gene_flow=None):
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
            flowScore, flow = self._calculate_flow_score_for_paths(flows, fgid, min_gene_flow=min_gene_flow)
            
            if not flowIDMod is None:
                fgid = flowIDMod(fgid)
            #print(fgid, flow, flowScore)
            
            outlist = list(flow)
            outlist.append(flowScore )
            weightSequence.append(
                (fgid, outlist)
            )
                                    
        return weightSequence

    def plot_genes_membership(self,genes, n_genes=30,use_flows=None,plot_hist=False, outfile=None):
        genes_flow=self.flows.filter(pl.col(self.symbol_column).is_in(genes) )
        
        
        if use_flows is None:
            use_flows = [x for x in self.flowid2flow]

        weightSequence = self._paths_to_weight_sequence( flows=genes_flow, use_flows=use_flows)

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

        flowScores_df=flowScores_df.with_columns(
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

    def plot_paths(self, use_flows = None, figsize=None, outfile=None, min_flow=None, min_gene_flow=None, transformCounts = lambda x: np.sqrt(x), verbose=False, seriesColors=None, colorMode="scaling"):

        if use_flows is None:
            use_flows = [x for x in self.flowid2flow]

        weightSequence = self._paths_to_weight_sequence( flows=self.flows, use_flows=use_flows)
        weightSequence = self._filter_weightSequence(weightSequence, min_flow)          
                    
        if verbose:
            for x in weightSequence:
                print(x)

        if not seriesColors is None:
            seriesColorMap = self._create_series_color_map(seriesColors, colorMode)
        else:
            seriesColorMap=None
            
        SankeyPlotter._make_plot(weightSequence, self.series2name, self.levelOrder, self.seriesOrder, transformCounts=transformCounts, fsize=figsize, outfile=outfile, verbose=verbose, seriesColorMap=seriesColorMap)


    def _calculate_flow_score_for_paths(self, flowDF, flowID, min_gene_flow=None):
        
        
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



    def _get_flow_columns(self, flowID):
        
        flow = self.flowid2flow[flowID]
        
        flowCols = ["{}{}{}".format(y, self.sep,x) for x,y in flow]
        
        return flowCols, flow
    
    
    def paths2edges(self, paths):
        """Converts a list of flowIDs or flows into an edgeID list

        Args:
            paths (list): list of flowIDs or flows

        Returns:
            list: edgeIDs
        """
        
        paths = list(paths)
        requiredEdges = set()
                
        if type(paths[0]) == int:
            # integer sequence => flowid2flow
            
            for pathid in paths:
                flow = self.flowid2flow[pathid]
                for i in range(0, len(self.series2name)-1):
                    edge = (flow[i], flow[i+1])
                    requiredEdges.add(edge)
                    
        else:
            
            for flow in paths:
                for i in range(0, len(self.series2name)-1):
                    edge = (flow[i], flow[i+1])
                    requiredEdges.add(edge)
                    
        edgeIDs = set()
        
        for edgeID in self.edgeid2flow:
            if tuple(self.edgeid2flow[edgeID]) in requiredEdges:
                edgeIDs.add(edgeID)
                
        return edgeIDs
    
    def edges2paths(self, edgeIDs):
        """Converts a list of edge IDs into paths

        Args:
            edgeIDs (list): edge IDs

        Returns:
            list: paths
        """

        series2edge = defaultdict(set)

        for edgeID in edgeIDs:

            edge = tuple(self.edgeid2flow[edgeID])
            src, tgt = edge

            srcSeries, srcLevel = src
            tgtSeries, tgtLevel = tgt

            series2edge[srcSeries].add(edge)

        currentPaths = series2edge[ self.seriesOrder[0] ]
        print(currentPaths)

        for series in self.seriesOrder[2:]:

            nextCurrentPaths = []
            for path in currentPaths:

                lastNode =  tuple(path[-1])
                
                for edge in series2edge[lastNode[0]]:

                    if edge[0] == (lastNode[0], lastNode[1]):
                        newpath = list(path) + [edge[1]]
                        nextCurrentPaths.append(newpath)

            currentPaths = nextCurrentPaths

        currentPaths = set( [tuple(x) for x in currentPaths] )
        return currentPaths
            
    def addpaths(self, paths):
        """Adds a list of paths as _flowid2flow (in case calculation of all paths is taking too long)

        Args:
            paths (list): list of paths
        """
        
        self._flowid2flow = {}
        for x in paths:
            self._flowid2flow[len(self._flowid2flow)] = x
            
                    
                        

    #
    ##
    ### ORA analysis
    ##
    #
    




    def get_pathways(self, pathways_file:str,to_upper=True):
        """Reads in gmt or gaf file with genesets

        Args:
            pathways_file (str): path to gmt/gaf-file with genesets

        Raises:
            ValueError: Invalid file format (if neither gmt nor gaf format)

        Returns:
            dict: Dictionary with term-id to tuple (term-name, geneset)
        """
        
        print("Loading pathways from", pathways_file)
        
        if pathways_file.endswith("gmt"):
            rp = self._read_gmt_file(pathways_file,to_upper=to_upper)
        elif pathways_file.endswith("gaf"):
            rp = self._read_gaf_file(pathways_file,to_upper=to_upper)
        else:
            raise ValueError("Invalid File Format")
        
        print("Identified", len(rp), "pathways")
        
        return rp




    def analyse_pathways(self, use_edges:None, genesets_file="ReactomePathways.gmt", additional_genesets=None, set_size_threshold=[ 1,2,3,4, 10, 50, 100], minSetSize=1, feature_modificator=None,to_upper=True):
        """Calculated geneset over-representation results for genesets provided in genesets_file and additional_genesets.

        Args:
            use_flows (set, optional): Set of edges/flows to consider for analysis (if None, all edges are considered)
            genesets_file (str, optional): Path to geneset file in gmt-format. Defaults to "ReactomePathways.gmt".
            additional_genesets (dict, optional): Dictionary (geneset name to geneset) of additional genesets. Defaults to None.
            set_size_threshold (list, optional): Borders of the bins for calculating p-values from z-scores. Defaults to [ 1,2,3,4, 10, 50, 100].

        Returns:
            _type_: _description_
        """

        pathways = self.get_pathways(genesets_file,to_upper=to_upper)
            
        if not additional_genesets is None:
            for pname, pgenes in additional_genesets:
                pathways[pname] = (pname, pgenes)

        # filter pathways
        pathways = {x: pathways[x] for x in pathways if len(pathways[x][1]) >= minSetSize}
        
        if len(pathways) == 0:
            return None


        flowMembership=self.calc_coarse_flow_memberships(use_edges=use_edges)
        
        if not feature_modificator is None:
            
            flowMembership=flowMembership.with_columns(pl.Series([feature_modificator(x) for x in flowMembership[self.symbol_column]]).alias(self.symbol_column))
            print(flowMembership.head())
            

        allPathwayGenes = set()
        for x in pathways:
            pwname, pwGenes = pathways[x]
            for gene in pwGenes:
                allPathwayGenes.add(gene)

        flowGenes = set(flowMembership.select(self.symbol_column).to_series())       
        systemMemberships=self.calc_coarse_flow_memberships(use_edges=None)
        
        if not feature_modificator is None:
            
            systemMemberships=systemMemberships.with_columns(pl.Series([feature_modificator(x) for x in systemMemberships[self.symbol_column]]).alias(self.symbol_column))

        allPWGeneMemberships = systemMemberships.filter(pl.col(self.symbol_column).is_in(list(allPathwayGenes)))["membership"].sum()
        totalFlowMembership=flowMembership["membership"].sum()
        outData = defaultdict(list)

        for pwID in pathways:

            pwName, pwGenes = pathways[pwID]
            pwGenes = list(pwGenes)
            inflow_inset = list(flowGenes.intersection(pwGenes))

            pathwayMembership = flowMembership.filter(pl.col(self.symbol_column).is_in(pwGenes))["membership"].sum()
            flowInPathwayScore = systemMemberships.filter(pl.col(self.symbol_column).is_in(pwGenes))["membership"].sum()
            
            genes_coverage = pathwayMembership / totalFlowMembership if totalFlowMembership > 0 else 0
            pathway_coverage = pathwayMembership / len(pwGenes) if len(pwGenes) > 0 else 0

            # population: all genes
            # condition: genes
            # subset: pathway

            outData["pwid"].append(pwID)
            outData["pwname"].append(pwName)
            outData["pwFlow"].append(pathwayMembership)
            outData["pwGenes"].append(len(pwGenes))
            outData["allPwFlow"].append(totalFlowMembership)
            outData["allPwGenes"].append(allPWGeneMemberships)

            outData["pw_gene_intersection"].append(len(inflow_inset))
            outData["pw_coverage"].append(pathway_coverage)
            outData["genes_coverage"].append(genes_coverage)

            outData["mean_coverage"].append(pathway_coverage*genes_coverage)

        outdf = pd.DataFrame.from_dict(outData)       
        allFGDFs = self._calculate_pvalues(outdf, set_size_threshold=set_size_threshold)
                
        return allFGDFs


    def plotORAresult( self, dfin, title, numResults=10, figsize=(10,10), outfile=None, sep=" ", entryformat="{x}{sep}({y}, pw_cov={z:.3f}/{s})", qvalueColumn="pw_coverage_adj_pval"):
        """ Plots dot-plot for overrepresentation analysis results

        Args:
            dfin (pandas.DataFrame): pandas dataframe containing ORA results
            title (str): title of plot
            numResults (int, optional): Number of results shown in plot. Defaults to 10.
            figsize (tuple, optional): Size of the figure. Defaults to (10,10).
            outfile (str, optional): path to store plot in. Defaults to None.
            sep (str, optional): separator between term name and term id. Defaults to " ".
            entryformat (str, optional): Format string for y-axis description of ORA results. Defaults to "{x}{sep}({y}, pw_cov={z:.3f}/{s})".
        """
        
        df_raw = dfin.copy()

        
        df_raw = df_raw[~np.isnan(df_raw[qvalueColumn])]
        
        if df_raw.shape[0] == 0:
            print("No input given", file=sys.stderr)

        
        def makeTitle(colDescr, colID, colSize, setSize):
            out = []
            for x,y,z, s in zip(colDescr, colID, colSize, setSize):
                out.append(entryformat.format(x=x, sep=sep, y=y, z=z, s=s))

            return out


        # Prepare Data
        #determine plot type

        termIDColumn = "pwid"
        termNameColumn = "pwname"
        #qvalueColumn = "adj_pval"
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
            rvb = self._custom_div_cmap(150, mincol=color2, maxcol=color3, midcol=None)
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
        
        maxNLog = max(-np.log10(df[qvalueColumn]))
        maxLine = ((maxNLog// 10)+1)*10       
        
        # Draw plot
        fig, ax = plt.subplots(figsize=figsize, dpi= 80)
        ax.hlines(y=df.termtitle, xmin=0, xmax=maxLine, color='gray', alpha=0.7, linewidth=1, linestyles='dashdot')
        ax.vlines(x=-np.log(0.05), ymin=0, ymax=numResults, color='red', alpha=0.7, linewidth=1, linestyles='dashdot')
        
        sizeFactor = 10    
        scatter = ax.scatter(y=df.termtitle, x=-np.log10(df[qvalueColumn]), s=df.pwGenes*sizeFactor, c=colorValues, alpha=0.7, )

        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6, func=lambda x: x/sizeFactor)
        labels = [x for x in labels]

        # Title, Label, Ticks and Ylim
        ax.set_title(title, fontdict={'size':12})
        ax.set_xlabel('Neg. Log10 Adj. p-Value', fontdict={'size':12})
        ax.set_yticks(df.termtitle)
        ax.set_yticklabels(df.termtitle, fontdict={'horizontalalignment': 'right'})
        
        plt.tick_params(axis='x', which="major", length=7, width=2)
        plt.tick_params(axis='x', which="minor", length=5, width=2)
        
        if (0-np.log10(df[qvalueColumn].min())) > 50:
            ax.set_xscale('log')
            
            
        ax.xaxis.set_major_formatter(ScalarFormatter())
        
        plt.grid(visible=False)
        plt.tight_layout()
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=8)
        plt.setp(ax.get_xmajorticklabels(), fontsize=8)
        plt.setp(ax.get_xminorticklabels(), fontsize=8)

        if not outfile is None:
            plt.savefig(outfile + ".png", bbox_inches='tight')
            plt.savefig(outfile + ".pdf", bbox_inches='tight')
        
        plt.show()








    #
    ##
    ### INTERNAL FUNCTIONS
    ##
    #
    
    def _create_series_color_map(self, seriesColors, mode="scaling", colormap=None):
        """Creates a colormap for given colors per entry of seriesColors

        Args:
            seriesColors (dict): Dictionary with series name to color
            mode (str, optional): Whether diverging or scaling ( continuous ) colormap should be generated. Defaults to     "scaling".

        Returns:
            dictionary(colormap): For each given seriesname a colormap is generated
        """
        if seriesColors is None:
            seriesColors = {}
            
            if colormap is None:
                colormap = plt.get_cmap("viridis")
            
            for si, x in enumerate(self.series2name):
                seriesColors[self.series2name[x]] = colormap(si/len(self.series2name))
        seriesColorMap = {}
        for x in seriesColors:
            scmap = SankeyPlotter.createColorMap(seriesColors[x], mode)
            seriesColorMap[x] = scmap
            
        return seriesColorMap
    
    def _custom_div_cmap(self, numcolors=11, name='custom_div_cmap',
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
    
    def _calculate_pvalues(self, df:pd.DataFrame, column:str="pw_coverage", set_size_threshold=None):
        """Adds p-values and adj. p-values to ORA data frame

        Args:
            df (pd.DataFrame): ORA data frame
            column (str): column for which to calculate zscores, p-values and adj. p-values
            set_size_threshold (int, optional): Borders of the bins for calculating p-values from z-scores. Defaults to None.

        Returns:
            pd.DataFrame: ORA data frame with p-values (columns: "{column}_adj_pval" and "{column}_pval")
        """
        
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
        
            pwC_mean = curDF[curDF[column] != 0][column].mean()
            pwC_std = curDF[curDF[column] != 0][column].std(ddof=0)
            
            curDF["{}_zscore".format(column)] = (curDF[column]-pwC_mean)/pwC_std
            curDF["{}_pval".format(column)] = norm.sf(abs(curDF["{}_zscore".format(column)]))

            curDF.loc[curDF[ "{}_zscore".format(column) ] < 0, "{}_pval".format(column)] = 1.0
            lastThreshold = t
            
            resultDFs.append(curDF)

        outDF = pd.concat(resultDFs, axis=0)

        _ , elemAdjPvals, _, _ = multipletests(outDF["{}_pval".format(column)], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        outDF["{}_adj_pval".format(column)] = elemAdjPvals
        
        return outDF
    
    def _read_gmt_file(self, filepath:str,to_upper=True):
        """Read in gmt geneset file

        Args:
            filepath (str): file to read in

        Returns:
            dict: dictionary geneset id -> (geneset name, genes)
        """

        geneset2genes = {}

        with open(filepath) as fin:

            for line in fin:

                line = line.strip().split("\t")
                pwName = line[0]
                pwID = line[1]
                pwGenes = line[2:]
                if(to_upper):
                    pwGenes = [x.upper() for x in pwGenes]

                geneset2genes[pwID] = (pwName, pwGenes)

        return geneset2genes
    
    def _read_gaf_file(self, filepath:str,to_upper=True):
        """Read in gaf geneset file

        Args:
            filepath (str): file to read in

        Returns:
            dict: dictionary geneset id -> (geneset name, genes)
        """
        
        geneset2genes = defaultdict(set)
        geneset2name = {}

        with open(filepath) as fin:

            for line in fin:
                
                if line.startswith("!"):
                    continue

                line = line.strip().split("\t")
                pwName = line[4]
                pwID = line[4]
                if(to_upper):
                    pwGene = line[2].upper()
                else:
                    pwGene = line[2]
                geneset2name[pwID] = pwName
                geneset2genes[pwID].add(pwGene)
                
        output = {}
        for x in geneset2name:
            output[x] = (x, geneset2genes[x])

        return output
    
    

    
    
    
    
    
    
    
    #
    ##
    ### OLD FUNCTIONS
    ##
    #
    
    # __deprecated_
    
    def __deprecated_analyse_genes_for_genesets(self, pathways, flowDF, bgFlowDF, considerFlows=None, populationSize=None):
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


    def __deprecated_analyse_pathways_grouped(self, use_flows, pathways_file="ReactomePathways.gmt", additional_pathways=None, set_size_threshold=[ 1,2,3,4, 10, 50, 100]):

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



    
    def __deprecated_analyse_pathways(self, pathways_file="ReactomePathways.gmt", additional_pathways=None, use_flows=None, parallel=True):

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
    

    







    def __deprecated_plot_flow_memberships(self,use_flows, genes=None,n_genes=30,color_genes=None, figsize=(2,5), outfile=None, plot_histogram=True,violin=False, labelsize=4,min_gene_flow=0.0001,parallel=True):
        flowScores_df=self.__deprecated_calc_flow_memberships(use_flows,genes=genes,parallel=parallel)
        return self.make_plot_flow_memberships(flowScores_df,n_genes=n_genes,color_genes=color_genes, figsize=figsize, outfile=outfile, plot_histogram=plot_histogram,violin=violin, labelsize=labelsize)


    
    def __deprecated_calc_flow_memberships(self,use_flows,genes=None, n_genes=30,color_genes=None, figsize=(2,5), outfile=None, plot_histogram=True,violin=False, labelsize=4,min_gene_flow=None,parallel=True):
        flowDF=self.flows.clone()
        

        if not genes is None:        
            flowDF = self.flows.filter( pl.col(self.symbol_column).is_in(genes) )
        
        flowScores = pl.DataFrame()

        
        def prepare_flow_calculation(self, flowDF,fgid):
            flowCols, flow = self._get_flow_columns(fgid)

            if not min_gene_flow is None:
                flowDF=flowDF.filter(pl.all(pl.col(flowCols) > min_gene_flow))
            flowScore_perGene = flowDF.select([self.symbol_column,
                            pl.struct(flowCols).apply(lambda x: np.prod(list(x.values()))).alias("pwscore")
                        ]       
                        )
            flowScore_perGene.columns=[self.symbol_column,str(fgid) ]

        
            return(flowScore_perGene)
        
        bar = makeProgressBar()

        
        if parallel:
            
            print("Starting Event Loop")
            from joblib import Parallel, delayed, parallel_backend
            with parallel_backend('loky', n_jobs=8):
                results = Parallel()(delayed(prepare_flow_calculation)(self,flowDF, fgid) for fgid in bar(use_flows))
            print("Event Loop Completed")
            
        else:
            results=[]
            for fgid in bar(use_flows):
                results.append(prepare_flow_calculation(self,flowDF,fgid))


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

        return flowScores_df    






    # __deprecated_