import ROOT 
from math import exp, pi
from HistoLab import CopyHist

def smooth2D(plot, h, norm= True):
    lenX = plot.GetNbinsX()
    lenY = plot.GetNbinsY()
    
    newplot = CopyHist(plot, True)

    for i in range(lenX):
        for j in range(lenY):
            thisPt = 0
            for reX in range(lenX):
                for reY in range(lenY):
                    thisPt += exp( -((i-reX))**2/2/h )*exp( -((j-reY))**2/2/h ) * plot.GetBinContent(reX+1, reY+1)

            newplot.SetBinContent(i+1, j+1, thisPt)
            newplot.SetBinError(i+1, j+1, plot.GetBinError(i+1, j+1))

    if norm == True:
        newplot.Scale(1/newplot.Integral())
    return newplot

def smooth1D(plot, h, norm= False):
    lenX = plot.GetNbinsX()
    
    newplot = CopyHist(plot, True)

    for i in range(lenX):
        thisPt = 0
        for reX in range(lenX):
            thisPt += exp( -((i-reX))**2/2/h )* plot.GetBinContent(reX+1)

        newplot.SetBinContent(i+1, thisPt)
        newplot.SetBinError(i+1, plot.GetBinError(i+1))

    if norm == True:
        newplot.Scale(1/newplot.Integral())
    return newplot
