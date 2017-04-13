# this import is required

from ROOT import gDirectory, TLorentzVector, TVector3, TH1D, TH2D, TFile, RooStats, TMath
from HistoLab import HistoTool, CopyHist, StackHists, GetPlot, ScaleHist, Project2D, CopyHists, project1D, NormHists, AddHistPairs, ResetBins, GetHistMedian, GetHistIQR, GetRatios, GetAveRatios, GetGausMean, CalSysVariation, CalSysDiff
from FitFunctions import BernsteinO2, BernsteinO3, BernsteinO4, BernsteinO5, ExpoBernsteinO2, ExpoBernsteinO3, ExpoBernsteinO4, ExpoBernsteinO5, Expo, Expo2, Expo3
from FitFunctions import BernsteinO2_raw, BernsteinO3_raw, BernsteinO4_raw, BernsteinO5_raw, ExpoBernsteinO2_raw, ExpoBernsteinO3_raw, ExpoBernsteinO4_raw, ExpoBernsteinO5_raw, Expo_raw, Expo2_raw, Expo3_raw
from FitFunctions import BernsteinO2Lin, BernsteinO3Lin, BernsteinO4Lin
from FitFunctions import bkgfit, bkgfit_2Region, fit_start, fit_end, fit_range
from vbfplot import mBB_Binning

#from rootpy.root2array import tree_to_ndarray
import root_numpy as rnp
import sys
from array import array
import numpy as np
import glob
from copy import deepcopy
from math import sqrt
from optparse import OptionParser
from copy import deepcopy
from MVA import * 
from Utils import *
import os

DrawTool_Diff = HistoTool()
DrawTool_Diff.lumi = "0"
DrawTool_Diff.sqrtS = "13"
DrawTool_Diff.shift = 0.11
DrawTool_Diff.doDiff = True
DrawTool_Diff.DrawRatio = "(Data-Fit) / Data"
DrawTool_Diff.SetCompData(True)


basedir = "../fitswisc/xmlAnaWSBuilder/config/vbf/"

DataFile       = "FitFiles_2cen/f_0tag-pseudo_data_2cen_CR_2tag_BDTI"
DataWeightFile = "FitFiles_2cen/f_0tag-pseudo_data_2cen_CR_2tag_weight_BDTI"


def TextToDataHist(histname, DataFile, DataWeightFile):
    hist = TH1D(histname, histname, len(mBB_Binning)-1, mBB_Binning)
    datalines = DataFile.readlines()
    dataweightlines = DataWeightFile.readlines()
    for i in range(len(datalines)):
        hist.Fill(float(datalines[i]), dataweightlines[i])

    return hist

def TextToBkgFit(POITxt, NPTxt, 
