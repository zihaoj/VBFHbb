# this import is required
import rootpy.tree
from rootpy.io import root_open
from ROOT import gDirectory, TLorentzVector, TVector3, TH1D, TH2D, TFile, RooStats, TMath
from HistoLab import HistoTool, CopyHist, StackHists, GetPlot, ScaleHist, Project2D, CopyHists, project1D, NormHists, AddHistPairs, ResetBins, GetHistMedian, GetHistIQR, GetRatios, GetAveRatios, GetGausMean, CalSysVariation, CalSysDiff
from FitFunctions import BernsteinO2, BernsteinO3, BernsteinO4, BernsteinO5, BernsteinO6, ExpoBernsteinO2, ExpoBernsteinO3, ExpoBernsteinO4, ExpoBernsteinO5, Expo, Expo2, Expo3
from FitFunctions import BernsteinO2_raw, BernsteinO3_raw, BernsteinO4_raw, BernsteinO5_raw, BernsteinO6_raw, ExpoBernsteinO2_raw, ExpoBernsteinO3_raw, ExpoBernsteinO4_raw, ExpoBernsteinO5_raw, Expo_raw, Expo2_raw, Expo3_raw, LinCustomRange_raw
from FitFunctions import BernsteinO2Lin, BernsteinO3Lin, BernsteinO4Lin, LinCustomRange
from FitFunctions import bkgfit, bkgfit_2Region, fit_start, fit_end, fit_range

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
from scipy import stats

p = OptionParser()
p.add_option('--NBDTReg', type = "string", default = '2',   dest = 'NBDTReg', help = 'number of bdt regions')
p.add_option('--BDTTotEff', type = "string", default = '0.4',   dest = 'BDTTotEff', help = 'total efficiency of BDT combined')
p.add_option('--btag', type = "string", default = 'tight',   dest = 'btag', help = 'btag working point')
p.add_option('--pTBBCut', action="store_true", default = False, dest = 'ApplypTBBCut', help = 'cut on pT(bb)')
p.add_option('--MVAList', type = "string", default = 'short',   dest = 'MVAList', help = 'MVA List')
p.add_option('--MVALeaveOut', type = "string", default = 'No',   dest = 'MVALeaveOut', help = 'MVA variable left out')
p.add_option('--test', action="store_true", default = False, dest = 'test', help = 'test mode')
p.add_option('--fit', action="store_true", default = False, dest = 'fit', help = 'fit mode')
p.add_option('--newmodel', type="string", default = 'n',  dest = 'newmodel', help = 'new model')
p.add_option('--channel', type="string", default = '4cen', dest = 'channel', help = 'channel')
p.add_option('--CR', type="string", default = '2tag', dest = 'CR', help = 'CR type')
p.add_option('--postfix', type="string", default = 'nocut', dest = 'postfix', help = 'post fix special option')
p.add_option('--period', type="string", default = 'combined', dest = 'period', help = 'combined/PreICHEP/PostICHEP')

(o,a) = p.parse_args()

lumi = "0"
if o.channel == "2cen":
    lumi = "24.5"

if o.channel == "4cen":
    if o.period == "combined":
        lumi = "24.5"
    if o.period == "PreICHEP":
        lumi = "9.4"
    if o.period == "PostICHEP":
        lumi = "13.3"
    if o.postfix == "4cencomp":
        lumi =  "0"

DrawTool = HistoTool()
DrawTool.lumi = "0"

DrawTool.sqrtS = "13"
DrawTool.SetCompData(False)
DrawTool.OutPlotDir = "Plots/"
DrawTool.lumi = lumi

DrawTool_Ratio = HistoTool()
DrawTool_Ratio.DrawRatio = "ratio"
DrawTool_Ratio.lumi = lumi
DrawTool_Ratio.sqrtS = "13"
DrawTool_Ratio.shift = 0.11
DrawTool_Ratio.SetCompData(True)

DrawTool_Diff = HistoTool()
DrawTool_Diff.lumi = lumi
DrawTool_Diff.sqrtS = "13"
DrawTool_Diff.shift = 0.11
DrawTool_Diff.doDiff = True
DrawTool_Diff.DrawRatio = "(Data-Fit) / Data"
DrawTool_Diff.SetCompData(True)

if o.ApplypTBBCut:
    os.system("rm f_*"+o.channel+"*")
    os.system("rm f_*"+o.channel+"*")
    os.system("rm workspace_*"+o.channel+"*")
    os.system("rm FitFiles_"+o.channel+"/*")
    os.system("rm FitFiles_"+o.channel+"/Poisson/*")

mBB_bin_width = 2

PDGID_Binning = array('d',range(16))
NJet_Binning = array('d',range(16))
NTrk_Binning = array('d',range(25))
NTrk_Binning_coarse = array('d',range(20))
JetpT_Binning = array('d',range(0, 400, 2))
pTJJ_Binning = array('d',range(0, 600, 20))
pTJJ_Binning_coarse = array('d',range(0, 800, 10))
MJJ_Binning = array('d',range(0, 3500, 100))
MJJ_Binning_coarse = array('d',range(0, 3500, 10))
Jeteta_Binning = array('d',[x / 10. for x in range(-56, 56, 2)])
MaxJeteta_Binning = array('d',[x / 10. for x in range(0, 56, 2)])
mBB_Binning = array('d', range(fit_start, fit_end+mBB_bin_width, mBB_bin_width))
mBB_Binning_sig = array('d', range(40, 150+mBB_bin_width, mBB_bin_width))
Eta_Binning = array('d',[x / 10. for x in range(-40, 40, 1)])
dR_Binning = array('d',[x / 10. for x in range(-15, 100, 1)])
dR_Binning_coarse = array('d',[x / 10. for x in range(0, 80, 5)])
dR_Binning_short = array('d',[x / 10. for x in range(0, 30, 1)])
JetWidth_Binning = array('d',[x / 100. for x in range(0, 40, 1)])
CosTheta_Binning = array('d',[x / 10. for x in range(-11, 11, 1)])
MV2_Binning = array('d', [x/100. for x in range(-100, 100, 1)])
BDT_Binning = array('d', [x/100. for x in range(-10, 10, 2)])
BDT_Binning_fine = array('d', [x/1000. for x in range(-200, 100, 4)])
NN_Binning = array('d', [x/100. for x in range(0, 100, 1)])
ptballance_Binning = array('d', [x/100. for x in range(0, 100, 5)])
mBB_sig_Binning = array('d', range(100, 140, mBB_bin_width))
empty_mBB_hist = TH1D("empty_mBB_hist", "empty_mBB_hist", len(mBB_Binning)-1, mBB_Binning)


### define b-tag systematics
BTAGSYSLIST = []
BTAGSYSNAMES = ["BTrig_J1_SF", "BTrig_J2_SF", "Eigen_B_0", "Eigen_B_1", "Eigen_B_2", "Eigen_B_3", "Eigen_B_4", 
                "Eigen_C_0", "Eigen_C_1", "Eigen_C_3", "Eigen_C_4",
                "Eigen_Light_0", "Eigen_Light_1", "Eigen_Light_2", "Eigen_Light_3", "Eigen_Light_4", "Eigen_Light_5",
                "Eigen_Light_6","Eigen_Light_7","Eigen_Light_8","Eigen_Light_9","Eigen_Light_10", "Eigen_Light_11","Eigen_Light_12","Eigen_Light_13", 
                "Eigen_Extrapolation", "Eigen_Extrapolation_Charm"]

for sys in BTAGSYSNAMES:
    BTAGSYSLIST.append(sys+"__1down")
    BTAGSYSLIST.append(sys+"__1up")

BTAGSYSMAP = {}
for ibsys in range(len(BTAGSYSLIST)):
    BTAGSYSMAP[ BTAGSYSLIST[ibsys] ] = [ibsys]

### define theory systematics
THSYSLIST = []
THSYSNAMES = ["TH_alpha_s", "TH_QCDScale", "TH_PDFVar"]
for sys in THSYSNAMES:
    THSYSLIST.append(sys+"__1down")
    THSYSLIST.append(sys+"__1up")


SYSLIST = ["JET_21NP_JET_EffectiveNP_1__1down", "JET_21NP_JET_EffectiveNP_1__1up",
           "JET_21NP_JET_EffectiveNP_2__1down", "JET_21NP_JET_EffectiveNP_2__1up",
           "JET_21NP_JET_EffectiveNP_3__1down", "JET_21NP_JET_EffectiveNP_3__1up",
           "JET_21NP_JET_EffectiveNP_4__1down", "JET_21NP_JET_EffectiveNP_4__1up",
           "JET_21NP_JET_EffectiveNP_5__1down", "JET_21NP_JET_EffectiveNP_5__1up",
           "JET_21NP_JET_EffectiveNP_6__1down", "JET_21NP_JET_EffectiveNP_6__1up",
           "JET_21NP_JET_EffectiveNP_7__1down", "JET_21NP_JET_EffectiveNP_7__1up",
           "JET_21NP_JET_EffectiveNP_8restTerm__1down", "JET_21NP_JET_EffectiveNP_8restTerm__1up",
           "JET_21NP_JET_EtaIntercalibration_Modelling__1down",  "JET_21NP_JET_EtaIntercalibration_Modelling__1up", 
           "JET_21NP_JET_EtaIntercalibration_TotalStat__1down",  "JET_21NP_JET_EtaIntercalibration_TotalStat__1up",
           "JET_21NP_JET_EtaIntercalibration_NonClosure__1down", "JET_21NP_JET_EtaIntercalibration_NonClosure__1up",
           "JET_21NP_JET_Pileup_OffsetMu__1down",                "JET_21NP_JET_Pileup_OffsetMu__1up",
           "JET_21NP_JET_Pileup_OffsetNPV__1down",               "JET_21NP_JET_Pileup_OffsetNPV__1up",
           "JET_21NP_JET_Pileup_PtTerm__1down",                  "JET_21NP_JET_Pileup_PtTerm__1up",
           "JET_21NP_JET_Pileup_RhoTopology__1down",             "JET_21NP_JET_Pileup_RhoTopology__1up",
           "JET_21NP_JET_Flavor_Composition__1down",             "JET_21NP_JET_Flavor_Composition__1up",
           "JET_21NP_JET_Flavor_Response__1down",                "JET_21NP_JET_Flavor_Response__1up",
           "JET_21NP_JET_BJES_Response__1down",                  "JET_21NP_JET_BJES_Response__1up",
           "JET_21NP_JET_PunchThrough_MC15__1down",              "JET_21NP_JET_PunchThrough_MC15__1up",
           "JET_21NP_JET_SingleParticle_HighPt__1down",          "JET_21NP_JET_SingleParticle_HighPt__1up",
           "JET_JER_SINGLE_NP__1up", 
           "JET_QG_trackEfficiency__1up",                             
           "JET_QG_trackFakes__1up", 
           "JET_QG_nchargedExp__1down",                          "JET_QG_nchargedExp__1up",                            
           "JET_QG_nchargedME__1down",                           "JET_QG_nchargedME__1up",
           "JET_QG_nchargedPDF__1down",                          "JET_QG_nchargedPDF__1up",
           "PRW_DATASF__1down",                                  "PRW_DATASF__1up"]


SYSNAMES = ["JET_21NP_JET_EffectiveNP_1", "JET_21NP_JET_EffectiveNP_2", 
            "JET_21NP_JET_EffectiveNP_3", "JET_21NP_JET_EffectiveNP_4",
            "JET_21NP_JET_EffectiveNP_5", "JET_21NP_JET_EffectiveNP_6",
            "JET_21NP_JET_EffectiveNP_7", "JET_21NP_JET_EffectiveNP_8restTerm", 
            "JET_21NP_JET_EtaIntercalibration_Modelling",   "JET_21NP_JET_EtaIntercalibration_TotalStat",
            "JET_21NP_JET_EtaIntercalibration_NonClosure",  "JET_21NP_JET_Pileup_OffsetMu",
            "JET_21NP_JET_Pileup_OffsetNPV",                "JET_21NP_JET_Pileup_PtTerm",
            "JET_21NP_JET_Pileup_RhoTopology",              "JET_21NP_JET_Flavor_Composition",
            "JET_21NP_JET_Flavor_Response",                 "JET_21NP_JET_BJES_Response",
            "JET_21NP_JET_PunchThrough_MC15",               "JET_21NP_JET_SingleParticle_HighPt",
            "JET_JER_SINGLE_NP",                            "JET_QG_trackEfficiency",
            "JET_QG_trackFakes",                            "JET_QG_nchargedExp",
            "JET_QG_nchargedME",                            "JET_QG_nchargedPDF",
            "PRW_DATASF"]



SAMPLES = {}
SYSSAMPLES = {}
FITSAMPLES= {}
FITRESULTS = {}


def CreateWeightSample(sample):

    weightsysarray = np.zeros([ SAMPLES[sample].var["weightSysts"].shape[0], len(SAMPLES[sample].var["weightSysts"][0]) ])

    print "create sample weights "
    for ievt in range(SAMPLES[sample].var["weightSysts"].shape[0]):
        weightsysarray[ievt] = SAMPLES[sample].var["weightSysts"][ievt]

    for ibsys in range(len(BTAGSYSLIST)):

        SAMPLES[sample+"_"+BTAGSYSLIST[ibsys]] = deepcopy(SAMPLES[sample])
        SAMPLES[sample+"_"+BTAGSYSLIST[ibsys]].name += " "
        SAMPLES[sample+"_"+BTAGSYSLIST[ibsys]].name += BTAGSYSLIST[ibsys]
        SAMPLES[sample+"_"+BTAGSYSLIST[ibsys]].isSys = True

        SAMPLES[sample+"_"+BTAGSYSLIST[ibsys]].var["eventWeight"] *=  weightsysarray[:, BTAGSYSMAP[ BTAGSYSLIST[ibsys]][0] ]
        SYSSAMPLES[sample+"_"+BTAGSYSLIST[ibsys]] = SAMPLES[sample+"_"+BTAGSYSLIST[ibsys]]

    
    if (sample != "VBF" and sample != "ggF"):
        return 

    ## alpha_s
    SAMPLES[sample+"_TH_alpha_s__1up"] = deepcopy(SAMPLES[sample])
    SAMPLES[sample+"_TH_alpha_s__1up"].name += " TH_alpha_s__1up"
    SAMPLES[sample+"_TH_alpha_s__1up"].isSys = True
    SAMPLES[sample+"_TH_alpha_s__1up"].var["eventWeight"] *=  SAMPLES[sample].var["TH_alpha_s_up"]
    SYSSAMPLES[sample+"_TH_alpha_s__1up"] = SAMPLES[sample+"_TH_alpha_s__1up"]

    SAMPLES[sample+"_TH_alpha_s__1down"] = deepcopy(SAMPLES[sample])
    SAMPLES[sample+"_TH_alpha_s__1down"].name += " TH_alpha_s__1down"
    SAMPLES[sample+"_TH_alpha_s__1down"].isSys = True
    SAMPLES[sample+"_TH_alpha_s__1down"].var["eventWeight"] *=  SAMPLES[sample].var["TH_alpha_s_dn"]
    SYSSAMPLES[sample+"_TH_alpha_s__1down"] = SAMPLES[sample+"_TH_alpha_s__1down"]

    ## NNPDF30 100NP
    THEORYSYSTS = ["weight_pdf4lhc", "weight_pdfnnpdf30", "weight_qcd_nnlops_2np", "weight_qcd_WG1"]
    THEORYSYSTS_Up = {}
    THEORYSYSTS_Dn = {}
    
    for sys in THEORYSYSTS:
        
        THEORYSYSTS_Up[sys] = np.zeros(SAMPLES[sample].var[sys].shape[0])
        THEORYSYSTS_Dn[sys] = np.zeros(SAMPLES[sample].var[sys].shape[0])

        for ievt in range(SAMPLES[sample].var[sys].shape[0]):
            THEORYSYSTS_Up[ievt] = np.max( SAMPLES[sample].var[sys][ievt, :] )
            THEORYSYSTS_Dn[ievt] = np.min( SAMPLES[sample].var[sys][ievt, :] )


    SAMPLES[sample+"_TH_QCDScale__1up"] = deepcopy(SAMPLES[sample])
    SAMPLES[sample+"_TH_QCDScale__1up"].name += " TH_QCDScale__1up"
    SAMPLES[sample+"_TH_QCDScale__1up"].isSys = True
    if sample == "VBF":
        SAMPLES[sample+"_TH_QCDScale__1up"].var["eventWeight"] *=  THEORYSYSTS_Up["weight_pdfnnpdf30"]
    if sample == "ggF":
        SAMPLES[sample+"_TH_QCDScale__1up"].var["eventWeight"] *=  THEORYSYSTS_Up["weight_pdf4lhc"]
    SYSSAMPLES[sample+"_TH_QCDScale__1up"] = SAMPLES[sample+"_TH_QCDScale__1up"]

    SAMPLES[sample+"_TH_QCDScale__1down"] = deepcopy(SAMPLES[sample])
    SAMPLES[sample+"_TH_QCDScale__1down"].name += " TH_QCDScale__1down"
    SAMPLES[sample+"_TH_QCDScale__1down"].isSys = True
    if sample == "VBF":
        SAMPLES[sample+"_TH_QCDScale__1down"].var["eventWeight"] *=  THEORYSYSTS_Dn["weight_pdfnnpdf30"]
    if sample == "ggF":
        SAMPLES[sample+"_TH_QCDScale__1down"].var["eventWeight"] *=  THEORYSYSTS_Dn["weight_pdf4lhc"]
    SYSSAMPLES[sample+"_TH_QCDScale__1down"] = SAMPLES[sample+"_TH_QCDScale__1down"]

    SAMPLES[sample+"_TH_PDFVar__1up"] = deepcopy(SAMPLES[sample])
    SAMPLES[sample+"_TH_PDFVar__1up"].name += " TH_PDFVar__1up"
    SAMPLES[sample+"_TH_PDFVar__1up"].isSys = True
    if sample == "VBF":
        SAMPLES[sample+"_TH_PDFVar__1up"].var["eventWeight"] *=  THEORYSYSTS_Up["weight_qcd_WG1"]
    if sample == "ggF":
        SAMPLES[sample+"_TH_PDFVar__1up"].var["eventWeight"] *=  THEORYSYSTS_Up["weight_qcd_nnlops_2np"]
    SYSSAMPLES[sample+"_TH_PDFVar__1up"] = SAMPLES[sample+"_TH_PDFVar__1up"]

    SAMPLES[sample+"_TH_PDFVar__1down"] = deepcopy(SAMPLES[sample])
    SAMPLES[sample+"_TH_PDFVar__1down"].name += " TH_PDFVar__1down"
    SAMPLES[sample+"_TH_PDFVar__1down"].isSys = True
    if sample == "VBF":
        SAMPLES[sample+"_TH_PDFVar__1down"].var["eventWeight"] *=  THEORYSYSTS_Dn["weight_qcd_WG1"]
    if sample == "ggF":
        SAMPLES[sample+"_TH_PDFVar__1down"].var["eventWeight"] *=  THEORYSYSTS_Dn["weight_qcd_nnlops_2np"]
    SYSSAMPLES[sample+"_TH_PDFVar__1down"] = SAMPLES[sample+"_TH_PDFVar__1down"]



def InitSamples():
    
    VBF = None
    ggF = None
    Zbb_QCD = None
    Zbb_EWK = None
    data_2tag = None
    data_0tag = None
    data_2tag_PostICHEP = None

    if o.channel == "4cen":

        if o.period == "combined":

            VBF = PhysicsProcess("VBF H->bb", "input/4cen/VBF_Reader_03_08_4cen_"+o.btag+".root")
            ggF = PhysicsProcess("ggF H->bb", "input/4cen/ggF_Reader_03_08_4cen_"+o.btag+".root")
            Zbb_QCD = PhysicsProcess("QCD Z->bb", "input/4cen/Zbb_QCD_Reader_03_08_4cen_"+o.btag+".root")
            Zbb_EWK = PhysicsProcess("EWK Z->bb", "input/4cen/Zbb_EWK_Reader_03_08_4cen_"+o.btag+".root")
            data_0tag = PhysicsProcess("data 0tag", "input/4cen/data_2016_Reader_03_08_0tag_4cen.root")
            data_2tag = PhysicsProcess("data 2tag", "input/4cen/data_2016_Reader_03_08_2tag_4cen_"+o.btag+".root")

            for sys in SYSLIST:
                systreename = sys
                if ("JET_QG_trackEfficiency" in sys or "JET_QG_trackFakes" in sys):
                    systreename = systreename.rstrip("__1up")
                
                SYSSAMPLES["VBF_"+sys] = PhysicsProcess("VBF H->bb "+sys, "input/4cen/VBF_Reader_03_08_4cen_"+o.btag+".root", systreename)
                SYSSAMPLES["ggF_"+sys] = PhysicsProcess("ggF H->bb "+sys, "input/4cen/ggF_Reader_03_08_4cen_"+o.btag+".root", systreename)
                SYSSAMPLES["Zbb_QCD_"+sys] = PhysicsProcess("QCD Z->bb "+sys, "input/4cen/Zbb_QCD_Reader_03_08_4cen_"+o.btag+".root", systreename)
                SYSSAMPLES["Zbb_EWK_"+sys] = PhysicsProcess("EWK Z->bb "+sys, "input/4cen/Zbb_EWK_Reader_03_08_4cen_"+o.btag+".root", systreename)

    if o.channel == "2cen":
        VBF = PhysicsProcess("VBF H->bb", "input/2cen/VBF_Reader_03_17_2cen.root")
        ggF = PhysicsProcess("ggF H->bb", "input/2cen/ggF_Reader_03_17_2cen.root")
        Zbb_QCD = PhysicsProcess("QCD Z->bb", "input/2cen/Zbb_QCD_Reader_03_17_2cen.root")
        Zbb_EWK = PhysicsProcess("EWK Z->bb", "input/2cen/Zbb_EWK_Reader_03_17_2cen.root")
        data_0tag = PhysicsProcess("data 0tag", "input/2cen/data_2016_Reader_03_17_0tag_2cen.root")
        data_2tag = PhysicsProcess("data 2tag", "input/2cen/data_2016_Reader_03_17_2tag_2cen.root")

        if o.postfix == "mbbcorrcheck":
            SAMPLES["VBF_OneMu"] = PhysicsProcess("VBF H->bb Muon Correction", "input/VBF_Reader_01_18_2cen_OneMu.root")
            SAMPLES["VBF_NoCorr"] = PhysicsProcess("VBF H->bb No Correction", "input/VBF_Reader_01_18_2cen_NoCorr.root")

            SAMPLES["Zbb_QCD_OneMu"] = PhysicsProcess("QCD Z->bb Muon Correction", "input/Zbb_QCD_Reader_01_20_2cen_OneMu.root")
            SAMPLES["Zbb_QCD_NoCorr"] = PhysicsProcess("QCD Z->bb No Correction", "input/Zbb_QCD_Reader_01_20_2cen_NoCorr.root")

            SAMPLES["Zbb_EWK_OneMu"] = PhysicsProcess("EWK Z->bb Muon Correction", "input/Zbb_EWK_Reader_01_20_2cen_OneMu.root")
            SAMPLES["Zbb_EWK_NoCorr"] = PhysicsProcess("EWK Z->bb No Correction", "input/Zbb_EWK_Reader_01_20_2cen_NoCorr.root")

            SAMPLES["Zbb_NoCorr"] = PhysicsProcess("QCD Z->bb No Correction", "input/Zbb_QCD_Reader_01_20_2cen_NoCorr.root")
            SAMPLES["Zbb_OneMu"] = PhysicsProcess("QCD Z->bb No Correction", "input/Zbb_QCD_Reader_01_20_2cen_NoCorr.root")

        if o.postfix == "NLOComp":
            SAMPLES["VBF_NLO"] = PhysicsProcess("VBF H->bb NLO", "input/2cen/VBF_Reader_03_17_2cen_Old.root")
            SAMPLES["ggF_NLO"] = PhysicsProcess("ggF H->bb EF4J", "input/2cen/ggF_Reader_03_17_2cen_Old.root")

        for sys in SYSLIST:

            systreename = sys
            if ("JET_QG_trackEfficiency" in sys or "JET_QG_trackFakes" in sys):
                systreename = systreename.rstrip("__1up")

            SYSSAMPLES["VBF_"+sys] = PhysicsProcess("VBF H->bb "+sys, "input/2cen/VBF_Reader_03_17_2cen.root", systreename)
            SYSSAMPLES["ggF_"+sys] = PhysicsProcess("ggF H->bb "+sys, "input/2cen/ggF_Reader_03_17_2cen.root", systreename)
            SYSSAMPLES["Zbb_QCD_"+sys] = PhysicsProcess("QCD Z->bb "+sys, "input/2cen/Zbb_QCD_Reader_03_17_2cen.root", systreename)
            SYSSAMPLES["Zbb_EWK_"+sys] = PhysicsProcess("EWK Z->bb "+sys, "input/2cen/Zbb_EWK_Reader_03_17_2cen.root", systreename)


    SAMPLES["VBF"] = VBF
    SAMPLES["ggF"] = ggF
    SAMPLES["Zbb_QCD"] = Zbb_QCD
    SAMPLES["Zbb_EWK"] = Zbb_EWK
    SAMPLES["data_0tag"] = data_0tag
    SAMPLES["data_2tag"] = data_2tag

    SAMPLES.update(SYSSAMPLES)

    FITSAMPLES["pseudo_bkg"] = PhysicsProcess("pseudo background", None)
    FITSAMPLES["pseudo_data"] = PhysicsProcess("pseudo data", None)
    FITSAMPLES["pseudo_bkg_ratio_fit"] = PhysicsProcess("pseudo background ratio fit", None)
    FITSAMPLES["pseudo_data_ratio_fit"] = PhysicsProcess("pseudo data ratio fit", None)
    FITSAMPLES["pseudo_bkg_lin_fit"] = PhysicsProcess("pseudo background lin fit", None)
    FITSAMPLES["pseudo_data_lin_fit"] = PhysicsProcess("pseudo data lin fit", None)
    FITSAMPLES["signal"] = PhysicsProcess("signal", None)

    for s in SAMPLES:

        ### Load event weight
        SAMPLES[s].AddEventVarFromTree("eventWeight", o.test)
        
        ### back up event weight used only for control region no blinding
        SAMPLES[s].var["eventWeight_noblind"]= SAMPLES[s].var["eventWeight"]

        SAMPLES[s].AddEventVarFromTree("weightSysts", o.test)

        if s == "VBF" or s == "ggF":

            SAMPLES[s].AddEventVarFromTree("weight_pdf4lhc", o.test)
            SAMPLES[s].AddEventVarFromTree("weight_pdfnnpdf30", o.test)
            SAMPLES[s].AddEventVarFromTree("weight_qcd_nnlops", o.test)
            SAMPLES[s].AddEventVarFromTree("weight_qcd_nnlops_2np", o.test)
            SAMPLES[s].AddEventVarFromTree("weight_qcd_WG1", o.test)
            SAMPLES[s].AddEventVarFromTree("weight_alternative_pdf", o.test)

        ### Load event variable
        SAMPLES[s].AddEventVarFromTree("nJets", o.test)
        SAMPLES[s].AddEventVarFromTree("pTJ1", o.test)
        SAMPLES[s].AddEventVarFromTree("pTJ2", o.test)
        SAMPLES[s].AddEventVarFromTree("mBB", o.test)
        SAMPLES[s].AddEventVarFromTree("mBB_no_corr", o.test)
        SAMPLES[s].AddEventVarFromTree("pTBB", o.test)
        SAMPLES[s].AddEventVarFromTree("mJJ", o.test)
        SAMPLES[s].AddEventVarFromTree("pTJJ", o.test)
        SAMPLES[s].AddEventVarFromTree("etaJ1", o.test)
        SAMPLES[s].AddEventVarFromTree("etaJ2", o.test)
        SAMPLES[s].AddEventVarFromTree("etaB1", o.test)
        SAMPLES[s].AddEventVarFromTree("etaB2", o.test)
        SAMPLES[s].AddEventVarFromTree("cosTheta_boost", o.test)

        SAMPLES[s].AddEventVarFromTree("mindRJ1_Ex", o.test)
        SAMPLES[s].AddEventVarFromTree("mindRJ2_Ex", o.test)
        SAMPLES[s].AddEventVarFromTree("pT_ballance", o.test)

        SAMPLES[s].AddEventVarFromTree("alpha_s_up", o.test)
        SAMPLES[s].AddEventVarFromTree("alpha_s_dn", o.test)

        SAMPLES[s].var["mindRJ1_Ex"][ SAMPLES[s].var["mindRJ1_Ex"]>10] = -1
        SAMPLES[s].var["mindRJ2_Ex"][ SAMPLES[s].var["mindRJ2_Ex"]>10] = -1

        SAMPLES[s].AddEventVarFromTree("MV2c10B1", o.test)
        SAMPLES[s].AddEventVarFromTree("MV2c10B2", o.test)
        SAMPLES[s].AddEventVarFromTree("MV2c10J1", o.test)
        SAMPLES[s].AddEventVarFromTree("MV2c10J2", o.test)

        SAMPLES[s].AddEventVarFromTree("dRBB", o.test)
        SAMPLES[s].AddEventVarFromTree("max_J1J2", o.test)
        SAMPLES[s].AddEventVarFromTree("eta_J_star", o.test)

        SAMPLES[s].AddEventVarFromTree("NTrk500PVJ1", o.test)
        SAMPLES[s].AddEventVarFromTree("NTrk500PVJ2", o.test)

        SAMPLES[s].AddEventVarFromTree("pTB1", o.test)
        SAMPLES[s].AddEventVarFromTree("pTB2", o.test)

        SAMPLES[s].AddNoCut()

        SAMPLES[s].AddCut("2b_mv2c10_70",  np.logical_and((SAMPLES[s].var["MV2c10B1"]>0.8244273),(SAMPLES[s].var["MV2c10B2"]>0.8244273))==1)
        SAMPLES[s].AddCut("0b_mv2c10_85",  np.logical_or( np.logical_or((SAMPLES[s].var["MV2c10B1"]>0.1758475),
                                                                        (SAMPLES[s].var["MV2c10B2"]>0.1758475)), 
                                                          np.logical_or((SAMPLES[s].var["MV2c10J1"]>0.1758475),
                                                                        (SAMPLES[s].var["MV2c10J2"]>0.1758475)))==0)


        SAMPLES[s].AddCut("pTBB>70", SAMPLES[s].var["pTBB"]>70)
        SAMPLES[s].AddCut("pTBB>100", SAMPLES[s].var["pTBB"]>100)
        SAMPLES[s].AddCut("pTBB>110", SAMPLES[s].var["pTBB"]>110)
        SAMPLES[s].AddCut("pTBB>130", SAMPLES[s].var["pTBB"]>130)
        SAMPLES[s].AddCut("pTBB>150", SAMPLES[s].var["pTBB"]>150)
        SAMPLES[s].AddCut("pTBB>160", SAMPLES[s].var["pTBB"]>160)
        SAMPLES[s].AddCut("pTBB>170", SAMPLES[s].var["pTBB"]>170)
        SAMPLES[s].AddCut("pTBB>190", SAMPLES[s].var["pTBB"]>190)
        SAMPLES[s].AddCut("pTBB>210", SAMPLES[s].var["pTBB"]>210)


        print ("sample", s)
        print ("0tag", sum(SAMPLES[s].cut["0b_mv2c10_85"].array))
        print ("2tag", sum(SAMPLES[s].cut["2b_mv2c10_70"].array))

        ## add sidebands cut
        blindarray_low = np.logical_and( SAMPLES[s].var["mBB"]<100, SAMPLES[s].var["mBB"]>80)
        blindarray_up = np.logical_and( SAMPLES[s].var["mBB"]<190, SAMPLES[s].var["mBB"]>150)

        fitarray_low = np.logical_and( SAMPLES[s].var["mBB"]<100, SAMPLES[s].var["mBB"]>fit_start)
        fitarray_up = np.logical_and( SAMPLES[s].var["mBB"]<fit_end, SAMPLES[s].var["mBB"]>140)

        SAMPLES[s].AddCut("sideband_low", blindarray_low)
        SAMPLES[s].AddCut("sideband_up", blindarray_up)
        SAMPLES[s].AddCut("sideband", np.logical_or(blindarray_low, blindarray_up))

        SAMPLES[s].AddCut("fitband_low", fitarray_low)
        SAMPLES[s].AddCut("fitband_up", fitarray_up)
        SAMPLES[s].AddCut("fitband", np.logical_or(fitarray_low, fitarray_up))

        ## add signal region cut
        signalarray = np.logical_and( SAMPLES[s].var["mBB"]>100, SAMPLES[s].var["mBB"]<140)
        SAMPLES[s].AddCut("signal", signalarray)
        SAMPLES[s].AddCut("blindsingalregion", np.logical_or( SAMPLES[s].var["mBB"]<100, SAMPLES[s].var["mBB"]>140) )
        SAMPLES[s].AddCut("blindsingalregion_no_corr", np.logical_or( SAMPLES[s].var["mBB_no_corr"]<100, SAMPLES[s].var["mBB_no_corr"]>140) )

def ScaleLumi():
    for s in SAMPLES:
        if "data" not in s:
            if o.channel == "2cen":
                SAMPLES[s].var["eventWeight"] *= 24.5
            if o.channel == "4cen":
                if o.period == "PreICHEP":
                    SAMPLES[s].var["eventWeight"] *= 1.0
                if o.period == "PostICHEP":
                    SAMPLES[s].var["eventWeight"] *= 1.0
    
def BlindData():
    for s in SAMPLES:
        if "data_2tag" in s:
            blindarray = np.logical_or( SAMPLES[s].var["mBB"]<100, SAMPLES[s].var["mBB"]>140)
            SAMPLES[s].var["eventWeight"] *= blindarray

def ApplyBtag():
    for s in SAMPLES:
        btagarray = None
        if "0tag" in s:
            btagarray = SAMPLES[s].cut["0b_mv2c10_85"].array
            SAMPLES[s].var["eventWeight"] *= btagarray

def ApplypTBBCut():
    if o.channel == "2cen":
        pTCut = 160
    if o.channel == "4cen":
        pTCut = 150
    for s in SAMPLES:
        SAMPLES[s].var["eventWeight"] = (SAMPLES[s].var["pTBB"]>pTCut) * SAMPLES[s].var["eventWeight"]
        if not SAMPLES[s].isSys:
            print "SAMPLE", s, " pT BB cutflow ", sum( SAMPLES[s].var["eventWeight"]>0 )

        SAMPLES[s].var["eventWeight_noblind"] = (SAMPLES[s].var["pTBB"]>pTCut) * SAMPLES[s].var["eventWeight_noblind"]
        SAMPLES[s].var["eventWeight_flat"] = (SAMPLES[s].var["pTBB"]>pTCut) * (SAMPLES[s].var["eventWeight"]>0)


def Make1DPlots(samples, histname, var, weight, cuts, binning, precuts=[]):

    for s in samples:

        if ("mBB" not in histname) and (samples[s].isSys):
            continue
            
        samples[s].Add1DHist(histname, binning, NoCut)
        samples[s].FillHist (histname, samples[s].var[var], samples[s].var[weight])

        for cut in cuts:
            samples[s].Add1DHist( histname, binning, samples[s].cut[cut])
            samples[s].FillHist ( histname, samples[s].var[var], samples[s].var[weight], samples[s].cut[cut])
            for precut in precuts:
                samples[s].Add1DHist( histname , binning, samples[s].cut[precut]+samples[s].cut[cut] )
                samples[s].FillHist ( histname , samples[s].var[var], samples[s].var[weight], samples[s].cut[precut] +samples[s].cut[cut])


def MakeCutFlow(histname, cuts):

    normedHists = []
    labels = []
    for s in SAMPLES:
        SAMPLES[s].Add1DHist(histname,  array('d', range( len(cuts)+1 )) )

        for icut in range(len(cuts)):
            SAMPLES[s].histograms[histname].GetXaxis().SetBinLabel(icut+1, cuts[icut])
            nEvts = sum(SAMPLES[s].var["eventWeight"]*SAMPLES[s].cut[cuts[icut]].array)
            SAMPLES[s].histograms[histname].SetBinContent(icut+1, nEvts)

        DrawTool.DrawHists(s+"_"+histname,  ["", "Events"], [SAMPLES[s].histograms[histname]], [s])

        SAMPLES[s].histograms[histname+"_Normed"] = CopyHist(SAMPLES[s].histograms[histname])
        SAMPLES[s].histograms[histname+"_Normed"].Scale(1.0/SAMPLES[s].histograms[histname].GetBinContent(1))

        normedHists.append( SAMPLES[s].histograms[histname+"_Normed"] )
        labels.append(s)

    DrawTool.DrawHists(histname +"_Fractional",   ["", "Fraction"], normedHists, labels)


def DrawMVAInputComparison(varlist, cuts):

    cuts = ["_"+x for x in cuts]
    cuts += [""]
    
    for var in varlist:
        for s in SAMPLES:
            if SAMPLES[s].isSys:
                continue

            SAMPLES[s].var[var+"_corr"] = np.zeros(len(varlist))
            ivar2 = 0
            for var2 in varlist:
                corrfac = CalWeightedCorrCoeff( SAMPLES[s].var[var], SAMPLES[s].var[var2], SAMPLES[s].var["eventWeight"] )
                SAMPLES[s].var[var+"_corr"][ivar2] = corrfac

                if o.channel == "2cen":
                    if var == "NTrk500PVJ1" or var2 == "NTrk500PVJ1":
                        SAMPLES[s].var[var+"_corr"][ivar2] = 0
                    if var == "NTrk500PVJ1" and var2 == "NTrk500PVJ1":
                        SAMPLES[s].var[var+"_corr"][ivar2] = 1

                ivar2 += 1

            for cut in cuts:
                SAMPLES[s].Norm1DHist(var+cut)
                SAMPLES[s].Norm1DHist(var+"_sideband_low"+cut)
                SAMPLES[s].Norm1DHist(var+"_sideband_up"+cut)

        for cut in cuts:
            DrawTool.DrawHists(var+"_mvacomp"+cut,  [var, "A.U."], [SAMPLES["VBF"].histograms[var+cut+"_Normed"], SAMPLES["data_2tag"].histograms[var+cut+"_Normed"]],
                               ["signal", "background"])

            DrawTool.DrawHists(var+"_mvacomp_sidband"+cut,  [var, "A.U."], [ SAMPLES["data_2tag"].histograms[var+"_sideband_low"+cut+"_Normed"], 
                                                                             SAMPLES["data_2tag"].histograms[var+"_sideband_up"+cut+"_Normed"], 
                                                                             SAMPLES["VBF"].histograms[var+cut+"_Normed"]],
                               ["data 2tag sideband low", "data 2tag sideband up", "VBF"])


    def drawBDTMVAInput2D(s, var, xbin, ybin, cut = "NoCut"):

        histname = "BDTScore_"+var+"_"+cut
        histname_norm = "BDTScore_RowNorm_"+var+"_"+cut

        if cut == "NoCut":
            histname = "BDTScore_"+var
            histname_norm = "BDTScore_RowNorm_"+var
            
        SAMPLES[s].Add2DHist(histname, xbin, ybin)
        SAMPLES[s].Add2DHist(histname_norm, xbin, ybin)

        var_array = np.dstack( ( MVAPRED["BDT_"+s+"_full"], SAMPLES[s].var[var] ))[0]
        #SAMPLES[s].FillHist( histname, var_array, SAMPLES[s].var["eventWeight"], SAMPLES[s].cut[cut])
        SAMPLES[s].FillHist( histname, var_array, SAMPLES[s].var["eventWeight"]* SAMPLES[s].cut[cut].array)
        
        nxbins = SAMPLES[s].histograms[histname].GetNbinsX()
        nybins = SAMPLES[s].histograms[histname].GetNbinsY()

        for ibiny in range(nybins):
            sumthisrow = 0

            for ibinx in range(nxbins):
                sumthisrow +=  SAMPLES[s].histograms[histname].GetBinContent( ibinx+1, ibiny+1)
            
            for ibinx in range(nxbins):
                thisbin = SAMPLES[s].histograms[histname].GetBinContent( ibinx+1, ibiny+1)
                if sumthisrow != 0:
                    SAMPLES[s].histograms[histname_norm].SetBinContent( ibinx+1, ibiny+1, thisbin/float(sumthisrow) )

        if var == "mBB":
            print "sample", s, "BDTScore Mbb correlation ", histname, SAMPLES[s].histograms[histname].GetCorrelationFactor()
            print "sample", s, "BDTScore Mbb correlation norm ", histname, SAMPLES[s].histograms[histname_norm].GetCorrelationFactor()

        DrawTool.DrawHists(s+histname,  ["BDT Score", var], [ SAMPLES[s].histograms[histname]],   [""])
        DrawTool.DrawHists(s+histname_norm,  ["BDT Score", var], [ SAMPLES[s].histograms[histname_norm]],   [""])

    
    for s in SAMPLES:
        if s != "VBF" and s != "data_2tag":
            continue

        SAMPLES[s].Add2DHist("MVAInput_Corr", np.array(range(len(varlist))), np.array(range(len(varlist))) )

        DrawTool.draw2Dtext = True
        
        for ivar in range(len(varlist)):

            SAMPLES[s].histograms["MVAInput_Corr"].GetXaxis().SetBinLabel(ivar+1, varlist[ivar])

            for jvar in range(len(varlist)):
                SAMPLES[s].histograms["MVAInput_Corr"].GetYaxis().SetBinLabel(jvar+1, varlist[jvar])
                SAMPLES[s].histograms["MVAInput_Corr"].SetBinContent(ivar+1, jvar+1, round(SAMPLES[s].var[varlist[ivar]+"_corr"][jvar], 4))
                
        SAMPLES[s].histograms["MVAInput_Corr"].GetZaxis().SetRangeUser(-1,1)
        DrawTool.DrawHists(s+"MVAInput_Corr_Canv",  ["", ""], [ SAMPLES[s].histograms["MVAInput_Corr"]],    [""])

        DrawTool.draw2Dtext = False
        if s == "VBF" or s == "data_2tag":
            drawBDTMVAInput2D(s, "mJJ", BDT_Binning, pTJJ_Binning_coarse)
            drawBDTMVAInput2D(s, "pTJJ", BDT_Binning, pTJJ_Binning_coarse)
            drawBDTMVAInput2D(s, "cosTheta_boost", BDT_Binning, CosTheta_Binning)
            drawBDTMVAInput2D(s, "mindRJ1_Ex", BDT_Binning, dR_Binning_coarse)
            drawBDTMVAInput2D(s, "mindRJ2_Ex", BDT_Binning, dR_Binning_coarse)
            drawBDTMVAInput2D(s, "max_J1J2", BDT_Binning, MaxJeteta_Binning)
            drawBDTMVAInput2D(s, "eta_J_star", BDT_Binning, Jeteta_Binning)
            drawBDTMVAInput2D(s, "NTrk500PVJ1", BDT_Binning, NTrk_Binning_coarse)
            drawBDTMVAInput2D(s, "NTrk500PVJ2", BDT_Binning, NTrk_Binning_coarse)

            SAMPLES[s].Add2DHist("pTBB_mBB", pTJJ_Binning_coarse, pTJJ_Binning_coarse)
            pTBB_mBB = np.dstack( (SAMPLES[s].var["pTBB"], SAMPLES[s].var["mBB"]))[0]
            SAMPLES[s].FillHist( "pTBB_mBB", pTBB_mBB, SAMPLES[s].var["eventWeight"])
            SAMPLES[s].histograms["pTBB_mBB"].Scale( 1/ SAMPLES[s].histograms["pTBB_mBB"].Integral() )
            DrawTool.DrawHists(s+"pTBB_mBB",  ["pTBB", "mBB"], [ SAMPLES[s].histograms["pTBB_mBB"]],    [""])

            drawBDTMVAInput2D(s, "mBB", BDT_Binning_fine, JetpT_Binning)
            drawBDTMVAInput2D(s, "mBB", BDT_Binning_fine, JetpT_Binning, "CR")
            drawBDTMVAInput2D(s, "mBB", BDT_Binning_fine, JetpT_Binning, "SRI")
            drawBDTMVAInput2D(s, "mBB", BDT_Binning_fine, JetpT_Binning, "SRII")
            if o.NBDTReg ==3:
                drawBDTMVAInput2D(s, "mBB", BDT_Binning_fine, JetpT_Binning, "SRIII")


def DrawMVAInputComparisonDifferentSamples(varlist):
    Make1DPlots(SAMPLES, "nJets",      "nJets",           "eventWeight", ["sideband"], NJet_Binning)
    Make1DPlots(SAMPLES, "mBB",        "mBB",             "eventWeight", ["sideband"], JetpT_Binning)
    Make1DPlots(SAMPLES, "mJJ",        "mJJ",             "eventWeight", ["sideband"], MJJ_Binning)
    Make1DPlots(SAMPLES, "pTJJ",       "pTJJ",            "eventWeight", ["sideband"], pTJJ_Binning)
    Make1DPlots(SAMPLES, "NTrk500PVJ1","NTrk500PVJ1",     "eventWeight", ["sideband"], NTrk_Binning)
    Make1DPlots(SAMPLES, "NTrk500PVJ2","NTrk500PVJ2",     "eventWeight", ["sideband"], NTrk_Binning)
    Make1DPlots(SAMPLES, "cosTheta_boost",   "cosTheta_boost",        "eventWeight", ["sideband"], CosTheta_Binning)
    Make1DPlots(SAMPLES, "mindRJ1_Ex", "mindRJ1_Ex",      "eventWeight", ["sideband"], dR_Binning)
    Make1DPlots(SAMPLES, "mindRJ2_Ex", "mindRJ2_Ex",      "eventWeight", ["sideband"], dR_Binning)
    Make1DPlots(SAMPLES, "max_J1J2",   "max_J1J2",        "eventWeight", ["sideband"], MaxJeteta_Binning)
    Make1DPlots(SAMPLES, "eta_J_star", "eta_J_star",      "eventWeight", ["sideband"], Jeteta_Binning)
    Make1DPlots(SAMPLES, "pT_balance", "pT_balance",    "eventWeight", ["sideband"], ptballance_Binning)
    Make1DPlots(SAMPLES, "mindRJ1_Ex_Noweight", "mindRJ1_Ex",    "eventWeight_flat", ["sideband"], dR_Binning)
    Make1DPlots(SAMPLES, "mindRJ2_Ex_Noweight", "mindRJ2_Ex",    "eventWeight_flat", ["sideband"], dR_Binning)

    varlist.append("mBB")
    varlist.append("nJets")

    for var in varlist:
        for s in SAMPLES:
            print "sample", s, "var", var
            if SAMPLES[s].isSys:
                continue
            SAMPLES[s].Norm1DHist(var)
            SAMPLES[s].Norm1DHist(var+"_sideband")
        
        if o.postfix == "4cencomp":
            DrawTool_Ratio.DrawHists(var+"_ICHEPComp",  [var, "A.U."], [SAMPLES["data_2tag_PreICHEP"].histograms[var+"_sideband_Normed"], SAMPLES["data_2tag_PostICHEP"].histograms[var+"_sideband_Normed"]],    ["Pre-ICHEP", "Post-ICHEP"])

        if o.postfix == "NLOComp":
            DrawTool_Ratio.DrawHists(var+"_sideband_VBFSamples",  [var, "A.U."], [SAMPLES["VBF"].histograms[var+"_sideband_Normed"], SAMPLES["VBF_NLO"].histograms[var+"_sideband_Normed"]], ["VBF NNPDF30", "VBF NLO"])
            DrawTool_Ratio.DrawHists(var+"_sideband_ggFSamples",  [var, "A.U."], [SAMPLES["ggF"].histograms[var+"_sideband_Normed"], SAMPLES["ggF_NLO"].histograms[var+"_sideband_Normed"]], ["ggF NNLO", "ggF NLO EF4J"])

            DrawTool_Ratio.DrawHists(var+"_VBFSamples",  [var, "A.U."], [SAMPLES["VBF"].histograms[var+"_Normed"], SAMPLES["VBF_NLO"].histograms[var+"_Normed"]], ["VBF NNPDF30", "VBF NLO"])
            DrawTool_Ratio.DrawHists(var+"_ggFSamples",  [var, "A.U."], [SAMPLES["ggF"].histograms[var+"_Normed"], SAMPLES["ggF_NLO"].histograms[var+"_Normed"]], ["ggF NNLO", "ggF NLO EF4J"])
            
    #SAMPLES["VBF"].Norm1DHist("mindRJ1_Ex_Noweight")
    #SAMPLES["VBF"].Norm1DHist("mindRJ2_Ex_Noweight")
    #DrawTool.DrawHists("mindRJ1_Ex_flatweight",  ["mindRJ1_Ex", "A.U."], [SAMPLES["VBF"].histograms["mindRJ1_Ex_Noweight_Normed"], SAMPLES["VBF"].histograms["mindRJ1_Ex_Normed"]],    ["MC W/O Weight", "MC W/ Weight"])
    #DrawTool.DrawHists("mindRJ2_Ex_flatweight",  ["mindRJ2_Ex", "A.U."], [SAMPLES["VBF"].histograms["mindRJ2_Ex_Noweight_Normed"], SAMPLES["VBF"].histograms["mindRJ2_Ex_Normed"]],    ["MC W/O Weight", "MC W/ Weight"])

        
def ComputeMVAWeights(varlist):

    MVAinputs = {}
    MVAinputs_sys = {}

    for s in SAMPLES:
        out_arr = SAMPLES[s].var["eventWeight"]
        if "data" in s: ## data apply side band cuts
            out_arr = SAMPLES[s].var["eventWeight"]*SAMPLES[s].cut["sideband"].array

        for var in varlist:
            out_arr = np.vstack( (out_arr, SAMPLES[s].var[var]) )

        if SAMPLES[s].isSys:
            MVAinputs_sys[s] = out_arr
        else:
            MVAinputs[s] = out_arr

    dataset_BDT = None
    dataset_BDT_sys= None

    dataset_BDT, dataset_BDT_sys = makeData(False, False, MVAinputs["data_2tag"], MVAinputs["data_0tag"], MVAinputs["VBF"], MVAinputs["ggF"], 
                                            MVAinputs["Zbb_QCD"], MVAinputs["Zbb_EWK"],MVAinputs_sys )

    model_BDT = None

    if o.newmodel == "y":
        print "rebuidling BDT"
        print dataset_BDT["data_2tag_full"].shape
        
        model_BDT  = buildBDT(dataset_BDT, o.channel, o.period, o.btag, doTest = o.test)

    else:
        filename = "BDT"+o.channel+o.period+o.btag+'.sav'
        if o.MVALeaveOut != "No":
            filename = "BDT"+o.channel+o.period+"_No"+o.MVALeaveOut+'.sav'
        model_BDT = pickle.load(open(filename, 'rb'))
    
    PredictModel(model_BDT, "BDT", dataset_BDT, dataset_BDT_sys)

    output.cd()

    DrawTool.DrawHists("BDTscorePlot", ["BDT score", "A.U."], 
                       NormHists([ MVAScoreHists["BDT_data_2tag"], MVAScoreHists["BDT_VBF"], MVAScoreHists["BDT_ggF"],
                                   MVAScoreHists["BDT_Zbb_QCD"], MVAScoreHists["BDT_Zbb_EWK"] ]),
                       ["Background", "VBF","ggF", "Zbb QCD", "Zbb EWK"])


    d, pval = stats.ks_2samp( MVAScoreHists["BDT_sig_train"], MVAScoreHists["BDT_sig_test"]  )
    print 'ks bkg', d, pval
    d, pval = stats.ks_2samp( MVAScoreHists["BDT_bkg_train"], MVAScoreHists["BDT_bkg_test"]  )
    print 'ks sig', d, pval

    MVAScoreHists["BDT_Z"] = deepcopy(MVAScoreHists["BDT_Zbb_QCD"])
    MVAScoreHists["BDT_Z"].Add(MVAScoreHists["BDT_Zbb_EWK"])
    
    DrawTool.DrawHists("BDTscorePlot_Short", ["BDT score", "A.U."], 
                       NormHists([ MVAScoreHists["BDT_sig_train"], MVAScoreHists["BDT_sig_test"], MVAScoreHists["BDT_bkg_train"], MVAScoreHists["BDT_bkg_test"]]),
                       ["Signal Training Sample", "Signal Test Sample","Background Training Sample", "Background Test Sample"])
    

    DrawTool.getROC( [ MVAPRED["BDT_sig_train"], MVAPRED["BDT_sig_test"]],
                     [ MVAPRED["BDT_bkg_train"], MVAPRED["BDT_bkg_test"]],
                     ["BDT train", "BDT test"] )

    BDT_Cut_00 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 1)
    BDT_Cut_01 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.01)
    BDT_Cut_015= ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.015)
    BDT_Cut_02 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.02)
    BDT_Cut_03 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.03)
    BDT_Cut_05 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.05)
    BDT_Cut_10 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.1)
    BDT_Cut_20 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.2)
    BDT_Cut_40 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.4)
    BDT_Cut_60 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.6)
    BDT_Cut_80 = ScanGetThreshold(MVAPRED["BDT_VBF_full"],  SAMPLES["VBF"].var["eventWeight"], 0.8)

    print MVAPRED["BDT_VBF_full"].shape
    print SAMPLES["VBF"].var["eventWeight"]

    BDT_Opt_Cuts = FindOptimalBDTCut(  MVAPRED["BDT_VBF_full"] , SAMPLES["VBF"].var["eventWeight"],MVAPRED["BDT_data_2tag_full"], 
                                       MVAPRED["BDT_data_0tag_full"], SAMPLES["data_0tag"].cut["fitband"].array, SAMPLES["data_0tag"].cut["signal"].array, float(o.BDTTotEff), int(o.NBDTReg))
    
    for s in SAMPLES:

        SAMPLES[s].AddCut("PassBDTCut00-01",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_01,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut00-015",   np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_015,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut00-02",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_02,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut00-03",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_03,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut00-05",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_05,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut00-10",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_10,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut00-20",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_20,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut20-40",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_40,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_20)  )
        SAMPLES[s].AddCut("PassBDTCut40-60",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_60,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_40)  )
        SAMPLES[s].AddCut("PassBDTCut60-80",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_80,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_60)  )
        SAMPLES[s].AddCut("PassBDTCut80-100",    MVAPRED["BDT_"+s+"_full"]>BDT_Cut_80 )


        if o.channel == "4cen":
            SAMPLES[s].AddCut("CR",      np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_02,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
            if s == "VBF":
                print "CR BDT CUT", BDT_Cut_02

        if o.channel == "2cen":
            SAMPLES[s].AddCut("CR",      np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_03,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
            if s == "VBF":
                print "CR BDT CUT", BDT_Cut_03

        
        SAMPLES[s].AddCut("SRI",    MVAPRED["BDT_"+s+"_full"]>BDT_Opt_Cuts[0] )
        SAMPLES[s].AddCut("SRII",   np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Opt_Cuts[0],  MVAPRED["BDT_"+s+"_full"]>BDT_Opt_Cuts[1]) )
        if int(o.NBDTReg) ==3:
            SAMPLES[s].AddCut("SRIII",  np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Opt_Cuts[1],  MVAPRED["BDT_"+s+"_full"]>BDT_Opt_Cuts[2]) )

        if SAMPLES[s].isSys == False:
            print s, "CR fraction",  sum( SAMPLES[s].cut["CR"].array*SAMPLES[s].var["eventWeight"])/ sum( SAMPLES[s].var["eventWeight"]) 


def CalculateSensitivity(MVAcuts):

    ## Calculate normalization factor
    for cut in MVAcuts:

        print ("------------------")
        print (cut)

        sideband_0tag_full_evts = sum(SAMPLES["data_0tag"].var["eventWeight"]* (SAMPLES["data_0tag"].cut["fitband"].array) )    
        sideband_0tag_evts = sum(SAMPLES["data_0tag"].var["eventWeight"]* SAMPLES["data_0tag"].cut[cut].array* (SAMPLES["data_0tag"].cut["fitband"].array) )
        sideband_2tag_evts = sum(SAMPLES["data_2tag"].var["eventWeight"]* SAMPLES["data_2tag"].cut[cut].array* (SAMPLES["data_2tag"].cut["fitband"].array) )

        scale_factor = sideband_2tag_evts / sideband_0tag_evts

        print ("2tag sidebands", sideband_2tag_evts)
        print ("0tag sidebands", sideband_0tag_evts)
        print ("0tag full sidebands", sideband_0tag_full_evts)
        print ("scale factor", sideband_2tag_evts / sideband_0tag_evts)
        
        nvbf = sum(SAMPLES["VBF"].var["eventWeight"] * SAMPLES["VBF"].cut["signal"].array * SAMPLES["VBF"].cut[cut].array)
        nggf = sum(SAMPLES["ggF"].var["eventWeight"] * SAMPLES["ggF"].cut["signal"].array * SAMPLES["ggF"].cut[cut].array)
        n0tag = sum(SAMPLES["data_0tag"].var["eventWeight"] * SAMPLES["data_0tag"].cut["signal"].array * SAMPLES["data_0tag"].cut[cut].array * scale_factor )
        nzqcd = sum(SAMPLES["Zbb_QCD"].var["eventWeight"] * SAMPLES["Zbb_QCD"].cut["signal"].array * SAMPLES["Zbb_QCD"].cut[cut].array)
        nzewk = sum(SAMPLES["Zbb_EWK"].var["eventWeight"] * SAMPLES["Zbb_EWK"].cut["signal"].array * SAMPLES["Zbb_EWK"].cut[cut].array)

        fithist = SAMPLES["data_2tag"].histograms["mBB_short_"+cut]
        fithist.Fit("BernsteinO3")
        fit = fithist.GetFunction("BernsteinO3")
        n2tag = 0

        for ibin in range(len(mBB_sig_Binning)):
            bincenter = mBB_sig_Binning[ibin]+mBB_bin_width/2.
            n2tag += fit.Eval(bincenter)

        print ("N 2tag sideband", sideband_2tag_evts)
        print ("N VBF",  nvbf)
        print ("N ggF",  nggf)
        print ("N zewk",  nzewk)
        print ("N zqcd",  nzqcd)
        print ("N 0tag extrapolated", n0tag)
        print ("S/SQRT(B)", (nvbf+nggf)/sqrt(n0tag))
        print ("S/SQRT(B*(1+0.0001B))", (nvbf+nggf)/sqrt(n0tag*(1+0.0001*n0tag)))
        print ("N 2tag fit", n2tag)
        print ("S/SQRT(B)", (nvbf+nggf)/sqrt(n2tag))
        print ("S/SQRT(B*(1+0.0001B))", (nvbf+nggf)/sqrt(n2tag*(1+0.0001*n2tag)))

        print ("------------------")


def DrawPlotsWithCutsForEachSample(var, canvname, xLabel, yLabel, cuts, dosys = False, norm=False):

    for s in SAMPLES:
        plots = []
        labels = []

        if not dosys:
            if "1down" in s:
                continue
            if "1up" in s:
                continue

        plots.append(SAMPLES[s].histograms[var])
        labels.append("pre-selection")

        for cut in cuts:
            if norm:
                SAMPLES[s].Norm1DHist(var+"_"+cut)    
                plots.append(SAMPLES[s].histograms[var+"_"+cut+"_Normed"])
            else:
                plots.append(SAMPLES[s].histograms[var+"_"+cut])
            labels.append(s+"_"+cut)

        DrawTool.DrawHists(var+"_"+s+"_"+canvname,  [xLabel, yLabel], plots, ["pre-selection"]+cuts)


def DrawPlotsPreCutMultiSample(var, canvname, xLabel, yLabel, norm = True):

    plots = []
    labels = []

    for s in SAMPLES:
        if s != "data_2tag" and s!= "VBF":
            continue

        labels.append(s)

        if norm:
            SAMPLES[s].Norm1DHist(var)
            plots.append(SAMPLES[s].histograms[var+"_Normed"])
        else:
            plots.append(SAMPLES[s].histograms[var])

    DrawTool.DrawHists(var+"_"+canvname,  [xLabel, yLabel], plots, labels)

    plots = []
    labels = []

    for s in SAMPLES:
        if s != "Zbb_QCD" and s!= "Zbb_EWK":
            continue

        labels.append(s)

        if norm:
            SAMPLES[s].Norm1DHist(var)
            plots.append(SAMPLES[s].histograms[var+"_Normed"])
        else:
            plots.append(SAMPLES[s].histograms[var])

    DrawTool.DrawHists(var+"_"+canvname+"_Z",  [xLabel, yLabel], plots, labels)


def GenerateSignalWS(process = "signal", cut = None, useFullDistribution = False, useSys = True):

    xmlfile = None
    if useSys:
        xmlfile = open("FitFiles_"+o.channel+"/"+process+"_"+cut+"_"+o.channel+".xml", "a")
    else:
        xmlfile = open("FitFiles_"+o.channel+"/"+process+"_"+cut+"_"+o.channel+"_nosys.xml", "a")

    btagging = open("FitFiles_"+o.channel+"/"+"btag_"+cut+"_"+o.channel+".txt", "a")

    xmlfile.write('''<!DOCTYPE Model SYSTEM '../AnaWSBuilder.dtd'>'''+"\n")
    if useSys:
        xmlfile.write('''    <Model Type="External" Input="config/vbf/FitFiles_{!s}/WS_combined_{!s}_{!s}_{!s}_model.root" WSName="combined" ModelName="C_{!s}_model" ObservableName="obs_x_C_{!s}">'''.format(o.channel, process, o.channel, cut, cut, cut)+"\n")
    else:
        xmlfile.write('''    <Model Type="External" Input="config/vbf/FitFiles_{!s}/WS_combined_{!s}_{!s}_{!s}_model_nosys.root" WSName="combined" ModelName="C_{!s}_model" ObservableName="obs_x_C_{!s}">'''.format(o.channel, process, o.channel, cut, cut, cut)+"\n")
    xmlfile.write('''    <Item Name="unit[1]"/>\n''')
    xmlfile.write('''    <Rename OldName="Lumi" NewName="unit"/>\n''')

    signal_hist = None
    if process == "signal":
        signal_hist = deepcopy(  SAMPLES["VBF"].histograms["mBB_short_"+cut] )
        signal_hist.Add       (  SAMPLES["ggF"].histograms["mBB_short_"+cut] )
        signal_hist.SetName("signal_"+cut)
        
    if process == "z":
        if useFullDistribution:
            signal_hist = deepcopy(  SAMPLES["Zbb_QCD"].histograms["mBB_short"] )
            signal_hist.Add       (  SAMPLES["Zbb_EWK"].histograms["mBB_short"] )
            signal_hist.Scale(  (SAMPLES["Zbb_QCD"].histograms["mBB_short_"+cut].Integral()+
                                 SAMPLES["Zbb_EWK"].histograms["mBB_short_"+cut].Integral())/ signal_hist.Integral() )
            signal_hist.SetName("z_"+cut)
            fit_results_txt = open("FitFiles_"+o.channel+"/f_"+cut+"_full_"+o.channel, "a")
            fit_results_txt.write( "Z "+str( signal_hist.Integral())+"\n")

        else:
            signal_hist = deepcopy(  SAMPLES["Zbb_QCD"].histograms["mBB_short_"+cut] )
            signal_hist.Add       (  SAMPLES["Zbb_EWK"].histograms["mBB_short_"+cut] )
            signal_hist.SetName("z_"+cut)
            fit_results_txt = open("FitFiles_"+o.channel+"/f_"+cut+"_"+o.channel, "a")
            fit_results_txt.write( "Z "+str( signal_hist.Integral())+"\n")

    ### generate workspace for signal
    meas_sig = RooStats.HistFactory.Measurement(process+"_"+o.channel+"_"+cut, process+"_"+o.channel+"_"+cut);
    meas_sig.SetOutputFilePrefix("WS")
    meas_sig.SetExportOnly(1);
    meas_sig.SetPOI("xs");
    meas_sig.SetLumi(1);
    meas_sig.SetLumiRelErr(0.05);

    R_Channel_cut = RooStats.HistFactory.Channel("C_"+cut);
    R_Sample_cut = RooStats.HistFactory.Sample("S_"+cut);
    R_Sample_cut.SetHisto(signal_hist);

    doStat = False
    doSyst = True
    if process == "signal":
        doStat = False
        doSyst = True

    if doStat:
        R_Channel_cut.SetStatErrorConfig(0.01, "Poisson");
        R_Sample_cut.ActivateStatError();
        
        for ibin in range(signal_hist.GetNbinsX()):
            if signal_hist.GetBinContent(ibin+1) ==0:
                continue
            xmlfile.write('''    <ExtSyst ConstrName="gamma_stat_C_{!s}_bin_{!s}_constraint" NPName="gamma_stat_C_{!s}_bin_{!s}" GOName="nom_gamma_stat_C_{!s}_bin_{!s}" />'''.format(cut, ibin, cut, ibin, cut, ibin)+"\n")

    if doSyst and useSys:
        for key in SYSNAMES:
            R_Histo_Sys = None
            
            R_Histo_Sys = RooStats.HistFactory.HistoSys("ATLAS_"+key)
            xmlfile.write('''    <ExtSyst ConstrName="alpha_ATLAS_{!s}Constraint" NPName="alpha_ATLAS_{!s}" GOName="nom_alpha_ATLAS_{!s}" />'''.format(key, key, key)+"\n")

            signal_hist_sys_up = None
            if process == "signal":
                signal_hist_sys_up=deepcopy(SAMPLES["VBF_"+key+"__1up"].histograms["mBB_short_"+cut])
                signal_hist_sys_up.Add     (SAMPLES["ggF_"+key+"__1up"].histograms["mBB_short_"+cut])
            if process == "z":
                if useFullDistribution:
                    signal_hist_sys_up=deepcopy(SAMPLES["Zbb_QCD_"+key+"__1up"].histograms["mBB_short"])
                    signal_hist_sys_up.Add     (SAMPLES["Zbb_EWK_"+key+"__1up"].histograms["mBB_short"])
                    signal_hist_sys_up.Scale( (SAMPLES["Zbb_QCD_"+key+"__1up"].histograms["mBB_short_"+cut].Integral()
                                                   +SAMPLES["Zbb_EWK_"+key+"__1up"].histograms["mBB_short_"+cut].Integral())/ signal_hist_sys_up.Integral() )
                else:
                    signal_hist_sys_up=deepcopy(SAMPLES["Zbb_QCD_"+key+"__1up"].histograms["mBB_short_"+cut])
                    signal_hist_sys_up.Add     (SAMPLES["Zbb_EWK_"+key+"__1up"].histograms["mBB_short_"+cut])

            if ( ("JER" not in key) and ("QG_trackEfficiency" not in key) and ("QG_trackFakes" not in key)):
                R_Histo_Sys.SetHistoHigh(signal_hist_sys_up)
                
            if ( ("JER" in key) or ("QG_trackEfficiency" in key) or ("QG_trackFakes" in key)):
                JERDn = deepcopy(signal_hist)
                JERUp = deepcopy(signal_hist)

                diff = deepcopy(signal_hist)
                diff.Add(signal_hist_sys_up, -1)
                diff.Scale(0.5)
                for ibin in range(diff.GetNbinsX()):
                    diff.SetBinContent( ibin+1, abs(diff.GetBinContent(ibin+1)))

                JERUp.Add(diff, -1)
                JERDn.Add(diff)
                R_Histo_Sys.SetHistoHigh(JERUp)
                R_Histo_Sys.SetHistoLow(JERDn)
                R_Sample_cut.AddHistoSys(R_Histo_Sys)


#            elif "Eigen" in key:
#                rate = signal_hist_sys_up.Integral()/ signal_hist.Integral()
#                if rate != 1:
#                    xml = ''' <Systematic Name="alpha_ATLAS_btagging_{!s}" Constr="gaus" CentralValue="1" Mag="{!s}" WhereTo="yield"/> '''.format(key, str(abs(rate-1)))+"\n"
#                    btagging.write(xml)


            else:
                signal_hist_sys_dn = None
                if process == "signal":
                    signal_hist_sys_dn=deepcopy(SAMPLES["VBF_"+key+"__1down"].histograms["mBB_short_"+cut])
                    signal_hist_sys_dn.Add     (SAMPLES["ggF_"+key+"__1down"].histograms["mBB_short_"+cut])
                if process == "z":
                    if useFullDistribution:
                        signal_hist_sys_dn=deepcopy(SAMPLES["Zbb_QCD_"+key+"__1down"].histograms["mBB_short"])
                        signal_hist_sys_dn.Add     (SAMPLES["Zbb_EWK_"+key+"__1down"].histograms["mBB_short"])
                        signal_hist_sys_dn.Scale( (SAMPLES["Zbb_QCD_"+key+"__1down"].histograms["mBB_short_"+cut].Integral()
                                                   +SAMPLES["Zbb_EWK_"+key+"__1down"].histograms["mBB_short_"+cut].Integral())/ signal_hist_sys_up.Integral() )

                    else:
                        signal_hist_sys_dn=deepcopy(SAMPLES["Zbb_QCD_"+key+"__1down"].histograms["mBB_short_"+cut])
                        signal_hist_sys_dn.Add     (SAMPLES["Zbb_EWK_"+key+"__1down"].histograms["mBB_short_"+cut])

                R_Histo_Sys.SetHistoLow(signal_hist_sys_dn)
                R_Sample_cut.AddHistoSys(R_Histo_Sys)

    ### Initialize and save work space for signal
    R_Channel_cut.AddSample(R_Sample_cut);
    meas_sig.AddChannel(R_Channel_cut)
    RooStats.HistFactory.MakeModelAndMeasurementFast(meas_sig)

    WSname = "WS_combined_"+process+"_"+o.channel+"_"+cut+"_model.root"
    WSname_nosys = "WS_combined_"+process+"_"+o.channel+"_"+cut+"_model_nosys.root"

    if useSys:
        os.system("mv "+WSname + " FitFiles_"+o.channel +"/.")
    else:
        os.system("mv "+WSname + " FitFiles_"+o.channel +"/" + WSname_nosys)

    os.system("rm WS_*"+o.channel+"*")

    xmlfile.write('''</Model>''')    



def StandaloneFit(data_hist, signal_hist, z_hist, linfit, norm, fitname, fit_function, cut):

    if "spurious" in cut:
        signal_hist = deepcopy(signal_hist)
        signal_hist.Scale(norm)

    output.cd()
    fit_results_txt = open("FitFiles_"+o.channel+"/f_"+fitname+"_"+cut+"_"+o.channel, "a")
    fit_results_txt.write( "S "+str( signal_hist.Integral()/norm  )+"\n")
    fit_results_txt.write( "norm "+str( norm  )+"\n")
    fit_results_txt.write( "S Scaled "+str( signal_hist.Integral()  )+"\n")
    if linfit != None:
        fit_results_txt.write( "p0 "+str( linfit.GetParameter(0)  )+"\n")
        fit_results_txt.write( "p1 "+str( linfit.GetParameter(1)  )+"\n")

    alpha_z = 0
    alpha_sig = 0

    if "spurious" in cut:
        alpha_z, alpha_sig = bkgfit(data_hist, fit_function, fitname, doFloatZ = False , signal_hist = signal_hist, z_hist = z_hist)

    else:
        alpha_z, alpha_sig = bkgfit(data_hist, fit_function, fitname, doFloatZ = False , signal_hist = None, z_hist = z_hist)

    print "FITTED alpha_z", cut,  fitname, alpha_z
    print "FITTED alpha_sig", cut, fitname, alpha_sig

    if z_hist != None:

        z_plus_bkg = CopyHist(z_hist)
        z_only = CopyHist(z_hist)
        z_only.Scale(alpha_z)
        z_plus_bkg.Scale(alpha_z)
        z_plus_bkg.Add(fit_function)
        DrawTool_Diff.DrawHists("Bkg_Z_Fit_"+fitname+"_"+cut, ["M(bb)", "Number of Events"], 
                           [  z_only, deepcopy(data_hist), deepcopy(z_plus_bkg)],
                           [  "Z", "data", "Z + fit "+ fitname])
        
        fit_results_txt.write( "Z "+str( z_hist.Integral())+"\n")


    if "spurious" in cut:
        signal_only = deepcopy(signal_hist)
        signal_plus_bkg = deepcopy(signal_hist)
        signal_only.Scale(alpha_sig)
        signal_plus_bkg.Scale(alpha_sig)
        signal_plus_bkg.Add(fit_function)
        DrawTool_Diff.DrawHists("Spurious_Fit_"+fitname+"_"+cut, ["M(bb)", "Number of Events"], 
                           [ signal_only, deepcopy(data_hist), deepcopy(signal_plus_bkg)],
                           [ "Spurious", "data", "Spurious fit "+ fitname])

    ## Write out hybrid pseudo data && Create Fit Workspace
    ndof = len(mBB_Binning)-len(mBB_sig_Binning)-fit_function.GetNumberFreeParameters()+1
    if "spurious" in cut:
        ndof = len(mBB_Binning)-fit_function.GetNumberFreeParameters()+1
        
    FITRESULTS[cut][fitname] = FitResults(fitname, cut, o.channel, fit_function, ndof)

    ##############################################
    ### make hybrid pseudo data in text files  ###
    ##############################################
    if "spurious" in cut:
        fit_results_txt.write( "B "+str( data_hist.Integral()/norm  )+"\n")
        fit_results_txt.write( "B Scaled "+str( data_hist.Integral()  )+"\n")
        return

    ## predict fit over full regime
    fit_mBB = []
    for ibin in range(len(mBB_Binning)):
        bincenter = mBB_Binning[ibin]+mBB_bin_width/2.
        fit_mBB.append(int(fit_function.Eval(bincenter)))
    fit_results_txt.write( "B "+str( sum(fit_mBB)  )+"\n")

    pseudo_data_mBB = np.array([])
    pseudo_data_hist = TH1D("pseudo_data_"+fitname+"_"+cut, "pseudo_data_"+fitname+"_"+cut, len(mBB_Binning)-1, mBB_Binning)
    pseudo_bkg_hist  = TH1D("pseudo_bkg_"+fitname+"_"+cut,  "pseudo_bkg_"+fitname+"_"+cut,  len(mBB_Binning)-1, mBB_Binning)

    for ibin in range(len(mBB_sig_Binning)):
        bincenter = mBB_sig_Binning[ibin]+mBB_bin_width/2.
        pseudo_data_mBB = np.concatenate( (bincenter*np.ones(int(fit_function.Eval(bincenter))), pseudo_data_mBB) )

    rnp.fill_hist( pseudo_bkg_hist, pseudo_data_mBB)
    rnp.fill_hist( pseudo_data_hist, pseudo_data_mBB)
    pseudo_bkg_hist.Add(SAMPLES["data_2tag"].histograms["mBB_short_"+cut])
    pseudo_data_hist.Add(SAMPLES["data_2tag"].histograms["mBB_short_"+cut])
    pseudo_data_hist.Add(SAMPLES["VBF"].histograms["mBB_short_"+cut])
    pseudo_data_hist.Add(SAMPLES["ggF"].histograms["mBB_short_"+cut])

    output.cd()
    pseudo_data_hist.Write()
    SAMPLES["data_2tag"].histograms["pseudo_bkg_"+fitname+"_"+cut] = pseudo_bkg_hist

    out_data_array   = pseudo_data_mBB
    out_data_array = np.concatenate( (out_data_array, SAMPLES["data_2tag"].var["mBB"]))
    out_data_array = np.concatenate( (out_data_array, SAMPLES["VBF"].var["mBB"]))
    out_data_array = np.concatenate( (out_data_array, SAMPLES["ggF"].var["mBB"]))
    
    out_weight_array = np.ones(len(pseudo_data_mBB))
    out_weight_array = np.concatenate( (out_weight_array, SAMPLES["data_2tag"].var["eventWeight"]*SAMPLES["data_2tag"].cut[cut].array))
    out_weight_array = np.concatenate( (out_weight_array, SAMPLES["VBF"].var["eventWeight"]*SAMPLES["VBF"].cut[cut].array))
    out_weight_array = np.concatenate( (out_weight_array, SAMPLES["ggF"].var["eventWeight"]*SAMPLES["ggF"].cut[cut].array))
    
    ## apply range cut
    out_weight_array  = out_weight_array* (out_data_array<mBB_Binning[-1])
    out_weight_array  = out_weight_array* (out_data_array>mBB_Binning[0])
    out_data_array = out_data_array[out_weight_array>0]
    out_weight_array = out_weight_array[out_weight_array>0]

    np.savetxt("FitFiles_"+o.channel+"/f_pseudo_data_"+o.channel+"_fit_"+fitname+"_"+cut, out_data_array, fmt='%1.6f')
    np.savetxt("FitFiles_"+o.channel+"/f_pseudo_data_"+o.channel+"_fit_"+fitname+"_weight_"+cut, out_weight_array, fmt='%1.6f')


def ConvertHistToTxt(histogram, datatype, cut, postfix = ""):
    out_data_array = np.array([])
    out_weight_array = np.array([])
    data_list = []
    weight_list = []
    
    for ibin in range(histogram.GetNbinsX()):
        bincenter = histogram.GetBinCenter(ibin+1)
        bincontent = histogram.GetBinContent(ibin+1)
        data_list.append(bincenter)
        weight_list.append(bincontent)

    out_data_array   = np.array( data_list)
    out_weight_array = np.array( weight_list)


    if postfix == "":
        np.savetxt("FitFiles_"+o.channel+"/f_"+datatype+"_data_"+o.channel+"_CR_"+o.CR+"_"+cut+postfix, out_data_array, fmt='%1.6f')
        np.savetxt("FitFiles_"+o.channel+"/f_"+datatype+"_data_"+o.channel+"_CR_"+o.CR+"_weight_"+cut+postfix, out_weight_array, fmt='%1.6f')
    else:
        np.savetxt("FitFiles_"+o.channel+"/Poisson/f_"+datatype+"_data_"+o.channel+"_CR_"+o.CR+"_"+cut+postfix, out_data_array, fmt='%1.6f')
        np.savetxt("FitFiles_"+o.channel+"/Poisson/f_"+datatype+"_data_"+o.channel+"_CR_"+o.CR+"_weight_"+cut+postfix, out_weight_array, fmt='%1.6f')


def FitmBB(cuts, WritePseudoData = False):

    for cut in cuts:
        output.cd()
        print "cut ", cut

        ## fit spurious data
        FITRESULTS["spurious_"+cut]= {}

        norm_CR = None
        data_CR_fitband = None
        if o.CR == "0tag":
            norm_CR = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short_fitband"])
            data_CR_fitband = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short_fitband"])
        if o.CR == "2tag":
            norm_CR = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_noblind_fitband_CR"])
            data_CR_fitband = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_noblind_fitband_CR"])
        if o.CR == "mix":
            norm_CR = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short_fitband"])
            norm_CR.Add       (SAMPLES["data_2tag"].histograms["mBB_short_noblind_fitband_CR"])

            data_CR_fitband = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short_fitband"])
            data_CR_fitband.Add       (SAMPLES["data_2tag"].histograms["mBB_short_noblind_fitband_CR"])

        n_CR_fitband= norm_CR.Integral()
        norm_CR.Scale(1.0/n_CR_fitband)

        norm_SR = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_fitband_"+cut])
        n_SR_fitband = norm_SR.Integral()
        norm_SR.Scale(1.0/n_SR_fitband)

        ratio = deepcopy(norm_SR)
        ratio.Divide(norm_CR)
        ratio.Fit("pol1")
        linfit = ratio.GetFunction("pol1")
        linfit.Write()

        data_CR_fitband.Multiply(linfit)
        n_CR_fitband_postcor= data_CR_fitband.Integral()
        
        data_hist = None
        if o.CR == "0tag":
            data_hist = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short"])
        if o.CR == "2tag":
            data_hist = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_noblind_CR"])
        if o.CR == "mix":
            data_hist = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short"])
            data_hist.Add(SAMPLES["data_2tag"].histograms["mBB_short_noblind_CR"])

        ### Make CR data
        ConvertHistToTxt(data_hist, "CR", cut)
        #for i in range(1000):
        #    poissonhist = PoisonFluctuate(data_hist)
        #    ConvertHistToTxt(poissonhist, "CR", cut, str(i))
            

        ### Make spurious data
        data_hist.Multiply(linfit)
        data_norm = n_CR_fitband/n_CR_fitband_postcor ## normalize data to pre-correction
        mc_norm = n_CR_fitband/n_SR_fitband ## normalize mc to acceptance in CR
        data_hist.Scale(data_norm)

        ConvertHistToTxt(data_hist, "spurious", cut)


        ### Make Pseudo data from CR 
        signal_hist = CopyHist(SAMPLES["VBF"].histograms["mBB_short_"+cut])
        signal_hist.Add(SAMPLES["ggF"].histograms["mBB_short_"+cut])

        pseudo_data_hist = deepcopy(data_hist)
        pseudo_data_hist.Scale(1.0/mc_norm)
        pseudo_data_hist.Add(signal_hist)

        ConvertHistToTxt(pseudo_data_hist, "pseudo", cut)


        ### Make Pseudo data from 0-tag
        bkg_fitband = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_"+cut])
        bkgnorm = bkg_fitband.Integral()
        cr_fitband = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_noblind_fitband_CR"])
        crnorm = cr_fitband.Integral()
        bkg_fitband.Scale(1.0/bkgnorm)
        cr_fitband.Scale(1.0/crnorm)

        ratio_fitband = deepcopy(bkg_fitband)
        ratio_fitband.Divide(cr_fitband)
        ratio_fitband.Fit("pol1")
        linfit_0tag = ratio_fitband.GetFunction("pol1")

        bkg_fitband.Scale(bkgnorm)
        cr_fitband.Scale(crnorm)

        cr_fitband.Multiply(linfit_0tag)
        crnorm_postcor = cr_fitband.Integral()

        pseudo_data_hist_0tag = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_noblind_signal_CR"])
        pseudo_data_hist_0tag.Multiply(linfit_0tag)
        pseudo_data_hist_0tag.Scale(bkgnorm/crnorm_postcor)

        pseudo_data_hist_0tag.Add(SAMPLES["data_2tag"].histograms["mBB_short_"+cut])
        pseudo_data_hist_0tag.Add(SAMPLES["VBF"].histograms["mBB_short_"+cut])
        pseudo_data_hist_0tag.Add(SAMPLES["ggF"].histograms["mBB_short_"+cut])

        pseudo_data_hist_0tag.Write()
        SAMPLES["data_2tag"].histograms["mBB_short_"+cut].Write()

        ConvertHistToTxt(pseudo_data_hist_0tag, "0tag-pseudo", cut)

        z_hist = CopyHist(SAMPLES["Zbb_QCD"].histograms["mBB_short_"+cut])
        z_hist.Add(SAMPLES["Zbb_EWK"].histograms["mBB_short_"+cut])

        StandaloneFit( data_hist, signal_hist, z_hist, linfit, mc_norm, "BernsteinO2", BernsteinO2, "spurious_"+cut)
        StandaloneFit( data_hist, signal_hist, z_hist, linfit, mc_norm, "BernsteinO3", BernsteinO3, "spurious_"+cut)
        StandaloneFit( data_hist, signal_hist, z_hist, linfit, mc_norm, "BernsteinO4", BernsteinO4, "spurious_"+cut)
        StandaloneFit( data_hist, signal_hist, z_hist, linfit, mc_norm, "BernsteinO5", BernsteinO5, "spurious_"+cut)

        ### background only fit 2tag data
        FITRESULTS[cut]= {}
        z_hist_fitband = CopyHist(SAMPLES["Zbb_QCD"].histograms["mBB_short_fitband_"+cut])
        z_hist_fitband.Add(SAMPLES["Zbb_EWK"].histograms["mBB_short_fitband_"+cut])
        
        data_hist = CopyHist(SAMPLES["data_2tag"].histograms["mBB_short_"+cut])
        StandaloneFit( data_hist, signal_hist, z_hist_fitband, None, 1, "BernsteinO2", BernsteinO2, cut)
        StandaloneFit( data_hist, signal_hist, z_hist_fitband, None, 1, "BernsteinO3", BernsteinO3, cut)
        StandaloneFit( data_hist, signal_hist, z_hist_fitband, None, 1, "BernsteinO4", BernsteinO4, cut)
        StandaloneFit( data_hist, signal_hist, z_hist_fitband, None, 1, "BernsteinO5", BernsteinO5, cut)

        GenerateSignalWS("signal", cut, False, True)
        GenerateSignalWS("signal", cut, False, False)
        GenerateSignalWS("z", cut, True, False)
        GenerateSignalWS("z", cut, True, True)


        ### build z injection test data 
        
        z_hist_full = deepcopy(  SAMPLES["Zbb_QCD"].histograms["mBB_short"] )
        z_hist_full.Add       (  SAMPLES["Zbb_EWK"].histograms["mBB_short"] )
        z_hist_full.Scale( z_hist.Integral() / z_hist_full.Integral())

        for strength in [0.5, 1, 2, 5, 10]:
            for fitname in ["BernsteinO2", "BernsteinO3", "BernsteinO4", "BernsteinO5"]:
                zinject_data_hist = deepcopy(empty_mBB_hist)
                zinject_data_hist.Add( FITRESULTS["CR"][fitname].fit)
                zinject_data_hist.Multiply(linfit_0tag)
                zinject_data_hist.Scale(bkgnorm/crnorm_postcor)
                zinject_data_hist.Add(z_hist_full, strength)
                if cut != "CR":
                    zinject_data_hist.Add(signal_hist)

                print "z inject data, signal norm", signal_hist.Integral()
                print "z inject data, z norm", z_hist_full.Integral()

                ConvertHistToTxt(zinject_data_hist, "Z_inject_"+str(strength)+"_"+fitname, cut)

        ### Unblind data
        ConvertHistToTxt(SAMPLES["data_2tag"].histograms["mBB_short_noblind_"+cut], "unblind", cut)
        output.cd()


def DrawMBBCorr():
    Make1DPlots(SAMPLES, "mBB_corrcheck",        "mBB",             "eventWeight", [], mBB_Binning_sig)

    def weighted_avg_and_std(values, weights):
        """
        Return the weighted average and standard deviation.
        
        values, weights -- Numpy ndarrays with the same shape.
        """
        average = np.average(values, weights=weights)
        variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
        return (average, sqrt(variance))

    mean_nocorr, var_nocorr=  weighted_avg_and_std(SAMPLES["VBF_NoCorr"].var["mBB"], weights=SAMPLES["VBF"].var["eventWeight"]*(SAMPLES["VBF_NoCorr"].var["mBB"]>40)*(SAMPLES["VBF_NoCorr"].var["mBB"]<150))
    mean_onemu, var_onemu =  weighted_avg_and_std(SAMPLES["VBF_OneMu"].var["mBB"], weights=SAMPLES["VBF_OneMu"].var["eventWeight"]*(SAMPLES["VBF_OneMu"].var["mBB"]>40)*(SAMPLES["VBF_OneMu"].var["mBB"]<150))
    mean_pt, var_pt =  weighted_avg_and_std(SAMPLES["VBF"].var["mBB"], weights=SAMPLES["VBF"].var["eventWeight"]*(SAMPLES["VBF"].var["mBB"]>40)*(SAMPLES["VBF"].var["mBB"]<150))

    DrawTool.AddTexts(0.18, 0.7, 1, 0.03, "No Correction: Mean="+ str(round( mean_nocorr,2))+" RMS="+str(round( var_nocorr,2)))
    DrawTool.AddTexts(0.18, 0.65, 1, 0.03, "Muon Correction: Mean="+ str(round( mean_onemu,2)) +" RMS="+str(round( var_onemu,2)))
    DrawTool.AddTexts(0.18, 0.6, 1, 0.03, "Muon + p_{T} Reco Correction: Mean="+ str(round( mean_pt,2))  +" RMS="+str(round( var_pt,2)))

    DrawTool.lumi = "0"
    SAMPLES["VBF_NoCorr"].histograms["mBB_corrcheck"].Scale( 1.0/SAMPLES["VBF_NoCorr"].histograms["mBB_corrcheck"].Integral())
    SAMPLES["VBF_OneMu"].histograms["mBB_corrcheck"].Scale(  1.0/SAMPLES["VBF_OneMu"].histograms["mBB_corrcheck"].Integral())
    SAMPLES["VBF"].histograms["mBB_corrcheck"].Scale( 1.0/SAMPLES["VBF"].histograms["mBB_corrcheck"].Integral())
    DrawTool.DrawHists("Mbb_Corrcheck", ["M(bb)", "A.U."], 
                       [ SAMPLES["VBF_NoCorr"].histograms["mBB_corrcheck"], SAMPLES["VBF_OneMu"].histograms["mBB_corrcheck"], SAMPLES["VBF"].histograms["mBB_corrcheck"]],
                       [ "Nominal", "One Mu", "One Mu + PtReco"])


    print "lower sideband mbb fraction", SAMPLES["VBF"].histograms["mBB_corrcheck"].Integral(1,30)
    print "uppder sideband mbb fraction", SAMPLES["VBF"].histograms["mBB_corrcheck"].Integral(50, 55)

    DrawTool.ClearTexts()


    Zbb_array_nocorr = np.concatenate( (SAMPLES["Zbb_QCD_NoCorr"].var["mBB"], SAMPLES["Zbb_EWK_NoCorr"].var["mBB"]))
    Zbb_array_onemu = np.concatenate( (SAMPLES["Zbb_QCD_OneMu"].var["mBB"], SAMPLES["Zbb_EWK_OneMu"].var["mBB"]))
    Zbb_array_pt = np.concatenate( (SAMPLES["Zbb_QCD"].var["mBB"], SAMPLES["Zbb_EWK"].var["mBB"]))
    evt_array = np.concatenate( (SAMPLES["Zbb_QCD"].var["eventWeight"], SAMPLES["Zbb_EWK"].var["eventWeight"]))

    mean_nocorr, var_nocorr=  weighted_avg_and_std( Zbb_array_nocorr, weights=evt_array*( (Zbb_array_nocorr>40)*(Zbb_array_nocorr<150) ) )
    mean_onemu, var_onemu =  weighted_avg_and_std( Zbb_array_onemu, weights=evt_array*( (Zbb_array_onemu>40)*(Zbb_array_onemu<150) ) )
    mean_pt, var_pt =  weighted_avg_and_std(  Zbb_array_pt, weights=evt_array*((Zbb_array_pt>40)*(Zbb_array_pt<150) ) )

    DrawTool.AddTexts(0.18, 0.7, 1, 0.03, "No Correction: Mean="+ str(round( mean_nocorr,2))+" RMS="+str(round( var_nocorr,2)))
    DrawTool.AddTexts(0.18, 0.65, 1, 0.03, "Muon Correction: Mean="+ str(round( mean_onemu,2)) +" RMS="+str(round( var_onemu,2)))
    DrawTool.AddTexts(0.18, 0.6, 1, 0.03, "Muon + p_{T} Reco Correction: Mean="+ str(round( mean_pt,2))  +" RMS="+str(round( var_pt,2)))


    Zbb_nocorr = deepcopy(SAMPLES["Zbb_QCD_NoCorr"].histograms["mBB_corrcheck"])
    Zbb_nocorr.Add(SAMPLES["Zbb_EWK_NoCorr"].histograms["mBB_corrcheck"])
    Zbb_nocorr.Scale( 1.0/Zbb_nocorr.Integral() )

    Zbb_OneMu = deepcopy(SAMPLES["Zbb_QCD_OneMu"].histograms["mBB_corrcheck"])
    Zbb_OneMu.Add(SAMPLES["Zbb_EWK_OneMu"].histograms["mBB_corrcheck"])
    Zbb_OneMu.Scale( 1.0/Zbb_OneMu.Integral())

    Zbb_pt = deepcopy(SAMPLES["Zbb_QCD"].histograms["mBB_corrcheck"])
    Zbb_pt.Add(SAMPLES["Zbb_EWK"].histograms["mBB_corrcheck"])
    Zbb_pt.Scale( 1.0/Zbb_pt.Integral())

    DrawTool.DrawHists("Mbb_Corrcheck_Z", ["M(bb)", "A.U."], [ Zbb_nocorr, Zbb_OneMu, Zbb_pt], [ "Nominal", "One Mu", "One Mu + PtReco"])


###### main code########            
wp = o.btag

btagcuts = ["2b_mv2c10_70", "0b_mv2c10_85"]
pTBBcuts = ["pTBB>70", "pTBB>100", "pTBB>130", "pTBB>150",  "pTBB>160", "pTBB>170", "pTBB>190", "pTBB>210"]
mJJcuts  = ["mJJ>50", "mJJ>100", "mJJ>150", "mJJ>200"]
sidebandcuts = ["sideband_up", "sideband_low", "signal"]
fitbandcuts = ["fitband_up", "fitband_low", "fitband"]
BDTcuts = ["CR", "SRII", "SRI"]
BDTcuts_More = ["SRI", "SRII", "CR", "PassBDTCut00-01","PassBDTCut00-015","PassBDTCut00-02", "PassBDTCut00-03","PassBDTCut00-05","PassBDTCut00-10","PassBDTCut00-20", "PassBDTCut20-40","PassBDTCut40-60","PassBDTCut60-80","PassBDTCut80-100"]

if o.NBDTReg== "3":
    BDTcuts = ["CR", "SRI", "SRII", "SRIII"]
    BDTcuts_More = ["SRI", "SRII", "SRIII", "CR", "PassBDTCut00-01","PassBDTCut00-015","PassBDTCut00-02", "PassBDTCut00-03","PassBDTCut00-05","PassBDTCut00-10","PassBDTCut00-20", "PassBDTCut20-40","PassBDTCut40-60","PassBDTCut60-80","PassBDTCut80-100"]

MVAVarList_Use= ["mJJ", "pTJJ", "cosTheta_boost","mindRJ1_Ex", "mindRJ2_Ex", "max_J1J2", "eta_J_star", "NTrk500PVJ1", "NTrk500PVJ2", "pT_balance"]

if o.channel == "2cen" and o.postfix != "looseVBF":
    #MVAVarList_Use= ["mJJ", "pTJJ", "dEtaJJ", "cosTheta_boost","mindRJ1_Ex", "mindRJ2_Ex", "max_J1J2", "eta_J_star", "NTrk500PVJ2", "pT_balance"]
    MVAVarList_Use= ["mJJ", "pTJJ", "cosTheta_boost","mindRJ1_Ex", "mindRJ2_Ex", "max_J1J2", "eta_J_star", "NTrk500PVJ2", "pT_balance"]

if o.MVALeaveOut != "No":
    ind = MVAVarList_Use.index(o.MVALeaveOut)
    del MVAVarList_Use[ind]

############################
####    Main Code       ####
############################


################################
#### Initialize and Apply Cut###
################################

InitSamples()

## creating b-tagging systematics
CreateWeightSample("VBF")
CreateWeightSample("ggF")
CreateWeightSample("Zbb_QCD")
CreateWeightSample("Zbb_EWK")

SYSLIST = SYSLIST + BTAGSYSLIST + THSYSLIST
SYSNAMES = SYSNAMES + BTAGSYSNAMES + THSYSNAMES

if o.ApplypTBBCut:
    ApplypTBBCut()

filename = "out_"+ o.channel+"_" + o.period + "_" +o.ApplypTBBCut*"_pTBBCut"+o.postfix+".root"
if o.MVALeaveOut != "No":
    filename = "out_"+ o.channel+"_" + o.period + "_" +o.ApplypTBBCut*"_pTBBCut"+o.postfix+"_LeaveOut_"+o.MVALeaveOut+".root"
output = TFile(filename, "recreate")
output.cd()

ScaleLumi()
BlindData()
ApplyBtag()

if o.postfix == "4cencomp" or o.postfix == "NLOComp":
    DrawMVAInputComparisonDifferentSamples(MVAVarList_Use)
    sys.exit(0)
if o.postfix == "mbbcorrcheck":
    DrawMBBCorr()
    sys.exit(0)

ComputeMVAWeights(MVAVarList_Use)

#Make1DPlots(SAMPLES, "mBB_no_corr","mBB_no_corr",     "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
#Make1DPlots(SAMPLES, "mBBsig_no_corr",  "mBB_no_corr","eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), mBB_Binning_sig)

Make1DPlots(SAMPLES, "mBBsig",     "mBB",             "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More + fitbandcuts), mBB_Binning_sig)
Make1DPlots(SAMPLES, "mBB",        "mBB",             "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More + fitbandcuts), JetpT_Binning,
            precuts = ["sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "mBB_short",  "mBB",             "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts + BDTcuts_More + fitbandcuts), mBB_Binning,
            precuts = ["fitband"])
Make1DPlots(SAMPLES, "mBB_short_noblind",  "mBB",     "eventWeight_noblind", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More + fitbandcuts), mBB_Binning,
            precuts = ["fitband", "signal"])

Make1DPlots(SAMPLES, "NJets",      "nJets",           "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), NJet_Binning)
Make1DPlots(SAMPLES, "pTJ1",       "pTJ1",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
Make1DPlots(SAMPLES, "pTJ2",       "pTJ2",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
#Make1DPlots(SAMPLES, "pTB1",       "pTB1",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
#Make1DPlots(SAMPLES, "pTB2",       "pTB2",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
Make1DPlots(SAMPLES, "etaJ1",      "etaJ1",           "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), Eta_Binning)
Make1DPlots(SAMPLES, "etaJ2",      "etaJ2",           "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), Eta_Binning)
Make1DPlots(SAMPLES, "pTBB",       "pTBB",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)

Make1DPlots(SAMPLES, "mJJ",        "mJJ",             "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), MJJ_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "pTJJ",       "pTJJ",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), pTJJ_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "NTrk500PVJ1","NTrk500PVJ1",     "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), NTrk_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "NTrk500PVJ2","NTrk500PVJ2",     "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), NTrk_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "pTB1",       "pTB1",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "pTB2",       "pTB2",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "cosTheta_boost",   "cosTheta_boost",        "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), CosTheta_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "mindRJ1_Ex", "mindRJ1_Ex",      "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), dR_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "mindRJ2_Ex", "mindRJ2_Ex",      "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), dR_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "max_J1J2",   "max_J1J2",        "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), MaxJeteta_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "eta_J_star", "eta_J_star",      "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), Jeteta_Binning,
            precuts = [ "sideband_up", "sideband_low"])
Make1DPlots(SAMPLES, "pT_balance", "pT_balance",    "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), ptballance_Binning,
            precuts = [ "sideband_up", "sideband_low"])

#DrawPlotsWithCutsForEachSample("mBB",  "BDTcuts", "M(bb)",  "Events", BDTcuts_More, dosys = False)
#DrawPlotsWithCutsForEachSample("pTBB", "pTBBCuts", "M(bb)",  "Events", pTBBcuts, dosys = False)

DrawPlotsWithCutsForEachSample("mBB",  "pTBBCuts", "M(bb)",  "Events", pTBBcuts, dosys = False)

DrawPlotsPreCutMultiSample("NJets", "NoCut",   "NJets",  "A.U.")
DrawPlotsPreCutMultiSample("pTJ1",  "NoCut",   "p_{T} J1",   "A.U.")
DrawPlotsPreCutMultiSample("pTJ2",  "NoCut",   "p_{T} J2",   "A.U.")
DrawPlotsPreCutMultiSample("etaJ1", "NoCut",   "#eta J1",  "A.U.")
DrawPlotsPreCutMultiSample("etaJ2", "NoCut",   "#eta J2",  "A.U.")
DrawPlotsPreCutMultiSample("pTBB",  "NoCut",    "p_{T}(bb)",  "A.U.")
DrawPlotsPreCutMultiSample("pTB1",  "NoCut",    "p_{T} B1",  "A.U.")
DrawPlotsPreCutMultiSample("pTB2",  "NoCut",    "p_{T} B2",  "A.U.")

if not o.ApplypTBBCut:
    sys.exit(0)

DrawMVAInputComparison(["mBB"]+MVAVarList_Use, BDTcuts)
CalculateSensitivity(BDTcuts_More)


#############################
####  Customized Plots  #####
#############################

#if o.channel== "2cen":
#    Make1DPlots(SAMPLES, "mBB_short",  "mBB",  "eventWeight", ["SecondBNotTrue"], mBB_Binning)
#    
#    DrawTool.DrawHists("mBB_Full_Vs_NotTrueB",  ["M(bb)", "Events"], 
#                       [SAMPLES["VBF"].histograms["mBB_short"], SAMPLES["VBF"].histograms["mBB_short_SecondBNotTrue"]],
#                       ["Full", "Second B Not True"])

DrawTool_Ratio.DrawHists("mBB_0tag_vs_2tag",  ["M(bb)", "Events"], 
                         [SAMPLES["data_0tag"].histograms["mBB_short"], SAMPLES["data_2tag"].histograms["mBB_short" ]],
                         ["data_0tag", "data_2tag"])


#### Draw M(bb) for 0-tag vs 2-tag with BDT cut
#### Draw M(bb) for 2-tag BDT cuts comparison
for icut in range(len(BDTcuts)):
    cut1 = BDTcuts[icut]

    data_0tag = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short_fitband"])
    data_0tag.Scale(1./data_0tag.Integral())

    data_2tag = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_"+cut1 ])
    data_2tag.Scale(1./data_2tag.Integral())

    DrawTool_Ratio.DrawHists("mBB_0tag_vs_2tag_"+cut1,  ["M(bb)", "Events"], 
                             [data_0tag, data_2tag],
                             ["data_0tag", "data_2tag_"+cut1])


    mBB_Z_nocut = deepcopy(SAMPLES["Zbb_EWK"].histograms["mBB_short"])
    mBB_Z_nocut.Add(SAMPLES["Zbb_QCD"].histograms["mBB_short"])

    mBB_Z_BDT = deepcopy(SAMPLES["Zbb_EWK"].histograms["mBB_short_"+cut1])
    mBB_Z_BDT.Add(SAMPLES["Zbb_QCD"].histograms["mBB_short_"+cut1])

    mBB_Z_nocut.Scale(1.0/mBB_Z_nocut.Integral())
    mBB_Z_BDT.Scale(1.0/mBB_Z_BDT.Integral())

    DrawTool_Ratio.DrawHists("Zbb_mBB_vs_"+cut1,  ["M(bb)", "Events"], 
                             [mBB_Z_BDT, mBB_Z_nocut], 
                             ["Zbb "+cut1, "Zbb Full"])

    for jcut in range(icut+1, len(BDTcuts)):
        cut2 = BDTcuts[jcut]

        tmp_cut1 = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_"+cut1])
        tmp_cut1.Add(SAMPLES["Zbb_QCD"].histograms["mBB_short_fitband_"+cut1], -1)
        tmp_cut1.Add(SAMPLES["Zbb_EWK"].histograms["mBB_short_fitband_"+cut1], -1)

        tmp_cut2 = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_"+cut2])
        tmp_cut2.Add(SAMPLES["Zbb_QCD"].histograms["mBB_short_fitband_"+cut2], -1)
        tmp_cut2.Add(SAMPLES["Zbb_EWK"].histograms["mBB_short_fitband_"+cut2], -1)

        tmp_cut1.Scale( 1./tmp_cut1.Integral())
        tmp_cut2.Scale( 1./tmp_cut2.Integral())

        DrawTool_Ratio.DrawHists("mBB_2tag_"+cut1+"_vs_"+cut2,  ["M(bb)", "A.U."], 
                                 [tmp_cut1, tmp_cut2],
                                 ["data_2tag_"+cut1, "data_2tag_"+cut2])

#########################
####  Cut Study  ########
#########################

#DrawPlotsWithCutsForEachSample("mBB_no_corr",  "pTBBCuts", "M(bb)",  "Events", pTBBcuts)
#DrawPlotsWithCutsForEachSample("pTBB", "pTBBCuts", "pT(bb)", "Events", pTBBcuts)
#DrawPlotsWithCutsForEachSample("mJJ",  "pTBBCuts", "M(JJ)",  "Events", pTBBcuts)

#DrawPlotsWithCutsForEachSample("mBB",  "BDTCuts", "M(bb)",     "Events", BDTcuts)
#DrawPlotsWithCutsForEachSample("mJJ",  "BDTCuts", "M(JJ)",     "Events", BDTcuts)
#DrawPlotsWithCutsForEachSample("pTBB", "BDTCuts", "pT(bb)",    "Events", BDTcuts)

#DrawPlotsWithCutsForEachSample("mBB",  "BDTCuts_Normed", "M(bb)",  "Events", BDTcuts, norm=True)
#DrawPlotsWithCutsForEachSample("mJJ",  "BDTCuts_Normed", "M(JJ)",  "Events", BDTcuts, norm=True)
#DrawPlotsWithCutsForEachSample("pTBB", "BDTCuts_Normed", "pT(bb)", "Events", BDTcuts, norm=True)

#DrawPlotsWithCutsForEachSample("pTB1",  "BDTCuts_Normed", "p_{T} B1",  "Events", BDTcuts, norm=True)
#DrawPlotsWithCutsForEachSample("pTB2",  "BDTCuts_Normed", "p_{T} B2",  "Events", BDTcuts, norm=True)
#
#DrawPlotsWithCutsForEachSample("pTB1",  "BDTCuts", "p_{T} B1",  "Events", BDTcuts, norm=True)
#DrawPlotsWithCutsForEachSample("pTB2",  "BDTCuts", "p_{T} B2",  "Events", BDTcuts, norm=True)


#########################
####  Fit M(bb)  ########
#########################

FitmBB(BDTcuts, True)

for cut in  FITRESULTS:
    for fit in FITRESULTS[cut]:
        FITRESULTS[cut][fit].GetxmlForm()

    ##### Calculate the F-test for the ratio of the fits
    if cut == "2region":
        continue

    CalFtest( FITRESULTS[cut]["BernsteinO2"], FITRESULTS[cut]["BernsteinO3"] )
    CalFtest( FITRESULTS[cut]["BernsteinO3"], FITRESULTS[cut]["BernsteinO4"] )
    CalFtest( FITRESULTS[cut]["BernsteinO4"], FITRESULTS[cut]["BernsteinO5"] )

output.Close()
