# this import is required

import rootpy.tree
from rootpy.io import root_open
from ROOT import gDirectory, TLorentzVector, TVector3, TH1D, TH2D, TFile, RooStats, TMath
from HistoLab import HistoTool, CopyHist, StackHists, GetPlot, ScaleHist, Project2D, CopyHists, project1D, NormHists, AddHistPairs, ResetBins, GetHistMedian, GetHistIQR, GetRatios, GetAveRatios, GetGausMean, CalSysVariation, CalSysDiff
from FitFunctions import BernsteinO2, BernsteinO3, BernsteinO4, BernsteinO5, ExpoBernsteinO2, ExpoBernsteinO3, ExpoBernsteinO4, ExpoBernsteinO5, Expo, Expo2, Expo3
from FitFunctions import BernsteinO2_raw, BernsteinO3_raw, BernsteinO4_raw, BernsteinO5_raw, ExpoBernsteinO2_raw, ExpoBernsteinO3_raw, ExpoBernsteinO4_raw, ExpoBernsteinO5_raw, Expo_raw, Expo2_raw, Expo3_raw
from FitFunctions import BernsteinO2Lin, BernsteinO3Lin, BernsteinO4Lin
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

DrawTool = HistoTool()
DrawTool.lumi = "0"
DrawTool.sqrtS = "13"
DrawTool.SetCompData(False)
DrawTool.OutPlotDir = "Plots/"

DrawTool_Ratio = HistoTool()
DrawTool_Ratio.DrawRatio = "ratio"
DrawTool_Ratio.lumi = "0"
DrawTool_Ratio.sqrtS = "13"
DrawTool_Ratio.shift = 0.11
DrawTool_Ratio.SetCompData(True)

DrawTool_Diff = HistoTool()
DrawTool_Diff.lumi = "0"
DrawTool_Diff.sqrtS = "13"
DrawTool_Diff.shift = 0.11
DrawTool_Diff.doDiff = True
DrawTool_Diff.DrawRatio = "(Data-Fit) / Data"
DrawTool_Diff.SetCompData(True)


p = OptionParser()
p.add_option('--btag', type = "string", default = 'tight',   dest = 'btag', help = 'btag working point')
p.add_option('--pTBBCut', action="store_true", default = False, dest = 'ApplypTBBCut', help = 'cut on pT(bb)')
p.add_option('--mJJCut', action="store_true", default = False, dest = 'ApplymJJCut', help = 'cut on m(jj)')
p.add_option('--MVAList', type = "string", default = 'short',   dest = 'MVAList', help = 'MVA List')
p.add_option('--test', action="store_true", default = False, dest = 'test', help = 'test mode')
p.add_option('--fit', action="store_true", default = False, dest = 'fit', help = 'fit mode')
p.add_option('--channel', type="string", default = '4cen', dest = 'channel', help = 'channel')
p.add_option('--spurious', type="string", default = '2tag', dest = 'spurious', help = 'spurious type')
(o,a) = p.parse_args()

os.system("rm f_*"+o.channel+"*")
os.system("rm workspace_*"+o.channel+"*")
os.system("rm FitFiles_"+o.channel+"/*")

mBB_bin_width = 2

PDGID_Binning = array('d',range(16))
NJet_Binning = array('d',range(16))
JetpT_Binning = array('d',range(0, 400, 10))
MJJ_Binning = array('d',range(0, 2500, 10))
Jeteta_Binning = array('d',[x / 10. for x in range(-46, 46, 2)])
mBB_Binning = array('d', range(fit_start, fit_end+mBB_bin_width, mBB_bin_width))
mBB_Binning_sig = array('d', range(70, 150+mBB_bin_width, mBB_bin_width))
dR_Binning = array('d',[x / 10. for x in range(0, 80, 1)])
dR_Binning_short = array('d',[x / 10. for x in range(0, 30, 1)])
JetWidth_Binning = array('d',[x / 100. for x in range(0, 40, 1)])
CosTheta_Binning = array('d',[x / 10. for x in range(0, 11, 1)])
MV2_Binning = array('d', [x/100. for x in range(-100, 100, 1)])
BDT_Binning = array('d', [x/100. for x in range(80, 120, 1)])
NN_Binning = array('d', [x/100. for x in range(0, 100, 1)])

mBB_sig_Binning = array('d', range(100, 140, mBB_bin_width))

SYSLIST = ["JET_GroupedNP_1__1up", "JET_GroupedNP_1__1down", "JET_GroupedNP_2__1up", "JET_GroupedNP_2__1down", "JET_GroupedNP_3__1up", "JET_GroupedNP_3__1down",
           "JET_EtaIntercalibration_NonClosure__1up", "JET_EtaIntercalibration_NonClosure__1down", "JET_JER_SINGLE_NP__1up"]
SYSNAMES = ["JET_GroupedNP_1", "JET_GroupedNP_2", "JET_GroupedNP_3", "JET_EtaIntercalibration_NonClosure", "JET_JER_SINGLE_NP"]
SAMPLES = {}
SYSSAMPLES = {}
FITSAMPLES= {}
FITRESULTS = {}

def InitSamples():
    
    #VBF = PhysicsProcess("VBF H->bb", "input/VBF_Reader_10_12_4cen_NoTrigMatch_BMED.root")
    #VBF_HYBRID = PhysicsProcess("VBF H->bb PC HLTTrig", "input/VBF_Reader_10_12_4cen_NoTrigMatch_HYBRID.root")
    #VBF_BMED   = PhysicsProcess("VBF H->bb BMED HLTTrig", "input/VBF_Reader_10_12_4cen_NoTrigMatch_BMED.root")
    #VBF_BMED_TrigMatch = PhysicsProcess("VBF H->bb BMED TrigMatch", "input/VBF_Reader_10_12_4cen_TrigMatch_BMED.root")
    #VBF_BMED_TrigMatch_NoWeight = PhysicsProcess("VBF H->bb BMED TrigMatch NoWeight", "input/VBF_Reader_10_12_4cen_TrigMatch_BMED_NoWeight.root")

    VBF = None
    ggF = None
    Zbb_QCD = None
    Zbb_EWK = None
    data_2tag = None
    data_0tag = None

    if o.channel == "4cen":
        VBF = PhysicsProcess("VBF H->bb", "input/VBF_Reader_10_12_4cen.root")
        ggF = PhysicsProcess("ggF H->bb", "input/ggF_Reader_10_12_4cen.root")
        data_0tag = PhysicsProcess("data 0tag", "input/data_0tag_10_12_4cen.root")
        data_2tag = PhysicsProcess("data 2tag", "input/data_ICHEP_Reader_10_16_2tag_4cen.root")
        Zbb_QCD = PhysicsProcess("QCD Z->bb", "input/Zbb_QCD_Reader_10_12_4cen.root")
        Zbb_EWK = PhysicsProcess("EWK Z->bb", "input/Zbb_EWK_Reader_10_12_4cen.root")

    if o.channel == "2cen":
        #VBF = PhysicsProcess("VBF H->bb", "input/VBF_Reader_10_19_2cen_tight.root")
        #VBF = PhysicsProcess("VBF H->bb", "input/VBF_Reader_10_17_2cen_MediumBTAG.root")

        VBF = PhysicsProcess("VBF H->bb", "input/VBF_Reader_10_19_2cen.root")
        ggF = PhysicsProcess("ggF H->bb", "input/ggF_Reader_10_27_2cen.root")
        data_2tag = PhysicsProcess("data 2tag", "input/data_2tag_10_12_2cen.root")

        #VBF = PhysicsProcess("VBF H->bb", "input/VBF_Reader_11_02_2cen_BTAG_Medium.root")
        #ggF = PhysicsProcess("ggF H->bb", "input/ggF_Reader_11_02_2cen_BTAG_Medium.root")
        #data_2tag = PhysicsProcess("data 2tag", "input/data_ICHEP_Reader_11_01_2tag_2cen_BTAG_Medium.root")

        data_0tag = PhysicsProcess("data 0tag", "input/data_ICHEP_Reader_10_28_0tag_2cen.root")
        Zbb_QCD = PhysicsProcess("QCD Z->bb", "input/Zbb_QCD_Reader_10_27_2cen.root")
        Zbb_EWK = PhysicsProcess("EWK Z->bb", "input/Zbb_EWK_Reader_10_27_2cen.root")

    Zbb_QCD_0tag = PhysicsProcess("QCD Z->bb 0tag", "input/Zbb_QCD_Reader_10_26_0tag_4cen.root")
    Zbb_EWK_0tag = PhysicsProcess("EWK Z->bb 0tag", "input/Zbb_EWK_Reader_10_26_0tag_4cen.root")


    SAMPLES["VBF"] = VBF
    SAMPLES["ggF"] = ggF
    SAMPLES["Zbb_QCD"] = Zbb_QCD
    SAMPLES["Zbb_EWK"] = Zbb_EWK
    SAMPLES["Zbb_QCD_0tag"] = Zbb_QCD_0tag
    SAMPLES["Zbb_EWK_0tag"] = Zbb_EWK_0tag
    SAMPLES["data_0tag"] = data_0tag
    SAMPLES["data_2tag"] = data_2tag

#    SAMPLES["VBF_HYBRID"] = VBF_HYBRID
#    SAMPLES["VBF_BMED"] = VBF_BMED
#    SAMPLES["VBF_BMED_TrigMatch"] = VBF_BMED_TrigMatch
#    SAMPLES["VBF_BMED_TrigMatch_NoWeight"] = VBF_BMED_TrigMatch_NoWeight

#    for sys in SYSLIST:
#        SYSSAMPLES["VBF_"+sys] = PhysicsProcess("VBF H->bb "+sys, "VBF_09_16.root", sys)
#        SYSSAMPLES["ggF_"+sys] = PhysicsProcess("ggF H->bb "+sys, "ggF_09_16.root", sys)
#        SYSSAMPLES["Zbb_QCD_"+sys] = PhysicsProcess("QCD Z->bb "+sys, "Zbb_QCD_09_16.root", sys)
#        SYSSAMPLES["Zbb_EWK_"+sys] = PhysicsProcess("EWK Z->bb "+sys, "Zbb_EWK_09_16.root", sys)


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
        SAMPLES[s].AddEventVarFromTree("weightSysts", o.test)

        ### back up event weight used only for control region no blinding
        SAMPLES[s].var["eventWeight_noblind"]= SAMPLES[s].var["eventWeight"]

        ### Load event variable
        SAMPLES[s].AddEventVarFromTree("nJets", o.test)
        SAMPLES[s].AddEventVarFromTree("mBB", o.test)
        SAMPLES[s].AddEventVarFromTree("mBB_no_corr", o.test)
        SAMPLES[s].AddEventVarFromTree("pTBB", o.test)
        SAMPLES[s].AddEventVarFromTree("mJJ", o.test)
        SAMPLES[s].AddEventVarFromTree("pTJJ", o.test)
        SAMPLES[s].AddEventVarFromTree("dRBB", o.test)
        SAMPLES[s].AddEventVarFromTree("dEtaJJ", o.test)
        SAMPLES[s].AddEventVarFromTree("dEtaBB", o.test)
        SAMPLES[s].AddEventVarFromTree("WidthJ1", o.test)
        SAMPLES[s].AddEventVarFromTree("WidthJ2", o.test)
        SAMPLES[s].AddEventVarFromTree("etaJ1", o.test)
        SAMPLES[s].AddEventVarFromTree("etaJ2", o.test)
        SAMPLES[s].AddEventVarFromTree("etaB1", o.test)
        SAMPLES[s].AddEventVarFromTree("etaB2", o.test)
        SAMPLES[s].AddEventVarFromTree("cosTheta_MVA", o.test)
        SAMPLES[s].AddEventVarFromTree("HT_soft", o.test)
        SAMPLES[s].AddEventVarFromTree("HT_MVA", o.test)
        SAMPLES[s].AddEventVarFromTree("BDT", o.test)
        SAMPLES[s].AddEventVarFromTree("whoIsB1", o.test)
        SAMPLES[s].AddEventVarFromTree("whoIsB2", o.test)
        SAMPLES[s].AddEventVarFromTree("whoIsJ1", o.test)
        SAMPLES[s].AddEventVarFromTree("whoIsJ2", o.test)

        SAMPLES[s].AddEventVarFromTree("mindRB1", o.test)
        SAMPLES[s].AddEventVarFromTree("mindRB2", o.test)
        SAMPLES[s].AddEventVarFromTree("mindRJ1", o.test)
        SAMPLES[s].AddEventVarFromTree("mindRJ2", o.test)

        SAMPLES[s].AddEventVarFromTree("MV2c10B1", o.test)
        SAMPLES[s].AddEventVarFromTree("MV2c10B2", o.test)
        SAMPLES[s].AddEventVarFromTree("MV2c10J1", o.test)
        SAMPLES[s].AddEventVarFromTree("MV2c10J2", o.test)

        SAMPLES[s].AddEventVarFromTree("nTightMv2c10", o.test)
        SAMPLES[s].AddEventVarFromTree("nMediumMv2c10", o.test)
        SAMPLES[s].AddEventVarFromTree("nLooseMv2c10", o.test)

        SAMPLES[s].AddEventVarFromTree("dRBB", o.test)
        SAMPLES[s].AddEventVarFromTree("max_J1J2", o.test)
        SAMPLES[s].AddEventVarFromTree("eta_J_star", o.test)

        SAMPLES[s].AddEventVarFromTree("NTrk1000PVJ1", o.test)
        SAMPLES[s].AddEventVarFromTree("NTrk1000PVJ2", o.test)
        SAMPLES[s].AddEventVarFromTree("NTrk500PVJ1", o.test)
        SAMPLES[s].AddEventVarFromTree("NTrk500PVJ2", o.test)

        SAMPLES[s].AddEventVarFromTree("pTB1", o.test)
        SAMPLES[s].AddEventVarFromTree("pTB2", o.test)

        if o.channel == "2cen":
            SAMPLES[s].AddEventVarFromTree("HadronConeExclTruthLabelB1", o.test)
            SAMPLES[s].AddEventVarFromTree("HadronConeExclTruthLabelB2", o.test)
            SAMPLES[s].AddCut("SecondBNotTrue",  SAMPLES[s].var["HadronConeExclTruthLabelB2"]<5)

        ### Load jet variables
#        for i in range(np.max(SAMPLES[s].var["nJets"])):
#            SAMPLES[s].AddJetVarsFromTree(["pT", "eta", "phi, "mv2c10"], i)
        ## make dR cuts
            #make_jet_dR_iso_Cuts(4, radius )
#        SAMPLES[s].var["mindRB1"]  = np.zeros(SAMPLES[s].var["nJets"].shape[0])
#        SAMPLES[s].var["mindRB2"]  = np.zeros(SAMPLES[s].var["nJets"].shape[0])
#        for i in range(np.max(SAMPLES[s].var["nJets"])):
#            ## cal mindR
#            thisJetmindR = cal_jet_min_dR(SAMPLES[s], i)
#

        SAMPLES[s].AddCut("2b_mv2c10_70",  np.logical_and((SAMPLES[s].var["MV2c10B1"]>0.8244273),(SAMPLES[s].var["MV2c10B2"]>0.8244273))==1)
        SAMPLES[s].AddCut("2b_mv2c10_77",  np.logical_and((SAMPLES[s].var["MV2c10B1"]>0.645925),(SAMPLES[s].var["MV2c10B2"]>0.645925))==1)
        SAMPLES[s].AddCut("2b_mv2c10_85",  np.logical_and((SAMPLES[s].var["MV2c10B1"]>0.1758475),(SAMPLES[s].var["MV2c10B2"]>0.1758475))==1)


        SAMPLES[s].AddCut("0b_mv2c10_70",  np.logical_or( np.logical_or((SAMPLES[s].var["MV2c10B1"]>0.8244273),
                                                                        (SAMPLES[s].var["MV2c10B2"]>0.8244273)), 
                                                          np.logical_or((SAMPLES[s].var["MV2c10J1"]>0.8244273),
                                                                        (SAMPLES[s].var["MV2c10J2"]>0.8244273)))==0)

        SAMPLES[s].AddCut("0b_mv2c10_77",  np.logical_or( np.logical_or((SAMPLES[s].var["MV2c10B1"]>0.645925),
                                                                        (SAMPLES[s].var["MV2c10B2"]>0.645925)), 
                                                          np.logical_or((SAMPLES[s].var["MV2c10J1"]>0.645925),
                                                                        (SAMPLES[s].var["MV2c10J2"]>0.645925)))==0)

        SAMPLES[s].AddCut("0b_mv2c10_85",  np.logical_or( np.logical_or((SAMPLES[s].var["MV2c10B1"]>0.1758475),
                                                                        (SAMPLES[s].var["MV2c10B2"]>0.1758475)), 
                                                          np.logical_or((SAMPLES[s].var["MV2c10J1"]>0.1758475),
                                                                        (SAMPLES[s].var["MV2c10J2"]>0.1758475)))==0)

        SAMPLES[s].AddCut("secondB_medium",  SAMPLES[s].var["MV2c10B2"]>0.8244273)

        SAMPLES[s].AddCut("pTBB>70", SAMPLES[s].var["pTBB"]>70)
        SAMPLES[s].AddCut("pTBB>90", SAMPLES[s].var["pTBB"]>90)
        SAMPLES[s].AddCut("pTBB>110", SAMPLES[s].var["pTBB"]>110)
        SAMPLES[s].AddCut("pTBB>130", SAMPLES[s].var["pTBB"]>130)
        SAMPLES[s].AddCut("pTBB>150", SAMPLES[s].var["pTBB"]>150)
        SAMPLES[s].AddCut("pTBB>170", SAMPLES[s].var["pTBB"]>170)
        SAMPLES[s].AddCut("pTBB>190", SAMPLES[s].var["pTBB"]>190)
        SAMPLES[s].AddCut("pTBB>210", SAMPLES[s].var["pTBB"]>210)
        SAMPLES[s].AddCut("pTBB>230", SAMPLES[s].var["pTBB"]>230)
        SAMPLES[s].AddCut("pTBB>250", SAMPLES[s].var["pTBB"]>250)
        SAMPLES[s].AddCut("pTBB>270", SAMPLES[s].var["pTBB"]>270)

        SAMPLES[s].AddCut("mJJ>50", SAMPLES[s].var["mJJ"]>50)
        SAMPLES[s].AddCut("mJJ>100", SAMPLES[s].var["mJJ"]>100)
        SAMPLES[s].AddCut("mJJ>150", SAMPLES[s].var["mJJ"]>150)
        SAMPLES[s].AddCut("mJJ>200", SAMPLES[s].var["mJJ"]>200)

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
        pTCut = 140
    for s in SAMPLES:
        SAMPLES[s].var["eventWeight"] = (SAMPLES[s].var["pTBB"]>pTCut) * SAMPLES[s].var["eventWeight"]
        SAMPLES[s].var["eventWeight_noblind"] = (SAMPLES[s].var["pTBB"]>pTCut) * SAMPLES[s].var["eventWeight_noblind"]


def ApplymJJCut():

    for s in SAMPLES:
        SAMPLES[s].var["eventWeight"] = (SAMPLES[s].var["mJJ"]>100) * SAMPLES[s].var["eventWeight"]


def Make1DPlots(samples, histname, var, weight, cuts, binning, precut=None):

    for s in samples:
        samples[s].Add1DHist(histname, binning, NoCut)
        if precut == None:
            samples[s].FillHist (histname, samples[s].var[var], samples[s].var[weight])
        else:
            samples[s].FillHist ( histname, samples[s].var[var], samples[s].var[weight], samples[s].cut[precut])
        
        for cut in cuts:
            samples[s].Add1DHist( histname + "_" + cut, binning, samples[s].cut[cut])
            if precut == None:
                samples[s].FillHist ( histname + "_" + cut, samples[s].var[var], samples[s].var[weight], samples[s].cut[cut])
            else:
                samples[s].FillHist ( histname + "_" + cut, samples[s].var[var], samples[s].var[weight], samples[s].cut[cut] + samples[s].cut[precut] )


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



def DrawMVAInputComparison(varlist):
    
    for var in varlist:
        for s in SAMPLES:
           SAMPLES[s].Norm1DHist(var)
           SAMPLES[s].Norm1DHist(var+"_sideband_low")
           SAMPLES[s].Norm1DHist(var+"_sideband_up")
           SAMPLES[s].Norm1DHist(var+"_PassBDTCut20")
           SAMPLES[s].Norm1DHist(var+"_PassNNCut20")
           SAMPLES[s].Norm1DHist(var+"_PassBDTCut40")
           SAMPLES[s].Norm1DHist(var+"_PassNNCut40")
           SAMPLES[s].Norm1DHist(var+"_PassBDTCut60")
           SAMPLES[s].Norm1DHist(var+"_PassNNCut60")
           SAMPLES[s].Norm1DHist(var+"_PassBDTCut80")
           SAMPLES[s].Norm1DHist(var+"_PassNNCut80")

           SAMPLES[s].var[var+"_corr"] = np.zeros(len(varlist))
           ivar2 = 0
           for var2 in varlist:
               SAMPLES[s].var[var+"_corr"][ivar2] = CalWeightedCorrCoeff( SAMPLES[s].var[var], SAMPLES[s].var[var2], SAMPLES[s].var["eventWeight"] )
               ivar2 += 1

        DrawTool.DrawHists(var+"_mvacomp",  ["", ""], [SAMPLES["VBF"].histograms[var+"_Normed"], SAMPLES["ggF"].histograms[var+"_Normed"], 
                                                       SAMPLES["data_0tag"].histograms[var+"_Normed"], SAMPLES["data_2tag"].histograms[var+"_Normed"]],
                                                       ["VBF", "ggF", "data 0tag", "data 2tag"])

        DrawTool.DrawHists(var+"_mvacomp_sidband",  ["", ""], [ SAMPLES["data_2tag"].histograms[var+"_sideband_low_Normed"], 
                                                                SAMPLES["data_2tag"].histograms[var+"_sideband_up_Normed"], 
                                                                SAMPLES["VBF"].histograms[var+"_Normed"]],
                           ["data 2tag sideband low", "data 2tag sideband up", "VBF"])

    for s in SAMPLES:
        SAMPLES[s].Add2DHist("MVAInput_Corr", np.array(range(len(varlist))), np.array(range(len(varlist))) )
        for ivar in range(len(varlist)):
            SAMPLES[s].histograms["MVAInput_Corr"].GetXaxis().SetBinLabel(ivar+1, varlist[ivar])

            for jvar in range(len(varlist)):
                SAMPLES[s].histograms["MVAInput_Corr"].GetYaxis().SetBinLabel(jvar+1, varlist[jvar])
                SAMPLES[s].histograms["MVAInput_Corr"].SetBinContent(ivar+1, jvar+1, round(SAMPLES[s].var[varlist[ivar]+"_corr"][jvar], 4))

        SAMPLES[s].histograms["MVAInput_Corr"].Write()
        SAMPLES[s].Add2DHist("NNScore_mBB",  NN_Binning,  JetpT_Binning)
        SAMPLES[s].Add2DHist("BDTScore_mBB", BDT_Binning, JetpT_Binning)

        NN_mBB  = np.dstack( (MVAPRED["NN_"+s+"_full"],  SAMPLES[s].var["mBB"]))[0]
        BDT_mBB = np.dstack( (MVAPRED["BDT_"+s+"_full"], SAMPLES[s].var["mBB"]))[0]

        SAMPLES[s].FillHist( "NNScore_mBB",  NN_mBB, SAMPLES[s].var["eventWeight"])
        SAMPLES[s].FillHist( "BDTScore_mBB", BDT_mBB, SAMPLES[s].var["eventWeight"])
        SAMPLES[s].histograms["NNScore_mBB"].Write()
        SAMPLES[s].histograms["BDTScore_mBB"].Write()

        print (s, "NN mBB corr", SAMPLES[s].histograms["NNScore_mBB"].GetCorrelationFactor(), "BDT mBB corr", SAMPLES[s].histograms["BDTScore_mBB"].GetCorrelationFactor())
        

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

    #dataset_NN, dataset_NN_sys  = makeData(True,  False, MVAinputs["data_2tag"], MVAinputs["data_0tag"], MVAinputs["VBF"], MVAinputs["ggF"], 
    #                                       MVAinputs["Zbb_QCD"], MVAinputs["Zbb_EWK"],MVAinputs["Zbb_QCD_0tag"], MVAinputs["Zbb_EWK_0tag"], MVAinputs_sys )

    dataset_BDT = None
    dataset_BDT_sys= None

    if o.channel == "4cen":
        dataset_BDT, dataset_BDT_sys = makeData(False, False, MVAinputs["data_2tag"], MVAinputs["data_0tag"], MVAinputs["VBF"], MVAinputs["ggF"], 
                                                MVAinputs["Zbb_QCD"], MVAinputs["Zbb_EWK"],MVAinputs["Zbb_QCD_0tag"], MVAinputs["Zbb_EWK_0tag"], MVAinputs_sys )
    if o.channel == "2cen":
        dataset_BDT, dataset_BDT_sys = makeData(False, False, MVAinputs["data_2tag"], MVAinputs["data_0tag"], MVAinputs["VBF"], MVAinputs["ggF"], 
                                                MVAinputs["Zbb_QCD"], MVAinputs["Zbb_EWK"],MVAinputs["Zbb_QCD_0tag"], MVAinputs["Zbb_EWK_0tag"], MVAinputs_sys )

    model_BDT  = buildBDT(dataset_BDT, doTest = o.test)
    PredictModel(model_BDT, "BDT", dataset_BDT, dataset_BDT_sys)

    if o.channel == "4cen":
        DrawTool.DrawHists("BDT score", ["BDT score", "A.U."], 
                           NormHists([ MVAScoreHists["BDT_sig_train"], MVAScoreHists["BDT_sig_test"], MVAScoreHists["BDT_bkg_train"], 
                                       MVAScoreHists["BDT_bkg_test"], MVAScoreHists["BDT_ggF"], MVAScoreHists["BDT_data_0tag"], 
                                       MVAScoreHists["BDT_Zbb_QCD"], MVAScoreHists["BDT_Zbb_EWK"]]),
                           ["Signal Training Sample", "Signal Test Sample","Background Training Sample", "Background Test Sample", "ggF", "BDT_data_0tag",
                            "Zbb QCD", "Zbb EWK"])
    else:
        DrawTool.DrawHists("BDT score", ["BDT score", "A.U."], 
                           NormHists([ MVAScoreHists["BDT_sig_train"], MVAScoreHists["BDT_sig_test"], MVAScoreHists["BDT_bkg_train"], 
                                       MVAScoreHists["BDT_bkg_test"], MVAScoreHists["BDT_ggF"], MVAScoreHists["BDT_data_0tag"]]),
                           ["Signal Training Sample", "Signal Test Sample","Background Training Sample", "Background Test Sample", "ggF", "BDT_data_0tag"])


    DrawTool.getROC( [ MVAPRED["BDT_sig_train"], MVAPRED["BDT_sig_test"]],
                     [ MVAPRED["BDT_bkg_train"], MVAPRED["BDT_bkg_test"]],
                     ["BDT train", "BDT test"] )

    BDT_Cut_00  = np.percentile( MVAPRED["BDT_VBF_full"], 0 )
    BDT_Cut_03  = np.percentile( MVAPRED["BDT_VBF_full"], 3 )
    BDT_Cut_05  = np.percentile( MVAPRED["BDT_VBF_full"], 5 )
    BDT_Cut_10  = np.percentile( MVAPRED["BDT_VBF_full"], 10 )
    BDT_Cut_20  = np.percentile( MVAPRED["BDT_VBF_full"], 20 )
    BDT_Cut_30  = np.percentile( MVAPRED["BDT_VBF_full"], 30 )
    BDT_Cut_40  = np.percentile( MVAPRED["BDT_VBF_full"], 40 )
    BDT_Cut_50  = np.median( MVAPRED["BDT_VBF_full"] )
    BDT_Cut_60  = np.percentile( MVAPRED["BDT_VBF_full"], 60 )
    BDT_Cut_70  = np.percentile( MVAPRED["BDT_VBF_full"], 70 )
    BDT_Cut_80  = np.percentile( MVAPRED["BDT_VBF_full"], 80 )
    BDT_Cut_90  = np.percentile( MVAPRED["BDT_VBF_full"], 90 )
    BDT_Cut_95  = np.percentile( MVAPRED["BDT_VBF_full"], 95 )
    
    for s in SAMPLES:

        SAMPLES[s].AddCut("PassBDTCut00",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00 )
        SAMPLES[s].AddCut("PassBDTCut10",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_10 )
        SAMPLES[s].AddCut("PassBDTCut20",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_20 )
        SAMPLES[s].AddCut("PassBDTCut30",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_30 )
        SAMPLES[s].AddCut("PassBDTCut40",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_40 )
        SAMPLES[s].AddCut("PassBDTCut50",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_50 )
        SAMPLES[s].AddCut("PassBDTCut60",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_60 )
        SAMPLES[s].AddCut("PassBDTCut70",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_70 )
        SAMPLES[s].AddCut("PassBDTCut80",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_80 )
        SAMPLES[s].AddCut("PassBDTCut90",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_90 )
        SAMPLES[s].AddCut("PassBDTCut95",   MVAPRED["BDT_"+s+"_full"]>BDT_Cut_95 )
        
        SAMPLES[s].AddCut("PassBDTCut00-03",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_03,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut00-05",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_05,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut00-10",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_10,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut00-20",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_20,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_00)  )
        SAMPLES[s].AddCut("PassBDTCut20-40",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_40,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_20)  )
        SAMPLES[s].AddCut("PassBDTCut40-60",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_60,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_40)  )
        SAMPLES[s].AddCut("PassBDTCut60-80",    np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_80,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_60)  )
        SAMPLES[s].AddCut("PassBDTCut80-100",    MVAPRED["BDT_"+s+"_full"]>BDT_Cut_80 )


        SAMPLES[s].AddCut("BDTI",    MVAPRED["BDT_"+s+"_full"]>BDT_Cut_80 )
        SAMPLES[s].AddCut("BDTII",   np.logical_and(MVAPRED["BDT_"+s+"_full"]<BDT_Cut_80,MVAPRED["BDT_"+s+"_full"]>BDT_Cut_60)  )


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
        
        nvbf = sum(SAMPLES["VBF"].var["eventWeight"] * SAMPLES["VBF"].cut["signal"].array * SAMPLES["VBF"].cut[cut].array / 10*9.9)
        nggf = sum(SAMPLES["ggF"].var["eventWeight"] * SAMPLES["ggF"].cut["signal"].array * SAMPLES["ggF"].cut[cut].array / 10*9.9)
        n0tag = sum(SAMPLES["data_0tag"].var["eventWeight"] * SAMPLES["data_0tag"].cut["signal"].array * SAMPLES["data_0tag"].cut[cut].array * scale_factor )

        print ("N 2tag sideband", sideband_2tag_evts)
        print ("N VBF",  nvbf)
        print ("N ggF",  nggf)
        print ("N 0tag", n0tag)
        print ("S/SQRT(B)", (nvbf+nggf)/sqrt(n0tag))
        print ("S/SQRT(B*(1+0.0001B))", (nvbf+nggf)/sqrt(n0tag*(1+0.0001*n0tag)))
        print ("------------------")


def DrawPlotsWithCutsForEachSample(var, canvname, xLabel, yLabel, cuts, norm=False):

    for s in SAMPLES:
        plots = []
        labels = []

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


def GenerateSignalWS(process = "signal", cut = None):

    xmlfile = open("FitFiles_"+o.channel+"/"+process+"_"+cut+"_"+o.channel+".xml", "a")

    xmlfile.write('''<!DOCTYPE Model SYSTEM 'AnaWSBuilder.dtd'>'''+"\n")
    xmlfile.write('''    <Model Type="External" Input="config/vbf/FitFiles_{!s}/WS_{!s}_{!s}_{!s}.root" WSName="combined" ModelName="C_{!s}_model" ObservableName="obs_x_C_{!s}">'''.format(o.channel, process, o.channel, cut, cut, cut)+"\n")
    xmlfile.write('''    <Item Name="unit[1]"/>\n''')
    xmlfile.write('''    <Rename OldName="Lumi" NewName="unit"/>\n''')


    signal_hist = None;
    if process == "signal":
        signal_hist = deepcopy(  SAMPLES["VBF"].histograms["mBB_short_"+cut] )
        signal_hist.Add       (  SAMPLES["ggF"].histograms["mBB_short_"+cut] )

    if process == "z":
        signal_hist = deepcopy(  SAMPLES["Zbb_QCD"].histograms["mBB_short_"+cut] )
        signal_hist.Add       (  SAMPLES["Zbb_EWK"].histograms["mBB_short_"+cut] )

    ### generate workspace for signal
    meas_sig = RooStats.HistFactory.Measurement(process+"_"+o.channel+"_"+cut, process+"_"+o.channel+"_"+cut);
    meas_sig.SetOutputFilePrefix("WS")
    meas_sig.SetExportOnly(1);
    meas_sig.SetPOI("xs");
    meas_sig.SetLumi(1);
    meas_sig.SetLumiRelErr(0.01);

    for ibin in range(signal_hist.GetNbinsX()):
        if  signal_hist.GetBinContent(ibin+1) ==0:
            signal_hist.SetBinContent(ibin+1, 0.000001)
            signal_hist.SetBinError(ibin+1,  0.000001)

    R_Channel_cut = RooStats.HistFactory.Channel("C_"+cut);
    R_Sample_cut = RooStats.HistFactory.Sample("S_"+cut);
    R_Sample_cut.SetHisto(signal_hist);

    doStat = True
    doSyst = False

    if doStat:
        R_Channel_cut.SetStatErrorConfig(0.0, "Poisson");
        R_Sample_cut.ActivateStatError();
        
        for ibin in range(signal_hist.GetNbinsX()):
            xmlfile.write('''    <ExtSyst ConstrName="gamma_stat_C_{!s}_bin_{!s}_constraint" NPName="gamma_stat_C_{!s}_bin_{!s}" GOName="nom_gamma_stat_C_{!s}_bin_{!s}" />'''.format(cut, ibin, cut, ibin, cut, ibin)+"\n")

    if doSyst:
        for key in SYSNAMES:
            
            xmlfile.write('''    <ExtSyst ConstrName="alpha_ATLAS_{!s}_Constraint" NPName="alpha_{!s}" GOName="nom_alpha_ATLAS_{!s}" />'''.format(key, key, key)+"\n")

            R_Sys = RooStats.HistFactory.HistoSys("ATLAS_"+key)

            signal_hist_sys_up = None
            if process == "signal":
                signal_hist_sys_up.deepcopy(SAMPLES["VBF_"+key+"__1up"].histograms["mBB_short_"+cut])
                signal_hist_sys_up.Add     (SAMPLES["ggF_"+key+"__1up"].histograms["mBB_short_"+cut])
            if process == "z":
                signal_hist_sys_up.deepcopy(SAMPLES["Zbb_QCD_"+key+"__1up"].histograms["mBB_short_"+cut])
                signal_hist_sys_up.Add     (SAMPLES["Zbb_EWK_"+key+"__1up"].histograms["mBB_short_"+cut])

            R_Sys.SetHistoHigh(signal_hist_sys_up)

            if "JER" not in key:
                signal_hist_sys_dn = None
                if process == "signal":
                    signal_hist_sys_dn.deepcopy(SAMPLES["VBF_"+key+"__1dn"].histograms["mBB_short_"+cut])
                    signal_hist_sys_dn.Add     (SAMPLES["ggF_"+key+"__1dn"].histograms["mBB_short_"+cut])
                if process == "z":
                    signal_hist_sys_dn.deepcopy(SAMPLES["Zbb_QCD_"+key+"__1dn"].histograms["mBB_short_"+cut])
                    signal_hist_sys_dn.Add     (SAMPLES["Zbb_EWK_"+key+"__1dn"].histograms["mBB_short_"+cut])

                R_Sys.SetHistoLow(signal_hist_sys_dn)
            else:
                R_Sys.SetHistoLow(CopyHist(signal_hist))
            R_Sample_cut.AddHistoSys(R_Sys)

    ### Initialize and save work space for signal
    R_Channel_cut.AddSample(R_Sample_cut);
    meas_sig.AddChannel(R_Channel_cut)
    RooStats.HistFactory.MakeModelAndMeasurementFast(meas_sig)
    os.system("rm WS_c*"+o.channel+"*")
    os.system("mv WS_*"+o.channel+"* FitFiles_"+o.channel)

    xmlfile.write('''</Model>''')    



def StandaloneFit(data_hist, z_hist, norm, fitname, fit_function, cut):
    
    signal_hist = TH1D("signal_"+fitname+"_"+cut, "signal_"+fitname+"_"+cut, len(mBB_Binning)-1, mBB_Binning)

    if "spurious" in cut:
        signal_hist.Add(SAMPLES["VBF"].histograms["mBB_short_"+cut.lstrip("spurious_")])
        signal_hist.Add(SAMPLES["ggF"].histograms["mBB_short_"+cut.lstrip("spurious_")])
        signal_hist.Scale(norm)
        
    else:
        signal_hist.Add(SAMPLES["VBF"].histograms["mBB_short_"+cut])
        signal_hist.Add(SAMPLES["ggF"].histograms["mBB_short_"+cut])

    output.cd()
    fit_results_txt = open("FitFiles_"+o.channel+"/f_"+fitname+"_"+cut+"_"+o.channel, "a")

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

    if "spurious" in cut:
        signal_only = CopyHist(signal_hist)
        signal_plus_bkg = CopyHist(signal_hist)
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
        fit_results_txt.write( "B "+str( data_hist.Integral()  )+"\n")
        out_data_array = np.array([])
        out_weight_array = np.array([])

        for ibin in range(data_hist.GetNbinsX()):
            bincenter = data_hist.GetBinCenter(ibin+1)
            bincontent = data_hist.GetBinContent(ibin+1)
            out_data_array = np.concatenate( (bincenter*np.ones(int(bincontent)), out_data_array) )

        out_weight_array = np.ones( len(out_data_array))

        np.savetxt("FitFiles_"+o.channel+"/f_pseudo_data_"+o.channel+"_fit_"+fitname+"_"+cut, out_data_array, fmt='%1.6f')
        np.savetxt("FitFiles_"+o.channel+"/f_pseudo_data_"+o.channel+"_fit_"+fitname+"_weight_"+cut, out_weight_array, fmt='%1.6f')

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



def FitmBB(cuts, WritePseudoData = False):
    #FITRESULTS["0tag"]= {}

#    data_hist_1 = CopyHist(SAMPLES["data_2tag"].histograms["mBB_short_PassBDTCut80-100"])
#    data_hist_2 = CopyHist(SAMPLES["data_2tag"].histograms["mBB_short_PassBDTCut60-80"])
#    z_hist_1 = CopyHist(SAMPLES["Zbb_QCD"].histograms["mBB_short_fitband_PassBDTCut80-100"])
#    z_hist_1.Add(SAMPLES["Zbb_EWK"].histograms["mBB_short_fitband_PassBDTCut80-100"])
#    z_hist_2 = CopyHist(SAMPLES["Zbb_QCD"].histograms["mBB_short_fitband_PassBDTCut60-80"])
#    z_hist_2.Add(SAMPLES["Zbb_EWK"].histograms["mBB_short_fitband_PassBDTCut60-80"])
#
#    GenerateFitWorkSpace2Region(data_hist_1, data_hist_2, z_hist_1, z_hist_2, "BernsteinO2Lin", BernsteinO2Lin)
#    GenerateFitWorkSpace2Region(data_hist_1, data_hist_2, z_hist_1, z_hist_2, "BernsteinO3Lin", BernsteinO3Lin)
#    GenerateFitWorkSpace2Region(data_hist_1, data_hist_2, z_hist_1, z_hist_2, "BernsteinO4Lin", BernsteinO4Lin)

    for cut in cuts:
        print "cut ", cut

        ## fit spurious data
        FITRESULTS["spurious_"+cut]= {}

        norm_0tag = None
        if o.spurious == "0tag":
            norm_0tag = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short_fitband"])
        if o.spurious == "2tag":
            norm_0tag = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_fitband_noblind_PassBDTCut00-03"])
        if o.spurious == "mix":
            norm_0tag = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short_fitband"])
            norm_0tag.Add(SAMPLES["data_2tag"].histograms["mBB_short_fitband_noblind_PassBDTCut00-03"])
            
        n_0tag_fitband= norm_0tag.Integral()
        norm_0tag.Scale(1.0/n_0tag_fitband)

        norm_2tag = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_fitband_"+cut])
        n_2tag_fitband = norm_2tag.Integral()
        norm_2tag.Scale(1.0/n_2tag_fitband)

        ratio = deepcopy(norm_2tag)
        ratio.Divide(norm_0tag)
        ratio.Fit("pol1")
        linfit = ratio.GetFunction("pol1")
        linfit.Write()

        data_0tag_fitband = None
        if o.spurious == "0tag":
            data_0tag_fitband = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short_fitband"])
        if o.spurious == "2tag":
            data_0tag_fitband = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_fitband_noblind_PassBDTCut00-03"])
        if o.spurious == "mix":
            data_0tag_fitband = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short_fitband"])
            data_0tag_fitband.Add(SAMPLES["data_2tag"].histograms["mBB_short_fitband_noblind_PassBDTCut00-03"])

        data_0tag_fitband.Multiply(linfit)
        n_0tag_fitband_postcor= data_0tag_fitband.Integral()
        
        data_hist = None
        if o.spurious == "0tag":
            data_hist = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short"])
        if o.spurious == "2tag":
            data_hist = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_noblind_PassBDTCut00-03"])
        if o.spurious == "mix":
            data_hist = deepcopy(SAMPLES["data_0tag"].histograms["mBB_short"])
            data_hist.Add(SAMPLES["data_2tag"].histograms["mBB_short_noblind_PassBDTCut00-03"])

        data_hist.Multiply(linfit)

        data_norm = n_0tag_fitband/n_0tag_fitband_postcor ## normalize data to pre-correction
        mc_norm = n_0tag_fitband/n_2tag_fitband ## normalize mc to acceptance in CR
        data_hist.Scale(data_norm)

        z_hist = CopyHist(SAMPLES["Zbb_QCD_0tag"].histograms["mBB_short"])
        z_hist.Add(SAMPLES["Zbb_EWK_0tag"].histograms["mBB_short"])

        StandaloneFit( data_hist, z_hist, mc_norm, "BernsteinO2", BernsteinO2, "spurious_"+cut)
        StandaloneFit( data_hist, z_hist, mc_norm, "BernsteinO3", BernsteinO3, "spurious_"+cut)
        StandaloneFit( data_hist, z_hist, mc_norm, "BernsteinO4", BernsteinO4, "spurious_"+cut)
        StandaloneFit( data_hist, z_hist, mc_norm, "BernsteinO5", BernsteinO5, "spurious_"+cut)

        StandaloneFit( data_hist, z_hist, mc_norm, "ExpoBernsteinO2", ExpoBernsteinO2, "spurious_"+cut)
        StandaloneFit( data_hist, z_hist, mc_norm, "ExpoBernsteinO3", ExpoBernsteinO3, "spurious_"+cut)
        StandaloneFit( data_hist, z_hist, mc_norm, "ExpoBernsteinO4", ExpoBernsteinO4, "spurious_"+cut)
        StandaloneFit( data_hist, z_hist, mc_norm, "ExpoBernsteinO5", ExpoBernsteinO5, "spurious_"+cut)

        ### fit 2tag data
        #FITRESULTS[cut]= {}
        #z_hist = CopyHist(SAMPLES["Zbb_QCD"].histograms["mBB_short_fitband_"+cut])
        #z_hist.Add(SAMPLES["Zbb_EWK"].histograms["mBB_short_fitband_"+cut])
        #
        #data_hist = CopyHist(SAMPLES["data_2tag"].histograms["mBB_short_"+cut])
        #StandaloneFit( data_hist, z_hist, "BernsteinO2", BernsteinO2, cut)
        #StandaloneFit( data_hist, z_hist, "BernsteinO3", BernsteinO3, cut)
        #StandaloneFit( data_hist, z_hist, "BernsteinO4", BernsteinO4, cut)
        #StandaloneFit( data_hist, z_hist, "BernsteinO5", BernsteinO5, cut)

        GenerateSignalWS("signal", cut)
        GenerateSignalWS("z", cut)

        output.cd()


###### main code########            
wp = o.btag

btagcuts = ["2b_mv2c10_70", "2b_mv2c10_77", "2b_mv2c10_85", "0b_mv2c10_70", "0b_mv2c10_77", "0b_mv2c10_85"]
channelcuts_2cen = ["secondB_medium"]
pTBBcuts = ["pTBB>70", "pTBB>90", "pTBB>110", "pTBB>130", "pTBB>150", "pTBB>170", "pTBB>190", "pTBB>210", "pTBB>230", "pTBB>250", "pTBB>270"]
mJJcuts  = ["mJJ>50", "mJJ>100", "mJJ>150", "mJJ>200"]
sidebandcuts = ["sideband_up", "sideband_low"]
fitbandcuts = ["fitband_up", "fitband_low", "fitband"]
BDTcuts  = ["PassBDTCut00", "PassBDTCut20", "PassBDTCut40", "PassBDTCut60", "PassBDTCut80", "PassBDTCut90", "PassBDTCut95"]
BDTcuts  = ["PassBDTCut00", "PassBDTCut20", "PassBDTCut40", "PassBDTCut60", "PassBDTCut80", "PassBDTCut90", "PassBDTCut20-40", "PassBDTCut40-60", "PassBDTCut60-80", "PassBDTCut80-100"]
BDTcuts_More = ["BDTI", "BDTII", "PassBDTCut00-03","PassBDTCut00-05","PassBDTCut00-10","PassBDTCut00-20", "PassBDTCut20-40","PassBDTCut40-60","PassBDTCut60-80","PassBDTCut80-100"]
BDTcuts = ["PassBDTCut60-80","PassBDTCut80-100"]
BDTcuts = ["BDTI", "BDTII"]

if o.test :
    BDTcuts = ["BDTI", "BDTII"]


MVAVarList_Check= ["mBB", "pTBB", "mJJ", "pTJJ", "dRBB", "dEtaJJ", "WidthJ1", "WidthJ2", "cosTheta_MVA", "max_J1J2", "eta_J_star", "HT_soft", "HT_MVA", "mindRB1", "mindRB2", "mindRJ1","mindRJ2", "MV2c10B1", "MV2c10B2"]
MVAVarList_Use= ["mJJ", "pTJJ", "dEtaJJ", "WidthJ1", "WidthJ2", "cosTheta_MVA", "HT_soft", "HT_MVA", "mindRJ1", "mindRJ2", "MV2c10B1", "MV2c10B2"]

if o.MVAList == "long":
    MVAVarList_Use= ["mJJ", "pTJJ", "dEtaJJ", "WidthJ1", "WidthJ2", "cosTheta_MVA", "HT_soft", "HT_MVA", "max_J1J2", "eta_J_star", "mindRJ1", "mindRJ2"]
if o.MVAList == "short":
    MVAVarList_Use= ["mJJ", "pTJJ", "dEtaJJ", "cosTheta_MVA", "HT_MVA", "mindRJ1", "mindRJ2", "max_J1J2", "eta_J_star"]
if o.MVAList == "shortandwidth":
    MVAVarList_Use= ["mJJ", "pTJJ", "dEtaJJ", "cosTheta_MVA", "HT_MVA", "mindRJ1", "mindRJ2", "max_J1J2", "eta_J_star", "WidthJ1", "WidthJ2"]
if o.MVAList == "shortandntrk500":
    MVAVarList_Use= ["mJJ", "pTJJ", "dEtaJJ", "cosTheta_MVA", "HT_MVA", "mindRJ1", "mindRJ2", "max_J1J2", "eta_J_star", "NTrk500PVJ1", "NTrk500PVJ2"]
if o.MVAList == "shortandntrk1000":
    MVAVarList_Use= ["mJJ", "pTJJ", "dEtaJJ", "cosTheta_MVA", "HT_MVA", "mindRJ1", "mindRJ2", "max_J1J2", "eta_J_star", "NTrk1000PVJ1", "NTrk1000PVJ2"]


#############################
####    Main Code       #####
#############################


################################
#### Initialize and Apply Cut###
################################

InitSamples()

if o.ApplypTBBCut:
    ApplypTBBCut()

output = TFile("out_"+ o.channel+"_" + o.ApplypTBBCut*"_pTBBCut" + "_" + o.MVAList + ".root", "recreate")
output.cd()

BlindData()
ApplyBtag()

ComputeMVAWeights(MVAVarList_Use)
CalculateSensitivity(BDTcuts_More)
#DrawMVAInputComparison(MVAVarList_Use)

Make1DPlots(SAMPLES, "NJets",      "nJets",           "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), NJet_Binning)
Make1DPlots(SAMPLES, "mBB",        "mBB",             "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
Make1DPlots(SAMPLES, "mBB_no_corr","mBB_no_corr",     "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
Make1DPlots(SAMPLES, "mBB_short",  "mBB",             "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), mBB_Binning)
Make1DPlots(SAMPLES, "mBB_short_fitband",  "mBB",     "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), mBB_Binning, precut = "fitband")
Make1DPlots(SAMPLES, "mBBsig",     "mBB",             "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), mBB_Binning_sig)
Make1DPlots(SAMPLES, "mBBsig_no_corr",  "mBB_no_corr","eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), mBB_Binning_sig)
Make1DPlots(SAMPLES, "pTBB",       "pTBB",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
Make1DPlots(SAMPLES, "mJJ",        "mJJ",             "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), MJJ_Binning)
Make1DPlots(SAMPLES, "pTJJ",       "pTJJ",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)

Make1DPlots(SAMPLES, "mBB_short_noblind",  "mBB",             "eventWeight_noblind", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), mBB_Binning)
Make1DPlots(SAMPLES, "mBB_short_fitband_noblind",  "mBB",     "eventWeight_noblind", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), mBB_Binning, precut = "fitband")

Make1DPlots(SAMPLES, "pTB1",       "pTB1",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
Make1DPlots(SAMPLES, "pTB2",       "pTB2",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
#Make1DPlots(SAMPLES, "dRBB",       "dRBB",            "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), dR_Binning)
#Make1DPlots(SAMPLES, "dEtaJJ",     "dEtaJJ",          "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), dR_Binning)
#Make1DPlots(SAMPLES, "dEtaBB",     "dEtaBB",          "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), dR_Binning)
#Make1DPlots(SAMPLES, "WidthJ1",    "WidthJ1",         "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetWidth_Binning)
#Make1DPlots(SAMPLES, "WidthJ2",    "WidthJ2",         "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetWidth_Binning)
#Make1DPlots(SAMPLES, "cosTheta_MVA", "cosTheta_MVA",  "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), CosTheta_Binning)
#Make1DPlots(SAMPLES, "HT_soft",    "HT_soft",         "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
#Make1DPlots(SAMPLES, "HT_MVA",     "HT_MVA",          "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), JetpT_Binning)
#Make1DPlots(SAMPLES, "mindRB1",    "mindRB1",         "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), dR_Binning)
#Make1DPlots(SAMPLES, "mindRB2",    "mindRB2",         "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), dR_Binning)
#Make1DPlots(SAMPLES, "mindRJ1",    "mindRJ1",         "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), dR_Binning)
#Make1DPlots(SAMPLES, "mindRJ2",    "mindRJ2",         "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), dR_Binning)
#Make1DPlots(SAMPLES, "max_J1J2",   "max_J1J2",        "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), Jeteta_Binning)
#Make1DPlots(SAMPLES, "eta_J_star", "eta_J_star",      "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), Jeteta_Binning)
#Make1DPlots(SAMPLES, "MV2c10B1",   "MV2c10B1",        "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), MV2_Binning)
#Make1DPlots(SAMPLES, "MV2c10B2",   "MV2c10B2",        "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts  +  BDTcuts_More ), MV2_Binning)
#Make1DPlots(SAMPLES, "HadronConeExclTruthLabelB1",    "HadronConeExclTruthLabelB1",         "eventWeight", UniqueList(channelcuts_2cen +  pTBBcuts ), PDGID_Binning)
#Make1DPlots(SAMPLES, "HadronConeExclTruthLabelB2",    "HadronConeExclTruthLabelB2",         "eventWeight", UniqueList(channelcuts_2cen +  pTBBcuts ), PDGID_Binning)

DrawPlotsWithCutsForEachSample("mBB",  "pTBBCuts", "M(bb)",  "Events", pTBBcuts)
DrawPlotsWithCutsForEachSample("pTBB", "pTBBCuts", "M(bb)",  "Events", pTBBcuts)


#############################
####  Customized Plots  #####
#############################

if o.channel== "2cen":
    Make1DPlots(SAMPLES, "mBB_short",  "mBB",  "eventWeight", ["SecondBNotTrue"], mBB_Binning)
    
    DrawTool.DrawHists("mBB_Full_Vs_NotTrueB",  ["M(bb)", "Events"], 
                       [SAMPLES["VBF"].histograms["mBB_short"], SAMPLES["VBF"].histograms["mBB_short_SecondBNotTrue"]],
                       ["Full", "Second B Not True"])

DrawTool_Ratio.DrawHists("mBB_0tag_vs_2tag",  ["M(bb)", "Events"], 
                         [SAMPLES["data_0tag"].histograms["mBB_short"], SAMPLES["data_2tag"].histograms["mBB_short" ]],
                         ["data_0tag", "data_2tag"])

#### Draw M(bb) w/ and w/o correction
#DrawTool_Ratio.DrawHists("mBB_corr_vs_nocorr",  ["M(bb)", "Events"], [SAMPLES["VBF"].histograms["mBBsig"],
#                                                                      SAMPLES["VBF"].histograms["mBBsig_no_corr"]], ["W/ correction", "W/O correction"])

### Trig and Jet assignment cross check
#Make1DPlots(SAMPLES, "mBB_short",        "mBB",             "eventWeight", UniqueList(btagcuts +  pTBBcuts + sidebandcuts), JetpT_Binning)
#DrawTool.DrawHists("mBB_Trig",  ["M(bb)", "Events"], [SAMPLES["VBF_HYBRID"].histograms["mBB_short"],
#                                                      SAMPLES["VBF_BMED"].histograms["mBB_short"],
#                                                      SAMPLES["VBF_BMED_TrigMatch_NoWeight"].histograms["mBB_short"],
#                                                      SAMPLES["VBF_BMED_TrigMatch"].histograms["mBB_short"]],
#                   ["TDT PC BTag ", "TDT Medium BTag", "Offline Emulation Medium BTag w/o weights ", "Offline Emulation Medium BTag w/ weights "])


#### Draw M(bb) for 0-tag vs 2-tag with BDT cut
#### Draw M(bb) for 2-tag BDT cuts comparison
for icut in range(len(BDTcuts_More)):
    cut1 = BDTcuts_More[icut]
    DrawTool_Ratio.DrawHists("mBB_0tag_vs_2tag_"+cut1,  ["M(bb)", "Events"], 
                             [SAMPLES["data_0tag"].histograms["mBB_short"], SAMPLES["data_2tag"].histograms["mBB_short_"+cut1 ]], 
                             ["data_0tag", "data_2tag_"+cut1])

    DrawTool_Ratio.DrawHists("mBB_2tag_3%BDT_vs_2tag_"+cut1,  ["M(bb)", "Events"], 
                             [SAMPLES["data_2tag"].histograms["mBB_short_PassBDTCut00-03"], SAMPLES["data_2tag"].histograms["mBB_short_"+cut1 ]], 
                             ["data_2tag_3%BDT", "data_2tag_"+cut1])

    for jcut in range(icut+1, len(BDTcuts_More)):
        cut2 = BDTcuts_More[jcut]

        tmp_cut1 = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_"+cut1])
        tmp_cut2 = deepcopy(SAMPLES["data_2tag"].histograms["mBB_short_"+cut2])
        tmp_cut1.Scale( 1./tmp_cut1.Integral())
        tmp_cut2.Scale( 1./tmp_cut2.Integral())

        DrawTool_Ratio.DrawHists("mBB_2tag_"+cut1+"_vs_"+cut2,  ["M(bb)", "Events"], 
                                 [tmp_cut1, tmp_cut2],
                                 ["data_2tag_"+cut1, "data_2tag_"+cut2])

#########################
####  Cut Study  ########
#########################

#DrawPlotsWithCutsForEachSample("HadronConeExclTruthLabelB1", "pTBBCuts", "pdgID",  "Events", pTBBcuts)
#DrawPlotsWithCutsForEachSample("HadronConeExclTruthLabelB1", "btag", "pdgID",  "Events", channelcuts_2cen)
#DrawPlotsWithCutsForEachSample("HadronConeExclTruthLabelB2", "pTBBcuts", "pdgID",  "Events", pTBBcuts)
#DrawPlotsWithCutsForEachSample("HadronConeExclTruthLabelB2", "btag", "pdgID",  "Events", channelcuts_2cen)

#DrawPlotsWithCutsForEachSample("mBB_no_corr",  "pTBBCuts", "M(bb)",  "Events", pTBBcuts)
#DrawPlotsWithCutsForEachSample("pTBB", "pTBBCuts", "pT(bb)", "Events", pTBBcuts)
#DrawPlotsWithCutsForEachSample("mJJ",  "pTBBCuts", "M(JJ)",  "Events", pTBBcuts)
#DrawPlotsWithCutsForEachSample("dRBB", "pTBBCuts", "dR(bb)",  "Events", pTBBcuts)
#DrawPlotsWithCutsForEachSample("dEtaBB", "pTBBCuts", "dEta(bb)",  "Events", pTBBcuts)

#DrawPlotsWithCutsForEachSample("mBB",  "BDTCuts", "M(bb)",     "Events", BDTcuts)
#DrawPlotsWithCutsForEachSample("mJJ",  "BDTCuts", "M(JJ)",     "Events", BDTcuts)
#DrawPlotsWithCutsForEachSample("pTBB", "BDTCuts", "pT(bb)",    "Events", BDTcuts)
#DrawPlotsWithCutsForEachSample("dRBB", "BDTCuts", "dR(bb)",  "Events", BDTcuts)
#DrawPlotsWithCutsForEachSample("dEtaBB", "BDTCuts", "dEta(bb)","Events", BDTcuts)

#DrawPlotsWithCutsForEachSample("mBB",  "BDTCuts_Normed", "M(bb)",  "Events", BDTcuts, norm=True)
#DrawPlotsWithCutsForEachSample("mJJ",  "BDTCuts_Normed", "M(JJ)",  "Events", BDTcuts, norm=True)
#DrawPlotsWithCutsForEachSample("pTBB", "BDTCuts_Normed", "pT(bb)", "Events", BDTcuts, norm=True)

DrawPlotsWithCutsForEachSample("pTB1",  "BDTCuts_Normed", "p_{T} B1",  "Events", BDTcuts, norm=True)
DrawPlotsWithCutsForEachSample("pTB2",  "BDTCuts_Normed", "p_{T} B2",  "Events", BDTcuts, norm=True)

DrawPlotsWithCutsForEachSample("pTB1",  "BDTCuts", "p_{T} B1",  "Events", BDTcuts, norm=True)
DrawPlotsWithCutsForEachSample("pTB2",  "BDTCuts", "p_{T} B2",  "Events", BDTcuts, norm=True)


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

    #CalFtest( FITRESULTS[cut]["BernsteinO2"], FITRESULTS[cut]["BernsteinO3"] )
    #CalFtest( FITRESULTS[cut]["BernsteinO3"], FITRESULTS[cut]["BernsteinO4"] )
    #CalFtest( FITRESULTS[cut]["BernsteinO4"], FITRESULTS[cut]["BernsteinO5"] )

    #CalFtest( FITRESULTS[cut]["ExpoBernsteinO2"], FITRESULTS[cut]["ExpoBernsteinO3"] )
    #CalFtest( FITRESULTS[cut]["ExpoBernsteinO3"], FITRESULTS[cut]["ExpoBernsteinO4"] )
    #CalFtest( FITRESULTS[cut]["ExpoBernsteinO4"], FITRESULTS[cut]["ExpoBernsteinO5"] )

output.Close()
