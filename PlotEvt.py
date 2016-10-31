from ROOT import *
from math import sqrt
from optparse import OptionParser
import AtlasStyle as Atlas
import array
from HistoLab import HistoTool, CopyHist, StackHists, GetPlot, ScaleHist, Project2D, CopyHists, project1D, NormHists, AddHistPairs, ResetBins, GetHistMedian, GetHistIQR, GetRatios, GetAveRatios, GetGausMean, CalSysVariation, CalSysDiff
import numpy as np
import sys as sys
from copy import deepcopy

p = OptionParser()
p.add_option('--input',  type = 'string', default = "GbbProcessEvt.root", dest = 'input', help = 'intput File' )
p.add_option('--output',  type = 'string', default = "PlotEvt.root", dest = 'output', help = 'intput File' )
p.add_option('--XHHoutput',  type = 'string', default = "HHbbbbProcessEvt_mc01_MuInJet.root", dest = 'HHoutput', help = 'intput File' )
p.add_option('--HerwigOutput',  type = 'string', default = "GbbProcessEvt_Herwig_mc104_data303_wegithpTEta_TrkJetpT10_FatJetpT200_dPhiCut2.5_CJVeto0.75.root", dest = 'HerwigOutput', help = 'herwig intput File' )
p.add_option('--PythiaTruth',  type = 'string',  default = "GbbProcessEvt_mc305_data303_wegithpTEta_flavcorr_TrkJetpT10_FatJetpT200_dPhiCut2.5_CJVeto0.75.root", dest = 'PythiaTruth', help = 'Pythia Truth input File' )
p.add_option('--isHerwig',  type = 'string', default = "n", dest = 'isHerwig', help = 'herwig?' )

(o,a) = p.parse_args()

if o.isHerwig == 'n':
    isHerwig = False
else:
    isHerwig = True

infile = TFile(o.input, "READ")
XHHfile = TFile(o.HHoutput, "READ")
outdir = TFile(o.output, "Recreate")
Herwigfile = TFile(o.HerwigOutput, "READ")
PythiaTruthfile = TFile(o.PythiaTruth, "READ")
NoReweightFile = TFile("GbbProcessEvt_mc309_data307_wegithpTEta_TrkJetpT10_FatJetpT200_dPhiCut0.root","READ")
#NoReweightFile = TFile("GbbProcessEvt_mc309_data307_wegithpTEta_TrkJetpT10_FatJetpT200_dPhiCut0_MuonpT15.root","READ")
#NoReweightFile = TFile("GbbProcessEvt_mc306_data303_wegithpTEta_TrkJetpT10_FatJetpT200_dPhiCut2.5.root","READ")
#NoReweightFile = TFile("GbbProcessEvt_mc308_data305_wegithpTEta_TrkJetpT10_FatJetpT200_dPhiCut0_MuonpT8.root", "READ")
#NoReweightFile = TFile("GbbProcessEvt_mc271_data182_wegithpTEta_TrkJetpT10.root", "READ")
#NoReweightFile = TFile("GbbProcessEvt_mc268_data250_wegithpTEta_unsmeared.root", "READ")


pTSysFiles = []
for ifile in range(12):
    pTSysFiles.append(TFile("FitUncert/GbbProcessEvt_mc309_data307_wegithpTEta_flavcorr_TrkJetpT10_FatJetpT200_dPhiCut0_pTBin"+str(ifile+1)+".root"))


DrawTool = HistoTool()
DrawTool.SetFillStyle(1001)
DrawTool.SetCompData(True)

DrawTool_Ratio = HistoTool()
DrawTool_Ratio.SetFillStyle(1001)
DrawTool_Ratio.SetCompData(True)


DrawTool_DoubleRatio = HistoTool()
DrawTool_DoubleRatio.SetFillStyle(1001)
DrawTool_DoubleRatio.SetCompData(True)
DrawTool_DoubleRatio.colorlist = DrawTool_DoubleRatio.boringcolor
DrawTool_DoubleRatio.studytype = "Preliminary"

#DrawTool_New = HistoTool()
#DrawTool_New.SetFillStyle(1001)
#DrawTool_New.SetCompData(True)

#DrawTool.DoPrintPlots()
#DrawTool.SetOutPlotDir("Plots/")
#DrawTool.SetOutPlotDir("Plots_0Mu/")
#
#

NoFlav = ["L"]
SingleFlav = ["B", "C", "L"]
global DoubleFlav
DoubleFlav = ["BB", "BC", "BL", "CC", "CL","LL"]
OrderFlav = ["BB", "BC", "BL", "CB", "CC", "CL", "LB", "LC", "LL"]
Sd0Flav = ["BB", "BL", "LL", "CL", "CC"]

def StackFlav(name, flav, doBtag=False, doStack = True, thisfile = infile):
    HistList = []
    #print name, flav, doBtag, doStack, thisfile
    for i in range(len(flav)):
        if doBtag:
            thisPlot = CopyHist(thisfile.Get( name + "_" + flav[i] +"_btag_MC"))
        else:
            thisPlot = CopyHist(thisfile.Get( name + "_" + flav[i] +"_MC"))

        if thisPlot !=None:
            #thisPlot.Sumw2()
            HistList.append(thisPlot)

    #print len(HistList)
    Stack, Leg, Sum = StackHists(HistList, flav)
    if doStack:
        return Stack, Leg, Sum
    else:
        return HistList, flav, Sum

def PlotFits(plots, PlotName, axis, labels):
    newplots = []
    for i in range(len(plots)):
        
        newplots.append(deepcopy(plots[i]))
        GetGausMean(newplots[i])
        newplots[i].GetFunction('gaus').SetLineColor(DrawTool.colorlist[i])

    DrawTool.DrawHists(PlotName, axis, newplots, labels)

def getCorMatrixAndFitUncert(filename):

    txtfile = open(filename, 'r')
    lines = txtfile.readlines()
    partot =int(lines[1].split()[-1])

    #print 'total n par', partot
    sigma = []
    matrix = []
    corfac = []

    for i in range(partot):
        sigma.append( float(lines[partot+4+i].split()[8])/ float(lines[partot+4+i].split()[7])  )
        corfac.append(float(lines[partot+4+i].split()[4]))
        #print 'sigma', sigma
        thiscov = []
        for j in range(partot):
            thiscov.append( float(lines[2+i].split()[2+j]) )
            matrix.append(thiscov)
    return matrix, sigma, corfac


def ComputeSys(plot, upsys, dnsys, sysname):
    ## load data first
    upratios = []
    dnratios = []
    for i in range(len(sysname)):
        upratio = CopyHist(upsys[i])
        dnratio = CopyHist(dnsys[i])
        upratio.Divide(plot)
        dnratio.Divide(plot)
        upratios.append(upratio)
        dnratios.append(dnratio)
        
    upsum = 0
    dnsum = 0

    if upratios != []:
        upsum = CopyHist(upratios[0], True)
        for ibin in range(upsum.GetNbinsX()+1):
            thisbin = 0
            for isys in range(len(upratios)):
                if plot.GetBinContent(ibin+1) <1:
                    corr = 0
                else:
                    corr = (upratios[isys].GetBinContent(ibin+1)-1)**2
                #corr = (upratios[isys].GetBinContent(ibin+1)-1)**2
                thisbin += corr
            upsum.SetBinContent(ibin, 1+sqrt(thisbin))

    if dnratios != []:
        dnsum = CopyHist(dnratios[0], True)
        for ibin in range(dnsum.GetNbinsX()+1):
            thisbin = 0
            for isys in range(len(dnratios)):
                if plot.GetBinContent(ibin+1) <1:
                    corr = 0
                else:
                    corr = (dnratios[isys].GetBinContent(ibin+1)-1)**2
                #corr = (upratios[isys].GetBinContent(ibin+1)-1)**2
                thisbin += corr
            dnsum.SetBinContent(ibin, 1-sqrt(thisbin))
    
    sysComb = CopyHist(upsum)
    for ibin in range(sysComb.GetNbinsX()+1):
        sysComb.SetBinContent(ibin+1, 1)
        sysComb.SetBinError(ibin+1, (upsum.GetBinContent(ibin+1)-dnsum.GetBinContent(ibin+1))/2)
#    for ibin in range(sysComb.GetNbinsX()+1):
#        sysComb.SetBinContent(ibin+2, 1)
#        sysComb.SetBinError(ibin+2, (upsum.GetBinContent(ibin+1)-dnsum.GetBinContent(ibin+1)))
        #sysComb.SetBinError(9, 0.08)

    return sysComb


def ComputeSysRatio(plot, upsys, dnsys, sysname):
    ## load data first
    upratios = []
    dnratios = []
    for i in range(len(sysname)):
        upratio = CopyHist(upsys[i])
        dnratio = CopyHist(dnsys[i])
        upratio.Divide(plot)
        dnratio.Divide(plot)
        upratios.append(upratio)
        dnratios.append(dnratio)
        
    upsum = 0
    dnsum = 0

    if upratios != []:
        upsum = CopyHist(upratios[0], True)
        for ibin in range(upsum.GetNbinsX()+1):
            thisbin = 0
            for isys in range(len(upratios)):
                corr = (upratios[isys].GetBinContent(ibin+1)-1)**2
                
                thisbin += corr
            upsum.SetBinContent(ibin+1, sqrt(thisbin))

    if dnratios != []:
        dnsum = CopyHist(dnratios[0], True)
        for ibin in range(dnsum.GetNbinsX()+1):
            thisbin = 0
            for isys in range(len(dnratios)):
                corr = (dnratios[isys].GetBinContent(ibin+1)-1)**2
                thisbin += corr
            dnsum.SetBinContent(ibin+1, sqrt(thisbin))
    
    sysComb = CopyHist(plot)
    for ibin in range(sysComb.GetNbinsX()+1):
        sysComb.SetBinError(ibin+1, sqrt( ((upsum.GetBinContent(ibin+1)*sysComb.GetBinContent(ibin+1)+dnsum.GetBinContent(ibin+1)*sysComb.GetBinContent(ibin+1))/2)**2   +sysComb.GetBinError(ibin+1)**2)  )
        
    return sysComb
    

def ComputeFitSysSingleFlav(plot, sysplots, flav, fitbinuncert):
    # sysplots is a MCData
    # plot is a plot
    
    nbins = plot.GetNbinsX()
    newplot = CopyHist(plot)
    newplot_dn = CopyHist(plot)

    for ibin in range(nbins+1):
        totaluncert = 0

        # pT bin
        for i in range(len(sysplots)):
            # read in uncertainty for each flavor
            ptbinuncert = 0
            thisplot = 0

            for j in range(len(sysplots[i].Leg)):
                if sysplots[i].Leg[j] == flav:
                    thisplot =  sysplots[i].Stack[j]
                    break
                else:
                    continue
            
            thisbin = thisplot.GetBinContent(ibin+1)

            if flav == "BB":
                ptbinuncert += (thisbin * fitbinuncert[i][0])**2
            if flav == "BL":
                ptbinuncert += (thisbin * fitbinuncert[i][1])**2

            if flav == "LL":
                ptbinuncert += (thisbin * fitbinuncert[i][2])**2

            if flav == "BC":
                ptbinuncert += (thisbin * fitbinuncert[i][2])**2

            if flav == "CL":
                ptbinuncert += (thisbin * fitbinuncert[i][3])**2

            if flav == "CC":
                ptbinuncert += (thisbin * fitbinuncert[i][4])**2
            
            totaluncert += abs(ptbinuncert)

        if newplot.GetBinContent(ibin+1) != 0:
            totaluncert = sqrt(abs(totaluncert))/newplot.GetBinContent(ibin+1)
        else:
            totaluncert = 0

        newplot.SetBinContent(ibin+1, newplot.GetBinContent(ibin+1)*(1+totaluncert))
        newplot_dn.SetBinContent(ibin+1, newplot_dn.GetBinContent(ibin+1)*(1-totaluncert))
    return newplot, newplot_dn


def ComputeFitSys(plot, sysplots, cormatrices, fitbinuncert):
    # sysplots is the list of MCData class, so is plot
    # cor matrix element BB, BL, LL, CL, CC

    nbins = plot.MC.GetNbinsX()
    newplot = CopyHist(plot.MC)
    newplot_dn = CopyHist(plot.MC)

    for ibin in range(nbins+1):
        totaluncert = 0

        # pT bin
        for i in range(len(sysplots)):
            # read in uncertainty for each flavor
            flavuncert = [0]*5
            ptbinuncert = 0

            for j in range(len(DoubleFlav)):
                #print sysplots[i].Leg
                index = sysplots[i].Leg.index(DoubleFlav[j])

                thisplot =  sysplots[i].Stack[index]
                thisbin = thisplot.GetBinContent(ibin+1)

                if DoubleFlav[j] == "BB":
                    flavuncert[0] += (thisbin * fitbinuncert[i][0])**2

                if DoubleFlav[j] == "BL":
                    flavuncert[1] += (thisbin * fitbinuncert[i][1])**2

                if DoubleFlav[j] == "LL":
                    flavuncert[2] += (thisbin * fitbinuncert[i][2])**2

                if DoubleFlav[j] == "BC":
                    flavuncert[2] += (thisbin * fitbinuncert[i][2])**2

                if DoubleFlav[j] == "CL":
                    flavuncert[3] += (thisbin * fitbinuncert[i][3])**2

                if DoubleFlav[j] == "CC":
                    flavuncert[4] += (thisbin * fitbinuncert[i][4])**2
            
            for f1 in range(5):
                for f2 in range(f1, 5):
                    ptbinuncert += sqrt(flavuncert[f1])* sqrt(flavuncert[f2]) * cormatrices[i][f1][f2] 

            totaluncert += abs(ptbinuncert)

        if newplot.GetBinContent(ibin+1) != 0:
            totaluncert = sqrt(abs(totaluncert))/newplot.GetBinContent(ibin+1)
        else:
            totaluncert = 0

        newplot.SetBinContent(ibin+1, newplot.GetBinContent(ibin+1)*(1+totaluncert))
        newplot_dn.SetBinContent(ibin+1, newplot_dn.GetBinContent(ibin+1))

    ## b sys
    bnewplot = CopyHist(plot.BMC)
    bnewplot_dn = CopyHist(plot.BMC)

    for ibin in range(nbins+1):
        totaluncert = 0

        # pT bin
        for i in range(len(sysplots)):
            # read in uncertainty for each flavor
            flavuncert = [0]*5
            ptbinuncert = 0

            for j in range(len(DoubleFlav)):
                #print sysplots[i].Leg
                index = sysplots[i].Leg.index(DoubleFlav[j])

                thisplot =  sysplots[i].BStack[index]
                thisbin = thisplot.GetBinContent(ibin+1)

                if DoubleFlav[j] == "BB":
                    flavuncert[0] += (thisbin * fitbinuncert[i][0]/2)**2

                if DoubleFlav[j] == "BL":
                    flavuncert[1] += (thisbin * fitbinuncert[i][1]/2)**2

                if DoubleFlav[j] == "LL":
                    flavuncert[2] += (thisbin * fitbinuncert[i][2]/2)**2

                if DoubleFlav[j] == "BC":
                    flavuncert[2] += (thisbin * fitbinuncert[i][2]/2)**2

                if DoubleFlav[j] == "CL":
                    flavuncert[3] += (thisbin * fitbinuncert[i][3]/2)**2

                if DoubleFlav[j] == "CC":
                    flavuncert[4] += (thisbin * fitbinuncert[i][4]/2)**2
            
            for f1 in range(5):
                for f2 in range(f1, 5):
                    ptbinuncert += sqrt(flavuncert[f1])* sqrt(flavuncert[f2]) * cormatrices[i][f1][f2] 

            totaluncert += abs(ptbinuncert)

        if bnewplot.GetBinContent(ibin+1) != 0:
            totaluncert = sqrt(abs(totaluncert))/bnewplot.GetBinContent(ibin+1)
        else:
            totaluncert = 0

        bnewplot.SetBinContent(ibin+1, bnewplot.GetBinContent(ibin+1)*(1+totaluncert))
        bnewplot_dn.SetBinContent(ibin+1, bnewplot_dn.GetBinContent(ibin+1))

    return newplot, newplot_dn, bnewplot, bnewplot_dn


class MCData:
    def __init__(self, HistName, FlavType, CanvName, XLabel, YLabel, doBtag, noStack = False, thisfile = infile):
      
        #print thisfile
        #print HistName

        self.CanvName = CanvName
        self.XLabel = XLabel
        self.YLabel = YLabel
        self.FlavType = FlavType
        self.Stack, self.Leg, self.Sum = StackFlav(HistName, FlavType, False, True, thisfile)
        if noStack:
            self.Stack, self.Leg, self.Sum = StackFlav(HistName, FlavType, False, False, thisfile)
        self.Data = thisfile.Get(HistName+"_data")
        self.doBtag = doBtag
        self.BStack = 0
        self.BLeg = 0
        self.BSum = 0
        self.MC = self.Stack[-1]

        if self.doBtag:
            self.BStack, self.BLeg, self.BSum = StackFlav(HistName, FlavType, True, True, thisfile)
            self.BData = thisfile.Get(HistName+"_btag_data")
            self.BMC = self.BStack[-1]


    def Draw(self):
        DrawTool.DrawHists(self.CanvName,  [self.XLabel, self.YLabel], [], [], self.Stack+[self.Data], self.Leg +["Data"])
        if self.doBtag:
            DrawTool.DrawHists(self.CanvName+"_btag",  [self.XLabel, self.YLabel], [], [], self.BStack+[self.BData], self.BLeg +["Data"])
        return 

    def DrawMC(self):
        DrawTool.DrawHists(self.CanvName,  [self.XLabel, self.YLabel], [], [], self.Stack, self.Leg)
        if self.doBtag:
            DrawTool.DrawHists(self.CanvName+"_btag",  [self.XLabel, self.YLabel], [], [], self.BStack, self.BLeg)
        return 

    def DrawSys(self, sysname, sys, bsys):
        DrawTool.DrawHists(self.CanvName+"_sys_"+sysname,  [self.XLabel, self.YLabel], [], [], self.Stack+[self.Data], self.Leg +["Data"], [sys])
        if self.doBtag:
            DrawTool.DrawHists(self.CanvName+"_btag_sys_"+sysname,  [self.XLabel, self.YLabel], [], [], self.BStack+[self.BData], self.BLeg +["Data"], [bsys])
        return 

    def DrawSysMultiple(self, sysname, sys, bsys, bsysonly=None):
        DrawTool.DrawHists(self.CanvName+"_sys_"+sysname,  [self.XLabel, self.YLabel], [], [], self.Stack+[self.Data], self.Leg +["Data"], sys)
        if self.doBtag:
            DrawTool.DrawHists(self.CanvName+"_btag_sys_"+sysname,  [self.XLabel, self.YLabel], [], [], self.BStack+[self.BData], self.BLeg +["Data"], bsys, bsysonly)
        return 
    
    def DrawMCTot(self):
        DrawTool.DrawHists(self.CanvName,  [self.XLabel, self.YLabel], [self.Stack[-1],self.Data],  ["MC","Data"])
        return 

    def DrawMCOnly(self):
        DrawTool.DrawHists(self.CanvName,  [self.XLabel, self.YLabel], [self.Stack[-1]],  ["MC"])
        return 


    def DrawSinglePlot(self, i, postfix, isdata):
        if isdata:
            DrawTool.DrawHists(self.CanvName+postfix+"_data",  [self.XLabel, self.YLabel], [self.Data]  ["Data"])
        else:
            DrawTool.DrawHists(self.CanvName+postfix,  [self.XLabel, self.YLabel], [self.Stack[i]],  ["MC"])
        return 

    def ScaleMC(self, factor):
        for i in range(len(self.Stack)):
            if self.Stack[i] == None:
                continue
            self.Stack[i].Scale(factor)
        for i in range(len(self.BStack)):
            if self.BStack[i] == None:
                continue
            self.BStack[i].Scale(factor)

        return

    def NormStacks(self):
        self.Stack = NormHists(self.Stack)
        self.BStack = NormHists(self.BStack)
        self.MC = self.Stack[-1]
        self.BMC = self.BStack[-1]
        
def SetHist(Hist, vals, errors=None):
    for i in range(len(vals)):
        Hist.SetBinContent(i+1, vals[i])
        if errors != None:
            Hist.SetBinError(i+1, errors[i])

def PlotStat(RawPlots, Legends, Bins, Title, XLabel, YLabel, StatType = 'mean'):
    Hist = []
    for i in range(len(Legends)):
        thisStat = []
        thisError = []
        thisHist = TH1D(Legends[i]+Title, Legends[i]+Title, len(Bins)-1, Bins)
        for j in range(len(RawPlots[i])):
            if StatType == 'mean':
                thisStat.append(RawPlots[i][j].GetMean(1))
                thisError.append(RawPlots[i][j].GetMeanError(1))
            if StatType == 'rms':
                thisStat.append(RawPlots[i][j].GetRMS(1))
                thisError.append(RawPlots[i][j].GetRMSError(1))

            if StatType == 'gausmean':
                mean, error = GetGausMean(RawPlots[i][j])
                thisStat.append(mean)
                thisError.append(RawPlots[i][j].GetMeanError(1))
            
            if StatType == 'median':
                thisStat.append(GetHistMedian(RawPlots[i][j]))
                thisError.append(RawPlots[i][j].GetMeanError(1))

        SetHist(thisHist, thisStat, thisError)
        thisHist.SetMarkerStyle(26)
        if "Data" in Legends[i]:
            thisHist.SetMarkerStyle(24)
        if ("Pythia" in Legends[i]) and ("W/O FC" in Legends[i]):
            thisHist.SetMarkerStyle(25)
        if "Herwig" in Legends[i]:
            thisHist.SetMarkerStyle(32)
        Hist.append(thisHist)

    DrawTool.DrawHists(Title, [XLabel, YLabel], Hist, Legends)
    return Hist    

outdir.cd()


#### Get Flav Correlation Matrices
fileprefix = 'FitErrorCal/'
fitversion = 'v62_pT'
CorMatrices = []
FitBinUncert = []
CorFacs = []
for i in range(12):
    matrix, uncert, corfac = getCorMatrixAndFitUncert(fileprefix+fitversion+str(i+1)+'.txt')
    CorMatrices.append( matrix )
    FitBinUncert.append( uncert )
    CorFacs.append( corfac )

## Event Level Plots
H_TrigJet_pT    =  MCData("TrigJet_pT",     NoFlav,   "TrigJet_pT",     "Trig Jet p_{T} [GeV]",        "Jets / Bin",      False)
H_TrigJet_Eta   =  MCData("TrigJet_Eta",    NoFlav,   "TrigJet_Eta",    "Trig Jet #eta",             "Jets / 0.5",      False)
H_TrigJet_Phi   =  MCData("TrigJet_Phi",    NoFlav,   "TrigJet_Phi",    "Trig Jet #phi",             "Jets / Bin",      False)
H_CloseJet_Eta   =   MCData("CloseJet_Eta", NoFlav,   "CloseJet_Eta",     "Close Jet #eta",        "Jets / 0.5",      False)
H_CloseJet_Phi   =   MCData("CloseJet_Phi", NoFlav,   "CloseJet_Phi",     "Close Jet #phi",        "Jets / Bin",      False)
H_CloseJet_pT    =   MCData("CloseJet_pT", NoFlav,   "CloseJet_pT",     "Close Jet p_{T} [GeV]",        "Jets / Bin",      False)
H_dR_CloseJet_TrigJet   =  MCData("dR_CloseJet_TrigJet", NoFlav,   "dR_CloseJet_TrigJet",     "dR",        "Jets / Bin",      False)

H_NJet          =  MCData("NJet",           NoFlav,   "NJet(AntiKt4)",  "NJet(AntiKt4)",            "Entries",   False)
H_dR_TrigJet_MuonJet  =   MCData("dR_TrigJet_MuonJet", NoFlav,   "dR_TrigJet_MuJet",    "#DeltaR(Trig Jet, Muon Jet)",       "nJet",      False)

H_MuJet_pT      =  MCData("MuJet_pT",     SingleFlav,   "MuJet_pT",     "p_{T} [GeV]",        "Jets / Bin",      True)
H_dR_Mu_MuJet   =  MCData("dR_Mu_MuJet",  SingleFlav,   "dR_Mu_MuJet",  "#DeltaR(Muon, Muon Jet)",      "Jets / 0.04",      True)
H_Muon_pT       =  MCData("Mu_pT",        SingleFlav,   "Mu_pT",        "p_{T} [GeV]",       "Muons / 5 GeV",      True)
H_NonMuJet_pT   =  MCData("NonMuJet_pT",  SingleFlav,   "NonMuJet_pT",  "p_{T} [GeV]",     "Jets / Bin",      True)
H_LeadB_pT      =  MCData("LeadB_pT",     SingleFlav,   "LeadB_pT",     "Muon Jet B p_{T} [GeV]",        "Jets / Bin",      True)
H_SecB_pT       =  MCData("SecB_pT",      SingleFlav,   "SecB_pT",      "Non-Muon Jet B p_{T} [GeV]",     "Jets / Bin",      True)
H_ThJet_pT      =  MCData("ThJet_pT",     SingleFlav,   "ThJet_pT",     "Third Jet p_{T} [GeV]",     "Jets / Bin",      True)

H_MuJet_Sd0     =  MCData("MuJet_Sd0",     Sd0Flav,   "MuJet_Sd0",     "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0  =  MCData("NonMuJet_Sd0",  Sd0Flav,   "NonMuJet_Sd0",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)

H_MuJet_Sd0_pT1 =  MCData("MuJet_Sd0_pT1", Sd0Flav,   "MuJet_Sd0_pT1", "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0_pT1  =  MCData("NonMuJet_Sd0_pT1",  Sd0Flav,   "NonMuJet_Sd0_pT1",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)
H_MuJet_Sd0_pT2 =  MCData("MuJet_Sd0_pT2", Sd0Flav,   "MuJet_Sd0_pT2", "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0_pT2  =  MCData("NonMuJet_Sd0_pT2",  Sd0Flav,   "NonMuJet_Sd0_pT2",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)
H_MuJet_Sd0_pT3 =  MCData("MuJet_Sd0_pT3", Sd0Flav,   "MuJet_Sd0_pT3", "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0_pT3  =  MCData("NonMuJet_Sd0_pT3",  Sd0Flav,   "NonMuJet_Sd0_pT3",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)
H_MuJet_Sd0_pT4 =  MCData("MuJet_Sd0_pT4", Sd0Flav,   "MuJet_Sd0_pT4", "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0_pT4  =  MCData("NonMuJet_Sd0_pT4",  Sd0Flav,   "NonMuJet_Sd0_pT4",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)
H_MuJet_Sd0_pT5 =  MCData("MuJet_Sd0_pT5", Sd0Flav,   "MuJet_Sd0_pT5", "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0_pT5  =  MCData("NonMuJet_Sd0_pT5",  Sd0Flav,   "NonMuJet_Sd0_pT5",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)
H_MuJet_Sd0_pT6 =  MCData("MuJet_Sd0_pT6", Sd0Flav,   "MuJet_Sd0_pT6", "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0_pT6  =  MCData("NonMuJet_Sd0_pT6",  Sd0Flav,   "NonMuJet_Sd0_pT6",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)
H_MuJet_Sd0_pT7 =  MCData("MuJet_Sd0_pT7", Sd0Flav,   "MuJet_Sd0_pT7", "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0_pT7  =  MCData("NonMuJet_Sd0_pT7",  Sd0Flav,   "NonMuJet_Sd0_pT7",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)
H_MuJet_Sd0_pT8 =  MCData("MuJet_Sd0_pT8", Sd0Flav,   "MuJet_Sd0_pT8", "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0_pT8  =  MCData("NonMuJet_Sd0_pT8",  Sd0Flav,   "NonMuJet_Sd0_pT8",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)
H_MuJet_Sd0_pT9 =  MCData("MuJet_Sd0_pT9", Sd0Flav,   "MuJet_Sd0_pT9", "Muon Jet S_{d0}",        "Jets / Bin",      True)
H_NonMuJet_Sd0_pT9  =  MCData("NonMuJet_Sd0_pT9",  Sd0Flav,   "NonMuJet_Sd0_pT9",  "NonMuon Jet S_{d0}",     "Jets / Bin",      True)

H_MuJet_LL_pT     =  MCData("MuJet_LL_pT",     SingleFlav,   "MuJet_LL_pT",     "Muon Jet p_{T} [GeV]",        "Jets / Bin",      True)
H_NonMuJet_LL_pT  =  MCData("NonMuJet_LL_pT",  SingleFlav,   "NonMuJet_LL_pT",  "NonMuon Jet p_{T} [GeV]",     "Jets / Bin",      True)
H_MuJet_NStack_pT     =  MCData("MuJet_pT",     SingleFlav,   "MuJet_pT",     "Muon Jet p_{T} [GeV]",        "Jets / Bin",  True, True)
H_NonMuJet_NStack_pT  =  MCData("NonMuJet_pT",  SingleFlav,   "NonMuJet_pT",  "NonMuon Jet p_{T} [GeV]",     "Jets / Bin",  True, True)

H_MuJet_Eta    =  MCData("MuJet_Eta",    SingleFlav,   "MuJet_Eta",    "Muon Jet #eta",        "Jets / 0.5",      True)
H_NonMuJet_Eta =  MCData("NonMuJet_Eta", SingleFlav,   "NonMuJet_Eta", "NonMuon Jet #eta",     "Jets / 0.5",      True)
H_MuJet_MV1    =  MCData("MuJet_MV1",    SingleFlav,   "MuJet_MV1",    "Muon Jet MV1",        "Jets / Bin",      True)
H_NonMuJet_MV1 =  MCData("NonMuJet_MV1", SingleFlav,   "NonMuJet_MV1", "NonMuon Jet MV1",     "Jets / Bin",      True)
H_MuJet_MV1_Tail =  MCData("MuJet_MV1_Tail",    SingleFlav,   "MuJet_MV1_Tail",    "Muon Jet MV1",        "Jets / Bin",      True)
H_NonMuJet_MV1_Tail =  MCData("NonMuJet_MV1_Tail", SingleFlav,   "NonMuJet_MV1_Tail", "NonMuon Jet MV1",     "Jets / Bin",      True)

H_dR           =  MCData("SubJet_dR",        DoubleFlav,   "SubJet_dR",    "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1",   True)
H_dR_JESUp     =  MCData("SubJet_dR_JESUp",  DoubleFlav,   "SubJet_dR_JESUp",    "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1",   True)
H_dR_JESDn     =  MCData("SubJet_dR_JESDn",  DoubleFlav,   "SubJet_dR_JESDn",    "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1",   True)
H_dR_BUp       =  MCData("SubJet_dR_BUp",    DoubleFlav,   "SubJet_dR_BUp",    "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1",   True)
H_dR_BDn       =  MCData("SubJet_dR_BDn",    DoubleFlav,   "SubJet_dR_BDn",    "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1",   True)
H_dR_CUp       =  MCData("SubJet_dR_CUp",    DoubleFlav,   "SubJet_dR_CUp",    "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1",   True)
H_dR_CDn       =  MCData("SubJet_dR_CDn",    DoubleFlav,   "SubJet_dR_CDn",    "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1",   True)
H_dR_LUp       =  MCData("SubJet_dR_LUp",    DoubleFlav,   "SubJet_dR_LUp",    "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1",   True)
H_dR_LDn       =  MCData("SubJet_dR_LDn",    DoubleFlav,   "SubJet_dR_LDn",    "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1",   True)
H_dR_13        =  MCData("SubJet_dR_13",     DoubleFlav,   "SubJet_dR_13",    "#DeltaR(Jet 1, Jet 3)", "Entries / 0.1",   True)
H_dR_23        =  MCData("SubJet_dR_23",     DoubleFlav,   "SubJet_dR_23",    "#DeltaR(Jet 2, Jet 3)", "Entries / 0.1",   True)
H_dR_LeadB_Jet3=  MCData("dR_LeadB_Jet3",    DoubleFlav,   "dR_LeadB_Jet3",   "#DeltaR(LeadB, Jet 3)", "Entries / 0.1",   True)
H_dR_SecB_Jet3 =  MCData("dR_SecB_Jet3",     DoubleFlav,   "dR_SecB_Jet3",    "#DeltaR(SecB, Jet 3)", "Entries / 0.1",   True)

H_FatJet_pT    =  MCData("FatJet_pT",    DoubleFlav,   "FatJet_pT",    "p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_FatJet_pT_NoStack    =  MCData("FatJet_pT",    DoubleFlav,   "FatJet_pT",    "p_{T} [GeV]",   "Jets / 50 GeV",   True, True)
H_FatJet_pT_NoReWeight   =  MCData("FatJet_pT",    DoubleFlav,   "FatJet_pT",    "p_{T} [GeV]",   "Jets / 50 GeV",   
                                   True, False, NoReweightFile)
H_FatJet_pT_NoReWeight_NoStack   =  MCData("FatJet_pT",    DoubleFlav,   "FatJet_pT",    "p_{T} [GeV]",   "Jets / 50 GeV",   
                                           True, True, NoReweightFile)

H_TrigJet_pT_Weight    =  MCData("TrigJet_pT_Weight",    DoubleFlav,   "TrigJet_pT_Weight",    "Trig Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_TrigJet_pT_Weight_JESUp  =  MCData("TrigJet_pT_Weight_JESUp",    DoubleFlav,   "TrigJet_pT_Weight_JESUp",    "Trig Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_TrigJet_pT_Weight_JESDn  =  MCData("TrigJet_pT_Weight_JESDn",    DoubleFlav,   "TrigJet_pT_Weight_JESDn",    "Trig Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_TrigJet_pT_Weight_BUp  =  MCData("TrigJet_pT_Weight_BUp",    DoubleFlav,   "TrigJet_pT_Weight_BUp",    "Trig Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_TrigJet_pT_Weight_BDn  =  MCData("TrigJet_pT_Weight_BDn",    DoubleFlav,   "TrigJet_pT_Weight_BDn",    "Trig Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)


H_FatJet_pT_JESUncert_Up  =  MCData("FatJet_pT_JESUncert_Up",    DoubleFlav,   "FatJet_pT_JESUncert_Up",    "Fat Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_FatJet_pT_JESUncert_Dn  =  MCData("FatJet_pT_JESUncert_Dn",    DoubleFlav,   "FatJet_pT_JESUncert_Dn",    "Fat Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_FatJet_pT_BtagUncert_Up =  MCData("FatJet_pT_BtagUncert_Up",    DoubleFlav,   "FatJet_pT_BtagUncert_Up",    "Fat Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_FatJet_pT_BtagUncert_Dn =  MCData("FatJet_pT_BtagUncert_Dn",    DoubleFlav,   "FatJet_pT_BtagUncert_Dn",    "Fat Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_FatJet_pT_BUp =  MCData("FatJet_pT_BUp",    DoubleFlav,   "FatJet_pT",    "Fat Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_FatJet_pT_BDn =  MCData("FatJet_pT_BDn",    DoubleFlav,   "FatJet_pT",    "Fat Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_FatJet_Eta_BUp =  MCData("FatJet_Eta_BUp",    DoubleFlav,   "FatJet_Eta",    "#eta",   "Jets / 0.5",   True)
H_FatJet_Eta_BDn =  MCData("FatJet_Eta_BDn",    DoubleFlav,   "FatJet_Eta",    "#eta",   "Jets / 0.5",   True)
H_FatJet_Eta_JESUp =  MCData("FatJet_Eta_JESUp",    DoubleFlav,   "FatJet_Eta",    "#eta",   "Jets / 0.5",   True)
H_FatJet_Eta_JESDn =  MCData("FatJet_Eta_JESDn",    DoubleFlav,   "FatJet_Eta",    "#eta",   "Jets / 0.5",   True)
H_FatJet_Eta   =  MCData("FatJet_Eta",   DoubleFlav,   "FatJet_Eta",   "#eta",  "Jets / 0.5",   True)
H_FatJet_Mass  =  MCData("FatJet_Mass",  DoubleFlav,   "FatJet_Mass",  "Mass [GeV]", "Jets / 25 GeV",   True)
H_FatJet_Mass_pT100  =  MCData("FatJet_Mass_pT100",  DoubleFlav,   "FatJet_Mass_pT100",  "Fat Jet Mass [GeV]", "Jets / 25 GeV",   True)
H_FatJet_Mass_pT300  =  MCData("FatJet_Mass_pT300",  DoubleFlav,   "FatJet_Mass_pT300",  "Fat Jet Mass [GeV]", "Jets / 25 GeV",   True)
H_FatJet_Mass_pT450  =  MCData("FatJet_Mass_pT450",  DoubleFlav,   "FatJet_Mass_pT450",  "Fat Jet Mass [GeV]", "Jets / 25 GeV",   True)
H_FatJet_Mass_BUp  =  MCData("FatJet_Mass_BUp",  DoubleFlav,   "FatJet_Mass",  "Fat Jet Mass [GeV]", "Jets / 25 GeV",   True)
H_FatJet_Mass_BDn  =  MCData("FatJet_Mass_BDn",  DoubleFlav,   "FatJet_Mass",  "Fat Jet Mass [GeV]", "Jets / 25 GeV",   True)
H_FatJet_Mass_JESUp  =  MCData("FatJet_Mass_JESUp",  DoubleFlav,   "FatJet_Mass",  "Fat Jet Mass [GeV]", "Jets / 25 GeV",   True)
H_FatJet_Mass_JESDn  =  MCData("FatJet_Mass_JESDn",  DoubleFlav,   "FatJet_Mass",  "Fat Jet Mass [GeV]", "Jets / 25 GeV",   True)
H_FatJet_JESUp =  MCData("FatJet_JESUp", DoubleFlav,   "FatJet_JESUp", "Fat Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_FatJet_JESDn =  MCData("FatJet_JESDn", DoubleFlav,   "FatJet_JESDn", "Fat Jet p_{T} [GeV]",   "Jets / 50 GeV",   True)
H_FatJet_JMSUp =  MCData("FatJet_JMSUp", DoubleFlav,   "FatJet_JMSUp", "Fat Jet Mass [GeV]",   "Jets / 25 GeV",   True)
H_FatJet_JMSDn =  MCData("FatJet_JMSDn", DoubleFlav,   "FatJet_JMSDn", "Fat Jet Mass [GeV]",   "Jets / 25 GeV",   True)

H_FatJet_MOPT         =  MCData("FatJet_MOPT", DoubleFlav,   "FatJet_MOPT", "Fat Jet Mass / p_{T}",   "Jets / Bin",   True)
H_FatJet_TrigJet_dPhi =  MCData("FatJet_TrigJet_dPhi", DoubleFlav, "FatJet_TrigJet_dPhi ", "#Delta#phi(FatJet, TrigJet)",   "Jets / Bin",   True)
H_FatJet_pTRef =   MCData("FatJet_pTRef", DoubleFlav, "FatJet_pTRef", "TrigJet p_{T}*Cos(#Detla R(FatJet, TrigJet))", "Jets / Bin", True)


H_FatJet_pTRef_Response  =  MCData("FatJet_pTRef_Response", DoubleFlav,   "FatJet_pTRef_Response", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_Herwig =  MCData("FatJet_pTRef_Response",  DoubleFlav,   "FatJet_pTRef_Response",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, Herwigfile)

H_FatJet_pTRef_Response_Trig  =  MCData("FatJet_pTRef_Response_Trig", DoubleFlav,   "FatJet_pTRef_Response_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_Herwig_Trig =  MCData("FatJet_pTRef_Response_Trig",  DoubleFlav,   "FatJet_pTRef_Response_Trig",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, Herwigfile)

H_FatJet_pTRef_Response_pT1  =  MCData("FatJet_pTRef_Response_pT1", DoubleFlav,   "FatJet_pTRef_Response_pT1", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2  =  MCData("FatJet_pTRef_Response_pT2", DoubleFlav,   "FatJet_pTRef_Response_pT2", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3  =  MCData("FatJet_pTRef_Response_pT3", DoubleFlav,   "FatJet_pTRef_Response_pT3", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_TJESUp  =  MCData("FatJet_pTRef_Response_pT1_TJESUp", DoubleFlav,   "FatJet_pTRef_Response_pT1_TJESUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_TJESUp  =  MCData("FatJet_pTRef_Response_pT2_TJESUp", DoubleFlav,   "FatJet_pTRef_Response_pT2_TJESUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_TJESUp  =  MCData("FatJet_pTRef_Response_pT3_TJESUp", DoubleFlav,   "FatJet_pTRef_Response_pT3_TJESUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_TJESDn  =  MCData("FatJet_pTRef_Response_pT1_TJESDn", DoubleFlav,   "FatJet_pTRef_Response_pT1_TJESDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_TJESDn  =  MCData("FatJet_pTRef_Response_pT2_TJESDn", DoubleFlav,   "FatJet_pTRef_Response_pT2_TJESDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_TJESDn  =  MCData("FatJet_pTRef_Response_pT3_TJESDn", DoubleFlav,   "FatJet_pTRef_Response_pT3_TJESDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_FJESUp  =  MCData("FatJet_pTRef_Response_pT1_FJESUp", DoubleFlav,   "FatJet_pTRef_Response_pT1_FJESUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_FJESUp  =  MCData("FatJet_pTRef_Response_pT2_FJESUp", DoubleFlav,   "FatJet_pTRef_Response_pT2_FJESUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_FJESUp  =  MCData("FatJet_pTRef_Response_pT3_FJESUp", DoubleFlav,   "FatJet_pTRef_Response_pT3_FJESUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_FJESDn  =  MCData("FatJet_pTRef_Response_pT1_FJESDn", DoubleFlav,   "FatJet_pTRef_Response_pT1_FJESDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_FJESDn  =  MCData("FatJet_pTRef_Response_pT2_FJESDn", DoubleFlav,   "FatJet_pTRef_Response_pT2_FJESDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_FJESDn  =  MCData("FatJet_pTRef_Response_pT3_FJESDn", DoubleFlav,   "FatJet_pTRef_Response_pT3_FJESDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_BUp  =  MCData("FatJet_pTRef_Response_pT1_BUp", DoubleFlav,   "FatJet_pTRef_Response_pT1_BUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_BUp  =  MCData("FatJet_pTRef_Response_pT2_BUp", DoubleFlav,   "FatJet_pTRef_Response_pT2_BUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_BUp  =  MCData("FatJet_pTRef_Response_pT3_BUp", DoubleFlav,   "FatJet_pTRef_Response_pT3_BUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_BDn  =  MCData("FatJet_pTRef_Response_pT1_BDn", DoubleFlav,   "FatJet_pTRef_Response_pT1_BDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_BDn  =  MCData("FatJet_pTRef_Response_pT2_BDn", DoubleFlav,   "FatJet_pTRef_Response_pT2_BDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_BDn  =  MCData("FatJet_pTRef_Response_pT3_BDn", DoubleFlav,   "FatJet_pTRef_Response_pT3_BDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)

H_FatJet_pTRef_Response_pT1_Trig  =  MCData("FatJet_pTRef_Response_pT1_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT1_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_Trig  =  MCData("FatJet_pTRef_Response_pT2_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT2_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_Trig  =  MCData("FatJet_pTRef_Response_pT3_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT3_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_TJESUp_Trig  =  MCData("FatJet_pTRef_Response_pT1_TJESUp_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT1_TJESUp_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_TJESUp_Trig  =  MCData("FatJet_pTRef_Response_pT2_TJESUp_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT2_TJESUp_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_TJESUp_Trig  =  MCData("FatJet_pTRef_Response_pT3_TJESUp_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT3_TJESUp_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_TJESDn_Trig  =  MCData("FatJet_pTRef_Response_pT1_TJESDn_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT1_TJESDn_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_TJESDn_Trig  =  MCData("FatJet_pTRef_Response_pT2_TJESDn_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT2_TJESDn_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_TJESDn_Trig  =  MCData("FatJet_pTRef_Response_pT3_TJESDn_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT3_TJESDn_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_FJESUp_Trig  =  MCData("FatJet_pTRef_Response_pT1_FJESUp_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT1_FJESUp_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_FJESUp_Trig  =  MCData("FatJet_pTRef_Response_pT2_FJESUp_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT2_FJESUp_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_FJESUp_Trig  =  MCData("FatJet_pTRef_Response_pT3_FJESUp_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT3_FJESUp_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_FJESDn_Trig  =  MCData("FatJet_pTRef_Response_pT1_FJESDn_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT1_FJESDn_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_FJESDn_Trig  =  MCData("FatJet_pTRef_Response_pT2_FJESDn_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT2_FJESDn_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_FJESDn_Trig  =  MCData("FatJet_pTRef_Response_pT3_FJESDn_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT3_FJESDn_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_BUp_Trig  =  MCData("FatJet_pTRef_Response_pT1_BUp_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT1_BUp_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_BUp_Trig  =  MCData("FatJet_pTRef_Response_pT2_BUp_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT2_BUp_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_BUp_Trig  =  MCData("FatJet_pTRef_Response_pT3_BUp_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT3_BUp_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT1_BDn_Trig  =  MCData("FatJet_pTRef_Response_pT1_BDn_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT1_BDn_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT2_BDn_Trig  =  MCData("FatJet_pTRef_Response_pT2_BDn_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT2_BDn_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_pT3_BDn_Trig  =  MCData("FatJet_pTRef_Response_pT3_BDn_Trig", DoubleFlav,   "FatJet_pTRef_Response_pT3_BDn_Trig", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)

H_FatJet_pTRef_Response_Herwig_pT1 =  MCData("FatJet_pTRef_Response_pT1",  DoubleFlav,   "FatJet_pTRef_Response_pT1",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, Herwigfile)
H_FatJet_pTRef_Response_Herwig_pT2 =  MCData("FatJet_pTRef_Response_pT2",  DoubleFlav,   "FatJet_pTRef_Response_pT2",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, Herwigfile)
H_FatJet_pTRef_Response_Herwig_pT3 =  MCData("FatJet_pTRef_Response_pT3",  DoubleFlav,   "FatJet_pTRef_Response_pT3",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, Herwigfile)

H_FatJet_pTRef_Response_Herwig_pT1_Trig =  MCData("FatJet_pTRef_Response_pT1_Trig",  DoubleFlav,   "FatJet_pTRef_Response_pT1_Trig",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, Herwigfile)
H_FatJet_pTRef_Response_Herwig_pT2_Trig =  MCData("FatJet_pTRef_Response_pT2_Trig",  DoubleFlav,   "FatJet_pTRef_Response_pT2_Trig",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, Herwigfile)
H_FatJet_pTRef_Response_Herwig_pT3_Trig =  MCData("FatJet_pTRef_Response_pT3_Trig",  DoubleFlav,   "FatJet_pTRef_Response_pT3_Trig",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, Herwigfile)


H_FatJet_pTRef_Response_Truth_pT1 =  MCData("FatJet_pTRef_Response_pT1",  DoubleFlav,   "FatJet_pTRef_Response_pT1",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, PythiaTruthfile)
H_FatJet_pTRef_Response_Truth_pT2 =  MCData("FatJet_pTRef_Response_pT2",  DoubleFlav,   "FatJet_pTRef_Response_pT2",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, PythiaTruthfile)
H_FatJet_pTRef_Response_Truth_pT3 =  MCData("FatJet_pTRef_Response_pT3",  DoubleFlav,   "FatJet_pTRef_Response_pT3",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, PythiaTruthfile)

H_FatJet_pTRef_Response_Truth_pT1_Trig =  MCData("FatJet_pTRef_Response_pT1_Trig",  DoubleFlav,   "FatJet_pTRef_Response_pT1_Trig",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, PythiaTruthfile)
H_FatJet_pTRef_Response_Truth_pT2_Trig =  MCData("FatJet_pTRef_Response_pT2_Trig",  DoubleFlav,   "FatJet_pTRef_Response_pT2_Trig",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, PythiaTruthfile)
H_FatJet_pTRef_Response_Truth_pT3_Trig =  MCData("FatJet_pTRef_Response_pT3_Trig",  DoubleFlav,   "FatJet_pTRef_Response_pT3_Trig",  "Fat Jet p_{T}/ Ref p_{T}", "Jets / Bin",   True, False, PythiaTruthfile)

H_FatJet_pTRef_Response_JESUp  =  MCData("FatJet_pTRef_Response_JESUp", DoubleFlav,   "FatJet_pTRef_Response_JESUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_JESDn  =  MCData("FatJet_pTRef_Response_JESDn", DoubleFlav,   "FatJet_pTRef_Response_JESDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_BUp  =  MCData("FatJet_pTRef_Response_BUp", DoubleFlav,   "FatJet_pTRef_Response_BUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_BDn  =  MCData("FatJet_pTRef_Response_BDn", DoubleFlav,   "FatJet_pTRef_Response_BDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)

H_FatJet_pTRef_Response_Trig_JESUp  =  MCData("FatJet_pTRef_Response_Trig_JESUp", DoubleFlav,   "FatJet_pTRef_Response_Trig_JESUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_Trig_JESDn  =  MCData("FatJet_pTRef_Response_Trig_JESDn", DoubleFlav,   "FatJet_pTRef_Response_Trig_JESDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_Trig_BUp  =  MCData("FatJet_pTRef_Response_Trig_BUp", DoubleFlav,   "FatJet_pTRef_Response_Trig_BUp", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)
H_FatJet_pTRef_Response_Trig_BDn  =  MCData("FatJet_pTRef_Response_Trig_BDn", DoubleFlav,   "FatJet_pTRef_Response_Trig_BDn", "Fat Jet p_{T}/ Ref p_{T}",   "Jets / Bin",   True)


H_FatJet_NJet      =  MCData("FatJet_NJet",  DoubleFlav,   "FatJet_NJet",  "NJet Matched", "Number of Events",   True)
H_FatJet_NCluster  =  MCData("FatJet_NCluster",  DoubleFlav,   "FatJet_NCluster",  "Number of Topo Cluster", "Number of Events",   True)

#### Mass Study
H_FatJet_TruthMass  =  MCData("FatJet_TruthMass",  DoubleFlav,   "FatJet_TruthMass",  "Fat Jet Truth Mass [GeV]", "Jets / 25 GeV",   True)
H_FatJet_TrackMass  =  MCData("FatJet_TrackMass",  DoubleFlav,   "FatJet_TrackMass",  "Fat Jet Track Mass [GeV]", "Jets / 25 GeV",   True)

### Response
H_FatJet_Mass_Response_NStack    =  MCData("FatJet_Mass_Response",  DoubleFlav,   "FatJet_Mass_Response_NStack",  "FatJet Response", "Jets / Bin",   True, True)
H_FatJet_Mass_Response_NStack_pT1_Eta1    =  MCData("FatJet_Mass_Response_pT1_Eta1",  DoubleFlav,   "FatJet_Mass_Response_NStack_pT1_Eta1",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_pT2_Eta1    =  MCData("FatJet_Mass_Response_pT2_Eta1",  DoubleFlav,   "FatJet_Mass_Response_NStack_pT2_Eta1",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_pT3_Eta1    =  MCData("FatJet_Mass_Response_pT3_Eta1",  DoubleFlav,   "FatJet_Mass_Response_NStack_pT3_Eta1",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_pT1_Eta2    =  MCData("FatJet_Mass_Response_pT1_Eta2",  DoubleFlav,   "FatJet_Mass_Response_NStack_pT1_Eta2",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_pT2_Eta2    =  MCData("FatJet_Mass_Response_pT2_Eta2",  DoubleFlav,   "FatJet_Mass_Response_NStack_pT2_Eta2",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_pT3_Eta2    =  MCData("FatJet_Mass_Response_pT3_Eta2",  DoubleFlav,   "FatJet_Mass_Response_NStack_pT3_Eta2",  "FatJet Response", "Jets / Bin",   True)

H_FatJet_Mass_Response_C_NStack_pT1_Eta1    =  MCData("FatJet_Mass_Response_C_pT1_Eta1",  DoubleFlav,   "FatJet_Mass_Response_C_NStack_pT1_Eta1",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_C_NStack_pT2_Eta1    =  MCData("FatJet_Mass_Response_C_pT2_Eta1",  DoubleFlav,   "FatJet_Mass_Response_C_NStack_pT2_Eta1",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_C_NStack_pT3_Eta1    =  MCData("FatJet_Mass_Response_C_pT3_Eta1",  DoubleFlav,   "FatJet_Mass_Response_C_NStack_pT3_Eta1",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_C_NStack_pT1_Eta2    =  MCData("FatJet_Mass_Response_C_pT1_Eta2",  DoubleFlav,   "FatJet_Mass_Response_C_NStack_pT1_Eta2",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_C_NStack_pT2_Eta2    =  MCData("FatJet_Mass_Response_C_pT2_Eta2",  DoubleFlav,   "FatJet_Mass_Response_C_NStack_pT2_Eta2",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_C_NStack_pT3_Eta2    =  MCData("FatJet_Mass_Response_C_pT3_Eta2",  DoubleFlav,   "FatJet_Mass_Response_C_NStack_pT3_Eta2",  "FatJet Response", "Jets / Bin",   True)

H_FatJet_Mass_Response_NStack_Eta1    =  MCData("FatJet_Mass_Response_Eta1",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta1",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_Eta2    =  MCData("FatJet_Mass_Response_Eta2",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta2",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_Eta3    =  MCData("FatJet_Mass_Response_Eta3",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta3",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_Eta4    =  MCData("FatJet_Mass_Response_Eta4",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta4",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_Eta5    =  MCData("FatJet_Mass_Response_Eta5",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta5",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_Eta6    =  MCData("FatJet_Mass_Response_Eta6",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta6",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_Eta7    =  MCData("FatJet_Mass_Response_Eta7",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta7",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_Eta8    =  MCData("FatJet_Mass_Response_Eta8",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta8",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_Eta9    =  MCData("FatJet_Mass_Response_Eta9",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta9",  "FatJet Response", "Jets / Bin",   True)
H_FatJet_Mass_Response_NStack_Eta10    =  MCData("FatJet_Mass_Response_Eta10",  DoubleFlav,   "FatJet_Mass_Response_NStack_Eta10",  "FatJet Response", "Jets / Bin",   True)

H_FatJet_Mass_Response_NStack.NormStacks()
H_FatJet_Mass_Response_NStack_pT1_Eta1.NormStacks()
H_FatJet_Mass_Response_NStack_pT2_Eta1.NormStacks()
H_FatJet_Mass_Response_NStack_pT3_Eta1.NormStacks()
H_FatJet_Mass_Response_NStack_pT1_Eta2.NormStacks()
H_FatJet_Mass_Response_NStack_pT2_Eta2.NormStacks()
H_FatJet_Mass_Response_NStack_pT3_Eta2.NormStacks()

H_FatJet_Mass_Response_C_NStack_pT1_Eta1.NormStacks()
H_FatJet_Mass_Response_C_NStack_pT2_Eta1.NormStacks()
H_FatJet_Mass_Response_C_NStack_pT3_Eta1.NormStacks()
H_FatJet_Mass_Response_C_NStack_pT1_Eta2.NormStacks()
H_FatJet_Mass_Response_C_NStack_pT2_Eta2.NormStacks()
H_FatJet_Mass_Response_C_NStack_pT3_Eta2.NormStacks()

H_FatJet_Mass_Response_NStack_Eta1.NormStacks()
H_FatJet_Mass_Response_NStack_Eta2.NormStacks()
H_FatJet_Mass_Response_NStack_Eta3.NormStacks()
H_FatJet_Mass_Response_NStack_Eta4.NormStacks()
H_FatJet_Mass_Response_NStack_Eta5.NormStacks()
H_FatJet_Mass_Response_NStack_Eta6.NormStacks()
H_FatJet_Mass_Response_NStack_Eta7.NormStacks()
H_FatJet_Mass_Response_NStack_Eta8.NormStacks()
H_FatJet_Mass_Response_NStack_Eta9.NormStacks()
H_FatJet_Mass_Response_NStack_Eta10.NormStacks()

H_FatJet_Mass_Response_HHbbbb_NStack    =  MCData("FatJet_Mass_Response",  DoubleFlav,   "FatJet_Mass_Response_HHbbbb_NStack",  "FatJet Response", "Jets / Bin",   True, True, XHHfile)
H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta1    =  MCData("FatJet_Mass_Response_pT1_Eta1",  DoubleFlav,   "FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta1",  "FatJet Response", "Jets / Bin",   True, False, XHHfile)
H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta1    =  MCData("FatJet_Mass_Response_pT2_Eta1",  DoubleFlav,   "FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta1",  "FatJet Response", "Jets / Bin",   True, False, XHHfile)
H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta1    =  MCData("FatJet_Mass_Response_pT3_Eta1",  DoubleFlav,   "FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta1",  "FatJet Response", "Jets / Bin",   True, False, XHHfile)
H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta2    =  MCData("FatJet_Mass_Response_pT1_Eta2",  DoubleFlav,   "FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta2",  "FatJet Response", "Jets / Bin",   True, False, XHHfile)
H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta2    =  MCData("FatJet_Mass_Response_pT2_Eta2",  DoubleFlav,   "FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta2",  "FatJet Response", "Jets / Bin",   True, False, XHHfile)
H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta2    =  MCData("FatJet_Mass_Response_pT3_Eta2",  DoubleFlav,   "FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta2",  "FatJet Response", "Jets / Bin",   True, False, XHHfile)

H_FatJet_Mass_Response_HHbbbb_NStack.NormStacks()
H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta1.NormStacks()
H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta1.NormStacks()
H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta1.NormStacks()
H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta2.NormStacks()
H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta2.NormStacks()
H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta2.NormStacks()


#### Track Mass Ratio
H_FatJet_TrackMass_Ratio  =  MCData("FatJet_TrackMass_Ratio",  DoubleFlav,   "FatJet_TrackMass_Ratio",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_BUp  =  MCData("FatJet_TrackMass_Ratio_BUp",  DoubleFlav,   "FatJet_TrackMass_Ratio_BUp",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_BDn  =  MCData("FatJet_TrackMass_Ratio_BDn",  DoubleFlav,   "FatJet_TrackMass_Ratio_BDn",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JMSUp  =  MCData("FatJet_TrackMass_Ratio_JMSUp",  DoubleFlav,   "FatJet_TrackMass_Ratio_JMSUp",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JMSDn  =  MCData("FatJet_TrackMass_Ratio_JMSDn",  DoubleFlav,   "FatJet_TrackMass_Ratio_JMSDn",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JESUp  =  MCData("FatJet_TrackMass_Ratio_JESUp",  DoubleFlav,   "FatJet_TrackMass_Ratio_JESUp",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JESDn  =  MCData("FatJet_TrackMass_Ratio_JESDn",  DoubleFlav,   "FatJet_TrackMass_Ratio_JESDn",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)

H_FatJet_TrackMass_Ratio_pT1_Eta1  =  MCData("FatJet_TrackMass_Ratio_pT1_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT1_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_pT2_Eta1  =  MCData("FatJet_TrackMass_Ratio_pT2_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT2_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_pT3_Eta1  =  MCData("FatJet_TrackMass_Ratio_pT3_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT3_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_pT1_Eta2  =  MCData("FatJet_TrackMass_Ratio_pT1_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT1_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_pT2_Eta2  =  MCData("FatJet_TrackMass_Ratio_pT2_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT2_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_pT3_Eta2  =  MCData("FatJet_TrackMass_Ratio_pT3_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT3_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)

H_FatJet_TrackMass_Ratio_BUp_pT1_Eta1  =  MCData("FatJet_TrackMass_Ratio_BUp_pT1_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_BUp_pT1_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_BUp_pT2_Eta1  =  MCData("FatJet_TrackMass_Ratio_BUp_pT2_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_BUp_pT2_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_BUp_pT3_Eta1  =  MCData("FatJet_TrackMass_Ratio_BUp_pT3_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_BUp_pT3_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_BUp_pT1_Eta2  =  MCData("FatJet_TrackMass_Ratio_BUp_pT1_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_BUp_pT1_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_BUp_pT2_Eta2  =  MCData("FatJet_TrackMass_Ratio_BUp_pT2_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_BUp_pT2_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_BUp_pT3_Eta2  =  MCData("FatJet_TrackMass_Ratio_BUp_pT3_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_BUp_pT3_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)

H_FatJet_TrackMass_Ratio_JESUp_pT1_Eta1  =  MCData("FatJet_TrackMass_Ratio_JESUp_pT1_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_JESUp_pT1_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JESUp_pT2_Eta1  =  MCData("FatJet_TrackMass_Ratio_JESUp_pT2_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_JESUp_pT2_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JESUp_pT3_Eta1  =  MCData("FatJet_TrackMass_Ratio_JESUp_pT3_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_JESUp_pT3_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JESUp_pT1_Eta2  =  MCData("FatJet_TrackMass_Ratio_JESUp_pT1_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_JESUp_pT1_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JESUp_pT2_Eta2  =  MCData("FatJet_TrackMass_Ratio_JESUp_pT2_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_JESUp_pT2_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JESUp_pT3_Eta2  =  MCData("FatJet_TrackMass_Ratio_JESUp_pT3_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_JESUp_pT3_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)

H_FatJet_TrackMass_Ratio_JMSUp_pT1_Eta1  =  MCData("FatJet_TrackMass_Ratio_JMSUp_pT1_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_JMSUp_pT1_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JMSUp_pT2_Eta1  =  MCData("FatJet_TrackMass_Ratio_JMSUp_pT2_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_JMSUp_pT2_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JMSUp_pT3_Eta1  =  MCData("FatJet_TrackMass_Ratio_JMSUp_pT3_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_JMSUp_pT3_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JMSUp_pT1_Eta2  =  MCData("FatJet_TrackMass_Ratio_JMSUp_pT1_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_JMSUp_pT1_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JMSUp_pT2_Eta2  =  MCData("FatJet_TrackMass_Ratio_JMSUp_pT2_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_JMSUp_pT2_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)
H_FatJet_TrackMass_Ratio_JMSUp_pT3_Eta2  =  MCData("FatJet_TrackMass_Ratio_JMSUp_pT3_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_JMSUp_pT3_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin",   True)

H_FatJet_TrackMass_Ratio_Herwig_pT1_Eta1  =  MCData("FatJet_TrackMass_Ratio_pT1_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT1_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin", True,  False, Herwigfile)
H_FatJet_TrackMass_Ratio_Herwig_pT2_Eta1  =  MCData("FatJet_TrackMass_Ratio_pT2_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT2_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin", True,   False, Herwigfile)
H_FatJet_TrackMass_Ratio_Herwig_pT3_Eta1  =  MCData("FatJet_TrackMass_Ratio_pT3_Eta1",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT3_Eta1",  "Calo Jet Mass/ Track Mass", "Jets / Bin", True,   False, Herwigfile)
H_FatJet_TrackMass_Ratio_Herwig_pT1_Eta2  =  MCData("FatJet_TrackMass_Ratio_pT1_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT1_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin", True,   False, Herwigfile)
H_FatJet_TrackMass_Ratio_Herwig_pT2_Eta2  =  MCData("FatJet_TrackMass_Ratio_pT2_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT2_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin", True,   False, Herwigfile)
H_FatJet_TrackMass_Ratio_Herwig_pT3_Eta2  =  MCData("FatJet_TrackMass_Ratio_pT3_Eta2",  DoubleFlav,   "FatJet_TrackMass_Ratio_pT3_Eta2",  "Calo Jet Mass/ Track Mass", "Jets / Bin", True,   False, Herwigfile)


#### Track Mass Ratio
H_FatJet_TrackpT_Ratio  =  MCData("FatJet_TrackpT_Ratio",  DoubleFlav,   "FatJet_TrackpT_Ratio",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio  =  MCData("FatJet_TrackD2_Ratio",  DoubleFlav,   "FatJet_TrackD2_Ratio",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_BUp  =  MCData("FatJet_TrackpT_Ratio_BUp",  DoubleFlav,   "FatJet_TrackpT_Ratio_BUp",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_BDn  =  MCData("FatJet_TrackpT_Ratio_BDn",  DoubleFlav,   "FatJet_TrackpT_Ratio_BDn",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_FJESUp  =  MCData("FatJet_TrackpT_Ratio_FJESUp",  DoubleFlav,   "FatJet_TrackpT_Ratio_FJESUp",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_FJESDn  =  MCData("FatJet_TrackpT_Ratio_FJESDn",  DoubleFlav,   "FatJet_TrackpT_Ratio_FJESDn",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)

H_FatJet_TrackpT_Ratio_pT1  =  MCData("FatJet_TrackpT_Ratio_pT1",  DoubleFlav,   "FatJet_TrackpT_Ratio_pT1",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_pT2  =  MCData("FatJet_TrackpT_Ratio_pT2",  DoubleFlav,   "FatJet_TrackpT_Ratio_pT2",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_pT3  =  MCData("FatJet_TrackpT_Ratio_pT3",  DoubleFlav,   "FatJet_TrackpT_Ratio_pT3",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)

H_FatJet_TrackD2_Ratio_pT1  =  MCData("FatJet_TrackD2_Ratio_pT1",  DoubleFlav,   "FatJet_TrackD2_Ratio_pT1",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_pT2  =  MCData("FatJet_TrackD2_Ratio_pT2",  DoubleFlav,   "FatJet_TrackD2_Ratio_pT2",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_pT3  =  MCData("FatJet_TrackD2_Ratio_pT3",  DoubleFlav,   "FatJet_TrackD2_Ratio_pT3",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)

H_FatJet_Herwig_TrackpT_Ratio_pT1  =  MCData("FatJet_TrackpT_Ratio_pT1",  DoubleFlav,   "FatJet_TrackpT_Ratio_pT1",  "Calo Jet pT/ Track pT", "Jets / Bin",   True, False, Herwigfile)
H_FatJet_Herwig_TrackpT_Ratio_pT2  =  MCData("FatJet_TrackpT_Ratio_pT2",  DoubleFlav,   "FatJet_TrackpT_Ratio_pT2",  "Calo Jet pT/ Track pT", "Jets / Bin",   True, False, Herwigfile)
H_FatJet_Herwig_TrackpT_Ratio_pT3  =  MCData("FatJet_TrackpT_Ratio_pT3",  DoubleFlav,   "FatJet_TrackpT_Ratio_pT3",  "Calo Jet pT/ Track pT", "Jets / Bin",   True, False, Herwigfile)

H_FatJet_Truth_TrackpT_Ratio_pT1  =  MCData("FatJet_TrackpT_Ratio_pT1",  DoubleFlav,   "FatJet_TrackpT_Ratio_pT1",  "Calo Jet pT/ Track pT", "Jets / Bin",   True, False, PythiaTruthfile)
H_FatJet_Truth_TrackpT_Ratio_pT2  =  MCData("FatJet_TrackpT_Ratio_pT2",  DoubleFlav,   "FatJet_TrackpT_Ratio_pT2",  "Calo Jet pT/ Track pT", "Jets / Bin",   True, False, PythiaTruthfile)
H_FatJet_Truth_TrackpT_Ratio_pT3  =  MCData("FatJet_TrackpT_Ratio_pT3",  DoubleFlav,   "FatJet_TrackpT_Ratio_pT3",  "Calo Jet pT/ Track pT", "Jets / Bin",   True, False, PythiaTruthfile)

H_FatJet_TrackpT_Ratio_BUp_pT1  =  MCData("FatJet_TrackpT_Ratio_BUp_pT1",  DoubleFlav,   "FatJet_TrackpT_Ratio_BUp_pT1",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_BUp_pT2  =  MCData("FatJet_TrackpT_Ratio_BUp_pT2",  DoubleFlav,   "FatJet_TrackpT_Ratio_BUp_pT2",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_BUp_pT3  =  MCData("FatJet_TrackpT_Ratio_BUp_pT3",  DoubleFlav,   "FatJet_TrackpT_Ratio_BUp_pT3",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_BDn_pT1  =  MCData("FatJet_TrackpT_Ratio_BDn_pT1",  DoubleFlav,   "FatJet_TrackpT_Ratio_BDn_pT1",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_BDn_pT2  =  MCData("FatJet_TrackpT_Ratio_BDn_pT2",  DoubleFlav,   "FatJet_TrackpT_Ratio_BDn_pT2",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_BDn_pT3  =  MCData("FatJet_TrackpT_Ratio_BDn_pT3",  DoubleFlav,   "FatJet_TrackpT_Ratio_BDn_pT3",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_FJESUp_pT1  =  MCData("FatJet_TrackpT_Ratio_FJESUp_pT1",  DoubleFlav,   "FatJet_TrackpT_Ratio_FJESUp_pT1",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_FJESUp_pT2  =  MCData("FatJet_TrackpT_Ratio_FJESUp_pT2",  DoubleFlav,   "FatJet_TrackpT_Ratio_FJESUp_pT2",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_FJESUp_pT3  =  MCData("FatJet_TrackpT_Ratio_FJESUp_pT3",  DoubleFlav,   "FatJet_TrackpT_Ratio_FJESUp_pT3",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_FJESDn_pT1  =  MCData("FatJet_TrackpT_Ratio_FJESDn_pT1",  DoubleFlav,   "FatJet_TrackpT_Ratio_FJESDn_pT1",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_FJESDn_pT2  =  MCData("FatJet_TrackpT_Ratio_FJESDn_pT2",  DoubleFlav,   "FatJet_TrackpT_Ratio_FJESDn_pT2",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)
H_FatJet_TrackpT_Ratio_FJESDn_pT3  =  MCData("FatJet_TrackpT_Ratio_FJESDn_pT3",  DoubleFlav,   "FatJet_TrackpT_Ratio_FJESDn_pT3",  "Calo Jet pT/ Track pT", "Jets / Bin",   True)

H_FatJet_TrackD2_Ratio_BUp_pT1  =  MCData("FatJet_TrackD2_Ratio_BUp_pT1",  DoubleFlav,   "FatJet_TrackD2_Ratio_BUp_pT1",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_BUp_pT2  =  MCData("FatJet_TrackD2_Ratio_BUp_pT2",  DoubleFlav,   "FatJet_TrackD2_Ratio_BUp_pT2",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_BUp_pT3  =  MCData("FatJet_TrackD2_Ratio_BUp_pT3",  DoubleFlav,   "FatJet_TrackD2_Ratio_BUp_pT3",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_D2Up_pT1  =  MCData("FatJet_TrackD2_Ratio_D2Up_pT1",  DoubleFlav,   "FatJet_TrackD2_Ratio_D2Up_pT1",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_D2Up_pT2  =  MCData("FatJet_TrackD2_Ratio_D2Up_pT2",  DoubleFlav,   "FatJet_TrackD2_Ratio_D2Up_pT2",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_D2Up_pT3  =  MCData("FatJet_TrackD2_Ratio_D2Up_pT3",  DoubleFlav,   "FatJet_TrackD2_Ratio_D2Up_pT3",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_BDn_pT1  =  MCData("FatJet_TrackD2_Ratio_BDn_pT1",  DoubleFlav,   "FatJet_TrackD2_Ratio_BDn_pT1",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_BDn_pT2  =  MCData("FatJet_TrackD2_Ratio_BDn_pT2",  DoubleFlav,   "FatJet_TrackD2_Ratio_BDn_pT2",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_BDn_pT3  =  MCData("FatJet_TrackD2_Ratio_BDn_pT3",  DoubleFlav,   "FatJet_TrackD2_Ratio_BDn_pT3",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_D2Dn_pT1  =  MCData("FatJet_TrackD2_Ratio_D2Dn_pT1",  DoubleFlav,   "FatJet_TrackD2_Ratio_D2Dn_pT1",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_D2Dn_pT2  =  MCData("FatJet_TrackD2_Ratio_D2Dn_pT2",  DoubleFlav,   "FatJet_TrackD2_Ratio_D2Dn_pT2",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_D2Dn_pT3  =  MCData("FatJet_TrackD2_Ratio_D2Dn_pT3",  DoubleFlav,   "FatJet_TrackD2_Ratio_D2Dn_pT3",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_FJESUp_pT1  =  MCData("FatJet_TrackD2_Ratio_FJESUp_pT1",  DoubleFlav,   "FatJet_TrackD2_Ratio_FJESUp_pT1",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_FJESUp_pT2  =  MCData("FatJet_TrackD2_Ratio_FJESUp_pT2",  DoubleFlav,   "FatJet_TrackD2_Ratio_FJESUp_pT2",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_FJESUp_pT3  =  MCData("FatJet_TrackD2_Ratio_FJESUp_pT3",  DoubleFlav,   "FatJet_TrackD2_Ratio_FJESUp_pT3",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_FJESDn_pT1  =  MCData("FatJet_TrackD2_Ratio_FJESDn_pT1",  DoubleFlav,   "FatJet_TrackD2_Ratio_FJESDn_pT1",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_FJESDn_pT2  =  MCData("FatJet_TrackD2_Ratio_FJESDn_pT2",  DoubleFlav,   "FatJet_TrackD2_Ratio_FJESDn_pT2",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)
H_FatJet_TrackD2_Ratio_FJESDn_pT3  =  MCData("FatJet_TrackD2_Ratio_FJESDn_pT3",  DoubleFlav,   "FatJet_TrackD2_Ratio_FJESDn_pT3",  "Calo Jet D2/ Track D2", "Jets / Bin",   True)


### nTrack Plots
H_FatJet_nTrack  =  MCData("FatJet_nTrack",  DoubleFlav,   "nTrack",  "nTrack", "",   True)
H_FatJet_nTrack_pT1  =  MCData("FatJet_nTrack_pT1",  DoubleFlav,   "nTrack_pT1",  "nTrack", "",   True)
H_FatJet_nTrack_pT2  =  MCData("FatJet_nTrack_pT2",  DoubleFlav,   "nTrack_pT2",  "nTrack", "",   True)
H_FatJet_nTrack_pT3  =  MCData("FatJet_nTrack_pT3",  DoubleFlav,   "nTrack_pT3",  "nTrack", "",   True)

### Mass Systemtics

H_FatJet_tau31 =  MCData("FatJet_tau31", DoubleFlav,   "FatJet_tau31", "#tau_{31}",   "Jets / Bin",   True)
H_FatJet_tau21 =  MCData("FatJet_tau21", DoubleFlav,   "FatJet_tau21", "#tau_{21}",   "Jets / Bin",   True)
H_FatJet_tau21WTA =  MCData("FatJet_tau21WTA", DoubleFlav,   "FatJet_tau21WTA", "#tau_{21}^{WTA}",   "Jets / 0.1",   True)
H_FatJet_tau32 =  MCData("FatJet_tau32", DoubleFlav,   "FatJet_tau32", "#tau_{32}",   "Jets / Bin",   True)
H_FatJet_d2 =  MCData("FatJet_d2", DoubleFlav,   "FatJet_d2", "D_{2}^{(#beta=1)}",   "Jets / 0.3",   True)
#H_FatJet_c2 =  MCData("FatJet_c2", DoubleFlav,   "FatJet_c2", "C2",   "Jets / Bin",   True)
H_FatJet_width =  MCData("FatJet_width", DoubleFlav,   "FatJet_width", "width",   "Jets / Bin",   True)

H_MuJet_pT_JESUp =  MCData("MuJet_pT_JESUp", SingleFlav,   "MuJet_pT_JESUp", "Muon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_MuJet_pT_JESDn =  MCData("MuJet_pT_JESDn", SingleFlav,   "MuJet_pT_JESDn", "Muon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_MuJet_pT_BUp =  MCData("MuJet_pT_BUp", SingleFlav,   "MuJet_pT_BUp", "Muon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_MuJet_pT_BDn =  MCData("MuJet_pT_BDn", SingleFlav,   "MuJet_pT_BDn", "Muon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_MuJet_pT_CUp =  MCData("MuJet_pT_CUp", SingleFlav,   "MuJet_pT_CUp", "Muon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_MuJet_pT_CDn =  MCData("MuJet_pT_CDn", SingleFlav,   "MuJet_pT_CDn", "Muon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_MuJet_pT_LUp =  MCData("MuJet_pT_LUp", SingleFlav,   "MuJet_pT_LUp", "Muon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_MuJet_pT_LDn =  MCData("MuJet_pT_LDn", SingleFlav,   "MuJet_pT_LDn", "Muon Jet p_{T} [GeV]", "Jets / Bin",   True)

H_MuJet_Eta_JESUp =  MCData("MuJet_Eta_JESUp", SingleFlav, "MuJet_Eta_JESUp", "Muon Jet #eta", "Jets / Bin",True)
H_MuJet_Eta_JESDn =  MCData("MuJet_Eta_JESDn", SingleFlav, "MuJet_Eta_JESDn", "Muon Jet #eta", "Jets / Bin",True)
H_MuJet_Eta_BUp =  MCData("MuJet_Eta_BUp", SingleFlav, "MuJet_Eta_BUp", "Muon Jet #eta", "Jets / Bin",True)
H_MuJet_Eta_BDn =  MCData("MuJet_Eta_BDn", SingleFlav, "MuJet_Eta_BDn", "Muon Jet #eta", "Jets / Bin",True)
H_MuJet_Eta_CUp =  MCData("MuJet_Eta_CUp", SingleFlav, "MuJet_Eta_CUp", "Muon Jet #eta", "Jets / Bin",True)
H_MuJet_Eta_CDn =  MCData("MuJet_Eta_CDn", SingleFlav, "MuJet_Eta_CDn", "Muon Jet #eta", "Jets / Bin",True)
H_MuJet_Eta_LUp =  MCData("MuJet_Eta_LUp", SingleFlav, "MuJet_Eta_LUp", "Muon Jet #eta", "Jets / Bin",True)
H_MuJet_Eta_LDn =  MCData("MuJet_Eta_LDn", SingleFlav, "MuJet_Eta_LDn", "Muon Jet #eta", "Jets / Bin",True)

H_Muon_pT_BUp  =  MCData("Mu_pT_BUp",    SingleFlav,   "Mu_pT_BUp",    "Muon p_{T} [GeV]", "Jets / Bin",   True)
H_Muon_pT_BDn  =  MCData("Mu_pT_BDn",    SingleFlav,   "Mu_pT_BDn",    "Muon p_{T} [GeV]", "Jets / Bin",   True)
H_Muon_pT_CUp  =  MCData("Mu_pT_CUp",    SingleFlav,   "Mu_pT_CUp",    "Muon p_{T} [GeV]", "Jets / Bin",   True)
H_Muon_pT_CDn  =  MCData("Mu_pT_CDn",    SingleFlav,   "Mu_pT_CDn",    "Muon p_{T} [GeV]", "Jets / Bin",   True)
H_Muon_pT_LUp  =  MCData("Mu_pT_LUp",    SingleFlav,   "Mu_pT_LUp",    "Muon p_{T} [GeV]", "Jets / Bin",   True)
H_Muon_pT_LDn  =  MCData("Mu_pT_LDn",    SingleFlav,   "Mu_pT_LDn",    "Muon p_{T} [GeV]", "Jets / Bin",   True)
H_Muon_pT_JESUp = MCData("Mu_pT_JESUp",    SingleFlav,   "Mu_pT_JESUp",    "Muon p_{T} [GeV]", "Jets / Bin",   True)
H_Muon_pT_JESDn = MCData("Mu_pT_JESDn",    SingleFlav,   "Mu_pT_JESDn",    "Muon p_{T} [GeV]", "Jets / Bin",   True)

H_dR_Mu_MuJet_BUp = MCData("dR_Mu_MuJet_BUp",  SingleFlav, "dR_Mu_MuJet", "#DeltaR(Muon, Muon Jet)", "Jets / Bin", True)
H_dR_Mu_MuJet_BDn = MCData("dR_Mu_MuJet_BDn",  SingleFlav, "dR_Mu_MuJet", "#DeltaR(Muon, Muon Jet)", "Jets / Bin", True)
H_dR_Mu_MuJet_CUp = MCData("dR_Mu_MuJet_CUp",  SingleFlav, "dR_Mu_MuJet", "#DeltaR(Muon, Muon Jet)", "Jets / Bin", True)
H_dR_Mu_MuJet_CDn = MCData("dR_Mu_MuJet_CDn",  SingleFlav, "dR_Mu_MuJet", "#DeltaR(Muon, Muon Jet)", "Jets / Bin", True)
H_dR_Mu_MuJet_LUp = MCData("dR_Mu_MuJet_LUp",  SingleFlav, "dR_Mu_MuJet", "#DeltaR(Muon, Muon Jet)", "Jets / Bin", True)
H_dR_Mu_MuJet_LDn = MCData("dR_Mu_MuJet_LDn",  SingleFlav, "dR_Mu_MuJet", "#DeltaR(Muon, Muon Jet)", "Jets / Bin", True)

H_NonMuJet_pT_BUp =  MCData("NonMuJet_pT_BUp", SingleFlav,   "NonMuJet_pT_BUp", "NonMuon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_NonMuJet_pT_BDn =  MCData("NonMuJet_pT_BDn", SingleFlav,   "NonMuJet_pT_BDn", "NonMuon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_NonMuJet_pT_CUp =  MCData("NonMuJet_pT_CUp", SingleFlav,   "NonMuJet_pT_CUp", "NonMuon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_NonMuJet_pT_CDn =  MCData("NonMuJet_pT_CDn", SingleFlav,   "NonMuJet_pT_CDn", "NonMuon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_NonMuJet_pT_LUp =  MCData("NonMuJet_pT_LUp", SingleFlav,   "NonMuJet_pT_LUp", "NonMuon Jet p_{T} [GeV]", "Jets / Bin",   True)
H_NonMuJet_pT_LDn =  MCData("NonMuJet_pT_LDn", SingleFlav,   "NonMuJet_pT_LDn", "NonMuon Jet p_{T} [GeV]", "Jets / Bin",   True)

H_NonMuJet_Eta_BUp =  MCData("NonMuJet_Eta_BUp", SingleFlav, "NonMuJet_Eta_BUp", "NonMuon Jet #eta", "Jets / Bin",True)
H_NonMuJet_Eta_BDn =  MCData("NonMuJet_Eta_BDn", SingleFlav, "NonMuJet_Eta_BDn", "NonMuon Jet #eta", "Jets / Bin",True)
H_NonMuJet_Eta_CUp =  MCData("NonMuJet_Eta_CUp", SingleFlav, "NonMuJet_Eta_CUp", "NonMuon Jet #eta", "Jets / Bin",True)
H_NonMuJet_Eta_CDn =  MCData("NonMuJet_Eta_CDn", SingleFlav, "NonMuJet_Eta_CDn", "NonMuon Jet #eta", "Jets / Bin",True)
H_NonMuJet_Eta_LUp =  MCData("NonMuJet_Eta_LUp", SingleFlav, "NonMuJet_Eta_LUp", "NonMuon Jet #eta", "Jets / Bin",True)
H_NonMuJet_Eta_LDn =  MCData("NonMuJet_Eta_LDn", SingleFlav, "NonMuJet_Eta_LDn", "NonMuon Jet #eta", "Jets / Bin",True)


H_FatJet_NJet_BUp    =  MCData("FatJet_NJet_BUp",    DoubleFlav,   "FatJet_NJet_BUp",    "NJet Matched", "Number of Events",   True)
H_FatJet_NJet_BDn    =  MCData("FatJet_NJet_BDn",    DoubleFlav,   "FatJet_NJet_BDn",    "NJet Matched", "Number of Events",   True)
H_FatJet_NJet_JESUp  =  MCData("FatJet_NJet_JESUp",    DoubleFlav,   "FatJet_NJet_JESUp",    "NJet Matched", "Number of Events",   True)
H_FatJet_NJet_JESDn  =  MCData("FatJet_NJet_JESDn",    DoubleFlav,   "FatJet_NJet_JESDn",    "NJet Matched", "Number of Events",   True)

H_FatJet_MOPT_BUp    =  MCData("FatJet_MOPT_BUp", DoubleFlav,   "FatJet_MOPT_BUp", "Fat Jet Mass / p_{T}",   "Jets / Bin",   True)
H_FatJet_MOPT_BDn    =  MCData("FatJet_MOPT_BDn", DoubleFlav,   "FatJet_MOPT_BUp", "Fat Jet Mass / p_{T}",   "Jets / Bin",   True)
H_FatJet_MOPT_JESUp  =  MCData("FatJet_MOPT_JESUp", DoubleFlav,   "FatJet_MOPT_JESUp", "Fat Jet Mass / p_{T}",   "Jets / Bin",   True)
H_FatJet_MOPT_JESDn  =  MCData("FatJet_MOPT_JESDn", DoubleFlav,   "FatJet_MOPT_JESDn", "Fat Jet Mass / p_{T}",   "Jets / Bin",   True)
H_FatJet_MOPT_JMSUp  =  MCData("FatJet_MOPT_JMSUp", DoubleFlav,   "FatJet_MOPT_JMSUp", "Fat Jet Mass / p_{T}",   "Jets / Bin",   True)
H_FatJet_MOPT_JMSDn  =  MCData("FatJet_MOPT_JMSDn", DoubleFlav,   "FatJet_MOPT_JMSDn", "Fat Jet Mass / p_{T}",   "Jets / Bin",   True)

H_FatJet_tau31_BUp   =  MCData("FatJet_tau31_BUp", DoubleFlav,   "FatJet_tau31_BUp", "#tau_{31}",   "Jets / Bin",   True)
H_FatJet_tau31_BDn   =  MCData("FatJet_tau31_BDn", DoubleFlav,   "FatJet_tau31_BDn", "#tau_{31}",   "Jets / Bin",   True)
H_FatJet_tau31_JESUp =  MCData("FatJet_tau31_JESUp", DoubleFlav,   "FatJet_tau31_JESUp", "#tau_{31}",   "Jets / Bin",   True)
H_FatJet_tau31_JESDn =  MCData("FatJet_tau31_JESDn", DoubleFlav,   "FatJet_tau31_JESDn", "#tau_{31}",   "Jets / Bin",   True)

H_FatJet_tau32_BUp   =  MCData("FatJet_tau32_BUp", DoubleFlav,   "FatJet_tau32_BUp", "#tau_{32}",   "Jets / Bin",   True)
H_FatJet_tau32_BDn   =  MCData("FatJet_tau32_BDn", DoubleFlav,   "FatJet_tau32_BDn", "#tau_{32}",   "Jets / Bin",   True)
H_FatJet_tau32_JESUp =  MCData("FatJet_tau32_JESUp", DoubleFlav,   "FatJet_tau32_JESUp", "#tau_{32}",   "Jets / Bin",   True)
H_FatJet_tau32_JESDn =  MCData("FatJet_tau32_JESDn", DoubleFlav,   "FatJet_tau32_JESDn", "#tau_{32}",   "Jets / Bin",   True)

H_FatJet_tau21_BUp   =  MCData("FatJet_tau21_BUp", DoubleFlav,   "FatJet_tau21_BUp", "#tau_{21}",   "Jets / Bin",   True)
H_FatJet_tau21_BDn   =  MCData("FatJet_tau21_BDn", DoubleFlav,   "FatJet_tau21_BDn", "#tau_{21}",   "Jets / Bin",   True)
H_FatJet_tau21_JESUp =  MCData("FatJet_tau21_JESUp", DoubleFlav,   "FatJet_tau21_JESUp", "#tau_{21}",   "Jets / Bin",   True)
H_FatJet_tau21_JESDn =  MCData("FatJet_tau21_JESDn", DoubleFlav,   "FatJet_tau21_JESDn", "#tau_{21}",   "Jets / Bin",   True)

H_FatJet_tau21WTA_BUp   =  MCData("FatJet_tau21WTA_BUp", DoubleFlav,   "FatJet_tau21WTA_BUp", "#tau_{21}",   "Jets / 0.1",   True)
H_FatJet_tau21WTA_BDn   =  MCData("FatJet_tau21WTA_BDn", DoubleFlav,   "FatJet_tau21WTA_BDn", "#tau_{21}",   "Jets / 0.1",   True)
H_FatJet_tau21WTA_JESUp =  MCData("FatJet_tau21WTA_JESUp", DoubleFlav,   "FatJet_tau21WTA_JESUp", "#tau_{21}",   "Jets / 0.1",   True)
H_FatJet_tau21WTA_JESDn =  MCData("FatJet_tau21WTA_JESDn", DoubleFlav,   "FatJet_tau21WTA_JESDn", "#tau_{21}",   "Jets / 0.1",   True)
H_FatJet_tau21WTA_tauUp =  MCData("FatJet_tau21WTA_tauUp", DoubleFlav,   "FatJet_tau21WTA_tauUp", "#tau_{21}",   "Jets / 0.1",   True)
H_FatJet_tau21WTA_tauDn =  MCData("FatJet_tau21WTA_tauDn", DoubleFlav,   "FatJet_tau21WTA_tauDn", "#tau_{21}",   "Jets / 0.1",   True)

H_FatJet_d2_BUp      =  MCData("FatJet_d2_BUp", DoubleFlav,   "FatJet_d2_BUp", "D_{2}^{(#beta=1)}",   "Jets / 0.3",   True)
H_FatJet_d2_BDn      =  MCData("FatJet_d2_BDn", DoubleFlav,   "FatJet_d2_BDn", "D_{2}^{(#beta=1)}",   "Jets / 0.3",   True)
H_FatJet_d2_JESUp    =  MCData("FatJet_d2_JESUp", DoubleFlav,   "FatJet_d2_JESUp", "D_{2}^{(#beta=1)}",   "Jets / 0.3",   True)
H_FatJet_d2_JESDn    =  MCData("FatJet_d2_JESDn", DoubleFlav,   "FatJet_d2_JESDn", "D_{2}^{(#beta=1)}",   "Jets / 0.3",   True)
H_FatJet_d2_d2Up    =  MCData("FatJet_d2_d2Up", DoubleFlav,   "FatJet_d2_d2Up", "D_{2}^{(#beta=1)}",   "Jets / 0.3",   True)
H_FatJet_d2_d2Dn    =  MCData("FatJet_d2_d2Dn", DoubleFlav,   "FatJet_d2_d2Dn", "D_{2}^{(#beta=1)}",   "Jets / 0.3",   True)


H_dR_FitSys = []
for ifile in range(12):
    H_dR_FitSys.append( MCData("SubJet_dR",  DoubleFlav,   "SubJet_dR", "#DeltaR(Jet_{1}, Jet_{2})", "Entries / 0.1", True, False, pTSysFiles[ifile]) )

H_dR_Mu_MuJet_FitSys = []
for ifile in range(12):
    H_dR_Mu_MuJet_FitSys.append( MCData("dR_Mu_MuJet",  DoubleFlav,   "dR_Mu_MuJet", "#DeltaR(Muon, Muon Jet)", "Entries", True, False, pTSysFiles[ifile]) )

H_Muon_pT_FitSys = []
for ifile in range(12):
    H_Muon_pT_FitSys.append( MCData("Mu_pT_FitSys",  DoubleFlav,   "Mu_pT_FitSys", "p_{T} [GeV]", "Jets / Bin",  True, False, pTSysFiles[ifile]) )

H_MuJet_pT_FitSys = []
for ifile in range(12):
    H_MuJet_pT_FitSys.append( MCData("MuJet_pT_FitSys",  DoubleFlav,   "MuJet_pT_FitSys", "p_{T} [GeV]", "Jets / Bin",  True, False, pTSysFiles[ifile]) )

H_NonMuJet_pT_FitSys = []
for ifile in range(12):
    H_NonMuJet_pT_FitSys.append( MCData("NonMuJet_pT_FitSys",  DoubleFlav,   "NonMuJet_pT_FitSys", "p_{T} [GeV]", "Jets / Bin",  True, False, pTSysFiles[ifile]) )

H_MuJet_Eta_FitSys = []
for ifile in range(12):
    H_MuJet_Eta_FitSys.append( MCData("MuJet_Eta_FitSys",  DoubleFlav,   "MuJet_Eta_FitSys", "#eta", "Jets / Bin",  True, False, pTSysFiles[ifile]) )

H_NonMuJet_Eta_FitSys = []
for ifile in range(12):
    H_NonMuJet_Eta_FitSys.append( MCData("NonMuJet_Eta_FitSys",  DoubleFlav,   "NonMuJet_Eta_FitSys", "#eta", "Jets / Bin",  True, False, pTSysFiles[ifile]) )

H_FatJet_pT_FitSys = []
for ifile in range(12):
    H_FatJet_pT_FitSys.append( MCData("FatJet_pT",  DoubleFlav,  "FatJet_pT", "pT [GeV]",   "Jets / 50 GeV", True, False, pTSysFiles[ifile]) )

H_FatJet_Eta_FitSys = []
for ifile in range(12):
    H_FatJet_Eta_FitSys.append( MCData("FatJet_Eta",  DoubleFlav,  "FatJet_Eta", "#eta",   "Jets / 0.5", True, False, pTSysFiles[ifile]) )

H_FatJet_Mass_FitSys = []
for ifile in range(12):
    H_FatJet_Mass_FitSys.append( MCData("FatJet_Mass",  DoubleFlav,  "FatJet_Mass", "Mass [GeV]",   "Jets / 50 GeV", True, False, pTSysFiles[ifile]) )

H_FatJet_d2_FitSys = []
for ifile in range(12):
    H_FatJet_d2_FitSys.append( MCData("FatJet_d2",  DoubleFlav,  "FatJet_d2", "D_{2}^{(#beta=1)}",   "Jets / 0.3" , True, False, pTSysFiles[ifile]) )

H_FatJet_tau21WTA_FitSys = []
for ifile in range(12):
    H_FatJet_tau21WTA_FitSys.append( MCData("FatJet_tau21WTA",  DoubleFlav,  "FatJet_tau21WTA", "#tau_{21}",   "Jets / 0.1", True, False, pTSysFiles[ifile]) )


######### b tagging rate

H_FatJet_pT_FitUp, H_FatJet_pT_FitDn, H_FatJet_pT_BFitUp, H_FatJet_pT_BFitDn   = ComputeFitSys(H_FatJet_pT, H_FatJet_pT_FitSys, CorMatrices, FitBinUncert)

DrawTool.SetCompData(False)

FatJet_pT_BtagRate_Data = CopyHist(H_FatJet_pT.BData)
FatJet_pT_BtagRate_Data.Divide(H_FatJet_pT.Data)

FatJet_pT_BtagRate_MC = CopyHist(H_FatJet_pT.BMC)
FatJet_pT_BtagRate_MC.Divide(H_FatJet_pT.MC)

FatJet_pT_BtagRate_MC_FitUp = CopyHist(H_FatJet_pT_BFitUp)
FatJet_pT_BtagRate_MC_FitUp.Divide(H_FatJet_pT_FitUp)

FatJet_pT_BtagRate_MC_FitDn = CopyHist(H_FatJet_pT_BFitDn)
FatJet_pT_BtagRate_MC_FitDn.Divide(H_FatJet_pT_FitDn)

FatJet_pT_BtagRate_MC_BFitUp = CopyHist(H_FatJet_pT_BFitUp)
FatJet_pT_BtagRate_MC_BFitUp.Divide(H_FatJet_pT.MC)

FatJet_pT_BtagRate_MC_BFitDn = CopyHist(H_FatJet_pT_BFitDn)
FatJet_pT_BtagRate_MC_BFitDn.Divide(H_FatJet_pT.MC)

FatJet_pT_BtagRate_MC_JESUp = CopyHist(H_FatJet_JESUp.BMC)
FatJet_pT_BtagRate_MC_JESUp.Divide(H_FatJet_pT.MC)

FatJet_pT_BtagRate_MC_JESDn = CopyHist(H_FatJet_JESDn.BMC)
FatJet_pT_BtagRate_MC_JESDn.Divide(H_FatJet_pT.MC)

FatJet_pT_BtagRate_MC_BUp = CopyHist(H_FatJet_pT_BUp.BMC)
FatJet_pT_BtagRate_MC_BUp.Divide(H_FatJet_pT.MC)

FatJet_pT_BtagRate_MC_BDn = CopyHist(H_FatJet_pT_BDn.BMC)
FatJet_pT_BtagRate_MC_BDn.Divide(H_FatJet_pT.MC)


FatJet_pT_BtagRate_MC_Uncert = ComputeSysRatio(FatJet_pT_BtagRate_MC,
                                               [FatJet_pT_BtagRate_MC_BUp, FatJet_pT_BtagRate_MC_JESDn],
                                               [FatJet_pT_BtagRate_MC_BDn, FatJet_pT_BtagRate_MC_JESUp], ["BRate", "FJJES"])

FatJet_pT_BtagRate_MC_Uncert_Tot = ComputeSysRatio(FatJet_pT_BtagRate_MC,
                                                   [FatJet_pT_BtagRate_MC_BUp, FatJet_pT_BtagRate_MC_JESUp, FatJet_pT_BtagRate_MC_FitUp],
                                                   [FatJet_pT_BtagRate_MC_BDn, FatJet_pT_BtagRate_MC_JESDn, FatJet_pT_BtagRate_MC_FitDn], ["BRate", "FJJES", "FJFit"])

FatJet_pT_BtagRate_TotSys         = ComputeSys(FatJet_pT_BtagRate_MC,          
                                            [FatJet_pT_BtagRate_MC_BUp, FatJet_pT_BtagRate_MC_JESUp, FatJet_pT_BtagRate_MC_FitUp],
                                            [FatJet_pT_BtagRate_MC_BDn, FatJet_pT_BtagRate_MC_JESDn, FatJet_pT_BtagRate_MC_FitDn], ["BRate", "FJJES", "FJFit"])

FatJet_pT_BtagRate_BSys         = ComputeSys(FatJet_pT_BtagRate_MC,          
                                             [FatJet_pT_BtagRate_MC_BUp],
                                             [FatJet_pT_BtagRate_MC_BDn], ["BRate"])

FatJet_pT_BtagRate_Fit         = ComputeSys(FatJet_pT_BtagRate_MC,          
                                            [FatJet_pT_BtagRate_MC_FitUp],
                                            [FatJet_pT_BtagRate_MC_FitDn], ["FJFit"])

FatJet_pT_BtagRate_MC.SetMarkerSize(1.3)
FatJet_pT_BtagRate_MC.SetLineStyle(2)
FatJet_pT_BtagRate_Data.SetMarkerSize(1.3)
FatJet_pT_BtagRate_MC.SetMarkerStyle(26)
DrawTool.SetCompData(True)

FatJet_pT_BtagRate_TotSys.GetXaxis().SetTitle("p_{T} [GeV]")
FatJet_pT_BtagRate_Fit.GetXaxis().SetTitle("p_{T} [GeV]")
FatJet_pT_BtagRate_BSys.GetXaxis().SetTitle("p_{T} [GeV]")
FatJet_pT_BtagRate_MC.GetXaxis().SetTitle("p_{T} [GeV]")
FatJet_pT_BtagRate_Data.GetXaxis().SetTitle("p_{T} [GeV]")

#DrawTool.shift = 0.13
DrawTool.colorlist = [2, 2, 1]
DrawTool.studytype = "Preliminary"
DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.DrawHists("BtagRate_FatJetpT", ["p_{T} [GeV]", "Double b-tagging Rate"], [FatJet_pT_BtagRate_MC, FatJet_pT_BtagRate_MC, FatJet_pT_BtagRate_Data], ["Pythia 8 MC", "Pythia 8 MC", "Data"],[],[], [FatJet_pT_BtagRate_TotSys, FatJet_pT_BtagRate_Fit], FatJet_pT_BtagRate_BSys)
FatJet_pT_BtagRate_TotSys.Write()
FatJet_pT_BtagRate_Fit.Write()
FatJet_pT_BtagRate_BSys.Write()

#DrawTool.DrawHists("BtagRate_FatJetpT_Uncert", ["p_{T} [GeV]", "Double b-tagging Rate"], [FatJet_pT_BtagRate_MC_Uncert, FatJet_pT_BtagRate_MC_Uncert, FatJet_pT_BtagRate_Data], ["Pythia 8 MC", "Pythia 8 MC", "Data"])
#DrawTool.DrawHists("BtagRate_FatJetpT_Uncert_Tot", ["p_{T} [GeV]", "Double b-tagging Rate"], [FatJet_pT_BtagRate_MC_Uncert_Tot, FatJet_pT_BtagRate_MC_Uncert_Tot, FatJet_pT_BtagRate_Data], ["Pythia 8 MC", "Pythia 8 MC", "Data"])
#DrawTool.DrawHists("BtagRate_FatJetpT_BFitUp", ["p_{T} [GeV]", "Double b-tagging Rate"], [FatJet_pT_BtagRate_MC_BFitUp, FatJet_pT_BtagRate_MC_BFitUp, FatJet_pT_BtagRate_Data], ["Pythia 8 MC", "Pythia 8 MC", "Data"])
DrawTool.ClearTexts()
DrawTool.colorlist = DrawTool.origcolorlist

######### b frac plot

DrawTool.SetCompData(False)
##### BB
FatJet_pT_Fit_BB_Up, FatJet_pT_Fit_BB_Dn = ComputeFitSysSingleFlav(H_FatJet_pT_NoStack.Stack[0], H_FatJet_pT_FitSys , "BB", FitBinUncert)

FatJet_pT_BFrac_Fit_BB = CopyHist(H_FatJet_pT_NoStack.Stack[0])
FatJet_pT_BFrac_Fit_BB.Divide(H_FatJet_pT.MC)

FatJet_pT_BFrac_Fit_BB_Up = CopyHist(FatJet_pT_Fit_BB_Up)
FatJet_pT_BFrac_Fit_BB_Up.Divide(H_FatJet_pT.MC)

FatJet_pT_BFrac_Fit_BB_Dn = CopyHist(FatJet_pT_Fit_BB_Dn)
FatJet_pT_BFrac_Fit_BB_Dn.Divide(H_FatJet_pT.MC)

FatJet_pT_BFrac_Fit_BB_Uncert = ComputeSysRatio(FatJet_pT_BFrac_Fit_BB,
                                                [FatJet_pT_BFrac_Fit_BB_Up],[FatJet_pT_BFrac_Fit_BB_Dn], ['Fit'])


FatJet_pT_BFrac_MC_BB = CopyHist(H_FatJet_pT_NoReWeight_NoStack.Stack[0])
FatJet_pT_BFrac_MC_BB.Divide(H_FatJet_pT_NoReWeight.MC)

#### BL 

FatJet_pT_Fit_BL_Up, FatJet_pT_Fit_BL_Dn = ComputeFitSysSingleFlav(H_FatJet_pT_NoStack.Stack[2], H_FatJet_pT_FitSys , "BL", FitBinUncert)

FatJet_pT_BFrac_Fit_BL = CopyHist(H_FatJet_pT_NoStack.Stack[2])
FatJet_pT_BFrac_Fit_BL.Divide(H_FatJet_pT.MC)

FatJet_pT_BFrac_Fit_BL_Up = CopyHist(FatJet_pT_Fit_BL_Up)
FatJet_pT_BFrac_Fit_BL_Up.Divide(H_FatJet_pT.MC)

FatJet_pT_BFrac_Fit_BL_Dn = CopyHist(FatJet_pT_Fit_BL_Dn)
FatJet_pT_BFrac_Fit_BL_Dn.Divide(H_FatJet_pT.MC)

FatJet_pT_BFrac_Fit_BL_Uncert = ComputeSysRatio(FatJet_pT_BFrac_Fit_BL,
                                                [FatJet_pT_BFrac_Fit_BL_Up],[FatJet_pT_BFrac_Fit_BL_Dn], ['Fit'])


FatJet_pT_BFrac_MC_BL = CopyHist(H_FatJet_pT_NoReWeight_NoStack.Stack[2])
FatJet_pT_BFrac_MC_BL.Divide(H_FatJet_pT_NoReWeight.MC)


#### LL 

FatJet_pT_Fit_LL_Up, FatJet_pT_Fit_LL_Dn = ComputeFitSysSingleFlav(H_FatJet_pT_NoStack.Stack[5], H_FatJet_pT_FitSys , "LL", FitBinUncert)

FatJet_pT_BFrac_Fit_LL = CopyHist(H_FatJet_pT_NoStack.Stack[5])
FatJet_pT_BFrac_Fit_LL.Divide(H_FatJet_pT.MC)

FatJet_pT_BFrac_Fit_LL_Up = CopyHist(FatJet_pT_Fit_LL_Up)
FatJet_pT_BFrac_Fit_LL_Up.Divide(H_FatJet_pT.MC)

FatJet_pT_BFrac_Fit_LL_Dn = CopyHist(FatJet_pT_Fit_LL_Dn)
FatJet_pT_BFrac_Fit_LL_Dn.Divide(H_FatJet_pT.MC)

FatJet_pT_BFrac_Fit_LL_Uncert = ComputeSysRatio(FatJet_pT_BFrac_Fit_LL,
                                                [FatJet_pT_BFrac_Fit_LL_Up],[FatJet_pT_BFrac_Fit_LL_Dn], ['Fit'])


FatJet_pT_BFrac_MC_LL = CopyHist(H_FatJet_pT_NoReWeight_NoStack.Stack[5])
FatJet_pT_BFrac_MC_LL.Divide(H_FatJet_pT_NoReWeight.MC)


FatJet_pT_BFrac_Fit_BB_Uncert.SetLineStyle(2)
FatJet_pT_BFrac_Fit_BL_Uncert.SetLineStyle(2)
FatJet_pT_BFrac_Fit_LL_Uncert.SetLineStyle(2)

FatJet_pT_BFrac_Fit_BB_Uncert.SetMarkerStyle(4)
FatJet_pT_BFrac_Fit_BL_Uncert.SetMarkerStyle(4)
FatJet_pT_BFrac_Fit_LL_Uncert.SetMarkerStyle(4)

FatJet_pT_BFrac_Fit_BB_Uncert.SetMarkerSize(0.9)
FatJet_pT_BFrac_Fit_BL_Uncert.SetMarkerSize(0.9)
FatJet_pT_BFrac_Fit_LL_Uncert.SetMarkerSize(0.9)
FatJet_pT_BFrac_MC_BB.SetMarkerSize(0.9)
FatJet_pT_BFrac_MC_BL.SetMarkerSize(0.9)
FatJet_pT_BFrac_MC_LL.SetMarkerSize(0.9)

#DrawTool.shift = 0.13
DrawTool.colorlist = [1, 13, 2, 46, 4, 38]
DrawTool.colorlist = [13, 1, 46, 2, 38, 38, 4]

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.DrawHists("FitFrac_FatJetpT", ["p_{T} [GeV]", "Flavor Fraction"], [FatJet_pT_BFrac_Fit_LL, FatJet_pT_BFrac_MC_LL, FatJet_pT_BFrac_Fit_BL, FatJet_pT_BFrac_MC_BL, FatJet_pT_BFrac_Fit_BB, FatJet_pT_BFrac_MC_BB, ], ["LL fraction from fit", "LL fraction from MC prediction","BL fraction from fit", "BL fraction from MC prediction", "BB fraction from fit", "BB fraction from MC prediction"])

DrawTool.DrawHists("FitFracUncert_FatJetpT", ["p_{T} [GeV]", "Flavor Fraction"], [FatJet_pT_BFrac_Fit_LL_Uncert, FatJet_pT_BFrac_MC_LL, FatJet_pT_BFrac_Fit_BL_Uncert, FatJet_pT_BFrac_MC_BL, FatJet_pT_BFrac_Fit_BB_Uncert, FatJet_pT_BFrac_Fit_BB_Uncert, FatJet_pT_BFrac_MC_BB ], ["Fit LL fraction", "MC LL fraction","Fit BL fraction", "MC BL fraction", "Fit BB fraction", "Fit BB fraction", "MC BB fraction" ])

DrawTool.DrawHists("FitFracUncert_FatJetpT_MC", ["p_{T} [GeV]", "Flavor Fraction"], [FatJet_pT_BFrac_MC_LL, FatJet_pT_BFrac_MC_BL, FatJet_pT_BFrac_MC_BL, FatJet_pT_BFrac_MC_BB ], ["MC LL fraction","MC BL fraction","MC BL fraction", "MC BB fraction" ])

DrawTool.ClearTexts()

#DrawTool.shift = 0.1
DrawTool.SetCompData(True)

DrawTool.colorlist = DrawTool.origcolorlist



### Draw Basic Plots
#H_TrigJet_pT.DrawMCTot()
#H_TrigJet_Eta.DrawMCTot()
#H_TrigJet_Phi.DrawMCTot()
#H_NJet.DrawMCTot()

DrawTool.SetCompData(False)
H_dR_TrigJet_MuonJet.DrawMCOnly()


#PlotFits([H_FatJet_pTRef_Response_pT1.Data, H_FatJet_pTRef_Response_pT1.MC],"pTRef_Fit1", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data", "MC"])
#PlotFits([H_FatJet_pTRef_Response_pT2.Data, H_FatJet_pTRef_Response_pT2.MC],"pTRef_Fit2", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data", "MC"])
#PlotFits([H_FatJet_pTRef_Response_pT3.Data, H_FatJet_pTRef_Response_pT3.MC],"pTRef_Fit3", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data", "MC"])
#
#PlotFits([H_FatJet_pTRef_Response_pT1_Trig.Data, H_FatJet_pTRef_Response_pT1_Trig.MC],"pTRef_Fit1_Trig", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data", "MC"])
#PlotFits([H_FatJet_pTRef_Response_pT2_Trig.Data, H_FatJet_pTRef_Response_pT2_Trig.MC],"pTRef_Fit2_Trig", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data", "MC"])
#PlotFits([H_FatJet_pTRef_Response_pT3_Trig.Data, H_FatJet_pTRef_Response_pT3_Trig.MC],"pTRef_Fit3_Trig", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data", "MC"])
#
#PlotFits([H_FatJet_pTRef_Response_pT1.BData, H_FatJet_pTRef_Response_pT1.BMC],"pTRef_B_Fit1", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data double B", "MC double B"])
#PlotFits([H_FatJet_pTRef_Response_pT2.BData, H_FatJet_pTRef_Response_pT2.BMC],"pTRef_B_Fit2", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data double B", "MC double B"])
#PlotFits([H_FatJet_pTRef_Response_pT3.BData, H_FatJet_pTRef_Response_pT3.BMC],"pTRef_B_Fit3", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data double B", "MC double B"])
#
#PlotFits([H_FatJet_pTRef_Response_pT1_Trig.BData, H_FatJet_pTRef_Response_pT1_Trig.BMC],"pTRef_B_Fit1_Trig", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data double B", "MC double B"])
#PlotFits([H_FatJet_pTRef_Response_pT2_Trig.BData, H_FatJet_pTRef_Response_pT2_Trig.BMC],"pTRef_B_Fit2_Trig", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data double B", "MC double B"])
#PlotFits([H_FatJet_pTRef_Response_pT3_Trig.BData, H_FatJet_pTRef_Response_pT3_Trig.BMC],"pTRef_B_Fit3_Trig", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], ["Data double B", "MC double B"])
#
DrawTool.SetCompData(True)

#H_LeadB_pT.DrawMCTot()
#H_SecB_pT.DrawMCTot()

H_CloseJet_Eta.Draw()
H_CloseJet_Phi.Draw()
H_CloseJet_pT.Draw()
H_dR_CloseJet_TrigJet.Draw()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=0.2 jet before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Muon matched")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction correction applied")
H_MuJet_pT.Draw()
DrawTool.ClearTexts()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=0.2 jet before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Non muon matched")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction correction applied")
H_NonMuJet_pT.Draw()
DrawTool.ClearTexts()

H_dR_Mu_MuJet.Draw()
H_Muon_pT.Draw()
H_ThJet_pT.Draw()
H_MuJet_Eta.Draw()
H_NonMuJet_Eta.Draw()
H_MuJet_MV1.Draw()
H_NonMuJet_MV1.Draw()
H_MuJet_MV1_Tail.Draw()
H_NonMuJet_MV1_Tail.Draw()


DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "After b-tagging")
#DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Flavor fraction correction applied")
#DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "flavor fraction correction applied")
#DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "b-tagging scale factor applied")
H_dR.Draw()
DrawTool.ClearTexts()


H_dR_13.Draw()
H_dR_23.Draw()
H_dR_LeadB_Jet3.Draw()
H_dR_SecB_Jet3.Draw()
H_FatJet_pT.Draw()
H_TrigJet_pT_Weight.Draw()
H_FatJet_Mass.Draw()
H_FatJet_Mass_pT100.Draw()
H_FatJet_Mass_pT300.Draw()
H_FatJet_Mass_pT450.Draw()
H_FatJet_NJet.Draw()
H_FatJet_NCluster.Draw()
H_FatJet_Eta.Draw()
H_FatJet_MOPT.Draw()
H_FatJet_TrigJet_dPhi.Draw()
H_FatJet_pTRef.Draw()
H_FatJet_pTRef_Response.Draw()
H_FatJet_pTRef_Response_pT1.Draw()
H_FatJet_pTRef_Response_pT2.Draw()
H_FatJet_pTRef_Response_pT3.Draw()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
H_FatJet_tau21WTA.Draw()
H_FatJet_tau32.Draw()
H_FatJet_tau21.Draw()
H_FatJet_tau31.Draw()
H_FatJet_d2.Draw()
DrawTool.ClearTexts()

#H_FatJet_c2.Draw()
H_FatJet_width.Draw()

H_FatJet_TrackMass.Draw()
H_FatJet_TrackMass_Ratio.Draw()
H_FatJet_nTrack.Draw()
H_FatJet_nTrack_pT1.Draw()
H_FatJet_nTrack_pT2.Draw()
H_FatJet_nTrack_pT3.Draw()

H_FatJet_TrackD2_Ratio_pT1.Draw()
H_FatJet_TrackD2_Ratio_pT2.Draw()
H_FatJet_TrackD2_Ratio_pT3.Draw()

DrawTool.DrawHists("pTRef_pT1", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], [H_FatJet_pTRef_Response_pT1.MC, H_FatJet_pTRef_Response_pT1.Data], ["MC", "Data"])
DrawTool.DrawHists("pTRef_pT2", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], [H_FatJet_pTRef_Response_pT2.MC, H_FatJet_pTRef_Response_pT2.Data], ["MC", "Data"])
DrawTool.DrawHists("pTRef_pT3", ["Fat Jet p_{T}/ Ref p_{T}", "NJet"], [H_FatJet_pTRef_Response_pT3.MC, H_FatJet_pTRef_Response_pT3.Data], ["MC", "Data"])

#Compute Systematics
dR_Sys        = ComputeSys(H_dR.MC,          
                          [H_dR_JESUp.MC],    
                          [H_dR_JESDn.MC], ["Sys"])
dR_BSys       = ComputeSys(H_dR.BMC,          
                          [H_dR_BUp.BMC, H_dR_JESUp.BMC],    
                          [H_dR_BDn.BMC, H_dR_JESDn.BMC], ["BSys", "JESSys"])
dR_BSys_Only  = ComputeSys(H_dR.BMC,          
                          [H_dR_BUp.BMC],    
                          [H_dR_BDn.BMC], ["BSys"])

H_dR_FitUp, H_dR_FitDn, H_dR_BFitUp, H_dR_BFitDn   = ComputeFitSys(H_dR, H_dR_FitSys, CorMatrices, FitBinUncert)

dR_FitSys     = ComputeSys(H_dR.MC, 
                           [H_dR_FitUp], 
                           [H_dR_FitDn], ["FitSys"])

dR_BFitSys     = ComputeSys(H_dR.BMC, 
                            [H_dR_BFitUp], 
                            [H_dR_BFitDn], ["BFitSys"])

dR_TotSys     = ComputeSys(H_dR.MC, 
                           [H_dR_FitUp, H_dR_JESUp.MC], 
                           [H_dR_FitDn, H_dR_JESDn.MC], ["FitSys", "JES"])

dR_BTotSys     = ComputeSys(H_dR.BMC, 
                            [H_dR_BFitUp, H_dR_BUp.BMC, H_dR_JESUp.BMC],    
                            [H_dR_BFitDn, H_dR_BDn.BMC, H_dR_JESDn.BMC], ["BFitSys", "JESSys", "BJES"])

### Fat Jet pT

FatJet_JESSys = ComputeSys(H_FatJet_pT.MC, 
                          [H_FatJet_JESUp.MC], 
                          [H_FatJet_JESDn.MC], ["JES"])

FatJet_pT_Sys   = ComputeSys(H_FatJet_pT.BMC,   
                            [H_FatJet_pT_BUp.BMC, H_FatJet_JESUp.BMC],    
                            [H_FatJet_pT_BDn.BMC, H_FatJet_JESDn.BMC], ["BSys", "JES"])

FatJet_pT_BSys_Only   = ComputeSys(H_FatJet_pT.BMC,   
                                   [H_FatJet_pT_BUp.BMC],    
                                   [H_FatJet_pT_BUp.BMC], ["BSys"])

H_FatJet_pT_FitUp, H_FatJet_pT_FitDn, H_FatJet_pT_BFitUp, H_FatJet_pT_BFitDn   = ComputeFitSys(H_FatJet_pT, H_FatJet_pT_FitSys, CorMatrices, FitBinUncert)


FatJet_pT_FitSys  = ComputeSys(H_FatJet_pT.MC, 
                               [H_FatJet_pT_FitUp], 
                               [H_FatJet_pT_FitDn], ["FitSys"])
FatJet_pT_BFitSys = ComputeSys(H_FatJet_pT.BMC, 
                               [H_FatJet_pT_BFitUp], 
                               [H_FatJet_pT_BFitDn], ["BFitSys"])

FatJet_pT_TotSys  = ComputeSys(H_FatJet_pT.MC, 
                               [H_FatJet_JESUp.MC, H_FatJet_pT_FitUp], 
                               [H_FatJet_JESDn.MC, H_FatJet_pT_FitDn], ["JES", "FitSys"])

FatJet_pT_BTotSys   = ComputeSys(H_FatJet_pT.BMC,   
                                 [H_FatJet_pT_BUp.BMC, H_FatJet_JESUp.BMC, H_FatJet_pT_BFitUp],
                                 [H_FatJet_pT_BDn.BMC, H_FatJet_JESDn.BMC, H_FatJet_pT_BFitDn], ["BSys", "JES", "BFitSys"])


### Fat Jet Mass

FatJet_JMSSys = ComputeSys(H_FatJet_Mass.MC, 
                          [H_FatJet_JMSUp.MC], 
                          [H_FatJet_JMSDn.MC], ["JMS"])

H_FatJet_Mass_FitUp, H_FatJet_Mass_FitDn, H_FatJet_Mass_BFitUp, H_FatJet_Mass_BFitDn  = ComputeFitSys(H_FatJet_Mass, H_FatJet_Mass_FitSys, CorMatrices, FitBinUncert)

FatJet_Mass_FitSys  = ComputeSys(H_FatJet_Mass.MC, 
                                 [H_FatJet_Mass_FitUp], 
                                 [H_FatJet_Mass_FitDn], ["FitSys"])

FatJet_Mass_BFitSys  = ComputeSys(H_FatJet_Mass.BMC, 
                                  [H_FatJet_Mass_BFitUp], 
                                  [H_FatJet_Mass_BFitDn], ["BFitSys"])

FatJet_Mass_TotSys  = ComputeSys(H_FatJet_Mass.MC, 
                                 [H_FatJet_JMSUp.MC, H_FatJet_Mass_FitUp], 
                                 [H_FatJet_JMSDn.MC, H_FatJet_Mass_FitDn], ["JMS", "FitSys"])

FatJet_Mass_BTotSys   = ComputeSys(H_FatJet_Mass.BMC,   
                                   [H_FatJet_Mass_BUp.BMC, H_FatJet_JMSUp.BMC, H_FatJet_Mass_BFitUp],
                                   [H_FatJet_Mass_BDn.BMC, H_FatJet_JMSDn.BMC, H_FatJet_Mass_BFitDn], ["BSys", "JMS", "BFitSys"])


FatJet_pTRef_Response_Sys = ComputeSys(H_FatJet_pTRef_Response.MC, 
                                      [H_FatJet_pTRef_Response_JESUp.MC],
                                      [H_FatJet_pTRef_Response_JESDn.MC], ["JES"])
FatJet_pTRef_Response_BSys = ComputeSys(H_FatJet_pTRef_Response.BMC, 
                                      [H_FatJet_pTRef_Response_JESUp.BMC, H_FatJet_pTRef_Response_BUp.BMC],
                                      [H_FatJet_pTRef_Response_JESDn.BMC, H_FatJet_pTRef_Response_BDn.BMC], ["JES, BSys"])
FatJet_pTRef_Response_pT1_Sys = ComputeSys(H_FatJet_pTRef_Response_pT1.MC, 
                                      [H_FatJet_pTRef_Response_pT1_TJESUp.MC],
                                      [H_FatJet_pTRef_Response_pT1_TJESDn.MC], ["JES"])
FatJet_pTRef_Response_pT1_BSys = ComputeSys(H_FatJet_pTRef_Response_pT1.BMC, 
                                      [H_FatJet_pTRef_Response_pT1_TJESUp.BMC, H_FatJet_pTRef_Response_pT1_BUp.BMC],
                                      [H_FatJet_pTRef_Response_pT1_TJESDn.BMC, H_FatJet_pTRef_Response_pT1_BDn.BMC], ["JES, BSys"])
FatJet_pTRef_Response_pT2_Sys = ComputeSys(H_FatJet_pTRef_Response_pT2.MC, 
                                      [H_FatJet_pTRef_Response_pT2_TJESUp.MC],
                                      [H_FatJet_pTRef_Response_pT2_TJESDn.MC], ["JES"])
FatJet_pTRef_Response_pT2_BSys = ComputeSys(H_FatJet_pTRef_Response_pT2.BMC, 
                                      [H_FatJet_pTRef_Response_pT2_TJESUp.BMC, H_FatJet_pTRef_Response_pT2_BUp.BMC],
                                      [H_FatJet_pTRef_Response_pT2_TJESDn.BMC, H_FatJet_pTRef_Response_pT2_BDn.BMC], ["JES, BSys"])
FatJet_pTRef_Response_pT3_Sys = ComputeSys(H_FatJet_pTRef_Response_pT3.MC, 
                                      [H_FatJet_pTRef_Response_pT3_TJESUp.MC],
                                      [H_FatJet_pTRef_Response_pT3_TJESDn.MC], ["JES"])

FatJet_pTRef_Response_pT3_BSys = ComputeSys(H_FatJet_pTRef_Response_pT3.BMC, 
                                      [H_FatJet_pTRef_Response_pT3_TJESUp.BMC, H_FatJet_pTRef_Response_pT3_BUp.BMC],
                                      [H_FatJet_pTRef_Response_pT3_TJESDn.BMC, H_FatJet_pTRef_Response_pT3_BDn.BMC], ["JES, BSys"])

FatJet_pTRef_Response_pT1_Trig_Sys = ComputeSys(H_FatJet_pTRef_Response_pT1_Trig.MC, 
                                      [H_FatJet_pTRef_Response_pT1_TJESUp_Trig.MC],
                                      [H_FatJet_pTRef_Response_pT1_TJESDn_Trig.MC], ["JES"])
FatJet_pTRef_Response_pT1_Trig_BSys = ComputeSys(H_FatJet_pTRef_Response_pT1_Trig.BMC, 
                                      [H_FatJet_pTRef_Response_pT1_TJESUp_Trig.BMC, H_FatJet_pTRef_Response_pT1_BUp_Trig.BMC],
                                      [H_FatJet_pTRef_Response_pT1_TJESDn_Trig.BMC, H_FatJet_pTRef_Response_pT1_BDn_Trig.BMC], ["JES, BSys"])
FatJet_pTRef_Response_pT2_Trig_Sys = ComputeSys(H_FatJet_pTRef_Response_pT2_Trig.MC, 
                                      [H_FatJet_pTRef_Response_pT2_TJESUp_Trig.MC],
                                      [H_FatJet_pTRef_Response_pT2_TJESDn_Trig.MC], ["JES"])
FatJet_pTRef_Response_pT2_Trig_BSys = ComputeSys(H_FatJet_pTRef_Response_pT2_Trig.BMC, 
                                      [H_FatJet_pTRef_Response_pT2_TJESUp_Trig.BMC, H_FatJet_pTRef_Response_pT2_BUp_Trig.BMC],
                                      [H_FatJet_pTRef_Response_pT2_TJESDn_Trig.BMC, H_FatJet_pTRef_Response_pT2_BDn_Trig.BMC], ["JES, BSys"])
FatJet_pTRef_Response_pT3_Trig_Sys = ComputeSys(H_FatJet_pTRef_Response_pT3_Trig.MC, 
                                      [H_FatJet_pTRef_Response_pT3_TJESUp_Trig.MC],
                                      [H_FatJet_pTRef_Response_pT3_TJESDn_Trig.MC], ["JES"])
FatJet_pTRef_Response_pT3_Trig_BSys = ComputeSys(H_FatJet_pTRef_Response_pT3_Trig.BMC, 
                                      [H_FatJet_pTRef_Response_pT3_TJESUp_Trig.BMC, H_FatJet_pTRef_Response_pT3_BUp_Trig.BMC],
                                      [H_FatJet_pTRef_Response_pT3_TJESDn_Trig.BMC, H_FatJet_pTRef_Response_pT3_BDn_Trig.BMC], ["JES, BSys"])


### Mu Jet pT
MuJet_pT_Sys    = ComputeSys(H_MuJet_pT.MC,    
                            [H_MuJet_pT_JESUp.MC],    
                            [H_MuJet_pT_JESDn.MC], ["JES"])

H_MuJet_pT_FitUp, H_MuJet_pT_FitDn, H_MuJet_pT_BFitUp, H_MuJet_pT_BFitDn = ComputeFitSys(H_MuJet_pT, H_MuJet_pT_FitSys, CorMatrices, FitBinUncert)

MuJet_pT_FitSys  = ComputeSys(H_MuJet_pT.MC, 
                               [H_MuJet_pT_FitUp], 
                               [H_MuJet_pT_FitDn], ["FitSys"])
MuJet_pT_BFitSys  = ComputeSys(H_MuJet_pT.BMC, 
                               [H_MuJet_pT_BFitUp], 
                               [H_MuJet_pT_BFitDn], ["BFitSys"])

MuJet_pT_BSys    = ComputeSys(H_MuJet_pT.BMC,    
                             [H_MuJet_pT_BUp.BMC, H_MuJet_pT_JESUp.BMC],    
                             [H_MuJet_pT_BDn.BMC, H_MuJet_pT_JESDn.BMC], ["BSys", "JES"])

MuJet_pT_BSys_Only = ComputeSys(H_MuJet_pT.BMC,    
                                [H_MuJet_pT_BUp.BMC],    
                                [H_MuJet_pT_BDn.BMC], ["BSys"])

MuJet_pT_TotSys  = ComputeSys(H_MuJet_pT.MC, 
                              [H_MuJet_pT_JESUp.MC, H_MuJet_pT_FitUp], 
                              [H_MuJet_pT_JESDn.MC, H_MuJet_pT_FitDn], ["JES", "FitSys"])

MuJet_pT_BTotSys  = ComputeSys(H_MuJet_pT.BMC,   
                               [H_MuJet_pT_BUp.BMC, H_MuJet_pT_JESUp.BMC, H_MuJet_pT_BFitUp],
                               [H_MuJet_pT_BDn.BMC, H_MuJet_pT_JESDn.BMC, H_MuJet_pT_BFitDn], ["BSys", "JES", "BFitSys"])


### Non Mu jet pT

NonMuJet_pT_BSys = ComputeSys(H_NonMuJet_pT.BMC, 
                              [H_NonMuJet_pT_BUp.BMC, H_NonMuJet_pT_CUp.BMC, H_NonMuJet_pT_LUp.BMC], 
                              [H_NonMuJet_pT_BDn.BMC, H_NonMuJet_pT_CDn.BMC, H_NonMuJet_pT_LDn.BMC], ["BSys", "CSys", "LSys"])

H_NonMuJet_pT_FitUp, H_NonMuJet_pT_FitDn, H_NonMuJet_pT_BFitUp, H_NonMuJet_pT_BFitDn   = ComputeFitSys(H_NonMuJet_pT, H_NonMuJet_pT_FitSys, CorMatrices, FitBinUncert)

NonMuJet_pT_FitSys  = ComputeSys(H_NonMuJet_pT.MC, 
                                 [H_NonMuJet_pT_FitUp], 
                                 [H_NonMuJet_pT_FitDn], ["FitSys"])
NonMuJet_pT_BFitSys  = ComputeSys(H_NonMuJet_pT.BMC, 
                                 [H_NonMuJet_pT_BFitUp], 
                                 [H_NonMuJet_pT_BFitDn], ["BFitSys"])

NonMuJet_pT_BSys_Only = ComputeSys(H_NonMuJet_pT.BMC,    
                                   [H_NonMuJet_pT_BUp.BMC],    
                                   [H_NonMuJet_pT_BDn.BMC], ["BSys"])


NonMuJet_pT_TotSys  = ComputeSys(H_NonMuJet_pT.MC, 
                                 [H_NonMuJet_pT_FitUp], 
                                 [H_NonMuJet_pT_FitDn], ["FitSys"])

NonMuJet_pT_BTotSys  = ComputeSys(H_NonMuJet_pT.BMC,   
                                  [H_NonMuJet_pT_BUp.BMC, H_NonMuJet_pT_BFitUp],
                                  [H_NonMuJet_pT_BDn.BMC, H_NonMuJet_pT_BFitDn], ["BSys", "BFitSys"])


### Muon


Muon_pT_BSys    = ComputeSys(H_Muon_pT.BMC,    
                            [H_Muon_pT_BUp.BMC],    
                            [H_Muon_pT_BDn.BMC], ["BSys"])

Muon_pT_BSys_Only    = ComputeSys(H_Muon_pT.BMC,    
                                  [H_Muon_pT_BUp.BMC],    
                                  [H_Muon_pT_BDn.BMC], ["BSys"])

H_Muon_pT_FitUp, H_Muon_pT_FitDn, H_Muon_pT_BFitUp, H_Muon_pT_BFitDn = ComputeFitSys(H_Muon_pT, H_Muon_pT_FitSys, CorMatrices, FitBinUncert)


Muon_pT_FitSys  = ComputeSys(H_Muon_pT.MC, 
                               [H_Muon_pT_FitUp], 
                               [H_Muon_pT_FitDn], ["FitSys"])
Muon_pT_BFitSys  = ComputeSys(H_Muon_pT.BMC, 
                               [H_Muon_pT_BFitUp], 
                               [H_Muon_pT_BFitDn], ["BFitSys"])

Muon_pT_TotSys  = ComputeSys(H_Muon_pT.MC, 
                              [H_Muon_pT_JESUp.MC, H_Muon_pT_FitUp], 
                              [H_Muon_pT_JESDn.MC, H_Muon_pT_FitDn], ["JES", "FitSys"])

Muon_pT_BTotSys  = ComputeSys(H_Muon_pT.BMC,   
                               [H_Muon_pT_BUp.BMC, H_Muon_pT_JESUp.BMC, H_Muon_pT_BFitUp],
                               [H_Muon_pT_BDn.BMC, H_Muon_pT_JESDn.BMC, H_Muon_pT_BFitDn], ["BSys", "JES", "BFitSys"])
####

MuJet_Eta_Sys    = ComputeSys(H_MuJet_Eta.BMC,    
                            [H_MuJet_Eta_BUp.BMC],    
                            [H_MuJet_Eta_BDn.BMC], ["BSys"])

NonMuJet_Eta_Sys = ComputeSys(H_NonMuJet_Eta.BMC, 
                            [H_NonMuJet_Eta_BUp.BMC], 
                            [H_NonMuJet_Eta_BDn.BMC], ["BSys"])

dR_Mu_MuJet_Sys = ComputeSys(H_dR_Mu_MuJet.BMC,    
                            [H_dR_Mu_MuJet_BUp.BMC],    
                            [H_dR_Mu_MuJet_BDn.BMC], ["BSys"])

FatJet_Mass_Sys = ComputeSys(H_FatJet_Mass.BMC,  
                            [H_FatJet_Mass_BUp.BMC, H_FatJet_JMSUp.BMC, H_FatJet_Mass_JESUp.BMC],    
                            [H_FatJet_Mass_BDn.BMC, H_FatJet_JMSDn.BMC, H_FatJet_Mass_JESDn.BMC], ["BSys", "JMS", "JES"])

FatJet_Mass_BSys_Only = ComputeSys(H_FatJet_Mass.BMC,  
                                   [H_FatJet_Mass_BUp.BMC],    
                                   [H_FatJet_Mass_BDn.BMC], ["BSys"])

FatJet_Eta_Sys  = ComputeSys(H_FatJet_Eta.MC,   
                            [H_FatJet_Eta_JESUp.MC],    
                            [H_FatJet_Eta_JESDn.MC], ["JESSys"])

### Fat Jet Eta
H_FatJet_Eta_FitUp, H_FatJet_Eta_FitDn, H_FatJet_Eta_BFitUp, H_FatJet_Eta_BFitDn   = ComputeFitSys(H_FatJet_Eta, H_FatJet_Eta_FitSys, CorMatrices, FitBinUncert)


FatJet_Eta_FitSys  = ComputeSys(H_FatJet_Eta.MC, 
                               [H_FatJet_Eta_FitUp], 
                               [H_FatJet_Eta_FitDn], ["FitSys"])
FatJet_Eta_BFitSys = ComputeSys(H_FatJet_Eta.BMC, 
                               [H_FatJet_Eta_BFitUp], 
                               [H_FatJet_Eta_BFitDn], ["BFitSys"])

FatJet_Eta_TotSys  = ComputeSys(H_FatJet_Eta.MC, 
                               [H_FatJet_Eta_JESUp.MC, H_FatJet_Eta_FitUp], 
                               [H_FatJet_Eta_JESDn.MC, H_FatJet_Eta_FitDn], ["JES", "FitSys"])

FatJet_Eta_BTotSys   = ComputeSys(H_FatJet_Eta.BMC,   
                                 [H_FatJet_Eta_BUp.BMC, H_FatJet_Eta_JESUp.BMC, H_FatJet_Eta_BFitUp],
                                 [H_FatJet_Eta_BDn.BMC, H_FatJet_Eta_JESDn.BMC, H_FatJet_Eta_BFitDn], ["BSys", "JES", "BFitSys"])


FatJet_TrackMass_Ratio_Sys     = ComputeSys(H_FatJet_TrackMass_Ratio.BMC,   
                                           [H_FatJet_TrackMass_Ratio_BUp.BMC, H_FatJet_TrackMass_Ratio_JMSUp.BMC],    
                                           [H_FatJet_TrackMass_Ratio_BDn.BMC, H_FatJet_TrackMass_Ratio_JMSDn.BMC], ["BSys", "JMS"])
FatJet_TrackMass_Ratio_JMSSys  = ComputeSys(H_FatJet_TrackMass_Ratio.MC, 
                                           [H_FatJet_TrackMass_Ratio_JMSUp.MC],  
                                           [H_FatJet_TrackMass_Ratio_JMSDn.MC], ["JMS"])

FatJet_TrackpT_Ratio_Sys     = ComputeSys(H_FatJet_TrackpT_Ratio.MC,   
                                         [H_FatJet_TrackpT_Ratio_FJESUp.MC],    
                                         [H_FatJet_TrackpT_Ratio_FJESDn.MC], ["JES"])
FatJet_TrackpT_Ratio_BSys     = ComputeSys(H_FatJet_TrackpT_Ratio.BMC,   
                                          [H_FatJet_TrackpT_Ratio_BUp.BMC, H_FatJet_TrackpT_Ratio_FJESUp.BMC],    
                                          [H_FatJet_TrackpT_Ratio_BDn.BMC, H_FatJet_TrackpT_Ratio_FJESDn.BMC], ["BSys", "JES"])

FatJet_TrackpT_Ratio_pT1_Sys     = ComputeSys(H_FatJet_TrackpT_Ratio_pT1.MC,   
                                             [H_FatJet_TrackpT_Ratio_FJESDn_pT1.MC],    
                                             [H_FatJet_TrackpT_Ratio_FJESDn_pT1.MC], ["JES"])
FatJet_TrackpT_Ratio_pT1_BSys     = ComputeSys(H_FatJet_TrackpT_Ratio_pT1.BMC,   
                                              [H_FatJet_TrackpT_Ratio_BUp_pT1.BMC, H_FatJet_TrackpT_Ratio_FJESDn_pT1.BMC],    
                                              [H_FatJet_TrackpT_Ratio_BDn_pT1.BMC, H_FatJet_TrackpT_Ratio_FJESDn_pT1.BMC], ["BSys", "JES"])
FatJet_TrackpT_Ratio_pT2_Sys     = ComputeSys(H_FatJet_TrackpT_Ratio_pT2.MC,   
                                             [H_FatJet_TrackpT_Ratio_FJESDn_pT2.MC],    
                                             [H_FatJet_TrackpT_Ratio_FJESDn_pT2.MC], ["JES"])
FatJet_TrackpT_Ratio_pT2_BSys     = ComputeSys(H_FatJet_TrackpT_Ratio_pT2.BMC,   
                                              [H_FatJet_TrackpT_Ratio_BUp_pT2.BMC, H_FatJet_TrackpT_Ratio_FJESDn_pT2.BMC],    
                                              [H_FatJet_TrackpT_Ratio_BDn_pT2.BMC, H_FatJet_TrackpT_Ratio_FJESDn_pT2.BMC], ["BSys", "JES"])
FatJet_TrackpT_Ratio_pT3_Sys     = ComputeSys(H_FatJet_TrackpT_Ratio_pT3.MC,   
                                             [H_FatJet_TrackpT_Ratio_FJESDn_pT3.MC],    
                                             [H_FatJet_TrackpT_Ratio_FJESDn_pT3.MC], ["JES"])
FatJet_TrackpT_Ratio_pT3_BSys     = ComputeSys(H_FatJet_TrackpT_Ratio_pT3.BMC,   
                                              [H_FatJet_TrackpT_Ratio_BUp_pT3.BMC, H_FatJet_TrackpT_Ratio_FJESDn_pT3.BMC],    
                                              [H_FatJet_TrackpT_Ratio_BDn_pT3.BMC, H_FatJet_TrackpT_Ratio_FJESDn_pT3.BMC], ["BSys", "JES"])


FatJet_NJet_Sys    = ComputeSys(H_FatJet_NJet.BMC,   
                               [H_FatJet_NJet_BUp.BMC, H_FatJet_NJet_JESUp.BMC],    
                               [H_FatJet_NJet_BDn.BMC, H_FatJet_NJet_JESDn.BMC], ["BSys", "JES"])

FatJet_tau31_Sys   = ComputeSys(H_FatJet_tau31.MC,   
                               [H_FatJet_tau31_JESUp.MC],    
                               [H_FatJet_tau31_JESDn.MC], ["JES"])

FatJet_tau31_BSys   = ComputeSys(H_FatJet_tau31.BMC,   
                                [H_FatJet_tau31_BUp.BMC, H_FatJet_tau31_JESUp.BMC],    
                                [H_FatJet_tau31_BDn.BMC, H_FatJet_tau31_JESDn.BMC], ["BSys", "JES"])

FatJet_tau32_Sys   = ComputeSys(H_FatJet_tau32.MC,   
                               [H_FatJet_tau32_JESUp.MC],    
                               [H_FatJet_tau32_JESDn.MC], ["JES"])
FatJet_tau32_BSys   = ComputeSys(H_FatJet_tau32.BMC,   
                               [H_FatJet_tau32_BUp.BMC, H_FatJet_tau32_JESUp.BMC],    
                               [H_FatJet_tau32_BDn.BMC, H_FatJet_tau32_JESDn.BMC], ["BSys", "JES"])

FatJet_tau21_Sys   = ComputeSys(H_FatJet_tau21.MC,   
                               [H_FatJet_tau21_JESUp.MC],    
                               [H_FatJet_tau21_JESDn.MC], ["JES"])
FatJet_tau21_BSys   = ComputeSys(H_FatJet_tau21.BMC,   
                               [H_FatJet_tau21_BUp.BMC, H_FatJet_tau21_JESUp.BMC],    
                               [H_FatJet_tau21_BDn.BMC, H_FatJet_tau21_JESDn.BMC], ["BSys", "JES"])


############# tau21WTA

FatJet_tau21WTA_Sys   = ComputeSys(H_FatJet_tau21WTA.MC,   
                                   [ H_FatJet_tau21WTA_JESUp.MC, H_FatJet_tau21WTA_tauUp.MC],    
                                   [ H_FatJet_tau21WTA_JESDn.MC, H_FatJet_tau21WTA_tauDn.MC], ["JES", "Tau21"])

FatJet_tau21WTA_BSys   = ComputeSys(H_FatJet_tau21WTA.BMC,   
                                   [H_FatJet_tau21WTA_BUp.BMC, H_FatJet_tau21WTA_JESUp.BMC, H_FatJet_tau21WTA_tauUp.BMC],    
                                   [H_FatJet_tau21WTA_BDn.BMC, H_FatJet_tau21WTA_JESDn.BMC, H_FatJet_tau21WTA_tauDn.BMC], ["BSys", "JES", "Tau21"])

FatJet_tau21WTA_BSys_Only   = ComputeSys(H_FatJet_tau21WTA.BMC,   
                                         [H_FatJet_tau21WTA_BUp.BMC],    
                                         [H_FatJet_tau21WTA_BDn.BMC], ["BSys"])

H_FatJet_tau21WTA_FitUp, H_FatJet_tau21WTA_FitDn, H_FatJet_tau21WTA_BFitUp, H_FatJet_tau21WTA_BFitDn   = ComputeFitSys(H_FatJet_tau21WTA, H_FatJet_tau21WTA_FitSys, CorMatrices, FitBinUncert)

FatJet_tau21WTA_FitSys  = ComputeSys(H_FatJet_tau21WTA.MC, 
                                     [H_FatJet_tau21WTA_FitUp], 
                                     [H_FatJet_tau21WTA_FitDn], ["FitSys"])
FatJet_tau21WTA_BFitSys  = ComputeSys(H_FatJet_tau21WTA.BMC, 
                                     [H_FatJet_tau21WTA_BFitUp], 
                                     [H_FatJet_tau21WTA_BFitDn], ["BFitSys"])

FatJet_tau21WTA_TotSys  = ComputeSys(H_FatJet_tau21WTA.MC, 
                                     [H_FatJet_tau21WTA_JESUp.MC, H_FatJet_tau21WTA_tauUp.MC,  H_FatJet_tau21WTA_FitUp], 
                                     [H_FatJet_tau21WTA_JESDn.MC, H_FatJet_tau21WTA_tauDn.MC,  H_FatJet_tau21WTA_FitDn], ["JES", "Tau21", "FitSys"])

FatJet_tau21WTA_BTotSys   = ComputeSys(H_FatJet_tau21WTA.BMC,   
                                       [H_FatJet_tau21WTA_BUp.BMC, H_FatJet_tau21WTA_JESUp.BMC, H_FatJet_tau21WTA_tauUp.BMC, H_FatJet_tau21WTA_BFitUp],    
                                       [H_FatJet_tau21WTA_BDn.BMC, H_FatJet_tau21WTA_JESDn.BMC, H_FatJet_tau21WTA_tauDn.BMC, H_FatJet_tau21WTA_BFitDn], ["BSys", "JES", "Tau21", "BFitSys"])

############# d2

FatJet_d2_BSys      = ComputeSys(H_FatJet_d2.BMC,   
                               [H_FatJet_d2_BUp.BMC, H_FatJet_d2_JESUp.BMC, H_FatJet_d2_d2Up.BMC],    
                               [H_FatJet_d2_BDn.BMC, H_FatJet_d2_JESDn.BMC, H_FatJet_d2_d2Dn.BMC], ["BSys", "JES", "d2"])

FatJet_d2_BSys_Only   = ComputeSys(H_FatJet_d2.BMC,   
                                   [H_FatJet_d2_BUp.BMC],   
                                   [H_FatJet_d2_BDn.BMC], ["BSys"])

FatJet_d2_Sys      = ComputeSys(H_FatJet_d2.MC,   
                               [ H_FatJet_d2_JESUp.MC, H_FatJet_d2_d2Up.MC],    
                               [ H_FatJet_d2_JESDn.MC, H_FatJet_d2_d2Dn.MC], [ "JES", "d2"])

H_FatJet_d2_FitUp, H_FatJet_d2_FitDn, H_FatJet_d2_BFitUp, H_FatJet_d2_BFitDn   = ComputeFitSys(H_FatJet_d2, H_FatJet_d2_FitSys, CorMatrices, FitBinUncert)

FatJet_d2_FitSys  = ComputeSys(H_FatJet_d2.MC, 
                               [H_FatJet_d2_FitUp], 
                               [H_FatJet_d2_FitDn], ["FitSys"])
FatJet_d2_BFitSys  = ComputeSys(H_FatJet_d2.BMC, 
                               [H_FatJet_d2_BFitUp], 
                               [H_FatJet_d2_BFitDn], ["BFitSys"])

FatJet_d2_TotSys  = ComputeSys(H_FatJet_d2.MC, 
                               [H_FatJet_d2_JESUp.MC, H_FatJet_d2_d2Up.MC,  H_FatJet_d2_FitUp], 
                               [H_FatJet_d2_JESDn.MC, H_FatJet_d2_d2Dn.MC,  H_FatJet_d2_FitDn], ["JES", "D2", "FitSys"])

FatJet_d2_BTotSys   = ComputeSys(H_FatJet_d2.BMC,   
                                 [H_FatJet_d2_BUp.BMC, H_FatJet_d2_JESUp.BMC, H_FatJet_d2_d2Up.BMC, H_FatJet_d2_BFitUp],    
                                 [H_FatJet_d2_BDn.BMC, H_FatJet_d2_JESDn.BMC, H_FatJet_d2_d2Dn.BMC, H_FatJet_d2_BFitDn], ["BSys", "JES", "D2", "BFitSys"])

FatJet_MOPT_Sys     = ComputeSys(H_FatJet_MOPT.MC,   
                                [ H_FatJet_MOPT_JESUp.MC, H_FatJet_MOPT_JMSUp.MC],
                                [ H_FatJet_MOPT_JESDn.MC, H_FatJet_MOPT_JMSDn.MC], [ "JES", "JMS"])
FatJet_MOPT_BSys     = ComputeSys(H_FatJet_MOPT.BMC,   
                                 [ H_FatJet_MOPT_JESUp.BMC, H_FatJet_MOPT_JMSUp.BMC, H_FatJet_MOPT_BUp.BMC],
                                 [ H_FatJet_MOPT_JESDn.BMC, H_FatJet_MOPT_JMSDn.BMC, H_FatJet_MOPT_BDn.BMC], [ "JES", "JMS", "BUp"])




H_TrigJet_pT_Weight    =  MCData("TrigJet_pT_Weight",    DoubleFlav,   "TrigJet_pT_Weight",    "Trig Jet p_{T} [GeV]",   "Jets / Bin",   True)
H_TrigJet_pT_Weight_JESUp  =  MCData("TrigJet_pT_Weight_JESUp",    DoubleFlav,   "TrigJet_pT_Weight_JESUp",    "Trig Jet p_{T} [GeV]",   "Jets / Bin",   True)
H_TrigJet_pT_Weight_JESDn  =  MCData("TrigJet_pT_Weight_JESDn",    DoubleFlav,   "TrigJet_pT_Weight_JESDn",    "Trig Jet p_{T} [GeV]",   "Jets / Bin",   True)
H_TrigJet_pT_Weight_BUp  =  MCData("TrigJet_pT_Weight_BUp",    DoubleFlav,   "TrigJet_pT_Weight_BUp",    "Trig Jet p_{T} [GeV]",   "Jets / Bin",   True)
H_TrigJet_pT_Weight_BDn  =  MCData("TrigJet_pT_Weight_BDn",    DoubleFlav,   "TrigJet_pT_Weight_BDn",    "Trig Jet p_{T} [GeV]",   "Jets / Bin",   True)

TrigJet_pT_Weight_Sys = ComputeSys(H_TrigJet_pT_Weight.MC,
                                  [ H_TrigJet_pT_Weight_JESUp.MC],
                                  [ H_TrigJet_pT_Weight_JESDn.MC], ["JES"])
TrigJet_pT_Weight_BSys = ComputeSys(H_TrigJet_pT_Weight.BMC,
                                   [ H_TrigJet_pT_Weight_JESUp.BMC, H_TrigJet_pT_Weight_BUp.BMC],
                                   [ H_TrigJet_pT_Weight_JESDn.BMC, H_TrigJet_pT_Weight_BDn.BMC], ["JES", "BSys"])


#Draw systematics ####
H_FatJet_TrackMass_Ratio.DrawSys("JMS", FatJet_TrackMass_Ratio_JMSSys,  FatJet_TrackMass_Ratio_JMSSys)
H_FatJet_TrackMass_Ratio.DrawSys("BSys", FatJet_TrackMass_Ratio_Sys,  FatJet_TrackMass_Ratio_Sys)

H_FatJet_MOPT.DrawSys("JMS", FatJet_MOPT_Sys,  FatJet_MOPT_Sys)
H_FatJet_MOPT.DrawSys("BSys", FatJet_MOPT_BSys,  FatJet_MOPT_BSys)

H_FatJet_TrackpT_Ratio.DrawSys("JMS", FatJet_TrackpT_Ratio_Sys,  FatJet_TrackpT_Ratio_Sys)
H_FatJet_TrackpT_Ratio.DrawSys("BSys", FatJet_TrackpT_Ratio_BSys,  FatJet_TrackpT_Ratio_BSys)

H_FatJet_TrackpT_Ratio_pT1.DrawSys("JMS", FatJet_TrackpT_Ratio_pT1_Sys,  FatJet_TrackpT_Ratio_pT1_Sys)
H_FatJet_TrackpT_Ratio_pT1.DrawSys("BSys", FatJet_TrackpT_Ratio_pT1_BSys,  FatJet_TrackpT_Ratio_pT1_BSys)
H_FatJet_TrackpT_Ratio_pT2.DrawSys("JMS", FatJet_TrackpT_Ratio_pT2_Sys,  FatJet_TrackpT_Ratio_pT2_Sys)
H_FatJet_TrackpT_Ratio_pT2.DrawSys("BSys", FatJet_TrackpT_Ratio_pT2_BSys,  FatJet_TrackpT_Ratio_pT2_BSys)
H_FatJet_TrackpT_Ratio_pT3.DrawSys("JMS", FatJet_TrackpT_Ratio_pT3_Sys,  FatJet_TrackpT_Ratio_pT3_Sys)
H_FatJet_TrackpT_Ratio_pT3.DrawSys("BSys", FatJet_TrackpT_Ratio_pT3_BSys,  FatJet_TrackpT_Ratio_pT3_BSys)

#### JES

#FatJet_pT_TotSys.SetBinError(1, 0)
#FatJet_pT_FitSys.SetBinError(1, 0)
#FatJet_pT_BTotSys.SetBinError(1, 0)
#FatJet_pT_BFitSys.SetBinError(1, 0)
#FatJet_pT_BSys_Only.SetBinError(1, 0)

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_FatJet_pT.DrawSys("JES", FatJet_JESSys, FatJet_JESSys)
H_FatJet_pT.DrawSysMultiple("SysTot", [FatJet_pT_TotSys, FatJet_pT_FitSys], [FatJet_pT_BTotSys, FatJet_pT_BFitSys])
DrawTool.ClearTexts()

#H_FatJet_pT.DrawSys("FitSys",FatJet_pT_FitSys, FatJet_pT_FitSys)
#H_FatJet_pT.DrawSys("BSys",FatJet_pT_Sys, FatJet_pT_Sys)

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_FatJet_pT.DrawSysMultiple("BSysTot", [FatJet_pT_TotSys, FatJet_pT_FitSys], [FatJet_pT_BTotSys, FatJet_pT_BFitSys], FatJet_pT_BSys_Only)
DrawTool.ClearTexts()



#### JMS
DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_FatJet_Mass.DrawSys("JMS",FatJet_JMSSys, FatJet_JMSSys)
H_FatJet_Mass.DrawSysMultiple("SysTot", [FatJet_Mass_TotSys, FatJet_Mass_FitSys], [FatJet_Mass_BTotSys, FatJet_Mass_BFitSys])
DrawTool.ClearTexts()

H_FatJet_Mass.DrawSys("FitSys",FatJet_Mass_FitSys, FatJet_Mass_FitSys)

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_FatJet_Mass.DrawSys("BSys",FatJet_Mass_Sys, FatJet_Mass_Sys)
H_FatJet_Mass.DrawSysMultiple("BSysTot", [FatJet_Mass_TotSys, FatJet_Mass_FitSys], [FatJet_Mass_BTotSys, FatJet_Mass_BFitSys], FatJet_Mass_BSys_Only)
DrawTool.ClearTexts()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_FatJet_Eta.DrawSysMultiple("SysTot", [FatJet_Eta_TotSys, FatJet_Eta_FitSys], [FatJet_Eta_TotSys, FatJet_Eta_FitSys])
H_FatJet_Eta.DrawSys("BSys",FatJet_Eta_Sys, FatJet_Eta_Sys)
FatJet_Eta_Sys.Write()
DrawTool.ClearTexts()



### mu jet pT
H_MuJet_pT.DrawSys("Sys",    MuJet_pT_Sys,    MuJet_pT_Sys)
H_MuJet_pT.DrawSys("FitSys", MuJet_pT_FitSys, MuJet_pT_FitSys)

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=0.2 jet before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Muon matched")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_MuJet_pT.DrawSysMultiple("SysTot", [MuJet_pT_TotSys, MuJet_pT_FitSys], [MuJet_pT_BTotSys, MuJet_pT_BFitSys])
DrawTool.ClearTexts()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=0.2 jet after b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Muon matched")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_MuJet_pT.DrawSys("BSys",    MuJet_pT_BSys,    MuJet_pT_BSys)
H_MuJet_pT.DrawSysMultiple("BSysTot", [MuJet_pT_TotSys, MuJet_pT_FitSys], [MuJet_pT_BTotSys, MuJet_pT_BFitSys], MuJet_pT_BSys_Only )
DrawTool.ClearTexts()


### non mu jet pT
H_NonMuJet_pT.DrawSys("FitSys", NonMuJet_pT_FitSys, NonMuJet_pT_FitSys)

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=0.2 jet before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Non muon matched")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_NonMuJet_pT.DrawSysMultiple("SysTot", [NonMuJet_pT_TotSys, NonMuJet_pT_FitSys], [NonMuJet_pT_BTotSys, NonMuJet_pT_BFitSys])
DrawTool.ClearTexts()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=0.2 jet after b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Non muon matched")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_NonMuJet_pT.DrawSys("BSys",    NonMuJet_pT_BSys,    NonMuJet_pT_BSys)
H_NonMuJet_pT.DrawSysMultiple("BSysTot", [NonMuJet_pT_TotSys, NonMuJet_pT_FitSys], [NonMuJet_pT_BTotSys, NonMuJet_pT_BFitSys], NonMuJet_pT_BSys_Only)
DrawTool.ClearTexts()


#### muon pT

H_Muon_pT.DrawSys("BSys",    Muon_pT_BSys,    Muon_pT_BSys)
H_Muon_pT.DrawSys("FitSys",  Muon_pT_FitSys, Muon_pT_FitSys)

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "Muon matched to anti-k_{t} R=0.2 jet before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Flavor fraction corrected")
H_Muon_pT.DrawSysMultiple("SysTot", [Muon_pT_TotSys, Muon_pT_FitSys], [Muon_pT_BTotSys, Muon_pT_BFitSys])
DrawTool.ClearTexts()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "Muon matched to anti-k_{t} R=0.2 jet after b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Flavor fraction corrected")
H_Muon_pT.DrawSys("BSys",    Muon_pT_BSys,    Muon_pT_BSys)
H_Muon_pT.DrawSysMultiple("BSysTot", [Muon_pT_TotSys, Muon_pT_FitSys], [Muon_pT_BTotSys, Muon_pT_BFitSys], Muon_pT_BSys_Only)
DrawTool.ClearTexts()

##########
H_MuJet_Eta.DrawSys("BSys",    MuJet_Eta_Sys,    MuJet_Eta_Sys)
H_NonMuJet_Eta.DrawSys("BSys", NonMuJet_Eta_Sys, NonMuJet_Eta_Sys)
H_dR_Mu_MuJet.DrawSys("BSys",    dR_Mu_MuJet_Sys, dR_Mu_MuJet_Sys)


## dR mu jet and non mu jet
H_dR.DrawSys("Sys", dR_Sys, dR_Sys)
H_dR.DrawSys("FitSys", dR_FitSys, dR_FitSys)
H_dR.DrawSys("BFitSys", dR_BFitSys, dR_BFitSys)

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "Before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Flavor fraction corrected")
H_dR.DrawSysMultiple("SysTot", [dR_TotSys, dR_FitSys], [dR_BTotSys, dR_BFitSys])
DrawTool.ClearTexts()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "After b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Flavor fraction corrected")
H_dR.DrawSysMultiple("BSysTot", [dR_TotSys, dR_FitSys], [dR_BTotSys, dR_BFitSys], dR_BSys_Only)
H_dR.DrawSys("BSys", dR_BSys, dR_BSys)
DrawTool.ClearTexts()


#H_FatJet_pTRef_Response.DrawSys("JES", FatJet_pTRef_Response_Sys, FatJet_pTRef_Response_Sys)
#H_FatJet_pTRef_Response.DrawSys("BSys", FatJet_pTRef_Response_BSys, FatJet_pTRef_Response_BSys)
#H_FatJet_pTRef_Response_pT1.DrawSys("JES", FatJet_pTRef_Response_pT1_Sys, FatJet_pTRef_Response_pT1_Sys)
#H_FatJet_pTRef_Response_pT1.DrawSys("BSys", FatJet_pTRef_Response_pT1_BSys, FatJet_pTRef_Response_pT1_BSys)
#H_FatJet_pTRef_Response_pT2.DrawSys("JES", FatJet_pTRef_Response_pT2_Sys, FatJet_pTRef_Response_pT2_Sys)
#H_FatJet_pTRef_Response_pT2.DrawSys("BSys", FatJet_pTRef_Response_pT2_BSys, FatJet_pTRef_Response_pT2_BSys)
#H_FatJet_pTRef_Response_pT3.DrawSys("JES", FatJet_pTRef_Response_pT3_Sys, FatJet_pTRef_Response_pT3_Sys)
#H_FatJet_pTRef_Response_pT3.DrawSys("BSys", FatJet_pTRef_Response_pT3_BSys, FatJet_pTRef_Response_pT3_BSys)
#
#H_FatJet_pTRef_Response_pT1_Trig.DrawSys("JES", FatJet_pTRef_Response_pT1_Trig_Sys, FatJet_pTRef_Response_pT1_Trig_Sys)
#H_FatJet_pTRef_Response_pT1_Trig.DrawSys("BSys", FatJet_pTRef_Response_pT1_Trig_BSys, FatJet_pTRef_Response_pT1_Trig_BSys)
#H_FatJet_pTRef_Response_pT2_Trig.DrawSys("JES", FatJet_pTRef_Response_pT2_Trig_Sys, FatJet_pTRef_Response_pT2_Trig_Sys)
#H_FatJet_pTRef_Response_pT2_Trig.DrawSys("BSys", FatJet_pTRef_Response_pT2_Trig_BSys, FatJet_pTRef_Response_pT2_Trig_BSys)
#H_FatJet_pTRef_Response_pT3_Trig.DrawSys("JES", FatJet_pTRef_Response_pT3_Trig_Sys, FatJet_pTRef_Response_pT3_Trig_Sys)
#H_FatJet_pTRef_Response_pT3_Trig.DrawSys("BSys", FatJet_pTRef_Response_pT3_Trig_BSys, FatJet_pTRef_Response_pT3_Trig_BSys)

H_FatJet_NJet.DrawSys("BSys",FatJet_NJet_Sys, FatJet_NJet_Sys)

H_FatJet_tau31.DrawSys("Sys",FatJet_tau31_Sys, FatJet_tau31_Sys)
H_FatJet_tau21.DrawSys("Sys",FatJet_tau21_Sys, FatJet_tau21_Sys)
H_FatJet_tau32.DrawSys("Sys",FatJet_tau32_Sys, FatJet_tau32_Sys)

H_FatJet_tau31.DrawSys("BSys",FatJet_tau31_BSys, FatJet_tau31_BSys)
H_FatJet_tau21.DrawSys("BSys",FatJet_tau21_BSys, FatJet_tau21_BSys)
H_FatJet_tau32.DrawSys("BSys",FatJet_tau32_BSys, FatJet_tau32_BSys)

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_FatJet_tau21WTA.DrawSys("Sys",FatJet_tau21WTA_Sys, FatJet_tau21WTA_Sys)
H_FatJet_tau21WTA.DrawSysMultiple("SysTot", [FatJet_tau21WTA_TotSys, FatJet_tau21WTA_FitSys], [FatJet_tau21WTA_BTotSys, FatJet_tau21WTA_BFitSys])
DrawTool.ClearTexts()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_FatJet_tau21WTA.DrawSys("BSys",FatJet_tau21WTA_BSys, FatJet_tau21WTA_BSys)
H_FatJet_tau21WTA.DrawSysMultiple("BSysTot", [FatJet_tau21WTA_TotSys, FatJet_tau21WTA_FitSys], [FatJet_tau21WTA_BTotSys, FatJet_tau21WTA_BFitSys], FatJet_tau21WTA_BSys_Only)
DrawTool.ClearTexts()

H_FatJet_tau21WTA.DrawSys("FitSys", FatJet_tau21WTA_FitSys, FatJet_tau21WTA_FitSys)
H_FatJet_tau21WTA.DrawSys("BFitSys", FatJet_tau21WTA_BFitSys, FatJet_tau21WTA_BFitSys)

########### d2 systematics
#DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet before b-tagging")
#DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet before b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
H_FatJet_d2.DrawSys("Sys",FatJet_d2_Sys, FatJet_d2_Sys)
H_FatJet_d2.DrawSysMultiple("SysTot", [FatJet_d2_TotSys, FatJet_d2_FitSys], [FatJet_d2_BTotSys, FatJet_d2_BFitSys])
DrawTool.ClearTexts()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.AddTexts(0.23, 0.59, 1, 0.048, "Flavor fraction corrected")
#DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging")
#DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
H_FatJet_d2.DrawSys("BSys",FatJet_d2_BSys, FatJet_d2_BSys)
H_FatJet_d2.DrawSysMultiple("BSysTot", [FatJet_d2_TotSys, FatJet_d2_FitSys], [FatJet_d2_BTotSys, FatJet_d2_BFitSys], FatJet_d2_BSys_Only)
DrawTool.ClearTexts()

H_FatJet_d2.DrawSys("FitSys", FatJet_d2_FitSys, FatJet_d2_FitSys )
H_FatJet_d2.DrawSys("BFitSys", FatJet_d2_BFitSys, FatJet_d2_BFitSys)


H_TrigJet_pT_Weight.DrawSys("Sys", TrigJet_pT_Weight_Sys, TrigJet_pT_Weight_Sys)
H_TrigJet_pT_Weight.DrawSys("BSys", TrigJet_pT_Weight_BSys, TrigJet_pT_Weight_BSys)


#### Draw Normalized Distributions #####

if (not isHerwig):
    H_FatJet_Mass_Response_NStack.MC.Rebin(10)
    H_FatJet_Mass_Response_NStack.BMC.Rebin(10)

DrawTool_Ratio.DrawHists("Mass_Response_Ratio", ["Response", "Arbitrary Units"], [H_FatJet_Mass_Response_NStack.MC, H_FatJet_Mass_Response_NStack.BMC], ["Inclusive", "Btag"])
DrawTool_Ratio.DrawHists("Mass_Response_C_Ratio_pT1_Eta1", ["Response", "Arbitrary Units"], [H_FatJet_Mass_Response_C_NStack_pT1_Eta1.MC, H_FatJet_Mass_Response_C_NStack_pT1_Eta1.BMC], ["Inclusive", "Btag"])
DrawTool_Ratio.DrawHists("Mass_Response_C_Ratio_pT1_Eta2", ["Response", "Arbitrary Units"], [H_FatJet_Mass_Response_C_NStack_pT1_Eta2.MC, H_FatJet_Mass_Response_C_NStack_pT1_Eta2.BMC], ["Inclusive", "Btag"])
DrawTool_Ratio.DrawHists("Mass_Response_C_Ratio_pT2_Eta1", ["", ""], [H_FatJet_Mass_Response_C_NStack_pT2_Eta1.MC, H_FatJet_Mass_Response_C_NStack_pT2_Eta1.BMC], ["Inclusive", "Btag"])
DrawTool_Ratio.DrawHists("Mass_Response_C_Ratio_pT2_Eta2", ["", ""], [H_FatJet_Mass_Response_C_NStack_pT2_Eta2.MC, H_FatJet_Mass_Response_C_NStack_pT2_Eta2.BMC], ["Inclusive", "Btag"])
DrawTool_Ratio.DrawHists("Mass_Response_C_Ratio_pT3_Eta1", ["", ""], [H_FatJet_Mass_Response_C_NStack_pT3_Eta1.MC, H_FatJet_Mass_Response_C_NStack_pT3_Eta1.BMC], ["Inclusive", "Btag"])
DrawTool_Ratio.DrawHists("Mass_Response_C_Ratio_pT3_Eta2", ["", ""], [H_FatJet_Mass_Response_C_NStack_pT3_Eta2.MC, H_FatJet_Mass_Response_C_NStack_pT3_Eta2.BMC], ["Inclusive", "Btag"])


DrawTool.SetCompData(False)

#Sd0Flav_s = ["BB", "BL", "BL", "LL", "CL", "CC"]
Sd0Flav_Rev = Sd0Flav
Sd0Flav_Rev.reverse()

H_MuJet_Sd0.Stack[0].SetMarkerStyle(20)
H_MuJet_Sd0.Stack[1].SetMarkerStyle(21)
H_MuJet_Sd0.Stack[2].SetMarkerStyle(22)
H_MuJet_Sd0.Stack[3].SetMarkerStyle(23)
H_MuJet_Sd0.Stack[4].SetMarkerStyle(33)
H_MuJet_Sd0.Stack[0].SetMarkerSize(1.5)
H_MuJet_Sd0.Stack[1].SetMarkerSize(1.5)
H_MuJet_Sd0.Stack[2].SetMarkerSize(1.5)
H_MuJet_Sd0.Stack[3].SetMarkerSize(1.5)
H_MuJet_Sd0.Stack[4].SetMarkerSize(1.5)
H_NonMuJet_Sd0.Stack[0].SetMarkerStyle(20)
H_NonMuJet_Sd0.Stack[1].SetMarkerStyle(21)
H_NonMuJet_Sd0.Stack[2].SetMarkerStyle(22)
H_NonMuJet_Sd0.Stack[3].SetMarkerStyle(23)
H_NonMuJet_Sd0.Stack[4].SetMarkerStyle(33)
H_NonMuJet_Sd0.Stack[0].SetMarkerSize(1.5)
H_NonMuJet_Sd0.Stack[1].SetMarkerSize(1.5)
H_NonMuJet_Sd0.Stack[2].SetMarkerSize(1.5)
H_NonMuJet_Sd0.Stack[3].SetMarkerSize(1.5)
H_NonMuJet_Sd0.Stack[4].SetMarkerSize(1.5)

for i in range(len(H_MuJet_Sd0.Stack)):
    H_MuJet_Sd0.Stack[i].GetYaxis().SetTitleOffset( H_MuJet_Sd0.Stack[i].GetYaxis().GetTitleOffset()*1.4)
    H_NonMuJet_Sd0.Stack[i].GetYaxis().SetTitleOffset( H_NonMuJet_Sd0.Stack[i].GetYaxis().GetTitleOffset()*1.4)

H_MuJet_Sd0_pT2.Stack[0].SetMarkerStyle(20)
H_MuJet_Sd0_pT2.Stack[1].SetMarkerStyle(21)
H_MuJet_Sd0_pT2.Stack[2].SetMarkerStyle(22)
H_MuJet_Sd0_pT2.Stack[3].SetMarkerStyle(23)
H_MuJet_Sd0_pT2.Stack[4].SetMarkerStyle(33)

H_NonMuJet_Sd0_pT2.Stack[0].SetMarkerStyle(20)
H_NonMuJet_Sd0_pT2.Stack[1].SetMarkerStyle(21)
H_NonMuJet_Sd0_pT2.Stack[2].SetMarkerStyle(22)
H_NonMuJet_Sd0_pT2.Stack[3].SetMarkerStyle(23)
H_NonMuJet_Sd0_pT2.Stack[4].SetMarkerStyle(33)

H_MuJet_Sd0_pT7.Stack[0].SetMarkerStyle(20)
H_MuJet_Sd0_pT7.Stack[1].SetMarkerStyle(21)
H_MuJet_Sd0_pT7.Stack[2].SetMarkerStyle(22)
H_MuJet_Sd0_pT7.Stack[3].SetMarkerStyle(23)
H_MuJet_Sd0_pT7.Stack[4].SetMarkerStyle(33)

H_NonMuJet_Sd0_pT7.Stack[0].SetMarkerStyle(20)
H_NonMuJet_Sd0_pT7.Stack[1].SetMarkerStyle(21)
H_NonMuJet_Sd0_pT7.Stack[2].SetMarkerStyle(22)
H_NonMuJet_Sd0_pT7.Stack[3].SetMarkerStyle(23)
H_NonMuJet_Sd0_pT7.Stack[4].SetMarkerStyle(33)


H_MuJet_Sd0_pT8.Stack[0].SetMarkerStyle(20)
H_MuJet_Sd0_pT8.Stack[1].SetMarkerStyle(21)
H_MuJet_Sd0_pT8.Stack[2].SetMarkerStyle(22)
H_MuJet_Sd0_pT8.Stack[3].SetMarkerStyle(23)
H_MuJet_Sd0_pT8.Stack[4].SetMarkerStyle(33)
H_NonMuJet_Sd0_pT8.Stack[0].SetMarkerStyle(20)
H_NonMuJet_Sd0_pT8.Stack[1].SetMarkerStyle(21)
H_NonMuJet_Sd0_pT8.Stack[2].SetMarkerStyle(22)
H_NonMuJet_Sd0_pT8.Stack[3].SetMarkerStyle(23)
H_NonMuJet_Sd0_pT8.Stack[4].SetMarkerStyle(33)

print 'marker size', H_NonMuJet_Sd0_pT8.Stack[4].GetMarkerSize()



DrawTool.lumi = "0"
DrawTool.colorlist = [DrawTool.origcolorlist[0], DrawTool.origcolorlist[1], DrawTool.origcolorlist[2], DrawTool.origcolorlist[3], DrawTool.origcolorlist[3], DrawTool.origcolorlist[5]]
#DrawTool.colorlist = DrawTool.boringcolor
DrawTool.shift = 0.18
DrawTool.studytype = "Simulation Preliminary"
DrawTool.AddTexts(0.41, 0.80, 1, 0.048, "Pythia 8")
DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=0.2 track jet")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Muon matched")
H_MuJet_Sd0.Stack = NormHists(H_MuJet_Sd0.Stack)
DrawTool.DrawHists("MuJet_Sd0",  ["S_{d_{0}}", "Arbitrary Unit"], [H_MuJet_Sd0.Stack[2], H_MuJet_Sd0.Stack[1], H_MuJet_Sd0.Stack[0], H_MuJet_Sd0.Stack[3], H_MuJet_Sd0.Stack[3], H_MuJet_Sd0.Stack[4]]
                   , ["LL+LC+LB", "CL+CB", "CC", "BL+BC", "BL+BC", "BB"])
DrawTool.ClearTexts()

DrawTool.AddTexts(0.41, 0.80, 1, 0.048, "Pythia 8")
DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=0.2 track jet")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Non muon matched")
H_NonMuJet_Sd0.Stack = NormHists(H_NonMuJet_Sd0.Stack)
DrawTool.DrawHists("NonMuJet_Sd0",  ["S_{d_{0}}", "Arbitrary Unit"], [H_NonMuJet_Sd0.Stack[2], H_NonMuJet_Sd0.Stack[1], H_NonMuJet_Sd0.Stack[0], H_NonMuJet_Sd0.Stack[3], H_NonMuJet_Sd0.Stack[3], 
                                                                      H_NonMuJet_Sd0.Stack[4]]
                   , ["LL+LC+LB", "CL+CB", "CC", "BL+BC", "BL+BC", "BB"])
DrawTool.ClearTexts()

DrawTool.colorlist = DrawTool.origcolorlist
DrawTool.lumi = "20.3"
DrawTool.shift = 0.135

H_MuJet_Sd0_pT1.Stack = NormHists(H_MuJet_Sd0_pT1.Stack)
DrawTool.DrawHists("MuJet_Sd0_pT1",  ["S_{d0}", "Arbitrary Unit"], H_MuJet_Sd0_pT1.Stack[0:5], Sd0Flav_Rev)
H_NonMuJet_Sd0_pT1.Stack = NormHists(H_NonMuJet_Sd0_pT1.Stack)
DrawTool.DrawHists("NonMuJet_Sd0_pT1",  ["S_{d0}", "Arbitrary Unit"], H_NonMuJet_Sd0_pT1.Stack[0:5], Sd0Flav_Rev)

H_MuJet_Sd0_pT2.Stack = NormHists(H_MuJet_Sd0_pT2.Stack)
DrawTool.DrawHists("MuJet_Sd0_pT2",  ["S_{d0}", "Arbitrary Unit"], H_MuJet_Sd0_pT2.Stack[0:5], Sd0Flav_Rev)
H_NonMuJet_Sd0_pT2.Stack = NormHists(H_NonMuJet_Sd0_pT2.Stack)
DrawTool.DrawHists("NonMuJet_Sd0_pT2",  ["S_{d0}", "Arbitrary Unit"], H_NonMuJet_Sd0_pT2.Stack[0:5], Sd0Flav_Rev)

H_MuJet_Sd0_pT3.Stack = NormHists(H_MuJet_Sd0_pT3.Stack)
DrawTool.DrawHists("MuJet_Sd0_pT3",  ["S_{d0}", "Arbitrary Unit"], H_MuJet_Sd0_pT3.Stack[0:5], Sd0Flav_Rev)
H_NonMuJet_Sd0_pT3.Stack = NormHists(H_NonMuJet_Sd0_pT3.Stack)
DrawTool.DrawHists("NonMuJet_Sd0_pT3",  ["S_{d0}", "Arbitrary Unit"], H_NonMuJet_Sd0_pT3.Stack[0:5], Sd0Flav_Rev)

H_MuJet_Sd0_pT4.Stack = NormHists(H_MuJet_Sd0_pT4.Stack)
DrawTool.DrawHists("MuJet_Sd0_pT4",  ["S_{d0}", "Arbitrary Unit"], H_MuJet_Sd0_pT4.Stack[0:5], Sd0Flav_Rev)
H_NonMuJet_Sd0_pT4.Stack = NormHists(H_NonMuJet_Sd0_pT4.Stack)
DrawTool.DrawHists("NonMuJet_Sd0_pT4",  ["S_{d0}", "Arbitrary Unit"], H_NonMuJet_Sd0_pT4.Stack[0:5], Sd0Flav_Rev)

H_MuJet_Sd0_pT5.Stack = NormHists(H_MuJet_Sd0_pT5.Stack)
DrawTool.DrawHists("MuJet_Sd0_pT5",  ["S_{d0}", "Arbitrary Unit"], H_MuJet_Sd0_pT5.Stack[0:5], Sd0Flav_Rev)
H_NonMuJet_Sd0_pT4.Stack = NormHists(H_NonMuJet_Sd0_pT5.Stack)
DrawTool.DrawHists("NonMuJet_Sd0_pT5",  ["S_{d0}", "Arbitrary Unit"], H_NonMuJet_Sd0_pT5.Stack[0:5], Sd0Flav_Rev)
 
H_MuJet_Sd0_pT6.Stack = NormHists(H_MuJet_Sd0_pT6.Stack)
DrawTool.DrawHists("MuJet_Sd0_pT6",  ["S_{d0}", "Arbitrary Unit"], H_MuJet_Sd0_pT6.Stack[0:5], Sd0Flav_Rev)
H_NonMuJet_Sd0_pT6.Stack = NormHists(H_NonMuJet_Sd0_pT6.Stack)
DrawTool.DrawHists("NonMuJet_Sd0_pT6",  ["S_{d0}", "Arbitrary Unit"], H_NonMuJet_Sd0_pT6.Stack[0:5], Sd0Flav_Rev)

H_MuJet_Sd0_pT7.Stack = NormHists(H_MuJet_Sd0_pT7.Stack)
DrawTool.DrawHists("MuJet_Sd0_pT7",  ["S_{d0}", "Arbitrary Unit"], H_MuJet_Sd0_pT7.Stack[0:5], Sd0Flav_Rev)
H_NonMuJet_Sd0_pT7.Stack = NormHists(H_NonMuJet_Sd0_pT7.Stack)
DrawTool.DrawHists("NonMuJet_Sd0_pT7",  ["S_{d0}", "Arbitrary Unit"], H_NonMuJet_Sd0_pT7.Stack[0:5], Sd0Flav_Rev)

H_MuJet_Sd0_pT8.Stack = NormHists(H_MuJet_Sd0_pT8.Stack)
DrawTool.DrawHists("MuJet_Sd0_pT8",  ["S_{d0}", "Arbitrary Unit"], H_MuJet_Sd0_pT8.Stack[0:5], Sd0Flav_Rev)
H_NonMuJet_Sd0_pT8.Stack = NormHists(H_NonMuJet_Sd0_pT8.Stack)
DrawTool.DrawHists("NonMuJet_Sd0_pT8",  ["S_{d0}", "Arbitrary Unit"], H_NonMuJet_Sd0_pT8.Stack[0:5], Sd0Flav_Rev)

H_MuJet_Sd0_pT9.Stack = NormHists(H_MuJet_Sd0_pT9.Stack)
DrawTool.DrawHists("MuJet_Sd0_pT9",  ["S_{d0}", "Arbitrary Unit"], H_MuJet_Sd0_pT9.Stack[0:5], Sd0Flav_Rev)
H_NonMuJet_Sd0_pT9.Stack = NormHists(H_NonMuJet_Sd0_pT9.Stack)
DrawTool.DrawHists("NonMuJet_Sd0_pT9",  ["S_{d0}", "Arbitrary Unit"], H_NonMuJet_Sd0_pT9.Stack[0:5], Sd0Flav_Rev)

sys.exit(0)
DrawTool.studytype = "Preliminary"


H_FatJet_TruthMass.DrawMC()


FatJetpTBin_L = [200,400,600,1000]
FatJetpTBin = array.array("d")
FatJetpTBin.fromlist(FatJetpTBin_L)

FatJetMOPTBin_L = [0, 0.25, 0.4, 1.5]
FatJetMOPTBin = array.array("d")
FatJetMOPTBin.fromlist(FatJetMOPTBin_L)

TrigJetpTBins_L = [100, 400, 550, 1000]
TrigJetpTBins = array.array("d")
TrigJetpTBins.fromlist(TrigJetpTBins_L)


FatJetResponseEtaBins_L = [-2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0]
FatJetResponseEtaBins = array.array("d")
FatJetResponseEtaBins.fromlist(FatJetResponseEtaBins_L)


H_FatJet_pTRef_Response_Herwig_pT1.MC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Herwig_pT2.MC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Herwig_pT3.MC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Herwig_pT1.BMC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Herwig_pT2.BMC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Herwig_pT3.BMC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Truth_pT1.MC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Truth_pT2.MC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Truth_pT3.MC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Truth_pT1.BMC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Truth_pT2.BMC.SetMarkerStyle(23)
H_FatJet_pTRef_Response_Truth_pT3.BMC.SetMarkerStyle(23)

#colorlist =[(150,45,62),(52,54,66),(151,156,156),(52,136,153),(242,235,199)]
colorlist = [2, 6, 64, 60, 92, 95, 8, 32]
DrawTool.SetColors(colorlist)

FatJet_pTRef_Ratio_Mean = PlotStat( [ [H_FatJet_pTRef_Response_pT1.MC, H_FatJet_pTRef_Response_pT2.MC, 
                                       H_FatJet_pTRef_Response_pT3.MC],
                                      [H_FatJet_pTRef_Response_pT1.BMC, H_FatJet_pTRef_Response_pT2.BMC, 
                                       H_FatJet_pTRef_Response_pT3.BMC], 
                                      [H_FatJet_pTRef_Response_pT1.Data, H_FatJet_pTRef_Response_pT2.Data, 
                                       H_FatJet_pTRef_Response_pT3.Data],
                                      [H_FatJet_pTRef_Response_pT1.BData, H_FatJet_pTRef_Response_pT2.BData, 
                                       H_FatJet_pTRef_Response_pT3.BData],
                                      [H_FatJet_pTRef_Response_Herwig_pT1.MC, H_FatJet_pTRef_Response_Herwig_pT2.MC, 
                                       H_FatJet_pTRef_Response_Herwig_pT3.MC],
                                      [H_FatJet_pTRef_Response_Herwig_pT1.BMC, H_FatJet_pTRef_Response_Herwig_pT2.BMC, 
                                       H_FatJet_pTRef_Response_Herwig_pT3.BMC],
                                      [H_FatJet_pTRef_Response_Truth_pT1.MC, H_FatJet_pTRef_Response_Truth_pT2.MC, 
                                       H_FatJet_pTRef_Response_Truth_pT3.MC],
                                      [H_FatJet_pTRef_Response_Truth_pT1.BMC, H_FatJet_pTRef_Response_Truth_pT2.BMC, 
                                       H_FatJet_pTRef_Response_Truth_pT3.BMC]],
                                    ["Pythia MC", "Pythia MC Double Btag", "Data", "Data Double Btag", "Herwig MC (W/O FC)", "Herwig MC Double B Truth (W/O FC)", "Pythia MC (W/O FC)", "Pythia MC Double B Truth (W/O FC)"  ],
                                    FatJetpTBin, "pTRef_ResponseMean", "Fat Jet p_{T}", "<r(pTRef)>", "mean" )

DrawTool.ReSetColors()


RefpT_Mean_DoubleRatios = GetRatios([FatJet_pTRef_Ratio_Mean[2], FatJet_pTRef_Ratio_Mean[3]],
                                    [FatJet_pTRef_Ratio_Mean[0], FatJet_pTRef_Ratio_Mean[1]])

for i in range(len(RefpT_Mean_DoubleRatios)):
    RefpT_Mean_DoubleRatios[i].SetMarkerStyle(20)

DrawTool.DrawHists("<R(pTRef)>", ["Trig Jet p_{T}", "<r(pTRef)>_{Data}/<r(pTRef)_{MC}>"], RefpT_Mean_DoubleRatios, ["Multi-jet inclusive", "Multi-jet b-tag"])



############## compute systematics for pT Ref
FatJet_pTRef_GenDep_Herwig_Mean = PlotStat( [ [H_FatJet_pTRef_Response_Truth_pT1.MC, H_FatJet_pTRef_Response_Truth_pT2.MC, 
                                              H_FatJet_pTRef_Response_Truth_pT3.MC],
                                             [H_FatJet_pTRef_Response_Truth_pT1.BMC, H_FatJet_pTRef_Response_Truth_pT2.BMC, 
                                              H_FatJet_pTRef_Response_Truth_pT3.BMC], 
                                             [H_FatJet_pTRef_Response_Herwig_pT1.MC, H_FatJet_pTRef_Response_Herwig_pT2.MC, 
                                              H_FatJet_pTRef_Response_Herwig_pT3.MC],
                                             [H_FatJet_pTRef_Response_Herwig_pT1.BMC, H_FatJet_pTRef_Response_Herwig_pT2.BMC, 
                                              H_FatJet_pTRef_Response_Herwig_pT3.BMC]],
                                           ["Pythia MC (W/O FC)", "Pythia MC B Truth (W/O FC)", "Herwig MC (W/O FC)", "Herwig MC B Truth (W/O FC)"],
                                           FatJetpTBin, "pTRef_ResponseMean_Herwig_Generator", "Fat Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_GenDep_Truth_Mean = PlotStat( [ [H_FatJet_pTRef_Response_pT1.MC, H_FatJet_pTRef_Response_pT2.MC, 
                                               H_FatJet_pTRef_Response_pT3.MC],
                                              [H_FatJet_pTRef_Response_pT1.BMC, H_FatJet_pTRef_Response_pT2.BMC, 
                                               H_FatJet_pTRef_Response_pT3.BMC], 
                                              [H_FatJet_pTRef_Response_Truth_pT1.MC, H_FatJet_pTRef_Response_Truth_pT2.MC, 
                                               H_FatJet_pTRef_Response_Truth_pT3.MC],
                                              [H_FatJet_pTRef_Response_Truth_pT1.BMC, H_FatJet_pTRef_Response_Truth_pT2.BMC, 
                                               H_FatJet_pTRef_Response_Truth_pT3.BMC]],
                                            ["Pythia MC", "Pythia MC Btag", "Pythia MC (W/O FC)", "Pythia MC B Truth (W/O FC)"],
                                            FatJetpTBin, "pTRef_ResponseMean_Truth_Generator", "Fat Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_TJESUp_Mean = PlotStat( [ [H_FatJet_pTRef_Response_pT1_TJESUp.MC, H_FatJet_pTRef_Response_pT2_TJESUp.MC, 
                                        H_FatJet_pTRef_Response_pT3_TJESUp.MC],
                                       [H_FatJet_pTRef_Response_pT1_TJESUp.BMC, H_FatJet_pTRef_Response_pT2_TJESUp.BMC, 
                                        H_FatJet_pTRef_Response_pT3_TJESUp.BMC], 
                                       [H_FatJet_pTRef_Response_pT1.Data, H_FatJet_pTRef_Response_pT2.Data, 
                                        H_FatJet_pTRef_Response_pT3.Data],
                                       [H_FatJet_pTRef_Response_pT1.BData, H_FatJet_pTRef_Response_pT2.BData, 
                                        H_FatJet_pTRef_Response_pT3.BData] ],
                                       ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                       FatJetpTBin, "pTRef_ResponseMean_TJESUp", "Fat Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_FJESUp_Mean = PlotStat( [ [H_FatJet_pTRef_Response_pT1_FJESUp.MC, H_FatJet_pTRef_Response_pT2_FJESUp.MC, 
                                        H_FatJet_pTRef_Response_pT3_FJESUp.MC],
                                       [H_FatJet_pTRef_Response_pT1_FJESUp.BMC, H_FatJet_pTRef_Response_pT2_FJESUp.BMC, 
                                        H_FatJet_pTRef_Response_pT3_FJESUp.BMC], 
                                       [H_FatJet_pTRef_Response_pT1.Data, H_FatJet_pTRef_Response_pT2.Data, 
                                        H_FatJet_pTRef_Response_pT3.Data],
                                       [H_FatJet_pTRef_Response_pT1.BData, H_FatJet_pTRef_Response_pT2.BData, 
                                        H_FatJet_pTRef_Response_pT3.BData] ],
                                       ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                       FatJetpTBin, "pTRef_ResponseMean_FJESUp", "Fat Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_BUp_Mean = PlotStat( [ [H_FatJet_pTRef_Response_pT1_BUp.MC, H_FatJet_pTRef_Response_pT2_BUp.MC, 
                                     H_FatJet_pTRef_Response_pT3_BUp.MC],
                                    [H_FatJet_pTRef_Response_pT1_BUp.BMC, H_FatJet_pTRef_Response_pT2_BUp.BMC, 
                                     H_FatJet_pTRef_Response_pT3_BUp.BMC], 
                                    [H_FatJet_pTRef_Response_pT1.Data, H_FatJet_pTRef_Response_pT2.Data, 
                                     H_FatJet_pTRef_Response_pT3.Data],
                                    [H_FatJet_pTRef_Response_pT1.BData, H_FatJet_pTRef_Response_pT2.BData, 
                                     H_FatJet_pTRef_Response_pT3.BData] ],
                                    ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                    FatJetpTBin, "pTRef_ResponseMean_BUp", "Fat Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_TJESDn_Mean = PlotStat( [ [H_FatJet_pTRef_Response_pT1_TJESDn.MC, H_FatJet_pTRef_Response_pT2_TJESDn.MC, 
                                        H_FatJet_pTRef_Response_pT3_TJESDn.MC],
                                       [H_FatJet_pTRef_Response_pT1_TJESDn.BMC, H_FatJet_pTRef_Response_pT2_TJESDn.BMC, 
                                        H_FatJet_pTRef_Response_pT3_TJESDn.BMC], 
                                       [H_FatJet_pTRef_Response_pT1.Data, H_FatJet_pTRef_Response_pT2.Data, 
                                        H_FatJet_pTRef_Response_pT3.Data],
                                       [H_FatJet_pTRef_Response_pT1.BData, H_FatJet_pTRef_Response_pT2.BData, 
                                        H_FatJet_pTRef_Response_pT3.BData] ],
                                       ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                       FatJetpTBin, "pTRef_ResponseMean_TJESDn", "Fat Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_FJESDn_Mean = PlotStat( [ [H_FatJet_pTRef_Response_pT1_FJESDn.MC, H_FatJet_pTRef_Response_pT2_FJESDn.MC, 
                                        H_FatJet_pTRef_Response_pT3_FJESDn.MC],
                                       [H_FatJet_pTRef_Response_pT1_FJESDn.BMC, H_FatJet_pTRef_Response_pT2_FJESDn.BMC, 
                                        H_FatJet_pTRef_Response_pT3_FJESDn.BMC], 
                                       [H_FatJet_pTRef_Response_pT1.Data, H_FatJet_pTRef_Response_pT2.Data, 
                                        H_FatJet_pTRef_Response_pT3.Data],
                                       [H_FatJet_pTRef_Response_pT1.BData, H_FatJet_pTRef_Response_pT2.BData, 
                                        H_FatJet_pTRef_Response_pT3.BData] ],
                                       ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                       FatJetpTBin, "pTRef_ResponseMean_FJESDn", "Fat Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_BDn_Mean = PlotStat( [ [H_FatJet_pTRef_Response_pT1_BDn.MC, H_FatJet_pTRef_Response_pT2_BDn.MC, 
                                     H_FatJet_pTRef_Response_pT3_BDn.MC],
                                    [H_FatJet_pTRef_Response_pT1_BDn.BMC, H_FatJet_pTRef_Response_pT2_BDn.BMC, 
                                     H_FatJet_pTRef_Response_pT3_BDn.BMC], 
                                    [H_FatJet_pTRef_Response_pT1.Data, H_FatJet_pTRef_Response_pT2.Data, 
                                     H_FatJet_pTRef_Response_pT3.Data],
                                    [H_FatJet_pTRef_Response_pT1.BData, H_FatJet_pTRef_Response_pT2.BData, 
                                     H_FatJet_pTRef_Response_pT3.BData] ],
                                    ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                    FatJetpTBin, "pTRef_ResponseMean_BDn", "Trig Jet p_{T}", "Response Mean", "gausmean" )

RefpT_Mean_GenDep_Truth_DoubleRatios = GetRatios([FatJet_pTRef_GenDep_Truth_Mean[2], FatJet_pTRef_GenDep_Truth_Mean[3]],
                                                 [FatJet_pTRef_GenDep_Truth_Mean[0], FatJet_pTRef_GenDep_Truth_Mean[1]])
for i in range(len(RefpT_Mean_GenDep_Truth_DoubleRatios)):
    RefpT_Mean_GenDep_Truth_DoubleRatios[i].SetMarkerStyle(20)
DrawTool.DrawHists("<R(pTRef)>GenDep_Truth", ["Fat Jet p_{T}", "Ratio of Mean"], RefpT_Mean_GenDep_Truth_DoubleRatios, ["Multi-jet inclusive", "Multi-jet b-tag"])


RefpT_Mean_GenDep_Herwig_DoubleRatios = GetRatios([FatJet_pTRef_GenDep_Herwig_Mean[2], FatJet_pTRef_GenDep_Herwig_Mean[3]],
                                                  [FatJet_pTRef_GenDep_Herwig_Mean[0], FatJet_pTRef_GenDep_Herwig_Mean[1]])
for i in range(len(RefpT_Mean_GenDep_Herwig_DoubleRatios)):
    RefpT_Mean_GenDep_Herwig_DoubleRatios[i].SetMarkerStyle(20)
DrawTool.DrawHists("<R(pTRef)>GenDep_Herwig", ["Fat Jet p_{T}", "Ratio of Mean"], RefpT_Mean_GenDep_Herwig_DoubleRatios, ["Multi-jet inclusive", "Multi-jet b-tag"])


RefpT_Mean_TJESUp_DoubleRatios = GetRatios([FatJet_pTRef_TJESUp_Mean[2], FatJet_pTRef_TJESUp_Mean[3]],
                                           [FatJet_pTRef_TJESUp_Mean[0], FatJet_pTRef_TJESUp_Mean[1]])
RefpT_Mean_FJESUp_DoubleRatios = GetRatios([FatJet_pTRef_FJESUp_Mean[2], FatJet_pTRef_FJESUp_Mean[3]],
                                           [FatJet_pTRef_FJESUp_Mean[0], FatJet_pTRef_FJESUp_Mean[1]])
RefpT_Mean_FJESDn_DoubleRatios = GetRatios([FatJet_pTRef_FJESDn_Mean[2], FatJet_pTRef_FJESDn_Mean[3]],
                                           [FatJet_pTRef_FJESDn_Mean[0], FatJet_pTRef_FJESDn_Mean[1]])
RefpT_Mean_FJES_DoubleRatios = GetAveRatios(RefpT_Mean_FJESUp_DoubleRatios, RefpT_Mean_FJESDn_DoubleRatios)

RefpT_Mean_BUp_DoubleRatios = GetRatios([FatJet_pTRef_BUp_Mean[2], FatJet_pTRef_BUp_Mean[3]],
                                        [FatJet_pTRef_BUp_Mean[0], FatJet_pTRef_BUp_Mean[1]])

#CalSysVariation(RefpT_Mean_DoubleRatios[0], [RefpT_Mean_GenDep_Truth_DoubleRatios[0]], True)
#CalSysVariation(RefpT_Mean_DoubleRatios[1], [RefpT_Mean_GenDep_Truth_DoubleRatios[1]], True)
#DrawTool.DrawHists("<R(pTRef)>Stat+Truth", ["Fat Jet p_{T}", "Ratio of Mean"], RefpT_Mean_DoubleRatios, ["Multi-jet inclusive", "Multi-jet b-tag"])
#
#CalSysVariation(RefpT_Mean_DoubleRatios[0], [RefpT_Mean_GenDep_Herwig_DoubleRatios[0]], True)
#CalSysVariation(RefpT_Mean_DoubleRatios[1], [RefpT_Mean_GenDep_Herwig_DoubleRatios[1]], True)
#DrawTool.DrawHists("<R(pTRef)>Stat+Truth+Gen", ["Fat Jet p_{T}", "Ratio of Mean"], RefpT_Mean_DoubleRatios, ["Multi-jet inclusive", "Multi-jet B"])

CalSysVariation(RefpT_Mean_DoubleRatios[0], [RefpT_Mean_TJESUp_DoubleRatios[0]])
CalSysVariation(RefpT_Mean_DoubleRatios[1], [RefpT_Mean_TJESUp_DoubleRatios[1]])
DrawTool.DrawHists("<R(pTRef)>Stat+Truth+Gen+TJES", ["Fat Jet p_{T}", "Ratio of Mean"], RefpT_Mean_DoubleRatios, ["Multi-jet inclusive", "Multi-jet b-tag"])

CalSysVariation(RefpT_Mean_DoubleRatios[1], [RefpT_Mean_BUp_DoubleRatios[1]])
DrawTool.DrawHists("<R(pTRef)>Stat+Truth+Gen+TJES+Btag", ["Fat Jet p_{T}", "Ratio of Mean"], RefpT_Mean_DoubleRatios, ["Multi-jet inclusive", "Multi-jet b-tag"])


JES = CalSysDiff(RefpT_Mean_DoubleRatios[0], RefpT_Mean_FJES_DoubleRatios[0] )
BJES = CalSysDiff(RefpT_Mean_DoubleRatios[1], RefpT_Mean_FJES_DoubleRatios[1] )

DrawTool.DrawHistsWithSys("JESUncert", ["Fat Jet pT", "<r(pTRef)>_{Data}/<r(pTRef)_{MC}>"], [RefpT_Mean_DoubleRatios[0]], ["Multi-jet inclusive"], [JES], [""])
DrawTool.DrawHistsWithSys("BJESUncert",   ["Fat Jet pT", "<r(pTRef)>_{Data}/<r(pTRef)_{MC}>"], [RefpT_Mean_DoubleRatios[1]], ["Multi-jet b-tag"],     [BJES], [""])


###################################trig pT  binned
FatJet_pTRef_Ratio_Mean_Trig = PlotStat( [ [H_FatJet_pTRef_Response_pT1_Trig.MC, H_FatJet_pTRef_Response_pT2_Trig.MC, 
                                            H_FatJet_pTRef_Response_pT3_Trig.MC],
                                           [H_FatJet_pTRef_Response_pT1_Trig.BMC, H_FatJet_pTRef_Response_pT2_Trig.BMC, 
                                            H_FatJet_pTRef_Response_pT3_Trig.BMC], 
                                           [H_FatJet_pTRef_Response_pT1_Trig.Data, H_FatJet_pTRef_Response_pT2_Trig.Data, 
                                            H_FatJet_pTRef_Response_pT3_Trig.Data],
                                           [H_FatJet_pTRef_Response_pT1_Trig.BData, H_FatJet_pTRef_Response_pT2_Trig.BData, 
                                            H_FatJet_pTRef_Response_pT3_Trig.BData],
                                           [H_FatJet_pTRef_Response_Herwig_pT1_Trig.MC, H_FatJet_pTRef_Response_Herwig_pT2_Trig.MC, 
                                            H_FatJet_pTRef_Response_Herwig_pT3_Trig.MC],
                                           [H_FatJet_pTRef_Response_Herwig_pT1_Trig.BMC, H_FatJet_pTRef_Response_Herwig_pT2_Trig.BMC, 
                                            H_FatJet_pTRef_Response_Herwig_pT3_Trig.BMC],
                                           [H_FatJet_pTRef_Response_Truth_pT1_Trig.MC, H_FatJet_pTRef_Response_Truth_pT2_Trig.MC, 
                                            H_FatJet_pTRef_Response_Truth_pT3_Trig.MC],
                                           [H_FatJet_pTRef_Response_Truth_pT1_Trig.BMC, H_FatJet_pTRef_Response_Truth_pT2_Trig.BMC, 
                                            H_FatJet_pTRef_Response_Truth_pT3_Trig.BMC]],
                                         ["Pythia MC", "Pythia MC Double Btag", "Data", "Data Double Btag", "Herwig MC (W/O FC)", "Herwig MC Double B Truth (W/O FC)", "Pythia MC (W/O FC)", "Pythia MC Double B Truth (W/O FC)"  ],
                                         TrigJetpTBins, "pTRef_ResponseMean_Trig", "Trig Jet p_{T}", "<r(pTRef)>", "mean" )

DrawTool.ReSetColors()

RefpT_Mean_DoubleRatios_Trig = GetRatios([FatJet_pTRef_Ratio_Mean_Trig[2], FatJet_pTRef_Ratio_Mean_Trig[3]],
                                         [FatJet_pTRef_Ratio_Mean_Trig[0], FatJet_pTRef_Ratio_Mean_Trig[1]])

for i in range(len(RefpT_Mean_DoubleRatios_Trig)):
    RefpT_Mean_DoubleRatios_Trig[i].SetMarkerStyle(20)

DrawTool.DrawHists("<R(pTRef)>_Trig", ["Trig Jet p_{T}", "<r(pTRef)}>_{Data}/<r(pTRef)_{MC}>"], RefpT_Mean_DoubleRatios_Trig, ["Multi-jet inclusive", "Multi-jet b-tag"])



############## compute systematics for pT Ref
FatJet_pTRef_GenDep_Herwig_Mean_Trig = PlotStat( [ [H_FatJet_pTRef_Response_Truth_pT1_Trig.MC, H_FatJet_pTRef_Response_Truth_pT2_Trig.MC, 
                                                    H_FatJet_pTRef_Response_Truth_pT3_Trig.MC],
                                                   [H_FatJet_pTRef_Response_Truth_pT1_Trig.BMC, H_FatJet_pTRef_Response_Truth_pT2_Trig.BMC, 
                                                    H_FatJet_pTRef_Response_Truth_pT3_Trig.BMC], 
                                                   [H_FatJet_pTRef_Response_Herwig_pT1_Trig.MC, H_FatJet_pTRef_Response_Herwig_pT2_Trig.MC, 
                                                    H_FatJet_pTRef_Response_Herwig_pT3_Trig.MC],
                                                   [H_FatJet_pTRef_Response_Herwig_pT1_Trig.BMC, H_FatJet_pTRef_Response_Herwig_pT2_Trig.BMC, 
                                                    H_FatJet_pTRef_Response_Herwig_pT3_Trig.BMC]],
                                                 ["Pythia MC (W/O FC)", "Pythia MC B Truth (W/O FC)", "Herwig MC (W/O FC)", "Herwig MC B Truth (W/O FC)"],
                                                 TrigJetpTBins, "pTRef_ResponseMean_Herwig_Generator_Trig", "Trig Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_GenDep_Truth_Mean_Trig = PlotStat( [ [H_FatJet_pTRef_Response_pT1_Trig.MC, H_FatJet_pTRef_Response_pT2_Trig.MC, 
                                               H_FatJet_pTRef_Response_pT3_Trig.MC],
                                              [H_FatJet_pTRef_Response_pT1_Trig.BMC, H_FatJet_pTRef_Response_pT2_Trig.BMC, 
                                               H_FatJet_pTRef_Response_pT3_Trig.BMC], 
                                              [H_FatJet_pTRef_Response_Truth_pT1_Trig.MC, H_FatJet_pTRef_Response_Truth_pT2_Trig.MC, 
                                               H_FatJet_pTRef_Response_Truth_pT3_Trig.MC],
                                              [H_FatJet_pTRef_Response_Truth_pT1_Trig.BMC, H_FatJet_pTRef_Response_Truth_pT2_Trig.BMC, 
                                               H_FatJet_pTRef_Response_Truth_pT3_Trig.BMC]],
                                            ["Pythia MC", "Pythia MC Btag", "Pythia MC (W/O FC)", "Pythia MC B Truth (W/O FC)"],
                                            TrigJetpTBins, "pTRef_ResponseMean_Truth_Generator_Trig", "Trig Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_TJESUp_Mean_Trig = PlotStat( [ [H_FatJet_pTRef_Response_pT1_TJESUp_Trig.MC, H_FatJet_pTRef_Response_pT2_TJESUp_Trig.MC, 
                                             H_FatJet_pTRef_Response_pT3_TJESUp_Trig.MC],
                                            [H_FatJet_pTRef_Response_pT1_TJESUp_Trig.BMC, H_FatJet_pTRef_Response_pT2_TJESUp_Trig.BMC, 
                                             H_FatJet_pTRef_Response_pT3_TJESUp_Trig.BMC], 
                                            [H_FatJet_pTRef_Response_pT1_Trig.Data, H_FatJet_pTRef_Response_pT2_Trig.Data, 
                                             H_FatJet_pTRef_Response_pT3_Trig.Data],
                                            [H_FatJet_pTRef_Response_pT1_Trig.BData, H_FatJet_pTRef_Response_pT2_Trig.BData, 
                                             H_FatJet_pTRef_Response_pT3_Trig.BData] ],
                                          ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                          TrigJetpTBins, "pTRef_ResponseMean_TJESUp_Trig", "Trig Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_FJESUp_Mean_Trig = PlotStat( [ [H_FatJet_pTRef_Response_pT1_FJESUp_Trig.MC, H_FatJet_pTRef_Response_pT2_FJESUp_Trig.MC, 
                                             H_FatJet_pTRef_Response_pT3_FJESUp_Trig.MC],
                                            [H_FatJet_pTRef_Response_pT1_FJESUp_Trig.BMC, H_FatJet_pTRef_Response_pT2_FJESUp_Trig.BMC, 
                                             H_FatJet_pTRef_Response_pT3_FJESUp_Trig.BMC], 
                                            [H_FatJet_pTRef_Response_pT1_Trig.Data, H_FatJet_pTRef_Response_pT2_Trig.Data, 
                                             H_FatJet_pTRef_Response_pT3_Trig.Data],
                                            [H_FatJet_pTRef_Response_pT1_Trig.BData, H_FatJet_pTRef_Response_pT2_Trig.BData, 
                                             H_FatJet_pTRef_Response_pT3_Trig.BData] ],
                                          ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                          TrigJetpTBins, "pTRef_ResponseMean_FJESUp_Trig", "Trig Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_BUp_Mean_Trig = PlotStat( [ [H_FatJet_pTRef_Response_pT1_BUp_Trig.MC, H_FatJet_pTRef_Response_pT2_BUp_Trig.MC, 
                                          H_FatJet_pTRef_Response_pT3_BUp_Trig.MC],
                                         [H_FatJet_pTRef_Response_pT1_BUp_Trig.BMC, H_FatJet_pTRef_Response_pT2_BUp_Trig.BMC, 
                                          H_FatJet_pTRef_Response_pT3_BUp_Trig.BMC], 
                                         [H_FatJet_pTRef_Response_pT1_Trig.Data, H_FatJet_pTRef_Response_pT2_Trig.Data, 
                                          H_FatJet_pTRef_Response_pT3_Trig.Data],
                                         [H_FatJet_pTRef_Response_pT1_Trig.BData, H_FatJet_pTRef_Response_pT2_Trig.BData, 
                                          H_FatJet_pTRef_Response_pT3_Trig.BData] ],
                                       ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                       TrigJetpTBins, "pTRef_ResponseMean_BUp_Trig", "Trig Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_TJESDn_Mean_Trig = PlotStat( [ [H_FatJet_pTRef_Response_pT1_TJESDn_Trig.MC, H_FatJet_pTRef_Response_pT2_TJESDn_Trig.MC, 
                                             H_FatJet_pTRef_Response_pT3_TJESDn_Trig.MC],
                                            [H_FatJet_pTRef_Response_pT1_TJESDn_Trig.BMC, H_FatJet_pTRef_Response_pT2_TJESDn_Trig.BMC, 
                                             H_FatJet_pTRef_Response_pT3_TJESDn_Trig.BMC], 
                                            [H_FatJet_pTRef_Response_pT1_Trig.Data, H_FatJet_pTRef_Response_pT2_Trig.Data, 
                                             H_FatJet_pTRef_Response_pT3_Trig.Data],
                                            [H_FatJet_pTRef_Response_pT1_Trig.BData, H_FatJet_pTRef_Response_pT2_Trig.BData, 
                                             H_FatJet_pTRef_Response_pT3_Trig.BData] ],
                                          ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                          TrigJetpTBins, "pTRef_ResponseMean_TJESDn_Trig", "Trig Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_FJESDn_Mean_Trig = PlotStat( [ [H_FatJet_pTRef_Response_pT1_FJESDn_Trig.MC, H_FatJet_pTRef_Response_pT2_FJESDn_Trig.MC, 
                                             H_FatJet_pTRef_Response_pT3_FJESDn_Trig.MC],
                                            [H_FatJet_pTRef_Response_pT1_FJESDn_Trig.BMC, H_FatJet_pTRef_Response_pT2_FJESDn_Trig.BMC, 
                                             H_FatJet_pTRef_Response_pT3_FJESDn_Trig.BMC], 
                                            [H_FatJet_pTRef_Response_pT1_Trig.Data, H_FatJet_pTRef_Response_pT2_Trig.Data, 
                                             H_FatJet_pTRef_Response_pT3_Trig.Data],
                                            [H_FatJet_pTRef_Response_pT1_Trig.BData, H_FatJet_pTRef_Response_pT2_Trig.BData, 
                                             H_FatJet_pTRef_Response_pT3_Trig.BData] ],
                                          ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                          TrigJetpTBins, "pTRef_ResponseMean_FJESDn_Trig", "Trig Jet p_{T}", "Response Mean", "gausmean" )

FatJet_pTRef_BDn_Mean_Trig = PlotStat( [ [H_FatJet_pTRef_Response_pT1_BDn_Trig.MC, H_FatJet_pTRef_Response_pT2_BDn_Trig.MC, 
                                          H_FatJet_pTRef_Response_pT3_BDn_Trig.MC],
                                         [H_FatJet_pTRef_Response_pT1_BDn_Trig.BMC, H_FatJet_pTRef_Response_pT2_BDn_Trig.BMC, 
                                          H_FatJet_pTRef_Response_pT3_BDn_Trig.BMC], 
                                         [H_FatJet_pTRef_Response_pT1_Trig.Data, H_FatJet_pTRef_Response_pT2_Trig.Data, 
                                          H_FatJet_pTRef_Response_pT3_Trig.Data],
                                         [H_FatJet_pTRef_Response_pT1_Trig.BData, H_FatJet_pTRef_Response_pT2_Trig.BData, 
                                          H_FatJet_pTRef_Response_pT3_Trig.BData] ],
                                       ["Pythia MC", "Pythia MC Btag", "Data", "Data Btag"],
                                       TrigJetpTBins, "pTRef_ResponseMean_BDn_Trig", "Trig Jet p_{T}", "Response Mean", "gausmean" )

RefpT_Mean_GenDep_Truth_DoubleRatios_Trig = GetRatios([FatJet_pTRef_GenDep_Truth_Mean_Trig[2], FatJet_pTRef_GenDep_Truth_Mean_Trig[3]],
                                                 [FatJet_pTRef_GenDep_Truth_Mean_Trig[0], FatJet_pTRef_GenDep_Truth_Mean_Trig[1]])
#for i in range(len(RefpT_Mean_GenDep_Truth_DoubleRatios)):
#    RefpT_Mean_GenDep_Truth_DoubleRatios[i].SetMarkerStyle(20)
DrawTool.DrawHists("<R(pTRef)>GenDep_Truth_Trig", ["Fat Jet p_{T}", "Ratio of Mean"], RefpT_Mean_GenDep_Truth_DoubleRatios_Trig, ["Multi-jet inclusive", "Multi-jet b-tag"])


RefpT_Mean_GenDep_Herwig_DoubleRatios_Trig = GetRatios([FatJet_pTRef_GenDep_Herwig_Mean_Trig[2], FatJet_pTRef_GenDep_Herwig_Mean_Trig[3]],
                                                       [FatJet_pTRef_GenDep_Herwig_Mean_Trig[0], FatJet_pTRef_GenDep_Herwig_Mean_Trig[1]])
for i in range(len(RefpT_Mean_GenDep_Herwig_DoubleRatios_Trig)):
    RefpT_Mean_GenDep_Herwig_DoubleRatios_Trig[i].SetMarkerStyle(20)
DrawTool.DrawHists("<R(pTRef)>GenDep_Herwig_Trig", ["Trig Jet p_{T}", "Ratio of Mean"], RefpT_Mean_GenDep_Herwig_DoubleRatios_Trig, ["Multi-jet inclusive", "Multi-jet b-tag"])


RefpT_Mean_TJESUp_DoubleRatios_Trig = GetRatios([FatJet_pTRef_TJESUp_Mean_Trig[2], FatJet_pTRef_TJESUp_Mean_Trig[3]],
                                                [FatJet_pTRef_TJESUp_Mean_Trig[0], FatJet_pTRef_TJESUp_Mean_Trig[1]])
RefpT_Mean_FJESUp_DoubleRatios_Trig = GetRatios([FatJet_pTRef_FJESUp_Mean_Trig[2], FatJet_pTRef_FJESUp_Mean_Trig[3]],
                                                [FatJet_pTRef_FJESUp_Mean_Trig[0], FatJet_pTRef_FJESUp_Mean_Trig[1]])
RefpT_Mean_FJESDn_DoubleRatios_Trig = GetRatios([FatJet_pTRef_FJESDn_Mean_Trig[2], FatJet_pTRef_FJESDn_Mean_Trig[3]],
                                                [FatJet_pTRef_FJESDn_Mean_Trig[0], FatJet_pTRef_FJESDn_Mean_Trig[1]])
RefpT_Mean_FJES_DoubleRatios_Trig = GetAveRatios(RefpT_Mean_FJESUp_DoubleRatios_Trig, RefpT_Mean_FJESDn_DoubleRatios_Trig)

RefpT_Mean_BUp_DoubleRatios_Trig = GetRatios([FatJet_pTRef_BUp_Mean_Trig[2], FatJet_pTRef_BUp_Mean_Trig[3]],
                                             [FatJet_pTRef_BUp_Mean_Trig[0], FatJet_pTRef_BUp_Mean_Trig[1]])


CalSysVariation(RefpT_Mean_DoubleRatios_Trig[0], [RefpT_Mean_GenDep_Truth_DoubleRatios[0]], True)
CalSysVariation(RefpT_Mean_DoubleRatios[1], [RefpT_Mean_GenDep_Truth_DoubleRatios[1]], True)
DrawTool.DrawHists("<R(pTRef)>Stat+Truth", ["Fat Jet p_{T}", "Ratio of Mean"], RefpT_Mean_DoubleRatios, ["Multi-jet inclusive", "Multi-jet b-tag"])

CalSysVariation(RefpT_Mean_DoubleRatios_Trig[0], [RefpT_Mean_GenDep_Herwig_DoubleRatios_Trig[0]], True)
CalSysVariation(RefpT_Mean_DoubleRatios_Trig[1], [RefpT_Mean_GenDep_Herwig_DoubleRatios_Trig[1]], True)
DrawTool.DrawHists("<R(pTRef)>Stat+Truth+Gen_Trig", ["Trig Jet p_{T}", "Ratio of Mean"], RefpT_Mean_DoubleRatios_Trig, ["Multi-jet inclusive", "Multi-jet b-tag"])

CalSysVariation(RefpT_Mean_DoubleRatios_Trig[0], [RefpT_Mean_TJESUp_DoubleRatios_Trig[0]])
CalSysVariation(RefpT_Mean_DoubleRatios_Trig[1], [RefpT_Mean_TJESUp_DoubleRatios_Trig[1]])
DrawTool.DrawHists("<R(pTRef)>Stat+Truth+Gen+TJES_Trig", ["Trig Jet p_{T}", "Ratio of Mean"], RefpT_Mean_DoubleRatios_Trig, ["Multi-jet inclusive", "Multi-jet b-tag"])

CalSysVariation(RefpT_Mean_DoubleRatios_Trig[1], [RefpT_Mean_BUp_DoubleRatios_Trig[1]])
DrawTool.DrawHists("<R(pTRef)>Stat+Truth+Gen+TJES+Btag_Trig", ["Trig Jet p_{T}", "Ratio of Mean"], RefpT_Mean_DoubleRatios_Trig, ["Multi-jet inclusive", "Multi-jet b-tag"])


JES_Trig = CalSysDiff(RefpT_Mean_DoubleRatios_Trig[0], RefpT_Mean_FJESUp_DoubleRatios_Trig[0] )
BJES_Trig = CalSysDiff(RefpT_Mean_DoubleRatios_Trig[1], RefpT_Mean_FJESUp_DoubleRatios_Trig[1] )

DrawTool.DrawHistsWithSys("JESUncert_Trig", ["Trig Jet pT", "<r(pTRef)>_{Data}/<r(pTRef)_{MC}>"], [RefpT_Mean_DoubleRatios_Trig[0]], ["Multi-jet inclusive"], [JES_Trig], [""])
DrawTool.DrawHistsWithSys("BJESUncert_Trig",   ["Trig Jet pT", "<r(pTRef)>_{Data}/<r(pTRef)_{MC}>"], [RefpT_Mean_DoubleRatios_Trig[1]], ["Multi-jet b-tag"],     [BJES_Trig], [""])


if (isHerwig):
    sys.exit(0)

####### Response Summary Plots ##################

DrawTool = HistoTool("", "Preliminary", "0", "8", True, True)
DrawTool.AddTexts(0.23, 0.76, 1, 0.048, "|#eta|<1.2")

BJMS_Eta1_S = PlotStat( [ [H_FatJet_Mass_Response_NStack_pT1_Eta1.MC, H_FatJet_Mass_Response_NStack_pT2_Eta1.MC, 
                           H_FatJet_Mass_Response_NStack_pT3_Eta1.MC],
                          [H_FatJet_Mass_Response_NStack_pT1_Eta1.BMC, H_FatJet_Mass_Response_NStack_pT2_Eta1.BMC, 
                           H_FatJet_Mass_Response_NStack_pT3_Eta1.BMC],
                          [H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta1.MC, H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta1.MC,
                           H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta1.MC],
                          [H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta1.BMC, H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta1.BMC,
                           H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta1.BMC]],
                        ["Multi-jet inclusive", "Multi-jet Btag", "HHbbbb", "HHbbbb Btag"],
                        FatJetMOPTBin, "ResponseGausMean_Eta1", "Fat Jet Mass / p_{T}", "Response Mean", "gausmean" )
DrawTool.ClearTexts()
BJMS_Eta1 = CopyHist(BJMS_Eta1_S[1])
BJMS_Eta1.Divide(BJMS_Eta1_S[0])
BJMS_Eta1.SetMarkerStyle(20)
DrawTool.DrawHists("BJMS_Eta1", ["Fat Jet Mass / p_{T}", "BJMS"], [BJMS_Eta1], ["|#eta|<1.2"])


DrawTool.AddTexts(0.23, 0.76, 1, 0.048, "|#eta|>1.2")
BJMS_Eta2_S = PlotStat( [ [H_FatJet_Mass_Response_NStack_pT1_Eta2.MC, H_FatJet_Mass_Response_NStack_pT2_Eta2.MC, 
                           H_FatJet_Mass_Response_NStack_pT3_Eta2.MC],
                          [H_FatJet_Mass_Response_NStack_pT1_Eta2.BMC, H_FatJet_Mass_Response_NStack_pT2_Eta2.BMC, 
                           H_FatJet_Mass_Response_NStack_pT3_Eta2.BMC],
                          [H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta2.MC, H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta2.MC,
                           H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta2.MC],
                          [H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta2.BMC, H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta2.BMC,
                           H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta2.BMC]],
                        ["Multi-jet inclusive", "Multi-jet Btag", "HHbbbb", "HHbbbb Btag"],
                        FatJetMOPTBin, "ResponseGausMean_Eta2", "Fat Jet Mass / p_{T}", "Response Mean", "gausmean" )
DrawTool.ClearTexts()
BJMS_Eta2 = CopyHist(BJMS_Eta2_S[1])
BJMS_Eta2.Divide(BJMS_Eta2_S[0])
BJMS_Eta2.SetMarkerStyle(20)
DrawTool.DrawHists("BJMS_Eta2", ["Fat Jet Mass / p_{T}", "BJMS"], [BJMS_Eta2], ["|#eta|>1.2"])


PlotStat( [ [H_FatJet_Mass_Response_NStack_Eta1.MC, H_FatJet_Mass_Response_NStack_Eta2.MC, 
             H_FatJet_Mass_Response_NStack_Eta3.MC, H_FatJet_Mass_Response_NStack_Eta4.MC, 
             H_FatJet_Mass_Response_NStack_Eta5.MC, H_FatJet_Mass_Response_NStack_Eta6.MC, 
             H_FatJet_Mass_Response_NStack_Eta7.MC, H_FatJet_Mass_Response_NStack_Eta8.MC],
            [H_FatJet_Mass_Response_NStack_Eta1.BMC, H_FatJet_Mass_Response_NStack_Eta2.BMC, 
             H_FatJet_Mass_Response_NStack_Eta3.BMC, H_FatJet_Mass_Response_NStack_Eta4.BMC, 
             H_FatJet_Mass_Response_NStack_Eta5.BMC, H_FatJet_Mass_Response_NStack_Eta6.BMC, 
             H_FatJet_Mass_Response_NStack_Eta7.BMC, H_FatJet_Mass_Response_NStack_Eta8.BMC] ],
          ["Multi-jet inclusive", "Multi-jet b-tag"],
          FatJetResponseEtaBins, "ResponseGausMean_Eta", "Eta", "Response Mean", "gausmean")


PlotStat( [ [H_FatJet_Mass_Response_NStack_pT1_Eta1.MC, H_FatJet_Mass_Response_NStack_pT2_Eta1.MC, 
             H_FatJet_Mass_Response_NStack_pT3_Eta1.MC],
            [H_FatJet_Mass_Response_NStack_pT1_Eta1.BMC, H_FatJet_Mass_Response_NStack_pT2_Eta1.BMC, 
             H_FatJet_Mass_Response_NStack_pT3_Eta1.BMC],
            [H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta1.MC, H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta1.MC,
             H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta1.MC],
            [H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta1.BMC, H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta1.BMC,
             H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta1.BMC]],
          ["Multi-jet inclusive", "Multi-jet b-tag", "HHbbbb", "HHbbbb b-tag"],
           FatJetMOPTBin, "ResponseRMS_Eta1", "Mass / p_{T}", "Response RMS", "rms" )

PlotStat( [ [H_FatJet_Mass_Response_NStack_pT1_Eta2.MC, H_FatJet_Mass_Response_NStack_pT2_Eta2.MC, 
             H_FatJet_Mass_Response_NStack_pT3_Eta2.MC],
            [H_FatJet_Mass_Response_NStack_pT1_Eta2.BMC, H_FatJet_Mass_Response_NStack_pT2_Eta2.BMC, 
             H_FatJet_Mass_Response_NStack_pT3_Eta2.BMC],
            [H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta2.MC, H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta2.MC,
             H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta2.MC],
            [H_FatJet_Mass_Response_HHbbbb_NStack_pT1_Eta2.BMC, H_FatJet_Mass_Response_HHbbbb_NStack_pT2_Eta2.BMC,
             H_FatJet_Mass_Response_HHbbbb_NStack_pT3_Eta2.BMC]],
          ["Multi-jet inclusive", "Multi-jet b-tag", "HHbbbb", "HHbbbb b-tag"],
           FatJetMOPTBin, "ResponseRMS_Eta2", "Mass / p_{T}", "Response RMS", "rms" )


######## nTrack Plots
nTrackMean        = [H_FatJet_nTrack_pT1.MC.GetMean(1), 
                     H_FatJet_nTrack_pT2.MC.GetMean(1), 
                     H_FatJet_nTrack_pT3.MC.GetMean(1)]
nTrackMean_Error  = [H_FatJet_nTrack_pT1.MC.GetMeanError(1), 
                     H_FatJet_nTrack_pT2.MC.GetMeanError(1), 
                     H_FatJet_nTrack_pT3.MC.GetMeanError(1)]
BnTrackMean       = [H_FatJet_nTrack_pT1.BMC.GetMean(1), 
                     H_FatJet_nTrack_pT2.BMC.GetMean(1), 
                     H_FatJet_nTrack_pT3.BMC.GetMean(1)]
BnTrackMean_Error = [H_FatJet_nTrack_pT1.BMC.GetMeanError(1), 
                     H_FatJet_nTrack_pT2.BMC.GetMeanError(1), 
                     H_FatJet_nTrack_pT3.BMC.GetMeanError(1)]

Hist_nTrackMean  = TH1D("nTrackMean", "nTrackMean", 3, FatJetMOPTBin)
Hist_BnTrackMean = TH1D("BnTrackMean", "BnTrackMean", 3, FatJetMOPTBin)
SetHist(Hist_nTrackMean, nTrackMean, nTrackMean_Error)
SetHist(Hist_BnTrackMean, BnTrackMean, BnTrackMean_Error)
DrawTool.DrawHists("nTrackMean", ["Mass / p_{T}", "nTrack Mean"], [Hist_nTrackMean, Hist_BnTrackMean], ["Inclusive", "b-tag"])



############# trackmass ratio
TrackMass_Ratios_Eta1 = PlotStat( [ [H_FatJet_TrackMass_Ratio_pT1_Eta1.MC, 
                                     H_FatJet_TrackMass_Ratio_pT2_Eta1.MC, 
                                     H_FatJet_TrackMass_Ratio_pT3_Eta1.MC],
                                    [H_FatJet_TrackMass_Ratio_pT1_Eta1.BMC, 
                                     H_FatJet_TrackMass_Ratio_pT2_Eta1.BMC, 
                                     H_FatJet_TrackMass_Ratio_pT3_Eta1.BMC],
                                    [H_FatJet_TrackMass_Ratio_pT1_Eta1.Data, 
                                     H_FatJet_TrackMass_Ratio_pT2_Eta1.Data, 
                                     H_FatJet_TrackMass_Ratio_pT3_Eta1.Data],
                                    [H_FatJet_TrackMass_Ratio_pT1_Eta1.BData, 
                                     H_FatJet_TrackMass_Ratio_pT2_Eta1.BData, 
                                     H_FatJet_TrackMass_Ratio_pT3_Eta1.BData] ],
           
                                  ["MC", "MC b-tag", "Data", "Data b-tag"],
                                  FatJetMOPTBin, "TrackMassRatio_Eta1", "Mass / p_{T}", "#LTr_{trk}^{Mass}#GT", "gausmean" )

TrackMass_Ratios_Gaus_Eta1 = PlotStat( [ [H_FatJet_TrackMass_Ratio_pT1_Eta1.MC, 
                                          H_FatJet_TrackMass_Ratio_pT2_Eta1.MC, 
                                          H_FatJet_TrackMass_Ratio_pT3_Eta1.MC],
                                         [H_FatJet_TrackMass_Ratio_pT1_Eta1.BMC, 
                                          H_FatJet_TrackMass_Ratio_pT2_Eta1.BMC, 
                                          H_FatJet_TrackMass_Ratio_pT3_Eta1.BMC],
                                         [H_FatJet_TrackMass_Ratio_pT1_Eta1.Data, 
                                          H_FatJet_TrackMass_Ratio_pT2_Eta1.Data, 
                                          H_FatJet_TrackMass_Ratio_pT3_Eta1.Data],
                                         [H_FatJet_TrackMass_Ratio_pT1_Eta1.BData, 
                                          H_FatJet_TrackMass_Ratio_pT2_Eta1.BData, 
                                          H_FatJet_TrackMass_Ratio_pT3_Eta1.BData] ],
           
                                       ["MC", "MC b-tag", "Data", "Data b-tag"],
                                       FatJetMOPTBin, "TrackMassRatio_Gaus_Eta1", "Mass / p_{T}", "<r(trk)^{Mass}>", "gausmean" )


TrackMass_Ratios_Eta2 = PlotStat( [ [H_FatJet_TrackMass_Ratio_pT1_Eta2.MC, 
                                     H_FatJet_TrackMass_Ratio_pT2_Eta2.MC, 
                                     H_FatJet_TrackMass_Ratio_pT3_Eta2.MC],
                                    [H_FatJet_TrackMass_Ratio_pT1_Eta2.BMC, 
                                     H_FatJet_TrackMass_Ratio_pT2_Eta2.BMC, 
                                     H_FatJet_TrackMass_Ratio_pT3_Eta2.BMC],
                                    [H_FatJet_TrackMass_Ratio_pT1_Eta2.Data, 
                                     H_FatJet_TrackMass_Ratio_pT2_Eta2.Data, 
                                     H_FatJet_TrackMass_Ratio_pT3_Eta2.Data],
                                    [H_FatJet_TrackMass_Ratio_pT1_Eta2.BData, 
                                     H_FatJet_TrackMass_Ratio_pT2_Eta2.BData, 
                                     H_FatJet_TrackMass_Ratio_pT3_Eta2.BData] ],
           
                                  ["MC", "MC b-tag", "Data", "Data b-tag"],
                                  FatJetMOPTBin, "TrackMassRatio_Eta2", "Mass / p_{T}", "#LTr_{trk}^{Mass}#GT", "gausmean" )

TrackMass_Ratios_Gaus_Eta2 = PlotStat( [ [H_FatJet_TrackMass_Ratio_pT1_Eta2.MC, 
                                          H_FatJet_TrackMass_Ratio_pT2_Eta2.MC, 
                                          H_FatJet_TrackMass_Ratio_pT3_Eta2.MC],
                                         [H_FatJet_TrackMass_Ratio_pT1_Eta2.BMC, 
                                          H_FatJet_TrackMass_Ratio_pT2_Eta2.BMC, 
                                          H_FatJet_TrackMass_Ratio_pT3_Eta2.BMC],
                                         [H_FatJet_TrackMass_Ratio_pT1_Eta2.Data, 
                                          H_FatJet_TrackMass_Ratio_pT2_Eta2.Data, 
                                          H_FatJet_TrackMass_Ratio_pT3_Eta2.Data],
                                         [H_FatJet_TrackMass_Ratio_pT1_Eta2.BData, 
                                          H_FatJet_TrackMass_Ratio_pT2_Eta2.BData, 
                                          H_FatJet_TrackMass_Ratio_pT3_Eta2.BData] ],
           
                                       ["MC", "MC b-tag", "Data", "Data b-tag"],
                                       FatJetMOPTBin, "TrackMassRatio_Gaus_Eta2", "Mass / p_{T}", "<r(trk)^{Mass}>", "gausmean" )


DoubleRatios = GetRatios([TrackMass_Ratios_Eta1[2], TrackMass_Ratios_Eta1[3],
                          TrackMass_Ratios_Eta2[2], TrackMass_Ratios_Eta2[3]],
                         [TrackMass_Ratios_Eta1[0], TrackMass_Ratios_Eta1[1],
                          TrackMass_Ratios_Eta2[0], TrackMass_Ratios_Eta2[1]])
for i in range(len(DoubleRatios)):
    DoubleRatios[i].SetMarkerStyle(20)


DrawTool.DrawHists("<R(trk)>", ["Mass / p_{T}", "R(trk) Mean"], DoubleRatios, ["|#eta|<1.2, inclusive", "|#eta|<1.2, btag", "|#eta|>1.2, inclusive", "|#eta|>1.2, btag"])

DoubleRatios_gaus = GetRatios([TrackMass_Ratios_Gaus_Eta1[2], TrackMass_Ratios_Gaus_Eta1[3],
                               TrackMass_Ratios_Gaus_Eta2[2], TrackMass_Ratios_Gaus_Eta2[3]],
                              [TrackMass_Ratios_Gaus_Eta1[0], TrackMass_Ratios_Gaus_Eta1[1],
                               TrackMass_Ratios_Gaus_Eta2[0], TrackMass_Ratios_Gaus_Eta2[1]])

DoubleRatios_gaus[0].SetMarkerStyle(26)
DoubleRatios_gaus[1].SetMarkerStyle(24)
DoubleRatios_gaus[2].SetMarkerStyle(25)
DoubleRatios_gaus[3].SetMarkerStyle(32)

DrawTool.DrawHists("<R(trk)>_gaus", ["Mass / p_{T}", "<r(trk)^{Mass}>_{Data}/<r(trk)^{Mass}>_{MC}"], DoubleRatios_gaus, ["|#eta|<1.2, inclusive", "|#eta|<1.2, btag", "|#eta|>1.2, inclusive", "|#eta|>1.2, btag"])




#####################JMS Constraints #####################

TrackMass_Ratios_BUp_Eta1 = PlotStat( [ [H_FatJet_TrackMass_Ratio_BUp_pT1_Eta1.MC, 
                                         H_FatJet_TrackMass_Ratio_BUp_pT2_Eta1.MC, 
                                         H_FatJet_TrackMass_Ratio_BUp_pT3_Eta1.MC],
                                        [H_FatJet_TrackMass_Ratio_BUp_pT1_Eta1.BMC, 
                                         H_FatJet_TrackMass_Ratio_BUp_pT2_Eta1.BMC, 
                                         H_FatJet_TrackMass_Ratio_BUp_pT3_Eta1.BMC],
                                        [H_FatJet_TrackMass_Ratio_pT1_Eta1.Data, 
                                         H_FatJet_TrackMass_Ratio_pT2_Eta1.Data, 
                                         H_FatJet_TrackMass_Ratio_pT3_Eta1.Data],
                                        [H_FatJet_TrackMass_Ratio_pT1_Eta1.BData, 
                                         H_FatJet_TrackMass_Ratio_pT2_Eta1.BData, 
                                         H_FatJet_TrackMass_Ratio_pT3_Eta1.BData] ],
                                      
                                      ["MC", "MC b-tag", "Data", "Data b-tag"],
                                      FatJetMOPTBin, "TrackMassRatio_Gaus_BUp_Eta1", "Mass / p_{T}", "#LTr_{trk}^{Mass}#GT", "gausmean" )

TrackMass_Ratios_Herwig_Eta1 = PlotStat( [ [H_FatJet_TrackMass_Ratio_Herwig_pT1_Eta1.MC, 
                                            H_FatJet_TrackMass_Ratio_Herwig_pT2_Eta1.MC, 
                                            H_FatJet_TrackMass_Ratio_Herwig_pT3_Eta1.MC],
                                           [H_FatJet_TrackMass_Ratio_Herwig_pT1_Eta1.BMC, 
                                            H_FatJet_TrackMass_Ratio_Herwig_pT2_Eta1.BMC, 
                                            H_FatJet_TrackMass_Ratio_Herwig_pT3_Eta1.BMC],
                                           [H_FatJet_TrackMass_Ratio_pT1_Eta1.Data, 
                                            H_FatJet_TrackMass_Ratio_pT2_Eta1.Data, 
                                            H_FatJet_TrackMass_Ratio_pT3_Eta1.Data],
                                           [H_FatJet_TrackMass_Ratio_pT1_Eta1.BData, 
                                            H_FatJet_TrackMass_Ratio_pT2_Eta1.BData, 
                                            H_FatJet_TrackMass_Ratio_pT3_Eta1.BData] ],
                                         
                                         ["MC", "MC b-tag", "Data", "Data b-tag"],
                                         FatJetMOPTBin, "TrackMassRatio_Gaus_Herwig_Eta1", "Mass / p_{T}", "r(trk) Mean", "gausmean" )


TrackMass_Ratios_BUp_Eta2 = PlotStat( [ [H_FatJet_TrackMass_Ratio_BUp_pT1_Eta2.MC, 
                                         H_FatJet_TrackMass_Ratio_BUp_pT2_Eta2.MC, 
                                         H_FatJet_TrackMass_Ratio_BUp_pT3_Eta2.MC],
                                        [H_FatJet_TrackMass_Ratio_BUp_pT1_Eta2.BMC, 
                                         H_FatJet_TrackMass_Ratio_BUp_pT2_Eta2.BMC, 
                                         H_FatJet_TrackMass_Ratio_BUp_pT3_Eta2.BMC],
                                        [H_FatJet_TrackMass_Ratio_pT1_Eta2.Data, 
                                         H_FatJet_TrackMass_Ratio_pT2_Eta2.Data, 
                                         H_FatJet_TrackMass_Ratio_pT3_Eta2.Data],
                                        [H_FatJet_TrackMass_Ratio_pT1_Eta2.BData, 
                                         H_FatJet_TrackMass_Ratio_pT2_Eta2.BData, 
                                         H_FatJet_TrackMass_Ratio_pT3_Eta2.BData] ],
                                      
                                      ["MC", "MC b-tag", "Data", "Data b-tag"],
                                      FatJetMOPTBin, "TrackMassRatio_Gaus_BUp_Eta2", "Mass / p_{T}", "r(trk) Mean", "gausmean" )


TrackMass_Ratios_Herwig_Eta2 = PlotStat( [ [H_FatJet_TrackMass_Ratio_Herwig_pT1_Eta2.MC, 
                                            H_FatJet_TrackMass_Ratio_Herwig_pT2_Eta2.MC, 
                                            H_FatJet_TrackMass_Ratio_Herwig_pT3_Eta2.MC],
                                           [H_FatJet_TrackMass_Ratio_Herwig_pT1_Eta2.BMC, 
                                            H_FatJet_TrackMass_Ratio_Herwig_pT2_Eta2.BMC, 
                                            H_FatJet_TrackMass_Ratio_Herwig_pT3_Eta2.BMC],
                                           [H_FatJet_TrackMass_Ratio_pT1_Eta2.Data, 
                                            H_FatJet_TrackMass_Ratio_pT2_Eta2.Data, 
                                            H_FatJet_TrackMass_Ratio_pT3_Eta2.Data],
                                           [H_FatJet_TrackMass_Ratio_pT1_Eta2.BData, 
                                            H_FatJet_TrackMass_Ratio_pT2_Eta2.BData, 
                                            H_FatJet_TrackMass_Ratio_pT3_Eta2.BData] ],
                                         
                                         ["MC", "MC b-tag", "Data", "Data b-tag"],
                                         FatJetMOPTBin, "TrackMassRatio_Gaus_Herwig_Eta2", "Mass / p_{T}", "r(trk) Mean", "gausmean" )


TrackMass_Ratios_JMSUp_Eta1 = PlotStat( [ [H_FatJet_TrackMass_Ratio_JMSUp_pT1_Eta1.MC, 
                                           H_FatJet_TrackMass_Ratio_JMSUp_pT2_Eta1.MC, 
                                           H_FatJet_TrackMass_Ratio_JMSUp_pT3_Eta1.MC],
                                          [H_FatJet_TrackMass_Ratio_JMSUp_pT1_Eta1.BMC, 
                                           H_FatJet_TrackMass_Ratio_JMSUp_pT2_Eta1.BMC, 
                                           H_FatJet_TrackMass_Ratio_JMSUp_pT3_Eta1.BMC],
                                          [H_FatJet_TrackMass_Ratio_pT1_Eta1.Data, 
                                           H_FatJet_TrackMass_Ratio_pT2_Eta1.Data, 
                                           H_FatJet_TrackMass_Ratio_pT3_Eta1.Data],
                                          [H_FatJet_TrackMass_Ratio_pT1_Eta1.BData, 
                                           H_FatJet_TrackMass_Ratio_pT2_Eta1.BData, 
                                           H_FatJet_TrackMass_Ratio_pT3_Eta1.BData] ],
                                        
                                        ["MC", "MC b-tag", "Data", "Data b-tag"],
                                        FatJetMOPTBin, "TrackMassRatio_Gaus_JMSUp_Eta1", "Mass / p_{T}", "r(trk) Mean", "gausmean" )

TrackMass_Ratios_JMSUp_Eta2 = PlotStat( [ [H_FatJet_TrackMass_Ratio_JMSUp_pT1_Eta2.MC, 
                                           H_FatJet_TrackMass_Ratio_JMSUp_pT2_Eta2.MC, 
                                           H_FatJet_TrackMass_Ratio_JMSUp_pT3_Eta2.MC],
                                          [H_FatJet_TrackMass_Ratio_JMSUp_pT1_Eta2.BMC, 
                                           H_FatJet_TrackMass_Ratio_JMSUp_pT2_Eta2.BMC, 
                                           H_FatJet_TrackMass_Ratio_JMSUp_pT3_Eta2.BMC],
                                          [H_FatJet_TrackMass_Ratio_pT1_Eta2.Data, 
                                           H_FatJet_TrackMass_Ratio_pT2_Eta2.Data, 
                                           H_FatJet_TrackMass_Ratio_pT3_Eta2.Data],
                                          [H_FatJet_TrackMass_Ratio_pT1_Eta2.BData, 
                                           H_FatJet_TrackMass_Ratio_pT2_Eta2.BData, 
                                           H_FatJet_TrackMass_Ratio_pT3_Eta2.BData] ],
                                        
                                        ["MC", "MC b-tag", "Data", "Data b-tag"],
                                        FatJetMOPTBin, "TrackMassRatio_Gaus_JMSUp_Eta2", "Mass / p_{T}", "r(trk) Mean", "gausmean" )


TrackMass_JMSUp_DoubleRatios_Eta1 = GetRatios([TrackMass_Ratios_JMSUp_Eta1[2], TrackMass_Ratios_JMSUp_Eta1[3]],
                                              [TrackMass_Ratios_JMSUp_Eta1[0], TrackMass_Ratios_JMSUp_Eta1[1]])
TrackMass_JMSUp_DoubleRatios_Eta2 = GetRatios([TrackMass_Ratios_JMSUp_Eta2[2], TrackMass_Ratios_JMSUp_Eta2[3]],
                                              [TrackMass_Ratios_JMSUp_Eta2[0], TrackMass_Ratios_JMSUp_Eta2[1]])
TrackMass_BUp_DoubleRatios_Eta1 = GetRatios([TrackMass_Ratios_BUp_Eta1[2], TrackMass_Ratios_BUp_Eta1[3]],
                                            [TrackMass_Ratios_BUp_Eta1[0], TrackMass_Ratios_BUp_Eta1[1]])
TrackMass_BUp_DoubleRatios_Eta2 = GetRatios([TrackMass_Ratios_BUp_Eta2[2], TrackMass_Ratios_BUp_Eta2[3]],
                                            [TrackMass_Ratios_BUp_Eta2[0], TrackMass_Ratios_BUp_Eta2[1]])

TrackMass_Herwig_DoubleRatios_Eta1 = GetRatios([TrackMass_Ratios_Herwig_Eta1[2], TrackMass_Ratios_Herwig_Eta1[3]],
                                            [TrackMass_Ratios_Herwig_Eta1[0], TrackMass_Ratios_Herwig_Eta1[1]])
TrackMass_Herwig_DoubleRatios_Eta2 = GetRatios([TrackMass_Ratios_Herwig_Eta2[2], TrackMass_Ratios_Herwig_Eta2[3]],
                                            [TrackMass_Ratios_Herwig_Eta2[0], TrackMass_Ratios_Herwig_Eta2[1]])

TrackMass_JMSUp_Sys_Eta1 = CalSysDiff(DoubleRatios_gaus[0], TrackMass_JMSUp_DoubleRatios_Eta1[0] )
TrackMass_JMSUp_Sys_Eta1_btag = CalSysDiff(DoubleRatios_gaus[1], TrackMass_JMSUp_DoubleRatios_Eta1[1] )
TrackMass_JMSUp_Sys_Eta2 = CalSysDiff(DoubleRatios_gaus[2], TrackMass_JMSUp_DoubleRatios_Eta2[0] )
TrackMass_JMSUp_Sys_Eta2_btag = CalSysDiff(DoubleRatios_gaus[3], TrackMass_JMSUp_DoubleRatios_Eta2[1] )


CalSysVariation(DoubleRatios_gaus[0], [TrackMass_Herwig_DoubleRatios_Eta1[0]], True)
#CalSysVariation(DoubleRatios_gaus[1], [TrackMass_Herwig_DoubleRatios_Eta1[1]], True)
CalSysVariation(DoubleRatios_gaus[2], [TrackMass_Herwig_DoubleRatios_Eta2[0]], True)
#CalSysVariation(DoubleRatios_gaus[3], [TrackMass_Herwig_DoubleRatios_Eta2[1]], True)

CalSysVariation(DoubleRatios_gaus[0], [TrackMass_BUp_DoubleRatios_Eta1[0]])
CalSysVariation(DoubleRatios_gaus[1], [TrackMass_BUp_DoubleRatios_Eta1[1]])
CalSysVariation(DoubleRatios_gaus[2], [TrackMass_BUp_DoubleRatios_Eta2[0]])
CalSysVariation(DoubleRatios_gaus[3], [TrackMass_BUp_DoubleRatios_Eta2[1]])

DoubleRatios_gaus[0].SetMarkerStyle(20)
DoubleRatios_gaus[1].SetMarkerStyle(20)
DoubleRatios_gaus[2].SetMarkerStyle(20)
DoubleRatios_gaus[3].SetMarkerStyle(20)


DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet, |#eta|<1.2")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.DrawHistsWithSys("Rtrk_inclusive_eta1", ["Mass / p_{T}", "#LTr_{trk}^{Mass}#GT_{Data}/#LTr_{trk}^{Mass}#GT_{MC}"], [DoubleRatios_gaus[0]], ["Multi-jet inclusive"], [TrackMass_JMSUp_Sys_Eta1], ["Existing JMS Uncertainty"] )
DrawTool.DrawHistsWithSys("Rtrk_btag_eta1",      ["Mass / p_{T}", "#LTr_{trk}^{Mass}#GT_{Data}/#LTr_{trk}^{Mass}#GT_{MC}"], [DoubleRatios_gaus[1]], ["Multi-jet b-tag"],      [TrackMass_JMSUp_Sys_Eta1_btag], ["Existing JMS Uncertainty"] )
DrawTool.ClearTexts()

DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet, |#eta|>1.2")
DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool.DrawHistsWithSys("Rtrk_inclusive_eta2", ["Mass / p_{T}", "#LTr_{trk}^{Mass}#GT_{Data}/#LTr_{trk}^{Mass}#GT_{MC}"], [DoubleRatios_gaus[2]], ["Multi-jet inclusive"], [TrackMass_JMSUp_Sys_Eta2], ["Existing JMS Uncertainty"] )
DrawTool.DrawHistsWithSys("Rtrk_btag_eta2",      ["Mass / p_{T}", "#LTr_{trk}^{Mass}#GT_{Data}/#LTr_{trk}^{Mass}#GT_{MC}"], [DoubleRatios_gaus[3]], ["Multi-jet b-tag"],      [TrackMass_JMSUp_Sys_Eta2_btag], ["Existing JMS Uncertainty"] )
DrawTool.ClearTexts()


DrawTool_DoubleRatio.SysLabels = ["Existing JMS Uncertainty"]
DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")

TrackMass_Ratios_Eta1[0].SetMarkerStyle(26)
TrackMass_Ratios_Eta1[1].SetMarkerStyle(26)
TrackMass_Ratios_Eta1[0].SetLineStyle(2)
TrackMass_Ratios_Eta1[1].SetLineStyle(2)
TrackMass_Ratios_Eta1[2].SetMarkerStyle(20)
TrackMass_Ratios_Eta1[3].SetMarkerStyle(20)

TrackMass_Ratios_Eta1[0].GetYaxis().SetTitleSize( TrackMass_Ratios_Eta1[0].GetYaxis().GetTitleSize()*1.1)
TrackMass_Ratios_Eta1[1].GetYaxis().SetTitleSize( TrackMass_Ratios_Eta1[0].GetYaxis().GetTitleSize()*1.1)
TrackMass_Ratios_Eta1[2].GetYaxis().SetTitleSize( TrackMass_Ratios_Eta1[0].GetYaxis().GetTitleSize()*1.1)
TrackMass_Ratios_Eta1[3].GetYaxis().SetTitleSize( TrackMass_Ratios_Eta1[0].GetYaxis().GetTitleSize()*1.1)


DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet before b-tagging, |#eta|<1.2")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool_DoubleRatio.PreCalRatio = DoubleRatios_gaus[0]
DrawTool_DoubleRatio.DrawHists("rtrk_mass_eta1", ["p_{T} [GeV]", "#LTr_{trk}^{Mass}#GT"],  [TrackMass_Ratios_Eta1[0], TrackMass_Ratios_Eta1[0], TrackMass_Ratios_Eta1[2]], ["Pythia 8 MC", "Pythia 8 MC", "Data"], [], [], [TrackMass_JMSUp_Sys_Eta1])
DrawTool_DoubleRatio.ClearTexts()

DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging, |#eta|<1.2")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool_DoubleRatio.PreCalRatio = DoubleRatios_gaus[1]
DrawTool_DoubleRatio.DrawHists("rtrk_mass_eta1_btag", ["p_{T} [GeV]", "#LTr_{trk}^{Mass}#GT"],  [TrackMass_Ratios_Eta1[1], TrackMass_Ratios_Eta1[1], TrackMass_Ratios_Eta1[3]], ["Pythia 8 MC", "Pythia 8 MC", "Data"], [], [], [TrackMass_JMSUp_Sys_Eta1_btag])
DrawTool_DoubleRatio.ClearTexts()

TrackMass_Ratios_Eta2[0].SetMarkerStyle(26)
TrackMass_Ratios_Eta2[1].SetMarkerStyle(26)
TrackMass_Ratios_Eta2[0].SetLineStyle(2)
TrackMass_Ratios_Eta2[1].SetLineStyle(2)
TrackMass_Ratios_Eta2[2].SetMarkerStyle(20)
TrackMass_Ratios_Eta2[3].SetMarkerStyle(20)

TrackMass_Ratios_Eta2[0].GetYaxis().SetTitleSize( TrackMass_Ratios_Eta2[0].GetYaxis().GetTitleSize()*1.1)
TrackMass_Ratios_Eta2[1].GetYaxis().SetTitleSize( TrackMass_Ratios_Eta2[0].GetYaxis().GetTitleSize()*1.1)
TrackMass_Ratios_Eta2[2].GetYaxis().SetTitleSize( TrackMass_Ratios_Eta2[0].GetYaxis().GetTitleSize()*1.1)
TrackMass_Ratios_Eta2[3].GetYaxis().SetTitleSize( TrackMass_Ratios_Eta2[0].GetYaxis().GetTitleSize()*1.1)

DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet before b-tagging, |#eta|>1.2")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool_DoubleRatio.PreCalRatio = DoubleRatios_gaus[2]
DrawTool_DoubleRatio.DrawHists("rtrk_mass_eta2", ["p_{T} [GeV]", "#LTr_{trk}^{Mass}#GT"],  [TrackMass_Ratios_Eta2[0], TrackMass_Ratios_Eta2[0], TrackMass_Ratios_Eta2[2]], ["Pythia 8 MC", "Pythia 8 MC", "Data"], [], [], [TrackMass_JMSUp_Sys_Eta2])
DrawTool_DoubleRatio.ClearTexts()

DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging, |#eta|>1.2")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool_DoubleRatio.PreCalRatio = DoubleRatios_gaus[3]
DrawTool_DoubleRatio.DrawHists("rtrk_mass_eta2_btag", ["p_{T} [GeV]", "#LTr_{trk}^{Mass}#GT"],  [TrackMass_Ratios_Eta2[1], TrackMass_Ratios_Eta2[1], TrackMass_Ratios_Eta2[3]], ["Pythia 8 MC", "Pythia 8 MC", "Data"], [], [], [TrackMass_JMSUp_Sys_Eta2_btag])
DrawTool_DoubleRatio.ClearTexts()


############# Track pT Ratio

TrackpT_Ratios_Gaus = PlotStat( [ [H_FatJet_TrackpT_Ratio_pT1.MC, 
                                   H_FatJet_TrackpT_Ratio_pT2.MC, 
                                   H_FatJet_TrackpT_Ratio_pT3.MC],
                                  [H_FatJet_TrackpT_Ratio_pT1.BMC, 
                                   H_FatJet_TrackpT_Ratio_pT2.BMC, 
                                   H_FatJet_TrackpT_Ratio_pT3.BMC],
                                  [H_FatJet_TrackpT_Ratio_pT1.Data, 
                                   H_FatJet_TrackpT_Ratio_pT2.Data, 
                                   H_FatJet_TrackpT_Ratio_pT3.Data],
                                  [H_FatJet_TrackpT_Ratio_pT1.BData, 
                                   H_FatJet_TrackpT_Ratio_pT2.BData, 
                                   H_FatJet_TrackpT_Ratio_pT3.BData] ],
                                
                                ["MC", "MC b-tag", "Data", "Data b-tag"],
                                FatJetpTBin, "TrackpTRatio_Gaus", "p_{T} [GeV]", "#LTr_{trk}^{p_{T}}#GT", "gausmean" )

DoubleRatios_gaus = GetRatios([TrackpT_Ratios_Gaus[2], TrackpT_Ratios_Gaus[3]],
                              [TrackpT_Ratios_Gaus[0], TrackpT_Ratios_Gaus[1]])

DoubleRatios_gaus[0].SetMarkerStyle(26)
DoubleRatios_gaus[1].SetMarkerStyle(24)

#DrawTool.DrawHists("<r_{trk}>_pT_gaus", ["Fat Jet p_{T}", "<r_{trk}^{p_{T}}>_{Data}/<r_{trk}^{p_{T}}>_{MC}"], DoubleRatios_gaus, ["Multi-jet inclusive", "Multi-jet b-tag"])



TrackpT_Ratios_GenDep_Herwig_Mean = PlotStat( [ [H_FatJet_Truth_TrackpT_Ratio_pT1.MC, H_FatJet_Truth_TrackpT_Ratio_pT2.MC, 
                                                 H_FatJet_Truth_TrackpT_Ratio_pT3.MC ],
                                                [H_FatJet_Truth_TrackpT_Ratio_pT1.BMC, H_FatJet_Truth_TrackpT_Ratio_pT2.BMC,  
                                                 H_FatJet_Truth_TrackpT_Ratio_pT3.BMC], 
                                                [H_FatJet_Herwig_TrackpT_Ratio_pT1.MC, H_FatJet_Herwig_TrackpT_Ratio_pT2.MC, 
                                                 H_FatJet_Herwig_TrackpT_Ratio_pT3.MC ],
                                                [H_FatJet_Herwig_TrackpT_Ratio_pT1.BMC, H_FatJet_Herwig_TrackpT_Ratio_pT2.BMC,  
                                                 H_FatJet_Herwig_TrackpT_Ratio_pT3.BMC]],
                                              
                                              ["Pythia MC (W/O FC)", "Pythia MC B Truth (W/O FC)", "Herwig MC (W/O FC)", "Herwig MC B Truth (W/O FC)"],
                                              FatJetpTBin, "TrackpTRatio_Gaus_Herwig", "Fat Jet p_{T}", "r_{trk} Mean", "gausmean" )

TrackpT_Ratios_GenDep_Truth_Mean = PlotStat( [ [H_FatJet_TrackpT_Ratio_pT1.MC, H_FatJet_TrackpT_Ratio_pT2.MC, 
                                                H_FatJet_TrackpT_Ratio_pT3.MC ],
                                               [H_FatJet_TrackpT_Ratio_pT1.BMC, H_FatJet_TrackpT_Ratio_pT2.BMC,  
                                                H_FatJet_TrackpT_Ratio_pT3.BMC], 
                                               [H_FatJet_Truth_TrackpT_Ratio_pT1.MC, H_FatJet_Truth_TrackpT_Ratio_pT2.MC, 
                                                H_FatJet_Truth_TrackpT_Ratio_pT3.MC ],
                                               [H_FatJet_Truth_TrackpT_Ratio_pT1.BMC, H_FatJet_Truth_TrackpT_Ratio_pT2.BMC,  
                                                H_FatJet_Truth_TrackpT_Ratio_pT3.BMC]],
                                              
                                             ["Pythia MC", "Pythia MC B", "Pythia MC (W/O FC)", "Pythia MC B Truth (W/O FC)"],
                                             FatJetpTBin, "TrackpTRatio_Gaus_Truth", "Fat Jet p_{T}", "r_{trk} Mean", "gausmean" )


TrackpT_Ratios_BUp = PlotStat( [ [H_FatJet_TrackpT_Ratio_BUp_pT1.MC, 
                                  H_FatJet_TrackpT_Ratio_BUp_pT2.MC, 
                                  H_FatJet_TrackpT_Ratio_BUp_pT3.MC],
                                 [H_FatJet_TrackpT_Ratio_BUp_pT1.BMC, 
                                  H_FatJet_TrackpT_Ratio_BUp_pT2.BMC, 
                                  H_FatJet_TrackpT_Ratio_BUp_pT3.BMC],
                                 [H_FatJet_TrackpT_Ratio_pT1.Data, 
                                  H_FatJet_TrackpT_Ratio_pT2.Data, 
                                  H_FatJet_TrackpT_Ratio_pT3.Data],
                                 [H_FatJet_TrackpT_Ratio_pT1.BData, 
                                  H_FatJet_TrackpT_Ratio_pT2.BData, 
                                  H_FatJet_TrackpT_Ratio_pT3.BData] ],
                               
                               ["MC", "MC b-tag", "Data", "Data b-tag"],
                               FatJetpTBin, "TrackpTRatio_Gaus_BUp", "Fat Jet p{T}", "r_{trk} Mean", "gausmean" )

TrackpT_Ratios_FJESUp = PlotStat( [ [H_FatJet_TrackpT_Ratio_FJESDn_pT1.MC, 
                                     H_FatJet_TrackpT_Ratio_FJESDn_pT2.MC, 
                                     H_FatJet_TrackpT_Ratio_FJESDn_pT3.MC],
                                    [H_FatJet_TrackpT_Ratio_FJESDn_pT1.BMC, 
                                     H_FatJet_TrackpT_Ratio_FJESDn_pT2.BMC, 
                                     H_FatJet_TrackpT_Ratio_FJESDn_pT3.BMC],
                                    [H_FatJet_TrackpT_Ratio_pT1.Data, 
                                     H_FatJet_TrackpT_Ratio_pT2.Data, 
                                     H_FatJet_TrackpT_Ratio_pT3.Data],
                                    [H_FatJet_TrackpT_Ratio_pT1.BData, 
                                     H_FatJet_TrackpT_Ratio_pT2.BData, 
                                     H_FatJet_TrackpT_Ratio_pT3.BData] ],
                                  
                                  ["MC", "MC b-tag", "Data", "Data b-tag"],
                                  FatJetpTBin, "TrackpTRatio_Gaus_FJESUp", "Fat Jet p_{T}", "r_{trk} Mean", "gausmean" )


TrackpT_GenDep_Herwig_Mean_DoubleRatios = GetRatios([TrackpT_Ratios_GenDep_Herwig_Mean[2], TrackpT_Ratios_GenDep_Herwig_Mean[3]],
                                                    [TrackpT_Ratios_GenDep_Herwig_Mean[0], TrackpT_Ratios_GenDep_Herwig_Mean[1]])

TrackpT_GenDep_Truth_Mean_DoubleRatios = GetRatios([TrackpT_Ratios_GenDep_Truth_Mean[2], TrackpT_Ratios_GenDep_Truth_Mean[3]],
                                                   [TrackpT_Ratios_GenDep_Truth_Mean[0], TrackpT_Ratios_GenDep_Truth_Mean[1]])

TrackpT_FJESUp_DoubleRatios = GetRatios([TrackpT_Ratios_FJESUp[2], TrackpT_Ratios_FJESUp[3]],
                                        [TrackpT_Ratios_FJESUp[0], TrackpT_Ratios_FJESUp[1]])
TrackpT_BUp_DoubleRatios = GetRatios([TrackpT_Ratios_BUp[2], TrackpT_Ratios_BUp[3]],
                                     [TrackpT_Ratios_BUp[0], TrackpT_Ratios_BUp[1]])

TrackpT_FJESUp_Sys = CalSysDiff(DoubleRatios_gaus[0], TrackpT_FJESUp_DoubleRatios[0] )
TrackpT_FJESUp_Sys_btag = CalSysDiff(DoubleRatios_gaus[1], TrackpT_FJESUp_DoubleRatios[1] )


#CalSysVariation(DoubleRatios_gaus[0], [TrackpT_GenDep_Truth_Mean_DoubleRatios[0]])

CalSysVariation(DoubleRatios_gaus[1], [TrackpT_GenDep_Truth_Mean_DoubleRatios[1]])
CalSysVariation(DoubleRatios_gaus[1], [TrackpT_GenDep_Truth_Mean_DoubleRatios[1]])
CalSysVariation(DoubleRatios_gaus[1], [TrackpT_GenDep_Truth_Mean_DoubleRatios[1]])
CalSysVariation(DoubleRatios_gaus[1], [TrackpT_GenDep_Truth_Mean_DoubleRatios[1]])

DrawTool.DrawHists("<R(TrackpT)>Stat+Truth+Gen", ["Fat Jet p_{T}", "Ratio of Mean"], DoubleRatios_gaus, ["Multi-jet inclusive", "Multi-jet b-tag"])
DrawTool.DrawHists("<R(TrackpT)>Stat+Truth+GenTot", ["Fat Jet p_{T}", "Ratio of Mean"], DoubleRatios_gaus, ["Multi-jet inclusive", "Multi-jet b-tag"])

CalSysVariation(DoubleRatios_gaus[1], [TrackpT_BUp_DoubleRatios[1]])

CalSysVariation(DoubleRatios_gaus[0], [TrackpT_GenDep_Truth_Mean_DoubleRatios[1]])
CalSysVariation(DoubleRatios_gaus[0], [TrackpT_GenDep_Truth_Mean_DoubleRatios[0]])

DoubleRatios_gaus[0].SetMarkerStyle(20)
DoubleRatios_gaus[1].SetMarkerStyle(20)

#DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet")
#DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
#DrawTool.DrawHistsWithSys("Rtrk_pT", ["p_{T} [GeV]", "#LTr_{trk}^{p_{T}}#GT_{Data}/#LTr_{trk}^{p_{T}}#GT_{MC}"], [DoubleRatios_gaus[0]], ["Multi-jet inclusive"], [TrackpT_FJESUp_Sys], ["Existing JES Uncertainty"] )
#DrawTool.DrawHistsWithSys("Rtrk_pT_btag",    ["p_{T} [GeV]", "#LTr_{trk}^{p_{T}}#GT_{Data}/#LTr_{trk}^{p_{T}}#GT_{MC}"], [DoubleRatios_gaus[1]], ["Multi-jet b-tag"],      [TrackpT_FJESUp_Sys_btag], ["Existing JES Uncertainty"] )
#DrawTool.ClearTexts()


DrawTool_DoubleRatio = HistoTool()
DrawTool_DoubleRatio.SetFillStyle(1001)
DrawTool_DoubleRatio.SetCompData(True)
DrawTool_DoubleRatio.colorlist = DrawTool_DoubleRatio.boringcolor
DrawTool_DoubleRatio.studytype = "Preliminary"

DrawTool_DoubleRatio.SysLabels = ["Existing JES Uncertainty"]
DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")

TrackpT_Ratios_Gaus[0].SetMarkerStyle(26)
TrackpT_Ratios_Gaus[1].SetMarkerStyle(26)
TrackpT_Ratios_Gaus[0].SetLineStyle(2)
TrackpT_Ratios_Gaus[1].SetLineStyle(2)
TrackpT_Ratios_Gaus[2].SetMarkerStyle(20)
TrackpT_Ratios_Gaus[3].SetMarkerStyle(20)

TrackpT_Ratios_Gaus[0].GetYaxis().SetTitleSize( TrackpT_Ratios_Gaus[0].GetYaxis().GetTitleSize()*1.1)
TrackpT_Ratios_Gaus[1].GetYaxis().SetTitleSize( TrackpT_Ratios_Gaus[0].GetYaxis().GetTitleSize()*1.1)
TrackpT_Ratios_Gaus[2].GetYaxis().SetTitleSize( TrackpT_Ratios_Gaus[0].GetYaxis().GetTitleSize()*1.1)
TrackpT_Ratios_Gaus[3].GetYaxis().SetTitleSize( TrackpT_Ratios_Gaus[0].GetYaxis().GetTitleSize()*1.1)

DrawTool_DoubleRatio.PreCalRatio = DoubleRatios_gaus[0]

DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet before b-tagging")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool_DoubleRatio.DrawHists("rtrk_pT", ["p_{T} [GeV]", "#LTr_{trk}^{p_{T}}#GT"],  [TrackpT_Ratios_Gaus[0], TrackpT_Ratios_Gaus[0], TrackpT_Ratios_Gaus[2]], ["Pythia 8 MC", "Pythia 8 MC", "Data"], [], [], [TrackpT_FJESUp_Sys])

DrawTool_DoubleRatio.ClearTexts()

DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")

DrawTool_DoubleRatio.PreCalRatio = DoubleRatios_gaus[1]
DrawTool_DoubleRatio.DrawHists("rtrk_pT_btag", ["p_{T} [GeV]", "#LTr_{trk}^{p_{T}}#GT"],  [TrackpT_Ratios_Gaus[1], TrackpT_Ratios_Gaus[1], TrackpT_Ratios_Gaus[3]], ["Pythia 8 MC", "Pythia 8 MC", "Data"], [], [], [TrackpT_FJESUp_Sys_btag])

DrawTool_DoubleRatio.ClearTexts()


############# Track D2 Ratio
TrackD2_Ratios_Gaus = PlotStat( [ [H_FatJet_TrackD2_Ratio_pT1.MC, 
                                   H_FatJet_TrackD2_Ratio_pT2.MC, 
                                   H_FatJet_TrackD2_Ratio_pT3.MC],
                                  [H_FatJet_TrackD2_Ratio_pT1.BMC, 
                                   H_FatJet_TrackD2_Ratio_pT2.BMC, 
                                   H_FatJet_TrackD2_Ratio_pT3.BMC],
                                  [H_FatJet_TrackD2_Ratio_pT1.Data, 
                                   H_FatJet_TrackD2_Ratio_pT2.Data, 
                                   H_FatJet_TrackD2_Ratio_pT3.Data],
                                  [H_FatJet_TrackD2_Ratio_pT1.BData, 
                                   H_FatJet_TrackD2_Ratio_pT2.BData, 
                                   H_FatJet_TrackD2_Ratio_pT3.BData] ],
                                
                                ["MC", "MC b-tag", "Data", "Data b-tag"],
                                FatJetpTBin, "TrackD2Ratio_Gaus", "p_{T} [GeV]", "#LTr_{trk}^{D_{2}}#GT",  "gausmean" )
 

DoubleRatios_gaus = GetRatios([TrackD2_Ratios_Gaus[2], TrackD2_Ratios_Gaus[3]],
                              [TrackD2_Ratios_Gaus[0], TrackD2_Ratios_Gaus[1]])

DoubleRatios_gaus[0].SetMarkerStyle(26)
DoubleRatios_gaus[1].SetMarkerStyle(24)

DrawTool.DrawHists("<r_{trk}>_D2_gaus", ["Fat Jet p_{T} [GeV]", "<r_{trk}^{p_{T}}>_{Data}/<r_{trk}^{p_{T}}>_{MC}"], DoubleRatios_gaus, ["Multi-jet inclusive", "Multi-jet b-tag"])


TrackD2_Ratios_BUp = PlotStat( [ [H_FatJet_TrackD2_Ratio_BUp_pT1.MC, 
                                  H_FatJet_TrackD2_Ratio_BUp_pT2.MC, 
                                  H_FatJet_TrackD2_Ratio_BUp_pT3.MC],
                                 [H_FatJet_TrackD2_Ratio_BUp_pT1.BMC, 
                                  H_FatJet_TrackD2_Ratio_BUp_pT2.BMC, 
                                  H_FatJet_TrackD2_Ratio_BUp_pT3.BMC],
                                 [H_FatJet_TrackD2_Ratio_pT1.Data, 
                                  H_FatJet_TrackD2_Ratio_pT2.Data, 
                                  H_FatJet_TrackD2_Ratio_pT3.Data],
                                 [H_FatJet_TrackD2_Ratio_pT1.BData, 
                                  H_FatJet_TrackD2_Ratio_pT2.BData, 
                                  H_FatJet_TrackD2_Ratio_pT3.BData] ],
                               
                               ["MC", "MC b-tag", "Data", "Data b-tag"],
                               FatJetpTBin, "TrackD2Ratio_Gaus_BUp", "Fat Jet p{T}", "r_{trk} Mean", "gausmean" )

TrackD2_Ratios_D2Up = PlotStat( [ [H_FatJet_TrackD2_Ratio_D2Up_pT1.MC, 
                                  H_FatJet_TrackD2_Ratio_D2Up_pT2.MC, 
                                  H_FatJet_TrackD2_Ratio_D2Up_pT3.MC],
                                 [H_FatJet_TrackD2_Ratio_D2Up_pT1.BMC, 
                                  H_FatJet_TrackD2_Ratio_D2Up_pT2.BMC, 
                                  H_FatJet_TrackD2_Ratio_D2Up_pT3.BMC],
                                 [H_FatJet_TrackD2_Ratio_pT1.Data, 
                                  H_FatJet_TrackD2_Ratio_pT2.Data, 
                                  H_FatJet_TrackD2_Ratio_pT3.Data],
                                 [H_FatJet_TrackD2_Ratio_pT1.BData, 
                                  H_FatJet_TrackD2_Ratio_pT2.BData, 
                                  H_FatJet_TrackD2_Ratio_pT3.BData] ],
                               
                               ["MC", "MC b-tag", "Data", "Data b-tag"],
                               FatJetpTBin, "TrackD2Ratio_Gaus_D2Up", "Fat Jet p{T}", "r_{trk} Mean", "gausmean" )

TrackD2_Ratios_D2Dn = PlotStat( [ [H_FatJet_TrackD2_Ratio_D2Dn_pT1.MC, 
                                  H_FatJet_TrackD2_Ratio_D2Dn_pT2.MC, 
                                  H_FatJet_TrackD2_Ratio_D2Dn_pT3.MC],
                                 [H_FatJet_TrackD2_Ratio_D2Dn_pT1.BMC, 
                                  H_FatJet_TrackD2_Ratio_D2Dn_pT2.BMC, 
                                  H_FatJet_TrackD2_Ratio_D2Dn_pT3.BMC],
                                 [H_FatJet_TrackD2_Ratio_pT1.Data, 
                                  H_FatJet_TrackD2_Ratio_pT2.Data, 
                                  H_FatJet_TrackD2_Ratio_pT3.Data],
                                 [H_FatJet_TrackD2_Ratio_pT1.BData, 
                                  H_FatJet_TrackD2_Ratio_pT2.BData, 
                                  H_FatJet_TrackD2_Ratio_pT3.BData] ],
                               
                               ["MC", "MC b-tag", "Data", "Data b-tag"],
                               FatJetpTBin, "TrackD2Ratio_Gaus_D2Dn", "Fat Jet p{T}", "r_{trk} Mean", "gausmean" )

TrackD2_Ratios_FJESUp = PlotStat( [ [H_FatJet_TrackD2_Ratio_FJESDn_pT1.MC, 
                                     H_FatJet_TrackD2_Ratio_FJESDn_pT2.MC, 
                                     H_FatJet_TrackD2_Ratio_FJESDn_pT3.MC],
                                    [H_FatJet_TrackD2_Ratio_FJESDn_pT1.BMC, 
                                     H_FatJet_TrackD2_Ratio_FJESDn_pT2.BMC, 
                                     H_FatJet_TrackD2_Ratio_FJESDn_pT3.BMC],
                                    [H_FatJet_TrackD2_Ratio_pT1.Data, 
                                     H_FatJet_TrackD2_Ratio_pT2.Data, 
                                     H_FatJet_TrackD2_Ratio_pT3.Data],
                                    [H_FatJet_TrackD2_Ratio_pT1.BData, 
                                     H_FatJet_TrackD2_Ratio_pT2.BData, 
                                     H_FatJet_TrackD2_Ratio_pT3.BData] ],
                                  
                                  ["MC", "MC b-tag", "Data", "Data b-tag"],
                                  FatJetpTBin, "TrackD2Ratio_Gaus_FJESUp", "Fat Jet p_{T} [GeV]", "r_{trk} Mean", "gausmean" )


#TrackD2_GenDep_Herwig_Mean_DoubleRatios = GetRatios([TrackD2_Ratios_GenDep_Herwig_Mean[2], TrackD2_Ratios_GenDep_Herwig_Mean[3]],
#                                                    [TrackD2_Ratios_GenDep_Herwig_Mean[0], TrackD2_Ratios_GenDep_Herwig_Mean[1]])

#TrackD2_GenDep_Truth_Mean_DoubleRatios = GetRatios([TrackD2_Ratios_GenDep_Truth_Mean[2], TrackD2_Ratios_GenDep_Truth_Mean[3]],
#                                                   [TrackD2_Ratios_GenDep_Truth_Mean[0], TrackD2_Ratios_GenDep_Truth_Mean[1]])

TrackD2_FJESUp_DoubleRatios = GetRatios([TrackD2_Ratios_FJESUp[2], TrackD2_Ratios_FJESUp[3]],
                                        [TrackD2_Ratios_FJESUp[0], TrackD2_Ratios_FJESUp[1]])
TrackD2_BUp_DoubleRatios = GetRatios([TrackD2_Ratios_BUp[2], TrackD2_Ratios_BUp[3]],
                                     [TrackD2_Ratios_BUp[0], TrackD2_Ratios_BUp[1]])
TrackD2_D2Up_DoubleRatios = GetRatios([TrackD2_Ratios_D2Up[2], TrackD2_Ratios_D2Up[3]],
                                      [TrackD2_Ratios_D2Up[0], TrackD2_Ratios_D2Up[1]])
TrackD2_D2Dn_DoubleRatios = GetRatios([TrackD2_Ratios_D2Dn[2], TrackD2_Ratios_D2Dn[3]],
                                      [TrackD2_Ratios_D2Dn[0], TrackD2_Ratios_D2Dn[1]])

TrackD2_D2Up_Sys = CalSysDiff(DoubleRatios_gaus[0], TrackD2_D2Up_DoubleRatios[0] )
TrackD2_D2Dn_Sys_btag = CalSysDiff(DoubleRatios_gaus[1], TrackD2_D2Dn_DoubleRatios[1] )


#CalSysVariation(DoubleRatios_gaus[0], [TrackD2_GenDep_Truth_Mean_DoubleRatios[0]])
#CalSysVariation(DoubleRatios_gaus[1], [TrackD2_GenDep_Truth_Mean_DoubleRatios[1]])
#CalSysVariation(DoubleRatios_gaus[1], [TrackD2_GenDep_Truth_Mean_DoubleRatios[1]])
#CalSysVariation(DoubleRatios_gaus[1], [TrackD2_GenDep_Truth_Mean_DoubleRatios[1]])
#CalSysVariation(DoubleRatios_gaus[1], [TrackD2_GenDep_Truth_Mean_DoubleRatios[1]])

CalSysVariation(DoubleRatios_gaus[0], [TrackD2_FJESUp_DoubleRatios[0]])
CalSysVariation(DoubleRatios_gaus[1], [TrackD2_FJESUp_DoubleRatios[1]])
CalSysVariation(DoubleRatios_gaus[1], [TrackD2_BUp_DoubleRatios[1]])
#CalSysVariation(DoubleRatios_gaus[1], [TrackD2_FJESUp_DoubleRatios[1]])

DoubleRatios_gaus[0].SetMarkerStyle(20)
DoubleRatios_gaus[1].SetMarkerStyle(20)

#DrawTool.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet")
#DrawTool.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
#DrawTool.DrawHistsWithSys("Rtrk_D2", ["p_{T} [GeV]", "#LTr_{trk}^{D_{2}^{(#beta=1)}}#GT_{Data}/#LTr_{trk}^{D_{2}^{(#beta=1)}}#GT_{MC}"], [DoubleRatios_gaus[0]], ["Multi-jet inclusive"], [TrackD2_D2Up_Sys], ["Existing D_{2}^{(#beta=1)} Uncertainty"] )
#
#DrawTool.DrawHistsWithSys("Rtrk_D2_btag",    ["p_{T} [GeV]", "#LTr_{trk}^{D_{2}^{(#beta=1)}}#GT_{Data}/#LTr_{trk}^{D_{2}^{(#beta=1)}}#GT_{MC}"], [DoubleRatios_gaus[1]], ["Multi-jet b-tag"],      [TrackD2_D2Up_Sys_btag], ["Existing D_{2}^{(#beta=1)} Uncertainty"] )
#
##TrackD2_D2Dn_Sys_btag.SetBinError(2, 0.08)
#DrawTool.DrawHistsWithSys("Rtrk_D2_btag_dn",    ["p_{T} [GeV]", "#LTr_{trk}^{D_{2}^{(#beta=1)}}#GT_{Data}/#LTr_{trk}^{D_{2}^{(#beta=1)}}#GT_{MC}"], [DoubleRatios_gaus[1]], ["Multi-jet b-tag"],      [TrackD2_D2Dn_Sys_btag], ["Existing D_{2}^{(#beta=1)} Uncertainty"] )
#DrawTool.ClearTexts()


DrawTool_DoubleRatio.SysLabels = ["Existing D_{2}^{(#beta=1)} Uncertainty"]
TrackD2_Ratios_Gaus[0].SetMarkerStyle(26)
TrackD2_Ratios_Gaus[1].SetMarkerStyle(26)
TrackD2_Ratios_Gaus[0].SetLineStyle(2)
TrackD2_Ratios_Gaus[1].SetLineStyle(2)
TrackD2_Ratios_Gaus[2].SetMarkerStyle(20)
TrackD2_Ratios_Gaus[3].SetMarkerStyle(20)

TrackD2_Ratios_Gaus[0].GetYaxis().SetTitleSize( TrackD2_Ratios_Gaus[0].GetYaxis().GetTitleSize()*1.1)
TrackD2_Ratios_Gaus[1].GetYaxis().SetTitleSize( TrackD2_Ratios_Gaus[0].GetYaxis().GetTitleSize()*1.1)
TrackD2_Ratios_Gaus[2].GetYaxis().SetTitleSize( TrackD2_Ratios_Gaus[0].GetYaxis().GetTitleSize()*1.1)
TrackD2_Ratios_Gaus[3].GetYaxis().SetTitleSize( TrackD2_Ratios_Gaus[0].GetYaxis().GetTitleSize()*1.1)

DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet before b-tagging")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool_DoubleRatio.PreCalRatio = DoubleRatios_gaus[0]
DrawTool_DoubleRatio.DrawHists("rtrk_D2", ["p_{T} [GeV]", "#LTr_{trk}^{D_{2}}#GT"],  [TrackD2_Ratios_Gaus[0], TrackD2_Ratios_Gaus[0], TrackD2_Ratios_Gaus[2]], ["Pythia 8 MC", "Pythia 8 MC", "Data"], [], [], [TrackD2_D2Up_Sys])
DrawTool_DoubleRatio.ClearTexts()

DrawTool_DoubleRatio.AddTexts(0.23, 0.73, 1, 0.048, "anti-k_{t} R=1.0 jet after b-tagging")
DrawTool_DoubleRatio.AddTexts(0.23, 0.66, 1, 0.048, "Trimmed (R_{sub}=0.3, f_{cut}=0.05)")
DrawTool_DoubleRatio.PreCalRatio = DoubleRatios_gaus[1]
DrawTool_DoubleRatio.DrawHists("rtrk_D2_btag", ["p_{T} [GeV]", "#LTr_{trk}^{D_{2}}#GT"],  [TrackD2_Ratios_Gaus[1], TrackD2_Ratios_Gaus[1], TrackD2_Ratios_Gaus[3]], ["Pythia 8 MC", "Pythia 8 MC", "Data"], [], [], [TrackD2_D2Dn_Sys_btag])
