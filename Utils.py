import rootpy.tree
from rootpy.io import root_open
from ROOT import gDirectory, TLorentzVector, TVector3, TH1D, TH2D, TFile, TMath, Math
from HistoLab import HistoTool, CopyHist, StackHists, GetPlot, ScaleHist, Project2D, CopyHists, project1D, NormHists, AddHistPairs, ResetBins, GetHistMedian, GetHistIQR, GetRatios, GetAveRatios, GetGausMean, CalSysVariation, CalSysDiff
from FitFunctions import BernsteinO3, BernsteinO4
import root_numpy as rnp
import numpy as np
from copy import deepcopy


def UniqueList(inList):
    out = []
    for i in inList:
        if i not in out:
            out.append(i)
    return out

def CalWeightedCorrCoeff(x, y, w):
    mu_x = np.average(x, weights= w)
    mu_y = np.average(y, weights= w)

    x_minus_mu = x-mu_x
    y_minus_mu = y-mu_y

    num = np.sum( x_minus_mu * y_minus_mu * w)

    denum = np.sqrt( np.sum(x_minus_mu*x_minus_mu*w) * np.sum(y_minus_mu*y_minus_mu*w) )
    return num/denum
    

def sort_jet_var(vararray, njet ):

    newarray = []
    for i in range(vararray.shape[0]):
        if vararray[i].shape[0]>njet:
            newarray.append(vararray[i][njet])
        else:
            newarray.append(-999)

    return np.array(newarray)


def sort_jet_vars(vararrays, njet ):

    newarrays = []

    for ivar in range(len(vararrays)):
        newarrays.append([])
    
    for ievt in range(vararrays[0].shape[0]):
        for ivar in range(len(vararrays)):

            if vararrays[ivar][ievt].shape[0]>njet:

                newarrays[ivar].append(vararrays[ivar][ievt][njet])

            else:
                newarrays[ivar].append(-999)

    for ivar in range(len(vararrays)):
        newarrays[ivar] = np.array(newarrays[ivar])

    return newarrays


def cal_jet_min_dR(sample, njet):

    newarray = []

    for ievt in range(sample.var["nJets"].shape[0]):

        if sample.var["nJets"][ievt]>njet:
            mindR = 999

            for ijet in range(sample.var["nJets"][ievt]):
                if ijet == njet:
                    continue
                else:
                    dR  = (sample.var["jet_eta_"+str(ijet)][ievt] - sample.var["jet_eta_"+str(njet)][ievt])**2
                    dR += (sample.var["jet_phi_"+str(ijet)][ievt] - sample.var["jet_phi_"+str(njet)][ievt])**2
                if mindR >dR:
                    mindR = sqrt(dR)

            newarray.append(mindR)

        else:
            newarray.append(0)

    sample.var["jet_mindR_"+str(ijet)] = newarray
    return np.array(newarray)


class Cut:
    def __init__(self, name, array):
        self.name = name
        self.array = array

    def __add__(self, cut2):
        return Cut( self.name+"_"+cut2.name, np.logical_and(self.array, cut2.array) )


NoCut = Cut("NoCut", np.array([0]))

class PhysicsProcess:
    def __init__(self, name, filename, sysName = ""):
        self.name = name
        self.filename = filename
        self.var = {}
        self.cut ={}
        self.histograms ={}
        self.file = None
        self.tree = None
        self.isSys = False
        if filename != None:
            self.file = root_open(filename) 
            if sysName =="":
                self.tree = self.file.Nominal
            else:
                self.tree = self.file.Get(sysName)
                self.isSys = True

    def Add1DHist(self, histname, bin, cut = NoCut):
        self.histograms[histname] = TH1D(self.name+"_"+histname +"_" +cut.name, self.name+"_"+histname+"_"+cut.name, len(bin)-1, bin)

    def Add2DHist(self, histname, binx, biny, cut = NoCut):
        self.histograms[histname] = TH2D(self.name+"_"+histname +"_" +cut.name, self.name+"_"+histname+"_"+cut.name, len(binx), 0, binx[-1], len(biny), 0, biny[-1])

    def Norm1DHist(self, histname):
        self.histograms[histname+"_Normed"] = NormHists([self.histograms[histname]])[0]

    def FillHist(self, histname, array, weight, cut = NoCut):
        if cut == NoCut:
            rnp.fill_hist( self.histograms[histname], array, weight)
        else:
            rnp.fill_hist( self.histograms[histname], array, weight*cut.array)

    def AddCut(self, name, array):
        self.cut[name] = Cut(name, array)

    def AddEventVarFromTree(self, varname, test=False):
        stop = None
        if test:
            stop = 3000000

        newarray = rnp.tree2array(self.tree, varname, stop=stop)
        if "mBB" in varname or "pTBB" in varname or "pTJJ" in varname or "mJJ" in varname or "HT" in varname:
            self.var[varname] = newarray/1000.
        else:
            self.var[varname] = newarray

    def AddJetVarFromTree(self, varname, njet):

        jetarray = rnp.tree2array(self.tree, varname)
        newarray = sort_jet_var(jetarray, njet)

        if "pT" in varname:
            newarray = newarray/1000.

        self.var["jet_"+varname+"_"+str(njet)] = newarray

    def AddJetVarsFromTree(self, varnames, njet):
        
        jetarrays = []
        for varname in varnames:
            jetarray = rnp.tree2array(self.tree, varname)
            jetarrays.append(jetarray)

        newarrays = sort_jet_vars(jetarrays, njet)

        for ivar in range(len(newarrays)):
            newarray = newarrays[ivar]

            if "pT" in varnames[ivar]:
                newarray = newarray/1000.

            self.var["jet_"+varnames[ivar]+"_"+str(njet)] = newarray


class FitResults:
    def __init__(self, fitname, region, fitfunction, ndof):
        self.fitname = fitname
        self.region = region
        self.fit = deepcopy(fitfunction)
        print fitfunction
        self.xmlform = ""
        self.chi = fitfunction.GetChisquare()
        self.ndof = ndof
        self.chindof = fitfunction.GetChisquare()/ndof
        self.chinprob = TMath.Prob( fitfunction.GetChisquare(), ndof)
        self.fitpar = []
        self.fitparup = []
        self.fitpardn = []

    def GetFitPar(self):
        for ipar in range(self.fit.GetNumberFreeParameters()):
            self.fitpar.append(self.fit.GetParameter(ipar))
            if self.fit.GetParameter(ipar)>0:
                self.fitparup.append( 10*self.fit.GetParameter(ipar))
                self.fitpardn.append( -10*self.fit.GetParameter(ipar))
            if self.fit.GetParameter(ipar)<0:
                self.fitparup.append( -10*self.fit.GetParameter(ipar))
                self.fitpardn.append( 10*self.fit.GetParameter(ipar))

    def GetxmlForm(self):

        self.GetFitPar()
        text = ""
        
        if self.fitname == "BernsteinO3":
            text = '''<ModelItem Name="EXPR::Bernstein3{!s}('@1*pow((1-@0), 3) + 3*@2*@0*pow((1-@0), 2) + 3*@3*(1-@0)*pow(@0, 2) + @4*pow(@0, 3)', x, x11[{!s}, {!s}, {!s}], x12[{!s}, {!s}, {!s}], x13[{!s}, {!s}, {!s}], x14[{!s}, {!s}, {!s}])"/>'''.format(self.region, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.fitpar[3], self.fitpardn[3],self.fitparup[3])

        if self.fitname == "BernsteinO4":
            text = '''<ModelItem Name="EXPR::Bernstein4{!s}('@1*pow((1-@0), 4) + 4*@2*@0*pow((1-@0), 3) + 6*@3*pow(@0, 2)*pow((1-@0),2) + 4*@4*(1-@0)*pow(@0, 3) + @5*pow(@0, 4)', x, x11[{!s}, {!s}, {!s}], x12[{!s}, {!s}, {!s}], x13[{!s}, {!s}, {!s}], x14[{!s}, {!s}, {!s}], x15[{!s}, {!s},{!s}])"/>'''.format(self.region, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.fitpar[4], self.fitpardn[4],self.fitparup[4])

        if self.fitname == "Expo":
            text = '''<ModelItem Name="EXPR::EXP{!s}('exp(@1+@0*@2)', :observable:, x11[{!s}, {!s}, {!s}], x12[{!s}, {!s}, {!s}])"/>'''.format(self.region, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.fitpar[1], self.fitpardn[1],self.fitparup[1])
            
        self.xmlform = text

        fit_results_txt = open("f_"+self.fitname+"_"+self.region, "a")

        for ipar in range(self.fit.GetNumberFreeParameters()):
            fit_results_txt.write(str(self.fit.GetParameter(ipar))+"\n")

        fit_results_txt.write("chi2 "+str(self.chi)+"\n")
        fit_results_txt.write("chi2 ndof "+str(self.chindof)+"\n")
        fit_results_txt.write("chi2 Prob "+str(self.chinprob)+"\n")
        fit_results_txt.write(self.xmlform+"\n")

        print self.region
        print "FitName", self.fitname
        print "Fit Chi2 prob", self.chinprob
        print "Fit Chi2 chi", self.chindof


def CalFtest(fit1, fit2):
    print fit1.fitname, fit2.fitname
    print 1.-Math.fdistribution_cdf( fit1.chindof/fit2.chindof, fit1.ndof, fit2.ndof)
