import rootpy.tree
from rootpy.io import root_open
from ROOT import gDirectory, TLorentzVector, TVector3, TH1D, TH2D, TFile, TMath, Math
from HistoLab import HistoTool, CopyHist, StackHists, GetPlot, ScaleHist, Project2D, CopyHists, project1D, NormHists, AddHistPairs, ResetBins, GetHistMedian, GetHistIQR, GetRatios, GetAveRatios, GetGausMean, CalSysVariation, CalSysDiff
from FitFunctions import BernsteinO3, BernsteinO4, fit_start, fit_end, fit_range
import root_numpy as rnp
import numpy as np
from copy import deepcopy
from math import sqrt

def UniqueList(inList):
    out = []
    for i in inList:
        if i not in out:
            out.append(i)
    return out

def PoisonFluctuate(hist, binstart =None, binend = None):
    newhist = deepcopy(hist)
    for ibin in range(hist.GetNbinsX()):
        bincontent = np.random.poisson(hist.GetBinContent(ibin+1 ))
        newhist.SetBinContent(ibin+1, bincontent)
    return newhist


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


def ScanGetThreshold(prediction, evtweight, cut):
    order = np.argsort(prediction)
    weight_sorted = evtweight[order]
    prediction_sorted = prediction[order]

    weighttot = sum(weight_sorted)
    weightsum = 0
    for i in range(len(weight_sorted)):
        weightsum += weight_sorted[i]
        if weightsum>cut*weighttot:
            return prediction_sorted[i]

def FindOptimalBDTCut(mcpred, evtweight, datapred, data0tagpred, data0tagsideband, data0tagsignal, datafraction, nregion):
    order = np.argsort(mcpred)
    weight_sorted = evtweight[order]
    mcpred_sorted = mcpred[order]

    weighttot = sum(weight_sorted)
    weightsum = 0
    bottomboundary =0

    #find botoom
    for i in range(weight_sorted.shape[0]):
        weightsum += weight_sorted[i]
        if weightsum>(1-datafraction)*weighttot:
            bottomboundary = mcpred_sorted[i]
            weight_sorted = weight_sorted[i:-1]
            mcpred_sorted = mcpred_sorted[i:-1]

            break

    if nregion ==2:
        maxSensitivity = 0
        cut = -100
        s_over_b_1 = 0
        s_over_b_2 = 0
        B1 = 0
        B2 = 0
        S1 = 0
        S2 = 0

        for i in range(weight_sorted.shape[0]):
            S1 = np.sum(weight_sorted[0:i])
            B1_2tag  = np.sum( np.logical_and(datapred< mcpred_sorted[i], datapred>bottomboundary))
            B1_0tag  = np.sum(np.logical_and(np.logical_and(  data0tagpred< mcpred_sorted[i] , data0tagsideband), data0tagpred>bottomboundary) )

            B1  = B1_2tag/float(B1_0tag) * np.sum(np.logical_and( np.logical_and( data0tagpred< mcpred_sorted[i] , data0tagsignal), data0tagpred>bottomboundary))
            
            S2 = np.sum(weight_sorted[i:-1])

            B2_2tag = np.sum( datapred> mcpred_sorted[i])
            B2_0tag  = np.sum(np.logical_and( data0tagpred> mcpred_sorted[i] , data0tagsideband) )

            B2  = B2_2tag/float(B2_0tag) * np.sum(np.logical_and( data0tagpred> mcpred_sorted[i], data0tagsignal) )

            if S1 < 10: #0.1 * np.sum(weight_sorted):
                continue

            if S2 < 10: #0.1 * np.sum(weight_sorted):
                continue

            if B1 == float('nan') or B2 == float('nan'):
                continue


            sensitivity = sqrt( S1**2/(S1+B1) + S2**2/(S2+B2) )

            if sensitivity>maxSensitivity:
                maxSensitivity = sensitivity
                cut = mcpred_sorted[i]
                s_over_b_1 = sqrt(S1**2/(B1+S1))
                s_over_b_2 = sqrt(S2**2/(B2+S2))

        print "region II", s_over_b_1
        print "region I", s_over_b_2
        print bottomboundary, cut

        return [cut, bottomboundary]

    if nregion ==3:
        maxSensitivity = 0
        cut = -100
        s_over_b_1 = 0
        s_over_b_2 = 0
        s_over_b_3 = 0
        B1 = 0
        B2 = 0
        B3 = 0
        S1 = 0
        S2 = 0
        S3 = 0

        searchgrid = []
        step = 100
        
        for i in range(step):
            searchgrid.append( (max(mcpred_sorted)-min(mcpred_sorted))/step*i+ min(mcpred_sorted) )

        for i in range( len(searchgrid) ):
            score_i = searchgrid[i]

            S1 = np.sum( weight_sorted[ mcpred_sorted<score_i]  )
            B1_2tag  = np.sum( np.logical_and(datapred< score_i, datapred>bottomboundary))
            B1_0tag  = np.sum(np.logical_and(np.logical_and(  data0tagpred< score_i , data0tagsideband), data0tagpred>bottomboundary) )
            B1  = B1_2tag/float(B1_0tag) * np.sum(np.logical_and( np.logical_and( data0tagpred< score_i , data0tagsignal), data0tagpred>bottomboundary))

            for j in range(i+1, len(searchgrid)):
                score_j = searchgrid[j]

                S2 = np.sum(  weight_sorted[ np.logical_and(mcpred_sorted>score_i, mcpred_sorted<score_j) ] )
                B2_2tag = np.sum( np.logical_and( datapred> score_i, datapred< score_j))
                B2_0tag  = np.sum(np.logical_and( np.logical_and( data0tagpred> score_i ,data0tagpred< score_j), data0tagsideband))
                B2  = B2_2tag/float(B2_0tag) * np.sum(np.logical_and( np.logical_and( data0tagpred> score_i ,data0tagpred< score_j), data0tagsignal))

                S3 = np.sum(weight_sorted[ mcpred_sorted>score_j]  )

                B3_2tag = np.sum(datapred> score_j)
                B3_0tag  = np.sum(np.logical_and( data0tagpred> score_j, data0tagsideband))
                B3  = B3_2tag/float(B3_0tag) * np.sum(np.logical_and( data0tagpred> score_j, data0tagsignal))

                if S1 < 0.2 * np.sum(weight_sorted) or S2 < 0.2 * np.sum(weight_sorted) or S3 < 0.2 * np.sum(weight_sorted):
                    continue

                if B1 < 10 or B2<10 or B3<10:
                    continue

                if B3_0tag < 20 or B2_0tag < 20 or B1_0tag < 20:
                    continue
                
                if B1 == float('nan') or B2 == float('nan') or B3 == float('nan'):
                    continue

                sensitivity = sqrt( S1**2/(B1+S1) + S2**2/(B2+S2) + S3**2/(B3+S3) )

                if sensitivity>maxSensitivity:
                    maxSensitivity = sensitivity
                    cut1 = score_i
                    cut2 = score_j
                    s_over_b_1 = sqrt(S1**2/(B1+S1))
                    s_over_b_2 = sqrt(S2**2/(B2+S2))
                    s_over_b_3 = sqrt(S3**2/(B3+S3))

        print "region III", s_over_b_1
        print "region II", s_over_b_2
        print "region I", s_over_b_3
        print bottomboundary, cut1, cut2

        return [cut2, cut1, bottomboundary]


class Cut:
    def __init__(self, name, array):
        self.name = name
        self.array = array

    def __add__(self, cut2):
        return Cut( self.name+"_"+cut2.name, np.logical_and(self.array, cut2.array) )


NoCut = Cut("NoCut", np.array([1]))

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
        if cut == NoCut:
            self.histograms[histname] = TH1D(self.name+"_"+histname , self.name+"_"+histname, len(bin)-1, bin)
        else:
            histname = histname +"_"+cut.name
            self.histograms[histname] = TH1D(self.name+"_"+histname, self.name+"_"+histname, len(bin)-1, bin)

    def Add2DHist(self, histname, binx, biny, cut = NoCut):
        self.histograms[histname] = TH2D(self.name+"_"+histname +"_" +cut.name, self.name+"_"+histname+"_"+cut.name, len(binx), binx[0], binx[-1], len(biny), biny[0], biny[-1])

    def Norm1DHist(self, histname):
        self.histograms[histname+"_Normed"] = NormHists([self.histograms[histname]])[0]

    def FillHist(self, histname, array, weight, cut = NoCut):
        if cut == NoCut:
            rnp.fill_hist( self.histograms[histname], array, weight)
        else:
            rnp.fill_hist( self.histograms[histname+"_"+cut.name], array, weight*cut.array)

    def AddCut(self, name, array):
        self.cut[name] = Cut(name, array)

    def AddNoCut(self):
        self.cut["NoCut"] = NoCut

    def AddEventVarFromTree(self, varname, test=False):
        stop = None
        if test:
            stop = 300000

        newarray = rnp.tree2array(self.tree, varname, stop=stop)

        if varname == "pT_ballance":
            varname = "pT_balance"

        VarInMeV = ["mBB", "pTBB", "pTJJ", "mJJ", "mJ1B1", "mJ1B2", "deltaMJJ", "HT_MVA", "HT_soft", "pTB1", "pTB2", "pTJ1", "pTJ2"]

        if varname in VarInMeV:
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
    def __init__(self, fitname, region, channel, fitfunction, ndof, period):
        self.fitname = fitname
        self.region = region
        self.channel = channel
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
        self.period = period

    def GetFitPar(self):
        for ipar in range(self.fit.GetNumberFreeParameters()):
            self.fitpar.append(self.fit.GetParameter(ipar))
            if self.fit.GetParameter(ipar)>0:
                self.fitparup.append( 10*self.fit.GetParameter(ipar))
                self.fitpardn.append( 0)
            if self.fit.GetParameter(ipar)<0:
                self.fitparup.append( 0)
                self.fitpardn.append( 10*self.fit.GetParameter(ipar))

    def GetxmlForm(self):

        self.GetFitPar()
        functiontext = ""
        
        if self.fitname == "BernsteinO3":
            functiontext = '''<ModelItem Name="EXPR::Bernstein3{!s}('@1*pow((1-@0), 3) + 3*@2*@0*pow((1-@0), 2) + 3*@3*(1-@0)*pow(@0, 2) + @4*pow(@0, 3)', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3])

            if "SR" in self.region:
                functiontext = '''<ModelItem Name="EXPR::Bernstein3{!s}('(@1*pow((1-@0), 3) + 3*@2*@0*pow((1-@0), 2) + 3*@3*(1-@0)*pow(@0, 2) + @4*pow(@0, 3))*(@5*@6+1)', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], :observable:, Lin0{!s}[-0.0001, -0.01, 1])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel+self.region, self.channel+self.region)


        if self.fitname == "BernsteinO4":
            functiontext = '''<ModelItem Name="EXPR::Bernstein4{!s}('@1*pow((1-@0), 4) + 4*@2*@0*pow((1-@0), 3) + 6*@3*pow(@0, 2)*pow((1-@0),2) + 4*@4*(1-@0)*pow(@0, 3) + @5*pow(@0, 4)', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s},{!s}])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4])


            if "SR" in self.region:
                functiontext = '''<ModelItem Name="EXPR::Bernstein4{!s}('(@1*pow((1-@0), 4) + 4*@2*@0*pow((1-@0), 3) + 6*@3*pow(@0, 2)*pow((1-@0),2) + 4*@4*(1-@0)*pow(@0, 3) + @5*pow(@0, 4))*(@6*@7+1)', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s},{!s}], :observable:, Lin0{!s}[-0.0001, -0.01, 1])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4], self.channel+self.region, self.channel+self.region)


        if self.fitname == "BernsteinO5":
            functiontext = '''<ModelItem Name="EXPR::Bernstein5{!s}('@1*pow((1-@0), 5) + 5*@2*@0*pow((1-@0), 4) + 10*@3*pow(@0, 2)*pow((1-@0),3) + 10*@4*pow(@0, 3)*pow((1-@0), 2) + 5*@5*(1-@0)*pow(@0, 4) + @6*pow(@0, 5)', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s},{!s}], x6{!s}[{!s}, {!s},{!s}])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4], self.channel, self.fitpar[5], self.fitpardn[5],self.fitparup[5])

            if "SR" in self.region:
                functiontext = '''<ModelItem Name="EXPR::Bernstein5{!s}('(@1*pow((1-@0), 5) + 5*@2*@0*pow((1-@0), 4) + 10*@3*pow(@0, 2)*pow((1-@0),3) + 10*@4*pow(@0, 3)*pow((1-@0), 2) + 5*@5*(1-@0)*pow(@0, 4) + @6*pow(@0, 5))*(@7*@8+1)', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s},{!s}], x6{!s}[{!s}, {!s},{!s}], :observable:, Lin0{!s}[-0.0001, -0.01, 1])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4], self.channel, self.fitpar[5], self.fitpardn[5],self.fitparup[5], self.channel+self.region, self.channel+self.region)


        if self.fitname == "ExpoBernsteinO2":
            functiontext = '''<ModelItem Name="EXPR::ExpoBernstein2{!s}('exp(@1*@0)*( @2*pow((1-@0), 2) + 2*@3*@0*(1-@0) + @4*(pow(@0, 2)) )', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3])

            if "SR" in self.region:
                functiontext = '''<ModelItem Name="EXPR::ExpoBernstein2{!s}('exp(@1*@0)*( @2*pow((1-@0), 2) + 2*@3*@0*(1-@0) + @4*(pow(@0, 2)))(@5*@6+1)', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], :observable:, Lin0{!s}[-0.0001, -0.01, 1])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel+self.region)


        if self.fitname == "ExpoBernsteinO3":
            functiontext = '''<ModelItem Name="EXPR::ExpoBernstein3{!s}('exp(@1*@0)*(@2*pow((1-@0), 3) + 3*@3*@0*pow((1-@0), 2) + 3*@4*(1-@0)*pow(@0, 2) + @5*pow(@0, 3))', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s}, {!s}])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4])
            
            if "SR" in self.region:
                functiontext = '''<ModelItem Name="EXPR::ExpoBernstein3{!s}('exp(@1*@0)*(@2*pow((1-@0), 3) + 3*@3*@0*pow((1-@0), 2) + 3*@4*(1-@0)*pow(@0, 2) + @5*pow(@0, 3))*(@6*@7+1)', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s}, {!s}], :observable:, Lin0{!s}[-0.0001, -0.1, 1])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4], self.channel+self.region)


        if self.fitname == "ExpoBernsteinO4":
            functiontext = '''<ModelItem Name="EXPR::ExpoBernstein4{!s}('exp(@1*@0)*(@2*pow((1-@0), 4) + 4*@3*@0*pow((1-@0), 3) + 6*@4*pow(@0, 2)*pow((1-@0),2) + 4*@5*(1-@0)*pow(@0, 3) + @6*pow(@0, 4))', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s},{!s}], x6{!s}[{!s}, {!s},{!s}])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4], self.channel, self.fitpar[5], self.fitpardn[5],self.fitparup[5])
            if "SR" in self.region:
                functiontext = '''<ModelItem Name="EXPR::ExpoBernstein4{!s}('exp(@1*@0)*(@2*pow((1-@0), 4) + 4*@3*@0*pow((1-@0), 3) + 6*@4*pow(@0, 2)*pow((1-@0),2) + 4*@5*(1-@0)*pow(@0, 3) + @6*pow(@0, 4))*(@7*@8+1)', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s},{!s}], x6{!s}[{!s}, {!s},{!s}], :observable:, Lin0{!s}[-0.0001, -0.1, 1])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4], self.channel, self.fitpar[5], self.fitpardn[5],self.fitparup[5], self.channel+self.region)


        if self.fitname == "ExpoPolO2":
            functiontext = '''<ModelItem Name="EXPR::ExpoPolO2{!s}('exp( -(@1 + @2*@0+ @3*pow(@0, 2) )', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2])

            if "SR" in self.region:
                functiontext = '''<ModelItem Name="EXPR::ExpoPolO2{!s}('exp( -(@1 + @2*@0+ @3*pow(@0, 2)*(@4*@5+1) )', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], :observable:, Lin0{!s}[-0.0001, -0.1, 1])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel+self.region)

        if self.fitname == "ExpoPolO3":
            functiontext = '''<ModelItem Name="EXPR::ExpoPolO3{!s}('exp( -(@1 + @2*@0+ @3*pow(@0, 2) +@4*pow(@0,3) )', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3])

            if "SR" in self.region:
                functiontext = '''<ModelItem Name="EXPR::ExpoPolO3{!s}('exp( -(@1 + @2*@0+ @3*pow(@0, 2) +@4*pow(@0,3)*(@5*@6+1) )', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], :observable:, Lin0{!s}[-0.0001, -0.1, 1])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel+self.region)

        if self.fitname == "ExpoPolO4":
            functiontext = '''<ModelItem Name="EXPR::ExpoPolO4{!s}('exp( -(@1 + @2*@0+ @3*pow(@0, 2) +@4*pow(@0,3) +@5*pow(@0,4) )', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s}, {!s}])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4])

            if "SR" in self.region:
                functiontext = '''<ModelItem Name="EXPR::ExpoPolO4{!s}('exp( -(@1 + @2*@0+ @3*pow(@0, 2) +@4*pow(@0,3) +@5*pow(@0,4)*(@6*@7+1) )', x, x1{!s}[{!s}, {!s}, {!s}], x2{!s}[{!s}, {!s}, {!s}], x3{!s}[{!s}, {!s}, {!s}], x4{!s}[{!s}, {!s}, {!s}], x5{!s}[{!s}, {!s}, {!s}], :observable:, Lin0{!s}[-0.0001, -0.1, 1])"/>'''.format(self.region, self.channel, self.fitpar[0], self.fitpardn[0],self.fitparup[0], self.channel, self.fitpar[1], self.fitpardn[1],self.fitparup[1], self.channel, self.fitpar[2], self.fitpardn[2],self.fitparup[2], self.channel, self.fitpar[3], self.fitpardn[3],self.fitparup[3], self.channel, self.fitpar[4], self.fitpardn[4],self.fitparup[4], self.channel+self.region)




        self.xmlform = functiontext

        fit_results_xml = open("FitFiles_"+self.channel+self.period+"/bkg_"+self.fitname+"_"+self.region+".xml", "a")

        fit_results_xml.write('''<!DOCTYPE Model SYSTEM '../AnaWSBuilder.dtd'>'''+"\n")
        fit_results_xml.write('''<Model Type="UserDef">'''+"\n")
        fit_results_xml.write('''  <Item Name="expr::x('((@0-@2)/@1)', :observable:, range[{!s}], shift[{!s}])"/>'''.format(fit_range, fit_start)+"\n")
        fit_results_xml.write("  "+self.xmlform+"\n")
        fit_results_xml.write('''</Model>''')

        print self.region
        print "FitName", self.fitname
        print "Fit Chi2 prob", self.chinprob
        print "Fit Chi2 chi", self.chindof



def CalFtest(fit1, fit2):
    print fit1.fitname, fit2.fitname
    print 1.-Math.fdistribution_cdf( fit1.chindof/fit2.chindof, fit1.ndof, fit2.ndof)



