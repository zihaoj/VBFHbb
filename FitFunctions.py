from ROOT import TF1, TH1D, TFile, Long, Double, TMinuit, RooRealVar, RooDataSet, RooGaussian, RooCBShape, RooAddPdf, RooExtendPdf, RooPlot, RooArgSet, RooArgList, TCanvas, RooFit
import random as rand
from math import sqrt, exp, log
import array
from cStringIO import StringIO
import root_numpy as rnp
import numpy as np

fit_start = 80
fit_end = 200
fit_range = fit_end-fit_start

################## Crystalball function
def CrystalBall_raw(x, par):
    if par[3]<0:
        return 0

    t = (x[0]-par[2])/par[3];
    if (par[0] < 0):
        t = -t

    abs_alpha = abs(par[0])

    if (t >= -abs_alpha):
        return par[4]*exp(-0.5*t*t)

    else:

        nDivAlpha = par[1]/abs_alpha
        AA = exp(-0.5*abs_alpha*abs_alpha)
        B = nDivAlpha-abs_alpha
        arg = nDivAlpha/(B-t)

        return par[4]*( arg**par[1])

CrystalBall = TF1("CrystalBall", CrystalBall_raw, 40, 150, 5)


def CrystalBallGaus_raw(x, par):

    if par[3]<0:
        return 0

    t = (x[0]-par[2])/par[3];
    if (par[0] < 0):
        t = -t

    abs_alpha = abs(par[0])

    if (t >= -abs_alpha):
        return par[4]*exp(-0.5*t*t  + exp(- (x[0]-par[5])**2/(2*par[6]**2) ))

    else:

        nDivAlpha = par[1]/abs_alpha
        AA = exp(-0.5*abs_alpha*abs_alpha)
        B = nDivAlpha-abs_alpha
        arg = nDivAlpha/(B-t)

        return par[4]*( arg**par[1]  + exp(- (x[0]-par[5])**2/(2*par[6]**2) ))


CrystalBallGaus = TF1("CrystalBallGaus", CrystalBallGaus_raw, 40, 150, 7)


################### Bernstein functions
def BernsteinO2_raw(x, par):
    return par[0]*(1-(x[0]-fit_start)/fit_range)**2+ 2*par[1]*(1-(x[0]-fit_start)/fit_range)*((x[0]-fit_start)/fit_range) +  par[2]*((x[0]-fit_start)/fit_range)**2

BernsteinO2 = TF1("BernsteinO2", BernsteinO2_raw, fit_start, fit_end, 3)


def BernsteinO3_raw(x, par):
    return par[0]*(1-((x[0]-fit_start)/fit_range))**3 + par[1]*(3*((x[0]-fit_start)/fit_range)*(1-((x[0]-fit_start)/fit_range))**2) + par[2]*(3*((x[0]-fit_start)/fit_range)**2*(1-((x[0]-fit_start)/fit_range))) + par[3]*((x[0]-fit_start)/fit_range)**3

BernsteinO3 = TF1("BernsteinO3", BernsteinO3_raw, fit_start, fit_end, 4)


def BernsteinO4_raw(x, par):
    return par[0]*(1-((x[0]-fit_start)/fit_range))**4 + par[1]*(4*((x[0]-fit_start)/fit_range)*(1-((x[0]-fit_start)/fit_range))**3) + par[2]*(6*((x[0]-fit_start)/fit_range)**2*(1-((x[0]-fit_start)/fit_range))**2) + par[3]*(4*((x[0]-fit_start)/fit_range)**3*(1-((x[0]-fit_start)/fit_range))) + par[4]*((x[0]-fit_start)/fit_range)**4

BernsteinO4 = TF1("BernsteinO4", BernsteinO4_raw, fit_start, fit_end, 5)


def BernsteinO5_raw(x, par):
    return par[0]*(1-((x[0]-fit_start)/fit_range))**5 + par[1]*(5*((x[0]-fit_start)/fit_range)*(1-((x[0]-fit_start)/fit_range))**4) + par[2]*(10*((x[0]-fit_start)/fit_range)**2*(1-((x[0]-fit_start)/fit_range))**3) + par[3]*(10*((x[0]-fit_start)/fit_range)**3*(1-((x[0]-fit_start)/fit_range))**2) + par[4]*(5*((x[0]-fit_start)/fit_range)**4*(1-((x[0]-fit_start)/fit_range))) + par[5]*((x[0]-fit_start)/fit_range)**5

BernsteinO5 = TF1("BernsteinO5", BernsteinO5_raw, fit_start, fit_end, 6)


def BernsteinO6_raw(x, par):
    return (par[0]*(1-((x[0]-fit_start)/fit_range))**6 + 
            par[1]*(6*((x[0]-fit_start)/fit_range)**1*(1-((x[0]-fit_start)/fit_range))**5) + 
            par[2]*(15*((x[0]-fit_start)/fit_range)**2*(1-((x[0]-fit_start)/fit_range))**4) + 
            par[3]*(20*((x[0]-fit_start)/fit_range)**3*(1-((x[0]-fit_start)/fit_range))**3) + 
            par[4]*(15*((x[0]-fit_start)/fit_range)**4*(1-((x[0]-fit_start)/fit_range))**2) + 
            par[5]*(6*((x[0]-fit_start)/fit_range)**5*(1-((x[0]-fit_start)/fit_range))**1) + 
            par[6]*((x[0]-fit_start)/fit_range)**6)

BernsteinO6 = TF1("BernsteinO6", BernsteinO6_raw, fit_start, fit_end, 7)


###################### Bernstein * exponentials
def ExpoBernsteinO2_raw(x, par):
    try:
        return exp(par[0]*(x[0]-fit_start)/fit_range)* (par[1]*(1-(x[0]-fit_start)/fit_range)**2+ 2*par[2]*(1-(x[0]-fit_start)/fit_range)*((x[0]-fit_start)/fit_range) +  par[3]*((x[0]-fit_start)/fit_range)**2)
    except OverflowError:
        return -1

ExpoBernsteinO2 = TF1("ExpoBernsteinO2", ExpoBernsteinO2_raw, fit_start, fit_end, 4)


def ExpoBernsteinO3_raw(x, par):
    try:
        return exp(par[0]*(x[0]-fit_start)/fit_range)* (par[1]*(1-((x[0]-fit_start)/fit_range))**3 + par[2]*(3*((x[0]-fit_start)/fit_range)*(1-((x[0]-fit_start)/fit_range))**2) + par[3]*(3*((x[0]-fit_start)/fit_range)**2*(1-((x[0]-fit_start)/fit_range))) + par[4]*((x[0]-fit_start)/fit_range)**3)
    except OverflowError:
        return -1

ExpoBernsteinO3 = TF1("ExpoBernsteinO3", ExpoBernsteinO3_raw, fit_start, fit_end, 5)


def ExpoBernsteinO4_raw(x, par):
    try:
        return exp(par[0]*(x[0]-fit_start)/fit_range)* (par[1]*(1-((x[0]-fit_start)/fit_range))**4 + par[2]*(4*((x[0]-fit_start)/fit_range)*(1-((x[0]-fit_start)/fit_range))**3) + par[3]*(6*((x[0]-fit_start)/fit_range)**2*(1-((x[0]-fit_start)/fit_range))**2) + par[4]*(4*((x[0]-fit_start)/fit_range)**3*(1-((x[0]-fit_start)/fit_range))) + par[5]*((x[0]-fit_start)/fit_range)**4)
    except OverflowError:
        return -1

ExpoBernsteinO4 = TF1("ExpoBernsteinO4", ExpoBernsteinO4_raw, fit_start, fit_end, 6)


def ExpoBernsteinO5_raw(x, par):
    try:
        return exp(par[0]*(x[0]-fit_start)/fit_range)* (par[1]*(1-((x[0]-fit_start)/fit_range))**5 + par[2]*(5*((x[0]-fit_start)/fit_range)*(1-((x[0]-fit_start)/fit_range))**4) + par[3]*(10*((x[0]-fit_start)/fit_range)**2*(1-((x[0]-fit_start)/fit_range))**3) + par[4]*(10*((x[0]-fit_start)/fit_range)**3*(1-((x[0]-fit_start)/fit_range))**2) + par[5]*(5*((x[0]-fit_start)/fit_range)**4*(1-((x[0]-fit_start)/fit_range))) + par[6]*((x[0]-fit_start)/fit_range)**5)
    except OverflowError:
        return -1

ExpoBernsteinO5 = TF1("ExpoBernsteinO5", ExpoBernsteinO5_raw, fit_start, fit_end, 7)


###################### Exponentials of polynomials
def ExpoPolO2_raw(x, par):
    return exp( - (par[0] + par[1]* ((x[0]-fit_start)/fit_range) + par[2]* ((x[0]-fit_start)/fit_range)**2 ) )

ExpoPolO2 = TF1("ExpoPolO2", ExpoPolO2_raw, fit_start, fit_end, 3)


def ExpoPolO3_raw(x, par):
    return exp( - (par[0] + par[1]* ((x[0]-fit_start)/fit_range) + par[2]* ((x[0]-fit_start)/fit_range)**2 + par[3]* ((x[0]-fit_start)/fit_range)**3 ) )

ExpoPolO3 = TF1("ExpoPolO3", ExpoPolO3_raw, fit_start, fit_end, 4)


def ExpoPolO4_raw(x, par):
    return exp( - (par[0] + par[1]* ((x[0]-fit_start)/fit_range) + par[2]* ((x[0]-fit_start)/fit_range)**2 + par[3]* ((x[0]-fit_start)/fit_range)**3 + par[4]* ((x[0]-fit_start)/fit_range)**4 ) )

ExpoPolO4 = TF1("ExpoPolO4", ExpoPolO4_raw, fit_start, fit_end, 5)


def LinCustomRange_raw(x, par):
    return par[0]+par[1]*x[0]

LinCustomRange = TF1("LinCustomRange", LinCustomRange_raw, fit_start, fit_end, 2)


####################### Fitting Utils

def HistToList (hist):
    hist_list = []
    for i in range(hist.GetNbinsX()):
        hist_list.append(hist.GetBinContent(i+1))
    return hist_list

def HistErrorList (hist):
    hist_list = []
    for i in range(hist.GetNbinsX()):
        hist_list.append(hist.GetBinError(i+1))
    return hist_list

def HistBinsToList (hist):
    bins_list = []
    for i in range(hist.GetNbinsX()):
        bins_list.append(hist.GetXaxis().GetBinCenter(i+1))
    return bins_list


#### main function used to run the minimization scheme
def bkgfit(data_hist, bkgfunction, bkgname, doFloatZ = False, signal_hist = None, z_hist = None):
    isBkgPlusZFit = False
    isSpuriousFit = False

    binning = HistBinsToList(data_hist)
    data_x = HistToList(data_hist)
    data_error = HistErrorList(data_hist)

    z_x = []
    signal_x = []

    if z_hist !=None:
        isBkgPlusZFit = True
        z_x = HistToList(z_hist)

    if signal_hist != None:
        isSpuriousFit = True
        signal_x = HistToList(signal_hist)

    parfunction = bkgfunction.GetNumberFreeParameters()
    partot = bkgfunction.GetNumberFreeParameters() + 2

    ### the fucntion used for TMinuit
    def fcn(npar, gin, f, par, ifag):
        L = 0
    
        # calculate likelihood, input par[0] is the N_B, par[1] is N_C, par[2] is N_L
        for ibin in range(len(binning)):
            if (data_x[ibin] <0.5):
                continue

            bincen = binning[ibin]
            
            bkg = 0
            data = data_x[ibin]

            if bkgname == "BernsteinO2":
                bkg = (par[0]*(1-(bincen-fit_start)/fit_range)**2
                       + 2*par[1]*(1-(bincen-fit_start)/fit_range)*((bincen-fit_start)/fit_range)
                       + par[2]*((bincen-fit_start)/fit_range)**2)

            if bkgname == "BernsteinO3":
                bkg = par[0]*(1-((bincen-fit_start)/fit_range))**3 + par[1]*(3*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**2) + par[2]*(3*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))) + par[3]*((bincen-fit_start)/fit_range)**3

            if bkgname == "BernsteinO4":
                bkg = par[0]*(1-((bincen-fit_start)/fit_range))**4 + par[1]*(4*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**3) + par[2]*(6*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))**2) + par[3]*(4*((bincen-fit_start)/fit_range)**3*(1-((bincen-fit_start)/fit_range))) + par[4]*((bincen-fit_start)/fit_range)**4

            if bkgname == "BernsteinO5":
                bkg = par[0]*(1-((bincen-fit_start)/fit_range))**5 + par[1]*(5*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**4) + par[2]*(10*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))**3) + par[3]*(10*((bincen-fit_start)/fit_range)**3*(1-((bincen-fit_start)/fit_range))**2) + par[4]*(5*((bincen-fit_start)/fit_range)**4*(1-((bincen-fit_start)/fit_range))) + par[5]*((bincen-fit_start)/fit_range)**5
            
            if bkgname ==  "BernsteinO6":
                bkg= (par[0]*(1-((bincen-fit_start)/fit_range))**6 + 
                      par[1]*(6*((bincen-fit_start)/fit_range)**1*(1-((bincen-fit_start)/fit_range))**5) + 
                      par[2]*(15*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))**4) + 
                      par[3]*(20*((bincen-fit_start)/fit_range)**3*(1-((bincen-fit_start)/fit_range))**3) + 
                      par[4]*(15*((bincen-fit_start)/fit_range)**4*(1-((bincen-fit_start)/fit_range))**2) + 
                      par[5]*(6*((bincen-fit_start)/fit_range)**5*(1-((bincen-fit_start)/fit_range))**1) + 
                      par[6]*((bincen-fit_start)/fit_range)**6)

            if bkgname == "ExpoBernsteinO2":
                try:
                    bkg = exp(par[0]*(bincen-fit_start)/fit_range)* (par[1]*(1-(bincen-fit_start)/fit_range)**2+ 2*par[2]*(1-(bincen-fit_start)/fit_range)*((bincen-fit_start)/fit_range) +  par[3]*((bincen-fit_start)/fit_range)**2)
                except OverflowError:
                    bkg = 0

            if bkgname == "ExpoBernsteinO3":
                try:
                    bkg = exp(par[0]*(bincen-fit_start)/fit_range)* (par[1]*(1-((bincen-fit_start)/fit_range))**3 + par[2]*(3*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**2) + par[3]*(3*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))) + par[4]*((bincen-fit_start)/fit_range)**3)
                except OverflowError:
                    bkg = 0

            if bkgname == "ExpoBernsteinO4":
                try:
                    bkg = exp(par[0]*(bincen-fit_start)/fit_range)* (par[1]*(1-((bincen-fit_start)/fit_range))**4 + par[2]*(4*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**3) + par[3]*(6*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))**2) + par[4]*(4*((bincen-fit_start)/fit_range)**3*(1-((bincen-fit_start)/fit_range))) + par[5]*((bincen-fit_start)/fit_range)**4)
                except OverflowError:
                    bkg = 0

            if bkgname == "ExpoBernsteinO5":
                try:
                    bkg = exp(par[0]*(bincen-fit_start)/fit_range)* (par[1]*(1-((bincen-fit_start)/fit_range))**5 + par[2]*(5*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**4) + par[3]*(10*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))**3) + par[4]*(10*((bincen-fit_start)/fit_range)**3*(1-((bincen-fit_start)/fit_range))**2) + par[5]*(5*((bincen-fit_start)/fit_range)**4*(1-((bincen-fit_start)/fit_range))) + par[6]*((bincen-fit_start)/fit_range)**5)
                except OverflowError:
                    bkg = 0


            if bkgname == "ExpoPolO2":
                bkg = exp( - (par[0] + par[1]* ((bincen-fit_start)/fit_range) + par[2]* ((bincen-fit_start)/fit_range)**2 ) )

            if bkgname == "ExpoPolO3":
                bkg = exp( - (par[0] + par[1]* ((bincen-fit_start)/fit_range) + par[2]* ((bincen-fit_start)/fit_range)**2 + par[3]* ((bincen-fit_start)/fit_range)**3 ) )

            if bkgname == "ExpoPolO4":
                bkg = exp( - (par[0] + par[1]* ((bincen-fit_start)/fit_range) + par[2]* ((bincen-fit_start)/fit_range)**2 + par[3]* ((bincen-fit_start)/fit_range)**3 + par[4]* ((bincen-fit_start)/fit_range)**4 ) )
            

            mu_x =  bkg
            
            #if isBkgPlusZFit:
            #    mu_x = mu_x + (par[partot-1] *z_x[ibin])

            if isSpuriousFit:
                mu_x =  mu_x + par[partot-2] *signal_x[ibin]

            #L = L + mu_x - data*log(mu_x)
                
            L = L + ((mu_x - data)/data_error[ibin])**2
            
        f[0] = L

    # initialize the TMinuit object
    arglist_p = 10 * [0]

    arglist = array.array('d')
    arglist.fromlist(arglist_p)
    ierflag = Long(0)
    maxiter = 1000000000
    
    arglist_p= [1]
    gMinuit = TMinuit(partot)
    gMinuit.mnexcm('SET PRIntout', arglist, 0, ierflag)

    gMinuit.SetPrintLevel(1)
    gMinuit.SetErrorDef(1.0)
    gMinuit.SetFCN(fcn)

    arglist_p = [2]
    arglist = array.array('d')
    arglist.fromlist(arglist_p)
    gMinuit.mnexcm('SET STRategy', arglist, 1, ierflag)

    arglist_p = [maxiter, 0.0000001]
    arglist = array.array('d')
    arglist.fromlist(arglist_p)
    gMinuit.mnexcm('MIGrad', arglist, 2, ierflag)
    gMinuit.SetMaxIterations(maxiter)
    
    # initialize fitting the variables
    vstart = [100.0] * partot
    # start alpha_z with 1
    vstart[partot-1] = 1.0
    vstart[partot-2] = 0

    step = [0.1] * partot
    upper = [100000]*partot
    lower = [0.1]*partot
    varname =[]

    if "ExpoPol" in bkgname:
        upper = [1000]*partot
        lower = [-1000]*partot

    if "ExpoBernstein" in bkgname:
        vstart[0] = -1
        upper[0] = 0
        lower[0] = -10

    for i in range( parfunction):
        varname.append("p"+str(i))

    varname.append("alpha_sig")
    varname.append("alpha_z")

    if doFloatZ:
        vstart[partot-1] = 1.0
        upper[partot-1] = 2
        lower[partot-1] = 0
        step[partot-1] = 0.01

    if isSpuriousFit:
        upper[partot-2] = 10.0
        lower[partot-2] = -10.0
        step[partot-2] = 0.1
        vstart[partot-2] = 1
    
    for i in range(partot):
        gMinuit.mnparm(i, varname[i], vstart[i], step[i], lower[i], upper[i], ierflag)

    if not isSpuriousFit:
        vstart[partot-2] =0
        gMinuit.FixParameter(partot-2)        
        
    if not doFloatZ:
        lower[partot-1] = 1
        upper[partot-1] = 1
        gMinuit.FixParameter(partot-1)

    if not isBkgPlusZFit:
        vstart[partot-1] = 0.0
        gMinuit.FixParameter(partot-1)

    
    # fitting procedure
    migradstat = gMinuit.Command('MIGrad ' + str(maxiter) + ' ' + str(0.001))
    improvestat = gMinuit.Command('IMProve ' + str(maxiter) + ' ' + str(0.01))

    for i in range(partot):
        arglist_p.append(i+1)

    arglist = array.array('d')
    arglist.fromlist(arglist_p)

    #gMinuit.mnmnos()

    # get fitting parameters
    fitval_p = [Double(0)] * partot
    fiterr_p = [Double(0)] * partot
    errup_p = [Double(0)] * partot
    errdown_p = [Double(0)] * partot
    eparab_p = [Double(0)] * partot
    gcc_p = [Double(0)] * partot

    fmin_p = [Double(0)]
    fedm_p = [Double(0)]
    errdef_p = [Double(0)]
    npari_p = Long(0)
    nparx_p = Long(0)
    istat_p = Long(0)

    fitval = array.array('d')
    fiterr = array.array('d')
    errup = array.array('d')
    errdown = array.array('d')
    eparab = array.array('d')
    gcc = array.array('d')
    
    for i in range(partot):
        gMinuit.GetParameter(i, fitval_p[i], fiterr_p[i])
        fitval.append(fitval_p[i])
        fiterr.append(fiterr_p[i])
        errup.append(errup_p[i])
        errdown.append(errdown_p[i])
        eparab.append(eparab_p[i])
        gcc.append(gcc_p[i])

    gMinuit.mnstat(fmin_p[0], fedm_p[0], errdef_p[0], npari_p, nparx_p, istat_p)
                        
    for p in range(bkgfunction.GetNumberFreeParameters()):
        bkgfunction.SetParameter(p, fitval[p])
        print "fit uncert",  fiterr_p[p]

    bkgfunction.SetChisquare(fmin_p[0])

    return  fitval[partot-1], fitval[partot-2]



def CalChi2(DataHist, function):

    chi2 = 0
    nbins = 0
    for i in range(DataHist.GetNbinsX()):
        databin = DataHist.GetBinContent(i+1)
        dataerr = DataHist.GetBinError(i+1)
        fval = function.Eval(DataHist.GetBinCenter(i+1))
          
        if dataerr !=0:
            nbins += 1
            chi2 += ((fval-databin)/dataerr)**2

            
    return chi2, nbins

#### main function used to run the minimization scheme
def signalfit(data_hist, signalfunction, signalname, process):

    binning = HistBinsToList(data_hist)
    data_x = HistToList(data_hist)
    data_error = HistErrorList(data_hist)

    parfunction = signalfunction.GetNumberFreeParameters()
    partot = signalfunction.GetNumberFreeParameters()
    print partot

    ### the fucntion used for TMinuit
    def fcn(npar, gin, f, par, ifag):
        L = 0
    
        # calculate likelihood, input par[0] is the N_B, par[1] is N_C, par[2] is N_L
        for ibin in range(len(binning)):
            #if (data_x[ibin] ==0):
            #    continue

            bincen = binning[ibin]
            
            mu_x = 0
            data = data_x[ibin]
            if data_error[ibin]==0:
                continue
            #if data<0.1:
            #    continue

            if signalname == "CrystalBall":

                if par[3]<0:
                    mu_x=0

                else:
                    t = (bincen-par[2])/(par[3]);
                    if (par[0] < 0):
                        t = -t

                    absAlpha = abs(par[0])

                    if (t >= -absAlpha):
                        mu_x = par[4]*exp(-0.5*t*t)

                    else:

                        nDivAlpha = par[1]/absAlpha
                        AA = exp(-0.5*absAlpha*absAlpha)
                        B = nDivAlpha-absAlpha
                        arg = nDivAlpha/(B-t)

                        mu_x = par[4]*(  arg**par[1])

            if signalname == "CrystalBallGaus":

                if par[3]<0:
                    mu_x=0

                else:
                    t = (bincen-par[2])/(par[3]);
                    if (par[0] < 0):
                        t = -t

                    absAlpha = abs(par[0])

                    if (t >= -absAlpha):
                        mu_x = par[4]*exp(-0.5*t*t + exp(- (bincen-par[5])**2/(2*par[6]**2) ))

                    else:

                        nDivAlpha = par[1]/absAlpha
                        AA = exp(-0.5*absAlpha*absAlpha)
                        B = nDivAlpha-absAlpha
                        arg = nDivAlpha/(B-t)

                        mu_x = par[4]*(  arg**par[1] + exp(- (bincen-par[5])**2/(2*par[6]**2) ))

                            
            #print mu_x, data, data_error[ibin]
            #L = L + mu_x - data*log(mu_x)
            L = L + ((mu_x - data)/data_error[ibin])**2
            
        f[0] = L

    # initialize the TMinuit object
    arglist_p = 10 * [0]

    arglist = array.array('d')
    arglist.fromlist(arglist_p)
    ierflag = Long(0)
    maxiter = 1000000000
    
    arglist_p= [1]
    gMinuit = TMinuit(partot)
    gMinuit.mnexcm('SET PRIntout', arglist, 0, ierflag)

    gMinuit.SetPrintLevel(1)
    gMinuit.SetErrorDef(1.0)
    gMinuit.SetFCN(fcn)

    arglist_p = [2]
    arglist = array.array('d')
    arglist.fromlist(arglist_p)
    gMinuit.mnexcm('SET STRategy', arglist, 1, ierflag)

    arglist_p = [maxiter, 0.0000001]
    arglist = array.array('d')
    arglist.fromlist(arglist_p)
    gMinuit.mnexcm('MIGrad', arglist, 2, ierflag)
    gMinuit.SetMaxIterations(maxiter)
    
    # initialize fitting the variables
    vstart = [125.0] * partot

    step = [0.1] * partot
    upper = [1000000]*partot
    lower = [-100]*partot
    varname =[]

    lower[3] = 0
    lower[4] = 0
    lower[1] = 0
    vstart[4] =data_hist.Integral()

    if process == "signal":
        vstart[2] = 125
        lower[2] = 110
        upper[2] = 140

        vstart[3] = 10
        lower[3] = 2
        upper[3] = 25

        if len(vstart)>5:
            vstart[5] = 125
            lower[5] = 110
            upper[5] = 140
            
            vstart[6] = 10
            lower[6] = 5
            upper[6] = 20

    if process == "z":
        vstart[2] = 90
        lower[2] = 70
        upper[2] = 110

        vstart[3] = 10
        lower[3] = 2
        upper[3] = 30

        if len(vstart)>5:
            vstart[5] = 90
            lower[5] = 70
            upper[5] = 110
            
            vstart[6] = 10
            lower[6] = 2
            upper[6] = 30


    for i in range( parfunction):
        varname.append("p"+str(i))

    for i in range(partot):
        gMinuit.mnparm(i, varname[i], vstart[i], step[i], lower[i], upper[i], ierflag)
    
    # fitting procedure
    migradstat = gMinuit.Command('MIGrad ' + str(maxiter) + ' ' + str(0.001))
    #improvestat = gMinuit.Command('IMProve ' + str(maxiter) + ' ' + str(0.01))

    for i in range(partot):
        arglist_p.append(i+1)

    arglist = array.array('d')
    arglist.fromlist(arglist_p)

    gMinuit.mnmnos()

    # get fitting parameters
    fitval_p = [Double(0)] * partot
    fiterr_p = [Double(0)] * partot
    errup_p = [Double(0)] * partot
    errdown_p = [Double(0)] * partot
    eparab_p = [Double(0)] * partot
    gcc_p = [Double(0)] * partot

    fmin_p = [Double(0)]
    fedm_p = [Double(0)]
    errdef_p = [Double(0)]
    npari_p = Long(0)
    nparx_p = Long(0)
    istat_p = Long(0)

    fitval = array.array('d')
    fiterr = array.array('d')
    errup = array.array('d')
    errdown = array.array('d')
    eparab = array.array('d')
    gcc = array.array('d')
    
    for i in range(partot):
        gMinuit.GetParameter(i, fitval_p[i], fiterr_p[i])
        fitval.append(fitval_p[i])
        fiterr.append(fiterr_p[i])
        errup.append(errup_p[i])
        errdown.append(errdown_p[i])
        eparab.append(eparab_p[i])
        gcc.append(gcc_p[i])

    gMinuit.mnstat(fmin_p[0], fedm_p[0], errdef_p[0], npari_p, nparx_p, istat_p)
                        
    for p in range(signalfunction.GetNumberFreeParameters()):
        signalfunction.SetParameter(p, fitval[p])
        print "fit uncert",  fiterr_p[p]

    signalfunction.SetChisquare(fmin_p[0])
    print fmin_p[0]
    return  fitval[partot-1], fitval[partot-2]



def RooFitSig(mbbarray, bdtarray, weightarray, TC_mass, binstart, binend):

    fitstart = 40
    fitend = 150

    mbbarray = range(200)
    bdtarray =  range(200)
    weightarray = range(200)

    mass = RooRealVar("X","m(bb)[GeV]",fitstart, fitend)
    BDT = RooRealVar("BDT","BDT",-1,100)
    weight = RooRealVar("weight","weight",-100,200)

    branchnames = ["X", "BDT", "weight"]

    dtype = np.dtype([ (branchnames[idx], np.float64) for idx in range(len(branchnames))])
    treearray = np.array([ (mbbarray[idx], bdtarray[idx], weightarray[idx]) for idx in range(len(mbbarray))], dtype)

    tree = rnp.array2tree(treearray)
    
    m0= RooRealVar("m0", "m0", TC_mass*1., TC_mass*1.-60., TC_mass*1.+60.)
    m02 = RooRealVar("m02", "m02", TC_mass*1., TC_mass*1.-60., TC_mass*1.+60.)
    alpha = RooRealVar("alpha", "alpha", 1.295, 1.0, 1.6)
    sigma2 = RooRealVar("sigma2", "sigma2", 35, 8., 100)
    n = RooRealVar("n", "n", 5,1,35)
    
    mean = RooRealVar("mean","mean of gaussian",750,0,6000)
    sigma = RooRealVar("sigma","width of gaussian",90,38,300)

    gauss= RooGaussian("gauss","gaussian PDF",mass,m0,sigma)   
    gauss2 =RooGaussian("gauss2","gaussian PDF",mass,m02,sigma2)   
    CBshape = RooCBShape("CBshape", "Crystal Ball PDF", mass, m0, sigma2, alpha, n)

    ##PDF normalization
    num1 = RooRealVar("num1","number of events",400,0,5000) 

    ##relative weight of 2 PDFs
    f = RooRealVar("f", "f", 0.95, 0.6, 1); 

    sigPdf = RooAddPdf("sigPdf", "Signal PDF", RooArgList(CBshape, gauss), RooArgList(f))
    extPdf = RooExtendPdf("extPdf", "extPdf", sigPdf, num1)
    data =  RooDataSet("data", "data", tree, RooArgSet(mass,BDT,weight), "BDT>0","weight");

    xframe = mass.frame()
    mass.setBins(20)
    data.plotOn(xframe) 
    extPdf.plotOn(xframe)#,Normalization(1.0,RooAbsReal.RelativeExpected),LineColor(1))

    hist = extPdf.createHistogram("X", fitend-fitstart)
    hist.SetAxisRange(binstart, binend)
    return deepcopy(hist)
