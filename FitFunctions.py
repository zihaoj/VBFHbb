from ROOT import TF1, TH1D, TFile, gMinuit, Long, Double, TMinuit
import random as rand
from math import sqrt, exp, log
import array
from cStringIO import StringIO

fit_start = 80
fit_end = 300
fit_range = fit_end-fit_start

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



def BernsteinO2Lin_raw(x, par):
    return (par[0]*(1-(x[0]-fit_start)/fit_range)**2+ 2*par[1]*(1-(x[0]-fit_start)/fit_range)*((x[0]-fit_start)/fit_range) +  par[2]*((x[0]-fit_start)/fit_range)**2)*(par[3]*(x[0]-fit_start)/fit_range+par[4])

BernsteinO2Lin = TF1("BernsteinO2Lin", BernsteinO2Lin_raw, fit_start, fit_end, 5)


def BernsteinO3Lin_raw(x, par):
    return (par[0]*(1-((x[0]-fit_start)/fit_range))**3 + par[1]*(3*((x[0]-fit_start)/fit_range)*(1-((x[0]-fit_start)/fit_range))**2) + par[2]*(3*((x[0]-fit_start)/fit_range)**2*(1-((x[0]-fit_start)/fit_range))) + par[3]*((x[0]-fit_start)/fit_range)**3)*(par[4]*(x[0]-fit_start)/fit_range+par[5])

BernsteinO3Lin = TF1("BernsteinO3Lin", BernsteinO3Lin_raw, fit_start, fit_end, 6)


def BernsteinO4Lin_raw(x, par):
    return (par[0]*(1-((x[0]-fit_start)/fit_range))**4 + par[1]*(4*((x[0]-fit_start)/fit_range)*(1-((x[0]-fit_start)/fit_range))**3) + par[2]*(6*((x[0]-fit_start)/fit_range)**2*(1-((x[0]-fit_start)/fit_range))**2) + par[3]*(4*((x[0]-fit_start)/fit_range)**3*(1-((x[0]-fit_start)/fit_range))) + par[4]*((x[0]-fit_start)/fit_range)**4)*(par[5]*(x[0]-fit_start)/fit_range+par[6])

BernsteinO4Lin = TF1("BernsteinO4Lin", BernsteinO4Lin_raw, fit_start, fit_end, 7)



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


###################### Sum of Exponentials
def Expo_raw(x, par):
    try:
        return exp(par[0]+par[1]*(x[0]-fit_start)/fit_range)
    except OverflowError:
        return -1

Expo = TF1("Expo", Expo_raw, fit_start, fit_end, 2)

def Expo2_raw(x, par):
    return (exp(par[0]+par[1]*(x[0]-fit_start)/fit_range) + exp(par[2]+par[3]*(x[0]-fit_start)/fit_range))

Expo2 = TF1("Expo2", Expo2_raw, fit_start, fit_end, 4)

def Expo3_raw(x, par):
    return (exp(par[0]+par[1]*(x[0]-fit_start)/fit_range) + exp(par[2]+par[3]*(x[0]-fit_start)/fit_range) + exp(par[4]+par[5]*(x[0]-fit_start)/fit_range))

Expo3 = TF1("Expo3", Expo3_raw, fit_start, fit_end, 6)


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
                bkg = par[0]*(1-(bincen-fit_start)/fit_range)**2+ 2*par[1]*(1-(bincen-fit_start)/fit_range)*((bincen-fit_start)/fit_range) +  par[2]*((bincen-fit_start)/fit_range)**2

            if bkgname == "BernsteinO3":
                bkg = par[0]*(1-((bincen-fit_start)/fit_range))**3 + par[1]*(3*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**2) + par[2]*(3*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))) + par[3]*((bincen-fit_start)/fit_range)**3

            if bkgname == "BernsteinO4":
                bkg = par[0]*(1-((bincen-fit_start)/fit_range))**4 + par[1]*(4*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**3) + par[2]*(6*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))**2) + par[3]*(4*((bincen-fit_start)/fit_range)**3*(1-((bincen-fit_start)/fit_range))) + par[4]*((bincen-fit_start)/fit_range)**4

            if bkgname == "BernsteinO5":
                bkg = par[0]*(1-((bincen-fit_start)/fit_range))**5 + par[1]*(5*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**4) + par[2]*(10*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))**3) + par[3]*(10*((bincen-fit_start)/fit_range)**3*(1-((bincen-fit_start)/fit_range))**2) + par[4]*(5*((bincen-fit_start)/fit_range)**4*(1-((bincen-fit_start)/fit_range))) + par[5]*((bincen-fit_start)/fit_range)**5


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

            if bkgname == "Expo2":
                try:
                    bkg = exp(par[0]+par[1]*(bincen-fit_start)/fit_range) + exp(par[2]+par[3]*(bincen-fit_start)/fit_range)
                    #bkg = exp(par[0]+par[1]*bincen) + exp(par[2]+par[3]*bincen)
                except OverflowError:
                    bkg = 0

            if bkgname == "Expo3":
                try:
                    bkg = exp(par[0]+par[1]*(bincen-fit_start)/fit_range) + exp(par[2]+par[3]*(bincen-fit_start)/fit_range) + exp(par[4]+par[5]*(bincen-fit_start)/fit_range)
                except OverflowError:
                    bkg = 0


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
    upper = [10000000]*partot
    lower = [0.1]*partot
    varname =[]

    #bkg = exp(par[0]+par[1]*(bincen-fit_start)/fit_range) + exp(par[2]+par[3]*(bincen-fit_start)/fit_range)

    if "ExpoBernstein" in bkgname:
        vstart[0] = -10
        upper[0] = 0
        lower[0] = -500

    elif "Expo2" in bkgname:
        vstart[1] = -5
        vstart[3] = 0
        upper[1] = 0
        upper[3] = 10000
        lower[1] = -10000
        lower[3] = -10000

    elif "Expo3" in bkgname:
        vstart[1] = -100
        vstart[3] = -10
        vstart[5] = -5
        upper[1] = 0
        upper[3] = 500
        upper[5] = 500
        lower[1] = -500
        lower[3] = -500
        lower[5] = -500

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



#### main function used to run the minimization scheme
def bkgfit_2Region(data_hist_1, data_hist_2, bkgfunction, bkgname, z_hist_1 = None, z_hist_2 = None):
    isBkgPlusZFit = True
    isSpuriousFit = False

    binning = HistBinsToList(data_hist_1)
    data_x_1 = HistToList(data_hist_1)
    data_error_1 = HistErrorList(data_hist_1)
    data_x_2 = HistToList(data_hist_2)
    data_error_2 = HistErrorList(data_hist_2)

    z_x_1 = []
    z_x_2 = []

    z_x_1 = HistToList(z_hist_1)
    z_x_2 = HistToList(z_hist_2)

    parfunction = bkgfunction.GetNumberFreeParameters()
    partot = bkgfunction.GetNumberFreeParameters()  ## norm and linear

    ### the fucntion used for TMinuit
    def fcn(npar, gin, f, par, ifag):
        L = 0
    
        # calculate likelihood, input par[0] is the N_B, par[1] is N_C, par[2] is N_L
        for ibin in range(len(binning)):
            if (data_x_1[ibin] <0.5):
                continue
            if (data_x_2[ibin] <0.5):
                continue

            bincen = binning[ibin]
            
            bkg = 0
            data_1 = data_x_1[ibin]
            data_2 = data_x_2[ibin]

            if bkgname == "BernsteinO2Lin":
                bkg_1 = par[0]*(1-(bincen-fit_start)/fit_range)**2+ 2*par[1]*(1-(bincen-fit_start)/fit_range)*((bincen-fit_start)/fit_range) +  par[2]*((bincen-fit_start)/fit_range)**2
                bkg_2 = bkg_1*(par[3]*bincen+par[4])

            if bkgname == "BernsteinO3Lin":
                bkg_1 = par[0]*(1-((bincen-fit_start)/fit_range))**3 + par[1]*(3*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**2) + par[2]*(3*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))) + par[3]*((bincen-fit_start)/fit_range)**3
                bkg_2 = bkg_1*(par[4]*bincen+par[5])

            if bkgname == "BernsteinO4Lin":
                bkg_1 = par[0]*(1-((bincen-fit_start)/fit_range))**4 + par[1]*(4*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**3) + par[2]*(6*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))**2) + par[3]*(4*((bincen-fit_start)/fit_range)**3*(1-((bincen-fit_start)/fit_range))) + par[4]*((bincen-fit_start)/fit_range)**4
                bkg_2 = bkg_1*(par[5]*bincen+par[6])

            if bkgname == "BernsteinO5":
                bkg_1 = par[0]*(1-((bincen-fit_start)/fit_range))**5 + par[1]*(5*((bincen-fit_start)/fit_range)*(1-((bincen-fit_start)/fit_range))**4) + par[2]*(10*((bincen-fit_start)/fit_range)**2*(1-((bincen-fit_start)/fit_range))**3) + par[3]*(10*((bincen-fit_start)/fit_range)**3*(1-((bincen-fit_start)/fit_range))**2) + par[4]*(5*((bincen-fit_start)/fit_range)**4*(1-((bincen-fit_start)/fit_range))) + par[5]*((bincen-fit_start)/fit_range)**5 
                bkg_2 = bkg_1*(par[6]*bincen+par[7])


            mu_x_1 =  bkg_1
            mu_x_2 =  bkg_2

            data_1 = data_x_1[ibin]
            data_2 = data_x_2[ibin]

            mu_x_1 = mu_x_1 + (z_x_1[ibin])
            mu_x_2 = mu_x_2 + (z_x_2[ibin])

            L = L + ((mu_x_1 - data_1)/data_error_1[ibin])**2 +((mu_x_2 - data_2)/data_error_2[ibin])**2
                
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
    vstart = [0.0] * partot

    step = [0.1] * partot
    upper = [1000000]*partot
    lower = [-1000000]*partot
    varname =[]

    for i in range( partot):
        varname.append("p"+str(i))

    varname.append("alpha_sig")
    varname.append("alpha_z")

    for i in range(partot):
        gMinuit.mnparm(i, varname[i], vstart[i], step[i], lower[i], upper[i], ierflag)

    
    # fitting procedure
    migradstat = gMinuit.Command('MIGrad ' + str(maxiter) + ' ' + str(0.001))
    improvestat = gMinuit.Command('IMProve ' + str(maxiter) + ' ' + str(0.01))

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
                        
    for p in range(bkgfunction.GetNumberFreeParameters()):
        bkgfunction.SetParameter(p, fitval[p])

    bkgfunction.SetChisquare(fmin_p[0])

    for i in range(partot):
        varname[i], fitval[i]
    print "2 regions chi^2", fmin_p[0]
