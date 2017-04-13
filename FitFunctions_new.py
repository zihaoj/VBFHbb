from ROOT import TF1, TH1D, TFile, gMinuit, Long, Double, TMinuit
import random as rand
from math import sqrt, exp, log
import array
from cStringIO import StringIO

fit_start = 80
fit_end = 200
fit_range = fit_end-fit_start


def HistToList (hist):
    hist_list = []
    for i in range(hist.GetNbinsX()):
        hist_list.append(hist.GetBinContent(i+1))
    return hist_list

def HistErrorList (hist):
    hist_list = []
    for i in range(hist.GetNbinsX()):
        err = hist.GetBinError(i+1)
        hist_list.append(err)
    return hist_list

def HistBinsToList (hist):
    bins_list = []
    for i in range(hist.GetNbinsX()):
        bins_list.append(hist.GetXaxis().GetBinCenter(i+1))
    return bins_list


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

CrystalBall = TF1("CrystalBall", CrystalBall_raw, fit_start, fit_end, 5)


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


CrystalBallGaus = TF1("CrystalBallGaus", CrystalBallGaus_raw, fit_start, fit_end, 7)


#### main function used to run the minimization scheme
def signalfit(data_hist, signalfunction, signalname):

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

    vstart[3] = data_hist.GetRMS()
    lower[3] = data_hist.GetRMS()/3.0
    upper[3] = data_hist.GetRMS()*3.0
    vstart[3] = 125
    lower[3] = 5
    upper[3] = 20

    vstart[4] =data_hist.Integral()

    vstart[2] = data_hist.GetMean()
    lower[2] = data_hist.GetMean()-10
    upper[2] = data_hist.GetMean()+10

    vstart[2] = 125
    lower[2] = 110
    upper[2] = 130

    if len(vstart)>5:
        vstart[5] = data_hist.GetMean()
        lower[5] = data_hist.GetMean()-10
        upper[5] = data_hist.GetMean()+10
        vstart[5] = 125
        lower[5] = 110
        upper[5] = 150

        #vstart[6] = data_hist.GetRMS()
        #lower[6] = data_hist.GetRMS()/3.0
        #upper[6] = data_hist.GetRMS()*3.0

        vstart[6] = 10
        lower[6] = 5
        upper[6] = 20


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
