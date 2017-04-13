from ROOT import gDirectory, TLorentzVector, TVector3, TH1D, TH2D, TFile, RooStats, TMath
from HistoLab import HistoTool, CopyHist, StackHists, GetPlot, ScaleHist, Project2D, CopyHists, project1D, NormHists, AddHistPairs, ResetBins, GetHistMedian, GetHistIQR, GetRatios, GetAveRatios, GetGausMean, CalSysVariation, CalSysDiff
import glob
import numpy as np

DrawTool = HistoTool()
DrawTool.lumi = "0"
DrawTool.sqrtS = "13"
DrawTool.SetCompData(False)
DrawTool.OutPlotDir = "Plots/"

output = TFile("MCToy.root", "recreate")
output.cd()

def PlotMuToy(mufile, muinject):

    Hist_mu = TH1D("Mu_"+str(muinject), "Mu_"+str(muinject), 20, muinject-7, muinject+7)
    Hist_pull = TH1D("Pull_"+str(muinject), "Pull_"+str(muinject), 20, -5, 5)

    muhats = []
    txt = open(mufile, 'r')
    lines = txt.readlines()
    for line in lines:
        linebreak = line.split(" ")
        mu = float( linebreak[1].lstrip("Mu") )
        err = float( linebreak[2].lstrip("Err") )
        muhats.append(mu)
        Hist_mu.Fill(mu)
        Hist_pull.Fill( (mu - muinject)/err/1.1  )

        muhats.append(mu)
    
    print sum(muhats)/len(muhats)
    print np.std(muhats)


    DrawTool.DrawHists("Mu_"+str(muinject),  ["#mu", "A.U."], [Hist_mu], ["Injected #mu = "+str(muinject)])

    Hist_pull.Fit("gaus")
    gausfit = Hist_pull.GetFunction("gaus")
    mean = round(gausfit.GetParameter(1), 3)
    mean_error = round(gausfit.GetParError(1), 3)
    width = round(gausfit.GetParameter(2), 3)
    width_error = round(gausfit.GetParError(2), 3)

    DrawTool.AddTexts(0.18, 0.7, 1, 0.04, "Mean = {!s}#pm{!s}".format(str(mean), str(mean_error)))
    DrawTool.AddTexts(0.18, 0.66, 1, 0.04, "Width = {!s}#pm{!s}".format(str(width), str(width_error)))

    DrawTool.DrawHists("Pull_"+str(muinject),  ["(#mu_{meas}-#mu_{inj})/#sigma", "A.U."], [Hist_pull], ["Injected #mu = "+str(muinject)])
    DrawTool.ClearTexts()

PlotMuToy("FitToy/toy_mc_batch_mu1.txt", 1)
PlotMuToy("FitToy/toy_mc_batch_mu2.txt", 2)
PlotMuToy("FitToy/toy_mc_batch_mu5.txt", 5)
PlotMuToy("FitToy/toy_mc_batch_mu05.txt", 0.5)

    
output.Close()
