from ROOT import *

file = TFile("data_cf_03_17.root", "r")

presel= file.Get("CutFlow/PreselectionCutFlow")

for ibin in range(presel.GetNbinsX()):
    print presel.GetBinContent(ibin+1)

sel= file.Get("CutFlow/Nominal/CutsNoWeight")

for ibin in range(presel.GetNbinsX()):
    print sel.GetBinContent(ibin+1)



    
