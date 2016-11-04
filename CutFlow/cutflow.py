from ROOT import *

file = TFile("cutflow_data_2tag_2cen.root", "r")

presel= file.Get("CutFlow/PreselectionCutFlow")

for ibin in range(presel.GetNbinsX()):
    print presel.GetBinContent(ibin+1)

sel= file.Get("CutFlow/Nominal/CutsNoWeight")

for ibin in range(presel.GetNbinsX()):
    print sel.GetBinContent(ibin+1)



    
