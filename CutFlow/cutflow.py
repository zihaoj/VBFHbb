from ROOT import *

file = TFile("Zbb_EWK_Reader_03_21_4cen_PostICHEP.root", "r")

presel= file.Get("CutFlow/PreselectionCutFlow")

for ibin in range(presel.GetNbinsX()):
    print presel.GetBinContent(ibin+1)

sel= file.Get("CutFlow/Nominal/CutsNoWeight")

for ibin in range(presel.GetNbinsX()):
    print sel.GetBinContent(ibin+1)



    
