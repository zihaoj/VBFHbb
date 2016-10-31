from ROOT import *

#file = TFile("cutflow_VBF_Reader_10_28_4cen.root", "r")
#file = TFile("cutflow_ggF_Reader_10_28_4cen.root", "r")
#file = TFile("cutflow_Zbb_QCD_Reader_10_28_4cen.root", "r")
#file = TFile("cutflow_Zbb_EWK_Reader_10_28_4cen.root", "r")

#file = TFile("cutflow_VBF_Reader_10_28_2cen.root", "r")
#file = TFile("cutflow_ggF_Reader_10_28_2cen.root", "r")
#file = TFile("cutflow_Zbb_QCD_Reader_10_28_2cen.root", "r")
file = TFile("cutflow_Zbb_EWK_Reader_10_28_2cen.root", "r")

presel= file.Get("CutFlow/PreselectionCutFlow")

for ibin in range(presel.GetNbinsX()):
    print presel.GetBinContent(ibin+1)

sel= file.Get("CutFlow/Nominal/CutsNoWeight")

for ibin in range(presel.GetNbinsX()):
    print sel.GetBinContent(ibin+1)



    
