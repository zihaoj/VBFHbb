from ROOT import *

orig_file = TFile("data_2tag_09_16.root", "read")
tree = orig_file.Get("Nominal")


new_file = TFile("data_2tag_09_16_copy.root", "recreate")
new_file.cd()
newtree = tree.CopyTree("mBB<100000 || mBB>140000")

newtree.Write()
new_file.Close()

