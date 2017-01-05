import ROOT
import os,sys

f = ROOT.TFile(sys.argv[1],'read')
h = f.Get("CutFlow")

for i in range(1, h.GetNbinsX()):
  print h.GetXaxis().GetBinLabel(i), '\t', h.GetBinContent(i)
