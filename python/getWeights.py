from ROOT import *

prefix = "/scratch/sdonato/VHbbRun2/V21_2/CMSSW_7_1_5/src/Xbb/env/ZvvHighPt_V21_"

def getWeight(fileInc, fileB, region):
    f = TFile(prefix+fileInc+".root")
    tree = f.Get("tree")
    countInc    = 1.* tree.Draw("",region)
    f.Close()

    f = TFile(prefix+fileB+".root")
    tree = f.Get("tree")
    countB      = 1.* tree.Draw("",region)
    f.Close()

    weight = countInc/(countB+countInc)
    return weight

WjetsHT100       = "WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
WjetsHT200       = "WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
WjetsHT400       = "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
WjetsHT600       = "WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"

ZjetsHT100       = "ZJetsToNuNu_HT-100To200_13TeV-madgraph"
ZjetsHT200       = "ZJetsToNuNu_HT-200To400_13TeV-madgraph"
ZjetsHT400       = "ZJetsToNuNu_HT-400To600_13TeV-madgraph"
ZjetsHT600       = "ZJetsToNuNu_HT-600ToInf_13TeV-madgraph"

ZLLjetsHT0     = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
ZLLjetsHT100     = "DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
ZLLjetsHT200     = "DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
ZLLjetsHT400     = "DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
ZLLjetsHT600     = "DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"

ZLLBjets        = "DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"

ZBjets          = "DYBJetsToNuNu_Zpt-40toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"

WBjets          = "WBJetsToLNu_Wpt-40toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
WjetsBgen       = "WJetsToLNu_BGenFilter_Wpt-40toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"

DYBJets             = "(lheNb>0 && lheV_pt>40)"
DYLightJets         = "!(lheNb>0 && lheV_pt>40)"

WBJets             = " (lheNb>0  && lheV_pt>40)"
WJetsBGen           = " (lheNb==0 && nGenStatus2bHad>0  && lheV_pt>40)"
WLightJets          = "((lheNb==0 && nGenStatus2bHad==0 && lheV_pt>40) || (lheV_pt<40))"

HT0            = "(lheHT<100)"
HT100          = "(lheHT>100&&lheHT<200)"
HT200          = "(lheHT>200&&lheHT<400)"
HT400          = "(lheHT>400&&lheHT<600)"
HT600          = "(lheHT>600)"

print "weightWBjetsHT100=\t%.2f" %getWeight(WjetsHT100,   WBjets, HT100+"&&"+WBJets)
print "weightWBjetsHT200=\t%.2f" %getWeight(WjetsHT200,   WBjets, HT200+"&&"+WBJets)
print "weightWBjetsHT400=\t%.2f" %getWeight(WjetsHT400,   WBjets, HT400+"&&"+WBJets)
print "weightWBjetsHT600=\t%.2f" %getWeight(WjetsHT600,   WBjets, HT600+"&&"+WBJets)
print ""
print "weightWjetsBgenHT100=\t%.2f" %getWeight(WjetsHT100,   WjetsBgen, HT100+"&&"+WJetsBGen)
print "weightWjetsBgenHT200=\t%.2f" %getWeight(WjetsHT200,   WjetsBgen, HT200+"&&"+WJetsBGen)
print "weightWjetsBgenHT400=\t%.2f" %getWeight(WjetsHT400,   WjetsBgen, HT400+"&&"+WJetsBGen)
print "weightWjetsBgenHT600=\t%.2f" %getWeight(WjetsHT600,   WjetsBgen, HT600+"&&"+WJetsBGen)
print ""
print "weightZBjetsHT100=\t%.2f" %getWeight(ZjetsHT100,   ZBjets, HT100+"&&"+DYBJets)
print "weightZBjetsHT200=\t%.2f" %getWeight(ZjetsHT200,   ZBjets, HT200+"&&"+DYBJets)
print "weightZBjetsHT400=\t%.2f" %getWeight(ZjetsHT400,   ZBjets, HT400+"&&"+DYBJets)
print "weightZBjetsHT600=\t%.2f" %getWeight(ZjetsHT600,   ZBjets, HT600+"&&"+DYBJets)
print ""
print "weightZLLBjetsHT0=\t%.2f" %getWeight(ZLLjetsHT0,   ZLLBjets, HT0+"&&"+DYBJets)
print "weightZLLBjetsHT100=\t%.2f" %getWeight(ZLLjetsHT100,   ZLLBjets, HT100+"&&"+DYBJets)
print "weightZLLBjetsHT200=\t%.2f" %getWeight(ZLLjetsHT200,   ZLLBjets, HT200+"&&"+DYBJets)
print "weightZLLBjetsHT400=\t%.2f" %getWeight(ZLLjetsHT400,   ZLLBjets, HT400+"&&"+DYBJets)
print "weightZLLBjetsHT600=\t%.2f" %getWeight(ZLLjetsHT600,   ZLLBjets, HT600+"&&"+DYBJets)


