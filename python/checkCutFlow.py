import ROOT
from myutils import util_funcs
import os,sys

nProcess = 5
isData = False

tree = ROOT.TChain("tree")

#for line in open('anafiles_tmp_1.txt').readlines():
#  tree.Add(line.replace('\n',''))

#tree.Add('/eos/uscms/store/user/lpcphys/noreplica/duong/VHbbTuples/Data_2016_merge/SingleMuon/ZC2016_11_Run2016_promptReco_3_SingleMuon__Run2016C-PromptReco-v2.root')
#tree.Add('/eos/uscms/store/user/lpcphys/noreplica/duong/VHbbTuples/Data_2016_merge/SingleMuon/ZC2016_11_Run2016_promptReco_3_SingleMuon__Run2016D-PromptReco-v2.root')
#tree.Add('/eos/uscms/store/user/yokugawa/Test/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root')
#tree.Add('/eos/uscms/store/user/lpcphys/noreplica/duong/VHbbTuples/MC_2016_merge/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0.root')
#for l in open('DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt').readlines():
#  tree.Add(l.replace('\n',''))
tree.Add('/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V24/Zmm_V24/syst//ZC2016_11_Run2016_promptReco_3_SingleMuon__Run2016C-PromptReco-v2.root')
tree.Add('/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V24/Zmm_V24/syst//ZC2016_11_Run2016_promptReco_3_SingleMuon__Run2016D-PromptReco-v2.root')


nentries = tree.GetEntries()
#nentries = 100000

nTmp = int(nentries/(nProcess))

evt_begin = []
evt_end = []
for i in range(0, nProcess):
  evt_begin.append(i*nTmp)
  evt_end.append((i+1)*nTmp-1)

evt_end[-1] = nentries - 1

print evt_begin
print evt_end



print '>>>>>>>>Total number of events are: ', nentries


def processEvents(tree, outputname, evtB, evtE, evtOut):

    h_cutFlow = ROOT.TH1D("CutFlow", "", 10, 0, 10)
    h_cutFlow.GetXaxis().SetBinLabel(1, "All")
    h_cutFlow.GetXaxis().SetBinLabel(2, "JSON")
    h_cutFlow.GetXaxis().SetBinLabel(3, "Zlepton")
    h_cutFlow.GetXaxis().SetBinLabel(4, "Zboson")
    h_cutFlow.GetXaxis().SetBinLabel(5, "Z+jet")
    h_cutFlow.GetXaxis().SetBinLabel(6, "MET")
    h_cutFlow.GetXaxis().SetBinLabel(7, "Z+bjet (pT sorted)")
    h_cutFlow.GetXaxis().SetBinLabel(8, "Z+bjet_1 (CSV sorted)")

    evtOut = []

    for i in range(evtB,evtE+1):

      if i/float(100000) == i/100000: print "Processing: ", i, " event"

      tree.GetEntry(i)
     
      #@@@@@@@@@@@lepton selection@@@@@@@@@@@@@@@@@ 
      h_cutFlow.Fill(0.5)
      
      #json (temp no trigger)
      if isData and not tree.json: continue
      h_cutFlow.Fill(1.5)

      if not tree.Vtype == 0 or tree.nvLeptons < 2: continue
      if tree.vLeptons_pt[0] < 20. or tree.vLeptons_pt[1] < 20. or abs(tree.vLeptons_eta[0]) > 2.4 or abs(tree.vLeptons_eta[1]) > 2.4: continue
      if (tree.vLeptons_relIso04[0] > 0.12 or tree.vLeptons_relIso04[1] > 0.12): continue
      if (tree.vLeptons_charge[0]*tree.vLeptons_charge[1] >= 0): continue
      h_cutFlow.Fill(2.5)

     
      if tree.V_mass < 70 or tree.V_mass > 110: continue
      h_cutFlow.Fill(3.5)
      
      nSelJet = 0
      sortJet_idx = []
      sortCSVjet_idx = []
      for iJ in range(0,tree.nJet):
        if tree.Jet_pt[iJ] > 30 and abs(tree.Jet_eta[iJ]) < 2.4:
          sortJet_idx.append(iJ)
          sortCSVjet_idx.append(iJ)
          nSelJet = nSelJet + 1

      if nSelJet < 1: continue
      h_cutFlow.Fill(4.5)
      util_funcs.sortPt(sortJet_idx, tree.Jet_pt) #sorted jet according to pT
      util_funcs.sortPtFatjet(sortCSVjet_idx,tree.Jet_btagCSV,nSelJet) #sorted jet according to CSV
      
      if tree.met_pt > 40: continue
      h_cutFlow.Fill(5.5)
      
      nSelJet_tag = 0
      for i in range(0,len(sortJet_idx)):
        iJ = sortJet_idx[i]
        if tree.Jet_btagCSV[iJ] > 0.935 and tree.Jet_vtxMass[iJ] > 0:
          nSelJet_tag = nSelJet_tag + 1

      nSelJet_tag_1 = 0
      if len(sortCSVjet_idx) > 0:
        iJ = sortCSVjet_idx[0]
        if tree.Jet_btagCSV[iJ] > 0.935 and tree.Jet_vtxMass[iJ] > 0:
          nSelJet_tag_1 = nSelJet_tag_1 + 1


      if nSelJet_tag >= 1: h_cutFlow.Fill(6.5)
      if nSelJet_tag_1 >= 1: h_cutFlow.Fill(7.5)


      #valid muon: Vtype = 0, two lepton passing pt and eta and isolation 
      #valid Z-boson: with mass
      #Z+jets: with one jet passing pT and eta
      #MET cut
      #jet with CSV and SVT

      #if tree.Vtype !=0: continue
      #sort jet by pt, if there atleas one jet passing pT and eta



       
    fOut = ROOT.TFile(outputname, "RECREATE")
    h_cutFlow.Write()
    fOut.cd()
    fOut.Close()
    
    #evtOut.append(nPassEvent)
    
    return evtOut
 
def processEventsOneInput(inputs):
  return processEvents(*inputs)

#evtOut = []
#processEvents(tree, "out1.root", 0, 1000, evtOut)

#print evtOut

inputs = []
for i in range(0, nProcess):
  outName = "out_" + str(i) + '.root'
  evtOut = []
  inputs.append((tree,outName,evt_begin[i],evt_end[i],evtOut))

print inputs

from multiprocessing import Pool
p = Pool(nProcess)
outputs = p.map(processEventsOneInput,inputs)

print outputs


if __name__ == '__main__':
   rep = ''
   while not rep in [ 'q', 'Q' ]:
     rep = raw_input( 'enter "q" to quit: ' )
     if 1 < len(rep):
       rep = rep[0]
  
