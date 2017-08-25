import os,sys
from ROOT import *
import ROOT
import array
import math

from myutils import util_funcs 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#btag SF
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
gSystem.Load('libCondFormatsBTauObjects') 
gSystem.Load('libCondToolsBTau')
bTagCalib = BTagCalibration('CSVv2', 'Data/CSVv2_Moriond17_B_H.csv')
v_sys = getattr(ROOT, 'vector<string>')()
#v_sys = ROOT.vector<string>()
v_sys.push_back('up')
v_sys.push_back('down')

# make a reader instance and load the sf data
btagSF_reader_M = BTagCalibrationReader(
    1,              # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
    "central",      # central systematic type
    v_sys,          # vector of other sys. types
)   
 
btagSF_reader_M.load(
    bTagCalib, 
    0,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

btagSF_reader_M.load(
    bTagCalib, 
    1,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

btagSF_reader_M.load(
    bTagCalib, 
    2,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "incl"      # measurement type
)


# make a reader instance and load the sf data
btagSF_reader_T = BTagCalibrationReader(
    2,              # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
    "central",      # central systematic type
    v_sys,          # vector of other sys. types
)    
btagSF_reader_T.load(
    bTagCalib, 
    0,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

btagSF_reader_T.load(
    bTagCalib, 
    1,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

btagSF_reader_T.load(
    bTagCalib, 
    2,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG  
    "incl"      # measurement type
)

#################################
#ctagSF
#################################
cTagCalib = BTagCalibration('cTag', 'Data/ctagger_Moriond17_B_H.csv')

# make a reader instance and load the sf data
ctagSF_reader_M = BTagCalibrationReader(
    1,              # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
    "central",      # central systematic type
    v_sys,          # vector of other sys. types
)   
 
ctagSF_reader_M.load(
    cTagCalib, 
    0,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "TnP"      # measurement type
)

ctagSF_reader_M.load(
    cTagCalib, 
    1,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

ctagSF_reader_M.load(
    cTagCalib, 
    2,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "incl"      # measurement type
)


# make a reader instance and load the sf data
ctagSF_reader_T = BTagCalibrationReader(
    2,              # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
    "central",      # central systematic type
    v_sys,          # vector of other sys. types
)

ctagSF_reader_T.load(
    cTagCalib, 
    0,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "TnP"      # measurement type
)

ctagSF_reader_T.load(
    cTagCalib, 
    1,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

ctagSF_reader_T.load(
    cTagCalib, 
    2,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG  
    "incl"      # measurement type
)

def fillHist(iJet, h, flav,tagName='CSVM'):
  if tr.Jet_pt[iJet] > 30 and abs(tr.Jet_eta[iJet]) < 2.4:
    if tr.Jet_hadronFlavour[iJet] == flav:
      h[0].Fill(tr.Jet_pt[iJet])
      passTag = False
      if tagName=='CSVM' and tr.Jet_btagCSV[iJet] > 0.8484:
        passTag = True
      if tagName=='CtagT' and tr.Jet_ctagVsL[iJet] > 0.69 and tr.Jet_ctagVsB[iJet] > -0.45:
        passTag = True
      if passTag: 
        h[1].Fill(tr.Jet_pt[iJet])
        #if tr.Jet_vtxMassCorr_IVF[iJet] > 0 and tr.Jet_vtxCat_IVF[iJet] == 0:
        if tr.Jet_vtxMass[iJet] > 0: 
          h[2].Fill(tr.Jet_pt[iJet])

def fillHist_8TeV(iJet, h, flav):
  if tr.allJet_pt[iJet] > 30 and abs(tr.allJet_eta[iJet]) < 2.4:
    passFlav = False
    if flav == 'bjet' and abs(tr.allJet_flavour[iJet]) == 5:
      passFlav = True
    if flav == 'cjet' and abs(tr.allJet_flavour[iJet]) == 4:
      passFlav = True
    if flav == 'ljet' and abs(tr.allJet_flavour[iJet]) != 4 and abs(tr.allJet_flavour[iJet]) != 5:
      passFlav = True
    if passFlav:
      h[0].Fill(tr.allJet_pt[iJet])
      if tr.allJet_csv[iJet] > 0.898:
        h[1].Fill(tr.allJet_pt[iJet])
        if tr.allJet_vtxMass[iJet] > 0: 
          h[2].Fill(tr.allJet_pt[iJet])


def calEff(h):
  hEffs = [h[1].Clone(h[0].GetName().replace('_all','') + '_eff_tag'),h[2].Clone(h[0].GetName().replace('_all','') + '_eff_tag_svt'),TH1D(h[0].GetName().replace('_all','') + '_eff_int','',2,0,1)] #int = average efficiency, bin 1 for tag, bin 2 for tag+SVT
  hEffs[0].Divide(h[0])
  hEffs[1].Divide(h[0])
  nEvt = [0,0,0] #deno, num 1 (for tag), num 2 (for tag+SVT)
  nEvt_err = [0,0,0]
  for i in range(3):
    nEvt[i] = h[i].Integral()
    for iBin in range(1,h[i].GetNbinsX()+1):
      nEvt_err[i] += h[i].GetBinError(iBin)*h[i].GetBinError(iBin)
  
  for i in range(1,3):
    r = 0
    err = 0
    if nEvt[0] != 0:
      r = nEvt[i]/nEvt[0]
      err = math.sqrt(nEvt_err[i]/(nEvt[i]*nEvt[i]) + nEvt_err[0]/(nEvt[0]*nEvt[0]))
    hEffs[2].SetBinContent(i, r)
    hEffs[2].SetBinError(i, r*err)

  return hEffs

#############################################################
#Main program
#############################################################
#tagName = 'CtagT' #CSVM
tagName = 'CSVM' 

f = TFile.Open('/uscms/home/duong/Scratch/Output_ZplusC_V25/Zll_inc_V25/syst_fromEOS/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root','read')

fOut = TFile.Open('Test/eff_13TeV_CSVM.root','recreate')

tr = f.Get('tree')

nEntries = tr.GetEntries()
#nEntries = 500000
debug = False 

xDiv = [30, 50, 70, 100, 140, 200, 300, 600, 1000]
if tagName == 'CtagT':
  xDiv = [30, 35, 50, 70, 200, 1000] 

#per jet efficiency
hB = [TH1D('hB_pt_all','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hB_pt_tag','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hB_pt_tag_svt','',len(xDiv)-1,array.array('f',xDiv))]
hC = [TH1D('hC_pt_all','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hC_pt_tag','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hC_pt_tag_svt','',len(xDiv)-1,array.array('f',xDiv))]
hL = [TH1D('hL_pt_all','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hL_pt_tag','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hL_pt_tag_svt','',len(xDiv)-1,array.array('f',xDiv))]

#per event efficiency
#central
hB_evt = [TH1D('hB_pt_evt_all','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hB_pt_evt_tag','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hB_pt_evt_tag_svt','',len(xDiv)-1,array.array('f',xDiv))]
hC_evt = [TH1D('hC_pt_evt_all','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hC_pt_evt_tag','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hC_pt_evt_tag_svt','',len(xDiv)-1,array.array('f',xDiv))]
hL_evt = [TH1D('hL_pt_evt_all','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hL_pt_evt_tag','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hL_pt_evt_tag_svt','',len(xDiv)-1,array.array('f',xDiv))]
#SF up and down
hB_evt_SFu = [TH1D('hB_pt_evt_all_SFu','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hB_pt_evt_tag_SFu','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hB_pt_evt_tag_svt_SFu','',len(xDiv)-1,array.array('f',xDiv))]
hB_evt_SFd = [TH1D('hB_pt_evt_all_SFd','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hB_pt_evt_tag_SFd','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hB_pt_evt_tag_svt_SFd','',len(xDiv)-1,array.array('f',xDiv))]
hC_evt_SFu = [TH1D('hC_pt_evt_all_SFu','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hC_pt_evt_tag_SFu','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hC_pt_evt_tag_svt_SFu','',len(xDiv)-1,array.array('f',xDiv))]
hC_evt_SFd = [TH1D('hC_pt_evt_all_SFd','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hC_pt_evt_tag_SFd','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hC_pt_evt_tag_svt_SFd','',len(xDiv)-1,array.array('f',xDiv))]
hL_evt_SFu = [TH1D('hL_pt_evt_all_SFu','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hL_pt_evt_tag_SFu','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hL_pt_evt_tag_svt_SFu','',len(xDiv)-1,array.array('f',xDiv))]
hL_evt_SFd = [TH1D('hL_pt_evt_all_SFd','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hL_pt_evt_tag_SFd','',len(xDiv)-1,array.array('f',xDiv)),TH1D('hL_pt_evt_tag_svt_SFd','',len(xDiv)-1,array.array('f',xDiv))]

for i in range(len(hB)):
  hB[i].Sumw2()
  hC[i].Sumw2()
  hL[i].Sumw2()
  hB_evt[i].Sumw2()
  hC_evt[i].Sumw2()
  hL_evt[i].Sumw2()

for iEvt in range(nEntries):
  tr.GetEntry(iEvt)
  if iEvt/10000 == iEvt/10000.: print 'Processing ', iEvt, ' in total ', nEntries, ', ', 100.*iEvt/nEntries, '%'
  if tr.Vtype != 0 and tr.Vtype != 1: continue
  if tr.vLeptons_pt[0] < 20 or tr.vLeptons_pt[1] < 20: continue
  if abs(tr.vLeptons_eta[0]) > 2.4 or abs(tr.vLeptons_eta[1]) > 2.4: continue
  for iJet in range(tr.nJet):
    fillHist(iJet, hB, 5, tagName)
    fillHist(iJet, hC, 4, tagName)
    fillHist(iJet, hL, 0, tagName)

  #calculate event weights
  selJet_idx = []
  for iJet in range(tr.nJet):
    if tr.Jet_pt[iJet] > 30 and abs(tr.Jet_eta[iJet]) < 2.4: selJet_idx.append(iJet)
  bcTagWeights = util_funcs.calBtagWeight(selJet_idx, tr, btagSF_reader_T, btagSF_reader_M, ctagSF_reader_T, ctagSF_reader_M)
  ws = [1,1,1]
  for i in range(3):
    if tagName == 'CtagT':
      ws[i] = bcTagWeights[2][i]
    if tagName == 'CSVM':
      ws[i] = bcTagWeights[1][i]
  if debug:
    print ws

  #loop over jets to see if it is Z+b, Z+c or Z+l events, also save the leading b-jet, c-jet, l-jet
  #check if events pass tagging and SVT requirement
  isZb = False
  iBj = -1
  isZc = False
  iCj = -1
  isZl = False
  iLj = -1
  isTagged = False
  isTagged_SVT = False
  for iJet in range(tr.nJet):
    if tr.Jet_pt[iJet] > 30 and abs(tr.Jet_eta[iJet]) < 2.4:
      if tr.Jet_hadronFlavour[iJet] == 5:
        isZb = True
        if iBj == -1: iBj = iJet
      if tr.Jet_hadronFlavour[iJet] == 4:
        isZc = True
        if iCj == -1: iCj = iJet
      if tr.Jet_hadronFlavour[iJet] == 0:
        isZl = True
        if iLj == -1: iLj = iJet
      if tagName == 'CtagT' and tr.Jet_ctagVsL[iJet] > 0.69 and tr.Jet_ctagVsB[iJet] > -0.45:
        isTagged = True
        if tr.Jet_vtxMass[iJet] > 0: isTagged_SVT = True
      if tagName == 'CSVM' and tr.Jet_btagCSV[iJet] > 0.8484:
        isTagged = True
        if tr.Jet_vtxMass[iJet] > 0: isTagged_SVT = True
 
  if isZb:
    hB_evt[0].Fill(tr.Jet_pt[iBj])
    hB_evt_SFu[0].Fill(tr.Jet_pt[iBj])
    hB_evt_SFd[0].Fill(tr.Jet_pt[iBj])
    if isTagged:
      hB_evt[1].Fill(tr.Jet_pt[iBj],ws[0])
      hB_evt_SFu[1].Fill(tr.Jet_pt[iBj],ws[1])
      hB_evt_SFd[1].Fill(tr.Jet_pt[iBj],ws[2])
    if isTagged_SVT:
      hB_evt[2].Fill(tr.Jet_pt[iBj],ws[0])
      hB_evt_SFu[2].Fill(tr.Jet_pt[iBj],ws[1])
      hB_evt_SFd[2].Fill(tr.Jet_pt[iBj],ws[2])
      
  if isZc and not isZb:
    hC_evt[0].Fill(tr.Jet_pt[iCj])
    hC_evt_SFu[0].Fill(tr.Jet_pt[iCj])
    hC_evt_SFd[0].Fill(tr.Jet_pt[iCj])
    if isTagged:
      hC_evt[1].Fill(tr.Jet_pt[iCj],ws[0])
      hC_evt_SFu[1].Fill(tr.Jet_pt[iCj],ws[1])
      hC_evt_SFd[1].Fill(tr.Jet_pt[iCj],ws[2])
    if isTagged_SVT:
      hC_evt[2].Fill(tr.Jet_pt[iCj],ws[0])
      hC_evt_SFu[2].Fill(tr.Jet_pt[iCj],ws[1])
      hC_evt_SFd[2].Fill(tr.Jet_pt[iCj],ws[2])

  if isZl and not isZb and not isZc:
    hL_evt[0].Fill(tr.Jet_pt[iLj])
    hL_evt_SFu[0].Fill(tr.Jet_pt[iLj])
    hL_evt_SFd[0].Fill(tr.Jet_pt[iLj])
    if isTagged:
      hL_evt[1].Fill(tr.Jet_pt[iLj],ws[0])
      hL_evt_SFu[1].Fill(tr.Jet_pt[iLj],ws[1])
      hL_evt_SFd[1].Fill(tr.Jet_pt[iLj],ws[2])
    if isTagged_SVT:
      hL_evt[2].Fill(tr.Jet_pt[iLj],ws[0])
      hL_evt_SFu[2].Fill(tr.Jet_pt[iLj],ws[1])
      hL_evt_SFd[2].Fill(tr.Jet_pt[iLj],ws[2])

hEffs_B = calEff(hB)
hEffs_C = calEff(hC)
hEffs_L = calEff(hL)

hEffs_B_evt = calEff(hB_evt)
hEffs_C_evt = calEff(hC_evt)
hEffs_L_evt = calEff(hL_evt)

hEffs_B_evt_SFu = calEff(hB_evt_SFu)
hEffs_C_evt_SFu = calEff(hC_evt_SFu)
hEffs_L_evt_SFu = calEff(hL_evt_SFu)

hEffs_B_evt_SFd = calEff(hB_evt_SFd)
hEffs_C_evt_SFd = calEff(hC_evt_SFd)
hEffs_L_evt_SFd = calEff(hL_evt_SFd)

fOut.cd()
for i in range(len(hB)):
  hB[i].Write()
  hC[i].Write()
  hL[i].Write()
  hB_evt[i].Write()
  hC_evt[i].Write()
  hL_evt[i].Write()
  hB_evt_SFu[i].Write()
  hC_evt_SFu[i].Write()
  hL_evt_SFu[i].Write()
  hB_evt_SFd[i].Write()
  hC_evt_SFd[i].Write()
  hL_evt_SFd[i].Write()
for i in range(len(hEffs_B)):
  hEffs_B[i].Write()
  hEffs_C[i].Write()
  hEffs_L[i].Write()
  hEffs_B_evt[i].Write()
  hEffs_C_evt[i].Write()
  hEffs_L_evt[i].Write()
  hEffs_B_evt_SFu[i].Write()
  hEffs_C_evt_SFu[i].Write()
  hEffs_L_evt_SFu[i].Write()
  hEffs_B_evt_SFd[i].Write()
  hEffs_C_evt_SFd[i].Write()
  hEffs_L_evt_SFd[i].Write()

fOut.Close()
