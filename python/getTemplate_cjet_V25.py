import ROOT

def scaleToLumi(fName, xSec, lumi):
  f = ROOT.TFile.Open(fName, 'read')
  #hTmp = f.Get('CountPosWeight')
  #hTmp1 = f.Get('CountNegWeight')
  #return lumi*xSec/(hTmp.GetBinContent(1)-hTmp1.GetBinContent(1))
  hTmp1 = f.Get('CountWeighted')
  return lumi*xSec/(hTmp1.GetBinContent(1))

def getHisto(fName, histname, axis, varName, cut, scale = '1'):
  debug = False 
  if fName == '':
    if len(varName) == 1:
      h = ROOT.TH1D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2])
    if len(varName) == 2:
      h = ROOT.TH2D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2], axis['y'][0], axis['y'][1], axis['y'][2])
    return h
  cut = '(' + cut + ')*(' + scale + ')'
  print '====================================================='
  print 'Input file: ', fName
  print 'Cut and scale: ', cut
  f = ROOT.TFile.Open(fName, 'read')
  tr = f.Get('tree')
  if len(varName) == 1:
    h = ROOT.TH1D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2])
  if len(varName) == 2:
    h = ROOT.TH2D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2], axis['y'][0], axis['y'][1], axis['y'][2])
  
  h.Sumw2()
  
  if len(varName) == 1:
    if not debug:
      tr.Draw(varName[0] + '>>' + histname, cut)
    else:
      tr.Draw(varName[0] + '>>' + histname, cut, '', 5000)
  if len(varName) == 2:
    if not debug:
      tr.Draw(varName[1] + ':' + varName[0] + '>>' + histname, cut)
    else:
      tr.Draw(varName[1] + ':' + varName[0] + '>>' + histname, cut, '', 5000)
  
  h.SetDirectory(0)

  return h


def PtReweight(h2D, hRat):
  
  h2D_reweight = h2D.Clone(h2D.GetName()+'_ptReweight')
  h2D_reweight.Reset()

  for iPt in range(1, hRat.GetNbinsX() + 1):
    pt = hRat.GetBinLowEdge(iPt)
    w = hRat.GetBinContent(iPt)
    pt_axis = h2D.GetYaxis()
    vtxMass_axis = h2D.GetXaxis()
    
    iPt1 = pt_axis.FindFixBin(pt)
    for iVtxMass in range(1, vtxMass_axis.GetNbins()+1):
      binContent = h2D.GetBinContent(iVtxMass,iPt1)
      binError = h2D.GetBinError(iVtxMass,iPt1)
      h2D_reweight.SetBinContent(iVtxMass,iPt1,binContent*w)
      h2D_reweight.SetBinError(iVtxMass,iPt1,binError*w)
  
  return h2D_reweight

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


lumi = 35862. 

csvCut = 'CSVM'
#csvCut = 'CtagT'

vtxMassType = 'vtxMass' #incVtxMass, vtxMassCorr_IVF


fOut = ROOT.TFile.Open('Test/template_' + vtxMassType + '_' + csvCut + '_ttsemi_allData_V25.root','recreate')


jetVtxMassVar = 'Jet_' + vtxMassType 

ROOT.gSystem.Load('../interface/VHbbNameSpace_h.so')

#make cut on lepton
#make cut on jets
#plot secondary mass for c-jet

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Cuts and corrections and scale factors
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
idxJet = 'idxJet_passCSV_SVT[0]'
if vtxMassType == 'incVtxMass':
  idxJet = 'idxJet_passCSV_SVT_1[0]'
if vtxMassType == 'vtxMassCorr_IVF':
  idxJet = 'idxJet_passCSV_SVT_2[0]'

if csvCut == 'CSVM':
  idxJet = idxJet.replace('[0]','[1]')

if csvCut == 'CtagT':
  idxJet = 'idxJet_passCtag_SVT[0]'

#jsonCut = '(json == 1)'
jsonCut = '(1)'
triggerCut = '((HLT_BIT_HLT_IsoMu22_v == 1) || (HLT_BIT_HLT_IsoTkMu22_v == 1) || (HLT_BIT_HLT_IsoMu24_v == 1) || (HLT_BIT_HLT_IsoTkMu24_v == 1))' 

#lepCut = '(is_ttsemi == 1) && (ttsemi_lep_veto[1] == 1) && (ttsemi_lep_pt > 25 && abs(ttsemi_lep_eta) < 2.4 && ttsemi_lep_mediumId == 1 && ttsemi_lep_iso < 0.25)'
lepCut = '((is_ttsemi == 1) && (abs(ttsemi_lep_pdgId) == 13 && ttsemi_lep_pt > 25 && abs(ttsemi_lep_eta) < 2.4))'
#metCut = '(met_pt < 40)'
metCut = '(1)'

jetCut1 = '(ttsemi_idxJet_sortWjetCSV[0] >= 0 && ttsemi_idxJet_sortWjetCSV[1] >= 0 && ttsemi_idxJet_sortWjetCSV[2] >= 0 && ttsemi_idxJet_sortWjetCSV[3] >= 0)'

jetCut2 = '(Jet_pt[ttsemi_idxJet_sortWjetCSV[0]] > 30 && abs(Jet_eta[ttsemi_idxJet_sortWjetCSV[0]]) < 2.4 && Jet_pt[ttsemi_idxJet_sortWjetCSV[1]] > 30 && abs(Jet_eta[ttsemi_idxJet_sortWjetCSV[1]]) < 2.4 &&Jet_pt[ttsemi_idxJet_sortWjetCSV[2]] > 30 && abs(Jet_eta[ttsemi_idxJet_sortWjetCSV[2]]) < 2.4 && Jet_pt[ttsemi_idxJet_sortWjetCSV[3]] > 30 && abs(Jet_eta[ttsemi_idxJet_sortWjetCSV[3]]) < 2.4)'

jetCut3 = '(Jet_btagCSV[ttsemi_idxJet_sortWjetCSV[0]] > 0.9535 && Jet_btagCSV[ttsemi_idxJet_sortWjetCSV[1]] > 0.9535 && Jet_btagCSV[ttsemi_idxJet_sortWjetCSV[2]] > 0.8484 && Jet_btagCSV[ttsemi_idxJet_sortWjetCSV[3]] < 0.546)'

jetCut4 = '(' + jetVtxMassVar + '[ttsemi_idxJet_sortWjetCSV[2]] > 0)'

if csvCut == 'CtagT':
  jetCut3 = '(Jet_btagCSV[ttsemi_idxJet_sortWjetCSV[0]] > 0.9535 && Jet_btagCSV[ttsemi_idxJet_sortWjetCSV[1]] > 0.9535 && Jet_ctagVsL[ttsemi_idxJet_sortWjetCSV[2]] > 0.69 && Jet_ctagVsB[ttsemi_idxJet_sortWjetCSV[2]] > -0.45 && Jet_btagCSV[ttsemi_idxJet_sortWjetCSV[3]] < 0.546)'

ttsemiCut = 'ttsemi_massChi2 < 1500'

#@@@@@@@@@@@@DY cut@@@@@@@@@@@@@@@@@@@@@@
#FIXME
jetCut_DY = '((' + idxJet + '>=0) && (Jet_pt[' + idxJet + '] > 30 && abs(Jet_eta[' + idxJet + ']) < 2.4 && ' + jetVtxMassVar + '[' + idxJet + '] > 0))'
lepCut_DY = '((Vtype_new==0) && (vLeptons_new_pt[0] > 25 && vLeptons_new_pt[1] > 25 && abs(vLeptons_new_eta[0]) < 2.4 && abs(vLeptons_new_eta[1]) < 2.4 && vLeptons_new_pfRelIso04[0] < 0.25 && vLeptons_new_pfRelIso04[1] < 0.25))'
massCut_DY = '((VHbb::HVMass(vLeptons_new_pt[0],vLeptons_new_eta[0],vLeptons_new_phi[0],vLeptons_new_mass[0],vLeptons_new_pt[1],vLeptons_new_eta[1],vLeptons_new_phi[1],vLeptons_new_mass[1]) > 70) && (VHbb::HVMass(vLeptons_new_pt[0],vLeptons_new_eta[0],vLeptons_new_phi[0],vLeptons_new_mass[0],vLeptons_new_pt[1],vLeptons_new_eta[1],vLeptons_new_phi[1],vLeptons_new_mass[1]) < 110))'
metCut_DY = '(met_pt < 40)'

#sf = 'sign(genWeight)*puWeight*ttsemiweight_trig[0]*ttsemiweight[0]*bTagWeight_CSVT[0]'
#sf_DY = 'sign(genWeight)*puWeight*muweight_trig[0]*muweight[0]*bTagWeight_CSVT[0]'

#sf = '(1.0)'
#sf_DY = '(1.0)'

sf = 'ttsemiweight[0]*bTagWeight_CSVT[0]'
sf_DY = 'muweight[0]*bTagWeight_CSVT[0]'
if csvCut == 'CSVM':
  sf = 'ttsemiweight[0]*bTagWeight_CSVM[0]'
  sf_DY = 'muweight[0]*bTagWeight_CSVM[0]'

if csvCut == 'CtagT':
  sf = 'ttsemiweight[0]*cTagWeight_CSVT[0]'
  sf_DY = 'muweight[0]*cTagWeight_CSVT[0]'

#@@@@@@@@@@@final cuts@@@@@@@@@@
cut = triggerCut + '&&' + lepCut + '&&' + ttsemiCut + '&&' + jetCut1 + '&&' + jetCut2 + '&&' + jetCut3 + '&&' + jetCut4 
bCut = cut + '&& (abs(Jet_hadronFlavour[ttsemi_idxJet_sortWjetCSV[2]]) == 5)'
cCut = cut + '&& (abs(Jet_hadronFlavour[ttsemi_idxJet_sortWjetCSV[2]]) == 4)'
lCut = cut + '&& (abs(Jet_hadronFlavour[ttsemi_idxJet_sortWjetCSV[2]]) != 4 && abs(Jet_hadronFlavour[ttsemi_idxJet_sortWjetCSV[2]]) != 5)'

cut_DY = lepCut_DY + '&&' + massCut_DY + '&&' + metCut_DY + '&&' + jetCut_DY 
cCut_DY = cut_DY + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) == 4)'


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_out = '/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V25/TTsemi_V25/syst_fromEOS/'
path_out_DY = '/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V25/Zll_V25/syst_fromEOS/'
samples = {'DataB':['ZC2017_V25_DATA_SINGLE_MUON_SingleMuon__Run2016B-23Sep2016-v3.root',1, 1, 1],
          'DataC':['ZC2017_V25_DATA_SINGLE_MUON_SingleMuon__Run2016C-23Sep2016-v1.root',1, 1, 1],
          'DataD':['ZC2017_V25_DATA_SINGLE_MUON_SingleMuon__Run2016D-23Sep2016-v1.root',1, 1, 1],
          'DataE':['ZC2017_V25_DATA_SINGLE_MUON_SingleMuon__Run2016E-23Sep2016-v1.root',1, 1, 1],
          'DataF':['ZC2017_V25_DATA_SINGLE_MUON_SingleMuon__Run2016F-23Sep2016-v1.root',1, 1, 1],
          'DataG':['ZC2017_V25_DATA_SINGLE_MUON_SingleMuon__Run2016G-23Sep2016-v1.root',1, 1, 1],
          'DataHv1':['ZC2017_V25_DATA_SINGLE_MUON_SingleMuon__Run2016H-PromptReco-v1.root',1, 1, 1],
          'DataHv2':['ZC2017_V25_DATA_SINGLE_MUON_SingleMuon__Run2016H-PromptReco-v2.root',1, 1, 1],
          'DataHv3':['ZC2017_V25_DATA_SINGLE_MUON_SingleMuon__Run2016H-PromptReco-v3.root',1, 1, 1],
          'DY':['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root', 0, 1921.8*3, 1], #6025.2*1.23
          'TT':['TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root', 0, 831.76, 1],
          'WW':['WW_TuneCUETP8M1_13TeV-pythia8.root', 0, 118.7, 1],
          'WZ':['WZ_TuneCUETP8M1_13TeV-pythia8.root', 0, 47.13, 1],
          'ZZ':['ZZ_TuneCUETP8M1_13TeV-pythia8.root', 0, 16.523, 1]}
#name,type=data/MC, xSection, scale

for sample,att in samples.items():
    
  fName = path_out + att[0]
  
  if att[0] == '': fName = '' 
  
  if att[1] != 1 and fName != '':
    att[3] = scaleToLumi(fName, att[2], lumi)
    samples[sample] = att

  print sample, samples[sample]


#@@@@@@@@@@get pt weight@@@@@@@@
h = {'DY':[], 'TT':[]}
axis = {'x':[20, 0, 200]}
var_DY = ['Jet_pt[' + idxJet + ']']
var_TT = ['Jet_pt[ttsemi_idxJet_sortWjetCSV[2]]']
for sample,att in samples.items():
    if sample == 'DY': 
      fName = path_out_DY + att[0]
      bName  = 'cjet_pt_' + sample
      h[sample].append(getHisto(fName, bName, axis, var_DY, cCut_DY, sf_DY))
    if sample == 'TT': 
      fName = path_out + att[0]
      cName  = 'cjet_pt_' + sample
      h[sample].append(getHisto(fName, cName, axis, var_TT, cCut, sf))


hRat = h['DY'][0].Clone('PtRatio_DYtoTT')
hTmp = h['TT'][0].Clone('Pt_tt_norm')
hRat.Scale(1./hRat.Integral())
hTmp.Scale(1./hTmp.Integral())
hRat.Divide(hTmp)

fOut.cd()
for sample,hists in h.items():
    for hist in hists:
        hist.Write()

hRat.Write()

#@@@@@@@@@@fill histograms@@@@@@@@@@@@
for syst in ['Central', 'JECUp', 'JECDown']: #actually not use syst.

  print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  print syst
  print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  
  ################
  #Cuts
  ################
  


  h = {'DataB':[], 'DataC':[], 'DataD':[], 'DataE':[], 'DataF':[], 'DataG':[], 'DataHv1':[], 'DataHv2':[], 'DataHv3':[], 'DY':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}
  
  for sample,att in samples.items():

    fName = path_out + att[0]
    
    if att[0] == '': fName = '' 
    
    ########################
    #modify cut and var according to syst
    ########################
    jetPt_syst = 'Jet_pt[ttsemi_idxJet_sortWjetCSV[2]]' 
    jer = '(1)'
    if 'Data' not in sample: 
      jer = '(Jet_corr_JER[ttsemi_idxJet_sortWjetCSV[2]]*(Jet_corr_JER[ttsemi_idxJet_sortWjetCSV[2]] > 0) + (1.0)*(Jet_corr_JER[ttsemi_idxJet_sortWjetCSV[2]] <= 0))'
    
    if syst == 'JECUp':
      jetPt_syst = 'Jet_rawPt[ttsemi_idxJet_sortWjetCSV[2]]*Jet_corr_JECUp[ttsemi_idxJet_sortWjetCSV[2]]*' + jer
    if syst == 'JECDown':
      jetPt_syst = 'Jet_rawPt[ttsemi_idxJet_sortWjetCSV[2]]*Jet_corr_JECDown[ttsemi_idxJet_sortWjetCSV[2]]*' + jer
    
    jetPt_syst_0 = jetPt_syst.replace('[2]','[0]')
    jetPt_syst_1 = jetPt_syst.replace('[2]','[1]')
    jetPt_syst_3 = jetPt_syst.replace('[2]','[3]')

    axis = {'x':[100, 0, 10],'y':[20,0,200]}
    var = ['Jet_' + vtxMassType + '[ttsemi_idxJet_sortWjetCSV[2]]','Jet_pt[ttsemi_idxJet_sortWjetCSV[2]]']
    if syst != 'Central':
      var[1] = var[1].replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[2]]',jetPt_syst)

    #@@@@@@@@@@file data histogram@@@@@@@
    if sample.find('Data') != -1:
      cutData = jsonCut + '&&' + cut.replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[2]]',jetPt_syst).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[0]]',jetPt_syst_0).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[1]]', jetPt_syst_1).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[3]]',jetPt_syst_3)
      
      print '==============================='
      print 'Cut for data ', cutData
      h[sample].append(getHisto(fName, 'AllJet_vtxMass_' + sample + '_semilep', axis, var, cutData))
      
      continue
    
    
    bName  = 'bjet_vtxMass_' + sample + '_semilep'
    cName  = 'cjet_vtxMass_' + sample + '_semilep'
    lName  = 'ljet_vtxMass_' + sample + '_semilep'
    bCutTmp = bCut.replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[2]]', jetPt_syst).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[0]]', jetPt_syst_0).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[1]]', jetPt_syst_1).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[3]]', jetPt_syst_3)
    cCutTmp = cCut.replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[2]]', jetPt_syst).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[0]]', jetPt_syst_0).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[1]]', jetPt_syst_1).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[3]]', jetPt_syst_3)
    lCutTmp = lCut.replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[2]]', jetPt_syst).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[0]]', jetPt_syst_0).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[1]]', jetPt_syst_1).replace('Jet_pt[ttsemi_idxJet_sortWjetCSV[3]]', jetPt_syst_3)
    
    h[sample].append(getHisto(fName, bName, axis, var, bCutTmp, str(samples[sample][3]) + '*' + sf)) #scale to lumi
    h[sample].append(getHisto(fName, cName, axis, var, cCutTmp, str(samples[sample][3]) + '*' + sf))
    h[sample].append(getHisto(fName, lName, axis, var, lCutTmp, str(samples[sample][3]) + '*' + sf))


  #@@@@@@@@@@@@@@@@@print all Data and MC@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  print 'Event yeilds: '
  for sample,hists in h.items():
    print sample,
    for iH in hists:
      print iH.Integral(0, 100000, 0, 100000),
    print '\n'

  #@@@@@@@@@@@@@@@@@@combine data and write@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  hData = None
  hBkgr = None
  hTT = None
  firstData = True
  firstTT = True
  firstBkgr = True

  for sample,hists in h.items():
    for iH in hists:
      if sample.find('Data') != -1:
        if firstData:
          name = 'AllJet_vtxMass_Data_semilep'
          if syst != 'Central':
            name = 'AllJet_vtxMass_Data_' + syst + '_semilep'
          hData = iH.Clone(name)
          firstData = False
        else: hData.Add(iH)
      if 'TT' in sample:
        if firstTT:
          name = 'AllJet_vtxMass_TT_semilep'
          if syst != 'Central':
            name = 'AllJet_vtxMass_TT_' + syst + '_semilep'
          hTT = iH.Clone(name)
          firstTT = False
        else: hTT.Add(iH)
      if 'TT' not in sample and 'Data' not in sample:
        if firstBkgr:
          name = 'AllJet_vtxMass_bkgr_semilep'
          if syst != 'Central':
            name = 'AllJet_vtxMass_bkgr_' + syst + '_semilep'
          hBkgr = iH.Clone(name)
          firstBkgr = False
        else: hBkgr.Add(iH)


  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  print '========================='
  print 'Observed data: ', hData.Integral(0, 10000, 0, 10000)
  print 'TT           : ', hTT.Integral(0,10000, 0, 10000)
  print 'Bkgr         : ', hBkgr.Integral(0, 1000, 0, 10000)

  name = 'AllJet_vtxMass_Data_bkgrSub_semilep'
  if syst != 'Central':
    name = 'AllJet_vtxMass_Data_bkgrSub_' + syst + '_semilep_'
  hData_bkgrSub = hData.Clone(name)
  hData_bkgrSub.Add(hBkgr, -1.0)

  scale_data_TT = hData_bkgrSub.Integral(0,10000, 0, 10000)/hTT.Integral(0, 10000, 0, 10000)
  print 'Scale factor Data_bkgr_subtracted/TT: ', scale_data_TT

  name = 'AllJet_vtxMass_Data_bkgrSub_TTsub_semilep'
  if syst != 'Central':
    name = 'AllJet_vtxMass_Data_bkgrSub_TTsub_' + syst + '_semilep'
  hData_bkgrSub_1 = hData_bkgrSub.Clone(name)
  hData_bkgrSub_1.Add(h['TT'][0], -1.*scale_data_TT) #bjet
  hData_bkgrSub_1.Add(h['TT'][2], -1.*scale_data_TT) #light jet

  hData_bkgrSub_1_ptReweight = PtReweight(hData_bkgrSub_1, hRat)

  fOut.cd()
  hData.Write()
  (hData.ProjectionX()).Write()
  hTT.Write()
  (hTT.ProjectionX()).Write()
  hBkgr.Write()
  (hBkgr.ProjectionX()).Write()

  hData_bkgrSub_1.Write()
  (hData_bkgrSub_1.ProjectionX()).Write()
  hData_bkgrSub_1_ptReweight.Write()
  (hData_bkgrSub_1_ptReweight.ProjectionX()).Write()

  #@@@@@@@@@@@@@@@@write MC@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  fOut.cd()
  for sample,hists in h.items():
    for hist in hists:
      hist.Write()
      (hist.ProjectionX()).Write()

  hCjet_tt_ptReweight = PtReweight(h['TT'][1], hRat)
  hCjet_tt_ptReweight_1D = hCjet_tt_ptReweight.ProjectionX()
  hCjet_tt_ptReweight.Write()
  hCjet_tt_ptReweight_1D.Write()



#@@@@@@@@@@@@@@@get DY vtx mass template to compare with TT semi pt reweighted@@@@@@
varTemp = []
varTemp.append('Jet_' + vtxMassType + '[' + idxJet + ']')
hVtxMass_cjet_DY = getHisto(path_out_DY + samples['DY'][0], 'cjet_vtxMass_DY', {'x':[100, 0, 10]}, varTemp, cCut_DY, str(samples['DY'][3]) + '*' + sf_DY)

fOut.cd()
hVtxMass_cjet_DY.Write()

fOut.Close()

  
