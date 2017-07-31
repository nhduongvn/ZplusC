import ROOT
from myutils.util_funcs import *
import copy

'''
def scaleToLumi(fName, xSec, lumi):
  f = ROOT.TFile.Open(fName, 'read')
  hTmp = f.Get('CountPosWeight')
  hTmp1 = f.Get('CountNegWeight')
  return lumi*xSec/(hTmp.GetBinContent(1)-hTmp1.GetBinContent(1))

def getHisto(fName, histname, axis, varName, cut, scale = '1'):
  if fName == '':
    if len(varName) == 1:
      h = ROOT.TH1D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2])
    if len(varName) == 2:
      h = ROOT.TH2D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2], axis['y'][0], axis['y'][1], axis['y'][2])
    return h
  cut = '(' + cut + ')*(' + scale + ')'
  print 'Cut and scale: ', cut
  #print 'Input file: ', fName
  f = ROOT.TFile.Open(fName, 'read')
  tr = f.Get('tree')
  if len(varName) == 1:
    h = ROOT.TH1D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2])
  if len(varName) == 2:
    h = ROOT.TH2D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2], axis['y'][0], axis['y'][1], axis['y'][2])
  
  h.Sumw2()
  
  if len(varName) == 1:
    tr.Draw(varName[0] + '>>' + histname, cut)
  if len(varName) == 2:
      tr.Draw(varName[1] + ':' + varName[0] + '>>' + histname, cut)
  
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

'''

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

lumi = 35862. 

csvCut = 'CSVM'

#csvCut = '_CSVT'

vtxMassType = 'incVtxMass' #vtxMass vtxMassCorr_IVF incVtxMass

ptReweightType = 'Jet_pt' #Jet_pt,vtxPt,vtxP,incVtxPt,incVtxP 

#fOut = ROOT.TFile.Open('Test/template_' + vtxMassType + '_' + ptReweightType + '_reweight_' + csvCut + '_emu_allData_V25.root','recreate')
fOut = ROOT.TFile.Open('Test/template_' + vtxMassType + '_' + csvCut + '_emu_allData_V25.root','recreate')

ROOT.gSystem.Load('../interface/VHbbNameSpace_h.so')

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

jetPt = 'Jet_pt[' + idxJet + ']'

jsonCut = '(json == 1)'
triggerCut = '(HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v == 1 || HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v == 1 || HLT_BIT_HLT_IsoMu24_v == 1 || HLT_BIT_HLT_IsoTkMu24_v == 1)'
lepCut = '((is_emu[0] == 1) && (emu_lep_pt[0] > 25 && emu_lep_pt[1] > 25 && abs(emu_lep_eta[0]) < 2.4 && abs(emu_lep_eta[1]) < 2.4 && emu_lep_iso[0] < 0.15 && emu_lep_iso[1] < 0.15))' #id and iso already applied when save leptons in write_regression...
emuMassCut = '(VHbb::HVMass(emu_lep_pt[0],emu_lep_eta[0],emu_lep_phi[0],emu_lep_mass[0],emu_lep_pt[1],emu_lep_eta[1],emu_lep_phi[1],emu_lep_mass[1]) > 50)'

jetCut = '((' + idxJet + '>=0) && (' + jetPt + ' > 30 && abs(Jet_eta[' + idxJet + ']) < 2.4 && Jet_' + vtxMassType + '[' + idxJet + '] > 0))'
if vtxMassType == 'vtxMassCorr_IVF':
  jetCut = '(' + jetCut + ' &&(Jet_vtxCat_IVF[' +idxJet + '] == 0))'
  
metCut = '(1)'

nJet = '(Sum$(Jet_pt > 30 && abs(Jet_eta) < 2.4) >= 2)'

#@@@@@@@@@@@@zbjet cut@@@@@@@@@@@@@@@@@@
lepCut_DY = '(Vtype==0) && (vLeptons_pt[0] > 20 && vLeptons_pt[1] > 20 && abs(vLeptons_eta[0]) < 2.4 && abs(vLeptons_eta[1]) < 2.4 && vLeptons_pfRelIso04[0] < 0.25 && vLeptons_pfRelIso04[1] < 0.25)'
massCut_DY = '(VHbb::HVMass(vLeptons_pt[0],vLeptons_eta[0],vLeptons_phi[0],vLeptons_mass[0],vLeptons_pt[1],vLeptons_eta[1],vLeptons_phi[1],vLeptons_mass[1]) > 70) && (VHbb::HVMass(vLeptons_pt[0],vLeptons_eta[0],vLeptons_phi[0],vLeptons_mass[0],vLeptons_pt[1],vLeptons_eta[1],vLeptons_phi[1],vLeptons_mass[1]) < 110)'
metCut_DY = '(met_pt < 40)'

#sf = 'sign(genWeight)*puWeight*emuweight_trig[0]*emuweight[0]*bTagWeight_CSVT[0]'
#sf_DY = 'sign(genWeight)*puWeight*muweight_trig[0]*muweight[0]*bTagWeight_CSVT[0]'
#TEMP
sf = 'emuweight[0]*bTagWeight_' + csvCut + '[0]'
sf_DY = 'muweight[0]*bTagWeight_' + csvCut + '[0]'
#sf = 'bTagWeight_' + csvCut + '[0]'
#sf_DY = 'bTagWeight_' + csvCut + '[0]'

#sf = '(1.)'
#sf_DY = '(1.)'

#@@@@@@@@@@@combining cuts@@@@@@@@@@
cut_DY = lepCut_DY + '&&' + massCut_DY + '&&' + metCut_DY + '&&' + jetCut 
bCut_DY = cut_DY + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) == 5)'
cCut_DY = cut_DY + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) == 4)'
lCut_DY = cut_DY + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) != 4 && abs(Jet_hadronFlavour[' + idxJet + ']) != 5)'

cut = triggerCut + '&&' + lepCut + '&&' + emuMassCut + '&&' + metCut + '&&' + jetCut + '&&' + nJet 
bCut = cut + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) == 5)'
cCut = cut + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) == 4)'
lCut = cut + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) != 4 && abs(Jet_hadronFlavour[' + idxJet + ']) != 5)'


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_out = '/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V25/TTemu_V25/syst_fromEOS/'
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
axis1 = {'x':[40, 0, 200]} #use for getting pT reweight correction function
var1 = ['Jet_pt[' + idxJet + ']'] #use for gettting pT reweight correction function

axis = {'x':[100, 0, 10],'y':[20,0,200]}
#var = ['Jet_' + vtxMassType + '[' + idxJet + ']','Jet_pt[' + idxJet + ']']
var = ['Jet_' + vtxMassType + '[' + idxJet + ']',var1[0]]

if ptReweightType == 'vtxPt':
  var1 = ['Jet_vtxPt[' + idxJet + ']']
if ptReweightType == 'vtxP':
  var1 = ['TMath::Sqrt(Jet_vtxPx[' + idxJet + ']*Jet_vtxPx[' + idxJet + ']+Jet_vtxPy[' + idxJet + ']*Jet_vtxPy[' + idxJet + ']+Jet_vtxPz[' + idxJet + ']*Jet_vtxPz[' + idxJet + '])']
if ptReweightType == 'incVtxPt':
  var1 = ['Jet_incVtxPt[' + idxJet + ']']
if ptReweightType == 'incVtxP':  
  var1 = ['TMath::Sqrt(Jet_incVtxPx[' + idxJet + ']*Jet_incVtxPx[' + idxJet + ']+Jet_incVtxPy[' + idxJet + ']*Jet_incVtxPy[' + idxJet + ']+Jet_incVtxPz[' + idxJet + ']*Jet_incVtxPz[' + idxJet + '])']

hTmp_2D = None
hTmp_2D_cjet = None
hTmp_2D_ljet = None
for sample,att in samples.items():
  if sample == 'DY' or sample == 'TT':
    if sample == 'TT': 
        fName = path_out + att[0]
        bName  = 'bjet_pt_' + sample + '_emu'
        h[sample].append(getHisto(fName, bName, axis1, var1, bCut, sf))
    if sample == 'DY': 
        fName = path_out_DY + att[0]
        bName  = 'bjet_pt_' + sample
        h[sample].append(getHisto(fName, bName, axis1, var1, bCut_DY, sf_DY))
        hTmp_2D = getHisto(fName, bName + '_2D', axis, var, bCut_DY, sf_DY)
        cName  = 'cjet_pt_' + sample
        hTmp_2D_cjet = getHisto(fName, cName + '_2D', axis, var, cCut_DY, sf_DY)
        lName  = 'ljet_pt_' + sample
        hTmp_2D_ljet = getHisto(fName, lName + '_2D', axis, var, lCut_DY, sf_DY)

hRat = h['DY'][0].Clone('PtRatio_DYtoTT')
hTmp = h['TT'][0].Clone('Pt_tt_norm')
nInt = hRat.Integral()
if nInt > 0:
  hRat.Scale(1./nInt)
nInt = hTmp.Integral()
if nInt > 0:
  hTmp.Scale(1./nInt)
hRat.Divide(hTmp)

fOut.cd()
for sample,hists in h.items():
    for hist in hists:
        hist.Write()

hRat.Write()
hTmp_2D.Write()
hTmp_2D_cjet.Write()
hTmp_2D_ljet.Write()

#@@@@@@@@@@fill histograms@@@@@@@@@@@@

for syst in ['Central', 'JECUp', 'JECDown']: # b-jet from TT not affected by gluon splitting so no gluon splitting unc.

  h = {'DataB':[], 'DataC':[], 'DataD':[], 'DataE':[], 'DataF':[], 'DataG':[], 'DataHv1':[], 'DataHv2':[], 'DataHv3':[], 'DY':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}

  print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  print syst
  print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

    
  for sample,att in samples.items():

    fName = path_out + att[0]
    
    if att[0] == '': fName = '' 
    
    ################
    #Cuts
    ################
    jetPt_syst = jetPt
    nJet_syst = nJet 
    jer = '(1)'
    jer1 = '(1)'
    if 'Data' not in sample: 
      jer = '(Jet_corr_JER[' + idxJet + ']*(Jet_corr_JER[' + idxJet + '] > 0) + (1.0)*(Jet_corr_JER[' + idxJet + '] <= 0))'
      jer1 = '(Jet_corr_JER*(Jet_corr_JER > 0) + (1.0)*(Jet_corr_JER <= 0))'
    
    if syst == 'JECUp':
      jetPt_syst = 'Jet_rawPt['+idxJet+']*Jet_corr_JECUp['+idxJet+']*' + jer
      nJet_syst = '(Sum$(Jet_rawPt*Jet_corr_JECUp*' + jer1 + ' > 30 && abs(Jet_eta) < 2.4) >= 2)'
    if syst == 'JECDown':
      jetPt_syst = 'Jet_rawPt['+idxJet+']*Jet_corr_JECDown['+idxJet+']*' + jer
      nJet_syst = '(Sum$(Jet_rawPt*Jet_corr_JECDown*' + jer1 + ' > 30 && abs(Jet_eta) < 2.4) >= 2)'
 
    #################
    #histo name
    #################
    bName  = 'bjet_vtxMass_' + sample + '_emu'
    cName  = 'cjet_vtxMass_' + sample + '_emu'
    lName  = 'ljet_vtxMass_' + sample + '_emu'
    dataHistName = 'AllJet_vtxMass_' + sample + '_emu'
    
    if syst != 'Central':
        bName  = 'bjet_vtxMass_' + sample + '_' + syst + '_emu'
        cName  = 'cjet_vtxMass_' + sample + '_' + syst + '_emu'
        lName  = 'ljet_vtxMass_' + sample + '_' + syst + '_emu'
        dataHistName = 'AllJet_vtxMass_' + sample + '_' + syst + '_emu'

    varTmp = copy.deepcopy(var)
    varTmp[1] = varTmp[1].replace(jetPt, jetPt_syst) 
    #@@@@@@@@@@file data histogram@@@@@@@
    if sample.find('Data') != -1:
      cutData = jsonCut + '&&' + triggerCut + '&&' + cut.replace(jetPt, jetPt_syst).replace(nJet, nJet_syst)
      print '==============================='
      print 'Cut for data ', cutData
          
      h[sample].append(getHisto(fName, dataHistName, axis, varTmp, cutData))
      #print '>>>>>>>>>>>>>>>>Integral: ', h[sample][-1].GetName(), ' ', h[sample][-1].Integral()
      
      continue
    
  
    h[sample].append(getHisto(fName, bName, axis, varTmp, bCut.replace(jetPt, jetPt_syst).replace(nJet, nJet_syst), str(samples[sample][3]) + '*' + sf)) #scale to lumi
    h[sample].append(getHisto(fName, cName, axis, varTmp, cCut.replace(jetPt, jetPt_syst).replace(nJet, nJet_syst), str(samples[sample][3]) + '*' + sf))
    h[sample].append(getHisto(fName, lName, axis, varTmp, lCut.replace(jetPt, jetPt_syst).replace(nJet, nJet_syst), str(samples[sample][3]) + '*' + sf))
  
  print h
  #@@@@@@@@@@@@@@@@@@combine data and write@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  dataHistName1 = 'AllJet_vtxMass_Data_emu'
  dataHistName2 = 'AllJet_vtxMass_Data_bkgrSub_emu'
  affixName = ''
  if syst != 'Central':
    dataHistName1 = 'AllJet_vtxMass_Data_' + syst + '_emu'
    dataHistName2 = 'AllJet_vtxMass_Data_bkgrSub_' + syst + '_emu'
    affixName = syst

  hDatas = []
  for sample,hists in h.items():
    if sample.find('Data') != -1:
      if len(hists) > 0: hDatas.append(hists[0])

  if len(hDatas) > 0:
    hData = hDatas[0].Clone(dataHistName1)
    for i in range(1, len(hDatas)):
      hData.Add(hDatas[i])
      
    hData_bkgrSub = hData.Clone(dataHistName2)

    bkgrs = ['DY', 'WW', 'WZ', 'ZZ']
    hTotBkgr = h['TT'][1].Clone('AllBkgr' + affixName) #c-jet from TT
    hTotBkgr.Add(h['TT'][2]) #light-jet from TT
    hBkgr_nonTT = h['DY'][0].Clone('Bkgr_nonTT' + affixName)
    hBkgr_nonTT.Reset()
    for bkgr in bkgrs:
        for iH in h[bkgr]:
            hTotBkgr.Add(iH)
            hBkgr_nonTT.Add(iH)

    hData_bkgrSub.Add(hTotBkgr, -1)
    hData_bkgrSub_1D = hData_bkgrSub.ProjectionX()
    
    hData_bkgrSub_ptReweight = PtReweight(hData_bkgrSub, hRat, vtxMassType)
    hData_bkgrSub_ptReweight_1D = hData_bkgrSub_ptReweight.ProjectionX()

    fOut.cd()
    hData.Write()
    hData_bkgrSub.Write()
    hData_bkgrSub_ptReweight.Write()
    hData_bkgrSub_1D.Write()
    hData_bkgrSub_ptReweight_1D.Write()
    
    hTotBkgr.Write()
    hTotBkgr_1D = hTotBkgr.ProjectionX()
    hTotBkgr_1D.Write()
    
    hBkgr_nonTT.Write()
    hBkgr_nonTT_1D = hBkgr_nonTT.ProjectionX()
    hBkgr_nonTT_1D.Write()

  #@@@@@@@@@@@@@@@@write MC@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  fOut.cd()
  for sample,hists in h.items():
    for hist in hists:
      hist.Write()
  
  hBjet_tt_1D = h['TT'][0].ProjectionX()
  hBjet_tt_ptReweight = PtReweight(h['TT'][0], hRat, vtxMassType)
  hBjet_tt_ptReweight_1D = hBjet_tt_ptReweight.ProjectionX()

  hBjet_tt_1D.Write()
  hBjet_tt_ptReweight.Write()
  hBjet_tt_ptReweight_1D.Write()


#@@@@@@@@@@@@@@@get DY vtx mass template to compare with TT emu pt reweighted@@@@@@
#get vtxMass variable
varTemp = []
varTemp.append(var[0])
hVtxMass_bjet_DY = getHisto(path_out_DY + samples['DY'][0], 'bjet_vtxMass_DY', {'x':[100, 0, 10]}, varTemp, bCut_DY, str(samples['DY'][3]) + '*' + sf_DY)

fOut.cd()
hVtxMass_bjet_DY.Write()

fOut.Close() 
