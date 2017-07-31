import ROOT
from myutils.util_funcs import *

'''
def scaleToLumi(fName, xSec, lumi):
  f = ROOT.TFile.Open(fName, 'read')
  hTmp = f.Get('CountPosWeight')
  hTmp1 = f.Get('CountNegWeight')
  return lumi*xSec/(hTmp.GetBinContent(1)-hTmp1.GetBinContent(1))

def getHisto(fName, histname, axis, varName, cut, scale = '1'):
  debug = False 
  if fName == '':
    if len(varName) == 1:
      h = ROOT.TH1D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2])
    if len(varName) == 2:
      h = ROOT.TH2D(histname, histname, axis['x'][0], axis['x'][1], axis['x'][2], axis['y'][0], axis['y'][1], axis['y'][2])
    return h
  cut = '(' + cut + ')*(' + scale + ')'
  print '=============================='
  print fName
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
    if not debug:
      tr.Draw(varName[0] + '>>' + histname, cut)
    else:
      tr.Draw(varName[0] + '>>' + histname, cut,'',10000)
    
  if len(varName) == 2:
    if not debug:
      tr.Draw(varName[1] + ':' + varName[0] + '>>' + histname, cut)
    else:
      tr.Draw(varName[1] + ':' + varName[0] + '>>' + histname, cut, '', 10000)
  
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
# Main program
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

lumi = 35862. 
#lumi = 1. 

csvCut = 'CSVM'
csvCut_val = '0.8484'

#csvCut = '_CSVT'

vtxMassType = 'vtxMass' # vtxMass, incVtxMass, vtxMassCorr_IVF

fOut = ROOT.TFile.Open('Test/template_CSVv2_deepCSV_' + vtxMassType + '_' + csvCut + '_iso0p05_Wjet_allData_V25.root','recreate')

ROOT.gSystem.Load('../interface/VHbbNameSpace_h.so')

jetVtxMassVar = 'Jet_' + vtxMassType 

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
#jsonCut = '(1)'
triggerCut = '(HLT_BIT_HLT_IsoTkMu24_v == 1 || HLT_BIT_HLT_IsoMu24_v == 1 || HLT_BIT_HLT_IsoTkMu27_v == 1 || HLT_BIT_HLT_IsoMu27_v == 1)' 
lepCut = '(Vtype_new == 2 && (vLeptons_new_pt[0] > 25 && abs(vLeptons_new_eta[0]) < 2.4 && vLeptons_new_pfRelIso04[0] < 0.05))' #previous cut is 0.01
#0.12
MTmassCut = '(VHbb::MTmass(vLeptons_new_pt[0],vLeptons_new_phi[0],met_pt, met_phi) > 50)'
#MTmassCut = '(1)'
jetCut = '((' + idxJet + '>= 0) && ' + jetPt + ' > 30 && abs(Jet_eta[' + idxJet + ']) < 2.4 && Jet_btagCSV[' + idxJet + '] > ' + csvCut_val + ' && Jet_btagDeepCSVdusg[' + idxJet + '] >= 0.05 && Jet_vtxCat_IVF[' +idxJet + '] == 0 && ' + jetVtxMassVar + '[' + idxJet + '] > 0)'

if vtxMassType == 'vtxMass' or vtxMassType == 'incVtxMass':
  jetCut = jetCut.replace('Jet_vtxCat_IVF[0] == 0','(1)')
  jetCut = jetCut.replace('Jet_vtxCat_IVF[0]==0','(1)')

#nJetCut = '(Sum$(Jet_pt > 25 && abs(Jet_eta) < 2.4) < 100)'
nJetCut = '(1)'
#metCut = '(met_pt > 40)'
metCut = '(1)'

#@@@@@@@@@@@@zbjet cut@@@@@@@@@@@@@@@@@@
lepCut_DY = '((Vtype_new==0) && (vLeptons_new_pt[0] > 25 && vLeptons_new_pt[1] > 25 && abs(vLeptons_new_eta[0]) < 2.4 && abs(vLeptons_new_eta[1]) < 2.4 && vLeptons_new_pfRelIso04[0] < 0.25 && vLeptons_new_pfRelIso04[1] < 0.25))'
jetCut_DY = '((' + idxJet + '>=0) && (' + jetPt + ' > 30 && abs(Jet_eta[' + idxJet + ']) < 2.4 && ' + jetVtxMassVar + '[' +idxJet + '] > 0 && Jet_vtxCat_IVF[' + idxJet + '] == 0))'
if vtxMassType == 'vtxMass':
  jetCut_DY = jetCut_DY.replace('Jet_vtxCat_IVF[' + idxJet + '] == 0','(1)')
  jetCut_DY = jetCut_DY.replace('Jet_vtxCat_IVF[' + idxJet + ']==0','(1)')

massCut_DY = '((VHbb::HVMass(vLeptons_new_pt[0],vLeptons_new_eta[0],vLeptons_new_phi[0],vLeptons_new_mass[0],vLeptons_new_pt[1],vLeptons_new_eta[1],vLeptons_new_phi[1],vLeptons_new_mass[1]) > 70) && (VHbb::HVMass(vLeptons_new_pt[0],vLeptons_new_eta[0],vLeptons_new_phi[0],vLeptons_new_mass[0],vLeptons_new_pt[1],vLeptons_new_eta[1],vLeptons_new_phi[1],vLeptons_new_mass[1]) < 110))'
metCut_DY = '(met_pt < 40)'

###########scale factor###############
sf = 'sign(genWeight)*puWeight*emuweight_trig[0]*emuweight[0]*bTagWeight_CSVT[0]'
sf_DY = 'sign(genWeight)*puWeight*muweight_trig[0]*muweight[0]*bTagWeight_CSVT[0]'

sf = 'bTagWeight_CSVT[0]'
sf_DY = 'muweight[0]*bTagWeight_CSVT[0]'
if csvCut == 'CSVM':
  sf = 'bTagWeight_CSVM[0]'
  sf_DY = 'muweight[0]*bTagWeight_CSVM[0]'

#@@@@@@@@@@@final cuts@@@@@@@@@@
cut = triggerCut + '&&' + lepCut + '&&' + jetCut + '&&' + metCut + '&&' + MTmassCut + '&&' + nJetCut
bCut = cut + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) == 5)'
cCut = cut + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) == 4)'
lCut = cut + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) != 4 && abs(Jet_hadronFlavour[' + idxJet + ']) != 5)'

cut_DY = lepCut_DY + '&&' + massCut_DY + '&&' + metCut_DY + '&&' + jetCut_DY 
lCut_DY = cut_DY + '&& (abs(Jet_hadronFlavour[' + idxJet + ']) != 4 && abs(Jet_hadronFlavour[' + idxJet + ']) != 5)'

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_out = '/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V25/Wjet_V25/syst_fromEOS/'
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
          #'QCD_15to30':['QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8.root',0, 1837410000,1],
          #'QCD_30to50':['QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8.root',0, 140932000,1],
          #'QCD_50to80':['QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8.root',0, 19204300,1],
          #'QCD_80to120':['QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8.root',0, 2762530,1],
          #'QCD_120to170':['QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8.root',0, 471100,1],
          #'QCD_170to300':['QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8.root',0, 117276,1],
          'Wjet':['WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root',0,20508.9*3,1],
          'DY':['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root', 0, 1921.8*3, 1], #6025.2*1.23
          'TT':['TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root', 0, 831.76, 1],
          'WW':['WW_TuneCUETP8M1_13TeV-pythia8.root', 0, 118.7, 1],
          'WZ':['WZ_TuneCUETP8M1_13TeV-pythia8.root', 0, 47.13, 1],
          'ZZ':['ZZ_TuneCUETP8M1_13TeV-pythia8.root', 0, 16.523, 1]
}
#name,type=data/MC, xSection, scale


for sample,att in samples.items():
    
  fName = path_out + att[0]
  
  if att[0] == '': fName = '' 
  
  if att[1] != 1 and fName != '':
    att[3] = scaleToLumi(fName, att[2], lumi)
    samples[sample] = att

  print sample, samples[sample]


#@@@@@@@@@@get pt weight@@@@@@@@
#h = {'DY':[], 'QCD_all':[], 'QCD_15to30':[], 'QCD_30to50':[], 'QCD_50to80':[], 'QCD_80to120':[], 'QCD_120to170':[],'QCD_170to300':[]}
h = {'DY':[], 'Wjet':[]}
axis = {'x':[20, 0, 200]}

for sample,att in samples.items():
    if sample.find('Wjet') != -1: 
        fName = path_out + att[0]
        lName  = 'ljet_pt_' + sample + '_Wjet'
        h[sample].append(getHisto(fName, lName, axis, [jetPt], lCut, str(samples[sample][3]) + '*' + sf)) # scale to luminosity and pthat weight
    if sample == 'DY': 
        fName = path_out_DY + att[0]
        lName  = 'ljet_pt_' + sample
        h[sample].append(getHisto(fName, lName, axis, [jetPt], lCut_DY, str(samples[sample][3]) + '*' + sf_DY))

#adding QCD together
#hTmp = h['QCD_15to30'][0].Clone('QCD_pt_all')
#hTmp.Reset()
#for sample,hists in h.items():
#  if sample.find('QCD') != -1 and sample.find('QCD_all') == -1:
#    if len(hists) != 1: print 'Warning in adding QCD to get QCD pt at len(hists) ', len(hists), '. Should be 1'
#    hTmp.Add(hists[0])
#h['QCD_all'].append(hTmp)

hRat = h['DY'][0].Clone('PtRatio_DYtoWjet')
hTmp = h['Wjet'][0].Clone('Pt_tt_norm')
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


for syst in ['Central', 'JECUp', 'JECDown']:
  #@@@@@@@@@@fill histograms@@@@@@@@@@@@
  #h = {'DataB':[], 'DataC':[], 'DataD':[], 'DataE':[], 'DataF':[], 'DataG':[], 'DataHv1':[], 'DataHv2':[], 'DataHv3':[], 'QCD_15to30':[], 'QCD_30to50':[], 'QCD_50to80':[], 'QCD_80to120':[], 'QCD_120to170':[],'QCD_170to300':[], 'DY':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}
  #h = {'DataB':[], 'DataC':[], 'DataD':[], 'DataE':[], 'DataF':[], 'DataG':[], 'DataHv1':[], 'DataHv2':[], 'DataHv3':[], 'Wjet':[], 'QCD_15to30':[], 'QCD_30to50':[], 'QCD_50to80':[], 'QCD_80to120':[], 'QCD_120to170':[],'QCD_170to300':[], 'TT':[], 'DY':[], 'WW':[], 'WZ':[], 'ZZ':[]}
  h = {'DataB':[], 'DataC':[], 'DataD':[], 'DataE':[], 'DataF':[], 'DataG':[], 'DataHv1':[], 'DataHv2':[], 'DataHv3':[], 'Wjet':[], 'DY':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}

  print '###########################'
  print syst
  print '###########################'

  for sample,att in samples.items():

    fName = path_out + att[0]
    
    if att[0] == '': fName = '' 
    
    ################
    #modify cuts according to systs.
    ################
    jetPt_syst = jetPt
    #nJet_syst = nJet 
    jer = '(1)'
    jer1 = '(1)'
    if 'Data' not in sample: 
      jer = '(Jet_corr_JER[' + idxJet + ']*(Jet_corr_JER[' + idxJet + '] > 0) + (1.0)*(Jet_corr_JER[' + idxJet + '] <= 0))'
      jer1 = '(Jet_corr_JER*(Jet_corr_JER > 0) + (1.0)*(Jet_corr_JER <= 0))'
    
    if syst == 'JECUp':
      jetPt_syst = 'Jet_rawPt['+idxJet+']*Jet_corr_JECUp['+idxJet+']*' + jer
      #nJet_syst = '(Sum$(Jet_rawPt*Jet_corr_JECUp*' + jer1 + ' > 30 && abs(Jet_eta) < 2.4) >= 2)'
    if syst == 'JECDown':
      jetPt_syst = 'Jet_rawPt['+idxJet+']*Jet_corr_JECDown['+idxJet+']*' + jer
      #nJet_syst = '(Sum$(Jet_rawPt*Jet_corr_JECDown*' + jer1 + ' > 30 && abs(Jet_eta) < 2.4) >= 2)'

    axis = {'x':[100, 0, 10],'y':[20,0,200]}
    var = [jetVtxMassVar + '[' + idxJet + ']',jetPt]    
    
    if syst != 'Central':
      var[1] = var[1].replace(jetPt,jetPt_syst)
    
    #####################################
    #Histo name
    #####################################
    dataHistName = 'AllJet_vtxMass_' + sample + '_Wjet'
    bName  = 'bjet_vtxMass_' + sample + '_Wjet'
    cName  = 'cjet_vtxMass_' + sample + '_Wjet'
    lName  = 'ljet_vtxMass_' + sample + '_Wjet'
    
    if syst != 'Central':
      dataHistName = 'AllJet_vtxMass_' + sample + '_' + syst + '_Wjet'
      bName  = 'bjet_vtxMass_' + sample + '_' + syst + '_Wjet'
      cName  = 'cjet_vtxMass_' + sample + '_' + syst + '_Wjet'
      lName  = 'ljet_vtxMass_' + sample + '_' + syst + '_Wjet'

    #@@@@@@@@@@file data histogram@@@@@@@
    if sample.find('Data') != -1:
      cutData = jsonCut + '&&' + triggerCut + '&&' + cut
      h[sample].append(getHisto(fName, dataHistName, axis, var, cutData.replace(jetPt, jetPt_syst)))
      
      continue
    
    
    h[sample].append(getHisto(fName, bName, axis, var, bCut.replace(jetPt, jetPt_syst), str(samples[sample][3]) + '*' + sf)) #scale to lumi
    h[sample].append(getHisto(fName, cName, axis, var, cCut.replace(jetPt, jetPt_syst), str(samples[sample][3]) + '*' + sf))
    h[sample].append(getHisto(fName, lName, axis, var, lCut.replace(jetPt, jetPt_syst), str(samples[sample][3]) + '*' + sf))

  #@@@@@@@@@@@@@@@@@@combine data and write@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  name = 'AllJet_vtxMass_Data_Wjet'
  if syst != 'Central':
    name =  'AllJet_vtxMass_Data_' + syst + '_Wjet'
  hData = h['DataB'][0].Clone(name)
  hData.Reset()

  for sample,hists in h.items():
    if sample.find('Data') != -1:
      if len(hists) > 0: hData.Add(hists[0])
      else: print sample, ' has ', len(hists), ' histogram. Should be 1'
      
  name = 'AllJet_vtxMass_Data_bkgrSub_Wjet'
  if syst != 'Central':
    name =  'AllJet_vtxMass_Data_bkgrSub_' + syst + '_Wjet'
  hData_bkgrSub = hData.Clone(name)
  #bkgrs_nonWjet = ['DY', 'TT', 'WW', 'WZ', 'ZZ','QCD_15to30','QCD_30to50','QCD_50to80','QCD_80to120','QCD_120to170','QCD_170to300']
  bkgrs_nonWjet = ['DY', 'TT', 'WW', 'WZ', 'ZZ']
  name = 'Bkgr_nonWjet'
  if syst != 'Central':
    name =  'Bkgr_nonWjet_' + syst
  hBkgr_nonWjet = h['DY'][0].Clone(name)
  hBkgr_nonWjet.Reset()
  for bkgr in bkgrs_nonWjet:
      for iH in h[bkgr]:
          hBkgr_nonWjet.Add(iH)

  name = 'AllMC'
  if syst != 'Central':
    name =  'AllMC_' + syst
  hMC = hBkgr_nonWjet.Clone(name)

  name = 'AllBkgr'
  if syst != 'Central':
    name =  'AllBkgr_' + syst
  hTotBkgr = hBkgr_nonWjet.Clone(name)
  if len(h['Wjet']) >= 3:
    hMC.Add(h['Wjet'][0])
    hMC.Add(h['Wjet'][1])
    hMC.Add(h['Wjet'][2])
    hTotBkgr.Add(h['Wjet'][0]) #b
    hTotBkgr.Add(h['Wjet'][1]) #c
  else:
    print '>>>hWjet should has three histograms'

  #scale_factor = hData.Integral()/hMC.Integral()
  #print 'Scale factor for background is: ', scale_factor
  #hBkgr_nonWjet.Scale(scale_factor)
  #hTotBkgr.Scale(scale_factor)

  hData_bkgrSub.Add(hTotBkgr, -1)
  hData_bkgrSub_1D = hData_bkgrSub.ProjectionX()

  hData_bkgrSub_ptReweight = PtReweight(hData_bkgrSub, hRat)
  hData_bkgrSub_ptReweight_1D = hData_bkgrSub_ptReweight.ProjectionX()

  fOut.cd()
  hData_1D = hData.ProjectionX()
  hData.Write()
  hData_1D.Write()
  hData_bkgrSub.Write()
  hData_bkgrSub_1D.Write()
  hData_bkgrSub_ptReweight.Write()
  hData_bkgrSub_ptReweight_1D.Write()

  hBkgr_nonWjet_1D = hBkgr_nonWjet.ProjectionX()
  hBkgr_nonWjet.Write()
  hBkgr_nonWjet_1D.Write()

  hTotBkgr_1D = hTotBkgr.ProjectionX()
  hTotBkgr.Write()
  hTotBkgr_1D.Write()

  #@@@@@@@@@@@@@@@@write MC@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  print 'Start to write the output'

  hBjet_Wjet_1D = h['Wjet'][0].ProjectionX()
  hCjet_Wjet_1D = h['Wjet'][1].ProjectionX()
  hLjet_Wjet_1D = h['Wjet'][2].ProjectionX()

  hLjet_Wjet_ptReweight = PtReweight(h['Wjet'][2], hRat)
  hLjet_Wjet_ptReweight_1D = hLjet_Wjet_ptReweight.ProjectionX()

  fOut.cd()
  for sample,hists in h.items():
    for hist in hists:
      hist.Write()

  hLjet_Wjet_1D.Write()
  hLjet_Wjet_ptReweight.Write()
  hLjet_Wjet_ptReweight_1D.Write()

  hCjet_Wjet_1D.Write()
  hBjet_Wjet_1D.Write()



#@@@@@@@@@@@@@@@get DY vtx mass template to compare with Wjet pt reweighted@@@@@@
hVtxMass_ljet_DY = None
hVtxMass_ljet_DY = getHisto(path_out_DY + samples['DY'][0], 'ljet_vtxMass_DY', {'x':[100, 0, 10]}, [jetVtxMassVar + '[' + idxJet + ']'], lCut_DY, str(samples['DY'][3]) + '*' + sf_DY)

fOut.cd()
hVtxMass_ljet_DY.Write()

print 'Closing output'
fOut.Close()

if __name__ == '__main__':
   rep = ''
   while not rep in [ 'q', 'Q' ]:
     rep = raw_input( 'enter "q" to quit: ' )
     if 1 < len(rep):
       rep = rep[0]

