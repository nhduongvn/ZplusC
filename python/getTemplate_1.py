import ROOT


def scaleToLumi(fName, xSec, lumi):
  f = ROOT.TFile.Open(fName, 'read')
  hTmp = f.Get('CountPosWeight')
  hTmp1 = f.Get('CountNegWeight')
  return lumi*xSec/(hTmp.GetBinContent(1)-hTmp1.GetBinContent(1))

def getHisto(fName, histname, cut, csvVal = '0.935', scale = '1'):
  if fName == '':
    h = ROOT.TH1D(histname, histname, 100, 0, 10)
    return h
  cut = '(' + cut + ')*(' + scale + ')'
  #print 'Cut and scale: ', cut
  #print 'Input file: ', fName
  f = ROOT.TFile.Open(fName, 'read')
  tr = f.Get('tree')
  h = ROOT.TH1D(histname, histname, 100, 0, 10)
  h.Sumw2()
  if csvVal == '0.935':
    tr.Draw('Jet_vtxMass[idxJet_passCSV_SVT[0]]>>' + histname, cut)
  if csvVal == '0.8':
    tr.Draw('Jet_vtxMass[idxJet_passCSV_SVT[1]]>>' + histname, cut)
  h.SetDirectory(0)
  return h

def combineDYtemplate(hDYinc_lheNj_1, hDYinc_lheNj_not1, hDY1Jet_lheNj_1, combineHist_name):
  hCom = hDY1Jet_lheNj_1.Clone(combineHist_name)
  hCom.Add(hDYinc_lheNj_1)
  scale = 1
  if hDYinc_lheNj_1.Integral(2,1000) > 0:
    scale = hCom.Integral(2, 1000)/hDYinc_lheNj_1.Integral(2, 1000)
  hTmp = hDYinc_lheNj_not1.Clone(hDYinc_lheNj_not1.GetName()+ '_tmp')
  hTmp.Scale(scale)
  hCom.Add(hTmp)
  return hCom





lumi = 7000. 

#csvCut = '_CSVM'
#csvVal = '0.8'

csvCut = '_CSVT'
csvVal = '0.935'

channels = ['Zmm', 'Zee']


fOut = ROOT.TFile.Open('template' + csvCut + '_allWeight.root','recreate')
#fOut = ROOT.TFile.Open('testTemplate.root','recreate')

ROOT.gSystem.Load('../interface/VHbbNameSpace_h.so')


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Cuts and corrections and scale factors
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
jsonCut = '(json == 1)'
triggerCut = '((HLT_BIT_HLT_IsoMu22_v == 1) || (HLT_BIT_HLT_IsoTkMu22_v == 1))' 
jetCut = '(idxJet_passCSV_SVT[0]>=0) && (Jet_pt[idxJet_passCSV_SVT[0]] > 30 && abs(Jet_eta[idxJet_passCSV_SVT[0]]) < 2.4)'
lepCut = '(vLeptons_pt[0] > 25 && vLeptons_pt[1] > 25 && abs(vLeptons_eta[0]) < 2.4 && abs(vLeptons_eta[1]) < 2.4 && vLeptons_pfRelIso04[0] < 0.25 && vLeptons_pfRelIso04[1] < 0.25)'
zjetMassCut = '(VHbb::HVMass(vLeptons_pt[0],vLeptons_eta[0],vLeptons_phi[0],vLeptons_mass[0],vLeptons_pt[1],vLeptons_eta[1],vLeptons_phi[1],vLeptons_mass[1]) > 70) && (VHbb::HVMass(vLeptons_pt[0],vLeptons_eta[0],vLeptons_phi[0],vLeptons_mass[0],vLeptons_pt[1],vLeptons_eta[1],vLeptons_phi[1],vLeptons_mass[1]) < 110)'
metCut = '(met_pt < 40)'
#metCut = '(1)'

sf_com = 'sign(genWeight)*puWeight'
sf_chan = [sf_com + '*muweight_trig[0]*muweight[0]*bTagWeight_CSVT[0]', sf_com + '*(1)']
#sf_chan = ['(1)', '(1)']

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_out = '/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V24/Zmm_V24/syst/'
samples = {'DataC':['ZC2016_11_Run2016_promptReco_3_SingleMuon__Run2016C-PromptReco-v2.root',1, 1, 1],
          'DataD':['ZC2016_11_Run2016_promptReco_3_SingleMuon__Run2016D-PromptReco-v2.root',1, 1, 1],
          'DY':['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root', 0, 1921.8*3, 1], #6025.2*1.23
#          'DY':['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root', 0, 1921.8*3, 1], #6025.2*1.23
          'DY1Jet':['', 0, 1, 1],
          'TT':['TT_TuneCUETP8M1_13TeV-powheg-pythia8.root', 0, 831.76, 1],
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



for chan in channels:
  #@@@@@@@@@@@final cuts@@@@@@@@@@
  Zll_cut = lepCut
  if chan == 'Zmm':
    Zll_cut = '(Vtype == 0) && ' + Zll_cut 
  if chan == 'Zee':
    Zll_cut = '(Vtype == 1) && ' + Zll_cut 
  cut = Zll_cut + '&&' + zjetMassCut + '&&' + metCut + '&&' + jetCut 
  bCut = cut + '&& (abs(Jet_hadronFlavour[idxJet_passCSV_SVT[0]]) == 5)'
  cCut = cut + '&& (abs(Jet_hadronFlavour[idxJet_passCSV_SVT[0]]) == 4)'
  lCut = cut + '&& (abs(Jet_hadronFlavour[idxJet_passCSV_SVT[0]]) != 4 && abs(Jet_hadronFlavour[idxJet_passCSV_SVT[0]]) != 5)'
  
  #@@@@@@@@@scale factor@@@@@@@@@
  sf = '(1)'
  if chan == 'Zmm':
    sf = sf_chan[0]
  if chan == 'Zee':
    sf = sf_chan[1] 


  #@@@@@@@@@@fill histograms@@@@@@@@@@@@
  h = {'DataC':[], 'DataD':[], 'DY':[], 'DY1Jet':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}
  hDY = {'DY':[], 'DY1Jet':[], 'DYcombine':[]} #b, c, l (lheNj==1,lheNj!=1) 

  for sample,att in samples.items():

    fName = path_out + att[0]
  
    if att[0] == '': fName = '' 

    #@@@@@@@@@@file data histogram@@@@@@@
    if sample.find('Data') != -1:
      cutData = jsonCut + '&&' + triggerCut + '&&' + cut
      print '==============================='
      print 'Cut for data ', cutData
      h[sample].append(getHisto(fName, 'AllJet_vtxMass_' + sample + '_' + chan, cutData))
    
      continue
  
  
    bName  = 'bjet_vtxMass_' + sample + '_' + chan
    cName  = 'cjet_vtxMass_' + sample + '_' + chan
    lName  = 'ljet_vtxMass_' + sample + '_' + chan
    h[sample].append(getHisto(fName, bName, bCut, csvVal, str(samples[sample][3]) + '*' + sf))
    h[sample].append(getHisto(fName, cName, cCut, csvVal, str(samples[sample][3]) + '*' + sf))
    h[sample].append(getHisto(fName, lName, lCut, csvVal, str(samples[sample][3]) + '*' + sf))

  
    if (sample == 'DY' or sample == 'DY1Jet'):
    
        bNameTmp = 'bjet_vtxMass_tmp_' + sample + '_' + chan
        cNameTmp = 'cjet_vtxMass_tmp_' + sample + '_' + chan
        lNameTmp = 'ljet_vtxMass_tmp_' + sample + '_' + chan
        
        h[sample].append(getHisto(fName, bNameTmp, bCut, csvVal, sf))
        h[sample].append(getHisto(fName, cNameTmp, cCut, csvVal, sf))
        h[sample].append(getHisto(fName, lNameTmp, lCut, csvVal, sf))
        
        b_lheNj_1_name = bNameTmp + '_lheNj_1' + '_' + chan
        c_lheNj_1_name = cNameTmp + '_lheNj_1' + '_' + chan
        l_lheNj_1_name = lNameTmp + '_lheNj_1' + '_' + chan
        bCut_lheNj_1 = bCut + ' && (lheNj == 1)' 
        cCut_lheNj_1 = cCut + ' && (lheNj == 1)' 
        lCut_lheNj_1 = lCut + ' && (lheNj == 1)'
     
        hDY[sample].append(getHisto(fName, b_lheNj_1_name, bCut_lheNj_1, csvVal, sf))
        hDY[sample].append(getHisto(fName, c_lheNj_1_name, cCut_lheNj_1, csvVal, sf))
        hDY[sample].append(getHisto(fName, l_lheNj_1_name, lCut_lheNj_1, csvVal, sf))

        b_lheNj_not1_name = bNameTmp + '_lheNj_not1' + '_' + chan
        c_lheNj_not1_name = cNameTmp + '_lheNj_not1' + '_' + chan
        l_lheNj_not1_name = lNameTmp + '_lheNj_not1' + '_' + chan
        bCut_lheNj_not1 = bCut + ' && (lheNj != 1)' 
        cCut_lheNj_not1 = cCut + ' && (lheNj != 1)' 
        lCut_lheNj_not1 = lCut + ' && (lheNj != 1)'
     
        hDY[sample].append(getHisto(fName, b_lheNj_not1_name, bCut_lheNj_not1, csvVal, sf))
        hDY[sample].append(getHisto(fName, c_lheNj_not1_name, cCut_lheNj_not1, csvVal, sf))
        hDY[sample].append(getHisto(fName, l_lheNj_not1_name, lCut_lheNj_not1, csvVal, sf))




  bNameTmp = 'bjet_vtxMass_tmp_DYcombine_' + chan
  hDY['DYcombine'].append(combineDYtemplate(hDY['DY'][0], hDY['DY'][3], hDY['DY1Jet'][0], bNameTmp))

  cNameTmp = 'cjet_vtxMass_tmp_DYcombine_' + chan
  hDY['DYcombine'].append(combineDYtemplate(hDY['DY'][1], hDY['DY'][4], hDY['DY1Jet'][1], cNameTmp))

  lNameTmp = 'ljet_vtxMass_tmp_DYcombine_' + chan
  hDY['DYcombine'].append(combineDYtemplate(hDY['DY'][2], hDY['DY'][5], hDY['DY1Jet'][2], lNameTmp))

  #@@@@@@@@@@@@@@@@@@combine data and write@@@@@@@@@@@@@@@@@@@@@

  hDatas = []
  for sample,hists in h.items():
    if sample.find('Data') != -1:
      if len(hists) > 0: hDatas.append(hists[0])

  if len(hDatas) > 0:
    hData = hDatas[0].Clone('AllJet_vtxMass_Data_' + chan)
    for i in range(1, len(hDatas)):
      hData.Add(hDatas[i])

    fOut.cd()
    hData.Write()

#@@@@@@@@@@@@@@@@write MC@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  fOut.cd()
  for sample,hists in h.items():
    for hist in hists:
      hist.Write()

  for sample,hists in hDY.items():
    for hist in hists:
      hist.Write()

fOut.Close()
  
