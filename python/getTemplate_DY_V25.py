import ROOT
from myutils.util_funcs import *

ROOT.gSystem.Load('../interface/VHbbNameSpace_h.so')

ROOT.gROOT.SetBatch(True)

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

def replace_for_list(l, st, st1):
  lOut = l[:]
  for i in range(len(l)):
    lOut[i] = l[i].replace(st, st1)
  return lOut

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Main program
##############################################################


dataPeriod = 'allData' #allData, B, C, D, E, F, G, H

csvCut = 'CSVM'
#csvCut = 'CtagT'

vtxMassType = 'vtxMass' # vtxMass, incVtxMass, vtxMassCorr_IVF

channels = ['Zmm', 'Zee']

systs = ['Central','JECUp','JECDown','JERUp','JERDown','gccUp', 'gccDown', 'gbbUp', 'gbbDown', 'puUp', 'puDown'] #always set central first
#systs = ['Central']

fOut = ROOT.TFile.Open('Test/template_' + vtxMassType + '_' + csvCut + '_pt20_noGapEle_DY_allData_allWeight_V25.root','recreate')
#fOut = ROOT.TFile.Open('Test/template_' + vtxMassType + '_' + csvCut + '_pt20_DY_' + dataPeriod + '_allWeight_V25.root','recreate')

###########################################
#Data
###########################################
#uncertainties: jec, jer, pu, gluon splitting, MC template statistic, BR B->C, ttbar, VV cross section unc, PDF.

lumi = 35862.  #all
if dataPeriod == 'B':
  lumi = 5783.74
if dataPeriod == 'C':
  lumi = 2573.399
if dataPeriod == 'D':
  lumi = 4248.384
if dataPeriod == 'E':
  lumi = 4009.132
if dataPeriod == 'F':
  lumi = 3101.618
if dataPeriod == 'G':
  lumi = 7540.488
if dataPeriod == 'H':
  lumi = 8605.689

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Cuts and corrections and scale factors
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
jetIdx = 'idxJet_passCSV_SVT[0]'
if vtxMassType == 'incVtxMass':
  jetIdx = 'idxJet_passCSV_SVT_1[0]'
if vtxMassType == 'vtxMassCorr_IVF':
  jetIdx = 'idxJet_passCSV_SVT_2[0]'

if csvCut == 'CSVM':
  jetIdx = jetIdx.replace('[0]','[1]')

if csvCut == 'CtagT':
  jetIdx = 'idxJet_passCtag_SVT[0]'

jsonCut = '(json == 1)'
triggerCut = {'Zmm': '((HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v == 1) || (HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v == 1) || (HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v == 1) || (HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v == 1))', 'Zee': '((HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v == 1) || (HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v == 1))'} 
jetSVTcut = '(Jet_' + vtxMassType + '[' + jetIdx + '] > 0)'
if vtxMassType == 'vtxMassCorr_IVF' or vtxMassType == 'incVtxMass':
  jetSVTcut = '(' + jetSVTcut + '&& (Jet_vtxCat_IVF[' + jetIdx + '] == 0))'
jetCut = '((' + jetIdx + ' >=0) && (Jet_pt[' + jetIdx + '] > 30 && abs(Jet_eta[' + jetIdx + ']) < 2.4) && ' + jetSVTcut + ')'
lepCut_Zee = '(vLeptons_new_pt[0] > 20 && vLeptons_new_pt[1] > 20 && abs(vLeptons_new_eta[0]) < 2.4 && abs(vLeptons_new_eta[1]) < 2.4 && vLeptons_new_pfRelIso03[0] < 0.25 && vLeptons_new_pfRelIso03[1] < 0.25 && (abs(vLeptons_new_etaSc[0]) < 1.4442 || abs(vLeptons_new_etaSc[0]) > 1.566) && (abs(vLeptons_new_etaSc[1]) < 1.4442 || abs(vLeptons_new_etaSc[1]) > 1.566))'
lepCut_Zmm = '(vLeptons_new_pt[0] > 20 && vLeptons_new_pt[1] > 20 && abs(vLeptons_new_eta[0]) < 2.4 && abs(vLeptons_new_eta[1]) < 2.4 && vLeptons_new_pfRelIso04[0] < 0.25 && vLeptons_new_pfRelIso04[1] < 0.25)'
zjetMassCut = '((VHbb::HVMass(vLeptons_new_pt[0],vLeptons_new_eta[0],vLeptons_new_phi[0],vLeptons_new_mass[0],vLeptons_new_pt[1],vLeptons_new_eta[1],vLeptons_new_phi[1],vLeptons_new_mass[1]) > 70) && (VHbb::HVMass(vLeptons_new_pt[0],vLeptons_new_eta[0],vLeptons_new_phi[0],vLeptons_new_mass[0],vLeptons_new_pt[1],vLeptons_new_eta[1],vLeptons_new_phi[1],vLeptons_new_mass[1]) < 110))'
metCut = '(met_pt < 40)'
#metCut = '(1)'

########################
#Scale
########################

sf_com = 'sign(genWeight)*puWeight'
#sf_chan = {'Zmm': sf_com + '*muweight_trig[0]*muweight[0]*bTagWeight_CSVT[0]', 'Zee': sf_com + '*(1)'}
sf_chan = {'Zmm': sf_com + '*muweight[0]*bTagWeight_CSVM[0]', 'Zee': sf_com + '*eleweight[0]*bTagWeight_CSVM[0]'}
if csvCut == 'CSVT': 
  sf_chan = {'Zmm': sf_com + '*muweight[0]*bTagWeight_CSVT[0]', 'Zee': sf_com + '*eleweight[0]*bTagWeight_CSVT[0]'}
if csvCut == 'CtagT': 
  sf_chan = {'Zmm': sf_com + '*muweight[0]*cTagWeight_CSVT[0]', 'Zee': sf_com + '*eleweight[0]*cTagWeight_CSVT[0]'}

#sf_chan = ['(1)', '(1)']

#########################################
#Sample
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_out = '/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V25/Zll_V25/syst_fromEOS/'
#TEMP
path_out_inc = '/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V25/Zll_inc_V25/syst_fromEOS/'
samples_data = {}
samples_data['Zmm'] = {
            'DataB':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016B-23Sep2016-v3.root',1, 1, 1],
            'DataC':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016C-23Sep2016-v1.root',1, 1, 1],
            'DataD':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016D-23Sep2016-v1.root',1, 1, 1],
            'DataE':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016E-23Sep2016-v1.root',1, 1, 1],
            'DataF':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016F-23Sep2016-v1.root',1, 1, 1],
            'DataG':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016G-23Sep2016-v1.root',1, 1, 1],
            'DataHv1':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016H-PromptReco-v1.root',1, 1, 1],
            'DataHv2':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016H-PromptReco-v2.root',1, 1, 1],
            'DataHv3':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016H-PromptReco-v3.root',1, 1, 1]
            }
samples_data['Zee'] = {
            'DataB':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016B-23Sep2016-v3.root',1, 1, 1],
            'DataC':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016C-23Sep2016-v1.root',1, 1, 1],
            'DataD':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016D-23Sep2016-v1.root',1, 1, 1],
            'DataE':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016E-23Sep2016-v1.root',1, 1, 1],
            'DataF':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016F-23Sep2016-v1.root',1, 1, 1],
            'DataG':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016G-23Sep2016-v1.root',1, 1, 1],
            'DataHv1':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016H-PromptReco-v1.root',1, 1, 1],
            'DataHv2':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016H-PromptReco-v2.root',1, 1, 1],
            'DataHv3':['ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016H-PromptReco-v3.root',1, 1, 1]
          }

print '>>>>>>>>>>>>>>>>>>>>>>>'
print samples_data['Zmm']
print samples_data['Zee']

samples_MC = {
          'DY':['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root', 0, 1921.8*3, 1], #6025.2*1.23
          'TT':['TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root', 0, 831.76, 1],
          'WW':['WW_TuneCUETP8M1_13TeV-pythia8.root', 0, 118.7, 1],
          'WZ':['WZ_TuneCUETP8M1_13TeV-pythia8.root', 0, 47.13, 1],
          'ZZ':['ZZ_TuneCUETP8M1_13TeV-pythia8.root', 0, 16.523, 1]
          }


###################################
#Now run
###################################

for sample,att in samples_MC.items():
    
  fName = path_out + att[0]
  
  if att[0] == '': fName = '' 
  
  if att[1] != 1 and fName != '':
    att[3] = scaleToLumi(fName, att[2], lumi)
    samples_MC[sample] = att

  print sample, samples_MC[sample]


for chan in channels:
#for chan in ['Zmm']:
  print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  print chan
  print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  
  #@@@@@@@@@@@final cuts@@@@@@@@@@
  Zll_cut = ""
  if chan == 'Zmm':
    Zll_cut = '((Vtype_new == 0) && ' + lepCut_Zmm + ')'
  if chan == 'Zee':
    Zll_cut = '((Vtype_new == 1) && ' + lepCut_Zee + ')'
  
  
  #hC = {'DY':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}
  #hTempC = {'DY':[]} #b, c, l 

  for syst in systs:

    print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    print syst
    print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

    jetCut_syst = ''
    jer = '(Jet_corr_JER[' + jetIdx + ']*(Jet_corr_JER[' + jetIdx + '] > 0) + (1.0)*(Jet_corr_JER[' + jetIdx + '] <= 0))'
    jerU = '(Jet_corr_JERUp[' + jetIdx + ']*(Jet_corr_JERUp[' + jetIdx + '] > 0) + (1.0)*(Jet_corr_JERUp[' + jetIdx + '] <= 0))'
    jerD = '(Jet_corr_JERDown[' + jetIdx + ']*(Jet_corr_JERDown[' + jetIdx + '] > 0) + (1.0)*(Jet_corr_JERDown[' + jetIdx + '] <= 0))'
    if syst == 'JECUp':
      jetCut_syst = jetCut.replace('Jet_pt['+jetIdx+']','Jet_rawPt['+jetIdx+']*Jet_corr_JECUp['+jetIdx+']*' + jer)
    if syst == 'JECDown':
      jetCut_syst = jetCut.replace('Jet_pt['+jetIdx+']','Jet_rawPt['+jetIdx+']*Jet_corr_JECDown['+jetIdx+']*'  + jer)
    if syst == 'JERUp':
      jetCut_syst = jetCut.replace('Jet_pt['+jetIdx+']','Jet_rawPt['+jetIdx+']*Jet_corr['+jetIdx+']*' + jerU)
    if syst == 'JERDown':
      jetCut_syst = jetCut.replace('Jet_pt['+jetIdx+']','Jet_rawPt['+jetIdx+']*Jet_corr['+jetIdx+']*' + jerD)
    
    cut = triggerCut[chan] + '&&' + Zll_cut + '&&' + zjetMassCut + '&&' + metCut
    if jetCut_syst != '': 
      cut = cut + '&&' + jetCut_syst
    else:
      cut = cut + '&&' + jetCut

    bCut = cut + '&& (abs(Jet_hadronFlavour[' + jetIdx + ']) == 5)'
    cCut = cut + '&& (abs(Jet_hadronFlavour[' + jetIdx + ']) == 4)'
    lCut = cut + '&& (abs(Jet_hadronFlavour[' + jetIdx + ']) != 4 && abs(Jet_hadronFlavour[' + jetIdx + ']) != 5)'
    
    #@@@@@@@@@scale factor@@@@@@@@@
    sf = sf_chan[chan]
    if syst == 'gccUp':
      sf1 = '(((' + jetIdx + ') >= 0 && Jet_gcc_weight[' + jetIdx + '] > 1)*(1.5) + ((' + jetIdx + ') >= 0 && Jet_gcc_weight[' + jetIdx + '] == 1))'
      sf = sf + '*' + sf1
    if syst == 'gccDown':
      sf1 = '(((' + jetIdx + ') >= 0 && Jet_gcc_weight[' + jetIdx + '] > 1)*(0.5) + ((' + jetIdx + ') >= 0 && Jet_gcc_weight[' + jetIdx + '] == 1))'
      sf = sf + '*' + sf1
    if syst == 'gbbUp':
      sf1 = '(((' + jetIdx + ') >= 0 && Jet_gbb_weight[' + jetIdx + '] > 1)*(1.5) + ((' + jetIdx + ') >= 0 && Jet_gbb_weight[' + jetIdx + '] == 1))'
      sf = sf + '*' + sf1
    if syst == 'gbbDown':
      sf1 = '(((' + jetIdx + ') >= 0 && Jet_gbb_weight[' + jetIdx + '] > 1)*(0.5) + ((' + jetIdx + ') >= 0 && Jet_gbb_weight[' + jetIdx + '] == 1))'
      sf = sf + '*' + sf1
    if syst == 'puUp':
      sf = sf.replace('puWeight','puWeightUp')
    if syst == 'puDown':
      sf = sf.replace('puWeight','puWeightDown')

    #@@@@@@@@putting togetter data and MC sample@@@@@@@@@@@@@
    samples = {}
    for l,v in samples_MC.items():
      samples[l] = v
    for l,v in samples_data[chan].items():
      if dataPeriod == 'B' and 'DataB' not in l: continue
      if dataPeriod == 'C' and 'DataC' not in l: continue
      if dataPeriod == 'D' and 'DataD' not in l: continue
      if dataPeriod == 'E' and 'DataE' not in l: continue
      if dataPeriod == 'F' and 'DataF' not in l: continue
      if dataPeriod == 'G' and 'DataG' not in l: continue
      if dataPeriod == 'H' and 'DataH' not in l: continue
      samples[l] = v

    print 'Samples are: ', samples
    print samples  

    #@@@@@@@@@@fill histograms@@@@@@@@@@@@
    h = {'DataB':[],'DataC':[], 'DataD':[], 'DataE':[], 'DataF':[], 'DataG': [], 'DataHv1': [], 'DataHv2':[], 'DataHv3':[], 'DY':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}
    hTemp = {'DY':[]} #b, c, l 
    
    h_inc = {'DataB':[],'DataC':[], 'DataD':[], 'DataE':[], 'DataF':[], 'DataG': [], 'DataHv1': [], 'DataHv2':[], 'DataHv3':[], 'DY':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}
    
    
    axis = {'x':[100, 0, 10]}
    var = ['Jet_' + vtxMassType + '[' + jetIdx + ']']
    
    for sample,att in samples.items():

      fName = path_out + att[0]
      fName_inc = path_out_inc + att[0]
    
      if att[0] == '':
        fName = ''
        fName_inc = ''

      #@@@@@@@@@@file data histogram@@@@@@@
      if sample.find('Data') != -1:
        #if syst == 'JERUp' or syst == 'JERDown': continue
        #cutData = jsonCut + '&&' + triggerCut[chan] + '&&' + cut.replace('*Jet_corr_JER['+jetIdx+']','')
        #hist_name = 'AllJet_vtxMass_' + sample + '_' + chan
        #if syst != 'Central':
        #  hist_name = 'AllJet_vtxMass_' + sample + '_' + syst + '_' + chan

        if syst != 'Central': continue # do not make syst shifts for data
        
        hist_name = 'AllJet_vtxMass_' + sample + '_' + chan
        cutData = jsonCut + '&&' + triggerCut[chan] + '&&' + cut
        h[sample].append(getHisto(fName, hist_name, axis, var, cutData))
        
        varTmp = replace_for_list(var,jetIdx,'0')
        cutDataTmp = cutData.replace(jetSVTcut,'(1)').replace(jetIdx,'0')
        print jetSVTcut
        print '>>>>>>>>>>>>>>>>>>Cut for inclusive is:>>>>>>>>>>>>>>>>>'
        print jetSVTcut
        print cutDataTmp
        h_inc[sample].append(getHisto(fName_inc, hist_name + '_inc', axis, varTmp, cutDataTmp))
      
        continue
    
    
      bName  = 'bjet_vtxMass_' + sample + '_' + chan
      cName  = 'cjet_vtxMass_' + sample + '_' + chan
      lName  = 'ljet_vtxMass_' + sample + '_' + chan
      if syst != 'Central':
        bName  = 'bjet_vtxMass_' + sample + '_' + syst + '_' + chan
        cName  = 'cjet_vtxMass_' + sample + '_' + syst + '_' + chan
        lName  = 'ljet_vtxMass_' + sample + '_' + syst + '_' + chan

      h[sample].append(getHisto(fName, bName, axis, var, bCut, str(samples[sample][3]) + '*' + sf))
      h[sample].append(getHisto(fName, cName, axis, var, cCut, str(samples[sample][3]) + '*' + sf))
      h[sample].append(getHisto(fName, lName, axis, var, lCut, str(samples[sample][3]) + '*' + sf))
      
      if syst == 'Central':
        sfTmp = sf.replace('*bTagWeight_CSVT[0]','').replace('*bTagWeight_CSVM[0]','').replace('*cTagWeight_CSVM[0]','').replace('*cTagWeight_CSVM[1]','') #use for inclusive, no tagging weights
        varTmp = replace_for_list(var,jetIdx,'0')
        bCutTmp = bCut.replace(jetSVTcut,'(1)').replace(jetIdx,'0')
        cCutTmp = cCut.replace(jetSVTcut,'(1)').replace(jetIdx,'0')
        lCutTmp = lCut.replace(jetSVTcut,'(1)').replace(jetIdx,'0')
        h_inc[sample].append(getHisto(fName_inc, bName + '_inc', axis, varTmp, bCutTmp, str(samples[sample][3]) + '*' + sfTmp))
        h_inc[sample].append(getHisto(fName_inc, cName + '_inc', axis, varTmp, cCutTmp, str(samples[sample][3]) + '*' + sfTmp))
        h_inc[sample].append(getHisto(fName_inc, lName + '_inc', axis, varTmp, lCutTmp, str(samples[sample][3]) + '*' + sfTmp))

      if (sample == 'DY'):
      
          bNameTmp = 'bjet_vtxMass_temp_' + sample + '_' + chan
          cNameTmp = 'cjet_vtxMass_temp_' + sample + '_' + chan
          lNameTmp = 'ljet_vtxMass_temp_' + sample + '_' + chan
          if syst != 'Central':
            bNameTmp = 'bjet_vtxMass_temp_' + sample + '_' + syst + '_' + chan
            cNameTmp = 'cjet_vtxMass_temp_' + sample + '_' + syst + '_' + chan
            lNameTmp = 'ljet_vtxMass_temp_' + sample + '_' + syst + '_' + chan
          
          hTemp[sample].append(getHisto(fName, bNameTmp, axis, var, bCut, sf))
          hTemp[sample].append(getHisto(fName, cNameTmp, axis, var, cCut, sf))
          hTemp[sample].append(getHisto(fName, lNameTmp, axis, var, lCut, sf))

    #####################Clone central value to use in elimination of large unc bins###########
    #if syst == 'Central':
    #  for sample,hists in h.items():
    #    for hist in hists:
    #      if 'Data' not in sample: hC[sample].append(hist.Clone(hist.GetName() + '_cloneForKillLargeUnc'))
    #  for sample,hists in hTemp.items():
    #    for hist in hists:
    #      if 'Data' not in sample: hTempC[sample].append(hist.Clone(hist.GetName() + '_cloneForKillLargeUnc'))

    #@@@@@@@@@@@@@@@@@@combine data and write@@@@@@@@@@@@@@@@@@@@@
    
    fOut.cd()
    
    if syst == 'Central': #there is no syst for data
      hDatas = []
      for sample,hists in h.items():
        if sample.find('Data') != -1:
          if len(hists) > 0: hDatas.append(hists[0])

      if len(hDatas) > 0:
        hist_name = 'AllJet_vtxMass_Data_' + chan
        hData = hDatas[0].Clone(hist_name)
        for i in range(1, len(hDatas)):
          hData.Add(hDatas[i])
        
        hData.Write()
      
      hDatas = []
      for sample,hists in h_inc.items():
        if sample.find('Data') != -1:
          if len(hists) > 0: hDatas.append(hists[0])

      if len(hDatas) > 0:
        hist_name = 'AllJet_vtxMass_Data_' + chan + '_inc'
        hData = hDatas[0].Clone(hist_name)
        for i in range(1, len(hDatas)):
          hData.Add(hDatas[i])

        hData.Write()

  #@@@@@@@@@@@@@@@@remove high error bins in JEC and JER and write MC@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   
    for sample,hists in h.items():
      #hTmps = []
      #if 'Data' not in sample: hTmps = hC[sample]
      for i in range(len(hists)):
        #if 'Data' not in sample:
        #  if 'JEC' in syst or 'JER' in syst:
        #    removeHighUncBins([hists[i],hTmps[i]])
        
        hists[i].Write()

    for sample,hists in hTemp.items():
      #hTmps = []
      #if 'Data' not in sample: hTmps = hTempC[sample]
      for i in range(len(hists)):
        #if 'Data' not in sample:
        #  if 'JEC' in syst or 'JER' in syst:
        #    removeHighUncBins([hists[i],hTmps[i]])
        
        hists[i].Write()
    
    if syst == 'Central':
      for sample,hists in h_inc.items():
        #hTmps = []
        #if 'Data' not in sample: hTmps = hC[sample]
        for i in range(len(hists)):
          #if 'Data' not in sample:
          #  if 'JEC' in syst or 'JER' in syst:
          #    removeHighUncBins([hists[i],hTmps[i]])
        
          hists[i].Write()


fOut.Close()
