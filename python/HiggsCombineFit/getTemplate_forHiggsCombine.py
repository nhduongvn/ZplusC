from ROOT import *
from math import sqrt
import sys
sys.path.append('/uscms_data/d3/duong/CMSSW/CMSSW_7_6_5/src/ZplusC/python/')
from myutils.util_funcs import *

import array

gROOT.SetBatch(False)

def customBin(h, xDiv):
  hOut = TH1D(h.GetName() + '_customBinning','',len(xDiv)-1,array.array('f',xDiv))
  for i in range(1,hOut.GetNbinsX()+1):
    eL = hOut.GetBinLowEdge(i) - 0.000001
    eH = hOut.GetBinLowEdge(i) + hOut.GetBinWidth(i) + 0.000001
    
    contS = 0
    errS = 0
    for j in range(1,h.GetNbinsX()+1):
      eL1 = h.GetBinLowEdge(j)
      eH1 = eL1 + h.GetBinWidth(j)
      if eL1 >= eL and eH1 <= eH:
        contS += h.GetBinContent(j)
        errS += h.GetBinError(j)*h.GetBinError(j)
    errS = TMath.Sqrt(errS)
    
    hOut.SetBinContent(i,contS)
    hOut.SetBinError(i,errS)
  
  hOut.SetBinError(0,h.GetBinError(0))
  hOut.SetBinContent(0,h.GetBinContent(0))
  hOut.SetBinError(hOut.GetNbinsX()+1,h.GetBinError(h.GetNbinsX()+1))
  hOut.SetBinContent(hOut.GetNbinsX()+1,h.GetBinContent(h.GetNbinsX()+1))

  return hOut


def makeName(pre,chan,syst):
  if syst != 'Central':
    if chan == '': return pre + syst
    return pre + syst + '_' + chan
  if chan != '':
    return pre + '_' + chan
  return pre


def getHisto(fIn, name, xRange=[-1,-1], xDiv=[-1,-1]):
  histIn = fIn.Get(name)
  #make a copy of histogram with the same bin width but with different range
  print '>>>>>>>Hist name: ', name
  hist = None
  if xRange[0] != -1 and xRange[1] != -1:
    nBin = int((xRange[1] - xRange[0])/histIn.GetBinWidth(1)) 
    print '>>>nBin: ', nBin
    hist = TH1D(name + '_clone', '', nBin, xRange[0], xRange[1])
    for i in range(1,hist.GetNbinsX() + 1):
      #get bin number
      iBin = histIn.FindFixBin(hist.GetBinLowEdge(i))
      hist.SetBinContent(i, histIn.GetBinContent(iBin))
      hist.SetBinError(i, histIn.GetBinError(iBin))
    
    #get overflow
    ofBin = histIn.FindFixBin(xRange[1])
    nOf = 0
    nErrOf = 0
    for i in range(ofBin, histIn.GetNbinsX() + 2):
      nOf += histIn.GetBinContent(i)
      nErrOf += histIn.GetBinError(i)*histIn.GetBinError(i)
    nErrOf = sqrt(nErrOf)
    hist.SetBinContent(hist.GetNbinsX() + 1, nOf)
    hist.SetBinError(hist.GetNbinsX() + 1, nErrOf)

    #get underflow
    ufBin = histIn.FindFixBin(xRange[0]-histIn.GetBinWidth(1))
    nUf = 0
    nErrUf = 0
    for i in range(0, ufBin + 1):
      nUf += histIn.GetBinContent(i)
      nErrUf += histIn.GetBinError(i)*histIn.GetBinError(i)
    nErrUf = sqrt(nErrUf)
    hist.SetBinContent(0, nUf)
    hist.SetBinError(0, nErrUf)
  
  else: hist = histIn.Clone(name + '_clone')
  
  hist1 = None
  if xDiv[0] != -1 and xDiv[1] != -1:
    hist1 = customBin(hist,xDiv)
    return hist1


  return hist

def makeline(words):
  tmp = '{0:30}'.format(words[0]) 
  for i in range(1,len(words)):
    tmp = tmp + '{0:20}'.format(words[i])
  return tmp
 
def makeShapeStatUnc(h,fOut):
  #loop over bin to find bins with 10% uncertainties GetBinError(i)/GetBinContent(i) > 0.01
  nBin = h.GetNbinsX()
  selBin = []
  selBinHigh = []
  highUncCutOff = 1.0 #100%
  for i in range(1,nBin + 1):
    #if h.GetBinError(i)/h.GetBinContent(i) > 0.01 and h.GetBinError(i)/h.GetBinContent(i) < 0.5:
    if h.GetBinContent(i) > 0 and h.GetBinError(i)/h.GetBinContent(i) > 0.0:
      selBin.append(i)
    if h.GetBinContent(i) > 0 and h.GetBinError(i)/h.GetBinContent(i) > highUncCutOff:
      selBinHigh.append(i)
      print 'Warning: high uncertainty ', h.GetName(), ' ', i, ' ', h.GetBinContent(i), ' ', h.GetBinError(i)
    if h.GetBinContent(i) <= 0:
      print 'Warning: negative bin content ', h.GetName(), ' ', i, ' ', h.GetBinContent(i), ' ', h.GetBinError(i)
    
  hOut = []
  for iBin in selBin:
    name = h.GetName() + '' + h.GetName() + '_stat_bin_' + str(iBin) + 'Up' 
    hTmp = h.Clone(name)
    if iBin in selBinHigh:
      hTmp.SetBinContent(iBin,hTmp.GetBinContent(iBin) + hTmp.GetBinContent(iBin)*highUncCutOff)
    else:
      hTmp.SetBinContent(iBin,hTmp.GetBinContent(iBin) + hTmp.GetBinError(iBin))
    hOut.append(hTmp)

    name = h.GetName() + '' + h.GetName() + '_stat_bin_' + str(iBin) + 'Down' 
    hTmp = h.Clone(name)
    delta = hTmp.GetBinContent(iBin) - hTmp.GetBinError(iBin)
    if iBin in selBinHigh:
      delta = hTmp.GetBinContent(iBin) - hTmp.GetBinContent(iBin)*highUncCutOff
    if delta > 0: hTmp.SetBinContent(iBin,delta)
    else:
      hTmp.SetBinContent(iBin,0)
      print 'Negative value found when fluctuate down bin statistical error, set to 0: ', h.GetName(), ' ', h.GetBinContent(iBin), ' ', h.GetBinError(iBin)
    hOut.append(hTmp)
  fOut.cd()
  for hTmp in hOut: hTmp.Write()

def write_mc_stat_line(proc, uncs, f_dc):
  for unc in uncs:
    ws = []
    if proc == 'Zb':
      ws = [unc,'shape','1.0','-','-','-','-']
    if proc == 'Zc':
      ws = [unc,'shape','-','1.0','-','-','-']
    if proc == 'Zl':
      ws = [unc,'shape','-','-','1.0','-','-']
    if proc == 'TT':
      ws = [unc,'shape','-','-','-','1.0','-']
    if proc == 'VV':
      ws = [unc,'shape','-','-','-','-','1.0']
    st = makeline(ws)
    f_dc.write(st + '\n')

def write_data_card(dc_name,fIn,systs, useDDbjet=True, useDDljet=True):
  f_dc = open(dc_name + '.txt','w')
  f_dc.write('imax 1 number of channels\n')
  f_dc.write('jmax 4 number of backgrounds (\'*\' = automatic)\n')
  f_dc.write('kmax * number of nuisance parameters (sources of systematical uncertainties)\n')
  st = 'shapes * * ' + fIn.GetName() + ' $PROCESS $PROCESS$SYSTEMATIC\n\n'
  f_dc.write(st)
  
  obs = (fIn.Get('data_obs')).Integral()
  ws = ['observation',str(obs)]
  st = makeline(ws)
  f_dc.write(st + '\n\n')
  
  ws = ['bin','SR']
  st = makeline(ws)
  f_dc.write(st + '\n\n')
  
  ws = ['bin','','SR','SR','SR','SR','SR']
  st = makeline(ws)
  f_dc.write(st + '\n')

  ws = ['process','','Zb','Zc','Zl','TT','VV']
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  ws = ['process','','-2','-1','0','1','2']
  st = makeline(ws)
  f_dc.write(st + '\n')
 
  #write rate
  zb = str((fIn.Get('Zb')).Integral())
  zc = str((fIn.Get('Zc')).Integral())
  zl = str((fIn.Get('Zl')).Integral())
  TT = str((fIn.Get('TT')).Integral())
  VV = str((fIn.Get('VV')).Integral())
  ws = ['rate','',zb,zc,zl,TT,VV]
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  #lumi
  #ws = ['lumi','lnN','1.026','1.026','1.026','1.026','1.026']
  #st = makeline(ws)
  #f_dc.write(st + '\n')

  #TTbar_bkgr_xSec
  ws = ['TTbar_bkgr_xSec','lnN','-','-','-','1.053','-']
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  #Shape sys
  for syst in systs:
    if syst != 'Central':
      ws = [syst,'shape','1.0','1.0','1.0','1.0','1.0']
      if syst == 'gbb':
        ws = [syst,'shape','1.0','-','-','1.0','1.0']
      if syst == 'gcc':
        ws = [syst,'shape','-','1.0','-','1.0','1.0']
      if syst == 'pu':
        ws = [syst,'shape','-','1.0','-','1.0','1.0']
      if syst == 'JER':
        if useDDbjet and useDDljet:
          ws = [syst,'shape','-','1.0','-','1.0','1.0']
        if useDDbjet and not useDDljet:
          ws = [syst,'shape','-','1.0','1.0','1.0','1.0']
        if not useDDbjet and useDDljet:
          ws = [syst,'shape','1.0','1.0','-','1.0','1.0']
      st = makeline(ws)
      f_dc.write(st + '\n')

  #Statistical of MC shape
  uncNames = {'Zb':[],'Zc':[],'Zl':[],'TT':[],'VV':[]}
  for k,v in uncNames.items(): 
    for key in fIn.GetListOfKeys():
      hist = key.ReadObj()
      histName = hist.GetName()
      if k in histName and 'stat_bin' in histName and 'Up' in histName:
        histName = histName.replace('Up','')
        histName = histName.replace(k,'')
        histName = k + histName
        uncNames[k].append(histName)
    write_mc_stat_line(k, uncNames[k], f_dc)

  
  #define unc. group
  MC_stat_systs = ''
  for k,vs in uncNames.items():
    for v in vs:
      MC_stat_systs += ' ' + v
  
  #ws = ['lumi_gr','group = lumi']
  #st = makeline(ws)
  #f_dc.write(st + '\n')
  
  ws = ['MC_stat_bin','group = ' + MC_stat_systs]
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  ws = ['JEC_gr','group = JEC']
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  ws = ['JER_gr','group = JER']
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  ws = ['gcc_gr','group = gcc']
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  ws = ['gbb_gr','group = gbb']
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  ws = ['pu_gr','group = pu']
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  ws = ['Bkgr_xSec','group = TTbar_bkgr_xSec']
  st = makeline(ws)
  f_dc.write(st + '\n')

  #ws = ['Exp','group = lumi JEC JER gcc gbb TTbar_bkgr_xSec']
  ws = ['Exp','group = JEC JER gcc gbb pu TTbar_bkgr_xSec']
  st = makeline(ws)
  f_dc.write(st + '\n')
  
  #ws = ['Exp','group = JEC JER']
  #st = makeline(ws)
  #f_dc.write(st + '\n')
  
  ws = ['Exp','group += ' + MC_stat_systs]
  st = makeline(ws)
  f_dc.write(st + '\n')
    
def removeHighUncBins(hIns): #[0]: syst, [1]: central
  hR = hIns[0].Clone('hRtmp')
  hR.Divide(hIns[1])
  for i in range(hR.GetNbinsX() + 1):
    if abs(hR.GetBinContent(i) - 1.0) > 0.2:
      hIns[0].SetBinContent(i, hIns[1].GetBinContent(i))
 


##############################################
# Main program
##############################################

#fIn = TFile.Open('Test/testTemplate_CSVM_G_V25.root','read')
#fOut = TFile.Open('Test/forRoofit_do_not_use.root','recreate')

#fIn = TFile.Open('Test/testTemplate_inverseCSVM_G_useLeadingJetPt_V25.root','read')
#fOut = TFile.Open('Test/forRoofit_inverseCSVM_useLeadingJetPt.root','recreate')
tagName = 'CSVM' #CSVM CtagT
useDDbjet = True 
useDDljet = False 
useDY_nlo = False 

fIn_name = '../Test/template_vtxMass_'+tagName+'_pt20_noGapEle_DY_allData_allWeight_V25.root'
if useDY_nlo:
    fIn_name = fIn_name.replace('DY','DY_nlo_1')

fIn_emu_name = '../Test/template_vtxMass_'+tagName+'_lep_pt_25_emu_allData_V25.root'
fIn_ttsemi_name = '../Test/template_vtxMass_'+tagName+'_ttsemi_allData_V25.root' #not use currently, TODO: data-driven cjet msv

fIn_Wjet_name = 'template_CSVv2_deepCSV_vtxMass_'+tagName+'_iso0p05_Wjet_allData_V25.root'
if tagName == 'CtagT':
  fIn_Wjet_name = 'template_deepCSV_vtxMass_'+tagName+'_iso0p05_Wjet_allData_V25.root'

fIn = TFile.Open(fIn_name,'read')
fIn_emu = TFile.Open(fIn_emu_name,'read')
fIn_ttsemi = TFile.Open(fIn_ttsemi_name,'read')
fIn_Wjet = TFile.Open('../Test/'+fIn_Wjet_name) #need to rerun to get JEC up

#fIn = TFile.Open('Test/forHiggsCombine_'+tagName+'_useLeadingJetPt.root','read')
#write_data_card('test',fIn)

temp_name = "vtxMass"
xRange = [0,6] #[-1,-1] for unchange range
#xDiv = [0.0, 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 2.6 , 2.8 , 3.0 , 3.2 , 3.4 , 3.6 , 3.8 , 4.0 , 4.2 , 4.4 , 4.6 , 4.8 , 5.0, 6.0] #[-1,-1] for unchange
#xDiv = [0.0, 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 2.6 , 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 6.0] #[-1,-1] for unchange
xDiv = [0.0, 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 2.6 , 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.4, 6.0] #[-1,-1] for unchange
nRebin = 1 #Set nRebin = 1 when using custom binning
#xDiv = [-1,-1]
#nRebin = 2
dc_name_pre = 'test_vtxMass_'+tagName+'_pt20_noGapEle_allData_binStatThresh_0'
fOut_name_pre = '../Test/forHiggsCombine_vtxMass_'+tagName+'_pt20_noGapEle_DY_allData_allWeight_V25_binStatThresh_0'

if useDY_nlo:
    dc_name_pre = dc_name_pre.replace('_allData_','_DY_nlo_1_allData_')
    fOut_name_pre = fOut_name_pre.replace('_DY_','_DY_nlo_1_')


##########################
#start to do something
##########################
print '#######################################'
print "Input histogram file: "
print "Data, bkgr and MC template(s): ", fIn.GetName()
if useDDbjet:
  print "b-jet data-driven template: ", fIn_emu.GetName()
if useDDljet:
  print "l-jet data-driven template: ", fIn_Wjet.GetName()
print '#######################################'

if useDDbjet:
  dc_name_pre = dc_name_pre + '_useDDbjet'
if useDDljet:
  dc_name_pre = dc_name_pre + '_useDDljet'

if xDiv[0] != -1 and xDiv[1] != -1:
  dc_name_pre = dc_name_pre + '_customBinning'
else:
  dc_name_pre = dc_name_pre + '_nRebin_' + str(nRebin)

if useDDbjet:
  fOut_name_pre = fOut_name_pre + '_useDDbjet'
if useDDljet:
  fOut_name_pre = fOut_name_pre + '_useDDljet'

if xDiv[0] != -1 and xDiv[1] != -1:
  fOut_name_pre = fOut_name_pre + '_customBinning'
else:
  fOut_name_pre = fOut_name_pre + '_nRebin_' + str(nRebin)


#fOut_name_pre = '../Test/forHiggsCombine_vtxMassCorr_IVF_'+tagName+'_DY_allData_allWeight_V25'

chans = ['Zee','Zmm']

systs = ['Central','JECUp', 'JECDown', 'JERUp', 'JERDown', 'gccUp', 'gccDown', 'gbbUp', 'gbbDown', 'puUp', 'puDown'] #always put central first
systs_forDc_writting = ['JEC', 'JER', 'gcc', 'gbb', 'pu']


for chan in chans:
  
  dc_name = dc_name_pre + '_' + chan
  fOut = TFile.Open(fOut_name_pre + '_' + chan + '.root','recreate')
  
  hC = {'DY':[], 'TT':[], 'VV':[]} #use in canceling high stat unc bin

  for syst in systs:  
    
    ######################
    #Get histogram to the list
    ######################
    h = {'Data':[], 'DY':[], 'DYtemp':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}
    h_inc = {'Data':[], 'DY':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]} #events before b-tagging

    for sample,hist in h.items():
      
      bName  = 'bjet_' + temp_name + '_' + sample + '_' + chan
      cName  = 'cjet_' + temp_name + '_' + sample + '_' + chan
      lName  = 'ljet_' + temp_name + '_' + sample + '_' + chan
      if syst != 'Central':
        bName  = 'bjet_' + temp_name + '_' + sample + '_' + syst + '_' + chan
        cName  = 'cjet_' + temp_name + '_' + sample + '_' + syst + '_' + chan
        lName  = 'ljet_' + temp_name + '_' + sample + '_' + syst + '_' + chan


      if sample == 'Data':
        if syst == 'Central':
          h[sample].append(getHisto(fIn, 'AllJet_' + temp_name + '_Data_' + chan,xRange))
          h_inc[sample].append(getHisto(fIn, 'AllJet_' + temp_name + '_Data_' + chan + '_inc',xRange))
        
        continue
      
      if sample != 'DYtemp':
        
        #!!!!!!!!!!do not swap the order b -> c -> l

        #get b-jet
        if not useDDbjet or sample != 'DY':
          h[sample].append(getHisto(fIn, bName, xRange))
        else:
          #the shape taken from DD and norm taken from MC DY
          hDY_bjetTmp = None
          if 'gcc' in syst or 'gbb' in syst or 'pu' in syst or 'JER' in syst: #use central for gcc, gbb, pu and JER unc consider no gcc, gbb, pu, JER effect on b-jet template from data.
            hDY_bjetTmp = getHisto(fIn, 'bjet_' + temp_name + '_' + sample + '_' + chan, xRange)
          else:
            hDY_bjetTmp = getHisto(fIn, bName, xRange)
          hTmp = SVT_mass_correction(getHisto(fIn_emu, 'AllJet_vtxMass_Data_bkgrSub_emu_px', xRange), syst) #other systs. rather than gbb and JEC are corrected to central 
          hTmp.Scale(hDY_bjetTmp.Integral()/hTmp.Integral())
          h[sample].append(hTmp)
        
        #get c-jet
        h[sample].append(getHisto(fIn, cName, xRange))
        
        #get l-jet
        if not useDDljet or sample != 'DY':
          h[sample].append(getHisto(fIn, lName, xRange))
        else:
          #the shape taken from DD and norm taken from MC DY
          hDY_ljetTmp = None
          if 'gcc' in syst or 'gbb' in syst or 'pu' in syst or 'JER' in syst: #use central for gcc, gbb, pu and JER unc; consider no gcc, gbb, pu, JER effect on l-jet templates from data.
            hDY_ljetTmp = getHisto(fIn, 'ljet_' + temp_name + '_' + sample + '_' + chan, xRange)
          else:
            hDY_ljetTmp = getHisto(fIn, lName, xRange) #the lName has syst effect
          
          name = 'AllJet_vtxMass_Data_bkgrSub'
          if 'JEC' in syst: #other syst use central value
            name = 'AllJet_vtxMass_Data_bkgrSub_' + syst

          hTmp = getHisto(fIn_Wjet, name + '_Wjet_px', xRange) 
          hTmp.Scale(hDY_ljetTmp.Integral()/hTmp.Integral())
          h[sample].append(hTmp)
        
        if syst == 'Central':
          h_inc[sample].append(getHisto(fIn, bName + '_inc', xRange))
          h_inc[sample].append(getHisto(fIn, cName + '_inc', xRange))
          h_inc[sample].append(getHisto(fIn, lName + '_inc', xRange))

      if sample == 'DYtemp':
        bNameTmp = 'bjet_' + temp_name + '_temp_DY_' + chan
        cNameTmp = 'cjet_' + temp_name + '_temp_DY_' + chan
        lNameTmp = 'ljet_' + temp_name + '_temp_DY_' + chan
        if syst != 'Central':
            bNameTmp = 'bjet_' + temp_name + '_temp_DY_' + syst + '_' + chan
            cNameTmp = 'cjet_' + temp_name + '_temp_DY_' + syst + '_' + chan
            lNameTmp = 'ljet_' + temp_name + '_temp_DY_' + syst + '_' + chan

        h[sample].append(getHisto(fIn, bNameTmp, xRange))
        h[sample].append(getHisto(fIn, cNameTmp, xRange))
        h[sample].append(getHisto(fIn, lNameTmp, xRange))

    #################
    #Adding TT and VV background
    #################
    hBkgr_TT = None
    hBkgr_VV = None
    first_TT = True
    first_VV = True
    nMC = 0
    for sample,hists in h.items():
      if sample.find('Data') != -1 or sample.find('DY') != -1 or sample.find('DYtemp') != -1:
        continue
      print '========================='
      print sample
      for hist in hists:
        print hist.GetName(), ' ', hist.Integral(), ' ', hist.GetBinContent(hist.GetNbinsX()+1)
        nMC = nMC + hist.Integral() #checking
        if sample.find('TT') != -1:
          if first_TT:
            nameTmp = makeName('TT','',syst) #chan = '' to consistent with datacard $PROCESS$SYSTEMATIC

            hBkgr_TT = hist.Clone(nameTmp)
            first_TT = False
          else:
            hBkgr_TT.Add(hist)
        if sample.find('WW') != -1 or sample.find('WZ') != -1 or sample.find('ZZ') != -1:
          if first_VV:
            nameTmp = makeName('VV','',syst) #chan = '' to consistent with datacard $PROCESS$SYSTEMATIC
            hBkgr_VV = hist.Clone(nameTmp)
            first_VV = False
          else:
            hBkgr_VV.Add(hist)
    
    nMC_inc = 0
    hBkgr_TT_inc = None
    hBkgr_VV_inc = None
    first_TT_inc = True
    first_VV_inc = True
    if syst == 'Central': 
      for sample,hists in h_inc.items():
        if sample.find('Data') != -1 or sample.find('DY') != -1 or sample.find('DYtemp') != -1:
          continue
        print '========================='
        print sample
        for hist in hists:
          print hist.GetName(), ' ', hist.Integral(), ' ', hist.GetBinContent(hist.GetNbinsX()+1)
          nMC_inc = nMC_inc + hist.Integral() #checking
          if sample.find('TT') != -1:
            if first_TT_inc:
              nameTmp = makeName('TT_inc','',syst) #chan = '' to consistent with datacard $PROCESS$SYSTEMATIC
              hBkgr_TT_inc = hist.Clone(nameTmp)
              first_TT_inc = False
            else:
              hBkgr_TT_inc.Add(hist)
          if sample.find('WW') != -1 or sample.find('WZ') != -1 or sample.find('ZZ') != -1:
            if first_VV_inc:
              nameTmp = makeName('VV_inc','',syst) #chan = '' to consistent with datacard $PROCESS$SYSTEMATIC
              hBkgr_VV_inc = hist.Clone(nameTmp)
              first_VV_inc = False
            else:
              hBkgr_VV_inc.Add(hist)
    
    ##########################################
    #Add DY to find total MC
    ##########################################
    nameTmp = makeName('AllJet_MC_'+ temp_name,chan,syst)
    hMC = hBkgr_TT.Clone(nameTmp)
    hMC.Add(hBkgr_VV)
    for hist in h['DY']:
      nMC = nMC + hist.Integral() #checking
      hMC.Add(hist)
    
    #########################################
    #Clone data
    #########################################
    hData = None
    hData_inc = None
    if syst == 'Central':
      hData = (h['Data'][0]).Clone("data_obs")
      hData_inc = (h_inc['Data'][0]).Clone("data_obs_inc")
      
    ###############################################
    #Scale MC by data/MC ratio to estimate background (bypass luminosity scaling)
    ###############################################
    #scale_norm = hData.Integral()/hMC.Integral()

    #print '========================='
    #print 'Scale for norm is: ', hData.Integral(), ' ', hMC.Integral(), ' ', scale_norm
    #print 'Will scale background by ', scale_norm, ' in background subtraction'

    #######################################
    #clone histograms for Roofit testing
    #######################################
    nameTmp = makeName('Zb','',syst)
    hZb = h['DY'][0].Clone(nameTmp)
    nameTmp = makeName('Zc','',syst)
    hZc = h['DY'][1].Clone(nameTmp)
    nameTmp = makeName('Zl','',syst)
    hZl = h['DY'][2].Clone(nameTmp)
    
    print '==============================='
    print hZb.GetName(), ' ', hZb.Integral(), ' ', hZb.GetBinContent(hZb.GetNbinsX()+1)
    print hZc.GetName(), ' ', hZc.Integral(), ' ', hZc.GetBinContent(hZc.GetNbinsX()+1)
    print hZl.GetName(), ' ', hZl.Integral(), ' ', hZl.GetBinContent(hZl.GetNbinsX()+1)
    

    #========print results===========
    nBkgr = []
    nBkgr_inc = []
    if syst == 'Central':
      nData = hData.Integral()
      nBkgr.append(hZb.Integral())
      nBkgr.append(hZc.Integral())
      nBkgr.append(hZl.Integral())
      nBkgr.append(hBkgr_TT.Integral())
      nBkgr.append(hBkgr_VV.Integral())
      print '>>>>>>>>>>>>>>>>>>>>>>'
      print '>>>>>', chan, '>>>>>>>'
      print 'Data: ', hData.Integral()  
      print 'DY + bjet: ', nBkgr[0] 
      print 'DY + cjet: ', nBkgr[1] 
      print 'DY + ljet: ', nBkgr[2]
      print 'TT: ', nBkgr[3]
      print 'VV: ', nBkgr[4]
      totBkgr = nBkgr[0] + nBkgr[1] + nBkgr[2] + nBkgr[3] + nBkgr[4]
      print 'Total MC: ', totBkgr 
      print 'Data/MC: ', nData/totBkgr
      
      nData_inc = hData_inc.Integral()
      nBkgr_inc.append(h_inc['DY'][0].Integral())
      nBkgr_inc.append(h_inc['DY'][1].Integral())
      nBkgr_inc.append(h_inc['DY'][2].Integral())
      nBkgr_inc.append(hBkgr_TT.Integral())
      nBkgr_inc.append(hBkgr_VV.Integral())
      print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      print '>>>>>', chan, ', inclusive >>>>>>>'
      print 'Data: ', hData_inc.Integral()  
      print 'DY + bjet: ', nBkgr_inc[0] 
      print 'DY + cjet: ', nBkgr_inc[1] 
      print 'DY + ljet: ', nBkgr_inc[2]
      print 'TT: ', nBkgr_inc[3]
      print 'VV: ', nBkgr_inc[4]
      totBkgr_inc = nBkgr_inc[0] + nBkgr_inc[1] + nBkgr_inc[2] + nBkgr_inc[3] + nBkgr_inc[4]
      print 'Total MC: ', totBkgr_inc 
      print 'Data/MC: ', nData_inc/totBkgr_inc

    ##############################
    #Write histogram
    ##############################

    fOut.cd()
    
    if syst == 'Central': 
      if xDiv[0] == -1 and xDiv[1] == -1: 
        hData.Rebin(nRebin)
      else:
        nameTmp = hData.GetName()
        hData = customBin(hData, xDiv)
        hData.SetName(nameTmp)
       
      hData.Write()
    
    #do something on Zb
    if xDiv[0] == -1 and xDiv[1] == -1: 
      hZb.Rebin(nRebin)
    else:
      nameTmp = hZb.GetName()
      hZb = customBin(hZb, xDiv)
      hZb.SetName(nameTmp)
    
    if syst == 'Central':
      makeShapeStatUnc(hZb,fOut)
      hC['DY'].append(hZb.Clone(hZb.GetName() + '_cloneToKillHighUnc'))

    if syst != 'Central':
      if 'JEC' in syst or 'JER' in syst:
        removeHighUncBins([hZb,hC['DY'][0]])

    hZb.Write()

    #do something on Zc
    if xDiv[0] == -1 and xDiv[1] == -1: 
      hZc.Rebin(nRebin)
    else:
      nameTmp = hZc.GetName()
      hZc = customBin(hZc, xDiv)
      hZc.SetName(nameTmp)

    if syst == 'Central':
      makeShapeStatUnc(hZc,fOut)
      hC['DY'].append(hZc.Clone(hZc.GetName() + '_cloneToKillHighUnc'))
    if syst != 'Central':
      if 'JEC' in syst or 'JER' in syst:
        removeHighUncBins([hZc,hC['DY'][1]])
    
    hZc.Write()

    #do something on Zl
    if xDiv[0] == -1 and xDiv[1] == -1: 
      hZl.Rebin(nRebin)
    else:
      nameTmp = hZl.GetName()
      hZl = customBin(hZl, xDiv)
      hZl.SetName(nameTmp)
    
    if syst == 'Central':
      makeShapeStatUnc(hZl,fOut)
      hC['DY'].append(hZl.Clone(hZl.GetName() + '_cloneToKillHighUnc'))
    if syst != 'Central':
      if 'JEC' in syst or 'JER' in syst:
        removeHighUncBins([hZl,hC['DY'][2]])

    hZl.Write()
    
    #do something on TT
    if xDiv[0] == -1 and xDiv[1] == -1: 
      hBkgr_TT.Rebin(nRebin)
    else:
      nameTmp = hBkgr_TT.GetName()
      hBkgr_TT = customBin(hBkgr_TT, xDiv)
      hBkgr_TT.SetName(nameTmp)

    if syst == 'Central':
      makeShapeStatUnc(hBkgr_TT,fOut)
      hC['TT'].append(hBkgr_TT.Clone(hBkgr_TT.GetName() + '_cloneToKillHighUnc'))
    if syst != 'Central':
      if 'JEC' in syst or 'JER' in syst:
        removeHighUncBins([hBkgr_TT,hC['TT'][0]])
    hBkgr_TT.Write()

    #do something on VV
    if xDiv[0] == -1 and xDiv[1] == -1: 
      hBkgr_VV.Rebin(nRebin)
    else:
      nameTmp = hBkgr_VV.GetName()
      hBkgr_VV = customBin(hBkgr_VV, xDiv)
      hBkgr_VV.SetName(nameTmp)

    if syst == 'Central':
      makeShapeStatUnc(hBkgr_VV,fOut)
      hC['VV'].append(hBkgr_VV.Clone(hBkgr_VV.GetName() + '_cloneToKillHighUnc'))
    if syst != 'Central':
      if 'JEC' in syst or 'JER' in syst:
        removeHighUncBins([hBkgr_VV,hC['VV'][0]])
    
    hBkgr_VV.Write()
      
    if syst == 'Central':
      if xDiv[0] == -1 and xDiv[1] == -1:
        hData_inc.Rebin(nRebin)
        h_inc['DY'][0].Rebin(nRebin)
        h_inc['DY'][1].Rebin(nRebin)
        h_inc['DY'][2].Rebin(nRebin)
        hBkgr_TT_inc.Rebin(nRebin)
        hBkgr_VV_inc.Rebin(nRebin)
      else:
        nameTmp = hData_inc.GetName()
        hData_inc = customBin(hData_inc, xDiv)
        hData_inc.SetName(nameTmp)
        nameTmp = h_inc['DY'][0].GetName()
        h_inc['DY'][0] = customBin(h_inc['DY'][0], xDiv)
        h_inc['DY'][0].SetName(nameTmp)
        nameTmp = h_inc['DY'][1].GetName()
        h_inc['DY'][1] = customBin(h_inc['DY'][1], xDiv)
        h_inc['DY'][1].SetName(nameTmp)
        nameTmp = h_inc['DY'][2].GetName()
        h_inc['DY'][2] = customBin(h_inc['DY'][2], xDiv)
        h_inc['DY'][2].SetName(nameTmp)
        nameTmp = hBkgr_TT_inc.GetName()
        hBkgr_TT_inc = customBin(hBkgr_TT_inc, xDiv)
        hBkgr_TT_inc.SetName(nameTmp)
        nameTmp = hBkgr_VV_inc.GetName()
        hBkgr_VV_inc = customBin(hBkgr_VV_inc, xDiv)
        hBkgr_VV_inc.SetName(nameTmp)
    
      hData_inc.Write()
      h_inc['DY'][0].Write()
      h_inc['DY'][1].Write()
      h_inc['DY'][2].Write()
      hBkgr_TT_inc.Write()
      hBkgr_VV_inc.Write()
  
  
  
  write_data_card(dc_name,fOut,systs_forDc_writting)

  fOut.Close()
