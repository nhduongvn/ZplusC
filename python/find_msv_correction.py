from ROOT import *

vtxMassType = 'vtxMass' #vtxMass, incVtxMass, vtxMassCorr_IVF
tagName = 'CSVM' #CtagT, CSVM

fIn = []
fIn.append(TFile.Open('Test/template_' + vtxMassType + '_' + tagName + '_emu_allData_V25.root','read'))
fIn.append(TFile.Open('Test/template_' + vtxMassType + '_' + tagName + '_ttsemi_allData_V25.root','read'))
fDY = TFile.Open('Test/template_' + vtxMassType + '_' + tagName + '_pt20_DY_allData_allWeight_V25.root', 'read')
#fIn.append(TFile.Open('Test/template_' + vtxMassType + '_' + tagName + '_emu_DY_nlo_allData_V25.root','read'))
#fIn.append(TFile.Open('Test/template_' + vtxMassType + '_' + tagName + '_ttsemi_DY_nlo_allData_V25.root','read'))
#fDY = TFile.Open('Test/template_' + vtxMassType + '_' + tagName + '_pt20_noGapEle_DY_nlo_allData_allWeight_V25.root', 'read')


systs = ['Central', 'JECDown', 'JECUp', 'gbbDown', 'gbbUp']

#hNames = {
#           'TT_bjet':["bjet_vtxMass_TT_emu_px","bjet_vtxMass_TT_JECDown_emu_px", "bjet_vtxMass_TT_JECUp_emu_px", 'bjet_vtxMass_TT_gbbDown_emu_px','bjet_vtxMass_TT_gbbUp_emu_px'],
#           'DY_bjet':['bjet_vtxMass_DY','bjet_vtxMass_DY_JECDown','bjet_vtxMass_DY_JECUp', '']
#         }

#fitFuncs = ["[0] + [1]*x + [2]*x*x + [3]*x^[4]","[0]+[1]*x"]
#fitFuncs = ["[0] + [1]*x + [2]*x*x + [3]*x*x*x","[0] + [1]*x"]
fitFuncs = { 'TT_bjet': ["[0] + [1]*x + [2]*x*x + [3]/x", "[0] + [1]*x + [2]*x*x + [3]/x", "[0] + [1]*x + [2]*x*x + [3]/x", "[0] + [1]*x + [2]*x*x + [3]/x", '[0] + [1]*x + [2]*x*x + [3]/x'],
           }

#################################
#get histograms
#################################
hist = {'DY_bjet':[], 'TT_bjet':[]} 

for syst in systs:
  name = 'bjet_vtxMass_DY'
  if syst != 'Central':
    name = name + '_' + syst
  hist['DY_bjet'].append(fDY.Get(name + '_Zee').Clone(name + '_All_for_'+syst))
  hTmp = fDY.Get(name + '_Zmm').Clone(name + '_Zmm_for_' + syst)
  hist['DY_bjet'][-1].Add(hTmp)
  
  name = 'bjet_vtxMass_TT_emu_px'
  if syst != 'Central' and 'gcc' not in syst and 'gbb' not in syst:
    name = 'bjet_vtxMass_TT_' + syst + '_emu_px'
  print name
  hist['TT_bjet'].append(fIn[0].Get(name).Clone(name + '_for_' + syst))


print hist

cans = []
print len(hist['DY_bjet'])
corr_funcs = {}

for iHist in range(len(hist['DY_bjet'])):
  hCorr = hist['DY_bjet'][iHist].Clone(hist['DY_bjet'][iHist].GetName() + '_bjet_corr')
  hTmp = hist['TT_bjet'][iHist].Clone(hist['TT_bjet'][iHist].GetName() + '_clone')
  nRebin = 2
  hCorr.Rebin(nRebin)
  hTmp.Rebin(nRebin)
  hCorr.Scale(1./hCorr.Integral())
  hTmp.Scale(1./hTmp.Integral())
  hCorr.Divide(hTmp)
  fitRange = [0.2,7]
  if tagName == 'CtagT': fitRange = [0.1,7]
  f1 = TF1("fit_bjet_" + str(iHist),fitFuncs['TT_bjet'][iHist],fitRange[0],fitRange[1])
  hCorr.Fit('fit_bjet_' + str(iHist),'R 0')
  nameTmp = 'msv_corr_bjet_' + systs[iHist] + '_' + tagName
  cans.append(TCanvas(nameTmp, nameTmp))
  hCorr.SetTitle("")
  hCorr.Draw()
  hCorr.SetStats(0)
  hCorr.GetXaxis().SetRangeUser(0,6)
  hCorr.GetYaxis().SetRangeUser(0.3,1.5)
  hCorr.SetMarkerStyle(20)
  hCorr.GetXaxis().SetTitle('M_{SV} [GeV]')
  hCorr.GetYaxis().SetTitle('Correction')
  f1.Draw("same")
  fitFuncs_pr = fitFuncs['TT_bjet'][iHist]
  fitFuncs_forText = fitFuncs['TT_bjet'][iHist]
  for iPar in range(f1.GetNpar()):
    if f1.GetParameter(iPar) >= 0:
      fitFuncs_pr = fitFuncs_pr.replace('['+str(iPar)+']',str(f1.GetParameter(iPar)))
      fitFuncs_forText = fitFuncs_forText.replace('['+str(iPar)+']',format(f1.GetParameter(iPar), '0.4f'))
    else:
      fitFuncs_pr = fitFuncs_pr.replace('+ ['+str(iPar)+']',str(f1.GetParameter(iPar)))
      fitFuncs_forText = fitFuncs_forText.replace('+ ['+str(iPar)+']',format(f1.GetParameter(iPar), '0.4f'))

  print fitFuncs_pr
  corr_funcs['TT_bjet_' + systs[iHist]] = fitFuncs_pr
  
  text = TLatex()
  text.DrawLatexNDC(0.15, 0.85, 'Correction function: ')
  text.DrawLatexNDC(0.15, 0.75, ' ' + fitFuncs_forText.replace('*x*x', '*x^{2}'))
  cans[-1].Update()
  cans[-1].Print('Test/' + cans[-1].GetName() + '.pdf')
  cans[-1].Print('Test/' + cans[-1].GetName() + '.C')
  cans[-1].Print('Test/' + cans[-1].GetName() + '.eps')
  cans[-1].Print('Test/' + cans[-1].GetName() + '.png')

print 'Number of canvas: ', len(cans)
print corr_funcs 

fTex = open('Test/msv_corr_bjet_' + tagName + '.txt','w')
c1 = TCanvas('corr_funcs_bjet_' + tagName)
leg = TLegend(0.6, 0.7, 0.9, 0.9)
f1 = []
i = 0
for k,v in corr_funcs.items():
  fTex.write('\n ' + k + ' ' + v)
  f1.append(TF1("fit_func_" + str(i),v,0.2,7))
  if i == 0:
    f1[-1].Draw()
    f1[-1].SetTitle("")
    f1[-1].GetYaxis().SetRangeUser(0,1.3)
    f1[-1].GetYaxis().SetTitle("Correction")
    f1[-1].GetXaxis().SetTitle("M_{SV} [GeV]")
    leg.AddEntry(f1[-1], k, 'l')
  else:
    f1[-1].Draw("same")
    f1[-1].SetLineColor(5+i)
    leg.AddEntry(f1[-1], k, 'l')

  if 'Central' not in k: f1[-1].SetLineStyle(2)
  
  i += 1

leg.Draw()

c1.Print('Test/' + c1.GetName() + '.pdf')

