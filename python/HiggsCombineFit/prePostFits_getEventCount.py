from ROOT import *
import os,sys
sys.path.append('/uscms_data/d3/duong/CMSSW/CMSSW_7_6_5/src/ZplusC/python/')
from myutils.util_funcs import *

from math import *

gROOT.SetBatch(False)

def getFitResults(tex):
  lines = open(tex,'read').readlines()
  tot_err = [0,0,0] #c, b, l
  stat_err = [0,0,0]
  mcStat_err = [0,0,0]
  rate = [0,0,0]
  err = [0,0,0]
  for l in lines:
   lsp = l.split('&')
   if 'Total' in l:
     tot_err[0] = float(lsp[1])/100.
     tot_err[1] = float(lsp[2])/100.
     tot_err[2] = float(lsp[3])/100.
   if 'MC temp. stat.' in l:
     mcStat_err[0] = float(lsp[1])/100.
     mcStat_err[1] = float(lsp[2])/100.
     mcStat_err[2] = float(lsp[3])/100.
   if 'Stat.' in l:
     stat_err[0] = float(lsp[1])/100.
     stat_err[1] = float(lsp[2])/100.
     stat_err[2] = float(lsp[3])/100.
   if 'SF &' in l:
     rate[0] = float(lsp[2].split('$\pm$')[0]) # 2 = rc ==>TODO fix print results in doFit.py so Zc, Zb, Zl not Zb, Zc, Zl as current
     rate[1] = float(lsp[1].split('$\pm$')[0]) 
     rate[2] = float(lsp[3].split('$\pm$')[0]) 

  for i in range(3):
    err[i] = sqrt(tot_err[i]*tot_err[i] - mcStat_err[i]*mcStat_err[i] - stat_err[i]*stat_err[i])

  return rate, err

def getFitResults1(tex):
  lines = open(tex,'read').readlines()
  rate = [0,0,0]
  err = [0,0,0]
  for l in lines:
   lsp = l.split('&')
 
   if 'SF$_c$' in l and 'stat' in l:
     print l
     rate[0] = float(lsp[1].split('$\pm$')[0]) # 2 = rc ==>TODO fix print results in doFit.py so Zc, Zb, Zl not Zb, Zc, Zl as current
     err[0]  = float(lsp[1].split('$\pm$')[2].replace('(syst.)','').replace('\\',''))
   if 'SF$_b$' in l and 'stat' in l:
     rate[1] = float(lsp[1].split('$\pm$')[0])
     err[1]  = float(lsp[1].split('$\pm$')[2].replace('(syst.)','').replace('\\',''))
   if 'SF$_l$' in l and 'stat' in l:
     rate[2] = float(lsp[1].split('$\pm$')[0])
     err[2]  = float(lsp[1].split('$\pm$')[2].replace('(syst.)','').replace('\\',''))

  return rate, err



def changeBinErrors(his,his1,syst_frac_err): #his after scale by r, his1 before scale
  nBin = his.GetNbinsX()
  for i in range(1, nBin + 1):
    stat = his1.GetBinError(i)
    syst = his1.GetBinContent(i)*syst_frac_err
    tot = sqrt(stat*stat + syst*syst)
    his.SetBinError(i, tot)

def getNum(his):
  num = [0,0,0,0] #fit range, fit range error, all, all error
  nBin = his.GetNbinsX()
  for i in range(1,nBin+1):
    num[0] = num[0] + his.GetBinContent(i)
    num[1] = num[1] + his.GetBinError(i)*his.GetBinError(i)
  num[2] = num[0] + his.GetBinContent(0) + his.GetBinContent(nBin + 1)
  num[3] = num[1] + his.GetBinError(0)*his.GetBinError(0) + his.GetBinError(nBin + 1)*his.GetBinError(nBin + 1)
  num[1] = sqrt(num[1])
  num[3] = sqrt(num[3])
  return num

def normHist(h, normBinWidth = 0.2): #normalize the bin content of histogram to bin width of 0.2
  hOut = h.Clone(h.GetName() + '_normBinWidth_' + str(normBinWidth))
  for iB in range(1, hOut.GetNbinsX()+1):
    scale = normBinWidth/hOut.GetBinWidth(iB)
    hOut.SetBinContent(iB,hOut.GetBinContent(iB)*scale)
    hOut.SetBinError(iB,hOut.GetBinError(iB)*scale)
  return hOut

########################
#Main program
########################

doMakePlot = True 

fIn_name_pre = '../Test/forHiggsCombine_vtxMass_CSVM_pt20_DY_allData_allWeight_V25_binStatThresh_0_useDDbjet_useDDljet_customBinning'

dcName_pre = 'test_vtxMass_pt20_allData_binStatThresh_0_useDDbjet_useDDljet_customBinning'

fLatex = open('../tables_eventCount_prefit_vtxMass_pt20_binStatThresh_0_useDDbjet_useDDljet_customBinning.txt','w')

fOut = open('../summary_eventCout_vtxMass_pt20_binStatThresh_0_useDDbjet_useDDljet_customBinning.txt','w')

plotName = 'vtxMass_pt20_binStatThresh_0_useDDbjet_useDDljet_customBinning'

chans = ['Zee','Zmm']
#chans = ['Zee']

nums = {'Zee':{}, 'Zmm':{}}
nums_inc = {'Zee':{}, 'Zmm':{}}

xAxisRange = [0,10] 
if '_vtxMass_' in dcName_pre:
  xAxisRange = [0,6]

for chan in chans:
  dcName = dcName_pre + '_' + chan
  tex = dcName + '/' + dcName + '_newErrorEstimation_latex.txt'
  print '>>>>>', tex

  fIn = TFile.Open(fIn_name_pre + '_' + chan + '.root','read')
  print fIn.GetName()
  h1 = []
  hName = []
  hName.append('Data')
  h1.append(fIn.Get('data_obs').Clone('Data1'))
  hName.append('VV')
  h1.append(fIn.Get('VV').Clone('VV1'))
  hName.append('TT')
  h1.append(fIn.Get('TT').Clone('TT1'))
  hName.append('Zl')
  h1.append(fIn.Get('Zl').Clone('Zl1'))
  hName.append('Zb')
  h1.append(fIn.Get('Zb').Clone('Zb1'))
  hName.append('Zc')
  h1.append(fIn.Get('Zc').Clone('Zc1'))
  h1_forPlot = []
  for i in range(len(h1)):
    h1_forPlot.append(normHist(h1[i]))

  cName = 'prefit_' + plotName + '_' + chan
  if doMakePlot: makeStackPlot(h1_forPlot,hName,cName,'Plots/','Jet M_{SV} [GeV]', xAxisRange, 'MC stat. unc.', False, False, 0.2)
  
  hMC = h1[-1].Clone('MC')
  for i in range(1, len(h1) - 1):
    hMC.Add(h1[i])

  for i in range(len(h1)):
    nums[chan][hName[i]] = getNum(h1[i])
  
  nums[chan]['MC'] = getNum(hMC)

  #get inclusive event count
  hTmp = (fIn.Get('data_obs_inc').Clone('Data1_inc'))
  hTmp.Add(fIn.Get('TT_inc'), -1)
  hTmp.Add(fIn.Get('VV_inc'), -1)
  nums_inc[chan]['Data_inc'] = getNum(hTmp)

  ##########################
  #after fit plot
  ##########################

  rate, r_err = getFitResults1(tex)
  print '>>>>>>>> rate, r_err: ', rate, ' ', r_err

  h2 = []
  h2.append(fIn.Get('data_obs').Clone('Data2'))
  h2.append(fIn.Get('VV').Clone('VV2'))
  h2.append(fIn.Get('TT').Clone('TT2'))
  h2.append(fIn.Get('Zl').Clone('Zl2'))
  h2.append(fIn.Get('Zb').Clone('Zb2'))
  h2.append(fIn.Get('Zc').Clone('Zc2'))

  for i in range(1,len(h2)):
    if 'Zc' in hName[i]:
      scale = rate[0]
      h2[i].Scale(scale)
      changeBinErrors(h2[i], h1[i], r_err[0])
    if 'Zb' in hName[i]:
      scale = rate[1]
      h2[i].Scale(scale)
      changeBinErrors(h2[i], h1[i], r_err[1])
    if 'Zl' in hName[i]:
      scale = rate[2]
      h2[i].Scale(scale)
      changeBinErrors(h2[i], h1[i], r_err[2])
  
  h2_forPlot = []
  for i in range(len(h2)):
    h2_forPlot.append(normHist(h2[i]))
  cName = 'postfit_' + plotName + '_' + chan
  if doMakePlot: makeStackPlot(h2_forPlot,hName,cName,'Plots/','Jet M_{SV} [GeV]', xAxisRange, 'Stat. #oplus syst.', False, False, 0.2)

  ########################
  #get inclusive events
  ########################

print nums
print nums_inc
##########################
#Prefit table
##########################
#fLatex = open('xxx.txt','w')
fLatex.write('\n')
fLatex.write('\\begin{table}\n')
fLatex.write('  \centering\n')
fLatex.write('  \\begin{tabular}{lcc}\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
l = ['', 'Electron ', 'Muon ']
st = makeLatexLine(l)
fLatex.write(st)
fLatex.write('  \\hline\n')
labels = ['Zc','Zb','Zl','TT','VV','MC','Data']
for k in labels:
  if k != 'Data':
    l = [k, format('%.0f' % nums['Zee'][k][0]) + ' $\pm$ '  + format('%.0f' % nums['Zee'][k][1]) + ' (' + format('%.0f' % nums['Zee'][k][2]) + ' $\pm$ '  + format('%.0f' % nums['Zee'][k][3]) + ')', format('%.0f' % nums['Zmm'][k][0]) + ' $\pm$ '  + format('%.0f' % nums['Zmm'][k][1]) + ' (' + format('%.0f' % nums['Zmm'][k][2]) + ' $\pm$ '  + format('%.0f' % nums['Zmm'][k][3]) + ')']
  else: 
    l = [k, format('%.0f' % nums['Zee'][k][0]) + ' (' + format('%.0f' % nums['Zee'][k][2]) + ')', format('%.0f' % nums['Zmm'][k][0]) + ' (' + format('%.0f' % nums['Zmm'][k][2]) + ')']
  if k == 'MC':
    fLatex.write('  \\hline\n')
  st = makeLatexLine(l)
  fLatex.write(st)

fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\end{tabular}\n')
fLatex.write('\\end{table}\n')

#######################
#Event cout
#######################
for ch in ['Zee','Zmm']:  
  for flav in ['Zc','Zb']:
    fOut.write(ch + ' ' + flav + ' ' + str(nums[ch][flav][2]) + ' ' + str(nums[ch][flav][3]) + '\n') # 2 = use overflow events

for ch in ['Zee','Zmm']:  
  fOut.write(ch + ' Data_inc ' + str(nums_inc[ch]['Data_inc'][2]) + ' ' + str(nums_inc[ch]['Data_inc'][3]) + '\n') # 2 = use overflow events
