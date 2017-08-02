import os,sys
import copy
from math import *
from ROOT import *
sys.path.append('/uscms_data/d3/duong/CMSSW/CMSSW_7_6_5/src/ZplusC/python/')
from myutils.util_funcs import *


def checkFit(log1, log2): #log1 multiDim, log2 mll
  lines = open(log2).readlines()
  print '>>>>>>>>Fit status: '
  for l in lines:
    if 'Hesse is valid - matrix is accurate' in l: print l #not working currently since log do not dump to log file
    if 'covariance matrix quality' in l: print l
    if 'covariance matrix status:' in l: print l #using Minuit convention : =0 not calculated, =1 approximated, =2 made pos def , =3 accurate (https://root.cern.ch/root/html/ROOT__Fit__FitResult.html)
    if 'Status : MINIMIZE=' in l: print l

  
def getMatrix(lines, pos, nNui): #pos: position of first row
  m = []
  for j in range(3): #loop over row
    l = lines[pos + j]
    row = []
    for i in range(3): #loop over column
      row.append(float(l.split()[i+nNui]))
    m.append(row)

  return m

def getResults(log1, log2, useMultiDimResult=True): #log1 multiDim, log2 mll
  debug = True
  rat = {'r1': [], 'r2': [], 'r3': []} #Zc, Zb, Zl: store fit value, error and global correlation getting from mll fit
  rat1 = {'r1': [], 'r2': [], 'r3': []} #Zc, Zb, Zl: store fit value, error and global correlation getting from multidim fit (for cross check)
  cov = [] #covariance matrix, from mll fit
  cor = [] #correlation matrix, from mll fit
  
  lines = open(log2).readlines()
  pos1 = -1 #line position to read r
  pos2 = -1 #line position to read cov, cor, glob cor
  nl = -1
  for l in lines:
    nl = nl + 1
    if 'Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.' in l: pos1 = nl
    if 'MnUserCovariance:' in l: pos2 = nl
  
  nNui = -1
  for i in range(pos1,len(lines)):
    #print lines[i]
    nNui = nNui + 1
    if 'r1' in lines[i]: break
  nNui = nNui - 2
     
  if debug:
    print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    print pos1, ' ', nNui
    print lines[pos1 + nNui + 2]
    print lines[pos1 + nNui + 3]
    print lines[pos1 + nNui + 4]
    print lines[pos2 + 3*nNui + 15]
    print lines[pos2 + 3*nNui + 16]
    print lines[pos2 + 3*nNui + 17]

  r1 = lines[pos1 + nNui+ 2].split()[2]
  r2 = lines[pos1 + nNui+ 3].split()[2]
  r3 = lines[pos1 + nNui+ 4].split()[2]
  r1_err = lines[pos1 + nNui+ 2].split()[4]
  r2_err = lines[pos1 + nNui+ 3].split()[4]
  r3_err = lines[pos1 + nNui+ 4].split()[4]

  r1_globCorr = lines[pos2 + 3*nNui+ 15].split()[0]
  r2_globCorr = lines[pos2 + 3*nNui+ 16].split()[0]
  r3_globCorr = lines[pos2 + 3*nNui+ 17].split()[0]

  if debug: print '>>> Found r1, r2, r3: ', r1, ' ', r2, ' ', r3, ' ', r1_err, ' ', r2_err, ' ', r3_err

  rat['r1'].append(float(r1))
  rat['r1'].append(float(r1_err))
  rat['r1'].append(float(r1_globCorr))
  rat['r2'].append(float(r2))
  rat['r2'].append(float(r2_err))
  rat['r2'].append(float(r2_globCorr))
  rat['r3'].append(float(r3))
  rat['r3'].append(float(r3_err))
  rat['r3'].append(float(r3_globCorr))

  cov = getMatrix(lines, pos2 + nNui + 2, nNui)
  cor = getMatrix(lines, pos2 + 2*nNui + 8, nNui)
  
  if debug:
    print rat
    print cov
    print cor
    print sqrt(cov[0][0]), ' ', sqrt(cov[1][1]), ' ', sqrt(cov[2][2]), ' ', r1_err, ' ', r2_err, ' ', r3_err
  
  #now read multiDimFit log to get r1, r2, r3 from multiDimFit
  lines = open(log1).readlines()
  for l in lines:
    if '(limited)' not in l: continue
    if 'r1' in l:
      rat1['r1'] = [float(l.split()[2]), float(l.split()[4])]
    if 'r2' in l:
      rat1['r2'] = [float(l.split()[2]), float(l.split()[4])]
    if 'r3' in l:
      rat1['r3'] = [float(l.split()[2]), float(l.split()[4])]

  if useMultiDimResult:
    for k,v in rat.items():
      rat[k][0] = rat1[k][0] 
      rat[k][1] = rat1[k][1] 

  if debug:
    print r1

  for k,v in rat.items():
    diff = abs(v[0] - rat1[k][0])/v[0]
    if diff > 1e-3:
      print '>>>Warning ', k, ' results are not the same ', v[0], ' ', rat1[k][0]
    diff = abs(v[1] - rat1[k][1])/v[1]
    if diff > 1e-3:
      print '>>>Warning ', k, ' error results are not the same ', v[1], ' ', rat1[k][1]
  
  ##########
  #return the result
  ##########
  return rat, cov, cor

def getPreFit(dcName):
  rate_prefit = {'Zb':[-1,-1], 'Zc':[-1,-1], 'Zl':[-1,-1]}
  lines = open(dcName + '.txt').readlines()
  nl = -1
  row = -1
  for l in lines:
    nl = nl + 1
    if 'Zb' in l and 'Zc' in l and 'Zl' in l and 'stat_bin' not in l:
      row = nl
  tmp = lines[row].split()
  columns = {'Zb':-1,'Zc':-1,'Zl':-1}
  columns['Zb'] = tmp.index('Zb')
  columns['Zc'] = tmp.index('Zc')
  columns['Zl'] = tmp.index('Zl')
  rate_prefit['Zb'][0] = float(lines[row+2].split()[columns['Zb']])
  rate_prefit['Zc'][0] = float(lines[row+2].split()[columns['Zc']])
  rate_prefit['Zl'][0] = float(lines[row+2].split()[columns['Zl']])
  for k,v in rate_prefit.items():
    rate_prefit[k][1] = sqrt(rate_prefit[k][0])

  return rate_prefit

def cal_error_f(r1, r2, r3, n1, n2, n3, V): #calculate the error of flavour fraction  r1*n1/(r1*n1+r2*r2+r3*r3)
  #cross check if V_00, V_11, V_22 is r1_err, r2_err, r3_err
  #print '>>> V_00, V_11, V_22: ', sqrt(V[0][0]), ' ', sqrt(V[1][1]), ' ', sqrt(V[2][2])
  N = r1*n1 + r2*n2 + r3*n3
  N1 = r1*n1
  N23 = r2*n2 + r3*n3
  dev_r1 = n1*N23/(N*N)
  dev_r2 = -N1*n2/(N*N)
  dev_r3 = -N1*n3/(N*N)
  sig = dev_r1*dev_r1*V[0][0] + dev_r2*dev_r2*V[1][1] + dev_r3*dev_r3*V[2][2] + 2*dev_r1*dev_r2*V[0][1] + 2*dev_r1*dev_r3*V[0][2] + 2*dev_r2*dev_r3*V[1][2]
  return sqrt(sig)
 
def cal_error_rcrb(r1, err_r1, r2, err_r2, V, mode = 1): #mode = 0,1,2 no correlation, useCovariant, useCorrelation y = r1/r2 r1 = rc, r2 = rb
  #print '>>> V_00, V_11, V_22: ', sqrt(V[0][0]), ' ', sqrt(V[1][1]), ' ', sqrt(V[2][2])
  #print '>>> V_01: ', V[0][1] 
  dev_r1 = 1/r2
  dev_r2 = -r1/(r2*r2)
  if mode == 0:
    sig = dev_r1*dev_r1*err_r1*err_r1 + dev_r2*dev_r2*err_r2*err_r2
  if mode == 1:
    sig = dev_r1*dev_r1*err_r1*err_r1 + dev_r2*dev_r2*err_r2*err_r2 + 2*dev_r1*dev_r2*V[0][1]
  if mode == 2:
    if abs(V[0][1]) > 1:
      print 'Warning: correlation is greater than 1, will use non-correlation'
      sig = dev_r1*dev_r1*err_r1*err_r1 + dev_r2*dev_r2*err_r2*err_r2
    else:
      print '>>>>>Correlation: ', V[0][1] 
      sig = dev_r1*dev_r1*err_r1*err_r1 + dev_r2*dev_r2*err_r2*err_r2 + 2*dev_r1*dev_r2*V[0][1]*err_r1*err_r2

  if sig < 0:
    print '>>>>>>>Warning: Negative found when using covariant for rc/rb error, use non-covariant or non-correlation', sig, ' ', V[0][1]
    sig = dev_r1*dev_r1*err_r1*err_r1 + dev_r2*dev_r2*err_r2*err_r2

  return sqrt(sig)

def cal_error_rcrb_mllFit(rc,rb,unc_c,unc_b,cov_c,cov_b, mode=1):
  zc_unc,zb_unc = 0,0
  if len(unc_c) == 1:
    zc_unc = abs(unc_c[0])
  else:
    zc_unc = max([abs(unc_c[0]),abs(unc_c[1])])
  if len(unc_b) == 1:
    zb_unc = abs(unc_b[0])
  else:
    zb_unc = max([abs(unc_b[0]),abs(unc_b[1])])
  uncTmp = cal_error_rcrb(rc,zc_unc,rb,zb_unc,cov_c,mode) 
  uncTmp1 = cal_error_rcrb(rc,zc_unc,rb,zb_unc,cov_b,mode)
  return max([uncTmp,uncTmp1])
  

def printResults(rate_prefit, rateIn, V, label, fIn): #rate_prefit = number of event, rate = correction, r, found after fitting --> number of events after fit = rate*rate_prefit
  
  #rate = rateIn[:] --> for array
  rate = copy.deepcopy(rateIn)

  #change the key to change label
  for k,v in label.items():
    rate[v] = rate.pop(k)

  frac_prefit = {'Zb': [-1, -1], 'Zc': [-1, -1], 'Zl': [-1, -1]}
  
  nTot_prefit = 0
  for k,v in rate_prefit.items(): nTot_prefit = nTot_prefit + v[0]
  for k,v in frac_prefit.items():
    frac_prefit[k][0] = rate_prefit[k][0]/nTot_prefit
  Vtmp = [[0,0,0],[0,0,0],[0,0,0]] #Zb, Zc, Zl
  Vtmp[0][0] = rate_prefit['Zb'][1]*rate_prefit['Zb'][1]
  Vtmp[1][1] = rate_prefit['Zc'][1]*rate_prefit['Zc'][1]
  Vtmp[2][2] = rate_prefit['Zl'][1]*rate_prefit['Zl'][1]
  frac_prefit['Zb'][1] = cal_error_f(rate_prefit['Zb'][0],rate_prefit['Zc'][0],rate_prefit['Zl'][0],1,1,1,Vtmp)
  frac_prefit['Zc'][1] = cal_error_f(rate_prefit['Zc'][0],rate_prefit['Zb'][0],rate_prefit['Zl'][0],1,1,1,Vtmp)
  frac_prefit['Zl'][1] = cal_error_f(rate_prefit['Zl'][0],rate_prefit['Zb'][0],rate_prefit['Zc'][0],1,1,1,Vtmp)
  
  rate_err = {'Zb': -1, 'Zc': -1, 'Zl': -1} #in percentage
  for k,v in rate.items(): rate_err[k] = rate[k][1]/rate[k][0]
  
  rate_afterFit = {'Zb':[-1,-1],'Zc':[-1,-1],'Zl':[-1,-1]}
  for k,v in rate_afterFit.items():
    rate_afterFit[k][0] = rate_prefit[k][0]*rate[k][0]
    rate_afterFit[k][1] = rate_prefit[k][0]*rate[k][0]*rate_err[k] #error

  nTot_afterFit = 0
  for k,v in rate_afterFit.items(): nTot_afterFit= nTot_afterFit + v[0]
  frac_afterFit = {'Zb':[-1,-1], 'Zc':[-1,-1], 'Zl':[-1,-1]}
  for k,v in frac_afterFit.items():
    frac_afterFit[k][0] = rate_afterFit[k][0]/nTot_afterFit
  frac_afterFit['Zb'][1] = cal_error_f(rate['Zb'][0],rate['Zc'][0],rate['Zl'][0],rate_prefit['Zb'][0], rate_prefit['Zc'][0], rate_prefit['Zl'][0], V)
  frac_afterFit['Zc'][1] = cal_error_f(rate['Zc'][0],rate['Zb'][0],rate['Zl'][0],rate_prefit['Zc'][0], rate_prefit['Zb'][0], rate_prefit['Zl'][0], V)
  frac_afterFit['Zl'][1] = cal_error_f(rate['Zl'][0],rate['Zb'][0],rate['Zc'][0],rate_prefit['Zl'][0], rate_prefit['Zb'][0], rate_prefit['Zc'][0], V)

  print ">>>>>>>>>>>>>>>", rate
  
  fIn.write('\\begin{table}\n')
  fIn.write('\\resizebox{\columnwidth}{!}{\n')
  fIn.write('\\begin{tabular}{c|c|c|c|c|c|c|c}\n')
  fIn.write('\\hline\n')
  fIn.write('\\hline\n')
  st = '{0:20}'.format('&') + '{0:20}'.format('Zb &') + '{0:20}'.format('Zc &') + '{0:20}'.format('Zl &') + '{0:20}'.format('Zc/Zb &') + '{0:20}'.format('fb &') + '{0:20}'.format('fc &') + '{0:20}'.format('fl\\\\\n')
  fIn.write(st)
  fIn.write('\\hline\n')
  st = '{0:20}'.format('Prefit &') + '{0:20}'.format('%.0f $\pm$ %.0f &' % (rate_prefit['Zb'][0], rate_prefit['Zb'][1])) + '{0:20}'.format('%.0f $\pm$ %.0f &' % (rate_prefit['Zc'][0], rate_prefit['Zc'][1])) + '{0:20}'.format('%.0f $\pm$ %.0f &' % (rate_prefit['Zl'][0], rate_prefit['Zl'][1])) + '{0:20}'.format('- & ') + '{0:20}'.format('%.3f $\pm$ %.3f &' % (frac_prefit['Zb'][0],frac_prefit['Zb'][1])) + '{0:20}'.format('%.3f $\pm$ %.3f &' % (frac_prefit['Zc'][0],frac_prefit['Zc'][1])) + '{0:20}'.format('%.3f $\pm$ %.3f \\\\\n' % (frac_prefit['Zl'][0], frac_prefit['Zl'][1]))
  fIn.write(st)
  st = '{0:20}'.format('SF &') + '{0:20}'.format('%.3f $\pm$ %.3f &' % (rate['Zb'][0], rate['Zb'][1])) + '{0:20}'.format('%.3f $\pm$ %.3f &' % (rate['Zc'][0], rate['Zc'][1])) + '{0:20}'.format('%.3f $\pm$ %.3f &' % (rate['Zl'][0], rate['Zl'][1])) + '{0:20}'.format('%.3f $\pm$ %.3f &' % (rate['Zc_over_Zb'][0], rate['Zc_over_Zb'][1])) + '{0:20}'.format('- & ') + '{0:20}'.format('- & ') + '{0:20}'.format('- \\\\\n')
  fIn.write(st)
  st = '{0:20}'.format('After fit &') + '{0:20}'.format('%.0f $\pm$ %.0f &' % (rate_afterFit['Zb'][0],rate_afterFit['Zb'][1])) + '{0:20}'.format('%.0f $\pm$ %.0f &' % (rate_afterFit['Zc'][0], rate_afterFit['Zc'][1])) + '{0:20}'.format('%.0f $\pm$ %.0f &' % (rate_afterFit['Zl'][0], rate_afterFit['Zl'][1])) + '{0:20}'.format('%.0f $\pm$ %.0f &' % (rate_afterFit['Zl'][0], rate_afterFit['Zl'][1])) + '{0:20}'.format('%.3f $\pm$ %.3f &' % (frac_afterFit['Zb'][0],frac_afterFit['Zb'][1])) + '{0:20}'.format('%.3f $\pm$ %.3f &' % (frac_afterFit['Zc'][0], frac_afterFit['Zc'][1])) + '{0:20}'.format('%.3f $\pm$ %.3f \\\\\n' % (frac_afterFit['Zl'][0], frac_afterFit['Zl'][1]))
  fIn.write(st)
  fIn.write('\\hline\n')
  fIn.write('\\hline\n')
  fIn.write('\\end{tabular}\n')
  fIn.write('}\n')
  fIn.write('\\end{table}\n')
  fIn.write('\n')


def plotGoodnessOfFit(fInName1, fInName2,dcName):
  doBatch = True 
  f = TFile('higgsCombine'+fInName1.replace('.root','.GoodnessOfFit.mH120.root'),'read')
  t = f.Get('limit')
  c = TCanvas('c_GoodnessOfFit_' + dcName)
  vals = []
  for e in range(t.GetEntries()):
        t.GetEntry(e)
        vals.append(t.limit)

  t.Draw('limit>>hist(200,0,200)')
  hist = gDirectory.Get('hist')
  hist.SetTitle('Test-statistic distribution')
  hist.Draw("HIST")

  fobs = TFile('higgsCombine'+fInName2.replace('.root','.GoodnessOfFit.mH120.123456.root'),'read')
  tobs = fobs.Get('limit')
  tobs.GetEntry(0)
  obs = tobs.limit

  pval = sum(1.0 for i in vals if i >= obs) / float(len(vals))

  arr = TArrow(obs, 0.001, obs, hist.GetMaximum()/4, 0.02, "<|")
  arr.SetLineColor(kBlue)
  arr.SetFillColor(kBlue)
  arr.SetFillStyle(1001)
  arr.SetLineWidth(6)
  arr.SetLineStyle(1)
  arr.SetAngle(60)
  arr.Draw("<|same")

  c.Draw()

  if not doBatch:
    pause()

  c.Print('goodnessOfFit_' + dcName + '.C')
  c.Print('goodnessOfFit_' + dcName + '.pdf')
  c.Print('goodnessOfFit_' + dcName + '.eps')
  c.Print('goodnessOfFit_' + dcName + '.png')

  print '======================'
  print dcName_pre
  print 'p-value for obs=%f is %f' % (obs, pval)

def getUnc(log):
  pipe = os.popen('grep "Best fit r1:" ' + log).readlines()
  up = pipe[0].split()[4].split('/')[1]
  down = pipe[0].split()[4].split('/')[0]
  
  pipe = os.popen('grep "r1" ' + log + ' | grep "<none>"').readlines()
  #if 'Zc' in log: print pipe
  if '-0.00' in pipe[-1] and '+0.00' in pipe[0]:
    print '>>>>Warning: no unc. found ', up, ' ', down
  elif '+0.00' in pipe[-1]:
    print '>>>>Warning: no UP unc. found ', up, ' ', down, ' . Setting up to down'
    up = down
  elif '-0.00' in pipe[-1]:
    print '>>>>Warning: no DOWN unc. found ', up, ' ', down, ' . Setting down to up'
    down = up
  
  return [abs(float(up)), -1.0*abs(float(down))]
  
def getCentral(log):
  pipe = os.popen('grep "Best fit r1:" ' + log).readlines()
  return float(pipe[0].split()[3])

def getCov(log): 
  debug = False
  cov = [] #covariance matrix, from mll fit
  cor = [] #correlation matrix, from mll fit
  
  lines = open(log).readlines()
  pos1 = -1 #line position to read r
  pos2 = -1 #line position to read cov, cor, glob cor
  nl = -1
  for l in lines:
    nl = nl + 1
    if 'Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.' in l: pos1 = nl
    if 'MnUserCovariance:' in l: pos2 = nl
  
  nNui = -1
  for i in range(pos1,len(lines)):
    #print lines[i]
    nNui = nNui + 1
    if 'r1' in lines[i]: break
  nNui = nNui - 2
     
  if debug:
    print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    print pos1, ' ', nNui
    print lines[pos1 + nNui + 2]
    print lines[pos1 + nNui + 3]
    print lines[pos1 + nNui + 4]
    print lines[pos2 + 3*nNui + 15]
    print lines[pos2 + 3*nNui + 16]
    print lines[pos2 + 3*nNui + 17]


  cov = getMatrix(lines, pos2 + nNui + 2, nNui)
  cor = getMatrix(lines, pos2 + 2*nNui + 8, nNui)

  return [cov, cor]

def getUnc_r(r,tot,stat):
  return [stat,sqrt(tot*tot-stat*stat)]


#############################
# Main program
#############################

dcName_pre = 'test_vtxMass_pt20_allData_binStatThresh_0_useDDbjet_useDDljet_customBinning'
chans = ['Zee','Zmm']
#chans = ['Zmm']
useRelUnc = True #use relative uncertainty in unc. breakdown


fOut = open('../summary_analyzeFit_vtxMass_pt20_binStatThresh_0_useDDbjet_useDDljet_customBinning.txt', 'w')
#fOut = open('testxxx.txt', 'w')

for chan in chans:
  print '#########################################################'
  print '#', chan
  print '#########################################################'
  
  dcName = dcName_pre + '_' + chan

  wDir = dcName
  os.chdir(wDir)

  fLatex = open(dcName + '_newErrorEstimation_latex.txt','w')
  
  #Get total unc.
  totUncs = {'Zc':[],'Zb':[],'Zl':[],'Zc_Zb':[]}
  rs = {'Zc':-1,'Zb':-1,'Zl':-1,'Zc_Zb':-1} #from best fit

  #Get unc from mll fit
  fitUncs = {'Zc':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
            'Zb':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
            'Zl':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
            'Zc_Zb':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]}}
             #syst when freezing that MC group [0,1] and not freezeing that MC group [0,1], Exp == freezing all nuisance syst --> this is stat uncertainty [0,1]
  rs_fits = {'Zc':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
            'Zb':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
            'Zl':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
            'Zc_Zb':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]}} #store the central value when doing syst breakdown 
 
  cov = {'Zc':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
      'Zb':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
      'Zl':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]}}

  cor = copy.deepcopy(cov)
  
  cor_central = {'Zc':[],'Zb':[],'Zl':[]}

  for flav,systs in fitUncs.items():
    if flav == 'Zc_Zb': continue
    log = 'Logs/log_' + flav + '_mllFit.txt'
    print '#################'
    print log
    totUncs[flav].append(getUnc(log)[0])
    totUncs[flav].append(getUnc(log)[1])
    rs[flav] = getCentral(log)
    cor_central[flav] = getCov(log)[1]

    for s,v in systs.items():
      log = []
      log.append('Logs/log_' + flav + '_syst_mllFit_freeze_' + s + '.txt') 
      log.append('Logs/log_' + flav + '_syst_mllFit_not_freeze_' + s + '.txt') 
      print '################'
      print log
      fitUncs[flav][s].append(getUnc(log[0])[0])
      fitUncs[flav][s].append(getUnc(log[0])[1])
      cov[flav][s].append(getCov(log[0])[0])
      cor[flav][s].append(getCov(log[0])[1])
      rs_fits[flav][s].append(getCentral(log[0]))

      if s != 'Exp':
        fitUncs[flav][s].append(getUnc(log[1])[0])
        fitUncs[flav][s].append(getUnc(log[1])[1])
        cov[flav][s].append(getCov(log[1])[0])
        cor[flav][s].append(getCov(log[1])[1])
        rs_fits[flav][s].append(getCentral(log[1]))
  
  #calculate ratio Zc_Zb and its unc. from best fit
  rs['Zc_Zb'] = rs['Zc']/rs['Zb']
  log_zc = 'Logs/log_Zc_mllFit.txt'
  log_zb = 'Logs/log_Zb_mllFit.txt'
  cov_zc = getCov(log_zc)[0] # 1 is for correlation
  cov_zb = getCov(log_zb)[0]
  zc_zb_unc = cal_error_rcrb_mllFit(rs['Zc'],rs['Zb'],totUncs['Zc'], totUncs['Zb'],cov_zc,cov_zb)
  totUncs['Zc_Zb'] = [zc_zb_unc, -1.0*zc_zb_unc]

  #calculate the fit uncertainty for ratio
  for s in fitUncs['Zc_Zb'].keys():
    mode = 1
    if s == 'Exp': mode = 0
    zc_unc = [fitUncs['Zc'][s][0],fitUncs['Zc'][s][1]]
    zb_unc = [fitUncs['Zb'][s][0],fitUncs['Zb'][s][1]]
    zc_zb_unc = cal_error_rcrb_mllFit(rs['Zc'],rs['Zb'],zc_unc,zb_unc,cov['Zc'][s][0],cov['Zb'][s][0],mode)
    fitUncs['Zc_Zb'][s].append(zc_zb_unc)
    fitUncs['Zc_Zb'][s].append(-1.0*zc_zb_unc)
    rs_fits['Zc_Zb'][s].append(rs_fits['Zc'][s][0]/rs_fits['Zb'][s][0])
    #print '>>>>Zc_Zb unc. ', s, ' ', fitUncs['Zc_Zb'][s]
    
    if s != 'Exp':
      zc_unc = [fitUncs['Zc'][s][2],fitUncs['Zc'][s][3]]
      zb_unc = [fitUncs['Zb'][s][2],fitUncs['Zb'][s][3]]
      zc_zb_unc = cal_error_rcrb_mllFit(rs['Zc'],rs['Zb'],zc_unc,zb_unc,cov['Zc'][s][1],cov['Zb'][s][1],mode)
      fitUncs['Zc_Zb'][s].append(zc_zb_unc)
      fitUncs['Zc_Zb'][s].append(-1.0*zc_zb_unc)
      rs_fits['Zc_Zb'][s].append(rs_fits['Zc'][s][1]/rs_fits['Zb'][s][1])
      #print '>>>>Zc_Zb unc. ', s, ' ', fitUncs['Zc_Zb'][s]

  print '>>>>>>>>>>>Unc. and central value from fit (read in log files) and total unc from best fit >>>>>>>>>>>'
  print '>> Central values from best fit >>'
  for k,v in rs.items():
    print '..', k, ' ', v
  print '>> Uncertainties from best fit >>'
  for k,v in totUncs.items():
    print '..', k, ' ', v
  print '>> Central values from uncertainty breakdown fits>>'
  for k,v in rs_fits.items():
    print '..', k
    for k1,v1 in v.items():
      print '....', k1, ' ', v1
  print '>> Unc. values from uncertainty breakdown fits>>'
  for k,v in fitUncs.items():
    print '..', k
    for k1,v1 in v.items():
      print '....', k1, ' ', v1

  print '#####Compare covariant for Zc and Zb #######'
  for s in ['MC_stat_bin','JEC_gr','JER_gr','gcc_gr','gbb_gr', 'pu_gr', 'Bkgr_xSec', 'Exp']:
    print s, ' ', cov['Zc'][s][0][0][1], ' ', cov['Zb'][s][0][0][1]
    if s != 'Exp':
      print s, ' ', cov['Zc'][s][1][0][1], ' ', cov['Zb'][s][1][0][1]
  
  print '#####Compare covariant for Zl and Zc and Zb #######'
  for s in ['MC_stat_bin','JEC_gr','JER_gr','gcc_gr','gbb_gr', 'pu_gr', 'Bkgr_xSec', 'Exp']:
    print s, ' Zl_Zc, Zl_Zb ', cov['Zl'][s][0][0][1], ' ', cov['Zl'][s][0][0][2] #Zl,Zc ; Zl, Zb
    print s, ' Zc_Zl, Zb_Zl ', cov['Zc'][s][0][0][2], ' ', cov['Zb'][s][0][0][2]
    if s != 'Exp':
      print s, ' Zl_Zc, Zl_Zb ', cov['Zl'][s][1][0][1], ' ', cov['Zl'][s][1][0][2] #Zl,Zc ; Zl, Zb
      print s, ' Zc_Zl, Zb_Zl ', cov['Zc'][s][1][0][2], ' ', cov['Zb'][s][1][0][2]

  print '######Correlation of Zl and Zc and Zb in free gcc_gr or gcc_gb fit#########'
  print 'gcc_gr ', cor['Zl']['gcc_gr']
  print 'gbb_gr ', cor['Zl']['gbb_gr']
  print 'pu_gr ', cor['Zl']['pu_gr']
  print 'Bkgr_xSec', cor['Zl']['Bkgr_xSec']


  #print cov['Zc']
  #print cov['Zb']
  
  ##########
  #calculate uncertainty
  ##########
  print '######################'
  print 'Calculate uncertainty'
  print '######################'
  
  uncs = {'Zc':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
          'Zb':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
          'Zl':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
          'Zc_Zb':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]}}
  
  print 'Use relative uncertainty: ', useRelUnc
  for flav,fitUnc in fitUncs.items():
    #if flav == 'Zc_Zb': continue
    for n,v in fitUnc.items():
      if n != 'Exp':
        #method 1: tot - freeze syst
        #down
        a,b = totUncs[flav][0],v[0]
        if useRelUnc:
          a,b = totUncs[flav][0]/rs[flav],v[0]/rs_fits[flav][n][0] #relative, 0 = freeze that unc. fit, 1 not freeze that unc. fit
        if a > b:
          if not useRelUnc: uncs[flav][n].append(sqrt(a*a-b*b))
          if useRelUnc:  uncs[flav][n].append(sqrt(a*a-b*b)*rs[flav])
        else:
          print '>>>>' + flav + ' Warning Tot unc. "smaller" than unc. ', n, ' ', a, ' ', b
          uncs[flav][n].append(0)
        #up
        a,b = totUncs[flav][1],v[1]
        if useRelUnc:
          a,b = totUncs[flav][1]/rs[flav],v[1]/rs_fits[flav][n][0] #relative
        if a < b:
          if not useRelUnc: uncs[flav][n].append(sqrt(a*a-b*b))
          if useRelUnc: uncs[flav][n].append(sqrt(a*a-b*b)*rs[flav])
        else:
          print '>>>>' + flav + ' Warning Tot unc. "smaller" than unc. ', n, ' ', a, ' ', b
          uncs[flav][n].append(0)
        
        #Method 2: not freeze syst - stat
        #down
        a,b = v[2],fitUnc['Exp'][0]
        if useRelUnc:
          a,b = v[2]/rs_fits[flav][n][1],fitUnc['Exp'][0]/rs_fits[flav]['Exp'][0]
        if a > b:
          if not useRelUnc: uncs[flav][n].append(sqrt(a*a-b*b))
          if useRelUnc: uncs[flav][n].append(sqrt(a*a-b*b)*rs[flav])
        else:
          print '>>>>' + flav + ' Warning Unc. "smaller" than stat. ', n, ' ', a, ' ', b
          uncs[flav][n].append(0)
        
        if abs(v[2]) > 0.9*rs[flav]:
          print '>>>Unreasonable large unc. when not freezing ', n, '. set unc. to 0 ', v[2], ' ', rs[flav] 
          uncs[flav][n][-1] = 0
        #up
        a,b = v[3],fitUnc['Exp'][1]
        if useRelUnc:
          a,b = v[3]/rs_fits[flav][n][1],fitUnc['Exp'][1]/rs_fits[flav]['Exp'][0]

        if a < b:
          if not useRelUnc: uncs[flav][n].append(sqrt(a*a-b*b))
          if useRelUnc: uncs[flav][n].append(sqrt(a*a-b*b)*rs[flav])
        else:
          print '>>>>' + flav + ' Warning Unc. "smaller" than stat. ', n, ' ', a, ' ', b
          uncs[flav][n].append(0)
        
        if abs(v[3]) > 0.9*rs[flav]:
          print '>>>Warning: ' + flav + ' Unreasonable large unc. when not freezing ', n, '. set unc. to 0 ', v[3], ' ', rs[flav]
          uncs[flav][n][-1] = 0

      else:
        if not useRelUnc: uncs[flav]['Exp'] = fitUncs[flav]['Exp'][:]
        if useRelUnc:
          tmp = fitUncs[flav]['Exp'][:]
          for iTmp in tmp:
            uncs[flav]['Exp'].append(iTmp/rs_fits[flav]['Exp'][0]*rs[flav])


  
  print '>>>>>>>> Final uncertainty >>>>>>>>>>'
  print uncs
  
  uncs_final = {'Zc':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
          'Zb':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
          'Zl':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]},
          'Zc_Zb':{'MC_stat_bin':[],'JEC_gr':[], 'JER_gr':[], 'gcc_gr':[], 'gbb_gr':[], 'pu_gr':[], 'Bkgr_xSec':[], 'Exp':[]}}
  
  for flav,unc in uncs.items():
    if flav == 'Zc_Zb': continue
    for s,v in unc.items():
      if s != 'Exp':
        if (v[0] != 0 and v[1] != 0):
          uncs_final[flav][s].append(max([v[0],v[1]]))
        elif (v[0] == 0 or v[1] == 0) and (v[2] != 0 and v[3] != 0):
          print '>>>Unc. from method 1 results 0, use method 2 ', flav, ' ', s, ' ', v 
          uncs_final[flav][s].append(max([v[2],v[3]]))
        else:
          print '>>>Unc. from both method 1,2 results 0, use max ', flav, ' ', s, ' ', v 
          uncs_final[flav][s].append(max(v))
      
      else: # if s != 'Exp'
        uncs_final[flav]['Exp'] = [max(uncs[flav]['Exp'])]
  
  #calculate Zc_Zb ratio by using the correlation from not freeze the syst 
  print '>>>>>>>>>Calculate the ratio of Zc_Zb using the correlation from not freeze the syst>>>>>'
  for s in uncs_final['Zc_Zb'].keys():
    print '>>>>>>>>>>>', s
    if s != 'Exp':
      #zc_zb_unc = cal_error_rcrb_mllFit(rs['Zc'],rs['Zb'],uncs_final['Zc'][s],uncs_final['Zb'][s],cov['Zc'][s][1],cov['Zb'][s][1])
      zc_zb_unc = cal_error_rcrb_mllFit(rs['Zc'],rs['Zb'],uncs_final['Zc'][s],uncs_final['Zb'][s],cor['Zc'][s][1],cor['Zb'][s][1],2)
      uncs_final['Zc_Zb'][s] = [zc_zb_unc]
    else: 
      zc_zb_unc = cal_error_rcrb_mllFit(rs['Zc'],rs['Zb'],uncs_final['Zc'][s],uncs_final['Zb'][s],0,0,0)
      uncs_final['Zc_Zb'][s] = [zc_zb_unc]
  
  print '>>>>>>>> Uncertainty final>>>>>>>>>>'
  sumUncs = {'Zc':0,'Zb':0,'Zl':0,'Zc_Zb':0}
  systUncs = {'Zc':0,'Zb':0,'Zl':0,'Zc_Zb':0}
  for flav,unc in uncs_final.items():
    for s,v in unc.items():
      sumUncs[flav] = sumUncs[flav] + v[0]*v[0]
      if s != 'Exp': 
        systUncs[flav] = systUncs[flav] + v[0]*v[0]

  for flav in sumUncs.keys():
    sumUncs[flav] = sqrt(sumUncs[flav])
    systUncs[flav] = sqrt(systUncs[flav])
  
  print '######## Uncertainty summary ########'
  uncNames = ['Exp','JEC_gr','JER_gr','MC_stat_bin','gcc_gr','gbb_gr', 'pu_gr', 'Bkgr_xSec']
  for flav in ['Zc','Zb','Zl','Zc_Zb']:
    print '>>>>>>>>>>>', flav, '>>>>>>>>>>>>'
    for n in uncNames:
      print n, ': ', uncs_final[flav][n][0]
    print 'Total unc: ', sumUncs[flav], ' ', totUncs[flav][0], ' ', totUncs[flav][1]

  print '>>>>>>> Total unc>>>>>>>>>>>>>>>>>'
  print totUncs
  print sumUncs

  print '>>>>>>> SF >>>>>>>>>>>>>>>>>>>>'
  print rs
  rs_unc = {}
  for flav in ['Zc','Zb','Zl','Zc_Zb']:
    #print rs[flav]
    #print totUncs[flav][0], totUncs[flav][1]
    #print uncs_final[flav]['Exp'][0]
    rs_unc[flav] = getUnc_r(rs[flav],max(abs(totUncs[flav][0]),abs(totUncs[flav][1])),uncs_final[flav]['Exp'][0]) # [stat, syst] 
  
  print '######## write the SF to latex file########'
  fLatex.write('\\begin{table}\n')
  fLatex.write('\\begin{tabular}{cc}\n')
  fLatex.write('\\hline\n')
  fLatex.write('\\hline\n')
  l = ['', chan]
  st = makeLatexLine(l)
  fLatex.write(st)
  fLatex.write('\hline\n')
  st = 'SF$_c$ & ' + '{0:20}'.format('%.2f' % rs['Zc']) + ' $\pm$ ' + '{0:20}'.format('%.2f (stat.)' % rs_unc['Zc'][0]) + ' $\pm$ ' + '{0:20}'.format('%.2f (syst.)' % rs_unc['Zc'][1]) + '\\\\\n' 
  fLatex.write(st)
  st = 'SF$_b$ & ' + '{0:20}'.format('%.2f' % rs['Zb']) + ' $\pm$ ' + '{0:20}'.format('%.2f (stat.)' % rs_unc['Zb'][0]) + ' $\pm$ ' + '{0:20}'.format('%.2f (syst.)' % rs_unc['Zb'][1]) + '\\\\\n'
  fLatex.write(st)
  st = 'SF$_l$ & ' + '{0:20}'.format('%.2f' % rs['Zl']) + ' $\pm$ ' + '{0:20}'.format('%.2f (stat.)' % rs_unc['Zl'][0]) + ' $\pm$ ' + '{0:20}'.format('%.2f (syst.)' % rs_unc['Zl'][1]) + '\\\\\n'
  fLatex.write(st)
  fLatex.write('\\hline\n')
  fLatex.write('\\hline\n')
  fLatex.write('\\end{tabular}\n')
  fLatex.write('\\end{table}\n')
  fLatex.write('\n')


  print '######## Write unc. to latex file##########'
  uncs_percent = copy.deepcopy(uncs_final)
  sumUncs_percent = copy.deepcopy(sumUncs)
  systUncs_percent = copy.deepcopy(systUncs)
  for flav,v in uncs_percent.items():
    sumUncs_percent[flav] = 100*sumUncs_percent[flav]/rs[flav]
    systUncs_percent[flav] = 100*systUncs_percent[flav]/rs[flav]
    for s in v.keys():
      uncs_percent[flav][s][0] = 100*uncs_percent[flav][s][0]/rs[flav]
  print uncs_percent
  print sumUncs_percent
  print systUncs_percent
  
  #err_names = [['JEC_gr','JEC'], ['JER_gr','JER'], ['MC_stat_bin','MC temp. stat.'], ['gcc_gr','gcc splitting'], ['gbb_gr','gbb splitting'], ['Bkgr_xSec', 'Bkgr. sub.'], ['Total_syst','Total syst.'], ['Exp', 'Stat.'], ['Total','Total']]
  err_names = [['JEC_gr','JEC'], ['JER_gr','JER'], ['MC_stat_bin','MC temp. stat.'], ['gcc_gr','gcc splitting'], ['gbb_gr','gbb splitting'], ['pu_gr','PU reweight'], ['Bkgr_xSec', 'Bkgr. sub.']]
  fLatex.write('\\begin{table}\n')
  #fLatex.write('\\begin{tabular}{ccccc}\n')
  fLatex.write('\\begin{tabular}{cccc}\n')
  fLatex.write('\\hline\n')
  fLatex.write('\\hline\n')
  #l = ['', 'SF$_c$ ', 'SF$_b$ ', 'SF$_l$ ', 'SF$_c$/SF$_b$']
  l = ['', 'SF$_c$ ', 'SF$_b$ ', 'SF$_l$ ']
  st = makeLatexLine(l)
  fLatex.write(st)
  fLatex.write('\hline\n')
  roundUp = 0.001 # to round up 0.049-> 0.1
  for err in err_names:
    l = [err[1], uncs_percent['Zc'][err[0]][0], uncs_percent['Zb'][err[0]][0] + roundUp, uncs_percent['Zl'][err[0]][0]]
    #l = []
    #if err[0] != 'Total_syst' and err[0] != 'Total':
    #  l = [err[1], uncs_percent['Zc'][err[0]][0], uncs_percent['Zb'][err[0]][0], uncs_percent['Zl'][err[0]][0], uncs_percent['Zc_Zb'][err[0]][0]]
    #if err[0] == 'Total_syst':
    #  fLatex.write('\hline\n')
    #  l = [err[1], systUncs_percent['Zc'], systUncs_percent['Zb'], systUncs_percent['Zl'], systUncs_percent['Zc_Zb']]
    #if err[0] == 'Total':
    #  l = [err[1], sumUncs_percent['Zc'], sumUncs_percent['Zb'], sumUncs_percent['Zl'], sumUncs_percent['Zc_Zb']]

    st = makeLatexLine(l,'%.1f')
    fLatex.write(st)

  fLatex.write('\\hline\n')
  fLatex.write('\\hline\n')
  fLatex.write('\\end{tabular}\n')
  fLatex.write('\\end{table}\n')
  fLatex.write('\n')

  #############write stuffs to estimate ratio###########
  for flav in ['Zc','Zb', 'Zc_Zb']:
      st = chan + ' ' + flav + ' ' + str(rs[flav]) + ' ' + str(rs_unc[flav][0]) + ' ' + str(rs_unc[flav][1]) + ' ' + str(uncs_final[flav]['MC_stat_bin'][0]) + '\n'
      fOut.write(st)

  for flav in ['Zc', 'Zb', 'Zl']:
      for i in range(3):
         st = chan + ' ' + flav + '_corr_' + str(i) + ' ' + str(cor_central[flav][i][0]) + ' ' + str(cor_central[flav][i][1]) + ' ' + str(cor_central[flav][i][2]) + '\n'
         fOut.write(st)
   
  for i in range(3):
    st = chan + ' ' + 'Zc_cov_' + str(i) + ' ' + str(cov_zc[i][0]) + ' ' + str(cov_zc[i][1]) + ' ' + str(cov_zc[i][2]) + '\n'
    fOut.write(st)

  for i in range(3):
    st = chan + ' ' + 'Zb_cov_' + str(i) + ' ' + str(cov_zb[i][0]) + ' ' + str(cov_zb[i][1]) + ' ' + str(cov_zb[i][2]) + '\n'
    fOut.write(st)

  ######################################
  print ">>>>>>>>>>>>>>>Check covariance>>>>>>>>>>>>>>"
  print cov_zc
  print cov_zb

  
  os.chdir('../')
