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
  print '>>> V_00, V_11, V_22: ', sqrt(V[0][0]), ' ', sqrt(V[1][1]), ' ', sqrt(V[2][2])
  N = r1*n1 + r2*n2 + r3*n3
  N1 = r1*n1
  N23 = r2*n2 + r3*n3
  dev_r1 = n1*N23/(N*N)
  dev_r2 = -N1*n2/(N*N)
  dev_r3 = -N1*n3/(N*N)
  sig = dev_r1*dev_r1*V[0][0] + dev_r2*dev_r2*V[1][1] + dev_r3*dev_r3*V[2][2] + 2*dev_r1*dev_r2*V[0][1] + 2*dev_r1*dev_r3*V[0][2] + 2*dev_r2*dev_r3*V[1][2]
  return sqrt(sig)
  
def cal_error_rcrb(r1, err_r1, r2, err_r2, V, useCovariant = True): #y = r1/r2 r1 = rc, r2 = rb
  print '>>> V_00, V_11, V_22: ', sqrt(V[0][0]), ' ', sqrt(V[1][1]), ' ', sqrt(V[2][2])
  dev_r1 = 1/r2
  dev_r2 = -r1/(r2*r2)
  if useCovariant:
    sig = dev_r1*dev_r1*err_r1*err_r1 + dev_r2*dev_r2*err_r2*err_r2 + 2*dev_r1*dev_r2*V[0][1]
  else:
    sig = dev_r1*dev_r1*err_r1*err_r1 + dev_r2*dev_r2*err_r2*err_r2
  if sig < 0:
    print '>>>>>>>Warning: Negative found when using covariant for rc/rb error, use non-covariant', sig, ' ', V[0][1]
    sig = dev_r1*dev_r1*err_r1*err_r1 + dev_r2*dev_r2*err_r2*err_r2

  return sqrt(sig)


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

def checkFitProblem(fIn, fitType):
  if fitType == 'mll':
    pipe = os.popen('grep "Fit failed" ' + fIn).readlines()
    if len(pipe) > 0:
      print '>>>>Waring: Mll fit fails ', fIn
      return False

  if fitType == 'md':
    lines = open(fIn).readlines()
    pos1 = -1
    pos2 = -1
    for i in range(len(lines)-1,-1,-1):
      if '--- MultiDimFit ---' in lines[i]: pos2 = i
      if 'r1' in lines[i] and '=' in lines[i] : pos1 = i
      if pos1 != -1 and pos2 != -1:
        break

    if pos2 - pos1 > 5:
      print '>>>>Warning: MultiDim fit do not converge'
      return False

  return True
 

#############################
# Main program
#############################

skipGoodnessOfFit = True
skipImpact = False 
skipUncBreakdown = False

dcName_pre = 'test_vtxMass_pt20_allData_binStatThresh_0_useDDbjet_useDDljet_customBinning'
chans = ['Zee','Zmm']
#chans = ['Zee']

#makeWS = 'text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:testModel'
makeWS_zc = 'text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO \'map=.*/Zc:r1[1,0,3]\' --PO \'map=.*/Zb:r2[1,0,3]\' --PO \'map=.*/Zl:r3[1,0,3]\''
makeWS_zb = 'text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO \'map=.*/Zb:r1[1,0,3]\' --PO \'map=.*/Zc:r2[1,0,3]\' --PO \'map=.*/Zl:r3[1,0,3]\''
makeWS_zl = 'text2workspace.py -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO \'map=.*/Zl:r1[1,0,3]\' --PO \'map=.*/Zc:r2[1,0,3]\' --PO \'map=.*/Zb:r3[1,0,3]\''

multiDimFit = 'combine -M MultiDimFit'
#multiDimFit_opt = '-v 3 -S 0' # verbose, -S 0 for no nuisance parameter (temporary for now)

#mllFit = 'combine -M MaxLikelihoodFit --robustFit 1 --minimizerTolerance 0.00000001'
mllFit = 'combine -M MaxLikelihoodFit'
#mllFit_opt = '-v 3 -S 0' # verbose, -S 0 for no nuisance parameter (temporary for now)

#--robustFit arg (=0)                  Search manually for 1 and 2 sigma bands instead of using Minos
#--stepSize (the default is 0.1, and is relative to the range of the parameter)
#--minimizerTolerance can be varied between 0.0001 and 1.0 or more (the default is 0.01); smaller tolerances can give a more accurate result and avoid false minima, while larger tolerances can allow a complex fit to converge to a satisfactory result rather than giving up.
#--minimizerStrategy can be set to 1 (default), 0 or 2; larger values correspond to using more derivatives in the fit, which is slower but can give better or more stable results for smooth pdfs (and can result in failures for less smooth pdfs)

fit_opt = '--robustFit 1 --stepSize 0.01 --minimizerStrategy 1 -v 3 --minimizerTolerance' # verbose, -S 0 for no nuisance parameter (temporary for now)
#fit_opt_setRange = {'Zee': '--setPhysicsModelParameterRanges r1=0.5,1.5:r2=0.5,1.5:r3=0.8,1.6', 'Zmm':'--setPhysicsModelParameterRanges r1=0.5,1.5:r2=0.5,1.5:r3=0.5,1.5', 'Com':'--setPhysicsModelParameterRanges r1=0.5,1.5:r2=0.5,1.5:r3=0.5,1.5'}
fit_opt_setRange = {'Zee': '--setPhysicsModelParameterRanges r1=0.7,1.5:r2=0.5,1.5:r3=0.6,1.8', 'Zmm':'--setPhysicsModelParameterRanges r1=0.5,1.5:r2=0.5,1.5:r3=0.5,1.5', 'Com':'--setPhysicsModelParameterRanges r1=0.5,1.5:r2=0.5,1.5:r3=0.5,1.5'}
#fit_opt1 = '--robustFit 1 --minimizerTolerance 0.00001 --minimizerStrategy 1 -v 3' # verbose, -S 0 for no nuisance parameter (temporary for now)
#tolerances = ['0.0000001','0.000001','0.00001','0.0001']
tolerances = ['0.0001','0.001','0.01']

label = {'r1':'Zc', 'r2':'Zb', 'r3':'Zl', 'r1_over_r2':'Zc_over_Zb'}

#rate, cov, cor = getResults('log_test_allData_multiDimFit.txt', 'log_test_allData_mllFit.txt')
#rate_prefit = getPreFit(dcName)

#for k,v in label.items():
#  rate[v] = rate.pop(k)

#print rate
#print cov
#print cor
#print rate_prefit

#printResults(rate_prefit, rate, cov)

#combine datacard
if 'Com' in chans:
  cmd = 'combineCards.py'
  chanTmps = ['Zee','Zmm','Com']
  for i in range(len(chanTmps)-1):
    dcName = dcName_pre + '_' + chanTmps[i] + '.txt'
    cmd = cmd + ' SR_' + chanTmps[i] + '=' + dcName
  cmd = cmd + ' > ' + dcName_pre + '_' + chanTmps[-1] + '.txt'
  print '>>>>>>>', cmd
  os.system(cmd)

for chan in chans:
  print '#########################################################'
  print '#', chan
  print '#########################################################'
  
  dcName = dcName_pre + '_' + chan

  ###########################
  #Make workspace
  ###########################

  if makeWS_zc != '':
    os.system(makeWS_zc + ' -o Zc_' + dcName + '_ws.root ' + dcName + '.txt')
  if makeWS_zb != '':
    os.system(makeWS_zb + ' -o Zb_' + dcName + '_ws.root ' + dcName + '.txt')
  if makeWS_zl != '':
    os.system(makeWS_zl + ' -o Zl_' + dcName + '_ws.root ' + dcName + '.txt')

  ###########################
  # create workdir, subdir ...
  ###########################
  wDir = dcName
  os.system('mkdir ' + wDir)
  os.system('rm ' + wDir + '/* -r')
  os.system('mkdir ' + wDir + '/Logs')
  os.system('cp ' + dcName + '.txt ' + wDir)
  os.system('cp ' + 'Zc_' + dcName + '_ws.root ' + wDir)
  os.system('cp ' + 'Zb_' + dcName + '_ws.root ' + wDir)
  os.system('cp ' + 'Zl_' + dcName + '_ws.root ' + wDir)

  os.chdir(wDir)
  fLatex = open(dcName + '_newErrorEstimation_latex.txt','w')

  ###########################
  #do fit
  ###########################
  
  print '####################'
  print '#', chan, ': Fitting'
  print '####################'
  #tor_Zc = '0.01'
  for z in ['Zc','Zb','Zl']:
    multiDimFit_log = 'Logs/log_' + z + '_multiDimFit.txt'
    mllFit_log = 'Logs/log_' + z + '_mllFit.txt'
    mllFit_log_tmp = 'Logs/log_' + z + '_mllFit_tmp.txt'
    
    torTmp = 0.01
    for tor in tolerances:
      opt = fit_opt + ' ' + tor
      cmd = multiDimFit + ' ' + opt + ' --saveWorkspace -n _' + z + '_bestfit ' + z + '_' + dcName + '_ws.root >& ' + multiDimFit_log
      print '>>>>>>>Best fit MultiDimFit: ', cmd #to get snapshot
      os.system('rm ' + multiDimFit_log)
      os.system(cmd)
      if checkFitProblem(multiDimFit_log,'md'): 
        print cmd + ': tolerant ' + str(tor)
        torTmp = tor
        break

    #do multiDimFit with narrow bands for systematic breakdown
    multiDimFit_log_1 = 'Logs/log_' + z + '_multiDimFit_forSystBreakdown.txt'
    opt = fit_opt + ' ' + torTmp + ' ' + fit_opt_setRange[chan]
    cmd = multiDimFit + ' ' + opt + ' --saveWorkspace -n _' + z + '_bestfit_forSystBreakdown ' + z + '_' + dcName + '_ws.root >& ' + multiDimFit_log_1
    print '>>>>>>>Best fit MultiDimFit for syst. breakdown: ', cmd
    os.system('rm ' + multiDimFit_log_1)
    os.system(cmd)
    if not checkFitProblem(multiDimFit_log_1,'md'): 
      print 'Warning: MultidimFit for syst. breakdown problem, check!!! ', cmd + ': tolerant ' + str(torTmp)
    
    for tor in tolerances:
      opt = fit_opt + ' ' + tor + ' ' + fit_opt_setRange[chan]
      cmd = mllFit + ' ' + opt + ' -n _' + z + '_bestfit ' + z + '_' + dcName + '_ws.root >& ' + mllFit_log_tmp
      print '>>>>>>>Best fit MllFit:      ', cmd #to get total errors
      os.system('rm ' + mllFit_log_tmp)
      os.system(cmd)
      if checkFitProblem(mllFit_log_tmp,'mll'):
        print cmd + ': tolerant ' + str(tor)
        #if z == 'Zc': tor_Zc = str(tor) #use for impact fit
        break

    cmd = 'grep -v "Iteration #" ' + mllFit_log_tmp + ' > ' + mllFit_log #remove long iteration message
    os.system(cmd)


  if not skipGoodnessOfFit: 
    print '####################'
    print '#', chan, ': GoodnessOfFit'
    print '####################'

    goodnessOfFit_log = 'Logs/log_' + dcName + '_GoodnessOfFit.txt'
    fName1 = 'GoodnessOfFit_' + dcName + '.root' 
    fName2 = 'GoodnessOfFit_toy_' + dcName + '.root' 
    os.system('combine -M GoodnessOfFit Zc_' + dcName + '_ws.root --algo=saturated -n ' + fName1.replace('.root','') + ' >& ' + goodnessOfFit_log)
    os.system('combine -M GoodnessOfFit Zc_' + dcName + '_ws.root --algo=saturated -t 500 -n ' + fName2.replace('.root','') + ' >> ' + goodnessOfFit_log)
    plotGoodnessOfFit(fName1,fName2,dcName)
    os.system('cp goodnessOfFit_*  ../Plots/')

  if not skipImpact: 
    print '####################'
    print '#', chan, ': Impact'
    print '####################'

    impact_log = 'Logs/log_' + dcName + '_impact.txt'
    #os.system('combineTool.py -M Impacts -d ' + dcName + '.root -m 125 --doInitialFit --robustFit 1 -n ' + dcName)
    #os.system('combineTool.py -M Impacts -d ' + dcName + '.root -m 125 --doFits --robustFit 1 --parallel 4 -n ' + dcName)
    #os.system('combineTool.py -M Impacts -d Zc_' + dcName + '_ws.root -m 125 --doInitialFit --robustFit 1 ')
    #os.system('combineTool.py -M Impacts -d Zc_' + dcName + '_ws.root -m 125 --doFits --robustFit 1 '
    os.system('combineTool.py -M Impacts -d Zc_' + dcName + '_ws.root -m 125 --doInitialFit --robustFit 1 --stepSize 0.01 ' + fit_opt_setRange[chan] + ' ')
    os.system('combineTool.py -M Impacts -d Zc_' + dcName + '_ws.root -m 125 --doFits --robustFit 1 --stepSize 0.01 --parallel 5 ' + fit_opt_setRange[chan] + ' ')
    
    #os.system('combineTool.py -M Impacts -d Zc_' + dcName + '_ws.root -m 125 --doInitialFit ' + fit_opt + ' ' + tor_Zc) 
    #os.system('combineTool.py -M Impacts -d Zc_' + dcName + '_ws.root -m 125 --doFits --parallel 5 ' + fit_opt + ' ' + tor_Zc)
    
    os.system('combineTool.py -M Impacts -d Zc_' + dcName + '_ws.root -m 125 -o ' + dcName + '.json')
  
  if not skipUncBreakdown:

    print '####################'
    print '#', chan, ': Uncertainty breakdown'
    print '####################'

    #redo MultiDimFit with the same settings as mll fit
    '''good_tor = {'Zc':[],'Zb':[],'Zl':[]}
    for z in ['Zc','Zb','Zl']:
      multiDimFit_log = 'Logs/log_' + z + '_multiDimFit_forSystBreakdown.txt'
      for tor in tolerances:
        opt = fit_opt + ' ' + tor + ' ' + fit_opt_setRange[chan]
        cmd = multiDimFit + ' ' + opt + ' --saveWorkspace -n _' + z + '_bestfit_forSystBreakdown ' + z + '_' + dcName + '_ws.root >& ' + multiDimFit_log
        print '>>>>>>>Best fit MultiDimFit for syst breakdown: ', cmd #to get snapshot
        os.system('rm ' + multiDimFit_log)
        os.system(cmd)
        if checkFitProblem(multiDimFit_log,'md'): 
          print cmd + ': tolerant ' + str(tor)
          good_tor[z].append(tor)
          break
    print '>> Found good tor: ', good_tor
    '''

    systs = ['MC_stat_bin','JEC_gr', 'JER_gr', 'gcc_gr', 'gbb_gr', 'pu_gr', 'Bkgr_xSec', 'Exp'] #syst when freezing all syst group except the one in question, Exp == freezing all nuisance syst --> this = stat uncertainty
    for k in systs:
      print '#############'
      print 'Breakdown syst: ', k
      print '#############'
      for z in ['Zc','Zb','Zl']:
        logPre = 'Logs/log_' + z + '_syst_mllFit'
        tmpS = ''
        freeze_syst_groups = ''
        if k != 'Exp':
          addComma = False
          for s in systs:
            if s != 'Exp' and s != k:
              if not addComma:
                addComma = True
                freeze_syst_groups += s
              else: 
                freeze_syst_groups += ',' + s
          
          #not_freeze fit
          tmpS = 'not_freeze_' + k
          log = logPre + '_' + tmpS + '.txt'
          logTmp = logPre + '_' + tmpS + '_tmp.txt'
          for tor in tolerances:
            opt = fit_opt + ' ' + tor + ' ' + fit_opt_setRange[chan]
            cmd = mllFit + ' ' + opt + ' --freezeNuisanceGroups ' + freeze_syst_groups + ' -d higgsCombine_' + z + '_bestfit_forSystBreakdown.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _' + z + '_' + tmpS + ' >& ' + logTmp
            print '>>>>> Syst breakdown ', cmd
            os.system('rm ' + logTmp)
            os.system(cmd)
            if checkFitProblem(logTmp,'mll'):
              print cmd + ': tolerant ' + str(tor)
              break
          os.system('grep -v "Iteration #" ' + logTmp + ' > ' + log)

          #freeze fit
          tmpS = 'freeze_' + k
          log = logPre + '_' + tmpS + '.txt'
          logTmp = logPre + '_' + tmpS + '_tmp.txt'
          for tor in tolerances:
            opt = fit_opt + ' ' + tor + ' ' + fit_opt_setRange[chan]
            cmd = mllFit + ' ' + opt + ' --freezeNuisanceGroups ' + k + ' -d higgsCombine_' + z + '_bestfit_forSystBreakdown.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _' + z + '_' + tmpS + ' >& ' + logTmp
            print '>>>>> Syst breakdown ', cmd
            os.system('rm ' + logTmp)
            os.system(cmd)
            if checkFitProblem(logTmp,'mll'):
              print cmd + ': tolerant ' + str(tor)
              break
          os.system('grep -v "Iteration #" ' + logTmp + ' > ' + log)
        
        #freeze Exp fit
        else:
          log = logPre + '_freeze_' + k + '.txt'
          logTmp = logPre + '_freeze_' + k + '_tmp.txt'
          tmpS = 'freeze_' + k
          for tor in tolerances:
            opt = fit_opt + ' ' + tor + ' ' + fit_opt_setRange[chan]
            cmd = mllFit + ' ' + opt + ' --freezeNuisanceGroups ' + k + ' -d higgsCombine_' + z + '_bestfit_forSystBreakdown.MultiDimFit.mH120.root -w w --snapshotName "MultiDimFit" -n _' + z + '_' + tmpS + ' >& ' + logTmp
            print '>>>>> Syst breakdown ', cmd
            os.system('rm ' + logTmp)
            os.system(cmd)
            if checkFitProblem(logTmp,'mll'):
              print cmd + ': tolerant ' + str(tor)
              break
          os.system('grep -v "Iteration #" ' + logTmp + ' > ' + log)


  os.chdir('../')
