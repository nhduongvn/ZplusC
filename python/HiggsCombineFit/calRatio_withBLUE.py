import os,sys,ROOT
from math import *
sys.path.append('/uscms_data/d3/duong/CMSSW/CMSSW_7_6_5/src/ZplusC/python/')
from myutils.util_funcs import *
from myutils.blue import *

def combineRatio(x,y,xerr,yerr,useAverage=False): #combine ratios x/y of electron and muon channels, for example event count ratio N_c/N_b
  #when error == 0 return the average
  if useAverage:
    s = 0
    for i in range(len(x)):
      s = s + x[i]/y[i]
    return [s/len(x),0]
  else:
    x_y = []
    x_y_err = []
    for i in range(len(x)):
      x_y.append(x[i]/y[i])
      x_y_err.append(x[i]/y[i]*sqrt(x_err[i]*x_err[i]/(x[i]*x[i]) + y_err[i]*y_err[i]/(y[i]*y[i])))
    
    #com = sum(xi/xi_err^2)/sum(1/xi_err^2)
    num = 0
    deno = 0
    for i in range(len(x)):
      num = num + x_y[i]/(x_y_err[i]*x_y_err[i])
      deno = deno + 1./(x_y_err[i]*x_y_err[i])
    
    return [num/deno, sqrt(1./deno)]

def getMultiplicationRelative(x,x_unc):
  s = 0
  for i in range(len(x)):
    s = s + x_unc[i]*x_unc[i]/(x[i]*x[i])
  return sqrt(s)

##################################################
#Main program
##################################################

#fR = '../summary_analyzeFit_pt20.txt'
#fN = '../summary_eventCout_pt20.txt'

fR = '../summary_analyzeFit_vtxMass_pt20_binStatThresh_0_useDDbjet_useDDljet_customBinning.txt'
fN = '../summary_eventCout_vtxMass_pt20_binStatThresh_0_useDDbjet_useDDljet_customBinning.txt'

fOutName = 'summary_all_vtxMass_pt20_binStatThresh_0_useDDbjet_useDDljet_customBinning.txt'

#eff_b = [0.61143, 0.00909]
#eff_c = [0.12162, 0.00242]

#eff_b = [0.585, 0.006, 0.00784]
#eff_c = [0.114, 0.00148, 0.00465]
eff_b = [0.5459, 0.00549, 0.00732]
eff_c = [0.0949, 0.00134, 0.00390]

A = {'Zee':{'Zb':[1,0.01],
            'Zc':[1,0.01],
            'Zinc':[1,0.01]
           },
     'Zmm':{'Zb':[1,0.01],
            'Zc':[1,0.01],
            'Zinc':[1,0.01]
           },
     'Com':{'Zb':[1,0.01],
            'Zc':[1,0.01],
            'Zinc':[1,0.01]}}

sf = {'Zee':{},'Zmm':{}, 'Com':{}}

N = {'Zee':{},'Zmm':{}}

chans = ['Zee','Zmm','Com']

##########
#Get SF and number of event
##########

lR = open(fR,'r').readlines()
lN = open(fN,'r').readlines()

for chan in chans:
  for flav in ['Zc','Zb','Zc_Zb']:
    for l in lR:
      lsp = l.split()
      if chan in lsp and flav in lsp:
        sf[chan][flav] = [float(lsp[2]),float(lsp[3]),float(lsp[4]),float(lsp[5])] #central, stat, syst, syst from MC_stat_bin 
    for l in lN:
      lsp = l.split()
      if flav == 'Zc_Zb' or chan == 'Com': continue
      if chan in lsp and flav in lsp:
        N[chan][flav] = [float(lsp[2]),float(lsp[3])]
      if chan in lsp and 'Data_inc' in lsp:
        N[chan]['Data_inc'] = [float(lsp[2]),float(lsp[3])]


print sf
print N

##############
#Estimate event ratio
##############

N_rat = {'Zee':{},'Zmm':{},'Com':{}} #chan,R(c/b),R(c/j),R(b/j), ratio with error

for chan in chans:
  if chan != 'Com':
    N_rat[chan]['Rcb'] = [N[chan]['Zc'][0]/N[chan]['Zb'][0],0] #error = 0 since error already taken into SF fit.
    errTmp = getRatioErrors([N[chan]['Zb'][0],0],[N[chan]['Data_inc'][0],N[chan]['Data_inc'][1]])
    N_rat[chan]['Rbj'] = [N[chan]['Zb'][0]/N[chan]['Data_inc'][0],errTmp]
    errTmp = getRatioErrors([N[chan]['Zc'][0],0],[N[chan]['Data_inc'][0],N[chan]['Data_inc'][1]])
    N_rat[chan]['Rcj'] = [N[chan]['Zc'][0]/N[chan]['Data_inc'][0],errTmp]
  if chan == 'Com':
    x = [N['Zee']['Zc'][0],N['Zmm']['Zc'][0]]
    x_err = [0,0]
    y = [N['Zee']['Zb'][0],N['Zmm']['Zb'][0]]
    y_err = [0,0]
    N_rat['Com']['Rcb'] = combineRatio(x,y,x_err,y_err,True)
    
    x = [N['Zee']['Zc'][0],N['Zmm']['Zc'][0]]
    x_err = [0,0]
    y = [N['Zee']['Data_inc'][0],N['Zmm']['Data_inc'][0]]
    y_err = [N['Zee']['Data_inc'][1],N['Zmm']['Data_inc'][1]]
    N_rat['Com']['Rcj'] = combineRatio(x,y,x_err,y_err,False)
    
    x = [N['Zee']['Zb'][0],N['Zmm']['Zb'][0]]
    x_err = [0,0]
    y = [N['Zee']['Data_inc'][0],N['Zmm']['Data_inc'][0]]
    y_err = [N['Zee']['Data_inc'][1],N['Zmm']['Data_inc'][1]]
    N_rat['Com']['Rbj'] = combineRatio(x,y,x_err,y_err,False)

print '>>>>>>>>>>>>N_rat>>>>>>>>>>>'
print N_rat

################
#Estimate acceptance ratio
################

A_rat = {'Zee':{},'Zmm':{},'Com':{}}
for chan in chans:
  A_rat[chan]['Rcb'] = [A[chan]['Zc'][0]/A[chan]['Zb'][0], getRatioErrors(A[chan]['Zc'],A[chan]['Zb'])] 
  A_rat[chan]['Rcj'] = [A[chan]['Zc'][0]/A[chan]['Zinc'][0], getRatioErrors(A[chan]['Zc'],A[chan]['Zinc'])] 
  A_rat[chan]['Rbj'] = [A[chan]['Zb'][0]/A[chan]['Zinc'][0], getRatioErrors(A[chan]['Zb'],A[chan]['Zinc'])]

print '>>>>>>>>>>A ratio>>>>>>>>>>>>'
print A_rat

###############
#Estimate tagging efficiency ratio
###############
btag_eff_rat = {'Zee':{}, 'Zmm':{}, 'Com':{}}
for chan in chans:
  btag_eff_rat[chan]['Rcb'] = [eff_b[0]/eff_c[0],getRatioErrors(eff_b,eff_c,1)]
  eff_b_tmp = [eff_b[0],eff_b[2]] #syst
  eff_c_tmp = [eff_c[0],eff_c[2]] #syst
  btag_eff_rat[chan]['Rcb'].append(getRatioErrors(eff_b_tmp,eff_c_tmp,1))
  tmp = [1,0]
  btag_eff_rat[chan]['Rcj'] = [1.0/eff_c[0],getRatioErrors(tmp,eff_c)]
  btag_eff_rat[chan]['Rcj'].append(getRatioErrors(tmp,eff_c_tmp,1))
  btag_eff_rat[chan]['Rbj'] = [1.0/eff_b[0],getRatioErrors(tmp,eff_b)]
  btag_eff_rat[chan]['Rbj'].append(getRatioErrors(tmp,eff_b_tmp,1))

print '>>>>>>>>>>>btag_eff>>>>>>>>>>>'
print btag_eff_rat

chans1 = ['Zee','Zmm']
###############
#Calculate ratio
###############
R = {'Zee':{},'Zmm':{},'Com':{}}
R_err = {'Zee':{},'Zmm':{},'Com':{}}
nm = {'Rcb':'Zc_Zb', 'Rcj': 'Zc', 'Rbj': 'Zb'}
syst_n = ['syst_sf','syst_btag_eff','syst_a']
for chan in chans1:
  for rN in ['Rcb','Rcj','Rbj']:
    R[chan][rN] = sf[chan][nm[rN]][0]*N_rat[chan][rN][0]*A_rat[chan][rN][0]*btag_eff_rat[chan][rN][0]

for chan in chans1: #no Com
  for rN in ['Rcb','Rcj','Rbj']:
    #combine stat error of sf and N_rat
    x = []
    x.append(sf[chan][nm[rN]][0])
    x.append(N_rat[chan][rN][0])
    x.append(btag_eff_rat[chan][rN][0])
    x_unc = []
    x_unc.append(sf[chan][nm[rN]][1])
    x_unc.append(N_rat[chan][rN][1]) 
    x_unc.append(btag_eff_rat[chan][rN][1]) 
    R_err[chan][rN+'_stat'] = R[chan][rN]*getMultiplicationRelative(x,x_unc)
    R_err[chan][rN+'_syst_sf'] = R[chan][rN]*sf[chan][nm[rN]][2]/sf[chan][nm[rN]][0]
    R_err[chan][rN+'_syst_a'] = R[chan][rN]*A_rat[chan][rN][1]/A_rat[chan][rN][0]
    R_err[chan][rN+'_syst_btag_eff'] = R[chan][rN]*btag_eff_rat[chan][rN][2]/btag_eff_rat[chan][rN][0]
    #further breakdown syst on SF for BLUE, since they have different correlation
    R_err[chan][rN + '_syst_sf_MC_stat_bin'] = R[chan][rN]*sf[chan][nm[rN]][3]/sf[chan][nm[rN]][0]
    R_err[chan][rN + '_syst_sf_other'] = R[chan][rN]*sqrt(sf[chan][nm[rN]][2]*sf[chan][nm[rN]][2] - sf[chan][nm[rN]][3]*sf[chan][nm[rN]][3])/sf[chan][nm[rN]][0]

#calculate total syst
for chan in chans1: #no Com
  for rN in ['Rcb','Rcj','Rbj']:
    tot = 0
    for s in syst_n:
      tot = tot + R_err[chan][rN + '_' + s]*R_err[chan][rN + '_' + s]
    R_err[chan][rN+'_syst'] = sqrt(tot)
#calculate total
for chan in chans1:
  for rN in ['Rcb','Rcj','Rbj']:
    R_err[chan][rN+'_tot'] = sqrt(R_err[chan][rN+'_stat']*R_err[chan][rN+'_stat'] + R_err[chan][rN+'_syst']*R_err[chan][rN+'_syst'])

print '>>>>>>>>>Ratio>>>>>>>>>'
print R
print R_err

#########################
#calculate ratio from BLUE
########################
weights_blue = {}
for rN in ['Rcb','Rcj','Rbj']:
#for rN in ['Rcb']:
  val_zee = R['Zee'][rN]
  val_zmm = R['Zmm'][rN]
  values = [val_zee,val_zmm]
  blue = BLUE(len(values))
  blue.AddMeasurement(values)
  
  values1 = [val_zee,val_zmm]
  blueStat = BLUE(len(values1))
  blueStat.AddMeasurement(values1)

  ### uncorrelated
  # statistical
  blue.AddUncertainty([R_err['Zee'][rN + '_stat'], R_err['Zmm'][rN + '_stat']], 0.0)
  blueStat.AddUncertainty([R_err['Zee'][rN + '_stat'], R_err['Zmm'][rN + '_stat']], 0.0)
  
  # MC statistical
  blue.AddUncertainty([R_err['Zee'][rN + '_syst_sf_MC_stat_bin'], R_err['Zmm'][rN + '_syst_sf_MC_stat_bin']], 0.0)

  ### correlated
  # other syst on SF
  blue.AddUncertainty([R_err['Zee'][rN + '_syst_sf_other'], R_err['Zmm'][rN + '_syst_sf_other']], 0.0)
  # syst on acceptance
  blue.AddUncertainty([R_err['Zee'][rN + '_syst_a'], R_err['Zmm'][rN + '_syst_a']], 0.0)
  # luminosity
  blue.AddUncertainty([R_err['Zee'][rN + '_syst_btag_eff'], R_err['Zmm'][rN + '_syst_btag_eff']], 0.0)
  #central + total
  tmp = blue.Simple()
  weights_blue['Zee_' + rN] = blue.weights.getA()[0][0]
  weights_blue['Zmm_' + rN] = blue.weights.getA()[1][0]

  print tmp
  R['Com'][rN] = tmp[0].getA()[0][0]
  R_err['Com'][rN + '_tot'] = tmp[1].getA()[0][0] 
  uncTmp = R_err['Com'][rN + '_tot']/R['Com'][rN]
  #stat
  tmpStat = blueStat.Simple() 
  uncTmp1 = tmpStat[1].getA()[0][0]/tmpStat[0].getA()[0][0] #combine value when only use stat
  R_err['Com'][rN + '_stat'] = uncTmp1*R['Com'][rN] 
  R_err['Com'][rN + '_syst']= sqrt(R_err['Com'][rN + '_tot']*R_err['Com'][rN + '_tot']-R_err['Com'][rN + '_stat']*R_err['Com'][rN + '_stat']) 
  print '>>>>>>Combined: ', R['Com'], ' ', R_err['Com']

print '>>>>>>Combine>>>>>>>>'
print R['Com']
print R_err['Com']
print weights_blue

#############
#Make latex tables
#############
fLatex = open(fOutName,'w')
fLatex.write('\\begin{table}\n')
fLatex.write('\\begin{tabular}{ccc}\n')
fLatex.write('\\hline\n')
fLatex.write('\\hline\n')
l = ['', 'Electron', 'Muon']
st = makeLatexLine(l)
fLatex.write(st)

fLatex.write('\\hline\n')

for fl in ['Zc','Zb','Zc_Zb']:
  st = ''
  if fl == 'Zc': st = 'SF$_c$'
  if fl == 'Zb': st = 'SF$_b$' 
  if fl == 'Zc_Zb': st = 'SF$_c$/SF$_b$' 
  for chan in chans1:
    st = st + ' & ' + format('%.2f' % sf[chan][fl][0]) + ' $\pm$ ' + format('%.2f (stat.)' % sf[chan][fl][1]) + ' $\pm$ ' + format('%.2f (syst.)' % sf[chan][fl][2])
    print fl, ' ', sf[chan][fl]
  st = st + '\\\\\n' 
  print st
  fLatex.write(st)

fLatex.write('\\hline\n')

st = 'N$_c^{MC}$'
for chan in chans1:
  if chan == 'Com':
    st = st + ' & - '
  else:
    st = st +  ' & ' + format('%.0f' % N[chan]['Zc'][0]) + ' $\pm$ ' + format('%.0f' % N[chan]['Zc'][1])
st = st + '\\\\\n'
fLatex.write(st)

st = 'N$_b^{MC}$'
for chan in chans1:
  if chan == 'Com':
    st = st + ' & - '
  else:
    st = st +  ' & ' + format('%.0f' % N[chan]['Zb'][0]) + ' $\pm$ ' + format('%.0f' % N[chan]['Zb'][1])
st = st + '\\\\\n'
fLatex.write(st)

st = 'N$_{inc}$'
for chan in chans1:
  if chan == 'Com':
    st = st + ' & - '
  else:
    st = st +  ' & ' + format('%.0f' % N[chan]['Data_inc'][0]) + ' $\pm$ ' + format('%.0f' % N[chan]['Data_inc'][1])
st = st + '\\\\\n'
fLatex.write(st)

for fl in ['Rcb','Rcj','Rbj']:
  st = ''
  if fl == 'Rcb': st = 'N$_c^{MC}$/N$_b^{MC}$'
  if fl == 'Rcj': st = 'N$_c^{MC}$/N$_{inc}$'
  if fl == 'Rbj': st = 'N$_b^{MC}$/N$_{inc}$'

  for chan in chans1:
    st = st + ' & ' + format('%.4f' % N_rat[chan][fl][0]) + ' $\pm$ ' + format('%.4f' % N_rat[chan][fl][1])
  st = st + '\\\\\n'
  fLatex.write(st)

fLatex.write('\\hline\n')

st = '$\epsilon^{c}_{tag}$'
for chan in chans1:
  st = st +  ' & ' + format('%.3f' % eff_c[0]) + ' $\pm$ ' + format('%.3f (stat.)' % eff_c[1]) + ' $\pm$ ' + format('%.3f (syst.)' % eff_c[2])
st = st + '\\\\\n'
fLatex.write(st)

st = '$\epsilon^{b}_{tag}$'
for chan in chans1:
  st = st +  ' & ' + format('%.3f' % eff_b[0]) + ' $\pm$ ' + format('%.3f' % eff_b[1]) + ' $\pm$ ' + format('%.3f (syst.)' % eff_b[2])
st = st + '\\\\\n'
fLatex.write(st)

st = '$\epsilon^{b}_{tag}/\epsilon^{c}_{tag}$'
for chan in chans1:
  st = st +  ' & ' + format('%.3f' % btag_eff_rat[chan]['Rcb'][0]) + ' $\pm$ ' + format('%.3f' % btag_eff_rat[chan]['Rcb'][1])
st = st + '\\\\\n'
fLatex.write(st)

fLatex.write('\\hline\n')

st = 'A$_c$'
for chan in chans1:
  st = st +  ' & ' + format('%.2f' % A[chan]['Zc'][0]) + ' $\pm$ ' + format('%.2f' % A[chan]['Zc'][1])
st = st + '\\\\\n'
fLatex.write(st)

st = 'A$_b$'
for chan in chans1:
  st = st +  ' & ' + format('%.2f' % A[chan]['Zb'][0]) + ' $\pm$ ' + format('%.2f' % A[chan]['Zb'][1])
st = st + '\\\\\n'
fLatex.write(st)

st = 'A$_{inc}$'
for chan in chans1:
  st = st +  ' & ' + format('%.2f' % A[chan]['Zinc'][0]) + ' $\pm$ ' + format('%.2f' % A[chan]['Zinc'][1])
st = st + '\\\\\n'
fLatex.write(st)

for fl in ['Rcb','Rcj','Rbj']:
  st = ''
  if fl == 'Rcb': st = 'A$_c$/A$_b$'
  if fl == 'Rcj': st = 'A$_c$/A$_{inc}$'
  if fl == 'Rbj': st = 'A$_b$/A$_{inc}$'

  for chan in chans1:
    st = st + ' & ' + format('%.2f' % A_rat[chan][fl][0]) + ' $\pm$ ' + format('%.2f' % A_rat[chan][fl][1])
  st = st + '\\\\\n'
  fLatex.write(st)

fLatex.write('\\hline\n')
fLatex.write('\\hline\n')
fLatex.write('\\end{tabular}\n')
fLatex.write('\\end{table}\n')
fLatex.write('\n')

#############
#Ratio table
#############
fLatex.write('\\begin{table}\n')
fLatex.write('\\resizebox{\columnwidth}{!}{\n')
fLatex.write('\\begin{tabular}{cccc}\n')
fLatex.write('\\hline\n')
fLatex.write('\\hline\n')
l = ['', 'Electron', 'Muon','Combined']
st = makeLatexLine(l)
fLatex.write(st)
fLatex.write('\\hline\n')
for fl in ['Rcb','Rcj','Rbj']:
  st = ''
  if fl == 'Rcb': st = 'R(c/b)'
  if fl == 'Rcj': st = 'R(c/j)'
  if fl == 'Rbj': st = 'R(b/j)'

  for chan in chans:
    st = st + ' & ' + format('%.3f' % R[chan][fl]) + '$\pm$' + format('%.3f (stat.)' % R_err[chan][fl + '_stat']) + ' $\pm$ ' + format('%.3f (syst.)' % R_err[chan][fl + '_syst'])
  st = st + '\\\\\n'
  fLatex.write(st)

fLatex.write('\\hline\n')
fLatex.write('\\hline\n')
fLatex.write('\\end{tabular}\n')
fLatex.write('}')
fLatex.write('\\end{table}\n')
fLatex.write('\n')

###########
#Systematic table
###########

err_names = [['syst_sf','SF syst. (\%)'],['syst_btag_eff','b-tagging (\%)'],['syst_a','Acceptance (\%)'],['syst','Total syst. (\%)'], ['stat','Stat. (\%)'], ['tot', 'Total (\%)']]
fLatex.write('\n')
fLatex.write('\\begin{table}\n')
fLatex.write('  \centering\n')
fLatex.write('  \\begin{tabular}{l|ccc|ccc}\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
l = ['', '\multicolumn{3}{c|}{Electron} ', '\multicolumn{3}{c}{Muon} ']
st = makeLatexLine(l)
fLatex.write(st)
l = ['','R(c/b)','R(c/j)','R(b/j)','R(c/b)','R(c/j)','R(b/j)']
st = makeLatexLine(l)
fLatex.write(st)
fLatex.write('  \\hline\n')
for err in err_names:
  if err[0] == 'syst':
    fLatex.write('\\hline\n')
  l = [err[1]]
  for chan in chans1:
    for fl in ['Rcb','Rcj','Rbj']:
      l.append(str(100*R_err[chan][fl + '_' + err[0]]/R[chan][fl]))
  st = makeLatexLine(l,'%.1f')
  fLatex.write(st)

fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\end{tabular}\n')
fLatex.write('\\end{table}\n')

###########
#Input for combine
###########
err_names = [['syst_sf_MC_stat_bin','SF MC temp. stat. (\%)'], ['syst_sf_other','SF other syst. (\%)'], ['syst_btag_eff','b-tagging (\%)'],['syst_a','Acceptance (\%)'], ['stat','Statistical (\%)']]
fLatex.write('\n')
fLatex.write('\\begin{table}\n')
fLatex.write('  \centering\n')
fLatex.write('  \\begin{tabular}{l|ccc|ccc|c}\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
l = ['', '\multicolumn{3}{c|}{Electron}', '\multicolumn{3}{c|}{Muon}',' ']
st = makeLatexLine(l)
fLatex.write(st)
l = ['','R(c/b)','R(c/j)','R(b/j)','R(c/b)','R(c/j)','R(b/j)', 'Correlation']
st = makeLatexLine(l)
fLatex.write(st)
fLatex.write('  \\hline\n')
st = makeLatexLine(l)
fLatex.write(st)
for err in err_names:
  l = [err[1]]
  for chan in chans1:
    for fl in ['Rcb','Rcj','Rbj']:
      l.append(str(100*R_err[chan][fl + '_' + err[0]]/R[chan][fl]))
  st = makeLatexLine(l,'%.1f')
  st = st.replace('\\\\','').replace('\n','')
  if err[0] == 'stat' or err[0] == 'syst_sf_MC_stat_bin':
    st = st + ' & 0 \\\\\n'
  else: 
    st = st + ' & 1 \\\\\n'

  fLatex.write(st)

fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\end{tabular}\n')
fLatex.write('\\end{table}\n')

###########
#Weights from combine
###########
fLatex.write('\n')
fLatex.write('\\begin{table}\n')
fLatex.write('  \centering\n')
fLatex.write('  \\begin{tabular}{l|ccc}\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
l = ['', 'Electron', 'Muon']
st = makeLatexLine(l)
fLatex.write(st)
fLatex.write('\\hline\n')
for fl in ['Rcb','Rcj','Rbj']:
  l = []
  if fl == 'Rcb': l.append('R(c/b)')
  if fl == 'Rcj': l.append('R(c/j)')
  if fl == 'Rbj': l.append('R(b/j)')
  l.append(100*weights_blue['Zee_' + fl])
  l.append(100*weights_blue['Zmm_' + fl])
  st = makeLatexLine(l,'%.1f')
  fLatex.write(st)
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\end{tabular}\n')
fLatex.write('\\end{table}\n')



