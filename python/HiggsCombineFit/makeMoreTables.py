import os,sys
sys.path.append('/uscms_data/d3/duong/CMSSW/CMSSW_7_6_5/src/ZplusC/python/')
from myutils.util_funcs import *

#make table of SF results

fLatex = open('tables_vtxMass_pt20_binStatThresh_0_useDDbjet_useDDljet_customBinning.tex','w')

#chans = ['Zee','Zmm','Com']
chans = ['Zee','Zmm']
preName = 'test_vtxMass_pt20_allData_binStatThresh_0_useDDbjet_useDDljet_customBinning'

#texts = {'Zee':{},'Zmm':{},'Com':{}}
texts = {'Zee':{},'Zmm':{}}
for chan in chans:
  lines = open(preName + '_' + chan + '/' + preName + '_' + chan + '_newErrorEstimation_latex.txt','r').readlines()
  for l in lines:
    lsp = l.split('&')
    print lsp
    if 'SF$_c$' in l and 'stat' in l:
      texts[chan]['SFc'] = lsp[1].replace('\\\\','').replace('\n','') 
    if 'SF$_b$' in l and 'stat' in l:
      texts[chan]['SFb'] = lsp[1].replace('\\\\','').replace('\n','')
    if 'SF$_l$' in l and 'stat' in l:
      texts[chan]['SFl'] = lsp[1].replace('\\\\','').replace('\n','')

print texts

fLatex.write('\n')
fLatex.write('\\begin{table}\n')
fLatex.write('  \centering\n')
fLatex.write('  \\begin{tabular}{lccc}\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
l = ['', 'Electron ', 'Muon ']
st = makeLatexLine(l)
fLatex.write(st)
fLatex.write('  \\hline\n')

labels = ['SFc','SFb','SFl']
label1s = ['SF$_c$','SF$_b$','SF$_l$']

for i in range(len(labels)):
  l = [label1s[i], texts['Zee'][labels[i]], texts['Zmm'][labels[i]]]
  st = makeLatexLine(l)
  fLatex.write(st)

fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\end{tabular}\n')
fLatex.write('\\end{table}\n')

#####################
#Uncertainty table
#####################

#texts = {'Zee':{},'Zmm':{},'Com':{}}
texts = {'Zee':{},'Zmm':{}}
tags = ['JEC &', 'JER &', 'MC temp. stat. &', 'gcc splitting &', 'gbb splitting &', 'PU reweight &', 'Bkgr. sub. &']
label1s = ['JEC (\%)', 'JER (\%)','MC temp. stat. (\%)', 'gcc splitting (\%)', 'gbb splitting (\%)', 'Pileup (\%)', 'Bkgr. subtraction (\%)']

for chan in chans:
  for tag in tags:
    texts[chan][tag] = {}

print texts

for chan in chans:
  lines = open(preName + '_' + chan + '/' + preName + '_' + chan + '_newErrorEstimation_latex.txt','r').readlines()
  for l in lines:
    lsp = l.split('&')
    for tag in tags:
      if tag in l:
        texts[chan][tag]['SFc'] = lsp[1].replace('\\\\','').replace('\n','') 
        texts[chan][tag]['SFb'] = lsp[2].replace('\\\\','').replace('\n','')
        texts[chan][tag]['SFl'] = lsp[3].replace('\\\\','').replace('\n','')

print texts

fLatex.write('\n')
fLatex.write('\\begin{table}\n')
fLatex.write('  \centering\n')
fLatex.write('  \\begin{tabular}{l|ccc|ccc}\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
l = ['', '\multicolumn{3}{c|}{Electron} ', '\multicolumn{3}{c}{Muon} ']
st = makeLatexLine(l)
fLatex.write(st)
l = ['','SF$_c$','SF$_b$','SF$_l$','SF$_c$','SF$_b$','SF$_l$']
st = makeLatexLine(l)
fLatex.write(st)

fLatex.write('  \\hline\n')

for i in range(len(label1s)):
  tag = tags[i]
  l = [label1s[i], texts['Zee'][tag]['SFc'], texts['Zee'][tag]['SFb'], texts['Zee'][tag]['SFl'], texts['Zmm'][tag]['SFc'], texts['Zmm'][tag]['SFb'], texts['Zmm'][tag]['SFl']]
  st = makeLatexLine(l)
  fLatex.write(st)

fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\end{tabular}\n')
fLatex.write('\\end{table}\n')

#################################
#make covariance matrix
#################################
cov = {'Zmm_Zc':[],'Zee_Zc':[]}
lines = open('../'+preName.replace('test_','summary_analyzeFit_').replace('allData_','') + '.txt').readlines()
for k,v in cov.items():
  for l in lines:
    if k.split('_')[0] in l and k.split('_')[1] + '_cov_' in l:
      cov[k].append([l.split()[2],l.split()[3],l.split()[4].replace('\n','')])
print cov['Zmm_Zc']
print cov['Zee_Zc']
sf_names = ['SF_c','SF_b','SF_l']
for k,v in cov.items():
  if 'Zmm' in k:
    fLatex.write('\n$Cov_{ij}^{muon} =\\bordermatrix{&SF_c&SF_b&SF_l\\cr')
  if 'Zee' in k:
    fLatex.write('\n$Cov_{ij}^{ele} =\\bordermatrix{&SF_c&SF_b&SF_l\\cr')
  for i in range(3):
    st = '\n ' + sf_names[i] + '&'
    for j in range(3):
      if j != 2:
        st = st + v[i][j] + '&'
      else:
        st = st + v[i][j] + '\\cr'
    if i != 2:
      fLatex.write(st)
    else:
      fLatex.write(st + '}$')

