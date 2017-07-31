from ROOT import *
from myutils.util_funcs import *

gStyle.SetOptStat(0)
gROOT.SetBatch(False)

vtxMassType = 'vtxMass' #vtxMass, incVtxMass, vtxMassCorr_IVF

fInNamePre = 'Test/template_' + vtxMassType + '_CSVM_ttsemi_allData_V25'

fLatex = open(fInNamePre + '_latex.txt','w')

#plot_prenames = ['Jet_vtxMass']
#data_groups = {'DataCD':['DataC', 'DataD']}

lumi = 1 
f = TFile.Open(fInNamePre + '.root','read')
plots = []
plots.append((f.Get('AllJet_vtxMass_Data_semilep_px')).Clone('Data1D'))
plots.append((f.Get('AllJet_vtxMass_bkgr_semilep_px')).Clone('Bkgr'))
plots.append((f.Get('ljet_vtxMass_TT_semilep_px')).Clone('ljet_TT'))
plots.append((f.Get('bjet_vtxMass_TT_semilep_px')).Clone('bjet_TT'))
plots.append((f.Get('cjet_vtxMass_TT_semilep_px')).Clone('cjet_TT'))
  
for i in range(1, len(plots)):
  plots[i].Scale(lumi)
    
print '=======Event yeild (MC scaled to ', lumi, ') ========='
for i in range(0, len(plots)):
  print plots[i].GetName(), ': ', plots[i].Integral(0, 10000)

#############################
#stack plot
#############################

nRebin = 5
xAxisRange = [0,10]
if vtxMassType == 'incVtxMass' or vtxMassType == 'vtxMass':
  nRebin = 2
  xAxisRange = [0,6]

print 'I am going to rebin: ', nRebin
for i in range(0, len(plots)):
  plots[i].Rebin(nRebin)
    
c = makeStackPlot(plots, ['Data', 'DY+VV bkgr', 'l-jet (TT)', 'b-jet (TT)', 'c-jet (TT)'], 'c_' + vtxMassType + '_CSVM_ttsemi_allData_V25', 'Test/', 'Jet M_{SV} [GeV]', xAxisRange, 'MC unc. (stat.)', False)

nums = {'bjet_TT':[],'cjet_TT':[],'ljet_TT':[],'Bkgr':[], 'MC':[], 'DMCratio':[], 'Data1D':[]}
getNumbers(nums,plots)
for k,vs in nums.items():
  vs_st = []
  for v in vs:
    if k != 'DMCratio': vs_st.append(format('%.0f' % v))
    else: vs_st.append(format('%.3f' % v))
  nums[k] = vs_st

#############################
#write the table
#############################

fLatex.write('\n')
fLatex.write('\\begin{table}\n')
fLatex.write('  \\caption{}\n')
fLatex.write('  \\label{tab:}\n')
fLatex.write('  \centering\n')
fLatex.write('  \\begin{tabular}{l|c}\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
l = ['', 'Numer of events']
st = makeLatexLine(l)
fLatex.write(st)
fLatex.write('  \\hline\n')
l = ['  c-jet ($t\\bar{t}$)', str(nums['cjet_TT'][0]) + ' $\pm$ ' + str(nums['cjet_TT'][1]) + ' (' + str(nums['cjet_TT'][2]) + ' $\pm$ ' + str(nums['cjet_TT'][3]) + ')']
fLatex.write(makeLatexLine(l))
l = ['  b-jet ($t\\bar{t}$)', str(nums['bjet_TT'][0]) + ' $\pm$ ' + str(nums['bjet_TT'][1]) + ' (' + str(nums['bjet_TT'][2]) + ' $\pm$ ' + str(nums['bjet_TT'][3]) + ')']
fLatex.write(makeLatexLine(l))
l = ['  l-jet ($t\\bar{t}$)', str(nums['ljet_TT'][0]) + ' $\pm$ ' + str(nums['ljet_TT'][1]) + ' (' + str(nums['ljet_TT'][2]) + ' $\pm$ ' + str(nums['ljet_TT'][3]) + ')']
fLatex.write(makeLatexLine(l))
l = ['  DY + VV', str(nums['Bkgr'][0]) + ' $\pm$ ' + str(nums['Bkgr'][1]) + ' (' + str(nums['Bkgr'][2]) + ' $\pm$ ' + str(nums['Bkgr'][3]) + ')']
fLatex.write(makeLatexLine(l))
fLatex.write('  \\hline\n')
l = ['  Total MC', str(nums['MC'][0]) + ' $\pm$ ' + str(nums['MC'][1]) + ' (' + str(nums['MC'][2]) + ' $\pm$ ' + str(nums['MC'][3]) + ')']
fLatex.write(makeLatexLine(l))
l = ['  Data (semileptonic $t\\bar{t}$)', str(nums['Data1D'][0]) + ' $\pm$ ' + str(nums['Data1D'][1]) + ' (' + str(nums['Data1D'][2]) + ' $\pm$ ' + str(nums['Data1D'][3]) + ')']
fLatex.write(makeLatexLine(l))
l = ['  Data/MC', str(nums['DMCratio'][0]) + ' $\pm$ ' + str(nums['DMCratio'][1]) + ' (' + str(nums['DMCratio'][2]) + ' $\pm$ ' + str(nums['DMCratio'][3]) + ')']
fLatex.write(makeLatexLine(l))
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
fLatex.write('\\end{tabular}\n')
fLatex.write('\\end{table}\n')

##############################
#plot ratio plots
##############################

plotPrefix = 'cJet_' + vtxMassType + '_CSVM_ttsemi_allData_V25'
h = []
h.append((f.Get('AllJet_vtxMass_Data_bkgrSub_TTsub_semilep_px')).Clone('Data_bkgrSub_TTsub'))
h.append((f.Get('cjet_vtxMass_TT_semilep_px')).Clone('cjet_TT_1'))

h[0].Rebin(nRebin)
h[1].Rebin(nRebin)
shapeComparePlot(h, ['Data (Semileptonic sample)', 'c-jet (TT)'], plotPrefix, 'Jet M_{SV} [GeV]', xAxisRange)

#newBinning = [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.5,5,6,7,8,10]
#newBinning = [0,0.5,1,1.5,2,2.5,3,4,5,10]
#shapeComparePlot(h, ['Data (Semileptonic TT)', 'c-jet (TT)'], plotPrefix, 'Jet M_{SV} [GeV]', [0,10], False, False, newBinning)

#############################
#plot ratio plot to compare Z+jets and TT 
#############################

plotPrefix = 'cJet_' + vtxMassType + '_CSVM_TT_DY_allData_V25'
h = []
h.append((f.Get('cjet_vtxMass_TT_semilep_px')).Clone('cjet_TT_2'))
h.append((f.Get('cjet_vtxMass_DY')).Clone('cjet_DY_1'))

h[0].Rebin(nRebin)
h[1].Rebin(nRebin)
shapeComparePlot(h, ['c-jet (TT)', 'c-jet (DY)'], plotPrefix, 'Jet M_{SV} [GeV]', xAxisRange)

############################
#plot ratio plots to compare Z+jets and not reweighted data
############################

plotPrefix = 'cJet_' + vtxMassType + '_CSVM_ttsemi_DY_allData_V25'
h = []
h.append((f.Get('AllJet_vtxMass_Data_bkgrSub_TTsub_semilep_px')).Clone('Data_bkgrSub_TTsub_1'))
h.append((f.Get('cjet_vtxMass_DY')).Clone('cjet_DY'))

h[0].Rebin(nRebin)
h[1].Rebin(nRebin)
shapeComparePlot(h, ['Data (p_{T} reweighted)', 'c-jet (DY)'], plotPrefix, 'Jet M_{SV} [GeV]', xAxisRange)


############################
#plot ratio plots to compare Z+jets and reweighted data
############################

plotPrefix = 'cJet_' + vtxMassType + '_CSVM_ttsemi_ptReweight_allData_V25'
h = []
h.append((f.Get('AllJet_vtxMass_Data_bkgrSub_TTsub_semilep_ptReweight_px')).Clone('Data_bkgrSub_TTsub_ptReweight'))
h.append((f.Get('cjet_vtxMass_DY')).Clone('cjet_DY'))

h[0].Rebin(nRebin)
h[1].Rebin(nRebin)
shapeComparePlot(h, ['Data (p_{T} reweighted)', 'c-jet (DY)'], plotPrefix, 'Jet M_{SV} [GeV]', xAxisRange)

#newBinning = [0,0.5,1,1.5,2,2.5,3,4,5,10]
#shapeComparePlot(h, ['Data (p_{T} reweighted)', 'c-jet (DY)'], plotPrefix, 'Jet M_{SV} [GeV]', [0,10], False, False, newBinning)

