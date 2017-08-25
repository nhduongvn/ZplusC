import os,sys
import ConfigParser
from myutils import BetterConfigParser,util_funcs
import ROOT
import copy

ROOT.gSystem.Load('../interface/VHbbNameSpace_h.so')

ROOT.gROOT.SetBatch(True)

def scalePlotToLumi(cfg,sam_name,his):
  #loop over samples to pick up the right process to scale, for example DY, TT ...
  for i in range(23):
    sam_name_tmp = cfg.get('Sample_'+str(i),'name')
    if sam_name == sam_name_tmp:
      xSecTmp = cfg.get('Sample_'+str(i),'xSec').split('*')
      xSec = 1.
      for x in xSecTmp:
          xSec = xSec*float(x)
      sf = util_funcs.scaleToLumi(cfg.get('General','path') + '/' + cfg.get('Sample_'+str(i),'file'), xSec,float(cfg.get('General','lumi')))
      his.Scale(sf)

def printEventCount(plots,nums,label_nums, label_translate, fLatex):
    #nums = {'bjet':[],'cjet':[],'ljet':[],'Bkgr':[], 'MC':[], 'DMCratio':[], 'Data1D':[]}
    util_funcs.getNumbers(nums,plots)
    for k,vs in nums.items():
      vs_st = []
      for v in vs:
        if k != 'DMCratio': vs_st.append(format('%.0f' % v))
        else: vs_st.append(format('%.3f' % v))
      nums[k] = vs_st

    fLatex.write(plots[0].GetName())
    fLatex.write('\n')
    fLatex.write('\\begin{table}\n')
    fLatex.write('  \\caption{}\n')
    fLatex.write('  \\label{tab:}\n')
    fLatex.write('  \\centering\n')
    fLatex.write('  \\begin{tabular}{l|c}\n')
    fLatex.write('  \\hline\n')
    fLatex.write('  \\hline\n')
    l = ['', 'Number of events']
    st = util_funcs.makeLatexLine(l)
    fLatex.write(st)
    fLatex.write('  \\hline\n')
    for label in label_nums:
      if 'MC' == label:
          fLatex.write('\\hline\n')
      l = [label_translate[label], str(nums[label][2]) + ' $\pm$ ' + str(nums[label][3])]
      fLatex.write(util_funcs.makeLatexLine(l))
    fLatex.write('  \\hline\n')
    fLatex.write('  \\hline\n')
    fLatex.write('\\end{tabular}\n')
    fLatex.write('\\end{table}\n')

def printEventCountTwoChannels(nums,label_nums,label_translate,fLatex):
    #print nums
    fLatex.write('\n')
    fLatex.write('\\begin{table}\n')
    fLatex.write('  \\caption{}\n')
    fLatex.write('  \\label{tab:}\n')
    fLatex.write('  \\centering\n')
    fLatex.write('  \\begin{tabular}{lcc}\n')
    fLatex.write('  \\hline\n')
    fLatex.write('  \\hline\n')
    l = ['', 'Electron','Muon']
    st = util_funcs.makeLatexLine(l)
    fLatex.write(st)
    fLatex.write('  \\hline\n')
    for label in label_nums:
      if 'MC' == label:
          fLatex.write('\\hline\n')
      l = [label_translate[label], str(nums['Zee'][label][2]) + ' $\pm$ ' + str(nums['Zee'][label][3]), str(nums['Zmm'][label][2]) + ' $\pm$ ' + str(nums['Zmm'][label][3])]
      fLatex.write(util_funcs.makeLatexLine(l))
    fLatex.write('  \\hline\n')
    fLatex.write('  \\hline\n')
    fLatex.write('\\end{tabular}\n')
    fLatex.write('\\end{table}\n')


#####################
#some settings
#####################
cfg = BetterConfigParser()
cfg.read('Config/config.ini')
regions = ['zjet','zHFjet'] #regions to make plots zjet = z + >=1jet, zHFjet = z + >=1 HF jet
#regions = ['zHFjet'] #regions to make plots zjet = z + >=1jet, zHFjet = z + >=1 HF jet
affix = 'zjet_zHFjet_with_sf_CSVM'  
removeAllPlot = False

fIn = ROOT.TFile.Open('Test/histos_' + affix + '.root','read')
plotFolder = 'Test/ControlPlots_' + affix
os.system('mkdir ' + plotFolder + '/')
if removeAllPlot:
    os.system('rm ' + plotFolder + '/*')

fLatexName = 'Test/' + affix + '_latex.txt'
os.system('rm ' + fLatexName)
fLatex = open(fLatexName,'w')

print '>>>>>>>>>>>>>>>>>'
print 'Input data: ', fIn.GetName()
print 'Output folder for plots: ', plotFolder
print 'Latex file containing event count: ', fLatexName

samples = ['Data','DY','TT','WW','WZ','ZZ']
#loop over channel
event_counts = {}
for chan in ['Zee','Zmm']:
  #loop over regions
  for r in regions:
    rName = r + '_' + chan
    plots = cfg.get('Plots',r + '_plot').split(',')
    cutCommon = cfg.get('Cuts', rName) #get cut for a region
    #loop over plots
    printNumber = True
    for plot in plots:
      #loop over samples 
      histos = {}
      for sam in samples:
        if sam != 'DY':
          plotName = plot + '_' + sam + '_' + r + '_' + chan
          print plotName
          histos[sam] = fIn.Get(chan + '/' + r + '/' + plotName)
          #scale plot to lumi
          if sam != 'Data': scalePlotToLumi(cfg,sam,histos[sam])
        else:
          for jetType in ['bjet','cjet','ljet']:
            plotName = plot + '_' + sam + '_' + jetType + '_' + r + '_' + chan
            print plotName
            histos[sam+'_'+jetType] = fIn.Get(chan + '/' + r + '/' + plotName)
            #scale plot to lumi
            scalePlotToLumi(cfg,sam,histos[sam+'_'+jetType])
      
      #rebin the histogram
      for k,v in histos.items():
          v.Rebin(int(cfg.get(plot,'rebin')))

      #adding VV together
      hVV = histos['WW'].Clone('VV')
      for sam in ['WZ','ZZ']:
          hVV.Add(histos[sam])
      histos_forPlotting = [histos['Data'],hVV,histos['TT'],histos['DY_ljet'],histos['DY_bjet'],histos['DY_cjet']]
      histos_forEventCount = [histos['Data'].Clone('Data1D'),hVV.Clone('VV'),histos['TT'].Clone('TT'),histos['DY_ljet'].Clone('DY_ljet'),histos['DY_bjet'].Clone('DY_bjet'),histos['DY_cjet'].Clone('DY_cjet')] #clone plot to change the name
      histosNames = ['Data (' + chan.replace('Zee','Electron').replace('Zmm','Muon') + ')', 'VV', 't#bar{t}','Z+l-jet','Z+b-jet','Z+c-jet']
      #print the event count
      if printNumber and ('eta' in plot):
          printNumber = False
          nums = {'Data1D':[],'MC':[],'DMCratio':[],'TT':[],'VV':[],'DY_ljet':[],'DY_bjet':[],'DY_cjet':[]}
          printEventCount(histos_forEventCount,nums,['DY_bjet','DY_cjet','DY_ljet','TT','VV','MC','Data1D','DMCratio'],{'DY_bjet':'Z+b-jets','DY_cjet':'Z+c-jets','DY_ljet':'Z+l-jets','TT':'$t\\bar{t}$','VV':'VV', 'MC':'Total MC', 'Data1D':'Data', 'DMCratio':'Data/MC'}, fLatex)
          event_counts[chan+'_'+r] = nums

      cName = plot + '_' + r + '_' + chan
      xAxisTitle = cfg.get(plot,'xAxisTitle')
      xAxisRangeTmp = cfg.get(plot,'xAxisRange').split(',')
      xAxisRange = []
      for tmp in xAxisRangeTmp:
        if '-TMath' in tmp: xAxisRange.append(-ROOT.TMath.Pi())
        elif 'TMath' in tmp and '-' not in tmp: xAxisRange.append(ROOT.TMath.Pi())
        else: xAxisRange.append(float(tmp))
      print '>>>>>>>>>>>', plot, ' ', xAxisRange
      xMin = 10
      if 'njet' in plot: xMin = 1 
      util_funcs.makeStackPlot(histos_forPlotting, histosNames, cName, plotFolder, xAxisTitle, xAxisRange, 'MC unc. (stat.)', False, False, -1, [5,2,4,1,3], xMin) #the list is the order to add plots to legend
      #util_funcs.pause()
      util_funcs.makeStackPlot(histos_forPlotting, histosNames, cName+'_log', plotFolder, xAxisTitle, xAxisRange, 'MC unc. (stat.)', False, True, -1, [5,2,4,1,3], xMin)
      #util_funcs.pause()

print event_counts

for r in regions:
  fLatex.write('\n>>>>>>>>>>>Region: ' + r + '>>>>>>>>>>>>>>>>')
  event_counts_tmp = {}
  event_counts_tmp['Zee'] = event_counts['Zee_'+r]
  event_counts_tmp['Zmm'] = event_counts['Zmm_'+r]
  printEventCountTwoChannels(event_counts_tmp,['DY_bjet','DY_cjet','DY_ljet','TT','VV','MC','Data1D','DMCratio'],{'DY_bjet':'Z+b-jets','DY_cjet':'Z+c-jets','DY_ljet':'Z+l-jets','TT':'$t\\bar{t}$','VV':'VV', 'MC':'Total MC', 'Data1D':'Data', 'DMCratio':'Data/MC'}, fLatex)
