import os,sys
import ConfigParser
from myutils import BetterConfigParser 
import ROOT

ROOT.gSystem.Load('../interface/VHbbNameSpace_h.so')

ROOT.gROOT.SetBatch(True)

#####################
#some settings
#####################
cfg = BetterConfigParser()
cfg.read('Config/config.ini')
regions = ['zjet','zHFjet'] #regions to make plots zjet = z + >=1jet, zHFjet = z + >=1 HF jet
#regions = ['zHFjet'] #regions to make plots zjet = z + >=1jet, zHFjet = z + >=1 HF jet
fOut = ROOT.TFile.Open('Test/histos_zjet_zHFjet_with_sf_CSVM.root','recreate')
######################
#get data and MC files and sample cat cuts
######################
chains = {'Data_Zee':ROOT.TChain('tree'), 'Data_Zmm':ROOT.TChain('tree'), 'DY':ROOT.TChain('tree'), 'TT':ROOT.TChain('tree'), 'WW':ROOT.TChain('tree'), 'WZ':ROOT.TChain('tree'), 'ZZ':ROOT.TChain('tree')}

path = cfg.get('General','path')
print '>>>>>>>Sample used: >>>>>>>>>>>>'
for i in range(23):
  sample = 'Sample_' + str(i)
  name = cfg.get(sample,'name')
  fileName = path + '/' + cfg.get(sample,'file')
  print fileName 
  if name == 'DoubleEG':
    chains['Data_Zee'].Add(fileName)
    print 'Current chain events: ', chains['Data_Zee'].GetEntries()
  if name == 'DoubleMuon':
    chains['Data_Zmm'].Add(fileName)
    print 'Current chain events: ', chains['Data_Zmm'].GetEntries()
  if name == 'DY':
    chains['DY'].Add(fileName)
    print 'Current chain events: ', chains['DY'].GetEntries()
  if name == 'TT':
    chains['TT'].Add(fileName)
    print 'Current chain events: ', chains['TT'].GetEntries()
  if name == 'WW':
    chains['WW'].Add(fileName)
    print 'Current chain events: ', chains['WW'].GetEntries()
  if name == 'WZ':
    chains['WZ'].Add(fileName)
    print 'Current chain events: ', chains['WZ'].GetEntries()
  if name == 'ZZ':
    chains['ZZ'].Add(fileName)
    print 'Current chain events: ', chains['ZZ'].GetEntries()

#####################
#Print number of entries
#####################
print '########################'
print '# Number of entry'
print '########################'
for k,v in chains.items():
  print k, ' ', v.GetEntries()

#TEMP
#sys.exit()
#loop over channel
for chan in ['Zee','Zmm']:
  hDir = fOut.mkdir(chan)
  hDir.cd()
  #loop over regions
  for r in regions:
    rName = r + '_' + chan
    plots = cfg.get('Plots',r + '_plot').split(',')
    cutCommon = cfg.get('Cuts', rName) #get cut for a region
    sf = cfg.get('SFs','sf_'+rName)
    cat_cuts = cfg.get('Cuts',r+'_cat').split(',')
    hDir1 = hDir.mkdir(r) #make 
    hDir1.cd()
    #loop over samples
    for k,v in chains.items():
      #Temp
      #if 'DY' not in k: continue
      if 'Data' in k and chan not in k: continue #for data only loop over corresponding data set 
      #add json requirement to cut
      cut = cutCommon
      if 'Data' in k: cut = cfg.get('Cuts','jsonCut') + '&&' + cut
      #if 'Data' not in k: cut = '(' + cut + ')' + '*(' + sf + ')'
      print '>>>>>>>>>>>>>>>>>>>>>>>'
      print 'Chan, region, sample: ', chan, ' ', r, ' ', k
      print 'Entries: ', v.GetEntries()
      print 'Cut: ', cut
      for plot in plots:
        print '>>>>', plot
        var = cfg.get(plot,'var')
        xAxisRange = cfg.get(plot,'range').split(',')
        xAxisRange1 = [] #make conversion to accommodate special case of Pi()
        for i in range(len(xAxisRange)):
          if 'TMath::Pi' in xAxisRange[i] and '-' in xAxisRange[i]:
            xAxisRange1.append(-ROOT.TMath.Pi())
          elif 'TMath::Pi' in xAxisRange[i] and '-' not in xAxisRange[i]:
            xAxisRange1.append(ROOT.TMath.Pi())
          else:
            xAxisRange1.append(xAxisRange[i])
        plotName = plot + '_' + k.replace('_' + chan,'') + '_' + rName
        print plotName
        h = ROOT.TH1D(plotName,'',int(xAxisRange1[0]),float(xAxisRange1[1]),float(xAxisRange1[2]))
        h.Sumw2()
        cutAndScale = '(' + cut + ')'
        if 'Data' not in k: cutAndScale = cutAndScale + '*(' + sf + ')'
        print 'Cut and scale: ', cutAndScale
        v.Draw(var + '>>' + plotName,cutAndScale)
        h.Write()
        #make plot for sample categories, for example DY+bjet, DY+cjet ...
        if 'DY' in k:
          for cat_cut in cat_cuts:
            plotName = plot + '_' + k.replace('_' + chan,'') + '_' + cat_cut.split(':')[0] + '_' + rName
            print plotName
            h = ROOT.TH1D(plotName,'',int(xAxisRange1[0]),float(xAxisRange1[1]),float(xAxisRange1[2]))
            h.Sumw2()
            cutTmp = cut+'&&('+cat_cut.split(':')[1]+')'
            cutAndScale = '(' + cutTmp + ')' + '*(' + sf + ')'
            print '>>>>>>Cut and scale tmp: >>>>>>>>', cutAndScale
            v.Draw(var + '>>' + plotName,cutAndScale)
            h.Write()

    hDir.cd()
  fOut.cd()
