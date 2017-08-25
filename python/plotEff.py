from ROOT import *
from myutils.util_funcs import *
gROOT.SetBatch(True)

def plotEff1(hs, h_label, canName, legHeader, hs_stat = []):
    c = TCanvas(canName)
    c.SetLogy()
    c.SetLogx()
    c.SetLeftMargin(0.0988539)
    c.SetRightMargin(0.0475759) 
    c.SetBottomMargin(0.157428)
    c.SetTopMargin(0.0199557)
    
    hs[0].Draw()
    hs[0].SetStats(0)
    hs[0].GetXaxis().SetTitle('p_{T} [GeV]')
    hs[0].GetXaxis().SetTitleOffset(1.4)
    hs[0].GetXaxis().SetTitleSize(0.05)
    hs[0].GetXaxis().SetLabelSize(0.05)
    hs[0].GetYaxis().SetRangeUser(0.001,70)
    hs[0].GetYaxis().SetTitle('Efficiency')
    hs[0].GetYaxis().SetTitleSize(0.05)
    hs[0].GetYaxis().SetLabelSize(0.05)
    hs[0].SetLineColor(633) #read
    hs[0].SetLineWidth(2)
    hs[0].SetMarkerColor(633)
    hs[0].SetMarkerStyle(20)
    hs[0].SetMarkerSize(1.5)
    if len(hs_stat) != 0:
      hs_stat[0].Draw("E2 same")
      hs_stat[0].SetFillColor(633)
      hs_stat[0].SetFillStyle(3004)
    
    hs[1].Draw('same')
    hs[1].SetLineColor(417) #green
    hs[1].SetLineWidth(2)
    hs[1].SetMarkerColor(417)
    hs[1].SetMarkerStyle(21)
    hs[1].SetMarkerSize(1.5)
    if len(hs_stat) != 0:
      hs_stat[1].Draw("E2 same")
      hs_stat[1].SetFillColor(417)
      hs_stat[1].SetFillStyle(3004)


    hs[2].Draw('same')
    hs[2].SetLineColor(601) #blue
    hs[2].SetLineWidth(2)
    hs[2].SetMarkerColor(601)
    hs[2].SetMarkerStyle(22)
    hs[2].SetMarkerSize(1.5)
    if len(hs_stat) != 0:
      hs_stat[2].Draw("E2 same")
      hs_stat[2].SetFillColor(601)
      hs_stat[2].SetFillStyle(3004)
    
    if len(hs_stat) != 0: 
      hDummy = hs_stat[0].Clone('Dummy')
      hDummy.SetFillColor(1)
    leg = TLegend(0.59,0.70,0.90,0.96)
    leg.SetHeader(legHeader)
    leg.SetLineColor(0)
    leg.SetFillColor(0)
    leg.SetTextSize(0.05)
    leg.AddEntry(hs[0],h_label[0],'lp')
    leg.AddEntry(hs[1],h_label[1],'lp')
    leg.AddEntry(hs[2],h_label[2],'lp')
    if len(hs_stat) != 0: leg.AddEntry(hDummy,'Stat. Unc.','f')
    leg.Draw()
    #pause()
    c.Print(canName + '.pdf')
    c.Print(canName + '.C')
    c.Print(canName + '.eps')
    c.Print(canName + '.png')

def getEfficiency(hist):
  eff = {}
  for iB in range(1,hist.GetNbinsX() + 1):
    binStr = str(hist.GetBinLowEdge(iB)) + '_' + str(hist.GetBinLowEdge(iB+1))
    eff[binStr] = [hist.GetBinContent(iB), hist.GetBinError(iB)]
  return eff
        
#########################################################
#Main program
#########################################################

f = TFile.Open('Test/eff_13TeV_CSVM.root')
fLatex = open(f.GetName().replace('.root','_latex.txt'),'w')


legHeaders = {'tag':'Efficiency (CSVM)', 'tag_svt':'Efficiency (CSVM + SVT)'}
effs = {'jet_tag':{'bjet':{},'cjet':{},'ljet':{}},'jet_tag_svt':{'bjet':{},'cjet':{},'ljet':{}},'event_tag':{'bjet':{},'cjet':{},'ljet':{}},'event_tag_svt':{'bjet':{},'cjet':{},'ljet':{}}}
effs_int = {'jet_tag':{'bjet':[],'cjet':[],'ljet':[]},'jet_tag_svt':{'bjet':[],'cjet':[],'ljet':[]},'event_tag':{'bjet':[],'cjet':[],'ljet':[]},'event_tag_svt':{'bjet':[],'cjet':[],'ljet':[]}}

ptBinLabel = []
getPtBinLabel = True

for selType in ['tag','tag_svt']:
    
    ##############
    #per jet efficiency
    ##############
    hs = [f.Get('hB_pt_eff_' + selType), f.Get('hC_pt_eff_' + selType), f.Get('hL_pt_eff_' + selType)]
    canName = "Test/Eff_plots/eff_pt_" + selType
    #TEMP turn off making plots
    plotEff1(hs,['b-jet','c-jet','l-jet'],canName,legHeaders[selType])
    #get ptBin label
    if getPtBinLabel:
      for iBin in range(1,hs[0].GetNbinsX()+1):
        ptBinLabel.append(str(hs[0].GetBinLowEdge(iBin)) + '_' + str(hs[0].GetBinLowEdge(iBin+1)))
      getPtBinLabel = False

    #get efficiency to make latex table
    effs['jet_' + selType]['bjet'] = getEfficiency(hs[0])
    effs['jet_' + selType]['cjet'] = getEfficiency(hs[1])
    effs['jet_' + selType]['ljet'] = getEfficiency(hs[2])

    #####################
    #event efficiency
    #####################
    hs = [f.Get('hB_pt_evt_eff_' + selType), f.Get('hC_pt_evt_eff_' + selType), f.Get('hL_pt_evt_eff_' + selType)]
    #statistical histogram
    hs_unc_stat = [hs[0].Clone(hs[0].GetName()+'_statUnc'), hs[1].Clone(hs[1].GetName()+'_statUnc'), hs[2].Clone(hs[2].GetName()+'_statUnc')]
    #get efficiency with statistical error to make latex table
    effs['event_' + selType]['bjet'] = getEfficiency(hs_unc_stat[0])
    effs['event_' + selType]['cjet'] = getEfficiency(hs_unc_stat[1])
    effs['event_' + selType]['ljet'] = getEfficiency(hs_unc_stat[2])
    #systematic uncertainty histogram
    hs_sf_up = [f.Get('hB_pt_evt_SFu_eff_' + selType), f.Get('hC_pt_evt_SFu_eff_' + selType), f.Get('hL_pt_evt_SFu_eff_' + selType)]
    hs_sf_down = [f.Get('hB_pt_evt_SFd_eff_' + selType), f.Get('hC_pt_evt_SFd_eff_' + selType), f.Get('hL_pt_evt_SFd_eff_' + selType)]
    for i in range(3): #loop over b,c,l
        for iBin in range(1, hs_unc_stat[i].GetNbinsX() + 1):
            systUnc = max(abs(hs_sf_up[i].GetBinContent(iBin)-hs[i].GetBinContent(iBin)),abs(hs_sf_down[i].GetBinContent(iBin)-hs[i].GetBinContent(iBin)))
            statUnc = hs_unc_stat[i].GetBinError(iBin)
            totUnc = sqrt(systUnc*systUnc+statUnc*statUnc)
            #set total uncertainty to hs
            hs[i].SetBinError(iBin,totUnc)
            #save systematic uncertainty to effs
            binStr = str(hs[i].GetBinLowEdge(iBin)) + '_' + str(hs[i].GetBinLowEdge(iBin+1))
            if i == 0:
                effs['event_' + selType]['bjet'][binStr].append(systUnc)
            if i == 1:
                effs['event_' + selType]['cjet'][binStr].append(systUnc)
            if i == 2:
                effs['event_' + selType]['ljet'][binStr].append(systUnc)
             
    
    canName = "Test/eff_pt_evt_" + selType
    #TEMP turn off making plots
    plotEff1(hs,['Z+b-jets','Z+c-jets','Z+l-jets'],canName,legHeaders[selType],hs_unc_stat)

    
#############
#Integral efficiency
#############
for effType in ['jet','event']: 
  hs = []
  hs_syst_up = []
  hs_syst_down = []
  if effType == 'jet':
    hs = [f.Get('hB_pt_eff_int'), f.Get('hC_pt_eff_int'), f.Get('hL_pt_eff_int')]
  if effType == 'event':
    hs = [f.Get('hB_pt_evt_eff_int'), f.Get('hC_pt_evt_eff_int'), f.Get('hL_pt_evt_eff_int')]
    hs_syst_up = [f.Get('hB_pt_evt_SFu_eff_int'), f.Get('hC_pt_evt_SFu_eff_int'), f.Get('hL_pt_evt_SFu_eff_int')]
    hs_syst_down = [f.Get('hB_pt_evt_SFd_eff_int'), f.Get('hC_pt_evt_SFd_eff_int'), f.Get('hL_pt_evt_SFd_eff_int')]
  
  for i in range(3):
    if i == 0: jetType = 'bjet'
    if i == 1: jetType = 'cjet'
    if i == 2: jetType = 'ljet'
    effs_int[effType + '_tag'][jetType] = [hs[i].GetBinContent(1), hs[i].GetBinError(1)]
    effs_int[effType + '_tag_svt'][jetType] = [hs[i].GetBinContent(2), hs[i].GetBinError(2)]
    if effType == 'event':
        syst_err = max(abs(hs_syst_up[i].GetBinContent(1)-hs[i].GetBinContent(1)), abs(hs_syst_down[i].GetBinContent(1) - hs[i].GetBinContent(1)))
        effs_int[effType + '_tag'][jetType].append(syst_err)
        syst_err = max(abs(hs_syst_up[i].GetBinContent(2)-hs[i].GetBinContent(2)), abs(hs_syst_down[i].GetBinContent(2) - hs[i].GetBinContent(2)))
        effs_int[effType + '_tag_svt'][jetType].append(syst_err)

#########################################
#make latex table for pt binned
#########################################
print effs
print ptBinLabel

#convert effs to percentage
for k,v in effs.items():
    for k1,v1 in v.items():
        for k2,v2 in v1.items():
            for i in range(len(v2)): v2[i] = 100*v2[i]
                

for effType in ['jet','event']:
  fLatex.write('\n')
  fLatex.write('\\begin{table}\n')
  fLatex.write('  \\caption{}\n')
  fLatex.write('  \\label{tab:}\n')
  fLatex.write('  \\centering\n')
  fLatex.write('  \\begin{tabular}{lcccccc}\n')
  fLatex.write('  \\hline\n')
  fLatex.write('  \\hline\n')
  l = ['', '\multicolumn{3}{c}{CSVM}','\multicolumn{3}{c}{CSVM + SVT}']
  st = makeLatexLine(l)
  fLatex.write(st)
  if effType == 'jet':
    l = ['', 'b-jets','c-jets','l-jets','b-jets','c-jets','l-jets']
  if effType == 'event':
    l = ['', 'Z+b-jets','Z+c-jets','Z+l-jets','Z+b-jets','Z+c-jets','Z+l-jets']
  st = makeLatexLine(l)
  fLatex.write(st)
  fLatex.write('  \\hline\n')
  for label in ptBinLabel:
      numFormat = '%.1f'
      l = []
      if effType == 'jet':
          l = [label.replace('.0','').replace('_','-')]
          for jetType in ['bjet','cjet','ljet']:
            if jetType == 'ljet': numFormat = '%.2f'
            else: numFormat = '%.1f'
            l.append(format(numFormat % effs[effType + '_tag'][jetType][label][0]) + ' $\pm$ ' +  format(numFormat % effs[effType + '_tag'][jetType][label][1]))
          for jetType in ['bjet','cjet','ljet']:
            if jetType == 'ljet': numFormat = '%.2f'
            else: numFormat = '%.1f'
            l.append(format(numFormat % effs[effType + '_tag_svt'][jetType][label][0]) + ' $\pm$ ' + format(numFormat % effs[effType + '_tag_svt'][jetType][label][1]))
      if effType == 'event':
          l = [label.replace('.0','').replace('_','-')]
          for jetType in ['bjet','cjet','ljet']:
            if jetType == 'ljet': numFormat = '%.2f'
            else: numFormat = '%.1f'
            l.append(format(numFormat % effs[effType + '_tag'][jetType][label][0]) + ' $\pm$ ' +  format(numFormat % effs[effType + '_tag'][jetType][label][1]) + ' $\pm$ ' + format(numFormat % effs[effType + '_tag'][jetType][label][2]))
          for jetType in ['bjet','cjet','ljet']:
            if jetType == 'ljet': numFormat = '%.2f'
            else: numFormat = '%.1f'
            l.append(format(numFormat % effs[effType + '_tag_svt'][jetType][label][0]) + ' $\pm$ ' + format(numFormat % effs[effType + '_tag_svt'][jetType][label][1]) + ' $\pm$ ' + format(numFormat % effs[effType + '_tag_svt'][jetType][label][2]))
      st = makeLatexLine(l)
      fLatex.write(st)
  fLatex.write('  \\hline\n')
  fLatex.write('  \\hline\n')
  fLatex.write('\\end{tabular}\n')
  fLatex.write('\\end{table}\n')


#########################################
#make latex table for integral efficiency 
#########################################
#convert efficiency to percentage
for k,v in effs_int.items():
    for k1,v1 in v.items():
        for i in range(len(v1)):
            v1[i] = 100*v1[i]

fLatex.write('\n')
fLatex.write('\\begin{table}\n')
fLatex.write('  \\caption{}\n')
fLatex.write('  \\label{tab:}\n')
fLatex.write('  \\centering\n')
fLatex.write('  \\begin{tabular}{lccc}\n')
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
l = ['', 'Z+b-jets','Z+c-jets','Z+l-jets']
st = makeLatexLine(l)
fLatex.write(st)
fLatex.write('  \\hline\n')
for selType in ['tag','tag_svt']:
  l = []
  if selType == 'tag':
      l=['CSVM']
  if selType == 'tag_svt':
      l=['CSVM + SVT']
  for jetType in ['bjet','cjet','ljet']:
      numFormat = '%.1f'
      if jetType != 'bjet': numFormat = '%.2f'
      l.append(format(numFormat % effs_int['event_'+selType][jetType][0]) + ' $\pm$ ' + format(numFormat % effs_int['event_'+selType][jetType][1]) + ' $\pm$ ' + format(numFormat % effs_int['event_'+selType][jetType][2]))
  st = makeLatexLine(l)
  fLatex.write(st)
fLatex.write('  \\hline\n')
fLatex.write('  \\hline\n')
fLatex.write('\\end{tabular}\n')
fLatex.write('\\end{table}\n')

