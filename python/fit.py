from ROOT import *

lumi = 7000
nFitTemplate = 3
templateChoice = 0 #0 = DYinc, 1 = DY1Jet, 2 = combine
useDDbjet = False 
nRebin = 2
xMax_fitRange = 6

#fIn = TFile.Open('template_CSVM.root','read')
fIn = TFile.Open('template_CSVT_allWeight.root','read')

fTmp = TFile('dd_bjet_temp.root', 'read')

gSystem.Load('../interface/VHbbNameSpace_h.so')

h = {'Data':[], 'DY':[], 'DY1Jet':[], 'DYcombine':[], 'TT':[], 'WW':[], 'WZ':[], 'ZZ':[]}

def getHisto(fIn, name):
  hist = fIn.Get(name)
  return hist

chan = 'Zmm'

for sample,val in h.items():
  
  bNameTmp = 'bjet_vtxMass_tmp_' + sample + '_' + chan
  cNameTmp = 'cjet_vtxMass_tmp_' + sample + '_' + chan
  lNameTmp = 'ljet_vtxMass_tmp_' + sample + '_' + chan
  bName  = 'bjet_vtxMass_' + sample + '_' + chan
  cName  = 'cjet_vtxMass_' + sample + '_' + chan
  lName  = 'ljet_vtxMass_' + sample + '_' + chan

  if sample == 'Data':
    h[sample].append(getHisto(fIn, 'AllJet_vtxMass_Data_' + chan))
    continue
  
  if sample != 'DYcombine':
    h[sample].append(getHisto(fIn, bName))
    h[sample].append(getHisto(fIn, cName))
    h[sample].append(getHisto(fIn, lName))

  
  if (sample == 'DY' or sample == 'DY1Jet' or sample == 'DYcombine'):
    h[sample].append(getHisto(fIn, bNameTmp))
    h[sample].append(getHisto(fIn, cNameTmp))
    h[sample].append(getHisto(fIn, lNameTmp))
     
#===Rebin histograms===
for sample,val in h.items():
  print sample, val
  for i in range(0, len(val)):
    val[i].Rebin(nRebin)


#===Get background====
hBkgr = []
for hTmp in h['TT']:
  name = hTmp.GetName()
  name = name.split('_')[0]
  hBkgr.append(hTmp.Clone('Bkgr_' + name))

for sample,hTmps in h.items():
  if sample == 'WW' or sample == 'WZ' or sample == 'ZZ':
    for i in range(0, len(hTmps)):
      name = hTmps[i].GetName()
      name = name.split('_')[0]
      #print '========================'
      #print name
      for j in range(0, len(hBkgr)):
        name1 = hBkgr[j].GetName()
        name1 = name1.split('_')[1]
        #print name1
        if name == name1:
          hBkgr[j].Add(hTmps[i])
          break

hTmp = hBkgr[0].Clone('Bkgr_all')
for i in range(1, len(hBkgr)):
  hTmp.Add(hBkgr[i])
  if nFitTemplate == 2:
    hTmp.Add(h['DY'][2])
hBkgr.append(hTmp)

#===Subtract background from data====
hData = h['Data'][0].Clone('Data_bkgrSubtracted')
hData.Add(hBkgr[len(hBkgr)-1], -1)

#===Get template=====================
hTemps = []
sample = 'DY'
nBegin = 3
if templateChoice == 1: sample = 'DY1Jet'
if templateChoice == 2: 
  sample = 'DYcombine'
  nBegin = 0

for i in range(nBegin, len(h[sample])):
  name = 'Clone_' + h[sample][i].GetName()
  hTemps.append(h[sample][i].Clone(name))

if useDDbjet:
 hTmp = fTmp.Get('jet_vtxMass_ptweight_emu')
 bNameTmp = 'bjet_vtxMass_tmp_DY' 
 hTmp1 = getHisto(fIn, bNameTmp)
 hTmp1.Reset()
 hTmp1.SetName('Data-drive_bjet_tmp')
 nBin = hTmp.GetNbinsX()
 for i in range(1, nBin+1):
   iBin = hTmp1.FindFixBin(hTmp.GetBinLowEdge(i))
   hTmp1.SetBinContent(iBin,hTmp.GetBinContent(i))
   hTmp1.SetBinError(iBin,hTmp.GetBinError(i))
 hTemps[0] = hTmp1


print '==============================='
for hTmp in hTemps:
  print hTmp.GetName(), hTmp.Integral()

#===Setting the first bin to zero====
hData.SetBinContent(1,0)
hData.SetBinError(1,0)
for i in range(0, len(hBkgr)):
  hBkgr[i].SetBinContent(1,0)
  hBkgr[i].SetBinError(1,0)

for i in range(0, len(hTemps)):
  hTemps[i].SetBinContent(1,0)
  hTemps[i].SetBinError(1,0)

#===plot=============================
c1 = TCanvas("c1")
c1.Divide(2,3)
c1.cd(1)
h['Data'][0].Draw()
c1.cd(2)
hBkgr[len(hBkgr)-1].Draw()
print '================='
print hBkgr[len(hBkgr)-1].Integral()
c1.cd(3)
hData.Draw()
c1.cd(4)
hTemps[0].Draw()
c1.cd(5)
hTemps[1].Draw()
c1.cd(6)
hTemps[2].Draw()
c1.Print('c1.png')

#======now do the fit================
temp_objs = TObjArray()
for i in range(0,3-(3-nFitTemplate)): temp_objs.Add(hTemps[i])
fit = TFractionFitter(hData,temp_objs)

for i in range(1, hData.GetNbinsX()+1):
  if hData.GetBinContent(i) < 0:
    fit.ExcludeBin(i)

nLastBin = hData.FindFixBin(xMax_fitRange)
fit.SetRangeX(2, nLastBin)
fit.Constrain(0,0.0,1.0)
fit.Constrain(1,0.0,1.0)
if nFitTemplate > 2:
  fit.Constrain(2,0.0,1.0)
fit.Fit()
fit.ErrorAnalysis(0.54)

frac = [0,0,0]
fracErr = [0,0,0]
for i in range(0,3):
 fracTmp = Double(0)
 fracErrTmp = Double(0)
 fit.GetResult(i, fracTmp, fracErrTmp)
 frac[i]= fracTmp
 fracErr[i] = fracErrTmp

chi2 = fit.GetChisquare() ;
prob = fit.GetProb() ;
npdf = fit.GetNDF() ;

#========print results===========
#===Print data and MC prediction===
nBkgr = [0, 0, 0, 0, 0] #DY_b, DY_c, DY_l, TT, VV
nBkgr_err = [0, 0, 0, 0, 0]
nData = 0

def getError(h, sBin):
  err = 0
  for i in range(sBin, h.GetNbinsX()+1):
    err = err + h.GetBinError(i)*h.GetBinError(i)
  return TMath.Sqrt(err)

for sample,hTmps in h.items():
  nTot = 0
  nErr = 0
  if sample != 'Data':
    for i in range(0,3):
      nTot = nTot + hTmps[i].Integral(2, 1000)
      nErr = nErr + getError(hTmps[i], 2)*getError(hTmps[i], 2)
    nErr = TMath.Sqrt(nErr)
  if sample == 'Data':
    nData = hTmps[0].Integral(2, 1000)
    #print 'Data: ',  nData     
  elif sample == 'DY':
    nBkgr[0] = hTmps[0].Integral(2, 1000)
    nBkgr_err[0] = getError(hTmps[0],2)
    nBkgr[1] = hTmps[1].Integral(2, 1000)
    nBkgr_err[1] = getError(hTmps[1],2)
    nBkgr[2] = hTmps[2].Integral(2, 1000)
    nBkgr_err[2] = getError(hTmps[2],2)
    #print 'DY bjet: ', nBkgr[0]
    #print 'DY cjet: ', nBkgr[1]
    #print 'DY ljet: ', nBkgr[2]
  else:
    #print sample + ': ', nTot
    if sample == 'TT':
      nBkgr[3] = nTot
      nBkgr_err[3] = nErr
    if sample == 'WW' or sample == 'WZ' or sample == 'ZZ':
      nBkgr[4] = nBkgr[4] + nTot
      nBkgr_err[4] = nBkgr_err[4] + nErr*nErr

nBkgr_err[4] = TMath.Sqrt(nBkgr_err[4])

print '======================='
print 'Data: ', nData
print 'DY + bjet: ', nBkgr[0], ' +/- ', nBkgr_err[0]
print 'DY + cjet: ', nBkgr[1], ' +/- ', nBkgr_err[1]
print 'DY + ljet: ', nBkgr[2], ' +/- ', nBkgr_err[2]
print 'TT: ', nBkgr[3], ' +/- ', nBkgr_err[3]
print 'VV: ', nBkgr[4], ' +/- ', nBkgr_err[4]
totBkgr = nBkgr[0] + nBkgr[1] + nBkgr[2] + nBkgr[3] + nBkgr[4]
print 'Total MC: ', totBkgr, ' +/- ', TMath.Sqrt(nBkgr_err[0]*nBkgr_err[0] + nBkgr_err[1]*nBkgr_err[1] + nBkgr_err[2]*nBkgr_err[2] + nBkgr_err[3]*nBkgr_err[3] + nBkgr_err[4]*nBkgr_err[4])
print 'Data/MC: ', nData/totBkgr

def ratioErr(a, b, a_err, b_err): #a/b
  return a/b*TMath.Sqrt(a_err*a_err/(a*a) + b_err*b_err/(b*b))


print '======================='
print 'MC-predicted fractions: '
print 'Two components (b and c): '
a = nBkgr[0]
a_err = nBkgr_err[0]
b = (nBkgr[0] + nBkgr[1])
b_err = TMath.Sqrt(nBkgr_err[0]*nBkgr_err[0] + nBkgr_err[1]*nBkgr_err[1])
print 'f_b: ', a/b, ' +/- ', ratioErr(a,b, a_err, b_err)

a = nBkgr[1]
a_err = nBkgr_err[1]
print 'f_c: ', a/b, ' +/- ', ratioErr(a, b, a_err, b_err)

print 'Three components (b, c and light): '
a = nBkgr[0]
a_err = nBkgr_err[0]
b = (nBkgr[0] + nBkgr[1] + nBkgr[2])
b_err = TMath.Sqrt(nBkgr_err[0]*nBkgr_err[0] + nBkgr_err[1]*nBkgr_err[1] + nBkgr_err[2]*nBkgr_err[2])
print 'f_b: ', a/b, ' +/- ', ratioErr(a,b, a_err, b_err)

a = nBkgr[1]
a_err = nBkgr_err[1]
print 'f_c: ', a/b, ' +/- ', ratioErr(a,b, a_err, b_err)

a = nBkgr[2]
a_err = nBkgr_err[2]
print 'f_l: ', a/b, ' +/- ', ratioErr(a,b, a_err, b_err)


print '========================'
print 'Fit resuts: '
print 'f_b: ', frac[0], ' +/- ', fracErr[0]
print 'f_c: ', frac[1], ' +/- ', fracErr[1]
print 'f_l: ', frac[2], ' +/- ', fracErr[2]

print 'Chi2: ', chi2
print 'Chi2/ndof: ', chi2/npdf
print 'Prob: ', prob

hFitPlot = fit.GetPlot() ;
hRat = hData.Clone("DataOverFit") ;
hRat.Divide(hFitPlot) ;
for i in range(1, hRat.GetNbinsX() + 1):
  if (hData.GetBinContent(i) != 0):
    hRat.SetBinError(i, hRat.GetBinContent(i)*hData.GetBinError(i)/hData.GetBinContent(i))
  else: hRat.SetBinError(i,0)

#=====scale the histogram==========
nData = hData.Integral(2, 1000)
for i in range(0,3):
  hTemps[i].Scale(frac[i]*nData/hTemps[i].Integral(2,1000))

#=====plot stack plot====
maxBin = hData.GetMaximumBin()
c2 = TCanvas("c2", "", 0, 22, 640, 700)
pad1 = TPad("pad1", "pad1", 0,0.3,1,1) 
pad1.SetTopMargin(0.05)
pad1.SetBottomMargin(0.019)
pad1.Draw()
pad1.cd()
  
st = THStack("st", "") ;
for i in range(0,3-(3-nFitTemplate)):
  hTemps[i].SetFillColor(i+2)
  st.Add(hTemps[i])
st.Draw('hist')
st.SetMaximum(1.2*hData.GetBinContent(maxBin))
st.GetYaxis().SetTitle("Events")
st.GetXaxis().SetRangeUser(0,6)
hData.SetMarkerStyle(20)
hData.Draw("same pe") 
hFitPlot.SetLineColor(46) 
hFitPlot.SetLineStyle(2) 
hFitPlot.SetLineWidth(2) 
hFitPlot.Draw("same hist") 
  
leg = TLegend(0.7, 0.6, 0.92, 0.92)   
leg.AddEntry(hData,"Data","lp") 
leg.AddEntry(hFitPlot,"Fit","lp") 
leg.AddEntry(hTemps[0],"b-jet","f") 
leg.AddEntry(hTemps[1],"c-jet","f") 
if nFitTemplate != 2:
  leg.AddEntry(hTemps[2],"l-jet","f") 
leg.SetShadowColor(0) 
leg.SetLineColor(0) 
leg.Draw() 

pad1.Update()

c2.cd()

pad2 = TPad("pad2","pad2",0,0,1,0.28)
pad2.SetTopMargin(0.006)
pad2.SetBottomMargin(0.32) 
pad2.Draw()

pad2.cd()

gStyle.SetOptStat(0) 
hRat.Draw("0 ep") 
hRat.SetTitle("") 
hRat.SetMarkerStyle(20)
hRat.GetXaxis().SetTitleSize(0.11) 
hRat.GetXaxis().SetLabelSize(0.115)
hRat.GetXaxis().SetRangeUser(0,6) 
hRat.GetYaxis().SetTitleSize(0.2) 
hRat.GetYaxis().SetLabelSize(0.115) 
hRat.GetYaxis().SetRangeUser(0,2.1) 
hRat.GetYaxis().SetNdivisions(505) 
hRat.GetXaxis().SetTitleFont(42) 
hRat.GetXaxis().SetTitleSize(0.1166667)
hRat.GetYaxis().SetTitle("DATA/Fit") 
#hRat.GetYaxis().SetTitleSize(0.1166667)
#hRat.GetYaxis().SetTitleOffset(0.5)
line = TLine(0,1,6,1) 
line.SetLineWidth(3) 
line.Draw("same") 
line.SetLineStyle(2)  

pad2.Update()

c2.Print('c2.png')


del fit


'''if __name__ == '__main__':
   rep = ''
   while not rep in [ 'q', 'Q' ]:
     rep = raw_input( 'enter "q" to quit: ' )
     if 1 < len(rep):
       rep = rep[0]
'''
