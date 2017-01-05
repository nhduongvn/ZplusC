import ROOT

ROOT.gStyle.SetOptStat(0)
f = []
f.append(ROOT.TFile.Open('template_CSVT_noWeight.root','read'))
f.append(ROOT.TFile.Open('template_CSVT_allWeight.root','read'))

template_name = ['bjet_vtxMass_tmp_DY_Zmm', 'cjet_vtxMass_tmp_DY_Zmm', 'ljet_vtxMass_tmp_DY_Zmm']

for temp in template_name:
  h = []
  h.append(f[0].Get(temp))
  h.append(f[1].Get(temp))

  for i in range(2):
    h[i].Rebin(2)
    h[i].Scale(1./h[i].Integral())
  hR = h[0].Clone('Ratio')
  hR.Divide(h[1])

  c = ROOT.TCanvas('c_' + temp,'', 600, 600)
  c.SetFillStyle(4000)
  c.SetFrameFillStyle(1000)
  c.SetFrameFillColor(0)

  oben = ROOT.TPad('oben','oben',0,0.3 ,1.0,1.0)
  oben.SetBottomMargin(0)
  oben.Draw()

  unten = ROOT.TPad('unten','unten',0,0.0,1.0,0.3)
  unten.SetTopMargin(0.)
  unten.SetBottomMargin(0.35)
  unten.Draw()

  oben.cd()
  l = ROOT.TLegend(0.5, 0.6,0.8,0.9)
  l.SetLineWidth(2)
  l.SetBorderSize(0)
  l.SetFillColor(0)
  l.SetTextFont(62)
  l.SetTextSize(0.035)

  xRange = 7

  h[1].SetMarkerStyle(21)
  h[1].SetMarkerSize(1.5)
  h[1].SetMarkerColor(4)
  h[1].SetLineColor(4)
  #h[1].SetFillColor(4)
  #h[1].Draw('hist')
  h[1].Draw('ep')
  #h[1].GetYaxis().SetRangeUser(0,0.2)
  h[1].GetYaxis().SetTitle("Events/bin")
  h[1].GetXaxis().SetRangeUser(0,xRange)


  h[0].Draw('ep same')
  h[0].SetMarkerStyle(20)
  h[0].SetMarkerSize(1.5)

  l.AddEntry(h[0], 'No weights', 'lp')
  l.AddEntry(h[1], 'With weights', 'lp')
  l.Draw()

  c.Update()

  unten.cd()
  ROOT.gPad.SetTicks(1,1)
  hR.Draw('E')
  hR.SetMarkerStyle(20)
  hR.SetMarkerSize(1.5)
  hR.GetYaxis().SetLabelSize(0.1)
  hR.GetXaxis().SetLabelSize(0.1)
  hR.GetXaxis().SetTitleSize(0.1)
  hR.GetYaxis().SetTitleSize(0.1)
  hR.GetYaxis().SetTitle('Ratio')
  hR.GetYaxis().SetTitleOffset(0.5)
  hR.GetXaxis().SetTitle('Jet vtxMass (GeV)')
  hR.GetXaxis().SetRangeUser(0,xRange)

  hR.GetYaxis().SetRangeUser(0.5,1.5)
  line = ROOT.TLine(0,1,xRange,1)
  line.Draw('same')

  c.Update()

  c.Print(temp + '.png')
