#include"TFile.h"
#include"TH2F.h"
#include"TLorentzVector.h"


TFile* ECALmap = 0;
TH2F* NewUnder = 0;
TH2F* NewOver = 0;
TH2F* NewUnderQCD = 0;
TH2F* NewOverQCD = 0;

bool Under(float phi, float eta) {
  if (ECALmap == 0) {
    ECALmap = TFile::Open("/scratch/sdonato/VHbbRun2/V21/CMSSW_7_1_5/src/Xbb/python/plot.root", "r");
  }

  if (NewUnder == 0) {
    NewUnder = (TH2F*) ECALmap->Get("NewUnder");
  }

  int bin_ = NewUnder->FindBin(phi, eta);
  int bit = NewUnder->GetBinContent(bin_);

  return (bit > 0) ? true : false;
}

bool Over(float phi, float eta) {
  if (ECALmap == 0) {
    ECALmap = TFile::Open("/scratch/sdonato/VHbbRun2/V21/CMSSW_7_1_5/src/Xbb/python/plot.root", "r");
  }

  if (NewOver == 0) {
    NewOver = (TH2F*) ECALmap->Get("NewOver");
  }

  int bin_ = NewOver->FindBin(phi, eta);
  int bit = NewOver->GetBinContent(bin_);

  return (bit > 0) ? true : false;
}

bool UnderQCD(float phi, float eta) {
  if (ECALmap == 0) {
    ECALmap = TFile::Open("/scratch/sdonato/VHbbRun2/V21/CMSSW_7_1_5/src/Xbb/python/plot.root", "r");
  }

  if (NewUnderQCD == 0) {
    NewUnderQCD = (TH2F*) ECALmap->Get("NewUnderQCD");
  }

  int bin_ = NewUnderQCD->FindBin(phi, eta);
  int bit = NewUnderQCD->GetBinContent(bin_);

  return (bit > 0) ? true : false;
}

bool OverQCD(float phi, float eta) {
  if (ECALmap == 0) {
    ECALmap = TFile::Open("/scratch/sdonato/VHbbRun2/V21/CMSSW_7_1_5/src/Xbb/python/plot.root", "r");
  }

  if (NewOverQCD == 0) {
    NewOverQCD = (TH2F*) ECALmap->Get("NewOverQCD");
  }

  int bin_ = NewOverQCD->FindBin(phi,eta);
  int bit = NewOverQCD->GetBinContent(bin_);

  return (bit > 0) ? true : false;
}

bool UnderAll(float phi, float eta) {
  return (Under(phi, eta) || UnderQCD(phi, eta));
}

bool OverAll(float phi, float eta) {
  return (Over(phi, eta) || OverQCD(phi, eta));
}

float M(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2){
  TLorentzVector p1, p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, 0);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, 0);
  return (p1+p2).M();
}

float M(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
  TLorentzVector p1, p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
  return (p1+p2).M();
}

float M3(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3) {
  TLorentzVector p1, p2, p3;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
  p3.SetPtEtaPhiM(pt3, eta3, phi3, m3);
  return (p1+p2+p3).M();
}

float Pt(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {
  TLorentzVector p1, p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, 0);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, 0);
  return (p1+p2).Pt();
}

float Phi(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {
  TLorentzVector p1, p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, 0);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, 0);
  return (p1+p2).Phi();
}

