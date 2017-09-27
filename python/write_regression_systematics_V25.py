#!/usr/bin/env python
####################################################################
#This code take the jet regression and systematic shift and apply them to some quanties and write those to the original tree
####################################################################

import sys,hashlib
import os,sys,subprocess
import ROOT 
import math
import shutil
from array import array
import warnings
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )
from optparse import OptionParser
from btag_reweight import *
from time import gmtime, strftime
from multiprocessing import Pool

from myutils import BetterConfigParser, ParseInfo, TreeCache, LeptonSF, util_funcs

import numpy as np

ROOT.gROOT.SetBatch(True)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#INPUT
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

argv = sys.argv
parser = OptionParser()
parser.add_option("-S", "--samples", dest="names", default="", 
                      help="samples you want to run on")
parser.add_option("-C", "--config", dest="config", default=[], action="append",
                      help="configuration defining the plots to make")

#don't use it anymore
#parser.add_option("-f", "--filelist", dest="filelist", default="",
#                              help="list of files you want to run on")

(opts, args) = parser.parse_args(argv)
if opts.config =="":
        opts.config = "config"

#print 'opts.filelist="'+opts.filelist+'"'
#filelist=filter(None,opts.filelist.replace(' ', '').split(';'))
#print filelist
#print "len(filelist)",len(filelist),
#if len(filelist)>0:
#    print "filelist[0]:",filelist[0];
#else:
#    print ''

print opts.config
config = BetterConfigParser()
config.read(opts.config)
debug = eval(config.get('General','Debug'))
useVtypeNew = eval(config.get('General','UseVtypeNew'))
anaTag = config.get("Analysis","tag")
TrainFlag = eval(config.get('Analysis','TrainFlag'))
btagLibrary = config.get('BTagReshaping','library')
samplesinfo=config.get('Directories','samplesinfo')
channel=config.get('Configuration','channel')
print 'channel is', channel

#VHbbNameSpace=config.get('VHbbNameSpace','library')
#ROOT.gSystem.Load(VHbbNameSpace)
AngLikeBkgs=eval(config.get('AngularLike','backgrounds'))
ang_yield=eval(config.get('AngularLike','yields'))

pathIN = config.get('Directories','SYSin')
pathOUT = config.get('Directories','SYSout')
#tmpDir = os.environ["TMPDIR"]
tmpDir = config.get("Directories",'tmpDir')

print 'INput samples:\t%s'%pathIN
print 'OUTput samples:\t%s'%pathOUT


#lhe_weight_map = False if not config.has_option('LHEWeights', 'weights_per_bin') else eval(config.get('LHEWeights', 'weights_per_bin'))


run_locally = config.get('Configuration', 'run_locally')

if debug:
  run_locally = 'True'
  pathOUT = '/uscmst1b_scratch/lpc1/lpctrig/duong/test/'


namelist=opts.names.split(',')
if run_locally == 'False':
  namelist = []
  a = opts.names.split('_')
  b = ''
  for i in range(0, len(a)-1):
    if i != len(a) - 2:
      b += a[i] + '_'
    else:
      b += a[i]
  namelist.append(b)


#load info
info = ParseInfo(samplesinfo,pathIN)

#@@@@@@@@@@@@@@@@muon sf@@@@@@@@@@@@@@@@@@
wdir = config.get('Directories','vhbbpath')

muIso_SF_reader_BtoF = LeptonSF(wdir+'/python/json/V25/SFs_BCDEF_muonIso.json', 'LooseISO_LooseID_pt_eta', 'abseta_pt_ratio')
muIso_SF_reader_GtoH = LeptonSF(wdir+'/python/json/V25/SFs_GH_muonIso.json', 'LooseISO_LooseID_pt_eta', 'abseta_pt_ratio')
muID_SF_reader_BtoF = LeptonSF(wdir+'/python/json/V25/EfficienciesAndSF_BCDEF_muonID.json', 'MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta', 'abseta_pt_ratio')
muID_SF_reader_GtoH = LeptonSF(wdir+'/python/json/V25/EfficienciesAndSF_GH_muonID.json', 'MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta', 'abseta_pt_ratio')
muTrk_SF_reader_BtoF = LeptonSF(wdir+'/python/json/V25/Tracking_EfficienciesAndSF_BCDEF_muon.json', 'Graph', 'ratio_eff_eta3_dr030e030_corr')
muTrk_SF_reader_GtoH = LeptonSF(wdir+'/python/json/V25/Tracking_EfficienciesAndSF_GH_muon.json', 'Graph', 'ratio_eff_eta3_dr030e030_corr')

muTrig_Eff_reader_1 = LeptonSF(wdir+'/python/json/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.json', 'IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093', 'abseta_pt_DATA')
muTrig_Eff_reader_2 = LeptonSF(wdir+'/python/json/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.json', 'IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097', 'abseta_pt_DATA')

#@@@@@@@@@@@@@@electronn sf@@@@@@@@@@@@@@@
eleID_SF_file = ROOT.TFile.Open('Data/egammaEffiAll.txt_SF2D.root', 'read')
eleID_SF_his = eleID_SF_file.Get('EGamma_SF2D_Medium')
eleRes_SF_file = ROOT.TFile.Open('Data/V25/egammaEffi_reconstruction_EGM2D.root', 'read')
eleRes_SF_his = eleRes_SF_file.Get('EGamma_SF2D')
eleMVA80_SF_file = ROOT.TFile.Open('Data/V25/egammaEffi_MVA80_EGM2D.root', 'read')
eleMVA80_SF_his = eleMVA80_SF_file.Get('EGamma_SF2D')
eleMVA90_SF_file = ROOT.TFile.Open('Data/V25/egammaEffi_MVA90_EGM2D.root', 'read')
eleMVA90_SF_his = eleMVA90_SF_file.Get('EGamma_SF2D')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#btag SF
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ROOT.gSystem.Load('libCondFormatsBTauObjects') 
ROOT.gSystem.Load('libCondToolsBTau')
bTagCalib = ROOT.BTagCalibration('CSVv2', 'Data/CSVv2_Moriond17_B_H.csv')
v_sys = getattr(ROOT, 'vector<string>')()
v_sys.push_back('up')
v_sys.push_back('down')

# make a reader instance and load the sf data
btagSF_reader_M = ROOT.BTagCalibrationReader(
    1,              # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
    "central",      # central systematic type
    v_sys,          # vector of other sys. types
)   
 
btagSF_reader_M.load(
    bTagCalib, 
    0,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

btagSF_reader_M.load(
    bTagCalib, 
    1,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

btagSF_reader_M.load(
    bTagCalib, 
    2,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "incl"      # measurement type
)


# make a reader instance and load the sf data
btagSF_reader_T = ROOT.BTagCalibrationReader(
    2,              # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
    "central",      # central systematic type
    v_sys,          # vector of other sys. types
)    
btagSF_reader_T.load(
    bTagCalib, 
    0,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

btagSF_reader_T.load(
    bTagCalib, 
    1,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

btagSF_reader_T.load(
    bTagCalib, 
    2,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG  
    "incl"      # measurement type
)

#################################
#ctagSF
#################################
cTagCalib = ROOT.BTagCalibration('cTag', 'Data/ctagger_Moriond17_B_H.csv')

# make a reader instance and load the sf data
ctagSF_reader_M = ROOT.BTagCalibrationReader(
    1,              # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
    "central",      # central systematic type
    v_sys,          # vector of other sys. types
)   
 
ctagSF_reader_M.load(
    cTagCalib, 
    0,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "TnP"      # measurement type
)

ctagSF_reader_M.load(
    cTagCalib, 
    1,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

ctagSF_reader_M.load(
    cTagCalib, 
    2,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "incl"      # measurement type
)


# make a reader instance and load the sf data
ctagSF_reader_T = ROOT.BTagCalibrationReader(
    2,              # 0 is for loose op, 1: medium, 2: tight, 3: discr. reshaping
    "central",      # central systematic type
    v_sys,          # vector of other sys. types
)

ctagSF_reader_T.load(
    cTagCalib, 
    0,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "TnP"      # measurement type
)

ctagSF_reader_T.load(
    cTagCalib, 
    1,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG 
    "comb"      # measurement type
)

ctagSF_reader_T.load(
    cTagCalib, 
    2,          # 0 is for b flavour, 1: FLAV_C, 2: FLAV_UDSG  
    "incl"      # measurement type
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Some functions
    
def computeSF(weight_SF,weight):
    weight_SF[0] = (weight[0][0]*weight[1][0])
    weight_SF[1] = ( (weight[0][0]-weight[0][1])*(weight[1][0]-weight[1][1]) )
    weight_SF[2] = ( (weight[0][0]+weight[0][1])*(weight[1][0]+weight[1][1]) )

def computeTrigWeight(weight_trig,weight):
  if weight[0][0] != 1.0 and weight[1][0] != 1.0:
    weight_trig[0] = 1-(1-weight[0][0])*(1-weight[1][0])
    weight_trig[1] = 1-(1-(weight[0][0]+weight[0][1]))*(1-(weight[1][0]+weight[1][1]))
    weight_trig[2] = 1-(1-(weight[0][0]-weight[0][1]))*(1-(weight[1][0]-weight[1][1]))
  if weight[0][0] == 1.0:
    weight_trig[0] = weight[1][0]
    weight_trig[1] = weight[1][0] + weight[1][1]
    weight_trig[2] = weight[1][0] - weight[1][1]
  if weight[1][0] == 1.0:
    weight_trig[0] = weight[0][0]
    weight_trig[1] = weight[0][0] + weight[0][1]
    weight_trig[2] = weight[0][0] - weight[0][1]

def getSF_2Dhis(eta, pt, h2D):
  weight = []
  iBin = h2D.FindFixBin(eta, pt)
  weight.append(h2D.GetBinContent(iBin))
  weight.append(h2D.GetBinError(iBin))
  if weight[0] <= 0.0: return [1.0,0.0] #case when eta, pt out of range
  return weight
  
def findMatchGenJet(jetEta,jetPhi,ch):
  minDr = 100
  iMatch = -1
  for iGenJet in range(0,ch.nGenJet):
    genJetEta = ch.GenJet_eta[iGenJet]
    genJetPhi = ch.GenJet_phi[iGenJet]
    genJetPt = ch.GenJet_pt[iGenJet]
    dRtmp = util_funcs.deltaR(jetEta,jetPhi,genJetEta,genJetPhi)
    if dRtmp < minDr:
      minDr = dRtmp
      iMatch = iGenJet
  if minDr < 0.3:
    return iMatch
  return -1


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Main functions
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


def fillTree(inputfile,outputfile,skimCut):

    input = ROOT.TFile.Open(inputfile,'read')
    output = ROOT.TFile.Open(outputfile,'recreate')
    print ''
    print 'inputfile',inputfile
    print "Writing: ",outputfile
    print ''

#write all stuffs from input (except tree, tree will be cloned later) to output

    input.cd()
    obj = ROOT.TObject
    for key in ROOT.gDirectory.GetListOfKeys():
        input.cd()
        obj = key.ReadObj()
        #skip writing old tree
        if obj.GetName() == job.tree:
            continue
        output.cd()
        obj.Write(key.GetName())

    input.cd()
    tree = input.Get(job.tree)
    nEntries = tree.GetEntries()

#now clone the tree

    output.cd()
    newtree = tree.CloneTree(0)

    # Add training Flag
    EventForTraining = array('i',[0])
    newtree.Branch('EventForTraining',EventForTraining,'EventForTraining/I')
    EventForTraining[0]=0

#divide MC sample by half to use for the MVA training 

    TFlag=ROOT.TTreeFormula("EventForTraining","evt%2",tree)
    specialWeight = ROOT.TTreeFormula('specialWeight',str(job.specialweight), tree)


    #Special weights
    DY_specialWeight= array('f',[0])
    DY_specialWeight[0] = 1
    newtree.Branch('DY_specialWeight',DY_specialWeight,'DY_specialWeight/F')

    Jet_vtxMassCorr1 = array('f',100*[-1])
    newtree.Branch('Jet_vtxMassCorr1',Jet_vtxMassCorr1,'Jet_vtxMassCorr1[tree.nJet]/F') #=1 for Vtype == 0, muon channel, applied for electron ID and trigger (HLT?)

    emu_id_iso_weight = array('f',[0])
    emu_id_iso_weight[0] = 1
    newtree.Branch('emu_id_iso_weight',emu_id_iso_weight,'emu_id_iso_weight/F') #=1 for Vtype == 0, muon channel, applied for electron ID and trigger (HLT?)
    
    emu_trig_weight = array('f',[0])
    emu_trig_weight[0] = 1
    newtree.Branch('emu_trig_weight',emu_trig_weight,'emu_trig_weight/F') #=1 for Vtype == 0, muon channel, applied for electron ID and trigger (HLT?)
    
    is_emu = array('i',[0])
    is_emu[0] = 0
    newtree.Branch('is_emu', is_emu, 'is_emu/I')

    emu_lep_pt = array('f',2*[-1])
    emu_lep_pt[0], emu_lep_pt[1] = -1, -1
    newtree.Branch('emu_lep_pt',emu_lep_pt,'emu_lep_pt[2]/F') #=0 for muon, 1 for electron
    
    emu_lep_phi = array('f',2*[-10])
    emu_lep_phi[0], emu_lep_phi[1] = -10, -10
    newtree.Branch('emu_lep_phi',emu_lep_phi,'emu_lep_phi[2]/F') #=0 for muon, 1 for electron

    emu_lep_eta = array('f',2*[-10])
    emu_lep_eta[0], emu_lep_eta[1] = -10, -10
    newtree.Branch('emu_lep_eta',emu_lep_eta,'emu_lep_eta[2]/F') #=0 for muon, 1 for electron
    
    emu_lep_etaSc = array('f',2*[-10])
    emu_lep_etaSc[0], emu_lep_etaSc[1] = -10, -10
    newtree.Branch('emu_lep_etaSc',emu_lep_etaSc,'emu_lep_etaSc[2]/F') #=0 for muon, 1 for electron
    
    emu_lep_mass = array('f',2*[-1])
    emu_lep_mass[0], emu_lep_mass[1] = -1, -1
    newtree.Branch('emu_lep_mass',emu_lep_mass,'emu_lep_mass[2]/F') #=0 for muon, 1 for electron
    
    emu_lep_charge = array('f',2*[0])
    emu_lep_charge[0], emu_lep_charge[1] = 0, 0
    newtree.Branch('emu_lep_charge',emu_lep_charge,'emu_lep_charge[2]/F') #=0 for muon, 1 for electron
    
    emu_lep_iso = array('f',2*[-1])
    emu_lep_iso[0], emu_lep_iso[1] = -1, -1
    newtree.Branch('emu_lep_iso',emu_lep_iso,'emu_lep_iso[2]/F') #=0 for muon, 1 for electron
    
    emu_lep_pdgId = array('i',2*[-1])
    emu_lep_pdgId[0], emu_lep_pdgId[1] = -1, -1
    newtree.Branch('emu_lep_pdgId',emu_lep_pdgId,'emu_lep_pdgId[2]/I') #=0 for muon, 1 for electron

    emu_lep_looseId = array('i',2*[-1])
    emu_lep_looseId[0], emu_lep_looseId[1] = -1, -1
    newtree.Branch('emu_lep_looseId',emu_lep_looseId,'emu_lep_looseId[2]/I') #=0 for muon, 1 for electron
    emu_lep_mediumId = array('i',2*[-1])
    emu_lep_mediumId[0], emu_lep_mediumId[1] = -1, -1
    newtree.Branch('emu_lep_mediumId',emu_lep_mediumId,'emu_lep_mediumId[2]/I') #=0 for muon, 1 for electron
    emu_lep_tightId = array('i',2*[-1])
    emu_lep_tightId[0], emu_lep_tightId[1] = -1, -1
    newtree.Branch('emu_lep_tightId',emu_lep_tightId,'emu_lep_tightId[2]/I') #=0 for muon, 1 for electron
    
    emu_lep_idEle = array('i',2*[-1])
    emu_lep_idEle[0], emu_lep_idEle[1] = -1, -1
    newtree.Branch('emu_lep_idEle',emu_lep_idEle,'emu_lep_idEle[2]/I') #=0 for muon, 1 for electron

    muweight = array('f', 3*[1])
    muweight[0], muweight[1], muweight[2] = 1., 1., 1.
    newtree.Branch('muweight',muweight,'muweight[3]/F')
    
    muweight_id = array('f', 3*[1])
    muweight_id[0], muweight_id[1], muweight_id[2] = 1., 1., 1.
    newtree.Branch('muweight_id',muweight_id,'muweight_id[3]/F')
    
    muweight_iso = array('f', 3*[1])
    muweight_iso[0], muweight_iso[1], muweight_iso[2] = 1., 1., 1.
    newtree.Branch('muweight_iso',muweight_iso,'muweight_iso[3]/F')
    
    muweight_trk = array('f', 3*[1])
    muweight_trk[0], muweight_trk[1], muweight_trk[2] = 1., 1., 1.
    newtree.Branch('muweight_trk',muweight_trk,'muweight_trk[3]/F')
    
    muweight_trig = array('f', 3*[1])
    muweight_trig[0], muweight_trig[1], muweight_trig[2] = 1., 1., 1.
    newtree.Branch('muweight_trig',muweight_trig,'muweight_trig[3]/F')
   
    eleweight = array('f', 3*[1])
    eleweight[0], eleweight[1], eleweight[2] = 1., 1., 1.
    newtree.Branch('eleweight',eleweight,'eleweight[3]/F')

    eleweight_res = array('f', 3*[1])
    eleweight_res[0], eleweight_res[1], eleweight_res[2] = 1., 1., 1.
    newtree.Branch('eleweight_res',eleweight_res,'eleweight_res[3]/F')
    
    eleweight_mva80 = array('f', 3*[1])
    eleweight_mva80[0], eleweight_mva80[1], eleweight_mva80[2] = 1., 1., 1.
    newtree.Branch('eleweight_mva80',eleweight_mva80,'eleweight_mva80[3]/F')
    
    eleweight_mva90 = array('f', 3*[1])
    eleweight_mva90[0], eleweight_mva90[1], eleweight_mva90[2] = 1., 1., 1.
    newtree.Branch('eleweight_mva90',eleweight_mva90,'eleweight_mva90[3]/F')

    emuweight = array('f', 3*[1])
    emuweight[0], emuweight[1], emuweight[2] = 1., 1., 1.
    newtree.Branch('emuweight',emuweight,'emuweight[3]/F')
    
    emuweight_muId = array('f', 3*[1])
    emuweight_muId[0], emuweight_muId[1], emuweight_muId[2] = 1., 1., 1.
    newtree.Branch('emuweight_muId',emuweight_muId,'emuweight_muId[3]/F')
    
    emuweight_muIso = array('f', 3*[1])
    emuweight_muIso[0], emuweight_muIso[1], emuweight_muIso[2] = 1., 1., 1.
    newtree.Branch('emuweight_muIso',emuweight_muIso,'emuweight_muIso[3]/F')
    
    emuweight_muTrk = array('f', 3*[1])
    emuweight_muTrk[0], emuweight_muTrk[1], emuweight_muTrk[2] = 1., 1., 1.
    newtree.Branch('emuweight_muTrk',emuweight_muTrk,'emuweight_muTrk[3]/F')
    
    emuweight_eleRes = array('f', 3*[1])
    emuweight_eleRes[0], emuweight_eleRes[1], emuweight_eleRes[2] = 1., 1., 1.
    newtree.Branch('emuweight_eleRes',emuweight_eleRes,'emuweight_eleRes[3]/F')
    
    emuweight_eleMVA90 = array('f', 3*[1])
    emuweight_eleMVA90[0], emuweight_eleMVA90[1], emuweight_eleMVA90[2] = 1., 1., 1.
    newtree.Branch('emuweight_eleMVA90',emuweight_eleMVA90,'emuweight_eleMVA90[3]/F')
    
    emuweight_trig = array('f', 3*[1])
    emuweight_trig[0], emuweight_trig[1], emuweight_trig[2] = 1., 1., 1.
    newtree.Branch('emuweight_trig',emuweight_trig,'emuweight_trig[3]/F')
   
    is_ttsemi = array('i',[0])
    is_ttsemi[0] = 0
    newtree.Branch('is_ttsemi', is_ttsemi, 'is_ttsemi/I')
   
    ttsemi_lep_pt = array('f',[-1])
    ttsemi_lep_pt[0] = -1
    newtree.Branch('ttsemi_lep_pt',ttsemi_lep_pt,'ttsemi_lep_pt/F') #=0 for muon, 1 for electron
    
    ttsemi_lep_phi = array('f',[-10])
    ttsemi_lep_phi[0] = -10
    newtree.Branch('ttsemi_lep_phi',ttsemi_lep_phi,'ttsemi_lep_phi/F') #=0 for muon, 1 for electron

    ttsemi_lep_eta = array('f',[-10])
    ttsemi_lep_eta[0] = -10
    newtree.Branch('ttsemi_lep_eta',ttsemi_lep_eta,'ttsemi_lep_eta/F') #=0 for muon, 1 for electron
    ttsemi_lep_etaSc = array('f',[-10])
    ttsemi_lep_etaSc[0] = -10
    newtree.Branch('ttsemi_lep_etaSc',ttsemi_lep_etaSc,'ttsemi_lep_etaSc/F') #=0 for muon, 1 for electron
    
    ttsemi_lep_mass = array('f',[-1])
    ttsemi_lep_mass[0] = -1
    newtree.Branch('ttsemi_lep_mass',ttsemi_lep_mass,'ttsemi_lep_mass/F') #=0 for muon, 1 for electron
    
    ttsemi_lep_iso = array('f',[-1])
    ttsemi_lep_iso[0] = -1
    newtree.Branch('ttsemi_lep_iso',ttsemi_lep_iso,'ttsemi_lep_iso/F') #=0 for muon, 1 for electron
    
    ttsemi_lep_pdgId = array('i',[-1])
    ttsemi_lep_pdgId[0] = -1
    newtree.Branch('ttsemi_lep_pdgId',ttsemi_lep_pdgId,'ttsemi_lep_pdgId/I') #=0 for muon, 1 for electron

    ttsemi_lep_looseId = array('i',[-1])
    ttsemi_lep_looseId[0] = -1
    newtree.Branch('ttsemi_lep_looseId',ttsemi_lep_looseId,'ttsemi_lep_looseId/I') #=0 for muon, 1 for electron
    ttsemi_lep_mediumId = array('i',[-1])
    ttsemi_lep_mediumId[0] = -1
    newtree.Branch('ttsemi_lep_mediumId',ttsemi_lep_mediumId,'ttsemi_lep_mediumId/I') #=0 for muon, 1 for electron
    ttsemi_lep_tightId = array('i',[-1])
    ttsemi_lep_tightId[0] = -1
    newtree.Branch('ttsemi_lep_tightId',ttsemi_lep_tightId,'ttsemi_lep_tightId/I') #=0 for muon, 1 for electron
    
    ttsemi_lep_veto = array('i', 6*[-1])
    ttsemi_lep_veto[0],ttsemi_lep_veto[1],ttsemi_lep_veto[2] = -1, -1, -1
    ttsemi_lep_veto[3],ttsemi_lep_veto[4],ttsemi_lep_veto[5] = -1, -1, -1
    newtree.Branch('ttsemi_lep_veto',ttsemi_lep_veto,'ttsemi_lep_veto[6]/I') #=0 loose veto, 1 medium veto, 3=tight veto

    ttsemi_massChi2 = array('f',[0])
    ttsemi_topMass = array('f',[0])
    ttsemi_wMass = array('f',[0])
    ttsemi_idxJet = array('i',4*[-1])
    ttsemi_idxJet_sortWjetCSV = array('i',4*[-1])
    ttsemi_idxJet_sortWjetCtag = array('i',4*[-1])

    newtree.Branch('ttsemi_massChi2', ttsemi_massChi2, 'ttsemi_massChi2/F')
    newtree.Branch('ttsemi_topMass', ttsemi_topMass, 'ttsemi_topMass/F')
    newtree.Branch('ttsemi_wMass', ttsemi_wMass, 'ttsemi_wMass/F')
    newtree.Branch('ttsemi_idxJet', ttsemi_idxJet, 'ttsemi_idxJet[4]/I')
    newtree.Branch('ttsemi_idxJet_sortWjetCSV', ttsemi_idxJet_sortWjetCSV, 'ttsemi_idxJet_sortWjetCSV[4]/I')
    newtree.Branch('ttsemi_idxJet_sortWjetCtag', ttsemi_idxJet_sortWjetCtag, 'ttsemi_idxJet_sortWjetCtag[4]/I')
   
    ttsemiweight = array('f', 3*[1])
    ttsemiweight[0], ttsemiweight[1], ttsemiweight[2] = 1., 1., 1.
    newtree.Branch('ttsemiweight',ttsemiweight,'ttsemiweight[3]/F')

    ttsemiweight_trig = array('f', 3*[1])
    ttsemiweight_trig[0], ttsemiweight_trig[1], ttsemiweight_trig[2] = 1., 1., 1.
    newtree.Branch('ttsemiweight_trig',ttsemiweight_trig,'ttsemiweight_trig[3]/F')
    
    ttsemiweight_muId = array('f', 3*[1])
    ttsemiweight_muId[0], ttsemiweight_muId[1], ttsemiweight_muId[2] = 1., 1., 1.
    newtree.Branch('ttsemiweight_muId',ttsemiweight_muId,'ttsemiweight_muId[3]/F')

    ttsemiweight_muIso = array('f', 3*[1])
    ttsemiweight_muIso[0], ttsemiweight_muIso[1], ttsemiweight_muIso[2] = 1., 1., 1.
    newtree.Branch('ttsemiweight_muIso',ttsemiweight_muIso,'ttsemiweight_muIso[3]/F')

    ttsemiweight_muTrk = array('f', 3*[1])
    ttsemiweight_muTrk[0], ttsemiweight_muTrk[1], ttsemiweight_muTrk[2] = 1., 1., 1.
    newtree.Branch('ttsemiweight_muTrk',ttsemiweight_muTrk,'ttsemiweight_muTrk[3]/F')

    ttsemiweight_eleRes = array('f', 3*[1])
    ttsemiweight_eleRes[0], ttsemiweight_eleRes[1], ttsemiweight_eleRes[2] = 1., 1., 1.
    newtree.Branch('ttsemiweight_eleRes',ttsemiweight_eleRes,'ttsemiweight_eleRes[3]/F')

    ttsemiweight_eleMVA90 = array('f', 3*[1])
    ttsemiweight_eleMVA90[0], ttsemiweight_eleMVA90[1], ttsemiweight_eleMVA90[2] = 1., 1., 1.
    newtree.Branch('ttsemiweight_eleMVA90',ttsemiweight_eleMVA90,'ttsemiweight_eleMVA90[3]/F')
   
    idxJet_passPtEta = array('i',10*[-1])
    for i in range(10):
      idxJet_passPtEta[i] = -1
    newtree.Branch('idxJet_passPtEta',idxJet_passPtEta,'idxJet_passPtEta[10]/I')  #index of jets passing pt and eta cut 

    idxJet_passCSV = array('i',2*[-1])
    idxJet_passCSV[0] = -1  
    idxJet_passCSV[1] = -1 
    newtree.Branch('idxJet_passCSV',idxJet_passCSV,'idxJet_passCSV[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium

    idxJet_passCSV_SVT = array('i',2*[-1])
    idxJet_passCSV_SVT[0] = -1  
    idxJet_passCSV_SVT[1] = -1 
    newtree.Branch('idxJet_passCSV_SVT',idxJet_passCSV_SVT,'idxJet_passCSV_SVT[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium
    
    idxJet_passCSV_SVT_1 = array('i',2*[-1])
    idxJet_passCSV_SVT_1[0] = -1  
    idxJet_passCSV_SVT_1[1] = -1 
    newtree.Branch('idxJet_passCSV_SVT_1',idxJet_passCSV_SVT_1,'idxJet_passCSV_SVT_1[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium

    idxJet_passCSV_SVT_2 = array('i',2*[-1])
    idxJet_passCSV_SVT_2[0] = -1  
    idxJet_passCSV_SVT_2[1] = -1 
    newtree.Branch('idxJet_passCSV_SVT_2',idxJet_passCSV_SVT_2,'idxJet_passCSV_SVT_2[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium
    
    idxJet_passCtag = array('i',2*[-1])
    idxJet_passCtag[0] = -1  
    idxJet_passCtag[1] = -1 
    newtree.Branch('idxJet_passCtag',idxJet_passCtag,'idxJet_passCtag[2]/I')  #index of highest pT jet passing Ctag cut: first entry is for tight, second entry is for medium
 
    idxJet_passCtag_SVT = array('i',2*[-1])
    idxJet_passCtag_SVT[0] = -1  
    idxJet_passCtag_SVT[1] = -1 
    newtree.Branch('idxJet_passCtag_SVT',idxJet_passCtag_SVT,'idxJet_passCtag_SVT[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium

    isCR_invertedCSV = array('i',2*[-1])
    isCR_invertedCSV[0] = 0
    isCR_invertedCSV[1] = 0
    newtree.Branch('isCR_invertedCSV',isCR_invertedCSV,'isCR_invertedCSV[2]/I') # flag to category the events in control regions by inverting CSV cut, if event doesn't pass CSV cuts ie idxJet_passCSV[0] or idxJet_passCSV[1] == -1, it is categoried as CR 
    idxJet_invertedCSV_SVT = array('i',2*[-1])
    idxJet_invertedCSV_SVT[0] = -1  
    idxJet_invertedCSV_SVT[1] = -1 
    newtree.Branch('idxJet_invertedCSV_SVT',idxJet_invertedCSV_SVT,'idxJet_invertedCSV_SVT[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium
    
    idxJet_invertedCSV_SVT_1 = array('i',2*[-1])
    idxJet_invertedCSV_SVT_1[0] = -1  
    idxJet_invertedCSV_SVT_1[1] = -1 
    newtree.Branch('idxJet_invertedCSV_SVT_1',idxJet_invertedCSV_SVT_1,'idxJet_invertedCSV_SVT_1[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium

    idxJet_invertedCSV_SVT_2 = array('i',2*[-1])
    idxJet_invertedCSV_SVT_2[0] = -1  
    idxJet_invertedCSV_SVT_2[1] = -1 
    newtree.Branch('idxJet_invertedCSV_SVT_2',idxJet_invertedCSV_SVT_2,'idxJet_invertedCSV_SVT_2[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium

    nCSVjet = array('i', [0])
    nCSVjet[0] = 0
    newtree.Branch('nCSVjet',nCSVjet,'nCSVjet/I')  #number of jet with CSV >= 0
    
    idxJet_sortedCSV = array('i', 20*[-1])
    newtree.Branch('idxJet_sortedCSV',idxJet_sortedCSV,'idxJet_sortedCSV[nCSVjet]/I')  #number of jet with CSV >= 0
    
    idxJet_sortedCSV_SVT = array('i',2*[-1])
    idxJet_sortedCSV_SVT[0] = -1  
    idxJet_sortedCSV_SVT[1] = -1 
    newtree.Branch('idxJet_sortedCSV_SVT',idxJet_sortedCSV_SVT,'idxJet_sortedCSV_SVT[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium
    
    bTagWeight_CSVM = array('f', 3*[1])
    newtree.Branch('bTagWeight_CSVM',bTagWeight_CSVM,'bTagWeight_CSVM[3]/F')
    bTagWeight_CSVT = array('f', 3*[1])
    newtree.Branch('bTagWeight_CSVT',bTagWeight_CSVT,'bTagWeight_CSVT[3]/F') 
    
    cTagWeight_CSVM = array('f', 3*[1])
    newtree.Branch('cTagWeight_CSVM',cTagWeight_CSVM,'cTagWeight_CSVM[3]/F')
    cTagWeight_CSVT = array('f', 3*[1])
    newtree.Branch('cTagWeight_CSVT',cTagWeight_CSVT,'cTagWeight_CSVT[3]/F') 
    ###############################
    #for gluon splitting
    ###############################
    Jet_gbb_weight = array('f',100*[-1])
    newtree.Branch('Jet_gbb_weight',Jet_gbb_weight,'Jet_gbb_weight[tree.nJet]/F')
    Jet_gcc_weight = array('f',100*[-1])
    newtree.Branch('Jet_gcc_weight',Jet_gcc_weight,'Jet_gcc_weight[tree.nJet]/F')

    ###############################
    #correct for Vtype stuffs
    ###############################
    Vtype_new = array('f',[0])
    newtree.Branch('Vtype_new',Vtype_new,'Vtype_new/F')

    vLeptonsBranches={}
    VBranches={}
    ##define Vleptons branch
    vLeptonsvar = ['pt', 'eta', 'etaSc', 'phi', 'mass', 'relIso03', 'relIso04', 'pdgId', 'looseIdPOG','mediumIdPOG_ICHEP2016','tightId','pfRelIso04','pfRelIso03']
    for var in vLeptonsvar:
        #vLeptonsBranches[var] = np.array([0]*2, dtype=float)
        vLeptonsBranches[var] = np.zeros(21, dtype=np.float32)
        obranch = newtree.Branch('vLeptons_new_%s'%var, vLeptonsBranches[var], 'vLeptons_new_%s[2]/F'%var)

    ##define Vleptons branch
    Vvar = ['pt', 'eta', 'phi', 'mass']
    LorentzDic = {'pt':'Pt','eta':'Eta','phi':'Phi','mass':'M'}
    for var in Vvar:
        #vLeptonsBranches[var] = np.array([0]*2, dtype=float)
        VBranches[var] = np.zeros(21, dtype=np.float32)
        obranch = newtree.Branch('V_new_%s'%var, VBranches[var], 'V_new_%s/F'%var)

#include the Vytpe reco here
    zEleSelection = lambda x : tree.selLeptons_pt[x] > 15 and tree.selLeptons_eleMVAIdSppring16GenPurp[x] >= 1
    zMuSelection = lambda x : tree.selLeptons_pt[x] > 15 and  tree.selLeptons_looseIdPOG[x] and tree.selLeptons_relIso04[x] < 0.25

    

    print 'starting event loop, processing',str(nEntries),'events'
    j_out=10000;

    #########################
    #Start event loop
    #########################
    #TEMP
    if debug: 
      nEntries = 5000
      print '>>>>>>>>>>>: ', useVtypeNew


    #nEntries = 100000
     

    nFilled = 0
    
    for entry in range(0,nEntries):
      # if entry>1000: break
      if ((entry%j_out)==0):
          if ((entry/j_out)==9 and j_out < 1e4): j_out*=10;
          print strftime("%Y-%m-%d %H:%M:%S", gmtime()),' - processing event',str(entry)+'/'+str(nEntries), '(cout every',j_out,'events)'
          sys.stdout.flush()
      
      #=====> Get entry
      tree.GetEntry(entry)
     
      ###############################
      #fixing Vtype
      ###############################

      #Variable to store Vtype and leptons info
      Vtype_new_ = -1
      V_mass_new = -1

      vLeptons_new = []
      #get all the lepton index
      lep_index = range(len(tree.selLeptons_pt))
      selectedElectrons = [i for i in  lep_index if abs(tree.selLeptons_pdgId[i]) == 11]
      selectedMuons =  [i for i in lep_index if abs(tree.selLeptons_pdgId[i]) == 13]

      zElectrons=[x for x in selectedElectrons if zEleSelection(x)]
      zMuons=[x for x in selectedMuons if zMuSelection(x) ]

      zMuons.sort(key=lambda x:tree.selLeptons_pt[x],reverse=True)
      zElectrons.sort(key=lambda x:tree.selLeptons_pt[x],reverse=True)

      if len(zMuons) >=  2 :
          if tree.selLeptons_pt[zMuons[0]] > 20:
              for i in zMuons[1:]:
                  if  tree.selLeptons_charge[zMuons[0]]*tree.selLeptons_charge[i] < 0:
                      #if tree.Vtype == 1:
                      Vtype_new_ = 0
                      for var in vLeptonsvar:
                          vLeptonsBranches[var][0] = getattr(tree,'selLeptons_%s'%var)[0]
                          vLeptonsBranches[var][1] = getattr(tree,'selLeptons_%s'%var)[i]
                      break
      elif len(zElectrons) >=  2 :
          if tree.selLeptons_pt[zElectrons[0]] > 20:
              for i in zElectrons[1:]:
                  if  tree.selLeptons_charge[zElectrons[0]]*tree.selLeptons_charge[i] < 0:
                      Vtype_new_ = 1
                      #if tree.Vtype == 0:
                      for var in vLeptonsvar:
                          vLeptonsBranches[var][0] = getattr(tree,'selLeptons_%s'%var)[0]
                          vLeptonsBranches[var][1] = getattr(tree,'selLeptons_%s'%var)[i]
                      break
      else:
          if tree.Vtype == 0 or tree.Vtype == 1:
              print '@ERROR: This is impossible, the new ele cut should be losser...'
              sys.exit(1)
          #add lepton if Vtype 2 or 3
          if tree.Vtype == 2 or tree.Vtype == 3:
              Vtype_new_ = tree.Vtype
              for var in vLeptonsvar:
                  vLeptonsBranches[var][0] = getattr(tree,'vLeptons_%s'%var)[0]
          #to handle missasigned Vtype 4 or -1 because of addtional electron cut
          elif (tree.Vtype == 4 or tree.Vtype == -1) and len(zElectrons) + len(zMuons) > 0:
              Vtype_new_ = 5
          #to handle missasigned Vtype 5 because of addtional electron cut
          elif tree.Vtype == 5 and len(zElectrons) + len(zMuons) == 0:
              if tree.met_pt < 80: 
                  Vtype_new_ = -1
              else: 
                  Vtype_new_ = 4 
          #if none of the exception above happen, it is save to copy the Vtype
          else:
              Vtype_new_ = tree.Vtype
           
      
      V = ROOT.TLorentzVector()
            
      if Vtype_new_ == 0 or Vtype_new_ == 1:
          lep1 = ROOT.TLorentzVector()
          lep2 = ROOT.TLorentzVector()
          lep1.SetPtEtaPhiM(vLeptonsBranches['pt'][0], vLeptonsBranches['eta'][0], vLeptonsBranches['phi'][0], vLeptonsBranches['mass'][0])
          lep2.SetPtEtaPhiM(vLeptonsBranches['pt'][1], vLeptonsBranches['eta'][1], vLeptonsBranches['phi'][1], vLeptonsBranches['mass'][1])
          V = lep1+lep2
          for var in Vvar:
              VBranches[var][0] = getattr(V,LorentzDic[var])()
      else: 
          for var in Vvar:
              VBranches[var][0] = getattr(tree,'V_%s'%var)

      Vtype_new[0] = Vtype_new_
     
      

      #################################
      #start other bussiness
      #################################

      if job.type != 'DATA':
          EventForTraining[0]=int(not TFlag.EvalInstance())
     
      #fill correctMass
      for i in range(0, tree.nJet):
        if tree.nprimaryVertices > 0:
          Jet_vtxMassCorr1[i] = util_funcs.VtxMassCorr(tree.primaryVertices_x[0], tree.primaryVertices_y[0], tree.primaryVertices_z[0], tree.Jet_vtxPosX[i], tree.Jet_vtxPosY[i], tree.Jet_vtxPosZ[i], tree.Jet_vtxPx[i], tree.Jet_vtxPy[i], tree.Jet_vtxPz[i], tree.Jet_vtxMass[i]) 
        else:
          Jet_vtxMassCorr1[i] = -1

      #find emu event
      #initial emu event
      is_emu[0] = 0
      emu_lep_pt[0], emu_lep_pt[1] = -1, -1
      emu_lep_eta[0], emu_lep_eta[1] = -10, -10
      emu_lep_etaSc[0], emu_lep_etaSc[1] = -10, -10
      emu_lep_phi[0], emu_lep_phi[1] = -10, -10
      emu_lep_mass[0], emu_lep_mass[1] = -1, -1
      emu_lep_charge[0], emu_lep_charge[1] = -1, -1
      emu_lep_pdgId[0], emu_lep_pdgId[1] = -1, -1
      emu_lep_looseId[0], emu_lep_looseId[1] = -1, -1
      emu_lep_mediumId[0], emu_lep_mediumId[1] = -1, -1
      emu_lep_tightId[0], emu_lep_tightId[1] = -1, -1
      emu_lep_iso[0], emu_lep_iso[1] = -1, -1

      #find tt semi lepton event
      #inital tt semi lepton event 
      is_ttsemi[0] = 0
      ttsemi_lep_pt[0] = -1
      ttsemi_lep_eta[0] = -10
      ttsemi_lep_etaSc[0] = -10
      ttsemi_lep_phi[0] = -10
      ttsemi_lep_mass[0] = -1
      ttsemi_lep_pdgId[0] = -1
      ttsemi_lep_looseId[0] = -1
      ttsemi_lep_mediumId[0] = -1
      ttsemi_lep_tightId[0] = -1
      ttsemi_lep_iso[0] = -1
      ttsemi_lep_veto[0],ttsemi_lep_veto[1],ttsemi_lep_veto[2] = 1, 1, 1
      ttsemi_lep_veto[3],ttsemi_lep_veto[4],ttsemi_lep_veto[5] = 1, 1, 1
      ttsemi_massChi2[0] = -1
      ttsemi_topMass[0] = -1
      ttsemi_wMass[0] = -1
      for i in range(0,4):
        ttsemi_idxJet[i] = -1
        ttsemi_idxJet_sortWjetCSV[i] = -1
        ttsemi_idxJet_sortWjetCtag[i] = -1
      
      lepTmp_idx = []
      lepTmp_pt = []
      lepTmp_eta = []
      lepTmp_etaSc = []
      lepTmp_phi = []
      lepTmp_mass = []
      lepTmp_charge = []
      lepTmp_pdgId = []
      lepTmp_looseId = []
      lepTmp_mediumId = []
      lepTmp_tightId = []
      lepTmp_idEle = []
      lepTmp_iso = []
      iLepTmp = -1
      for i in range(tree.nselLeptons):
        passId = False
        passIso = False
        if abs(tree.selLeptons_pdgId[i]) == 13:
          passId = tree.selLeptons_looseIdPOG[i] >= 1
          passIso = tree.selLeptons_pfRelIso04[i] < 0.25 #loose Iso
        if abs(tree.selLeptons_pdgId[i]) == 11:
          passId = (tree.selLeptons_eleMVAIdSppring16GenPurp[i] >= 1) #WP90, see VHbbAnalysis/Heppy/python/vhbbobj.py
          passIso = tree.selLeptons_pfRelIso03[i] < 0.25
          #passIso = tree.selLeptons_dr03TkSumPt[i] < 0.18 #see VHbbAnalysis/Heppy/test/vhbb.py

        if tree.selLeptons_pt[i] > 25 and abs(tree.selLeptons_eta[i]) < 2.4 and passId and passIso: 
          iLepTmp = iLepTmp + 1
          lepTmp_idx.append(iLepTmp)
          lepTmp_pt.append(tree.selLeptons_pt[i])
          lepTmp_eta.append(tree.selLeptons_eta[i])
          lepTmp_etaSc.append(tree.selLeptons_etaSc[i])
          lepTmp_phi.append(tree.selLeptons_phi[i])
          lepTmp_mass.append(tree.selLeptons_mass[i])
          lepTmp_charge.append(tree.selLeptons_charge[i])
          lepTmp_pdgId.append(tree.selLeptons_pdgId[i])
          lepTmp_looseId.append(tree.selLeptons_looseIdPOG[i]) #only valid for muon, for electron always = 1
          lepTmp_mediumId.append(tree.selLeptons_mediumIdPOG_ICHEP2016[i]) #only valid for muon, for electron always = 1 
          lepTmp_tightId.append(tree.selLeptons_tightId[i]) #only valid for muon, for electron always = 1
          if abs(tree.selLeptons_pdgId[i]) == 13:
            lepTmp_iso.append(tree.selLeptons_pfRelIso04[i])
            lepTmp_idEle.append(1)
          elif abs(tree.selLeptons_pdgId[i]) == 11:
            lepTmp_iso.append(tree.selLeptons_pfRelIso03[i])
            lepTmp_idEle.append(tree.selLeptons_eleMVAIdSppring16GenPurp[i])
          else:
            lepTmp_iso.append(-1)
            lepTmp_idEle.append(-1)
      
      util_funcs.sortPt(lepTmp_idx, lepTmp_pt)

      if len(lepTmp_idx) >= 2:
        iLepTmp1 = lepTmp_idx[0]
        iLepTmp2 = lepTmp_idx[1]
        muIdx = -1
        eIdx = -1
        if (abs(lepTmp_pdgId[iLepTmp1]) == 11 and abs(lepTmp_pdgId[iLepTmp2]) == 13):
            is_emu[0] = 1
            muIdx = iLepTmp2
            eIdx = iLepTmp1
        if (abs(lepTmp_pdgId[iLepTmp1]) == 13 and abs(lepTmp_pdgId[iLepTmp2]) == 11):
            is_emu[0] = 1
            muIdx = iLepTmp1
            eIdx = iLepTmp2
         
        if is_emu[0] == 1:
           emu_lep_pt[0], emu_lep_pt[1] = lepTmp_pt[muIdx], lepTmp_pt[eIdx]
           emu_lep_eta[0], emu_lep_eta[1] = lepTmp_eta[muIdx], lepTmp_eta[eIdx]
           emu_lep_etaSc[0], emu_lep_etaSc[1] = -10, lepTmp_etaSc[eIdx]
           emu_lep_phi[0], emu_lep_phi[1] = lepTmp_phi[muIdx], lepTmp_phi[eIdx]
           emu_lep_mass[0], emu_lep_mass[1] = lepTmp_mass[muIdx], lepTmp_mass[eIdx]
           emu_lep_charge[0], emu_lep_charge[1] = lepTmp_charge[muIdx], lepTmp_charge[eIdx]
           emu_lep_iso[0], emu_lep_iso[1] = lepTmp_iso[muIdx], lepTmp_iso[eIdx]
           emu_lep_pdgId[0], emu_lep_pdgId[1] = lepTmp_pdgId[muIdx], lepTmp_pdgId[eIdx]
           emu_lep_looseId[0], emu_lep_looseId[1] = lepTmp_looseId[muIdx], lepTmp_looseId[eIdx]
           emu_lep_mediumId[0], emu_lep_mediumId[1] = lepTmp_mediumId[muIdx], lepTmp_mediumId[eIdx]
           emu_lep_tightId[0], emu_lep_tightId[1] = lepTmp_tightId[muIdx], lepTmp_tightId[eIdx]
           emu_lep_idEle[0], emu_lep_idEle[1] = lepTmp_idEle[muIdx], lepTmp_idEle[eIdx]
      
      passVtype = (tree.Vtype == 2 or tree.Vtype == 3)
      if useVtypeNew: 
        passVtype = (Vtype_new[0] == 2 or Vtype_new[0] == 3)

      if passVtype: #2 = muon chan, 3 = electron chan
          
        #loop over jets and select jets with pT > 25 GeV and eta < 2.4
        idxSelJets = []
        for iTmp in range(0,tree.nJet):
          if tree.Jet_pt[iTmp] > 25 and abs(tree.Jet_eta[iTmp]) < 2.4:
            idxSelJets.append(iTmp)

        if len(idxSelJets) >= 4:
          
          is_ttsemi[0] = 1
          #fill leading lepton
          if not useVtypeNew:
            ttsemi_lep_pt[0]= tree.vLeptons_pt[0]
            ttsemi_lep_eta[0]= tree.vLeptons_eta[0]
            ttsemi_lep_etaSc[0]= tree.vLeptons_etaSc[0]
            ttsemi_lep_phi[0]= tree.vLeptons_phi[0]
            ttsemi_lep_mass[0]= tree.vLeptons_mass[0]
            ttsemi_lep_pdgId[0]= tree.vLeptons_pdgId[0]
            ttsemi_lep_looseId[0]= tree.vLeptons_looseIdPOG[0]
            ttsemi_lep_mediumId[0]= tree.vLeptons_mediumIdPOG_ICHEP2016[0]
            ttsemi_lep_tightId[0]= tree.vLeptons_tightId[0]
            #ttsemi_lep_iso[0]= tree.vLeptons_pfRelIso04[0]
            if abs(tree.vLeptons_pdgId[0]) == 13:
              ttsemi_lep_iso[0] = tree.vLeptons_pfRelIso04[0]
            elif abs(tree.vLeptons_pdgId[0]) == 11:
              ttsemi_lep_iso[0] = tree.vLeptons_pfRelIso03[0]
            else:
              ttsemi_lep_iso[0] = -1
          
          #Vtype_new
          if useVtypeNew:
            ttsemi_lep_pt[0]= vLeptonsBranches['pt'][0]
            ttsemi_lep_eta[0]= vLeptonsBranches['eta'][0]
            ttsemi_lep_etaSc[0]= vLeptonsBranches['etaSc'][0]
            ttsemi_lep_phi[0]= vLeptonsBranches['phi'][0]
            ttsemi_lep_mass[0]= vLeptonsBranches['mass'][0]
            ttsemi_lep_pdgId[0]= vLeptonsBranches['pdgId'][0]
            ttsemi_lep_looseId[0]= vLeptonsBranches['looseIdPOG'][0]
            ttsemi_lep_mediumId[0]= vLeptonsBranches['mediumIdPOG_ICHEP2016'][0]
            ttsemi_lep_tightId[0]= vLeptonsBranches['tightId'][0]
            if abs(vLeptonsBranches['pdgId'][0]) == 13:
              ttsemi_lep_iso[0] = vLeptonsBranches['pfRelIso04'][0]
            elif abs(vLeptonsBranches['pdgId'][0]) == 11:
              ttsemi_lep_iso[0] = vLeptonsBranches['pfRelIso03'][0]
            else:
              ttsemi_lep_iso[0] = -1


          maxJetInd = min(len(idxSelJets),6)
          for i in range(0,maxJetInd): #b_t_1
            idxJet_i = idxSelJets[i]
            p_i = ROOT.TLorentzVector()
            p_i.SetPtEtaPhiM(tree.Jet_pt[idxJet_i], tree.Jet_eta[idxJet_i], tree.Jet_phi[idxJet_i], tree.Jet_mass[idxJet_i])
            for j in range(0,maxJetInd): #b_t_2 comes together with W
              idxJet_j = idxSelJets[j]
              if idxJet_j != idxJet_i:
                p_j = ROOT.TLorentzVector()
                p_j.SetPtEtaPhiM(tree.Jet_pt[idxJet_j], tree.Jet_eta[idxJet_j], tree.Jet_phi[idxJet_j], tree.Jet_mass[idxJet_j])
                
                tmp = []
                for zx in range(0, maxJetInd):
                  tmp.append(idxSelJets[zx])
                tmp.remove(idxJet_i)
                tmp.remove(idxJet_j)
                
                k = min(tmp) #leading W jet
                l = max(tmp) #subleading Wjet
                p_k = ROOT.TLorentzVector()
                p_k.SetPtEtaPhiM(tree.Jet_pt[k], tree.Jet_eta[k], tree.Jet_phi[k], tree.Jet_mass[k])
                p_l = ROOT.TLorentzVector()
                p_l.SetPtEtaPhiM(tree.Jet_pt[l], tree.Jet_eta[l], tree.Jet_phi[l], tree.Jet_mass[l])
                w_mass = (p_k + p_l).M()
                t_mass = (p_j + p_k + p_l).M()
                chi2 = pow((t_mass - 173.2),2) + pow((w_mass - 80.4),2)
                if ttsemi_massChi2[0] < 0 or chi2 < ttsemi_massChi2[0] :
                    ttsemi_massChi2[0] = chi2
                    ttsemi_topMass[0] = t_mass
                    ttsemi_wMass[0] = w_mass
                    ttsemi_idxJet[0] = idxJet_i
                    ttsemi_idxJet[1] = idxJet_j
                    ttsemi_idxJet[2] = k
                    ttsemi_idxJet[3] = l

                    ttsemi_idxJet_sortWjetCSV[0] = idxJet_i
                    ttsemi_idxJet_sortWjetCSV[1] = idxJet_j
                    if tree.Jet_btagCSV[k] > tree.Jet_btagCSV[l]:
                      ttsemi_idxJet_sortWjetCSV[2] = k
                      ttsemi_idxJet_sortWjetCSV[3] = l
                    else:
                      ttsemi_idxJet_sortWjetCSV[2] = l
                      ttsemi_idxJet_sortWjetCSV[3] = k
                    
                    ttsemi_idxJet_sortWjetCtag[0] = idxJet_i
                    ttsemi_idxJet_sortWjetCtag[1] = idxJet_j
                    if tree.Jet_ctagVsL[k] > tree.Jet_ctagVsL[l]:
                      ttsemi_idxJet_sortWjetCtag[2] = k
                      ttsemi_idxJet_sortWjetCtag[3] = l
                    else:
                      ttsemi_idxJet_sortWjetCtag[2] = l
                      ttsemi_idxJet_sortWjetCtag[3] = k



      #find emu event
      #idx_emu[0] = -1
      #idx_emu[1] = -1


      #if tree.nselLeptons >= 2:
      #    if abs(tree.selLeptons_pdgId[0]) == 13 and abs(tree.selLeptons_pdgId[1]) == 11:
      #        idx_emu[0] = 1
      #        idx_emu[1] = 0
      #    if abs(tree.selLeptons_pdgId[1]) == 13 and abs(tree.selLeptons_pdgId[0]) == 11:
      #        idx_emu[0] = 0
      #        idx_emu[1] = 1
      
      nCSVjet[0] = 0
      selJet_idx = []
      selJet_idx1 = [] #with higher pT cut = 30
      for i in range(0, tree.nJet):
          if tree.Jet_pt[i] < float(config.get('General', 'jetPt')) or abs(tree.Jet_eta[i]) > float(config.get('General', 'absJetEta')): continue
          selJet_idx.append(i)
          if tree.Jet_pt[i] > 30: selJet_idx1.append(i)
          if tree.Jet_btagCSV[i] >= 0:
              idxJet_sortedCSV[nCSVjet[0]] = i
              nCSVjet[0] += 1
      
      #sort jet according to pt
      util_funcs.sortPt(selJet_idx, tree.Jet_pt)

      #store index of jet passing Pt and Eta requirement
      for i in range(len(idxJet_passPtEta)): idxJet_passPtEta[i] = -1
      for i in range(len(selJet_idx1)):
        if i < len(idxJet_passPtEta): idxJet_passPtEta[i] = selJet_idx1[i]

      idxJet_passCSV[0] = -1
      idxJet_passCSV[1] = -1
      #for AVR 
      idxJet_passCSV_SVT[0] = -1
      idxJet_passCSV_SVT[1] = -1
      #for uncorrected IVF
      idxJet_passCSV_SVT_1[0] = -1 
      idxJet_passCSV_SVT_1[1] = -1
      #for corrected IVF
      idxJet_passCSV_SVT_2[0] = -1 
      idxJet_passCSV_SVT_2[1] = -1
      
      idxJet_passCtag[0] = -1
      idxJet_passCtag[1] = -1
      idxJet_passCtag_SVT[0] = -1
      idxJet_passCtag_SVT[1] = -1

      for iJ in selJet_idx:
        if tree.Jet_btagCSV[iJ] >= float(config.get('General', 'CSV_Tight')):
          if idxJet_passCSV[0] == -1: idxJet_passCSV[0] = iJ
          if tree.Jet_vtxMass[iJ] > 0 and idxJet_passCSV_SVT[0] == -1:
            idxJet_passCSV_SVT[0] = iJ
          if tree.Jet_incVtxMass[iJ] > 0 and idxJet_passCSV_SVT_1[0] == -1:
            idxJet_passCSV_SVT_1[0] = iJ
          if tree.Jet_vtxMassCorr_IVF[iJ] > 0 and idxJet_passCSV_SVT_2[0] == -1:
            idxJet_passCSV_SVT_2[0] = iJ
        if tree.Jet_btagCSV[iJ] >= float(config.get('General', 'CSV_Medium')):
          if idxJet_passCSV[1] == -1: idxJet_passCSV[1] = iJ
          if tree.Jet_vtxMass[iJ] > 0 and idxJet_passCSV_SVT[1] == -1:
            idxJet_passCSV_SVT[1] = iJ
          if tree.Jet_incVtxMass[iJ] > 0 and idxJet_passCSV_SVT_1[1] == -1:
            idxJet_passCSV_SVT_1[1] = iJ
          if tree.Jet_vtxMassCorr_IVF[iJ] > 0 and idxJet_passCSV_SVT_2[1] == -1:
            idxJet_passCSV_SVT_2[1] = iJ
        if tree.Jet_ctagVsL[iJ] >= 0.69 and tree.Jet_ctagVsB[iJ] >= -0.45: #tight
          if idxJet_passCtag[0] == -1: idxJet_passCtag[0] = iJ
          if tree.Jet_vtxMass[iJ] > 0 and idxJet_passCtag_SVT[0] == -1:
            idxJet_passCtag_SVT[0] = iJ
        if tree.Jet_ctagVsL[iJ] >= -0.1 and tree.Jet_ctagVsB[iJ] >= 0.08: #medium
          if idxJet_passCtag[1] == -1: idxJet_passCtag[1] = iJ
          if tree.Jet_vtxMass[iJ] > 0 and idxJet_passCtag_SVT[1] == -1:
            idxJet_passCtag_SVT[1] = iJ
      
      #sort tagged jet according to CSV
      util_funcs.sortPtFatjet(idxJet_sortedCSV,tree.Jet_btagCSV,nCSVjet[0]) #sorting csv is similar to sort for jet pt., sorting applied to the first nCSVjet
      
      idxJet_sortedCSV_SVT[0] = -1
      idxJet_sortedCSV_SVT[1] = -1
      
      for iTmp in range(0, nCSVjet[0]):
        iJ = idxJet_sortedCSV[iTmp]
        if tree.Jet_btagCSV[iJ] >= float(config.get('General', 'CSV_Tight')) and tree.Jet_vtxMass[iJ] > 0 and idxJet_sortedCSV_SVT[0] == -1:
          idxJet_sortedCSV_SVT[0] = iJ
        if tree.Jet_btagCSV[iJ] >= float(config.get('General', 'CSV_Medium')) and tree.Jet_vtxMass[iJ] > 0 and idxJet_sortedCSV_SVT[1] == -1:
          idxJet_sortedCSV_SVT[1] = iJ
      
      #flag for control region using inverted CSV cut
      isCR_invertedCSV[0] = 0 if idxJet_passCSV[0] != -1 else 1
      isCR_invertedCSV[1] = 0 if idxJet_passCSV[1] != -1 else 1
      #loop over selected jets to find the first on with positive vtxMass 
      #for AVR 
      idxJet_invertedCSV_SVT[0] = -1
      idxJet_invertedCSV_SVT[1] = -1
      #for uncorrected IVF
      idxJet_invertedCSV_SVT_1[0] = -1 
      idxJet_invertedCSV_SVT_1[1] = -1
      #for corrected IVF
      idxJet_invertedCSV_SVT_2[0] = -1 
      idxJet_invertedCSV_SVT_2[1] = -1
      for iJ in selJet_idx:
        if isCR_invertedCSV[0] == 1:
          if tree.Jet_vtxMass[iJ] > 0 and idxJet_invertedCSV_SVT[0] == -1:
            idxJet_invertedCSV_SVT[0] = iJ
          if tree.Jet_incVtxMass[iJ] > 0 and idxJet_invertedCSV_SVT_1[0] == -1:
            idxJet_invertedCSV_SVT_1[0] = iJ
          if tree.Jet_vtxMassCorr_IVF[iJ] > 0 and idxJet_invertedCSV_SVT_2[0] == -1:
            idxJet_invertedCSV_SVT_2[0] = iJ
        if isCR_invertedCSV[1] == 1:
          if tree.Jet_vtxMass[iJ] > 0 and idxJet_invertedCSV_SVT[1] == -1:
            idxJet_invertedCSV_SVT[1] = iJ
          if tree.Jet_incVtxMass[iJ] > 0 and idxJet_invertedCSV_SVT_1[1] == -1:
            idxJet_invertedCSV_SVT_1[1] = iJ
          if tree.Jet_vtxMassCorr_IVF[iJ] > 0 and idxJet_invertedCSV_SVT_2[1] == -1:
            idxJet_invertedCSV_SVT_2[1] = iJ

      for iJet in range(tree.nJet):
        Jet_gbb_weight[iJet] = 1 
        Jet_gcc_weight[iJet] = 1 

      if job.type != 'DATA':
        for iJet in range(tree.nJet):
          iGenMatch = findMatchGenJet(tree.Jet_eta[iJet], tree.Jet_phi[iJet], tree)
          if iGenMatch >=0:
            if tree.GenJet_numBHadrons[iGenMatch] >= 2: Jet_gbb_weight[iJet] = 1.5
            if tree.GenJet_numCHadrons[iGenMatch] >= 2: Jet_gcc_weight[iJet] = 1.5

      #bTag scale factor
      bcTagWeights = [[1,1,1], [1,1,1],[1,1,1],[1,1,1]] #CSVT,CSVM,CtagT,CtagM 
      if not tree.isData:
        bcTagWeights = util_funcs.calBtagWeight(selJet_idx, tree, btagSF_reader_T, btagSF_reader_M, ctagSF_reader_T, ctagSF_reader_M,debug)
      for iTmp in range(0, 3):
        bTagWeight_CSVT[iTmp] = bcTagWeights[0][iTmp]
        bTagWeight_CSVM[iTmp] = bcTagWeights[1][iTmp]
        cTagWeight_CSVT[iTmp] = bcTagWeights[2][iTmp]
        cTagWeight_CSVM[iTmp] = bcTagWeights[3][iTmp]
       

      ##########################
      # Adding mu SFs LOOK mu SF use specialweight?
      ##########################
      if not job.specialweight:
          DY_specialWeight[0] = 1
      else :
          #specialWeight = ROOT.TTreeFormula('specialWeight',str(job.specialweight), tree) --> move outside of the loop otherwise memory problem
          specialWeight_ = specialWeight.EvalInstance()
          if entry%10000 == 0: print specialWeight_
          DY_specialWeight[0] = specialWeight_
      
      #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      # lepton scale factor
      #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
      #Muon ID + iso + trk, 0: central, 1: up, 2: down 
      muweight[0], muweight[1], muweight[2] = 1., 1., 1. 
      muweight_id[0],  muweight_id[1],  muweight_id[2] = 1., 1., 1.
      muweight_iso[0],  muweight_iso[1],  muweight_iso[2] = 1., 1., 1.
      muweight_trk[0],  muweight_trk[1],  muweight_trk[2] = 1., 1., 1.
      #Muon in heppy use looseID  
      muweight_trig[0],  muweight_trig[1],  muweight_trig[2] = 1., 1., 1.
      
      #electron channel
      eleweight[0], eleweight[1], eleweight[2] = 1., 1., 1. 
      eleweight_res[0], eleweight_res[1], eleweight_res[2] = 1., 1., 1. 
      eleweight_mva80[0], eleweight_mva80[1], eleweight_mva80[2] = 1., 1., 1. 
      eleweight_mva90[0], eleweight_mva90[1], eleweight_mva90[2] = 1., 1., 1. 
      
      #id and iso weight for emu sample
      emuweight[0], emuweight[1], emuweight[2] = 1., 1., 1.
      emuweight_muId[0], emuweight_muId[1], emuweight_muId[2] = 1., 1., 1.
      emuweight_muIso[0], emuweight_muIso[1], emuweight_muIso[2] = 1., 1., 1.
      emuweight_muTrk[0], emuweight_muTrk[1], emuweight_muTrk[2] = 1., 1., 1.
      emuweight_eleRes[0], emuweight_eleRes[1], emuweight_eleRes[2] = 1., 1., 1.
      emuweight_eleMVA90[0], emuweight_eleMVA90[1], emuweight_eleMVA90[2] = 1., 1., 1.
      emuweight_trig[0], emuweight_trig[1], emuweight_trig[2] = 1., 1., 1.
      
      ttsemiweight[0], ttsemiweight[1], ttsemiweight[2] = 1., 1., 1.
      ttsemiweight_trig[0], ttsemiweight_trig[1], ttsemiweight_trig[2] = 1., 1., 1.
      ttsemiweight_muId[0], ttsemiweight_muId[1], ttsemiweight_muId[2] = 1., 1., 1.
      ttsemiweight_muIso[0], ttsemiweight_muIso[1], ttsemiweight_muIso[2] = 1., 1., 1.
      ttsemiweight_muTrk[0], ttsemiweight_muTrk[1], ttsemiweight_muTrk[2] = 1., 1., 1.
      ttsemiweight_eleRes[0], ttsemiweight_eleRes[1], ttsemiweight_eleRes[2] = 1., 1., 1.
      ttsemiweight_eleMVA90[0], ttsemiweight_eleMVA90[1], ttsemiweight_eleMVA90[2] = 1., 1., 1.

      if job.type != 'DATA':
 
        passVtype = (tree.Vtype == 0)
        if useVtypeNew: 
          passVtype = (Vtype_new[0] == 0)
        
        if passVtype:
          #recalculate SF for muon iso
          weight_id = []
          weight_iso = []
          weight_trk = []
          if tree.run < 278820:
            if not useVtypeNew:
              weight_id.append(muID_SF_reader_BtoF.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0]))
              weight_id.append(muID_SF_reader_BtoF.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
              weight_iso.append(muIso_SF_reader_BtoF.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0]))
              weight_iso.append(muIso_SF_reader_BtoF.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
              weight_trk.append(muTrk_SF_reader_BtoF.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0]))
              weight_trk.append(muTrk_SF_reader_BtoF.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
            else:
              weight_id.append(muID_SF_reader_BtoF.get_2D( vLeptonsBranches['pt'][0], vLeptonsBranches['eta'][0]))
              weight_id.append(muID_SF_reader_BtoF.get_2D( vLeptonsBranches['pt'][1], vLeptonsBranches['eta'][1]))
              weight_iso.append(muIso_SF_reader_BtoF.get_2D( vLeptonsBranches['pt'][0], vLeptonsBranches['eta'][0]))
              weight_iso.append(muIso_SF_reader_BtoF.get_2D( vLeptonsBranches['pt'][1], vLeptonsBranches['eta'][1]))
              weight_trk.append(muTrk_SF_reader_BtoF.get_2D( vLeptonsBranches['pt'][0], vLeptonsBranches['eta'][0]))
              weight_trk.append(muTrk_SF_reader_BtoF.get_2D( vLeptonsBranches['pt'][1], vLeptonsBranches['eta'][1]))

          else:
            if not useVtypeNew:
              weight_id.append(muID_SF_reader_GtoH.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0]))
              weight_id.append(muID_SF_reader_GtoH.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
              weight_iso.append(muIso_SF_reader_GtoH.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0]))
              weight_iso.append(muIso_SF_reader_GtoH.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
              weight_trk.append(muTrk_SF_reader_GtoH.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0]))
              weight_trk.append(muTrk_SF_reader_GtoH.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
            else:
              weight_id.append(muID_SF_reader_GtoH.get_2D( vLeptonsBranches['pt'][0], vLeptonsBranches['eta'][0]))
              weight_id.append(muID_SF_reader_GtoH.get_2D( vLeptonsBranches['pt'][1], vLeptonsBranches['eta'][1]))
              weight_iso.append(muIso_SF_reader_GtoH.get_2D( vLeptonsBranches['pt'][0], vLeptonsBranches['eta'][0]))
              weight_iso.append(muIso_SF_reader_GtoH.get_2D( vLeptonsBranches['pt'][1], vLeptonsBranches['eta'][1]))
              weight_trk.append(muTrk_SF_reader_GtoH.get_2D( vLeptonsBranches['pt'][0], vLeptonsBranches['eta'][0]))
              weight_trk.append(muTrk_SF_reader_GtoH.get_2D( vLeptonsBranches['pt'][1], vLeptonsBranches['eta'][1]))


          computeSF(muweight_id, weight_id)
          computeSF(muweight_iso, weight_iso)
          computeSF(muweight_trk, weight_trk)
          
          if debug:
            print "================"
            if not useVtypeNew:
              print "mu1 pt, eta: ", tree.vLeptons_pt[0], ' ', tree.vLeptons_eta[0]
              print "mu2 pt, eta: ", tree.vLeptons_pt[1], ' ', tree.vLeptons_eta[1]
            else:
              print "mu1 pt, eta: ", vLeptonsBranches['pt'][0], ' ', vLeptonsBranches['eta'][0]
              print "mu2 pt, eta: ", vLeptonsBranches['pt'][1], ' ', vLeptonsBranches['eta'][1]
            print "Iso weights: ", weight_iso
            print 'weight_SF_LooseIso: ', muweight_iso

          #muweight[0] = tree.vLeptons_SF_IdCutLoose[0]*tree.vLeptons_SF_IdCutLoose[1]*muweight_iso[0]*tree.vLeptons_SF_trk_eta[0]*tree.vLeptons_SF_trk_eta[1]
          #muweight[1] = (tree.vLeptons_SF_IdCutLoose[0]-tree.vLeptons_SFerr_IdCutLoose[0])*(tree.vLeptons_SF_IdCutLoose[1]-tree.vLeptons_SFerr_IdCutLoose[1])*muweight_iso[1]*(tree.vLeptons_SF_trk_eta[0]-tree.vLeptons_SFerr_trk_eta[0])*(tree.vLeptons_SF_trk_eta[1]-tree.vLeptons_SFerr_trk_eta[1])
          #muweight[2] = (tree.vLeptons_SF_IdCutLoose[0]+tree.vLeptons_SFerr_IdCutLoose[0])*(tree.vLeptons_SF_IdCutLoose[1]+tree.vLeptons_SFerr_IdCutLoose[1])*muweight_iso[2]*(tree.vLeptons_SF_trk_eta[0]+tree.vLeptons_SFerr_trk_eta[0])*(tree.vLeptons_SF_trk_eta[1]+tree.vLeptons_SFerr_trk_eta[1])   
          muweight[0] = muweight_id[0]*muweight_iso[0]*muweight_trk[0]
          muweight[1] = muweight_id[1]*muweight_iso[1]*muweight_trk[1]
          muweight[2] = muweight_id[2]*muweight_iso[2]*muweight_trk[2]

          weight_trig = []
          if not useVtypeNew:
            weight_trig.append(muTrig_Eff_reader_2.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0]))
            weight_trig.append(muTrig_Eff_reader_2.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
          else: 
            weight_trig.append(muTrig_Eff_reader_2.get_2D( vLeptonsBranches['pt'][0], vLeptonsBranches['eta'][0]))
            weight_trig.append(muTrig_Eff_reader_2.get_2D( vLeptonsBranches['pt'][1], vLeptonsBranches['eta'][1]))
          computeTrigWeight(muweight_trig,weight_trig)
          
          if debug:
            print "Trigger eff: ", weight_trig
            print 'muweight_trig: ', muweight_trig

        passVtype = (tree.Vtype == 1)
        if useVtypeNew: 
          passVtype = (Vtype_new[0] == 1)
        
        if passVtype:
          if not useVtypeNew:
            weight_eleRes = [getSF_2Dhis(tree.vLeptons_etaSc[0],tree.vLeptons_pt[0],eleRes_SF_his),getSF_2Dhis(tree.vLeptons_etaSc[1],tree.vLeptons_pt[1],eleRes_SF_his)]
            weight_ele_mva80 = [getSF_2Dhis(tree.vLeptons_etaSc[0],tree.vLeptons_pt[0],eleMVA80_SF_his),getSF_2Dhis(tree.vLeptons_etaSc[1],tree.vLeptons_pt[1],eleMVA80_SF_his)]
            weight_ele_mva90 = [getSF_2Dhis(tree.vLeptons_etaSc[0],tree.vLeptons_pt[0],eleMVA90_SF_his),getSF_2Dhis(tree.vLeptons_etaSc[1],tree.vLeptons_pt[1],eleMVA90_SF_his)]
          else: 
            weight_eleRes = [getSF_2Dhis(vLeptonsBranches['etaSc'][0],vLeptonsBranches['pt'][0],eleRes_SF_his),getSF_2Dhis(vLeptonsBranches['etaSc'][1],vLeptonsBranches['pt'][1],eleRes_SF_his)]
            weight_ele_mva80 = [getSF_2Dhis(vLeptonsBranches['etaSc'][0],vLeptonsBranches['pt'][0],eleMVA80_SF_his),getSF_2Dhis(vLeptonsBranches['etaSc'][1],vLeptonsBranches['pt'][1],eleMVA80_SF_his)]
            weight_ele_mva90 = [getSF_2Dhis(vLeptonsBranches['etaSc'][0],vLeptonsBranches['pt'][0],eleMVA90_SF_his),getSF_2Dhis(vLeptonsBranches['etaSc'][1],vLeptonsBranches['pt'][1],eleMVA90_SF_his)]
          
          computeSF(eleweight_res,weight_eleRes)
          computeSF(eleweight_mva80,weight_ele_mva80)
          computeSF(eleweight_mva90,weight_ele_mva90)
          
          eleweight[0] = eleweight_res[0]*eleweight_mva90[0]
          eleweight[1] = eleweight_res[1]*eleweight_mva90[1]
          eleweight[2] = eleweight_res[2]*eleweight_mva90[2]
           
        #@@@@@@@@@@@@@@@@@@@@@@@e-mu weight@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #e-mu weight = ID+iso+trk for muon and ID+iso+trk for electron

        if is_emu[0] == 1:
          #Iso
          sfTmp = (muIso_SF_reader_BtoF.get_2D(emu_lep_pt[0], emu_lep_eta[0]))
          if tree.run >= 278820:
            sfTmp = (muIso_SF_reader_GtoH.get_2D(emu_lep_pt[0], emu_lep_eta[0]))
          if sfTmp[0] > 0:
            emuweight[0] = emuweight[0]*sfTmp[0]
            emuweight_muIso[0], emuweight_muIso[1], emuweight_muIso[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])
          
          #ID
          sfTmp = (muID_SF_reader_BtoF.get_2D(emu_lep_pt[0], emu_lep_eta[0]))
          if tree.run >= 278820:
            sfTmp = (muID_SF_reader_GtoH.get_2D(emu_lep_pt[0], emu_lep_eta[0]))
          if sfTmp[0] > 0:
            emuweight[0] = emuweight[0]*sfTmp[0]
            emuweight_muId[0], emuweight_muId[1], emuweight_muId[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])
          
          #trk
          sfTmp = (muTrk_SF_reader_BtoF.get_2D(emu_lep_pt[0], emu_lep_eta[0]))
          if tree.run >= 278820:
            sfTmp = (muTrk_SF_reader_GtoH.get_2D(emu_lep_pt[0], emu_lep_eta[0]))
          if sfTmp[0] > 0: 
            emuweight[0] = emuweight[0]*sfTmp[0]
            emuweight_muTrk[0], emuweight_muTrk[1], emuweight_muTrk[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])

          #electron: res 
          sfTmp = getSF_2Dhis(emu_lep_etaSc[1],emu_lep_pt[1],eleMVA90_SF_his)
          if sfTmp[0] > 0:
            emuweight[0] = emuweight[0]*sfTmp[0]
            emuweight_eleRes[0], emuweight_eleRes[1], emuweight_eleRes[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])
          
          #electron: ID
          sfTmp = getSF_2Dhis(emu_lep_etaSc[1],emu_lep_pt[1],eleRes_SF_his)
          if sfTmp[0] > 0:
            emuweight[0] = emuweight[0]*sfTmp[0]
            emuweight_eleMVA90[0], emuweight_eleMVA90[1], emuweight_eleMVA90[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])
          
          #trigger weight
          sfTmp = (muTrig_Eff_reader_2.get_2D(emu_lep_pt[0], emu_lep_eta[0]))
          if sfTmp[0] > 0:
            emuweight_trig[0],emuweight_trig[1],emuweight_trig[2] = sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])
        
        #@@@@@@@@@@@@@@@@@@@@@@@ttsemi weight@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if is_ttsemi[0] == 1:
          if abs(ttsemi_lep_pdgId[0]) == 13: 
            #ID
            sfTmp = (muID_SF_reader_BtoF.get_2D(ttsemi_lep_pt[0], ttsemi_lep_eta[0]))
            if tree.run >= 278820:
              sfTmp = (muID_SF_reader_GtoH.get_2D(ttsemi_lep_pt[0], ttsemi_lep_eta[0]))
            if sfTmp[0] > 0:
              ttsemiweight[0] = ttsemiweight[0]*sfTmp[0]
              ttsemiweight_muId[0], ttsemiweight_muId[1], ttsemiweight_muId[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])
            
            #iso
            sfTmp = (muIso_SF_reader_BtoF.get_2D(ttsemi_lep_pt[0], ttsemi_lep_eta[0]))
            if tree.run >= 278820:
              sfTmp = (muIso_SF_reader_GtoH.get_2D(ttsemi_lep_pt[0], ttsemi_lep_eta[0]))
            if sfTmp[0] > 0:
              ttsemiweight[0] = ttsemiweight[0]*sfTmp[0]
              ttsemiweight_muIso[0], ttsemiweight_muIso[1], ttsemiweight_muIso[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])

            #trk
            sfTmp = (muTrk_SF_reader_BtoF.get_2D(ttsemi_lep_pt[0], ttsemi_lep_eta[0]))
            if tree.run >= 278820:
              sfTmp = (muTrk_SF_reader_GtoH.get_2D(ttsemi_lep_pt[0], ttsemi_lep_eta[0]))
            if sfTmp[0] > 0:
              ttsemiweight[0] = ttsemiweight[0]*sfTmp[0]
              ttsemiweight_muTrk[0], ttsemiweight_muTrk[1], ttsemiweight_muTrk[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])
            
            #trigger weight
            sfTmp = muTrig_Eff_reader_2.get_2D(ttsemi_lep_pt[0], ttsemi_lep_eta[0])
            if sfTmp[0] > 0:
              ttsemiweight_trig[0],ttsemiweight_trig[1],ttsemiweight_trig[2] = sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])
          
          if abs(ttsemi_lep_pdgId[0]) == 11:
            #Res
            sfTmp = getSF_2Dhis(ttsemi_lep_etaSc[0],ttsemi_lep_pt[0], eleRes_SF_his)
            if sfTmp[0] > 0:
              ttsemiweight[0] = ttsemiweight[0]*sfTmp[0]
              ttsemiweight_eleRes[0], ttsemiweight_eleRes[1], ttsemiweight_eleRes[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])
            #Id
            sfTmp = getSF_2Dhis(ttsemi_lep_etaSc[0],ttsemi_lep_pt[0], eleMVA90_SF_his)
            if sfTmp[0] > 0:
              ttsemiweight[0] = ttsemiweight[0]*sfTmp[0]
              ttsemiweight_eleMVA90[0], ttsemiweight_eleMVA90[1], ttsemiweight_eleMVA90[2]= sfTmp[0], (sfTmp[0] - sfTmp[1]), (sfTmp[0] + sfTmp[1])



        '''
        emu_id_iso_weight[0] = 1
        emu_trig_weight[0] = 1 #apply for muon in the leading two muon
        pTcut = 20

        if tree.Vtype == 0 or tree.Vtype == 1: 
            #=========> Add HLT trigger SF for muon=======
            debug = False

            # dR matching between lepton and trigger object
            #FIXME: need to update for other hlt object
            SF1 = tree.vLeptons_SF_HLT_RunD4p2[0]*0.1801911165 + tree.vLeptons_SF_HLT_RunD4p3[0]*0.8198088835 #SF_HLT == 1 for electron in 76X tree
            SF2 = tree.vLeptons_SF_HLT_RunD4p2[1]*0.1801911165 + tree.vLeptons_SF_HLT_RunD4p3[1]*0.8198088835 
            eff1 = tree.vLeptons_Eff_HLT_RunD4p2[0]*0.1801911165 + tree.vLeptons_Eff_HLT_RunD4p3[0]*0.8198088835 #Eff_HLT == 1 for electron in 76X tree
            eff2 = tree.vLeptons_Eff_HLT_RunD4p2[1]*0.1801911165 + tree.vLeptons_Eff_HLT_RunD4p3[1]*0.8198088835

            DR = [999, 999]
            for k in range(0,2):
                for l in range(0,len(tree.trgObjects_hltIsoMu18_eta)):
                    dr_ = util_funcs.deltaR(tree.vLeptons_eta[k], tree.vLeptons_phi[k], tree.trgObjects_hltIsoMu18_eta[l], tree.trgObjects_hltIsoMu18_phi[l])
                    if dr_ < DR[k] and tree.vLeptons_pt[k] > pTcut:
                        DR[k] = dr_

            Mu1pass = DR[0] < 0.5
            Mu2pass = DR[1] < 0.5
             
            if tree.Vtype == 0:
                if not Mu1pass and not Mu2pass:
                    vLeptons_SFweight_HLT[0] = 0
                elif Mu1pass and not Mu2pass:
                    vLeptons_SFweight_HLT[0] = SF1
                elif not Mu1pass and Mu2pass:
                    vLeptons_SFweight_HLT[0] = SF2
                elif Mu1pass and Mu2pass:
                    effdata = 1 - (1-SF1*eff1)*(1-SF2*eff2);
                    effmc = 1 - (1-eff1)*(1-eff2);
                    vLeptons_SFweight_HLT[0] = effdata/effmc
             
            DR[0] = 999
            DR[1] = 999

            for k in range(0,2):
              for l in range(0,len(tree.trgObjects_hltEle22eta2p1WPLoose_eta)):
                dr_ = util_funcs.deltaR(tree.vLeptons_eta[k], tree.vLeptons_phi[k], tree.trgObjects_hltEle22eta2p1WPLoose_eta[l], tree.trgObjects_hltEle22eta2p1WPLoose_eta[l])
                if dr_ < DR[k] and tree.vLeptons_pt[k] > 24:
                  DR[k] = dr_

            Ele1pass = DR[0] < 0.5
            Ele2pass = DR[0] < 0.5

            if tree.Vtype == 1:
                #vLeptons_SFweight_HLT[0] = 1
                if not Ele1pass and not Ele2pass:
                    vLeptons_SFweight_HLT[0] = 0
                elif Ele1pass and not Ele2pass:
                    vLeptons_SFweight_HLT[0] = SF1
                elif not Ele1pass and Ele2pass:
                    vLeptons_SFweight_HLT[0] = SF2
                elif Ele1pass and Ele2pass:
                    effdata = 1 - (1-SF1*eff1)*(1-SF2*eff2);
                    effmc = 1 - (1-eff1)*(1-eff2);
                    vLeptons_SFweight_HLT[0] = effdata/effmc

                        #print 'vLeptSFw afer fill is', vLeptons_SFweight_HLT[0]

        '''
        #calculate id, iso and trigger scale factor for emu event
        #if idx_emu[0] >= 0 and idx_emu[1] >= 0:
        #  SF1 = tree.selLeptons_SF_HLT_RunD4p2[idx_emu[1]]*0.1801911165 + tree.selLeptons_SF_HLT_RunD4p3[idx_emu[1]]*0.8198088835 #SF_HLT == 1 for electron in 76X tree ignore RunC
        #  dRtmp = 100
        #  for l in range(0,len(tree.trgObjects_hltIsoMu18_eta)):
        #        dr_ = util_funcs.deltaR(tree.selLeptons_eta[idx_emu[1]], tree.selLeptons_phi[idx_emu[1]], tree.trgObjects_hltIsoMu18_eta[l], tree.trgObjects_hltIsoMu18_phi[l])
        #        if dr_ < dRtmp and tree.selLeptons_pt[idx_emu[1]] > pTcut:
        #            dRtmp = dr_

        #  MuPass = dRtmp < 0.5
        #  if not MuPass: emu_trig_weight[0] = 0
        #  if MuPass: emu_trig_weight[0] = SF1

        #  emu_id_iso_weight[0] = tree.selLeptons_SF_IsoLoose[idx_emu[0]]*tree.selLeptons_SF_IdCutLoose[idx_emu[0]]*tree.selLeptons_SF_IsoLoose[idx_emu[1]]*tree.selLeptons_SF_IdCutLoose[idx_emu[1]]
        
        #btag scale factor: not available need MC b-tagging efficiency
      

#@@@@@@@@@@@@@Fill tree here@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      newtree.Fill()
      nFilled += 1      


#@@@@@@@@@@@@@Apply the skims, selection cuts@@@@@@@@@@@@@@@@@@@@@@@@@@@
    output.cd()
    print "============================================================="
    print "Apply skim cuts to new tree: ", skimCut 

    newtree1 = newtree.CopyTree(skimCut) 

    print 'Total filled event: ', nFilled
    print 'Total events after skim: ', newtree1.GetEntries()

    #newtree.FlushBaskets()
    #newtree.AutoSave()
    newtree1.FlushBaskets()
    newtree1.AutoSave()

    print 'Save'
    output.Close()
    print 'Close'
  
    return outputfile

def fillTreeSingleInput(inputs):
  return fillTree(*inputs)



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#if namelist[0] == All, only process indentifier. the info contain all the sub sample with the same identifier. we only need to process once for a identifier
#processed_identifiers = []
for job in info:
    if namelist[0] != 'All' and not job.name in namelist and len([x for x in namelist if x==job.identifier])==0:
        print 'job.name',job.name,'and job.identifier',job.identifier,'not in namelist',namelist
        continue
#    if namelist[0] == 'All' and job.indentifier in processed_identifiers: continue
#    processed_identifiers.append(job.identifier)

 #skip if it is subsample
    if job.subsample:
        continue
    
    print '\t match - %s' %(job.name)

#define input and output file list    
#create file list
    #inFile_folder = pathIN + '/' + job.prefix+job.identifier
    #tmpD =  config.get('Directories','samplefiles').replace('/','') + '_afterPrepStep/'
    #os.system('mkdir ' + tmpD)
    #tmpFileList =  tmpD + job.prefix + job.identifier + '.txt'
    #os.system('rm -f ' + tmpFileList)
    #util_funcs.findSubFolders(inFile_folder,tmpFileList,False)
   
    tmpFileList = ''
    if run_locally == 'True':
      inFile_folder = pathIN + '/' + job.prefix+job.identifier
      if pathIN.find('.txt') != -1 or pathIN.find('.tex') != -1: #use a file as list of folder, in this case samples located in different folder
        inFile_folder = util_funcs.getInputFolder(pathIN, job.identifier)
        #lines = open(pathIN).readlines()
        #for l in lines:
        #  if '#' in l: continue
        #  if job.identifier in l:
        #    inFile_folder = l.split()[0] + '/' + job.identifier
        #    if 'root://cmseos.fnal.gov' in inFile_folder:
        #      inFile_folder = inFile_folder.split('root://cmseos.fnal.gov')[1]
        #    break
      print '>>>>>>>>>>>>>>>>>>FOLDER of root file: ', inFile_folder
 
      tmpD =  config.get('Directories','samplefiles').replace('/','') + '_afterPrepStep/'
      os.system('mkdir ' + tmpD)
      tmpFileList =  tmpD + job.prefix + job.identifier + '.txt'
      os.system('rm -f ' + tmpFileList)
      if '/store' in inFile_folder:
        util_funcs.findSubFolders(inFile_folder,tmpFileList,True)
      else:
        util_funcs.findSubFolders(inFile_folder,tmpFileList,False)

    else:
      tmpFileList = opts.names


    outputFolder = "%s/%s/" %(pathOUT,job.identifier)
    if run_locally == 'True':
      print "Remove output folder: ", outputFolder
      os.system('rm -rf ' + outputFolder)
     
    try:
        if run_locally == 'True': os.mkdir(outputFolder)
    except:
        pass
    

    input_processings = []
    fileList = []
    nFile_for_processing = int(config.get('Configuration','nFile_for_processing'))
    if run_locally == 'False':
      nFile_for_processing = -1
    if debug:
      nFile_for_processing = 1

    print '>>>>> nFile_for_processing: ', nFile_for_processing
    nFile = 0
    for line in open(tmpFileList).readlines():
      line = line.replace("\n","")
      line_tmp = line.split('/')
      tmpStr = job.identifier
      pos = -1
      for tmp in line_tmp:
        pos = pos + 1
        if tmp.find(job.identifier) != -1: break
      if pos != -1 and pos != len(line_tmp):
        pos = pos + 1
        tmpStr = line_tmp[pos]
      outTmp = line_tmp[len(line_tmp)-1].replace('.root', '_moreInfo.root')
      if tmpStr != '':
        outTmp = line_tmp[len(line_tmp)-1].replace('.root', '_' + tmpStr + '_moreInfo.root')
      if debug:
        print '>>>>>>>>>>>', line_tmp
        print '>>>>>>>>>>>', job.name, ' ', job.identifier
        print '>>>>>>>>>>>', tmpStr
        print '>>>>>>>>>>>', outTmp
      if run_locally == 'True': outTmp = outputFolder + '/' + outTmp
      input_processings.append((line,outTmp,job.addtreecut))
      fileList.append(line)
      nFile += 1
      if nFile_for_processing > 0 and nFile >= nFile_for_processing: break


    #print input_processings

 ## process the input list (using multiprocess)#######
    multiprocess=int(config.get('Configuration','nprocesses'))
    if run_locally == 'False': multiprocess = 1
    outputs = []
    if multiprocess>1:
      from multiprocessing import Pool
      p = Pool(multiprocess)
      outputs = p.map(fillTreeSingleInput,input_processings)
      p.close()
      #p.terminate()
      p.join()

    else:
      for input_ in input_processings:
              output = fillTreeSingleInput(input_)
              outputs.append(output)



#TEMP for LPC
    if run_locally == 'True':
      outfileName = pathOUT + '/' + job.identifier + '.root'
      print 'Remove old file before hadd: ', outfileName 
      os.system('rm -f ' + outfileName)
      command = "hadd -f -k " + outfileName + " " + outputFolder + '/*'
      print command
      os.system(command)


    '''
    input_processings = []
    for line in open(tmpFileList).readlines():
      line = line.replace("\n","")
      line_tmp = line.split('/')
      outTmp = line_tmp[len(line_tmp)-1].replace('.root', '_moreInfo.root')
      outTmp = outputFolder + '/' + outTmp
      input_processings.append((line,outTmp))


    print input_processings

 ## process the input list (using multiprocess)#######
    multiprocess=int(config.get('Configuration','nprocesses'))
    outputs = []
    if multiprocess>1:
        from multiprocessing import Pool
        if len(input_processings) < multiprocess: multiprocess = len(input_processings)
        print 'Number of multiprocess is: ', multiprocess
        p = Pool(multiprocess)
        outputs = p.map(fillTreeSingleInput,input_processings)
        p.close()
        p.join()

    else:
        for input_ in input_processings:
                output = fillTreeSingleInput(input_)
                outputs.append(output)



#TEMP for LPC
    
    whereToLaunch = config.get('Configuration','whereToLaunch')
    if ('LPC' in whereToLaunch or 'pisa' in whereToLaunch):
      if run_locally != 'False':
        outfileName = pathOUT + '/' + job.identifier + '.root'
        print 'Remove old file before hadd: ', outfileName 
        os.system('rm -f ' + outfileName)
        command = "hadd -f -k " + outfileName + " " + outputFolder + '/*'
        print command
        os.system(command)
    
    '''

    




    '''inputfiles = []
    outputfiles = []
    tmpfiles = [] #temporary store the output before copy output to final destination
    if len(filelist) == 0:
        inputfiles.append(pathIN+'/'+job.prefix+job.identifier+'.root')
        print('opening '+pathIN+'/'+job.prefix+job.identifier+'.root')
        tmpfiles.append(tmpDir+'/'+job.prefix+job.identifier+'.root')
        outputfiles.append("%s/%s/%s" %(pathOUT,job.prefix,job.identifier+'.root'))
    else:
        for inputFile in filelist:
            subfolder = inputFile.split('/')[-4]
            filename = inputFile.split('/')[-1]
            filename = filename.split('_')[0]+'_'+subfolder+'_'+filename.split('_')[1]
            hash = hashlib.sha224(filename).hexdigest()
            inputFile = "%s/%s/%s" %(pathIN,job.identifier,filename.replace('.root','')+'_'+str(hash)+'.root')
            if not os.path.isfile(inputFile): continue
            outputFile = "%s/%s/%s" %(pathOUT,job.identifier,filename.replace('.root','')+'_'+str(hash)+'.root')
            tmpfile = "%s/%s" %(tmpDir,filename.replace('.root','')+'_'+str(hash)+'.root')
            if inputFile in inputfiles: continue

#check if outputFile exists and is good, if it does don't add to the inputfiles which is the list of file to process            

            if os.path.isfile(outputFile):
                f = ROOT.TFile.Open(outputFile,'read')
                if not f:
                  print 'file is null, adding to input'
                  inputfiles.append(inputFile)
                  outputfiles.append(outputFile)
                  tmpfiles.append(tmpfile)
                  continue
                # f.Print()
                if f.GetNkeys() == 0 or f.TestBit(ROOT.TFile.kRecovered) or f.IsZombie():
                    print 'f.GetNkeys()',f.GetNkeys(),'f.TestBit(ROOT.TFile.kRecovered)',f.TestBit(ROOT.TFile.kRecovered),'f.IsZombie()',f.IsZombie()
                    print 'File', outputFile, 'already exists but is buggy, gonna delete and rewrite it.'
                    command = 'rm %s' %(outputFile)
                    subprocess.call([command], shell=True)
                    print(command)
                else: continue

            inputfiles.append(inputFile)
            outputfiles.append(outputFile)
            tmpfiles.append(tmpfile)

        print 'inputfiles',inputfiles,'tmpfiles',tmpfiles

#loop over inputfile and do something. The loop will be for each inputfile having one tmpfile and one outputFile. They are all go together     
    
    for inputfile,tmpfile,outputFile in zip(inputfiles,tmpfiles,outputfiles):
    '''
