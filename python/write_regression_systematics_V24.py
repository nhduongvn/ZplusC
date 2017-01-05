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
ROOT.gROOT.SetBatch(True)
from optparse import OptionParser
from btag_reweight import *
from time import gmtime, strftime
from multiprocessing import Pool

from myutils import BetterConfigParser, ParseInfo, TreeCache, LeptonSF, util_funcs


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#INPUT
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

debug = False 

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
anaTag = config.get("Analysis","tag")
TrainFlag = eval(config.get('Analysis','TrainFlag'))
btagLibrary = config.get('BTagReshaping','library')
samplesinfo=config.get('Directories','samplesinfo')
channel=config.get('Configuration','channel')
print 'channel is', channel

VHbbNameSpace=config.get('VHbbNameSpace','library')
ROOT.gSystem.Load(VHbbNameSpace)
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
#jsons = {
#  wdir+'/python/json/ScaleFactor_GsfElectronToRECO_passingTrigWP80.json' : ['ScaleFactor_GsfElectronToRECO_passingTrigWP80', 'eta_pt_ratio'],
#            #Not available for the moment
#  wdir+'/python/json/ScaleFactor_HLT_Ele23_WPLoose_Gsf_v.json' : ['ScaleFactor_HLT_Ele23_WPLoose_Gsf_v', 'eta_pt_ratio'],
#}
jsons = {
    wdir+'/python/json/MuonIso_Z_RunBCD_prompt80X_7p65.json' : ['MC_NUM_LooseRelIso_DEN_LooseID_PAR_pt_spliteta_bin1', 'abseta_pt_ratio'], #muonIso, recalculate muonIso since heppy ntuple store Iso SF with respect to tightID while vLepton is looseID
    #  wdir+'/python/json/HLT_Ele23_WPLoose.json' : ['ScaleFactor_egammaEff_WP80', 'eta_pt_ratio'],
    wdir+'/python/json/ScaleFactor_egammaEff_WP80.json' : ['ScaleFactor_egammaEff_WP80', 'pt_eta_ratio'],
    wdir+'/python/json/ScaleFactor_egammaEff_WP90.json' : ['ScaleFactor_egammaEff_WP90', 'eta_pt_ratio'],
    #  wdir+'/python/json/eff_Ele27_WPLoose_Eta2p1_RunBtoF.json' : ['Trigger_Eff', 'eta_pt_ratio'],
    #  wdir+'/python/json/egammaEffi_tracker.json' : ['egammaEffi_tracker', 'eta_pt_ratio'],
    #  wdir+'/python/json/SingleMuonTrigger_LooseMuons_beforeL2fix_Z_RunBCD_prompt80X_7p65.json' : ['MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_beforeL2Fix', 'abseta_pt_MC'],
    #  wdir+'/python/json/SingleMuonTrigger_LooseMuons_afterL2fix_Z_RunBCD_prompt80X_7p65.json' : ['MuonTrigger_data_all_IsoMu22_OR_IsoTkMu22_pteta_Run2016B_afterL2Fix', 'abseta_pt_MC'],
    #ID+ISO
    wdir+'/python/json/WP90PlusIso_BCD.json' : ['WP90PlusIso_BCD', 'eta_pt_ratio'],
    #  wdir+'/python/json/WP90PlusIso_BCDEF.json' : ['WP90PlusIso_BCDEF', 'eta_pt_ratio'],
    #trigg
    #  wdir+'/python/json/WP90_BCD_withRelIso.json' : ['electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCD', 'eta_pt_ratio'],
    #  wdir+'/python/json/WP90_BCDEF_withRelIso.json' : ['electronTriggerEfficiencyHLT_Ele27_WPLoose_eta2p1_WP90_BCDEF', 'eta_pt_ratio']
}

muIso_SF_reader = LeptonSF(wdir+'/python/json/MuonIso_Z_RunBCD_prompt80X_7p65.json', 'MC_NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1', 'abseta_pt_ratio')
muID_SF_reader = LeptonSF(wdir+'/python/json/MuonID_Z_RunBCD_prompt80X_7p65.json', 'MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1', 'abseta_pt_ratio')
muTrig_Eff_reader_1 = LeptonSF(wdir+'/python/json/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.json', 'IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093', 'abseta_pt_DATA')
muTrig_Eff_reader_2 = LeptonSF(wdir+'/python/json/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.json', 'IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097', 'abseta_pt_DATA')

#@@@@@@@@@@@@@@electronn sf@@@@@@@@@@@@@@@
eleID_SF_file = ROOT.TFile.Open('Data/egammaEffiAll.txt_SF2D.root', 'read')
eleID_SF_his = eleID_SF_file.Get('EGamma_SF2D_Medium')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#btag SF
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ROOT.gSystem.Load('libCondFormatsBTauObjects') 
ROOT.gSystem.Load('libCondToolsBTau')
bTagCalib = ROOT.BTagCalibration('CSVv2', 'Data/CSVv2_ichep.csv')
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


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Some functions
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def calBtagWeight(selJet_idx, tree, btagSF_reader_T, btagSF_reader_M):
  #FIXME use constant efficiency for jet
  eff_b = 0.475 #consistent with AN-16-251 page 28
  eff_c = 0.0315
  eff_l = 0.000726 #AN-16-251, tab 14 page 43, 80-120: 0.00147 +/- 0.00017 
  pMC = [[1,1,1], [1, 1, 1]] #central, up, down, tight, loose
  pData = [[1,1,1], [1, 1, 1]] #central, up, down, tight, loose
  bTagWeights = [[1,1,1], [1, 1, 1]] #central, up, down, tight, loose
  w = [[1,1,1], [1, 1, 1]] #central, up, down, tight, loose
  #has_lightJet_passCSVT = False
  for i in range(0, len(selJet_idx)):
    iJ = selJet_idx[i]
    btagTmp = tree.Jet_btagCSV[iJ]
    jetFlavTmp = abs(tree.Jet_hadronFlavour[iJ])
    jetFlav = 2
    eff = eff_l
    if jetFlavTmp == 5:
      jetFlav = 0
      eff = eff_b
    if jetFlavTmp == 4:
      jetFlav = 1
      eff = eff_c
    jetEta = tree.Jet_eta[iJ]
    jetPt = tree.Jet_pt[iJ]
    #get scale factor
    passT = (btagTmp > float(config.get('General', 'CSV_Tight')))
    passM = (btagTmp > float(config.get('General', 'CSV_Medium')))
    sf_Ts = []
    sf_Ts.append(btagSF_reader_T.eval_auto_bounds('central', jetFlav, jetEta, jetPt))
    sf_Ts.append(btagSF_reader_T.eval_auto_bounds('up', jetFlav, jetEta, jetPt))
    sf_Ts.append(btagSF_reader_T.eval_auto_bounds('down', jetFlav, jetEta, jetPt))
    sf_Ms = []
    sf_Ms.append(btagSF_reader_M.eval_auto_bounds('central', jetFlav, jetEta, jetPt))
    sf_Ms.append(btagSF_reader_M.eval_auto_bounds('up', jetFlav, jetEta, jetPt))
    sf_Ms.append(btagSF_reader_M.eval_auto_bounds('down', jetFlav, jetEta, jetPt))
    
    if passT:
     pMC[0][0] = pMC[0][0]*eff #T, central
     pMC[0][1] = pMC[0][1]*eff #T, up
     pMC[0][2] = pMC[0][2]*eff #T, down
     pData[0][0] = pData[0][0]*eff*sf_Ts[0] #T, central
     pData[0][1] = pData[0][1]*eff*sf_Ts[1] #T, up
     pData[0][2] = pData[0][2]*eff*sf_Ts[2] #T, down
     #if jetFlav == 2:
     #  has_lightJet_passCSVT = True
     #  print 'Light jet: ', jetPt, ' ', jetEta, ' ', eff, ' ', sf_Ts[0], sf_Ts[1], sf_Ts[2]
    else:
     pMC[0][0] = pMC[0][0]*(1-eff)
     pMC[0][1] = pMC[0][1]*(1-eff)
     pMC[0][2] = pMC[0][2]*(1-eff)
     pData[0][0] = pData[0][0]*(1-sf_Ts[0]*eff)
     pData[0][1] = pData[0][1]*(1-sf_Ts[1]*eff)
     pData[0][2] = pData[0][2]*(1-sf_Ts[2]*eff)

    if passM:
     pMC[1][0] = pMC[1][0]*eff #M, central
     pMC[1][1] = pMC[1][1]*eff #M, up
     pMC[1][2] = pMC[1][2]*eff #M, down
     pData[1][0] = pData[1][0]*eff*sf_Ms[0] #M, central
     pData[1][1] = pData[1][1]*eff*sf_Ms[1] #M, up
     pData[1][2] = pData[1][2]*eff*sf_Ms[2] #M, down
    else:
     pMC[1][0] = pMC[1][0]*(1-eff)
     pMC[1][1] = pMC[1][1]*(1-eff)
     pMC[1][2] = pMC[1][2]*(1-eff)
     pData[1][0] = pData[1][0]*(1-sf_Ms[0]*eff)
     pData[1][1] = pData[1][1]*(1-sf_Ms[1]*eff)
     pData[1][2] = pData[1][2]*(1-sf_Ms[2]*eff)

  for i in range(0,len(pMC)):
    for j in range(0,len(pMC[i])):
      w[i][j] = pData[i][j]/pMC[i][j]
  #if has_lightJet_passCSVT:
  #  print 'Final weight with light jet passing CSVT: ', w

  return w  
    
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


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Main functions
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


def fillTree(inputfile,outputfile):

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

    emu_lep_eta = array('f',2*[-10])
    emu_lep_eta[0], emu_lep_eta[1] = -10, -10
    newtree.Branch('emu_lep_eta',emu_lep_eta,'emu_lep_eta[2]/F') #=0 for muon, 1 for electron
    
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

    muweight = array('f', 3*[1])
    muweight[0], muweight[1], muweight[2] = 1., 1., 1.
    newtree.Branch('muweight',muweight,'muweight[3]/F')
    
    muweight_trig = array('f', 3*[1])
    muweight_trig[0], muweight_trig[1], muweight_trig[2] = 1., 1., 1.
    newtree.Branch('muweight_trig',muweight_trig,'muweight_trig[3]/F')
    
    emuweight = array('f', 3*[1])
    emuweight[0], emuweight[1], emuweight[2] = 1., 1., 1.
    newtree.Branch('emuweight',emuweight,'emuweight[3]/F')

    emuweight_trig = array('f', 3*[1])
    emuweight_trig[0], emuweight_trig[1], emuweight_trig[2] = 1., 1., 1.
    newtree.Branch('emuweight_trig',emuweight_trig,'emuweight_trig[3]/F')


    idxJet_passCSV = array('i',2*[-1])
    idxJet_passCSV[0] = -1  
    idxJet_passCSV[1] = -1 
    newtree.Branch('idxJet_passCSV',idxJet_passCSV,'idxJet_passCSV[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium

    idxJet_passCSV_SVT = array('i',2*[-1])
    idxJet_passCSV_SVT[0] = -1  
    idxJet_passCSV_SVT[1] = -1 
    newtree.Branch('idxJet_passCSV_SVT',idxJet_passCSV_SVT,'idxJet_passCSV_SVT[2]/I')  #index of highest pT jet passing CSV cut: first entry is for tight, second entry is for medium
    
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
    

    print 'starting event loop, processing',str(nEntries),'events'
    j_out=10000;

    #########################
    #Start event loop
    #########################
    #TEMP
    #nEntries = 5000
    
    nFilled = 0
    
    for entry in range(0,nEntries):
      # if entry>1000: break
      if ((entry%j_out)==0):
          if ((entry/j_out)==9 and j_out < 1e4): j_out*=10;
          print strftime("%Y-%m-%d %H:%M:%S", gmtime()),' - processing event',str(entry)+'/'+str(nEntries), '(cout every',j_out,'events)'
          sys.stdout.flush()
      
      #=====> Get entry
      tree.GetEntry(entry)
      
      
      if job.type != 'DATA':
          EventForTraining[0]=int(not TFlag.EvalInstance())
     
      #find emu event
      #initial emu event
      is_emu[0] = 0
      emu_lep_pt[0], emu_lep_pt[1] = -1, -1
      emu_lep_eta[0], emu_lep_eta[1] = -10, -10
      emu_lep_pdgId[0], emu_lep_pdgId[1] = -1, -1
      emu_lep_looseId[0], emu_lep_looseId[1] = -1, -1
      emu_lep_mediumId[0], emu_lep_mediumId[1] = -1, -1
      emu_lep_tightId[0], emu_lep_tightId[1] = -1, -1
      emu_lep_iso[0], emu_lep_iso[1] = -1, -1
      
      if tree.Vtype != 0 and tree.Vtype != 1:
        lepTmp_idx = []
        lepTmp_pt = []
        lepTmp_eta = []
        lepTmp_pdgId = []
        lepTmp_looseId = []
        lepTmp_mediumId = []
        lepTmp_tightId = []
        lepTmp_iso = []
        iLepTmp = -1
        for i in range(tree.nvLeptons):
          iLepTmp = iLepTmp + 1
          lepTmp_idx.append(iLepTmp)
          lepTmp_pt.append(tree.vLeptons_pt[i])
          lepTmp_eta.append(tree.vLeptons_eta[i])
          lepTmp_pdgId.append(tree.vLeptons_pdgId[i])
          lepTmp_looseId.append(tree.vLeptons_looseIdPOG[i])
          lepTmp_mediumId.append(tree.vLeptons_mediumIdPOG_ICHEP2016[i])
          lepTmp_tightId.append(tree.vLeptons_tightId[i])
          lepTmp_iso.append(tree.vLeptons_pfRelIso04[i])
        for i in range(tree.naLeptons):
          iLepTmp = iLepTmp + 1
          lepTmp_idx.append(iLepTmp)
          lepTmp_pt.append(tree.aLeptons_pt[i])
          lepTmp_eta.append(tree.aLeptons_eta[i])
          lepTmp_pdgId.append(tree.aLeptons_pdgId[i])
          lepTmp_looseId.append(tree.aLeptons_looseIdPOG[i])
          lepTmp_mediumId.append(tree.aLeptons_mediumIdPOG_ICHEP2016[i])
          lepTmp_tightId.append(tree.aLeptons_tightId[i])
          lepTmp_iso.append(tree.aLeptons_pfRelIso04[i])
        
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
             emu_lep_iso[0], emu_lep_iso[1] = lepTmp_iso[muIdx], lepTmp_iso[eIdx]
             emu_lep_pdgId[0], emu_lep_pdgId[1] = lepTmp_pdgId[muIdx], lepTmp_pdgId[eIdx]
             emu_lep_looseId[0], emu_lep_looseId[1] = lepTmp_looseId[muIdx], lepTmp_looseId[eIdx]
             emu_lep_mediumId[0], emu_lep_mediumId[1] = lepTmp_mediumId[muIdx], lepTmp_mediumId[eIdx]
             emu_lep_tightId[0], emu_lep_tightId[1] = lepTmp_tightId[muIdx], lepTmp_tightId[eIdx]
              
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
      for i in range(0, tree.nJet):
          if tree.Jet_pt[i] < float(config.get('General', 'jetPt')) or abs(tree.Jet_eta[i]) > float(config.get('General', 'absJetEta')): continue
          selJet_idx.append(i)
          if tree.Jet_btagCSV[i] >= 0:
              idxJet_sortedCSV[nCSVjet[0]] = i
              nCSVjet[0] += 1
      
      #sort jet according to pt
      util_funcs.sortPt(selJet_idx, tree.Jet_pt)

      idxJet_passCSV[0] = -1
      idxJet_passCSV[1] = -1
      idxJet_passCSV_SVT[0] = -1
      idxJet_passCSV_SVT[1] = -1
      
      for iJ in selJet_idx:
        if tree.Jet_btagCSV[iJ] >= float(config.get('General', 'CSV_Tight')):
          if idxJet_passCSV[0] == -1: idxJet_passCSV[0] = iJ
          if tree.Jet_vtxMass[iJ] > 0 and idxJet_passCSV_SVT[0] == -1:
            idxJet_passCSV_SVT[0] = iJ
        if tree.Jet_btagCSV[iJ] >= float(config.get('General', 'CSV_Medium')):
          if idxJet_passCSV[1] == -1: idxJet_passCSV[1] = iJ
          if tree.Jet_vtxMass[iJ] > 0 and idxJet_passCSV_SVT[1] == -1:
            idxJet_passCSV_SVT[1] = iJ


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

      #bTag scale factor
      bTagWeights = [[1,1,1], [1,1,1]] 
      if not tree.isData:
        bTagWeights = calBtagWeight(selJet_idx, tree, btagSF_reader_T, btagSF_reader_M)
      for iTmp in range(0, 3):
        bTagWeight_CSVT[iTmp] = bTagWeights[0][iTmp]
        bTagWeight_CSVM[iTmp] = bTagWeights[1][iTmp]
       

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
      #Muon in heppy use looseID  
      muweight_trig[0],  muweight_trig[1],  muweight_trig[2] = 1., 1., 1.
      #id and iso weight for emu sample
      emuweight[0], emuweight[1], emuweight[2] = 1., 1., 1.
      emuweight_trig[0], emuweight_trig[1], emuweight_trig[2] = 1., 1., 1.

      if job.type != 'DATA':
 
        if tree.Vtype == 0:
          #recalculate SF for muon iso
          weight_SF_LooseISO = [1, 1, 1] #central, up, down
          weight = []
          weight.append(muIso_SF_reader.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0]))
          weight.append(muIso_SF_reader.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
          computeSF(weight_SF_LooseISO, weight)
          if debug:
            print "================"
            print "mu1 pt, eta: ", tree.vLeptons_pt[0], ' ', tree.vLeptons_eta[0]
            print "mu2 pt, eta: ", tree.vLeptons_pt[1], ' ', tree.vLeptons_eta[1]
            print "Iso weights: ", weight
            print 'weight_SF_LooseIso: ', weight_SF_LooseISO

          muweight[0] = tree.vLeptons_SF_IdCutLoose[0]*tree.vLeptons_SF_IdCutLoose[1]*weight_SF_LooseISO[0]*tree.vLeptons_SF_trk_eta[0]*tree.vLeptons_SF_trk_eta[1]
          muweight[1] = (tree.vLeptons_SF_IdCutLoose[0]-tree.vLeptons_SFerr_IdCutLoose[0])*(tree.vLeptons_SF_IdCutLoose[1]-tree.vLeptons_SFerr_IdCutLoose[1])*weight_SF_LooseISO[1]*(tree.vLeptons_SF_trk_eta[0]-tree.vLeptons_SFerr_trk_eta[0])*(tree.vLeptons_SF_trk_eta[1]-tree.vLeptons_SFerr_trk_eta[1])
          muweight[2] = (tree.vLeptons_SF_IdCutLoose[0]+tree.vLeptons_SFerr_IdCutLoose[0])*(tree.vLeptons_SF_IdCutLoose[1]+tree.vLeptons_SFerr_IdCutLoose[1])*weight_SF_LooseISO[2]*(tree.vLeptons_SF_trk_eta[0]+tree.vLeptons_SFerr_trk_eta[0])*(tree.vLeptons_SF_trk_eta[1]+tree.vLeptons_SFerr_trk_eta[1])   

          weight = []
          weight.append(muTrig_Eff_reader_2.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0]))
          weight.append(muTrig_Eff_reader_2.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
          computeTrigWeight(muweight_trig,weight)
          
          if debug:
            print "Trigger eff: ", weight
            print 'muweight_trig: ', muweight_trig

        
        #@@@@@@@@@@@@@@@@@@@@@@@e-mu weight@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #e-mu weight = ID+iso+trk for muon and ID+iso+trk for electron

        if is_emu[0] == 1:
          sfTmp = (muIso_SF_reader.get_2D(emu_lep_pt[0], emu_lep_eta[0]))[0]
          if sfTmp > 0: emuweight[0] = emuweight[0]*sfTmp
          sfTmp = (muID_SF_reader.get_2D(emu_lep_pt[0], emu_lep_eta[0]))[0]
          if sfTmp > 0: emuweight[0] = emuweight[0]*sfTmp
          #now doing for electron, only ID SF available
          iBin = eleID_SF_his.FindFixBin(emu_lep_eta[1], emu_lep_pt[1])
          sfTmp = eleID_SF_his.GetBinContent(iBin)
          if sfTmp > 0: emuweight[0] = emuweight[0]*sfTmp
          #trigger weight
          sfTmp = (muTrig_Eff_reader_2.get_2D(emu_lep_pt[0], emu_lep_eta[0]))[0]
          if sfTmp > 0: emuweight_trig[0] = sfTmp

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

    print 'Exit loop. Total filled event: ', nFilled
    newtree.FlushBaskets()
    newtree.AutoSave()
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
      tmpD =  config.get('Directories','samplefiles').replace('/','') + '_afterPrepStep/'
      os.system('mkdir ' + tmpD)
      tmpFileList =  tmpD + job.prefix + job.identifier + '.txt'
      os.system('rm -f ' + tmpFileList)
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
    print '>>>>> nFile_for_processing: ', nFile_for_processing
    nFile = 0
    for line in open(tmpFileList).readlines():
      line = line.replace("\n","")
      line_tmp = line.split('/')
      outTmp = line_tmp[len(line_tmp)-1].replace('.root', '_moreInfo.root')
      if run_locally == 'True': outTmp = outputFolder + '/' + outTmp
      input_processings.append((line,outTmp))
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
