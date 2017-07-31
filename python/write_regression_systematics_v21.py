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

from myutils import BetterConfigParser, ParseInfo, TreeCache, LeptonSF, util_funcs

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

applyBTagweights=eval(config.get('Regression','applyBTagweights'))
print 'applyBTagweights is', applyBTagweights
csv_rwt_hf=config.get('BTagHFweights','file')
csv_rwt_lf=config.get('BTagLFweights','file')
applyRegression=eval(config.get('Regression','applyRegression'))
print 'applyRegression is', applyRegression
print "csv_rwt_hf",csv_rwt_hf,"csv_rwt_lf",csv_rwt_lf
bweightcalc = BTagWeightCalculator(
    csv_rwt_hf,
    csv_rwt_lf
)
bweightcalc.btag = "btag"

#lhe_weight_map = False if not config.has_option('LHEWeights', 'weights_per_bin') else eval(config.get('LHEWeights', 'weights_per_bin'))


namelist=opts.names.split(',')

#load info
info = ParseInfo(samplesinfo,pathIN)

def isInside(map_,eta,phi):
    bin_ = map_.FindBin(phi,eta)
    bit = map_.GetBinContent(bin_)
    if bit>0:
        return True
    else:
        return False

#def deltaPhi(phi1, phi2): 
#    result = phi1 - phi2
#    while (result > math.pi): result -= 2*math.pi
#    while (result <= -math.pi): result += 2*math.pi
#    return result

#def deltaR(phi1, eta1, phi2, eta2):
#    deta = eta1-eta2
#    dphi = deltaPhi(phi1, phi2)
#    result = math.sqrt(deta*deta + dphi*dphi)
#    return result

def addAdditionalJets(H, tree):
    for i in range(tree.nhjidxaddJetsdR08):
        idx = tree.hjidxaddJetsdR08[i]
        if (idx == tree.hJCidx[0]) or (idx == tree.hJCidx[1]): continue
        addjet = ROOT.TLorentzVector()
        if idx<tree.nJet:
		addjet.SetPtEtaPhiM(tree.Jet_pt[idx],tree.Jet_eta[idx],tree.Jet_phi[idx],tree.Jet_mass[idx])
        H = H + addjet
    return H

def resolutionBias(eta):
    if(eta< 0.5): return 0.052
    if(eta< 1.1): return 0.057
    if(eta< 1.7): return 0.096
    if(eta< 2.3): return 0.134
    if(eta< 5): return 0.28
    return 0

def corrPt(pt,eta,mcPt):
    return 1 ##FIXME
    # return (pt+resolutionBias(math.fabs(eta))*(pt-mcPt))/pt

def corrCSV(btag,  csv, flav):
    if(csv < 0.): return csv
    if(csv > 1.): return csv;
    if(flav == 0): return csv;
    if(math.fabs(flav) == 5): return  btag.ib.Eval(csv)
    if(math.fabs(flav) == 4): return  btag.ic.Eval(csv)
    if(math.fabs(flav) != 4  and math.fabs(flav) != 5): return  btag.il.Eval(csv)
    return -10000


def csvReshape(sh, pt, eta, csv, flav):
    return sh.reshape(float(eta), float(pt), float(csv), int(flav))

#=====> main program start here

ROOT.gROOT.ProcessLine( #not in use currently
        "struct H {\
        int         HiggsFlag;\
        float         mass;\
        float         pt;\
        float         eta;\
        float         phi;\
        float         dR;\
        float         dPhi;\
        float         dEta;\
        } ;"
    )

def findFatJetAndSubjet(fatjet_eta, fatjet_phi, subjet_eta, subjet_phi, subjet_btag, v_eta, v_phi, dRcut=1.5):
  
  
  #select fat jet not overlap with V
  fatjetIdxs = []
  for i in range(0, len(fatjet_eta)):
    if util_funcs.deltaR(fatjet_eta[i], fatjet_phi[i], v_eta, v_phi) > dRcut: fatjetIdxs.append(i)
 
  #select fat jet with at least two matched subjets
  fatjetIdx_twoSubjet = [] #fat jet with at lest two subjet
  subjetIdx_match_sortSubjetBtag = [] # index of two subjet with highest btag that match with fatjet
  for i in range(0, len(fatjetIdxs)):
    idxTmp = fatjetIdxs[i]
    subjetIdxs = []
    for j in range(0, len(subjet_eta)):
      if util_funcs.deltaR(fatjet_eta[idxTmp], fatjet_phi[idxTmp], subjet_eta[j], subjet_phi[j]) < dRcut: subjetIdxs.append(j)
    if len(subjetIdxs) >= 2:
      util_funcs.sortPt(subjetIdxs, subjet_btag)
      fatjetIdx_twoSubjet.append(idxTmp)
      subjetIdx_match_sortSubjetBtag.append([subjetIdxs[0], subjetIdxs[1]])
   
  #finally sort fat jet according to minimum btag of two subjets
  for i in range(0, len(fatjetIdx_twoSubjet)):
    minBtag1 = min(subjet_btag[subjetIdx_match_sortSubjetBtag[i][0]], subjet_btag[subjetIdx_match_sortSubjetBtag[i][1]])
    for j in range(i + 1, len(fatjetIdx_twoSubjet)):
      minBtag2 = min(subjet_btag[subjetIdx_match_sortSubjetBtag[j][0]], subjet_btag[subjetIdx_match_sortSubjetBtag[j][1]])
      if minBtag1 < minBtag2:
        tmp = fatjetIdx_twoSubjet[i]
        fatjetIdx_twoSubjet[i] = fatjetIdx_twoSubjet[j]
        fatjetIdx_twoSubjet[j] = tmp
        subjetTmp = subjetIdx_match_sortSubjetBtag[i]
        subjetIdx_match_sortSubjetBtag[i] = subjetIdx_match_sortSubjetBtag[j]
        subjetIdx_match_sortSubjetBtag[j] = subjetTmp
  if len(fatjetIdx_twoSubjet) >= 1: return [fatjetIdx_twoSubjet[0],subjetIdx_match_sortSubjetBtag[0][0],subjetIdx_match_sortSubjetBtag[0][1]] 
  else: return []
      
def findFatJetTest(njet, fatjet_eta, fatjet_phi, v_eta, v_phi, dRcut=1.5):
  if njet != len(fatjet_eta): print 'Warning: nfatjet and jet_eta array not equal: ', njet, len(fatjet_eta)
  #select fat jet not overlap with V
  fatjetIdxs = []
  for i in range(0, len(fatjet_eta)):
    if util_funcs.deltaR(fatjet_eta[i], fatjet_phi[i], v_eta, v_phi) > dRcut: fatjetIdxs.append(i)
  return fatjetIdxs 




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

    # Add bTag weights
    if applyBTagweights:
        bTagWeight_new = array('f',[0])
        bTagWeightJESUp_new = array('f',[0])
        bTagWeightJESDown_new = array('f',[0])
        bTagWeightLFUp_new = array('f',[0])
        bTagWeightLFDown_new = array('f',[0])
        bTagWeightHFUp_new = array('f',[0])
        bTagWeightHFDown_new = array('f',[0])
        bTagWeightLFStats1Up_new = array('f',[0])
        bTagWeightLFStats1Down_new = array('f',[0])
        bTagWeightLFStats2Up_new = array('f',[0])
        bTagWeightLFStats2Down_new = array('f',[0])
        bTagWeightHFStats1Up_new = array('f',[0])
        bTagWeightHFStats1Down_new = array('f',[0])
        bTagWeightHFStats2Up_new = array('f',[0])
        bTagWeightHFStats2Down_new = array('f',[0])
        bTagWeightcErr1Up_new = array('f',[0])
        bTagWeightcErr1Down_new = array('f',[0])
        bTagWeightcErr2Up_new = array('f',[0])
        bTagWeightcErr2Down_new = array('f',[0])


        bTagWeight_new[0] = 1
        bTagWeightJESUp_new[0] = 1
        bTagWeightJESDown_new[0] = 1
        bTagWeightLFUp_new[0] = 1
        bTagWeightLFDown_new[0] = 1
        bTagWeightHFUp_new[0] = 1
        bTagWeightHFDown_new[0] = 1
        bTagWeightLFStats1Up_new[0] = 1
        bTagWeightLFStats1Down_new[0] = 1
        bTagWeightLFStats2Up_new[0] = 1
        bTagWeightLFStats2Down_new[0] = 1
        bTagWeightHFStats1Up_new[0] = 1
        bTagWeightHFStats1Down_new[0] = 1
        bTagWeightHFStats2Up_new[0] = 1
        bTagWeightHFStats2Down_new[0] = 1
        bTagWeightcErr1Up_new[0] = 1
        bTagWeightcErr1Down_new[0] = 1
        bTagWeightcErr2Up_new[0] = 1
        bTagWeightcErr2Down_new[0] = 1

        newtree.Branch('bTagWeight_new',bTagWeight_new,'bTagWeight_new/F')
        newtree.Branch('bTagWeightJESUp_new',bTagWeightJESUp_new,'bTagWeightJESUp_new/F')
        newtree.Branch('bTagWeightJESDown_new',bTagWeightJESDown_new,'bTagWeightJESDown_new/F')
        newtree.Branch('bTagWeightLFUp_new',bTagWeightLFUp_new,'bTagWeightLFUp_new/F')
        newtree.Branch('bTagWeightLFDown_new',bTagWeightLFDown_new,'bTagWeightLFDown_new/F')
        newtree.Branch('bTagWeightHFUp_new',bTagWeightHFUp_new,'bTagWeightHFUp_new/F')
        newtree.Branch('bTagWeightHFDown_new',bTagWeightHFDown_new,'bTagWeightHFDown_new/F')
        newtree.Branch('bTagWeightLFStats1Up_new',bTagWeightLFStats1Up_new,'bTagWeightLFStats1Up_new/F')
        newtree.Branch('bTagWeightLFStats1Down_new',bTagWeightLFStats1Down_new,'bTagWeightLFStats1Down_new/F')
        newtree.Branch('bTagWeightLFStats2Up_new',bTagWeightLFStats2Up_new,'bTagWeightLFStats2Up_new/F')
        newtree.Branch('bTagWeightLFStats2Down_new',bTagWeightLFStats2Down_new,'bTagWeightLFStats2Down_new/F')
        newtree.Branch('bTagWeightHFStats1Up_new',bTagWeightHFStats1Up_new,'bTagWeightHFStats1Up_new/F')
        newtree.Branch('bTagWeightHFStats1Down_new',bTagWeightHFStats1Down_new,'bTagWeightHFStats1Down_new/F')
        newtree.Branch('bTagWeightHFStats2Up_new',bTagWeightHFStats2Up_new,'bTagWeightHFStats2Up_new/F')
        newtree.Branch('bTagWeightHFStats2Down_new',bTagWeightHFStats2Down_new,'bTagWeightHFStats2Down_new/F')
        newtree.Branch('bTagWeightcErr1Up_new',bTagWeightcErr1Up_new,'bTagWeightcErr1Up_new/F')
        newtree.Branch('bTagWeightcErr1Down_new',bTagWeightcErr1Down_new,'bTagWeightcErr1Down_new/F')
        newtree.Branch('bTagWeightcErr2Up_new',bTagWeightcErr2Up_new,'bTagWeightcErr2Up_new/F')
        newtree.Branch('bTagWeightcErr2Down_new',bTagWeightcErr2Down_new,'bTagWeightcErr2Down_new/F')


        bTagWeightSubjetCA15pruned = array('f',[0])
        bTagWeightSubjetAK08softdrop = array('f',[0])
       
        bTagWeightSubjetCA15pruned[0] = 1
        bTagWeightSubjetAK08softdrop[0] = 1

        newtree.Branch('bTagWeightSubjetCA15pruned',bTagWeightSubjetCA15pruned,'bTagWeightSubjetCA15pruned/F')
        newtree.Branch('bTagWeightSubjetAK08softdrop',bTagWeightSubjetAK08softdrop,'bTagWeightSubjetAK08softdrop/F')





    if channel == "Zmm" or channel == 'Zee':
    #Special weights

        DY_specialWeight= array('f',[0])
        DY_specialWeight[0] = 1
        newtree.Branch('DY_specialWeight',DY_specialWeight,'DY_specialWeight/F')


    #Add reg VHDphi
        HVdPhi_reg = array('f',[0])
        HVdPhi_reg[0] = 300
        newtree.Branch('HVdPhi_reg',HVdPhi_reg,'HVdPhi_reg/F')

    # Add muon SF
        vLeptons_SFweight_HLT = array('f',[0])
        lepton_EvtWeight = array('f',[0])

        # vLeptons_SFweight_IdLoose[0] = 1
        # vLeptons_SFweight_IsoLoose[0] = 1
        vLeptons_SFweight_HLT[0] = 1 
        lepton_EvtWeight[0] = 1

        # newtree.Branch('vLeptons_SFweight_IdLoose',vLeptons_SFweight_IdLoose,'vLeptons_SFweight_IdLoose/F')
        # newtree.Branch('vLeptons_SFweight_IsoLoose',vLeptons_SFweight_IsoLoose,'vLeptons_SFweight_IsoLoose/F')
        newtree.Branch('vLeptons_SFweight_HLT',vLeptons_SFweight_HLT,'vLeptons_SFweight_HLT/F') #=1 if Vtype == 1, electron channel, applied for muon
        newtree.Branch('lepton_EvtWeight',lepton_EvtWeight,'lepton_EvtWeight/F') #=1 for Vtype == 0, muon channel, applied for electron ID and trigger (HLT?)

        
    ### Adding new variable from configuration ###
    newVariableNames = []
    try:
        writeNewVariables = eval(config.get("Regression","writeNewVariables"))

        ## remove MC variables in data ##
        if job.type == 'DATA':
            for idx in dict(writeNewVariables):
                formula = writeNewVariables[idx]
                if 'gen' in formula or 'Gen' in formula or 'True' in formula or 'true' in formula or 'mc' in formula or 'Mc' in formula:
                    print "Removing: ",idx," with ",formula
                    del writeNewVariables[idx]

        newVariableNames = writeNewVariables.keys()
        newVariables = {}
        newVariableFormulas = {}
        for variableName in newVariableNames:
            formula = writeNewVariables[variableName]
            newVariables[variableName] = array('f',[0])
            newtree.Branch(variableName,newVariables[variableName],variableName+'/F')
            newVariableFormulas[variableName] =ROOT.TTreeFormula(variableName,formula,tree)
            print "adding variable ",variableName,", using formula",writeNewVariables[variableName]," ."
    except:
        pass

    lepCorrs = {}
    wdir = config.get('Directories','vhbbpath')
    jsons = {
      wdir+'/python/json/ScaleFactor_GsfElectronToRECO_passingTrigWP80.json' : ['ScaleFactor_GsfElectronToRECO_passingTrigWP80', 'eta_pt_ratio'],
                #Not available for the moment
      wdir+'/python/json/ScaleFactor_HLT_Ele23_WPLoose_Gsf_v.json' : ['ScaleFactor_HLT_Ele23_WPLoose_Gsf_v', 'eta_pt_ratio'],
    }
    
    for j, name in jsons.iteritems():
      lepCorr = LeptonSF(j , name[0], name[1])
      lepCorrs[j] = lepCorr
    
    print lepCorrs

#@@@@@@@@@@@@@@@@@@@initialization to fill Fatjet that not overlap with V@@@@@@@@@@@@@@@@@@@@@@@@@@@
    MAX_FATJET = 10
    nFatjetCA15ungroomed_notV = array('i',[0])
    nFatjetCA15pruned_notV= array('i',[0])
    nFatjetCA15trimmed_notV= array('i',[0])
    nFatjetCA15softdrop_notV= array('i',[0])
    nFatjetCA15softdropz2b1_notV= array('i',[0])
    nFatjetAK08ungroomed_notV= array('i',[0])
    idx_FatjetCA15ungroomed_notV = array('i',MAX_FATJET*[-1])
    idx_FatjetCA15pruned_notV= array('i',MAX_FATJET*[-1])
    idx_FatjetCA15trimmed_notV = array('i',MAX_FATJET*[-1])
    idx_FatjetCA15softdrop_notV = array('i',MAX_FATJET*[-1])
    idx_FatjetCA15softdropz2b1_notV = array('i',MAX_FATJET*[-1])
    idx_FatjetAK08ungroomed_notV = array('i',MAX_FATJET*[-1])
    newtree.Branch('nFatjetCA15ungroomed_notV', nFatjetCA15ungroomed_notV, 'nFatjetCA15ungroomed_notV/I')
    newtree.Branch('idx_FatjetCA15ungroomed_notV', idx_FatjetCA15ungroomed_notV, 'idx_FatjetCA15ungroomed_notV[nFatjetCA15ungroomed_notV]/I')
    newtree.Branch('nFatjetCA15pruned_notV', nFatjetCA15pruned_notV, 'nFatjetCA15pruned_notV/I')
    newtree.Branch('idx_FatjetCA15pruned_notV', idx_FatjetCA15pruned_notV, 'idx_FatjetCA15pruned_notV[nFatjetCA15pruned_notV]/I')
    newtree.Branch('nFatjetCA15trimmed_notV', nFatjetCA15trimmed_notV, 'nFatjetCA15trimmed_notV/I')
    newtree.Branch('idx_FatjetCA15trimmed_notV', idx_FatjetCA15trimmed_notV, 'idx_FatjetCA15trimmed_notV[nFatjetCA15trimmed_notV]/I')
    newtree.Branch('nFatjetCA15softdrop_notV', nFatjetCA15softdrop_notV, 'nFatjetCA15softdrop_notV/I')
    newtree.Branch('idx_FatjetCA15softdrop_notV', idx_FatjetCA15softdrop_notV, 'idx_FatjetCA15softdrop_notV[nFatjetCA15softdrop_notV]/I')
    newtree.Branch('nFatjetCA15softdropz2b1_notV', nFatjetCA15softdropz2b1_notV, 'nFatjetCA15softdropz2b1_notV/I')
    newtree.Branch('idx_FatjetCA15softdropz2b1_notV', idx_FatjetCA15softdropz2b1_notV, 'idx_FatjetCA15softdropz2b1_notV[nFatjetCA15softdropz2b1_notV]/I')
    newtree.Branch('nFatjetAK08ungroomed_notV', nFatjetAK08ungroomed_notV, 'nFatjetAK08ungroomed_notV/I')
    newtree.Branch('idx_FatjetAK08ungroomed_notV', idx_FatjetAK08ungroomed_notV, 'idx_FatjetAK08ungroomed_notV[nFatjetAK08ungroomed_notV]/I')

    
#@@@@@@@@@@@@@@@@@@@initialization to fill subjects matched with leading fatjet@@@@@@@@@@@@@@@@@@@@@@@@@@@
    MAX_NSUBJET = 10
    nSubjet_matchCA15ungroomed = array('i',[0])
    idx_subjet_matchCA15ungroomed = array('i',MAX_NSUBJET*[-1])
    nSubjet_matchCA15pruned = array('i',[0])
    idx_subjet_matchCA15pruned = array('i',MAX_NSUBJET*[-1])
    nSubjet_matchCA15trimmed = array('i',[0])
    idx_subjet_matchCA15trimmed = array('i',MAX_NSUBJET*[-1])
    nSubjet_matchCA15softdrop = array('i',[0])
    idx_subjet_matchCA15softdrop = array('i',MAX_NSUBJET*[-1])
    nSubjet_matchCA15softdropz2b1 = array('i',[0])
    idx_subjet_matchCA15softdropz2b1 = array('i',MAX_NSUBJET*[-1])
    nSubjet_matchAK08ungroomed = array('i',[0])
    idx_subjet_matchAK08ungroomed = array('i',MAX_NSUBJET*[-1])
    
    newtree.Branch('nSubjet_matchCA15ungroomed', nSubjet_matchCA15ungroomed, 'nSubjet_matchCA15ungroomed/I')
    newtree.Branch('idx_subjet_matchCA15ungroomed', idx_subjet_matchCA15ungroomed, 'idx_subjet_matchCA15ungroomed[nSubjet_matchCA15ungroomed]/I')
    newtree.Branch('nSubjet_matchCA15pruned', nSubjet_matchCA15pruned, 'nSubjet_matchCA15pruned/I')
    newtree.Branch('idx_subjet_matchCA15pruned', idx_subjet_matchCA15pruned, 'idx_subjet_matchCA15pruned[nSubjet_matchCA15pruned]/I')
    newtree.Branch('nSubjet_matchCA15trimmed', nSubjet_matchCA15trimmed, 'nSubjet_matchCA15trimmed/I')
    newtree.Branch('idx_subjet_matchCA15trimmed', idx_subjet_matchCA15trimmed, 'idx_subjet_matchCA15trimmed[nSubjet_matchCA15trimmed]/I')
    newtree.Branch('nSubjet_matchCA15softdrop', nSubjet_matchCA15softdrop, 'nSubjet_matchCA15softdrop/I')
    newtree.Branch('idx_subjet_matchCA15softdrop', idx_subjet_matchCA15softdrop, 'idx_subjet_matchCA15softdrop[nSubjet_matchCA15softdrop]/I')
    newtree.Branch('nSubjet_matchCA15softdropz2b1', nSubjet_matchCA15softdropz2b1, 'nSubjet_matchCA15softdropz2b1/I')
    newtree.Branch('idx_subjet_matchCA15softdropz2b1', idx_subjet_matchCA15softdropz2b1, 'idx_subjet_matchCA15softdropz2b1[nSubjet_matchCA15softdropz2b1]/I')
    newtree.Branch('nSubjet_matchAK08ungroomed', nSubjet_matchAK08ungroomed, 'nSubjet_matchAK08ungroomed/I')
    newtree.Branch('idx_subjet_matchAK08ungroomed', idx_subjet_matchAK08ungroomed, 'idx_subjet_matchAK08ungroomed[nSubjet_matchAK08ungroomed]/I')
   
    

    print 'starting event loop, processing',str(nEntries),'events'
    j_out=10000;

    #########################
    #Start event loop
    #########################
    #TEMP
    #nEntries = 1000
    nFilled = 0
    for entry in range(0,nEntries):
      # if entry>1000: break
      if ((entry%j_out)==0):
          if ((entry/j_out)==9 and j_out < 1e4): j_out*=10;
          print strftime("%Y-%m-%d %H:%M:%S", gmtime()),' - processing event',str(entry)+'/'+str(nEntries), '(cout every',j_out,'events)'
          sys.stdout.flush()
      
      #=====> Get entry
      tree.GetEntry(entry)
      
      ### Fill new variable from configuration ###
      for variableName in newVariableNames:
          newVariableFormulas[variableName].GetNdata()
          newVariables[variableName][0] = newVariableFormulas[variableName].EvalInstance()
      
      if job.type != 'DATA':
          EventForTraining[0]=int(not TFlag.EvalInstance())

      ##########################
      # Loop to fill bTag weights variables
      ##########################
      if applyBTagweights:
          if not job.type == 'DATA':
              jetsForBtagWeight = []
              for i in range(tree.nJet):
                  jetsForBtagWeight.append(Jet(tree.Jet_pt[i], tree.Jet_eta[i], tree.Jet_hadronFlavour[i], tree.Jet_btagCSV[i], bweightcalc.btag))

              bTagWeight_new[0] = bweightcalc.calcEventWeight(
                  jetsForBtagWeight, kind="final", systematic="nominal",
              )
              weights = {}
              for syst in ["JES", "LF", "HF", "LFStats1", "LFStats2", "HFStats1", "HFStats2", "cErr1", "cErr2"]:
                  for sdir in ["Up", "Down"]:
                      weights[syst+sdir] = bweightcalc.calcEventWeight(
                          jetsForBtagWeight, kind="final", systematic=syst+sdir
                          )
              bTagWeightJESUp_new[0] = weights["JESUp"]
              bTagWeightJESDown_new[0] = weights["JESDown"]
              bTagWeightLFUp_new[0] = weights["LFUp"]
              bTagWeightLFDown_new[0] = weights["LFDown"]
              bTagWeightHFUp_new[0] = weights["HFUp"]
              bTagWeightHFDown_new[0] = weights["HFDown"]
              bTagWeightLFStats1Up_new[0] = weights["LFStats1Up"]
              bTagWeightLFStats1Down_new[0] = weights["LFStats1Down"]
              bTagWeightLFStats2Up_new[0] = weights["LFStats2Up"]
              bTagWeightLFStats2Down_new[0] = weights["LFStats2Down"]
              bTagWeightHFStats1Up_new[0] = weights["HFStats1Up"]
              bTagWeightHFStats1Down_new[0] = weights["HFStats1Down"]
              bTagWeightHFStats2Up_new[0] = weights["HFStats2Up"]
              bTagWeightHFStats2Down_new[0] = weights["HFStats2Down"]
              bTagWeightcErr1Up_new[0] = weights["cErr1Up"]
              bTagWeightcErr1Down_new[0] = weights["cErr1Down"]
              bTagWeightcErr2Up_new[0] = weights["cErr2Up"]
              bTagWeightcErr2Down_new[0] = weights["cErr2Down"]

              #@@@@@btag weight for subjet@@@@@@@@@@@@

              bTagWeightSubjetCA15pruned[0] = 1
              bTagWeightSubjetAK08softdrop[0] = 1
              #match subjet with GenJet to identify subjet flavour
              SubjetCA15pruned_hadronFlavour = []
              for iSj in range(0, tree.nSubjetCA15pruned):
                dRmin = 999
                iGj_match = -1
                for iGj in range(0, tree.nGenJet):
                  dRtmp = util_funcs.deltaR(tree.SubjetCA15pruned_eta[iSj],tree.SubjetCA15pruned_phi[iSj], tree.GenJet_eta[iGj], tree.GenJet_phi[iGj])
                  if dRtmp < dRmin:
                    dRmin = dRtmp
                    iGj_match = iGj
                
                if dRmin < 0.3:
                  flavour = 0
                  if tree.GenJet_numBHadrons[iGj_match] >= 1: flavour = 5
                  if tree.GenJet_numBHadrons[iGj_match] == 0 and tree.GenJet_numCHadrons[iGj_match] >=1: flavour = 4
                  SubjetCA15pruned_hadronFlavour.append(flavour)

                else: SubjetCA15pruned_hadronFlavour.append(0)
              
              if tree.nSubjetCA15pruned != len(SubjetCA15pruned_hadronFlavour): print 'Warning: not equal subjet number and hadronFlavour', tree.nSubjetCA15pruned, len(SubjetCA15pruned_hadronFlavour) 

              #find subjets match with Fatjet that is not V
              #idx_subjetTmp = []
              subjet_ptTmp = []
              subjet_etaTmp = []
              subjet_flavourTmp = []
              subjet_btagTmp = []
              for i in range(0,tree.nFatjetCA15pruned):
                  if util_funcs.deltaR(tree.V_eta,tree.V_phi,tree.FatjetCA15pruned_eta[i],tree.FatjetCA15pruned_phi[i]) > 1.5:
                      for j in range(0,tree.nSubjetCA15pruned):
                          if util_funcs.deltaR(tree.SubjetCA15pruned_eta[j],tree.SubjetCA15pruned_phi[j],tree.FatjetCA15pruned_eta[i],tree.FatjetCA15pruned_phi[i]) < 1.5:
                              #idx_subjetTmp.append(j)
                              subjet_ptTmp.append(tree.SubjetCA15pruned_pt[j])
                              subjet_etaTmp.append(tree.SubjetCA15pruned_eta[j])
                              subjet_flavourTmp.append(SubjetCA15pruned_hadronFlavour[j])
                              subjet_btagTmp.append(tree.SubjetCA15pruned_btag[j])

              
              jetsForBtagWeight = []
              for i in range(0,len(subjet_ptTmp)):
                  jetsForBtagWeight.append(Jet(subjet_ptTmp[i], subjet_etaTmp[i], subjet_flavourTmp[i], subjet_btagTmp[i], bweightcalc.btag)) #CHECK what is calc.btag?
              #CHECK do we need to reset bweightcalc?
              bTagWeightSubjetCA15pruned[0] = bweightcalc.calcEventWeight(
                  jetsForBtagWeight, kind="final", systematic="nominal",
              )

              weights = {}
              for syst in ["JES", "LF", "HF", "LFStats1", "LFStats2", "HFStats1", "HFStats2", "cErr1", "cErr2"]:
                for sdir in ["Up", "Down"]:
                  weights[syst+sdir] = bweightcalc.calcEventWeight(
                          jetsForBtagWeight, kind="final", systematic=syst+sdir
                          )
              

#========> ================ Lepton Scale Factors =================
              # For custom made form own JSON files



      # End JSON loop ====================================
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
      
      HVdPhi_reg[0] = 300
      dphi = tree.HCSV_reg_phi - tree.V_phi
      if dphi > math.pi:
          dphi = dphi - 2*math.pi
      elif dphi <= -math.pi:
          dphi = dphi + 2*math.pi
      HVdPhi_reg[0] = dphi
      
      if job.type != 'DATA':

      #========> electron scale factor is here                    
      
        eTrigSFWeight = 1
        eIDLooseSFWeight = 1
      
        for j,lepCorr in lepCorrs.iteritems():
          weight = []
          weight.append(lepCorr.get_2D( tree.vLeptons_pt[0], tree.vLeptons_eta[0])) #return [value, error]
          weight.append(lepCorr.get_2D( tree.vLeptons_pt[1], tree.vLeptons_eta[1]))
          if j.find('WP80') != -1:
            eIDLooseSFWeight = weight[0][0]*weight[1][0]  #second index = 0 -> value
          elif j.find('ScaleFactor_HLT_Ele23_WPLoose_Gsf_v') != -1:
            #eTrigSFWeight = weight[0][0]*weight[1][0]
            SF1 = weight[0][0]
            SF2 = weight[1][0]
            eff1 = tree.vLeptons_Eff_HLT_RunD4p2[0]*0.1801911165 + tree.vLeptons_Eff_HLT_RunD4p3[0]*0.8198088835 #Eff_HLT == 1 for electron in 76X tree
            eff2 = tree.vLeptons_Eff_HLT_RunD4p2[1]*0.1801911165 + tree.vLeptons_Eff_HLT_RunD4p3[1]*0.8198088835
            DR = [999, 999]
            for k in range(0,2):
              for l in range(0,len(tree.trgObjects_hltEle23WPLoose_eta)):
                  dr_ = util_funcs.deltaR(tree.vLeptons_eta[k], tree.vLeptons_phi[k], tree.trgObjects_hltEle23WPLoose_eta[l], tree.trgObjects_hltEle23WPLoose_phi[l])
                  if dr_ < DR[k] and tree.vLeptons_pt[k] > 23:
                      DR[k] = dr_
            
            Ele1pass = DR[0] < 0.5
            Ele2pass = DR[0] < 0.5
 
            if tree.Vtype == 1:
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

          else:
            sys.exit('@ERROR: SF list doesn\'t match json files. Abort')
          
          #=====> ID and trigger for electron need to recalculate since not available in the ntuple==== 
          #FIXME: need to fix trigger scale factor, what ID applied for vLepton Cut or MVA?, need iso scale factor as well 
          if tree.Vtype == 1:
              #lepton_EvtWeight[0] = eIDLooseSFWeight*eTrigSFWeight
              lepton_EvtWeight[0] = eIDLooseSFWeight #only ID and iso. trigger scale store in vLeptons_SFweight_HLT
          
          #=========> Add HLT trigger SF for muon=======
          pTcut = 22;
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
                  if dr_ < DR[k] and tree.vLeptons_pt[k] > 22:
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

          '''DR[0] = 999
          DR[1] = 999

          for k in range(0,2):
            for l in range(0,len(tree.trgObjects_hltEle22eta2p1WPLoose_eta)):
              dr_ = util_funcs.deltaR(tree.vLeptons_eta[k], tree.vLeptons_phi[k], tree.trgObjects_hltEle22eta2p1WPLoose_eta[l], tree.trgObjects_hltEle22eta2p1WPLoose_eta[l])
              if dr_ < DR[k] and tree.vLeptons_pt[k] > 24:
                DR[k] = dr_

          Ele1pass = DR[0] < 0.5
          Ele2pass = DR[0] < 0.5
      
          #print 'vLeptSFw is', vLeptons_SFweight_HLT[0]
          #print 'Vtype is', tree.Vtype
          
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
      
#@@@@@@@@@@@@@fill the index of Fatjet that doesn't overlap with V@@@@@@
      nFatjetCA15ungroomed_notV[0] = 0
      nFatjetCA15pruned_notV[0] = 0
      nFatjetCA15trimmed_notV[0] = 0
      nFatjetCA15softdrop_notV[0] = 0
      nFatjetCA15softdropz2b1_notV[0] = 0
      nFatjetAK08ungroomed_notV[0] = 0
      '''for i in range(0,tree.nFatjetCA15ungroomed):
        dRtmp = deltaR(tree.FatjetCA15ungroomed_eta[i],tree.FatjetCA15ungroomed_phi[i],tree.V_eta,tree.V_phi)
        if dRtmp > 1.5:
          idx_FatjetCA15ungroomed_notV[nFatjetCA15ungroomed_notV[0]] = i
          nFatjetCA15ungroomed_notV[0] += 1
      for i in range(0,tree.nFatjetCA15pruned):
        dRtmp = deltaR(tree.FatjetCA15pruned_eta[i],tree.FatjetCA15pruned_phi[i],tree.V_eta,tree.V_phi)
        if dRtmp > 1.5:
          idx_FatjetCA15pruned_notV[nFatjetCA15pruned_notV[0]] = i
          nFatjetCA15pruned_notV[0] += 1
      for i in range(0,tree.nFatjetCA15trimmed):
        dRtmp = deltaR(tree.FatjetCA15trimmed_eta[i],tree.FatjetCA15trimmed_phi[i],tree.V_eta,tree.V_phi)
        if dRtmp > 1.5:
          idx_FatjetCA15trimmed_notV[nFatjetCA15trimmed_notV[0]] = i
          nFatjetCA15trimmed_notV[0] += 1
      for i in range(0,tree.nFatjetCA15softdrop):
        dRtmp = deltaR(tree.FatjetCA15softdrop_eta[i],tree.FatjetCA15softdrop_phi[i],tree.V_eta,tree.V_phi)
        if dRtmp > 1.5:
          idx_FatjetCA15softdrop_notV[nFatjetCA15softdrop_notV[0]] = i
          nFatjetCA15softdrop_notV[0] += 1
      for i in range(0,tree.nFatjetCA15softdropz2b1):
        dRtmp = deltaR(tree.FatjetCA15softdropz2b1_eta[i],tree.FatjetCA15softdropz2b1_phi[i],tree.V_eta,tree.V_phi)
        if dRtmp > 1.5:
          idx_FatjetCA15softdropz2b1_notV[nFatjetCA15softdropz2b1_notV[0]] = i
          nFatjetCA15softdropz2b1_notV[0] += 1
      for i in range(0,tree.nFatjetAK08ungroomed):
        dRtmp = deltaR(tree.FatjetAK08ungroomed_eta[i],tree.FatjetAK08ungroomed_phi[i],tree.V_eta,tree.V_phi)
        if dRtmp > 0.8:
          idx_FatjetAK08ungroomed_notV[nFatjetAK08ungroomed_notV[0]] = i
          nFatjetAK08ungroomed_notV[0] += 1
     
#@@@@@@@@@@@@sort idx according to pt@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      util_funcs.sortPtFatjet(idx_FatjetCA15ungroomed_notV,tree.FatjetCA15ungroomed_pt,nFatjetCA15ungroomed_notV[0])
      util_funcs.sortPtFatjet(idx_FatjetCA15pruned_notV,tree.FatjetCA15pruned_pt,nFatjetCA15pruned_notV[0])
      util_funcs.sortPtFatjet(idx_FatjetCA15trimmed_notV,tree.FatjetCA15trimmed_pt,nFatjetCA15trimmed_notV[0])
      util_funcs.sortPtFatjet(idx_FatjetCA15softdrop_notV,tree.FatjetCA15softdrop_pt,nFatjetCA15softdrop_notV[0])
      util_funcs.sortPtFatjet(idx_FatjetCA15softdropz2b1_notV,tree.FatjetCA15softdropz2b1_pt,nFatjetCA15softdropz2b1_notV[0])
      util_funcs.sortPtFatjet(idx_FatjetAK08ungroomed_notV,tree.FatjetAK08ungroomed_pt,nFatjetAK08ungroomed_notV[0])
      '''
#@@@@@@@@@@@@@fill matched subjet to leading fat jet@@@@@@@@@@@@@@@@@@@@

      nSubjet_matchCA15ungroomed[0] = 0 
      nSubjet_matchCA15pruned[0] = 0        
      nSubjet_matchCA15trimmed[0] = 0       
      nSubjet_matchCA15softdrop[0] = 0
      nSubjet_matchCA15softdropz2b1[0] = 0
      nSubjet_matchAK08ungroomed[0] = 0
      
      tmp = findFatJetAndSubjet(tree.FatjetCA15ungroomed_eta, tree.FatjetCA15ungroomed_phi, tree.SubjetCA15pruned_eta, tree.SubjetCA15pruned_phi, tree.SubjetCA15pruned_btag, tree.V_eta, tree.V_phi, 1.5)
      if len(tmp) == 3: 
        nFatjetCA15ungroomed_notV[0] = 1
        idx_FatjetCA15ungroomed_notV[0] = tmp[0]
        nSubjet_matchCA15ungroomed[0] = 2
        idx_subjet_matchCA15ungroomed[0] = tmp[1]
        idx_subjet_matchCA15ungroomed[1] = tmp[2]
      if len(tmp) > 0 and len(tmp) != 3: print 'Warning not consistent jet and subjet finding ', len(tmp)
      
      tmp = findFatJetAndSubjet(tree.FatjetCA15pruned_eta, tree.FatjetCA15pruned_phi, tree.SubjetCA15pruned_eta, tree.SubjetCA15pruned_phi, tree.SubjetCA15pruned_btag, tree.V_eta, tree.V_phi, 1.5)
      if len(tmp) == 3: 
        nFatjetCA15pruned_notV[0] = 1
        idx_FatjetCA15pruned_notV[0] = tmp[0]
        nSubjet_matchCA15pruned[0] = 2
        idx_subjet_matchCA15pruned[0] = tmp[1]
        idx_subjet_matchCA15pruned[1] = tmp[2]
      if len(tmp) > 0 and len(tmp) != 3: print 'Warning not consistent jet and subjet finding ', len(tmp)
      
      '''tmp1 = findFatJetTest(tree.nFatjetCA15pruned, tree.FatjetCA15pruned_eta,tree.FatjetCA15pruned_phi,tree.V_eta,tree.V_phi, 1.5)
      if len(tmp1) >= 1:
        nFatjetCA15pruned_notV[0] = 1
        idx_FatjetCA15pruned_notV[0] = tmp1[0]
      '''


      tmp = findFatJetAndSubjet(tree.FatjetCA15trimmed_eta, tree.FatjetCA15trimmed_phi, tree.SubjetCA15pruned_eta, tree.SubjetCA15pruned_phi, tree.SubjetCA15pruned_btag, tree.V_eta, tree.V_phi, 1.5)
      if len(tmp) == 3: 
        nFatjetCA15trimmed_notV[0] = 1
        idx_FatjetCA15trimmed_notV[0] = tmp[0]
        nSubjet_matchCA15trimmed[0] = 2
        idx_subjet_matchCA15trimmed[0] = tmp[1]
        idx_subjet_matchCA15trimmed[1] = tmp[2]
      if len(tmp) > 0 and len(tmp) != 3: print 'Warning not consistent jet and subjet finding ', len(tmp)

      tmp = findFatJetAndSubjet(tree.FatjetCA15softdrop_eta, tree.FatjetCA15softdrop_phi, tree.SubjetCA15pruned_eta, tree.SubjetCA15pruned_phi, tree.SubjetCA15pruned_btag, tree.V_eta, tree.V_phi, 1.5)
      if len(tmp) == 3: 
        nFatjetCA15softdrop_notV[0] = 1
        idx_FatjetCA15softdrop_notV[0] = tmp[0]
        nSubjet_matchCA15softdrop[0] = 2
        idx_subjet_matchCA15softdrop[0] = tmp[1]
        idx_subjet_matchCA15softdrop[1] = tmp[2]
      if len(tmp) > 0 and len(tmp) != 3: print 'Warning not consistent jet and subjet finding ', len(tmp)
      
      tmp = findFatJetAndSubjet(tree.FatjetCA15softdropz2b1_eta, tree.FatjetCA15softdropz2b1_phi, tree.SubjetCA15pruned_eta, tree.SubjetCA15pruned_phi, tree.SubjetCA15pruned_btag, tree.V_eta, tree.V_phi, 1.5)
      if len(tmp) == 3: 
        nFatjetCA15softdropz2b1_notV[0] = 1
        idx_FatjetCA15softdropz2b1_notV[0] = tmp[0]
        nSubjet_matchCA15softdropz2b1[0] = 2
        idx_subjet_matchCA15softdropz2b1[0] = tmp[1]
        idx_subjet_matchCA15softdropz2b1[1] = tmp[2]
      if len(tmp) > 0 and len(tmp) != 3: print 'Warning not consistent jet and subjet finding ', len(tmp)
      
      tmp = findFatJetAndSubjet(tree.FatjetAK08ungroomed_eta, tree.FatjetAK08ungroomed_phi, tree.SubjetAK08softdrop_eta, tree.SubjetAK08softdrop_phi, tree.SubjetAK08softdrop_btag, tree.V_eta, tree.V_phi, 0.8)
      if len(tmp) == 3: 
        nFatjetAK08ungroomed_notV[0] = 1
        idx_FatjetAK08ungroomed_notV[0] = tmp[0]
        nSubjet_matchAK08ungroomed[0] = 2
        idx_subjet_matchAK08ungroomed[0] = tmp[1]
        idx_subjet_matchAK08ungroomed[1] = tmp[2]
      if len(tmp) > 0 and len(tmp) != 3: print 'Warning not consistent jet and subjet finding ', len(tmp)


      '''for i in range(0,tree.nSubjetCA15pruned):
                   
        if nFatjetCA15ungroomed_notV[0] > 0:
          idxTmp = idx_FatjetCA15ungroomed_notV[0]
          #print '>>>>>>>>>>>>>>', idxTmp, len(tree.FatjetCA15ungroomed_eta)
          dRtmp = deltaR(tree.SubjetCA15pruned_eta[i],tree.SubjetCA15pruned_phi[i], tree.FatjetCA15ungroomed_eta[idxTmp], tree.FatjetCA15ungroomed_phi[idxTmp])
          if dRtmp < 1.5:
            idx_subjet_matchCA15ungroomed[nSubjet_matchCA15ungroomed[0]]=i
            nSubjet_matchCA15ungroomed[0] += 1
        
        if nFatjetCA15pruned_notV[0] > 0:
          idxTmp = idx_FatjetCA15pruned_notV[0]
          dRtmp = deltaR(tree.SubjetCA15pruned_eta[i],tree.SubjetCA15pruned_phi[i], tree.FatjetCA15pruned_eta[idxTmp], tree.FatjetCA15pruned_phi[idxTmp])
          if dRtmp < 1.5:
            idx_subjet_matchCA15pruned[nSubjet_matchCA15pruned[0]]=i
            nSubjet_matchCA15pruned[0] += 1

        if nFatjetCA15trimmed_notV[0] > 0:
          idxTmp = idx_FatjetCA15trimmed_notV[0]
          dRtmp = deltaR(tree.SubjetCA15pruned_eta[i],tree.SubjetCA15pruned_phi[i], tree.FatjetCA15trimmed_eta[idxTmp], tree.FatjetCA15trimmed_phi[idxTmp])
          if dRtmp < 1.5:
            idx_subjet_matchCA15trimmed[nSubjet_matchCA15trimmed[0]]=i
            nSubjet_matchCA15trimmed[0] += 1

        if nFatjetCA15softdrop_notV[0] > 0:
          idxTmp = idx_FatjetCA15softdrop_notV[0]
          dRtmp = deltaR(tree.SubjetCA15pruned_eta[i],tree.SubjetCA15pruned_phi[i], tree.FatjetCA15softdrop_eta[idxTmp], tree.FatjetCA15softdrop_phi[idxTmp])
          if dRtmp < 1.5:
            idx_subjet_matchCA15softdrop[nSubjet_matchCA15softdrop[0]]=i
            nSubjet_matchCA15softdrop[0] += 1
         
        if nFatjetCA15softdropz2b1_notV[0] > 0:
          idxTmp = idx_FatjetCA15softdropz2b1_notV[0]
          dRtmp = deltaR(tree.SubjetCA15pruned_eta[i],tree.SubjetCA15pruned_phi[i], tree.FatjetCA15softdropz2b1_eta[idxTmp], tree.FatjetCA15softdropz2b1_phi[idxTmp])
          if dRtmp < 1.5:
            idx_subjet_matchCA15softdrop_z2b1[nSubjet_matchCA15softdrop_z2b1[0]]=i
            nSubjet_matchCA15softdrop_z2b1[0] += 1

      for i in range(0,tree.nSubjetAK08softdrop):
        if  nFatjetAK08ungroomed_notV[0] > 0:
          idxTmp = idx_FatjetAK08ungroomed_notV[0]
          dRtmp = deltaR(tree.SubjetAK08softdrop_eta[i],tree.SubjetAK08softdrop_phi[i], tree.FatjetAK08ungroomed_eta[idxTmp], tree.FatjetAK08ungroomed_phi[idxTmp])
          if dRtmp < 0.8:
            idx_subjet_matchAK08ungroomed[nSubjet_matchAK08ungroomed[0]]=i
            nSubjet_matchAK08ungroomed[0] += 1
      '''

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
    inFile_folder = pathIN + '/' + job.prefix+job.identifier
    tmpD =  config.get('Directories','samplefiles').replace('/','') + '_afterPrepStep/'
    os.system('mkdir ' + tmpD)
    tmpFileList =  tmpD + job.prefix + job.identifier + '.txt'
    os.system('rm -f ' + tmpFileList)
    util_funcs.findSubFolders(inFile_folder,tmpFileList,False)
    
    outputFolder = "%s/%s/" %(pathOUT,job.identifier)
    print "Remove output folder: ", outputFolder
    os.system('rm -rf ' + outputFolder)
     
    try:
        os.mkdir(outputFolder)
    except:
        pass

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
        p = Pool(multiprocess)
        outputs = p.map(fillTreeSingleInput,input_processings)

    else:
        for input_ in input_processings:
                output = fillTreeSingleInput(input_)
                outputs.append(output)



#TEMP for LPC
    
    whereToLaunch = config.get('Configuration','whereToLaunch')
    if ('LPC' in whereToLaunch or 'pisa' in whereToLaunch):
      run_locally = config.get('Configuration', 'run_locally')
      if run_locally != 'False':
        outfileName = pathOUT + '/' + job.identifier + '.root'
        print 'Remove old file before hadd: ', outfileName 
        os.system('rm -f ' + outfileName)
        command = "hadd -f -k " + outfileName + " " + outputFolder + '/*'
        print command
        os.system(command)


    




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
