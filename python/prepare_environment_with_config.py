#! /usr/bin/env python
import os, pickle, sys, ROOT
ROOT.gROOT.SetBatch(True)
from optparse import OptionParser
from myutils import BetterConfigParser, copytree, copytreePSI, ParseInfo
import utils
import pdb

print 'start prepare_environment_with_config.py'

import os
if os.path.exists("../interface/DrawFunctions_C.so"):
    print 'ROOT.gROOT.LoadMacro("../interface/DrawFunctions_C.so")'
    ROOT.gROOT.LoadMacro("../interface/DrawFunctions_C.so")

argv = sys.argv

#get files info from config
parser = OptionParser()
parser.add_option("-C", "--config", dest="config", default=[], action="append",
                      help="directory config")
parser.add_option("-S", "--samples", dest="names", default="",
                              help="samples you want to run on")
parser.add_option("-f", "--filelist", dest="filelist", default="",
                              help="list of files you want to run on")

(opts, args) = parser.parse_args(argv)

config = BetterConfigParser()
config.read(opts.config)

namelist=opts.names.split(',')
filelist=opts.filelist.split(';')
print "namelist:",namelist
print "len(filelist)",len(filelist),"filelist[0]:",filelist[0]

run_locally = config.get('Configuration', 'run_locally')
pathIN = config.get('Directories','PREPin')
pathOUT = config.get('Directories','PREPout')
if run_locally=='False':
  pathOUT = './'

samplesinfo=config.get('Directories','samplesinfo')
sampleconf = BetterConfigParser()
sampleconf.read(samplesinfo)

whereToLaunch = config.get('Configuration','whereToLaunch')
TreeCopierPSI = config.get('Configuration','TreeCopierPSI')
prefix=config.get('General','prefix')

#It list stuff in pathIN and put the name of root file '.root' into 'name', LOOK into myutils/sample_parser.py ~line 108 to see how ParseInfo does 
info = ParseInfo(samplesinfo,pathIN)
#print "samplesinfo:",samplesinfo

#for job in info:
#  print job.name, ' ', job.identifier, '  ', job.addtreecut

for job in info:
    #print ">>>>>>>>job.name: ", job.name
    #print ">>>>>>>>job.identifier: ", job.identifier
    #job.name = samplename or subnames (if subsamples = true ) in samples_nosplit.ini
    #job.identifier = identifier for that sample
    #[DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8] -->identifier
    #  sampleName = HT0to100DYJetsToLL_M-50_TuneCUETP8M1
    #  subnames =['HT0to100ZJets_udscg','HT0to100ZJets_1b','HT0to100ZJets_2b'] -->name

    if namelist[0] != 'All' and not job.name in namelist and not job.identifier in namelist:
        continue

    #skip if it is subsample
    if job.subsample:
        continue
    
    #print '>>>>>>>>>>>>>>>>', job.name, job.identifier, job.subsample
    
    if('lxplus' in whereToLaunch):
        utils.TreeCopier(pathIN, pathOUT, job.identifier, job.prefix, job.addtreecut)
    else:
        if TreeCopierPSI == 'True':
            samplefiles = config.get('Directories','samplefiles')
            copytreePSI(samplefiles,pathOUT,prefix,job.prefix,job.identifier,'',job.addtreecut, config, filelist)
        else:
            #print ">>>>>>>>>>>>>>: ", job.prefix
            #print ">>>>>>>>>>>>>>: ", job.addtreecut
            #job.prefix = <General|newprefix>
            #job.addtreecut = "cut" in sample_noplit.ini
            copytree(pathIN,pathOUT,prefix,job.prefix,job.identifier,'',job.addtreecut, config, filelist)
            pass
    
print 'end prepare_environment_with_config.py'
