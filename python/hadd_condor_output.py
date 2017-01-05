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

from myutils import utils_condor

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

from myutils import BetterConfigParser, ParseInfo

print opts.config
config = BetterConfigParser()
config.read(opts.config)
samplesinfo=config.get('Directories','samplesinfo')

pathIN = config.get('Directories','SYSin')
pathOUT = config.get('Directories','SYSout')
#tmpDir = os.environ["TMPDIR"]
tmpDir = config.get("Directories",'tmpDir')

print 'INput samples:\t%s'%pathIN
print 'OUTput samples:\t%s'%pathOUT

namelist=opts.names.split(',')

#load info
info = ParseInfo(samplesinfo,pathIN)



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#if namelist[0] == All, only process indentifier. the info contain all the sub sample with the same identifier. we only need to process once for a identifier
#processed_identifiers = []

scratchDir = config.get('Directories','scratch')
submitDir = scratchDir + '/SubmitToCondor/'

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
    outputFolder = "%s/%s/" %(pathOUT,job.identifier)
    print "Remove output folder: ", outputFolder
    os.system('rm -rf ' + outputFolder)
    
    condor_folder = submitDir + '/' + job.identifier + '/'

    try:
        os.mkdir(outputFolder)
    except:
        pass

    os.system('rm tmp.txt')
    os.system('ls ' + condor_folder + '/*.stderr > tmp.txt')
    
    for line in open('tmp.txt').readlines():
        line = line.replace('\n','')
        if 'Break' in open(line).read():
            print '>>>Break found in job: ', line
        if 'Error' in open(line).read():
            print '>>>Error found in job: ', line

    tmpD =  config.get('Directories','samplefiles').replace('/','') + '_afterPrepStep/'
    input_fileList =  tmpD + job.prefix + job.identifier + '.txt'
    nInputFiles = sum(1 for line in open(input_fileList))
    os.system('rm tmp.txt')
    os.system('ls ' + condor_folder + '/*.root > tmp.txt')
    nOutputFiles = sum(1 for line in open('tmp.txt'))
    if nInputFiles != nOutputFiles: print '>>>>Missing output files in processing ' + job.name + '. Input = ', nInputFiles, ', Output = ', nOutputFiles
    
    
    outfileName = pathOUT + '/' + job.identifier + '.root'
    print 'Remove old file before hadd: ', outfileName 
    os.system('rm -f ' + outfileName)
    command = "hadd -f -k " + outfileName + " " + condor_folder + '/*.root'
    print ">>>>>>>>", command
    os.system(command)
    command = 'mv ' + condor_folder + '/*.root ' + pathOUT + '/' + job.identifier
    print ">>>>>>>>", command
    os.system(command)
