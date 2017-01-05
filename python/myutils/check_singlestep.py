#! /usr/bin/env python
import ROOT,sys,os,subprocess,random,string,hashlib
ROOT.gROOT.SetBatch(True)
from printcolor import printc
import pickle
import sys
from optparse import OptionParser
from BetterConfigParser import BetterConfigParser
from sample_parser import ParseInfo
import StackMaker, HistoMaker

if __name__ == "__main__":

    print 'start check_singlestep.py'

    argv = sys.argv

    #get files info from config
    parser = OptionParser()
    parser.add_option("-C", "--config", dest="config", default=[], action="append",
                          help="directory config")
    parser.add_option("-S", "--samples", dest="names", default="",
                                  help="samples you want to run on")
    parser.add_option("-R", "--region", dest="region", default="",
                                  help="region to plot")

    parser.add_option("-t", "--task", dest="task", default="",
                                  help="task (prep, sys, eval, plot) to be checked")

    (opts, args) = parser.parse_args(argv)

    vhbbPlotDef=opts.config[0].split('/')[0]+'/vhbbPlotDef.ini'
    opts.config.append(vhbbPlotDef)#adds it to the config list

    config = BetterConfigParser()
    config.read(opts.config)

    samplesinfo=config.get('Directories','samplesinfo')
    sampleconf = BetterConfigParser()
    sampleconf.read(samplesinfo)

    if opts.task == 'checksingleprep':
        pathOUT_orig = config.get('Directories','PREPout')
    elif opts.task == 'checksinglesys':
        pathOUT_orig = config.get('Directories','SYSout')
    elif opts.task == 'checksingleeval':
        pathOUT_orig = config.get('Directories','MVAout')
    elif opts.task == 'checksingleplot':
        pathOUT_orig = config.get('Directories','tmpSamples')

    info = ParseInfo(samplesinfo,pathOUT_orig)
    pathOUT_orig = pathOUT_orig.replace('gsidcap://t3se01.psi.ch:22128/','').replace('dcap://t3se01.psi.ch:22125/','').replace('root://t3dcachedb03.psi.ch:1094/','')
    print "opts.task",opts.task

    hash = ''
    if opts.region:
        print 'opts.region',opts.region
        region = opts.region
        print 'evaluating cuts for region',region
        section='Plot:%s'%region
        vars = (config.get(section, 'vars')).split(',')#get the variables to be ploted in each region
        print 'vars',vars
        SignalRegion = False
        if config.has_option(section,'Signal'):
            # mc.append(config.get(section,'Signal'))
            SignalRegion = True
        # # #GETALL AT ONCE
        options = []
        Stacks = []
        #print "Start Loop over the list of variables(to fill the StackMaker )" print "==============================================================\n"
        for i in range(len(vars)):# loop over the list of variables to be ploted in each reagion
            # print "The variable is ", vars[i], "\n"
            Stacks.append(StackMaker.StackMaker(config,vars[i],region,SignalRegion))# defined in myutils DoubleStackMaker. The StackMaker merge together all the informations necessary to perform the plot (plot region, variables, samples and signal region ). "options" contains the variables information, including the cuts.
            options.append(Stacks[i].options)
        cuts = []
        for option in options:
            cuts.append(option['cut'])
        _cutList = []
        #! Make the cut lists from inputs
        for cut in cuts:
            _cutList.append('(%s)'%cut.replace(' ',''))
        effective_cuts = []
        for cut in _cutList:
            if not cut in effective_cuts:
                effective_cuts.append(cut)
        _cutList = effective_cuts
        minCut = '||'.join(_cutList)

        hash = hashlib.sha224(minCut).hexdigest()
        print 'pathOUT_orig', pathOUT_orig
        print hash

    mc_dataset_missing_files = []
    data_dataset_missing_files = []
    dataset_identifiers = []
    # print "info:",info
    for job in info:
        dataset_identifiers.append(job.identifier)
    dataset_identifiers = set(dataset_identifiers)
    for identifier in dataset_identifiers:
        sampleType = config.get(identifier,'sampleType')
        print 'sampleType',sampleType

        # print identifier
        # identifier=opts.names.split(',')[0]
        print "identifier:",identifier
        pathOUT = pathOUT_orig+'/'+identifier
        samplefiles = config.get('Directories','samplefiles')
        filenames = open(samplefiles+'/'+identifier+'.txt').readlines()
        print 'number of files on DAS:',len(filenames),#filenames[0]

        # print 'pathOUT', pathOUT
        # print('ls '+pathOUT+'/*.root |wc -l')
        grep_hash = ''
        if hash != '':
            grep_hash = ' |grep '+hash
        nFilesInPathOut = int(os.popen('ls '+pathOUT+'/*.root '+grep_hash+'|wc -l').read())
        firstfileInPathOut = os.popen('ls '+pathOUT+'/*.root '+grep_hash+'|head -n 1').read()
        print '    --->>    number of files produced by '+opts.task+':',nFilesInPathOut#,firstfileInPathOut[0]

        if len(filenames) == int(nFilesInPathOut):
            print 'ALL FILES CORRECTLY PROCESSED BY '+opts.task+'\n'
            # sys.exit(0)
        else:
            print '-------->>>>>> MISSING '+str(len(filenames)-int(nFilesInPathOut))+'/'+str(len(filenames))+' FILES IN THE '+opts.task+' TASK FOR THE SAMPLE '+identifier+'\n'
            # sys.exit(1)
            if sampleType == 'DATA':
                data_dataset_missing_files.append(identifier)
            else:
                mc_dataset_missing_files.append(identifier)


    print '\n\nFINAL RECAP: \n\nmissing files for the following MC datasets:'
    print mc_dataset_missing_files
    print '\nmissing files for the following DATA datasets (VERY IMPORTANT!!!):'
    print data_dataset_missing_files

