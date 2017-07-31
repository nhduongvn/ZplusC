from __future__ import print_function
import os,sys,subprocess,hashlib
import ROOT
from samplesclass import Sample
import time

class TreeCache:
    def __init__(self, cutList, sampleList, path, config,filelist=[],mergeplot=False):
        ROOT.gROOT.SetBatch(True)
        self.path = path
        self.config = config
        print("Init path",path," sampleList",sampleList)
        self._cutList = []
        #! Make the cut lists from inputs
        for cut in cutList:
            self._cutList.append('(%s)'%cut.replace(' ',''))
        try:
            #self.__tmpPath = os.environ["TMPDIR"]
            #print('The TMPDIR is ', os.environ["TMPDIR"])
            self.__tmpPath = config.get('Directories','tmpDir')
            print('The TMPDIR is ', self.__tmpPath)
            
        except KeyError:
            print("\x1b[32;5m %s \x1b[0m" %open('%s/data/vhbb.txt' %config.get('Directories','vhbbpath')).read())
            print("\x1b[31;5;1m\n\t>>> %s: Please set your TMPDIR and try again... <<<\n\x1b[0m" %os.getlogin())
            sys.exit(-1)

        if filelist:
            print('len(filelist)',len(filelist))

        self.__doCache = True
        if config.has_option('Directories','tmpSamples'):
            self.__cachedPath = config.get('Directories','tmpSamples')
        self.__hashDict = {}
        self.minCut = None
        self.__find_min_cut()# store the cut list as one string in minCut, using ROOT syntax (i.e. || to separate between each cut) 
        self.__sampleList = sampleList
        print('\n\t>>> Caching FILES <<<\n')
        self.__cache_samples(filelist,mergeplot)
        
    def putOptions(self):
        return (self.__sampleList,self.__doCache,self.__tmpPath,self._cutList,self.__hashDict,self.minCut,self.path)
    
    def _mkdir_recursive(self, path):
        sub_path = os.path.dirname(path)
        if not os.path.exists(sub_path):
            self._mkdir_recursive(sub_path)
        if not os.path.exists(path):
            os.mkdir(path)

    def __find_min_cut(self):
        effective_cuts = []
        for cut in self._cutList:
            if not cut in effective_cuts:
                effective_cuts.append(cut)
        self._cutList = effective_cuts
        self.minCut = '||'.join(self._cutList)

    def _trim_tree(self, sample, filelist = [], mergeplot = False, forceReDo = False):

        start_time = time.time()
        print("Caching the sample")
        print("==================\n")

        ''' Create temporary file for each sample '''
        theName = sample.name
        print('Reading sample <<<< %s' %sample)
        # print('path',self.path,'get_path',sample.get_path)
        # print('filelist',filelist)
        inputfiles = []
        outputfiles = []
        tmpfiles = []
        # prefix
        tmpCache = ''
        tmpSource = ''
        theHash = ''

        if len(filelist) == 0:
            source = '%s/%s' %(self.path,sample.get_path)
            inputfiles.append(source)
            checksum = self.get_checksum(source)
            theHash = hashlib.sha224('%s_s%s_%s' %(sample,checksum,self.minCut)).hexdigest()
            self.__hashDict[theName] = theHash
            tmpSource = '%s/tmp_%s.root'%(self.__tmpPath,theHash)
            print('opening '+source)
            tmpfiles.append(tmpSource)
            tmpCache = '%s/tmp_%s.root'%(self.__cachedPath,theHash)
            outputfiles.append(tmpCache)
            print('tmpSource is', tmpSource)
            print('tmpCache is', tmpCache)
            if mergeplot:
                 # mergetreePSI(pathIN, pathOUT,           prefix,  newprefix, folderName,        Aprefix, Acut, config):
                from mergetreePSI import mergetreePSI_def
                mergetreePSI_def(self.path, self.__cachedPath, theHash, "tmp_",    sample.identifier, hashlib.sha224(self.minCut).hexdigest(),      "",   "")
                sys.exit(1)

        else:
            for inputFile in filelist:
                print('inputFile',inputFile)
                subfolder = inputFile.split('/')[-4]
                print('subfolder',subfolder)
                filename = inputFile.split('/')[-1]
                print('filename1',filename)
                filename = filename.split('_')[0]+'_'+subfolder+'_'+filename.split('_')[1]
                print('filename2',filename)
                hash = hashlib.sha224(filename).hexdigest()
                inputFile = "%s/%s/%s" %(self.path,sample.identifier,filename.replace('.root','')+'_'+str(hash)+'.root')
                print('inputFile',inputFile)
                if not os.path.isfile(inputFile.replace('root://t3dcachedb03.psi.ch:1094/','')): continue
                hash = hashlib.sha224(self.minCut).hexdigest()
                outputFile = "%s/%s/%s" %(self.__cachedPath,sample.identifier,filename.replace('.root','')+'_'+str(hash)+'.root')
                print('outputFile',outputFile)
                tmpfile = "%s/%s" %(self.__tmpPath,filename.replace('.root','')+'_'+str(hash)+'.root')
                print('tmpfile',tmpfile)
                if inputFile in inputfiles: continue
                del_protocol = outputFile
                del_protocol = del_protocol.replace('gsidcap://t3se01.psi.ch:22128/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=')
                del_protocol = del_protocol.replace('dcap://t3se01.psi.ch:22125/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=')
                del_protocol = del_protocol.replace('root://t3dcachedb03.psi.ch:1094/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=')
                if os.path.isfile(del_protocol.replace('srm://t3se01.psi.ch:8443/srm/managerv2?SFN=','')):
                    f = ROOT.TFile.Open(outputFile,'read')
                    if not f:
                      print('file is null, adding to input')
                      inputfiles.append(inputFile)
                      outputfiles.append(outputFile)
                      tmpfiles.append(tmpfile)
                      continue
                    # f.Print()
                    if f.GetNkeys() == 0 or f.TestBit(ROOT.TFile.kRecovered) or f.IsZombie():
                        print('f.GetNkeys()',f.GetNkeys(),'f.TestBit(ROOT.TFile.kRecovered)',f.TestBit(ROOT.TFile.kRecovered),'f.IsZombie()',f.IsZombie())
                        print('File', del_protocol.replace('srm://t3se01.psi.ch:8443/srm/managerv2?SFN=',''), 'already exists but is buggy, gonna delete and rewrite it.')
                        #command = 'rm %s' %(outputFile)
                        command = 'srmrm %s' %(del_protocol)
                        subprocess.call([command], shell=True)
                        print(command)
                    else: continue
                inputfiles.append(inputFile)
                outputfiles.append(outputFile)
                tmpfiles.append(tmpfile)
        print('inputfiles',inputfiles,'tmpfiles',tmpfiles)

        ######################################################################
        for inputfile,tmpfile,outputFile in zip(inputfiles,tmpfiles,outputfiles):

            #print('the tmp source is ', tmpSource)
            #print ('self.__doCache',self.__doCache,'self.file_exists(tmpSource)',self.file_exists(tmpSource))
            print("==================================================================")
            print ('The cut is ', self.minCut)
            print("==================================================================\n")
            if self.__doCache and self.file_exists(outputFile) and not forceReDo and len(filelist) == 0:
                print('sample',theName,'skipped, filename=',outputFile)
                return (theName,theHash)
            else:
                if self.__doCache and self.file_exists(tmpfile) and not forceReDo:
                    print ('File exists in TMPDIR, proceeding to the copy')
                    print('sample',theName,'skipped, filename=',tmpfile)
                    command = 'srmcp -2 -globus_tcp_port_range 20000,25000 file:///'+tmpfile+' '+ outputFile.replace('root://t3dcachedb03.psi.ch:1094','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/')
                    #if whereToLaunch == 'LPC': TODO make switch to run on LPC
                    command = 'cp ' + tmpfile + '  ' + outputFile
                    subprocess.call([command], shell=True)
                    if len(filelist) == 0: return (theName,theHash)

            print ('trying to create',tmpfile)
            print ('self.__tmpPath',self.__tmpPath)
            if self.__tmpPath.find('root://t3dcachedb03.psi.ch:1094/') != -1:
                mkdir_command = self.__tmpPath.replace('root://t3dcachedb03.psi.ch:1094/','srm://t3se01.psi.ch/')
                print('mkdir_command',mkdir_command)
                # RECURSIVELY CREATE REMOTE FOLDER ON PSI SE, but only up to 3 new levels
                mkdir_command1 = mkdir_command.rsplit('/',1)[0]
                mkdir_command2 = mkdir_command1.rsplit('/',1)[0]
                mkdir_command3 = mkdir_command2.rsplit('/',1)[0]
                my_user = os.popen("whoami").read().strip('\n').strip('\r')+'/'
                if my_user in mkdir_command3:
                  print ('mkdir_command3',mkdir_command3)
                  subprocess.call(['srmmkdir '+mkdir_command3], shell=True)# delete the files already created ?
                if my_user in mkdir_command2:
                  print ('mkdir_command2',mkdir_command2)
                  subprocess.call(['srmmkdir '+mkdir_command2], shell=True)# delete the files already created ?
                if my_user in mkdir_command1:
                  print ('mkdir_command1',mkdir_command1)
                  subprocess.call(['srmmkdir '+mkdir_command1], shell=True)# delete the files already created ?
                if my_user in mkdir_command:
                  print ('mkdir_command',mkdir_command)
                  subprocess.call(['srmmkdir '+mkdir_command], shell=True)# delete the files already created ?
            else:
                mkdir_command = self.__tmpPath
                print('mkdir_command',mkdir_command)
                # RECURSIVELY CREATE REMOTE FOLDER ON PSI SE, but only up to 3 new levels
                mkdir_command1 = mkdir_command.rsplit('/',1)[0]
                mkdir_command2 = mkdir_command1.rsplit('/',1)[0]
                mkdir_command3 = mkdir_command2.rsplit('/',1)[0]
                my_user = os.popen("whoami").read().strip('\n').strip('\r')+'/'
                if my_user in mkdir_command and not os.path.exists(mkdir_command):
                  print ('mkdir_command',mkdir_command)
                  subprocess.call(['mkdir '+mkdir_command], shell=True)# delete the files already created ?

            try:
                print("Writing: ",tmpfile)
                #! read the tree from the input
                if forceReDo:
                    output = ROOT.TFile.Open(tmpfile,'recreate')
                else:
                    output = ROOT.TFile.Open(tmpfile,'create')
                output.cd()
            except:
                ## in case there are problems go to the next dataset [probably another process is working on this dataset]
                if len(filelist) == 0: return (theName,theHash)
                else: print('PROBLEM WITH FILE!!',tmpfile); continue
            print ('reading inputfile',inputfile)
            print ("I am reading")
            input = ROOT.TFile.Open(inputfile,'read')
            input.Print()
            print ("I read the tree")
            tree = input.Get(sample.tree)
            assert type(tree) is ROOT.TTree

            input.cd()
            obj = ROOT.TObject
            for key in ROOT.gDirectory.GetListOfKeys():
                input.cd()
                obj = key.ReadObj()
                if obj.GetName() == 'tree':
                    continue
                output.cd()
                obj.Write(key.GetName())
            output.cd()
            theCut = self.minCut
            if sample.subsample:
                theCut += '& (%s)' %(sample.subcut)
            print ("the cut is", theCut)
            #Problem here: not working when empty tree
            cuttedTree=tree.CopyTree(theCut)
            cuttedTree.Write()
            output.Write()
            input.Close()
            del input
            output.Close()
    #        tmpSourceFile = ROOT.TFile.Open(tmpSource,'read')
    #        if tmpSourceFile.IsZombie():
    #            print("@ERROR: Zombie file")
            del output
            print ("debug4")
            print ("I've done " + theName + " in " + str(time.time() - start_time) + " s.")
            print ("Copying file to the tmp folder.")
            print('outputFile',outputFile)
            command = 'srmcp -2 -globus_tcp_port_range 20000,25000 file:///'+tmpfile+' '+ outputFile.replace('root://t3dcachedb03.psi.ch:1094','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/')
            #if whereToLaunch == 'LPC': TODO make switch to run on LPC
            command = 'cp ' + tmpfile + ' ' + outputFile  
            subprocess.call([command], shell=True)
            if '/scratch/' in  tmpfile: command = 'rm %s' %(tmpfile)
            else: command = 'gfal-rm %s' %(tmpfile)
            if len(filelist) == 0: return (theName,theHash)

            # ######################################################################
            # print('Exit loop')
            # newtree.AutoSave()
            # print('Save')
            # output.Close()
            # print('Close')
            # targetStorage = pathOUT.replace('gsidcap://t3se01.psi.ch:22128/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=')+'/'+job.prefix+job.identifier+'.root'
            # if('pisa' in config.get('Configuration','whereToLaunch')):
                # command = 'cp %s %s' %(tmpDir+'/'+job.prefix+job.identifier+'.root',targetStorage)
                # print(command)
                # subprocess.call([command], shell=True)
            # else:
                # command = 'srmmkdir %s' %(pathOUT.replace('gsidcap://t3se01.psi.ch:22128/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=')+'/'+job.identifier).replace('root://t3dcachedb03.psi.ch:1094/','srm://t3se01.psi.ch/')
                # print(command)
                # subprocess.call([command], shell=True)
                # if len(filelist) == 0:
                    # command = 'srmrm %s' %(targetStorage.replace('root://t3dcachedb03.psi.ch:1094/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/'))
                    # print(command)
                    # os.system(command)
                    # command = 'env -i X509_USER_PROXY=/shome/$USER/.x509up_u`id -u` gfal-copy file:////%s %s' %(tmpDir.replace('/mnt/t3nfs01/data01','')+'/'+job.prefix+job.identifier+'.root',targetStorage.replace('root://t3dcachedb03.psi.ch:1094/','srm://t3se01.psi.ch/'))
                    # print(command)
                    # os.system(command)
                # else:
                    # srmpathOUT = pathOUT.replace('gsidcap://t3se01.psi.ch:22128/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=').replace('dcap://t3se01.psi.ch:22125/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=').replace('root://t3dcachedb03.psi.ch:1094/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=')
                    # command = 'srmcp -2 -globus_tcp_port_range 20000,25000 file:///'+tmpfile+' '+outputFile.replace('root://t3dcachedb03.psi.ch:1094/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=')
                    # print(command)
                    # subprocess.call([command], shell=True)

                    # print('checking output file',outputFile)
                    # f = ROOT.TFile.Open(outputFile,'read')
                    # if not f or f.GetNkeys() == 0 or f.TestBit(ROOT.TFile.kRecovered) or f.IsZombie():
                        # print('TERREMOTO AND TRAGEDIA: THE MERGED FILE IS CORRUPTED!!! ERROR: exiting'
                        # sys.exit(1)

                    # command = 'rm '+tmpfile
                    # print(command)
                    # subprocess.call([command], shell=True)
            # ######################################################################

    def __cache_samples(self,filelist=[],mergeplot=False):
        inputs=[]
        for job in self.__sampleList:
            # test1 = [method for method in dir(job) if callable(getattr(job, method))]
            # print(test1)
            # print(dir(job))
            # print('__doc__',job.__doc__,'\n__eq__',job.__eq__,'\n__init__',job.__init__,'\n__module__',job.__module__,'\n__str__',job.__str__,'\nactive',job.active,'\naddtreecut',job.addtreecut,'\ncount',job.count,'\nget_path',job.get_path,'\ngroup',job.group,'\nidentifier',job.identifier,'\nlumi',job.lumi,'\nname',job.name,'\nprefix',job.prefix,'\nsf',job.sf,'\nspecialweight',job.specialweight,'\nsubcut',job.subcut,'\nsubsample',job.subsample,'\ntree',job.tree,'\ntreecut',job.treecut,'\ntype',job.type,'\nweightexpression',job.weightexpression,'\nxsec',job.xsec)
            # sys.exit()

            samplematch = False
            if filelist: samplematch = '/'+job.identifier+'/' in filelist[0]
            if samplematch: print('job.name',job.name,'samplematch',samplematch)#,'filelist',filelist)
            if not filelist or samplematch:
                inputs.append((self,"_trim_tree",(job),(filelist),(mergeplot)))
        multiprocess=0
        # if('pisa' in self.config.get('Configuration','whereToLaunch')):
        multiprocess=int(self.config.get('Configuration','nprocesses_treecache'))
        outputs = []
        print('launching __cache_samples with ',multiprocess,' processes')
        if multiprocess>1:
            from multiprocessing import Pool
            from myutils import GlobalFunction
            p = Pool(multiprocess)
            outputs = p.map(GlobalFunction, inputs)
        else:
            for input_ in inputs:
                #CACHING OF THE FILES HAPPENS HERE!!!!! THIS IS EQUIVALENT TO
                # self._trim_tree(job,filenames,mergeplot)
                outputs.append((getattr(input_[0],input_[1])(input_[2],input_[3],input_[4])))
        if len(filelist)==0:
        # for output in outputs:
            (theName,theHash) = outputs[0]
            self.__hashDict[theName]=theHash

    def get_tree(self, sample, cut):
        print('input file %s/tmp_%s.root'%(self.__cachedPath,self.__hashDict[sample.name]))
        # print ('Opening %s/tmp_%s.root'%(self.__tmpPath,self.__hashDict[sample.name]))
        input = ROOT.TFile.Open('%s/tmp_%s.root'%(self.__cachedPath,self.__hashDict[sample.name]),'read')
        print ('Opening %s/tmp_%s.root'%(self.__cachedPath,self.__hashDict[sample.name]))
        try:
            tree = input.Get(sample.tree)
            assert type(tree) is ROOT.TTree
        except:
            print ("%s/tmp_%s.root is corrupted. I'm relaunching _trim_tree"%(self.__cachedPath,self.__hashDict[sample.name]))
            self._trim_tree(sample, [], False,forceReDo=True)
            input = ROOT.TFile.Open('%s/tmp_%s.root'%(self.__cachedPath,self.__hashDict[sample.name]),'read')
            tree = input.Get(sample.tree)
            print("Type of sample.tree ROOT.TTree? (again) ", type(tree) is ROOT.TTree)


        #fill all Count* histos as lists, like self.CountWeighted = [123.23]
        for obj in input.GetListOfKeys():
            name = obj.GetName()
            if "Count" in name:
                obj = obj.ReadObj()
                assert(type(obj) is ROOT.TH1F)
                counts = []
                for i in range(obj.GetNbinsX()):
                    value = obj.GetBinContent(i+1)
                    if value<=0:
                        print("WARNING: bin ",i+1," of ",name," is ",value,". I'm forcing it to be 1.")
                        value=1
                    counts.append(value)

                setattr(self,name,counts)
        
        #VHbbNameSpace=config.get('VHbbNameSpace','library')
        #ROOT.gSystem.Load(VHbbNameSpace)
        if sample.subsample:
            cut += '& (%s)' %(sample.subcut)
        print('cut is', cut)
        ROOT.gROOT.cd()
        print('getting the tree after applying cuts')
        cuttedTree=tree.CopyTree(cut)
        cuttedTree.SetDirectory(0)
        input.Close()
        del input
        del tree
        return cuttedTree

    @staticmethod
    def get_slc_version():
        command = 'uname -a'
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=True)
        lines = p.stdout.readlines()
        line = lines[0].split()[2]
        if 'el5' in line:
            return 'SLC5'
        elif 'el6' in line:
            return 'SLC6'
        else:
            sys.exit(-1)

    @staticmethod
    def get_scale(sample, config, lumi = None, count=1):
#        print float(sample.lumi)
        try: sample.xsec = sample.xsec[0]
        except: pass
        anaTag=config.get('Analysis','tag')
        theScale = 1.
        lumi = float(sample.lumi)
        theScale = lumi*sample.xsec*sample.sf/(count)
        print("sample: ",sample,"lumi: ",lumi,"xsec: ",sample.xsec,"sample.sf: ",sample.sf,"count: ",count," ---> using scale: ", theScale)
        return theScale

    @staticmethod
    def get_checksum(file):
        if 'gsidcap://t3se01.psi.ch:22128' in file:
            srmPath = 'srm://t3se01.psi.ch:8443/srm/managerv2?SFN='
            if TreeCache.get_slc_version() == 'SLC5':
                command = 'lcg-ls -b -D srmv2 -l %s' %file.replace('gsidcap://t3se01.psi.ch:22128/','%s/'%srmPath)
                p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=True)
                lines = p.stdout.readlines()
                if any('No such' in line for line in lines):
                    print('File not found')
                    print(command)
                line = lines[1].replace('\t* Checksum: ','')
                checksum = line.replace(' (adler32)\n','')
            elif TreeCache.get_slc_version() == 'SLC6':
                command = 'srmls -l %s' %file.replace('gsidcap://t3se01.psi.ch:22128/','%s/'%srmPath)
                p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=True)
                lines = p.stdout.readlines()
                if any('does not exist' in line for line in lines):
                    print('File not found')
                    print(command)
                checksum = lines[6].replace('- Checksum value:','')
                checksum = checksum.strip()
                #srmPath = 'srm://t3se01.psi.ch'
                #command = 'gfal-sum %s ADLER32' %file.replace('gsidcap://t3se01.psi.ch:22128/','%s/'%srmPath)
                #print(command)
                #p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=True)
                #lines = p.stdout.readlines()
                #if any('No such' in line for line in lines):
                #    print('File not found')
                #    print(command)
                #checksum = lines[0].split()[1]

        else:
            command = 'md5sum %s' %file
            p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=True)
            lines = p.stdout.readlines()
            checksum = lines[0]
        return checksum
    
    @staticmethod
    def file_exists(file):
        print ('Will now check if the file exists')
        print ('=================================\n')

        file_exists = False

        file_dummy = file
        #srmPath = 'srm://t3se01.psi.ch:8443/srm/managerv2?SFN='
        file_dummy = file_dummy.replace('root://t3dcachedb03.psi.ch:1094/','')
        file_dummy = file_dummy.replace('srm://t3se01.psi.ch:8443/srm/managerv2?SFN=','')

        print ('The command is', 'os.path.isfile(',file_dummy,')', os.path.isfile(file_dummy))
        if os.path.isfile(file_dummy):
            print(file_dummy, 'exists.')
            f = ROOT.TFile.Open(file,'read')
            if not f:
                print ('File is null. Gonna redo it.')
                del_protocol = file.replace('gsidcap://t3se01.psi.ch:22128/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=').replace('dcap://t3se01.psi.ch:22125/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=').replace('root://t3dcachedb03.psi.ch:1094/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=')
                if '/scratch/' in  del_protocol: command = 'rm %s' %(del_protocol)
                else: command = 'gfal-rm %s' %(del_protocol)
                subprocess.call([command], shell=True)
                print(command)
            elif f.GetNkeys() == 0 or f.TestBit(ROOT.TFile.kRecovered) or f.IsZombie():
                print ('f.GetNkeys()',f.GetNkeys(),'f.TestBit(ROOT.TFile.kRecovered)',f.TestBit(ROOT.TFile.kRecovered),'f.IsZombie()',f.IsZombie())
                print ('File', file_dummy, 'already exists but is buggy, gonna delete and rewrite it.')
                del_protocol = file.replace('gsidcap://t3se01.psi.ch:22128/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=').replace('dcap://t3se01.psi.ch:22125/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=').replace('root://t3dcachedb03.psi.ch:1094/','srm://t3se01.psi.ch:8443/srm/managerv2?SFN=')
                # print('del_protocol',del_protocol)
                if '/scratch/' in  del_protocol: command = 'rm %s' %(del_protocol)
                else: command = 'gfal-rm %s' %(del_protocol)
                subprocess.call([command], shell=True)
                print(command)
            else:
                file_exists = True

        #exist = os.path.exists(file_dummy)
        #print('os.path.exists(',file_dummy,')',exist)
        #return os.path.exists(file_dummy)
        return file_exists


