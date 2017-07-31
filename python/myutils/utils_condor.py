def make_input_file_list(nFile, outDir_file_list, file_list_name, prefix_file_list_for_condor):
  
  lines = open(file_list_name).readlines()
  tmp = len(lines)
  nJob = tmp/nFile + 1 
  if tmp%nFile == 0: nJob = tmp/nFile
  print '>>>>>>>>Total file, number of file per job, nJobs: ', tmp, nFile, nJob

  iFile = 0 
  for line in range(0, len(lines), nFile):
    tmp = file_list_name.split('/')
    #newFile = open(outDir_file_list + '/' + tmp[len(tmp)-1] + '_' + str(iFile) + '.txt', 'w')  
    newFile = open(outDir_file_list + '/' + prefix_file_list_for_condor + '_' + str(iFile) + '.txt', 'w')  #sampleList_1.txt for example
        
    tmpFiles = lines[line:line+nFile]
    for i in range(0, len(tmpFiles)):
      newFile.write(tmpFiles[i])
    
    iFile = iFile + 1 
  
  return nJob


def write_condor_config(workDir, name_script_to_run, argv_of_script_to_run, name_output_subdir, nJob, cmssw='7_6_5'):
  f = open(workDir + '/condor_config.script', 'w')
  f.write('universe = vanilla \n')
  f.write('Executable = condor_executable.sh \n')
  f.write('Arguments = $(Process) ' + name_output_subdir + ' \n') 
  f.write('Should_Transfer_Files = YES \n')
  f.write('WhenToTransferOutput = ON_EXIT \n')
  f.write('Transfer_Input_Files = input.tar \n')
  f.write('Output = ctagana_$(Cluster)_$(Process).stdout \n')
  f.write('Error = ctagana_$(Cluster)_$(Process).stderr \n')
  f.write('Log = ctagana_$(Cluster)_$(Process).log \n')
  f.write('notify_user = ${LOGNAME}@FNAL.GOV \n')
  f.write('+LENGTH="SHORT" \n')
  #f.write('x509userproxy = /tmp/x509up_u12772 \n')
  f.write('x509userproxy = $ENV(X509_USER_PROXY) \n')
  tmp = 'Queue ' + str(nJob) + '\n'
  f.write(tmp)
  f.close()

  f = open(workDir + '/condor_executable.sh', 'w')
  f.write('#!/bin/bash \n')
  f.write('echo $SHELL \n')
  f.write('source /cvmfs/cms.cern.ch/cmsset_default.sh \n')
  f.write('cd /uscms_data/d3/duong/CMSSW/CMSSW_' + cmssw + '/src \n')
  f.write('eval `scramv1 runtime -sh` \n')
  f.write('cd ${_CONDOR_SCRATCH_DIR} \n')
  f.write('tar -xvf input.tar \n')
  f.write('python ' + name_script_to_run + ' ' + argv_of_script_to_run + '\n')
  #f.write('xrdcp *.root $2 \n')
  #f.write('rm *.root \n')
  f.close()
  
def write_ana_macro(workDir, nJob, treename):
  for i in range(0, nJob):
    f = open(workDir + '/run_' + str(i) + '.py', 'w')
    f.write('import os,sys \n')
    f.write('from myutils import addMoreInfo \n')
    f.write('lines = open("sampleList_' + str(i) + '.txt").readlines() \n')
    f.write('for line in lines: \n')
    f.write('  line = line.replace("\\n","") \n')
    f.write('  outFile = line.split("/")[-1] \n')
    f.write('  addMoreInfo.fillTree(line,outFile,"' + treename + '") \n')
    f.close()
