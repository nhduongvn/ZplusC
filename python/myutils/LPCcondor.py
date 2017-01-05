def write_condor_config(workDir, job_id, option, dir_root_file_on_condor, dir_root_file_out):
  f = open(workDir + '/condor_config_' + job_id + '.script', 'w')
  f.write('universe = vanilla \n')
  f.write('Executable = condor_executable_' + job_id + '.sh \n')
  f.write('Arguments = ' + option + ' \n') 
  f.write('Should_Transfer_Files = YES \n')
  f.write('WhenToTransferOutput = ON_EXIT \n')
  f.write('Transfer_Input_Files = input.tar \n')
  f.write('Output = log_$(Cluster)_$(Process).stdout \n')
  f.write('Error = log_$(Cluster)_$(Process).stderr \n')
  f.write('Log = log_$(Cluster)_$(Process).log \n')
  f.write('notify_user = ${LOGNAME}@FNAL.GOV \n')
  f.write('+LENGTH="SHORT" \n')
  f.write('x509userproxy = /tmp/x509up_u12772 \n')
  tmp = 'Queue 1 \n'
  f.write(tmp)
  f.close()

  f = open(workDir + '/condor_executable_' + job_id + '.sh', 'w')
  f.write('#!/bin/bash \n')
  f.write('echo $SHELL \n')
  f.write('source /cvmfs/cms.cern.ch/cmsset_default.sh \n')
  f.write('cd /uscms_data/d3/duong/CMSSW/CMSSW_7_6_5/src/ \n')
  f.write('eval `scramv1 runtime -sh` \n')
  f.write('cd ${_CONDOR_SCRATCH_DIR} \n')
  f.write('tar -xvf input.tar \n')
  f.write('./runAll.sh $1 $2 $3 $4 $5 $6 $7 \n')
#  f.write('ls -lhtr ' + dir_root_file_on_condor + '/ \n')
  f.write('mv ' + dir_root_file_on_condor + '/* . \n')
  f.write('ls -lhtr \n')
#  f.write('xrdcp ' + dir_root_file_on_condor + '/*.root ' + 'root://cmseos.fnal.gov/' + dir_root_file_out + '/ \n')
  f.write('xrdcp *.root ' + 'root://cmseos.fnal.gov/' + dir_root_file_out + '/ \n')
  #f.write('rm *.root \n')
  f.close()

#def hadd_skim_tree_afterCondor(sample, dir_in, dir_out):
#hadd dir_out `xrdfsls -u /store/user/foo/bar/ | grep ".root"`

