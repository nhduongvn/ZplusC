#!/bin/bash

#====================================================================
#
#        FILE: runAll.sh
#
#       USAGE: runAll.sh sample energy task
#
# DESCRIPTION: Script to be launched in the batch system.
#              Can also be used, with some care, to run locally.
#
#      AUTHOR: VHbb team
#              ETH Zurich
#
#====================================================================

# Fix Python escape sequence bug.
export TERM=""

#-------------------------------------------------
# Parse Input Arguments

sample=$1     # The sample to run on. It must match a sampleName in samples_nosplit.ini.
tag=$2        # The analysis configuration tag, e.g. 13TeV.
task=$3       # The task to perform.
nprocesses=$4 # Dummy variable used to shift the other parameters by +1. FIXME: Remove this argument?
job_id=$5     # Needed for split step and train optimisation. FIXME: It does not have a unique meaning.
bdt_params=$6 # The set of hyperparameters for BDT optimisation.
filelist=$7   # Needed to run the prep and sys step with a limited number of files per job, could be text file with a list of file in the format of xxx.txt.
region=$8     # plotting region

# Verify the number of input arguments.
if [ $# -lt 3 ]; then
    echo "RuntimeError: At least 3 arguments are required, e.g."
    echo "./runAll.sh sample tag task"
    exit
elif [ $task = "mva_opt" -a $# -lt 5 ]; then
    echo "RuntimeError: At least 5 arguments are required for BDT hyperparameter optimisation, e.g."
    echo "./runAll.sh sample tag mva_opt job_id bdt_params"
    exit
fi

echo "Task: $task"
echo

#-------------------------------------------------
# Setup Environment

# Change to the job submission directory when using lxbatch.
if [ -n "${LS_SUBCWD-}" ]; then
    cd $LS_SUBCWD
fi

echo "Parsing files in ${tag}config..."
echo

# The job submission environment.
whereToLaunch=$(
python - << END
import myutils
parser = myutils.BetterConfigParser()
parser.read('${tag}config/general.ini')
print parser.get('Configuration', 'whereToLaunch')
END
)

echo "whereToLaunch: $whereToLaunch"
echo

if [ $whereToLaunch == "pisa" ]; then
    export TMPDIR=$CMSSW_BASE/src/tmp
fi

# The configuration file names, formatted as arguments to the task scripts.
config_filenames=( $(
python - << END
import myutils
parser = myutils.BetterConfigParser()
parser.read('${tag}config/general.ini')
print parser.get('Configuration', 'List')
END
) )

echo "Configuration Files: ${config_filenames[@]}"
echo

for (( i=0; i<${#config_filenames[@]}; i++ )); do
    config_filenames[$i]="--config ${tag}config/${config_filenames[$i]}"
done

# The log file path.
logpath=$(
python - << END
import myutils
parser = myutils.BetterConfigParser()
parser.read('${tag}config/paths.ini')
print parser.get('Directories', 'logpath')
END
)

if [ ! -d $logpath ]; then
    mkdir -p $logpath
fi

# The MVA list.
MVAList=$(
python << END
import myutils
parser = myutils.BetterConfigParser()
parser.read('${tag}config/training.ini')
print parser.get('MVALists', 'List_for_submitscript')
END
)

#-------------------------------------------------
# Run Task

if [ $task = "prep" ]; then
    echo "./prepare_environment_with_config.py --samples $sample ${config_filenames[@]}"
    ./prepare_environment_with_config.py --samples $sample ${config_filenames[@]}

elif [ $task = "singleprep" ]; then
    echo "./prepare_environment_with_config.py --samples $sample ${config_filenames[@]} --filelist $filelist"
    ./prepare_environment_with_config.py --samples $sample ${config_filenames[@]} --filelist $filelist

elif [ $task = "mergesingleprep" ]; then
    echo "./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]}"
    ./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]}

elif [ $task = "trainReg" ]; then
    echo "./trainRegression.py --config ${tag}config/regression.ini ${config_filenames[@]}"
    ./trainRegression.py --config ${tag}config/regression.ini ${config_filenames[@]}

elif [ $task = "sys" ]; then
    echo "./write_regression_systematics_V24.py --samples $sample ${config_filenames[@]}"
    ./write_regression_systematics_V24.py --samples $sample ${config_filenames[@]}

elif [ $task = "sys_condor" ]; then
    echo "./write_regression_systematics_condor.py --samples $sample ${config_filenames[@]}"
    ./write_regression_systematics_condor.py --samples $sample ${config_filenames[@]}

elif [ $task = "sys_condor_hadd" ]; then
    echo "./hadd_condor_output.py --samples $sample ${config_filenames[@]}"
    ./hadd_condor_output.py --samples $sample ${config_filenames[@]}

elif [ $task = "singlesys" ]; then
    echo "./write_regression_systematics.py --samples $sample ${config_filenames[@]} --filelist $filelist"
    ./write_regression_systematics.py --samples $sample ${config_filenames[@]} --filelist $filelist

elif [ $task = "mergesinglesys" ]; then
    echo "./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]}  --mergesys True"
    ./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]} --mergesys True

elif [ $task = "reg" ]; then
    echo "./only_regression.py --samples $sample ${config_filenames[@]}"
    ./only_regression.py --samples $sample ${config_filenames[@]}

elif [ $task = "eval" ]; then
    echo "./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]}"
    ./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]}

elif [ $task = "singleeval" ]; then
    echo "./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]} --filelist $filelist"
    ./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]} --filelist $filelist

elif [ $task = "mergesingleeval" ]; then
    echo "./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]}  --mergeeval True"
    ./myutils/mergetreePSI.py --samples $sample ${config_filenames[@]} --mergeeval True

elif [ $task = "syseval" ]; then
    echo "./write_regression_systematics.py --samples $sample ${config_filenames[@]}"
    ./write_regression_systematics.py --samples $sample ${config_filenames[@]}
    echo "./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]}"
    ./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]}

elif [ $task = "train" ]; then
    echo "./train.py --training $sample ${config_filenames[@]} --local True"
    ./train.py --training $sample ${config_filenames[@]} --local True

elif [ $task = "singleplot" ]; then
    echo "./tree_stack.py --region $region ${config_filenames[@]} --filelist $filelist"
    ./tree_stack.py --region $region ${config_filenames[@]} --filelist $filelist

elif [ $task = "mergesingleplot" ]; then
    echo "./tree_stack.py --region $region ${config_filenames[@]} --mergeplot True"
    ./tree_stack.py --region $region ${config_filenames[@]} --mergeplot True

elif [ $task = "plot" ]; then
    echo "./tree_stack.py --region $sample ${config_filenames[@]}"
    ./tree_stack.py --region $sample ${config_filenames[@]}

elif [ $task = "checksingleprep" ] || [ $task = "checksinglesys" ] || [ $task = "checksingleeval" ] || [ $task = "checksingleplot" ]; then
    if [[ $region ]]; then
        echo "./myutils/check_singlestep.py --region $region ${config_filenames[@]} --task $task"
        ./myutils/check_singlestep.py --region $region ${config_filenames[@]} --task $task
    else
        echo "./myutils/check_singlestep.py ${config_filenames[@]} --task $task"
        ./myutils/check_singlestep.py ${config_filenames[@]} --task $task
    fi
elif [ $task = "dc" ]; then
    echo "./workspace_datacard_working.py --variable $sample ${config_filenames[@]}"
    ./workspace_datacard_working.py --variable $sample ${config_filenames[@]}

elif [ $task = "split" ]; then
    echo "./split_tree.py --samples $sample ${config_filenames[@]} --max-events $job_id"
    ./split_tree.py --samples $sample ${config_filenames[@]} --max-events $job_id

elif [ $task = "stack" ]; then
    echo "./manualStack.py --config ${config_filenames[@]}"
    ./manualStack.py ${config_filenames[@]}

elif [ $task = "plot_sys" ]; then
    echo "./plot_systematics.py ${config_filenames[@]}"
    ./plot_systematics.py ${config_filenames[@]}

elif [ $task = "mva_opt" ]; then
    echo "BDT Hyperparameters: $bdt_params"
    echo "./train.py --name ${sample} --training ${job_id} ${config_filenames[@]} --setting $bdt_params  --local True"
    ./train.py --name ${sample} --training ${job_id} ${config_filenames[@]} --setting $bdt_params --local True

elif [ $task = "mva_opt_eval" ]; then
    echo "./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]} --weight $bdt_params"
    ./evaluateMVA.py --discr $MVAList --samples $sample ${config_filenames[@]} --weight $bdt_params

# WORK IN PROGRESS
elif [ $task = "mva_opt_dc" ]; then
    echo "./workspace_datacard.py --variable $sample ${config_filenames[@]} --optimisation $bdt_params"
    ./workspace_datacard.py --variable $sample ${config_filenames[@]} --optimisation $bdt_params

fi

echo
echo "Exiting runAll.sh"
echo
