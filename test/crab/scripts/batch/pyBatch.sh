#!/bin/bash -f
#
# PhysQC processing script
#
# Expects output directory and full input file path as input. 
#

# Avoid wild cards
set noglob

##### CONFIGURATION ##############################################
# Where to store local results (will be copied over at the end)
TOPWORKDIR=/scratch/`whoami`
##################################################################
JOBDIR=sgejob-$JOB_ID-$DIR

############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N p_run

# Run time soft and hard limits hh:mm:ss
# soft=CPU time; hard=Wallclock time
# #$ -l s_rt=40:00:00,h_rt=42:00:00

### Specify the queue on which to run
#$ -q all.q

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err

##### MONITORING/DEBUG INFORMATION ###############################
DATE_START=`date +%s`
echo "################################################################"
echo "Job $JOB_ID started at " `date`
echo "################################################################"

##### RUNNING ####################################################

echo "Arguments: $EXE, $ARGUMENTS"
CMSSW_DIR=/shome/predragm/CMSSW_3_8_6/
RUNDIR=`pwd`
echo "Pwd: " `pwd`

# Set the CMSSW environment (for ROOT, mainly...)
source $VO_CMS_SW_DIR/cmsset_default.sh
#source /swshare/psit3/etc/profile.d/cms_ui_env.sh
#cmsenv

cd $CMSSW_DIR/src
eval `scramv1 runtime -sh`
if test $? -ne 0; then
   echo "ERROR: Failed to source scram environment" >&2
   exit 1
fi

# Fix problem with DCache and ROOT
LD_LIBRARY_PATH="/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH"

# Run the job
cd $RUNDIR
echo "./$EXE $ARGUMENTS"
python `pwd`/$EXE $ARGUMENTS

###########################################################################
DATE_END=`date +%s`
RUNTIME=$((DATE_END-DATE_START))
echo "################################################################"
echo "Job $JOB_ID finished at " `date`
echo "Wallclock running time: $RUNTIME s"
echo "################################################################"
echo " "
exit 0
