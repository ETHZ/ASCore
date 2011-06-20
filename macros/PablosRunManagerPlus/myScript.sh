#!/bin/bash
#################################
# PSI Tier-3 example batch Job  #
#################################

# 64bit & all queue version!

##### CONFIGURATION ##############################################
# Output files to be copied back to the User Interface
# (the file path must be given relative to the working directory)
OUTFILES="myout_job_$1.txt myerr_job_$1.txt"

# Output files to be copied to the SE
# (as above the file path must be given relative to the working directory)
SEOUTFILES="output_$1.root"

#SOURCEFILES="dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/susy/ntuples/mc/V01-03-01/ZJets-madgraph-V01-03-01.root"

SOURCEFILES=sourcefile

#
# By default, the files will be copied to $USER_SRM_HOME/$username/$JOBDIR,
# but here you can define the subdirectory under your SE storage path
# to which files will be copied (uncomment line)
#SEUSERSUBDIR="mytestsubdir/somedir"
#
# User's CMS hypernews name (needed for user's SE storage home path
# USER_SRM_HOME below)

HN_NAME=hpname

# set DBG=1 for additional debug output in the job report files
# DBG=2 will also give detailed output on SRM operations
DBG=0

#### The following configurations you should not need to change
# The SE's user home area (SRMv2 URL)
USER_SRM_HOME="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/"

# Top working directory on worker node's local disk. The batch
# job working directory will be created below this
TOPWORKDIR=/scratch/`whoami`

# Basename of job sandbox (job workdir will be $TOPWORKDIR/$JOBDIR)
JOBDIR=jobdir

SEUSERSUBDIR=srmdir

#exe=/shome/pablom/CMSSW_3_8_4_patch2/src/DiLeptonAnalysis/NTupleProducer/macros/RunJZBAnalyzer
exe=exefile
##################################################################

ARGUMENTS=argumentsname

RELEASE=therelease

############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N example_job

# Run time soft and hard limits hh:mm:ss
# soft=CPU time; hard=Wallclock time
#$ -l s_rt=12:00:00,h_rt=14:00:00

### Specify the queue on which to run
#$ -q all.q

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
#  #$ -o /shome/username/mydir/
#  #$ -e /shome/username/mydir/

##################################################################



##### MONITORING/DEBUG INFORMATION ###############################
DATE_START=`date +%s`
echo "Job started at " `date`
cat <<EOF
################################################################
## QUEUEING SYSTEM SETTINGS:
HOME=$HOME
USER=$USER
JOB_ID=$JOB_ID
JOB_NAME=$JOB_NAME
HOSTNAME=$HOSTNAME
TASK_ID=$TASK_ID
QUEUE=$QUEUE

EOF

echo "################################################################"

if test 0"$DBG" -gt 0; then
   echo "######## Environment Variables ##########"
   env
   echo "################################################################"
fi


##### SET UP WORKDIR AND ENVIRONMENT ######################################
STARTDIR=`pwd`
WORKDIR=$TOPWORKDIR/$JOBDIR
RESULTDIR=$STARTDIR/$JOBDIR
if test x"$SEUSERSUBDIR" = x; then
   SERESULTDIR=$USER_SRM_HOME/$HN_NAME/$JOBDIR
else
   SERESULTDIR=$USER_SRM_HOME/$HN_NAME/$SEUSERSUBDIR
fi
if test -e "$WORKDIR"; then
   echo "ERROR: WORKDIR ($WORKDIR) already exists! Aborting..." >&2
   exit 1
fi
mkdir -p $WORKDIR
if test ! -d "$WORKDIR"; then
   echo "ERROR: Failed to create workdir ($WORKDIR)! Aborting..." >&2
   exit 1
fi

cd $WORKDIR
cat <<EOF
################################################################
## JOB SETTINGS:
STARTDIR=$STARTDIR
WORKDIR=$WORKDIR
RESULTDIR=$RESULTDIR
SERESULTDIR=$SERESULTDIR
EOF

###########################################################################
## YOUR FUNCTIONALITY CODE GOES HERE
# set up CMS environment
export SCRAM_ARCH=slc5_amd64_gcc434
source $VO_CMS_SW_DIR/cmsset_default.sh

cd $RELEASE

eval `scramv1 runtime -sh`

cd $WORKDIR

export LD_LIBRARY_PATH=/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH

# Output goes to myout, errors to myerr
$exe -o $SEOUTFILES $ARGUMENTS $SOURCEFILES > myout_job_$1.txt 2>myerr_job_$1.txt

# Here we produce some additional output for the files that are copied back to
#  our shared home
scramv1 list >> myout_job_$1.txt 2>>myerr_job_$1.txt
srmls srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/ >> myout_job_$1.txt 2>>myerr_job_$1.txt

# create a dummy file for copying back to the SE
dd if=/dev/urandom of=mybigfile count=100 &>/dev/null





#### RETRIEVAL OF OUTPUT FILES AND CLEANING UP ############################
cd $WORKDIR
if test 0"$DBG" -gt 0; then
    echo "########################################################"
    echo "############# Working directory contents ###############"
    echo "pwd: " `pwd`
    ls -Rl
    echo "########################################################"
    echo "YOUR OUTPUT WILL BE MOVED TO $RESULTDIR"
    echo "########################################################"
fi

if test x"$OUTFILES" != x; then
   mkdir -p $RESULTDIR
   if test ! -e "$RESULTDIR"; then
          echo "ERROR: Failed to create $RESULTDIR ...Aborting..." >&2
          exit 1
   fi
   for n in $OUTFILES; do
       if test ! -e $WORKDIR/$n; then
          echo "WARNING: Cannot find output file $WORKDIR/$n. Ignoring it" >&2
       else
          cp -a $WORKDIR/$n $RESULTDIR/$n
          if test $? -ne 0; then
             echo "ERROR: Failed to copy $WORDFIR/$n to $RESULTDIR/$n" >&2
          fi
   fi
   done
fi

if test x"$SEOUTFILES" != x; then
   if test 0"$DBG" -ge 2; then
      srmdebug="-debug"
   fi
   for n in $SEOUTFILES; do
       if test ! -e $WORKDIR/$n; then
          echo "WARNING: Cannot find output file $WORKDIR/$n. Ignoring it" >&2
       else
          srmcp $srmdebug -2 file:///$WORKDIR/$n $SERESULTDIR/$n
          if test $? -ne 0; then
             echo "ERROR: Failed to copy $WORKDIR/$n to $SERESULTDIR/$n" >&2
          fi
   fi
   done
fi

echo "Cleaning up $WORKDIR"
rm -rf $WORKDIR

###########################################################################
DATE_END=`date +%s`
RUNTIME=$((DATE_END-DATE_START))
echo "################################################################"
echo "Job finished at " `date`
echo "Wallclock running time: $RUNTIME s"
exit 0

