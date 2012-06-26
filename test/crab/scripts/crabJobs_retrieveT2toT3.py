#!/usr/bin/python
import sys, os, pwd, commands, shlex, glob

# check the command line arguments
if (len(sys.argv) != 2):      # if no argument is specified,
	print 'Usage:' 
	print ' python crabJobs_retrieveT2toT3.py <task directory>'
	print ' e.g.'
	print ' 	python crabJobs_retrieveT2toT3.py DYJetsToLL-Summer12-PU_S7_START52_V9-v1'
	sys.exit()
else:
   jobName = sys.argv[1]


# We need to extract the output datasetname from the crab directory
crabLogFile = jobName+'/log/crab.log'
if (not os.path.exists(crabLogFile)):
   print '*** ERROR: crab log file <'+crabLogFile+'> not found! Stopping here.'
   sys.exit(-1)
else:
   for line in open(crabLogFile).readlines():
      if "User Dataset Name" in line:
         oDatasetPath = line.split('=')[1].splitlines()[0]
         break

# We NEED this information!
if ( not len(oDatasetPath)>0 ):
   print 'Couldn\'t find dataset name! Stopping here.'
   sys.exit(-1)

# Remove leading white spaces
oDatasetPath = oDatasetPath.lstrip()

command = '/swshare/psit3/bin/dbs_transferRegister.py --dbs=ph02 --to-site=T3_CH_PSI --retransfer '+oDatasetPath

print 'Running',command
return_value = os.system(command)

if return_value != 0:
   print 'Problem in replicating data!'
   sys.exit(-1)

   
