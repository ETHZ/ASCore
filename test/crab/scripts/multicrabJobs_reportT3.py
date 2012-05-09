#!/usr/bin/python
import sys, os, commands, shlex

# check the command line arguments
if (len(sys.argv[1:]) > 3) or len(sys.argv[1:]) < 2:        # if no number or filename is specified,
	print 'Usage:' 
	print ' python multicrabJobs_reportT3.py datamc ntupleVersion [multicrabConfig]'
	print ' interactively as e.g.:'
	print ' 	python multicrabJobs_reportT3.py mc V0X-0Y-0Z'
	print ' 	(the default multicrab config file is multicrab.cfg)'
	print ' 	python multicrabJobs_reportT3.py mc V0X-0Y-0Z multicrab.cfg'
	print ' or submit it to the qsub as e.g.:'
	print ' 	../scripts/pyT3 ../scripts/multicrabJobs_reportT3.py "mc V0X-0Y-0Z"'
	print ' 	(the default multicrab config file is multicrab.cfg)'
	print ' 	../scripts/pyT3 ../scripts/multicrabJobs_reportT3.py "mc V0X-0Y-0Z multicrab.cfg"'
	sys.exit()
else:
	if len(sys.argv[1:]) == 3:                              # if only T2Name is specified, use the default multicrab.cfg file
		datamc = sys.argv[1]
		ntupleVersion = sys.argv[2]
		multicrabConfigName = sys.argv[3]
	if len(sys.argv[1:]) == 2:                              # otherwise, use the default T2_CH_CSCS and multicrab.cfg file
		datamc = sys.argv[1]
		ntupleVersion = sys.argv[2]
		multicrabConfigName = "multicrab.cfg"

# check if multicrab config file exists
if (not os.path.exists(multicrabConfigName)):
	print 'Multicrab config file <'+multicrabConfigName+'> does not exist. Exiting...'
	sys.exit()

# get the list of all multicrab tasks that need to be processed
command_getMulticrabJobs = "grep \"\[\" "+multicrabConfigName+" | awk '{if (!index($1,\"#\")) {split($1,strarr,\"[\"); split(strarr[2],strarr,\"]\"); if (strarr[1]!=\"COMMON\" && strarr[1]!=\"MULTICRAB\") print strarr[1]}}'"
status, output = commands.getstatusoutput(command_getMulticrabJobs)
jobNames = shlex.split(output)								# get list of tasks
if status != 0:
	print 'Problem to determine which multicrab tasks need to be processed. Exiting...'
	sys.exit()

# process all tasks with crabJobs_retrieveT2toT3.py command
for jobName in jobNames:
	command = "python ../scripts/crabJobs_reportT3.py "+jobName+" "+datamc+" "+ntupleVersion
	return_value = os.system(command)
	if return_value != 0:
		print '\nProblem in processing report for the task '+jobName+'. Continuing...'
