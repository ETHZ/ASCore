#!/usr/bin/python
import sys, os, commands, shlex

print sys.argv[1:]
# check the command line arguments
if (len(sys.argv[1:]) > 4) or len(sys.argv[1:]) < 2:        # if no number or filename is specified,
	print 'Usage:' 
	print ' python multicrabJobs_retrieveT2toT3.py DATAMC ntupleVersion [T2Name [multicrabConfig]]'
	print ' e.g.'
	print ' 	python multicrabJobs_retrieveT2toT3.py MC V0X-0Y-0Z'
	print ' 	(the default T2 is T2_CH_CSCS)'
	print ' 	python multicrabJobs_retrieveT2toT3.py MC V0X-0Y-0Z T2_DE_RWTH'
	print ' 	(the default multicrab config file is multicrab.cfg)'
	print ' 	python multicrabJobs_retrieveT2toT3.py MC V0X-0Y-0Z T2_DE_RWTH multicrab.cfg'
	sys.exit()
else:
	if len(sys.argv[1:]) == 4:                              # if both T2Name and multicrab file are specified,
		datamc = sys.argv[1]
		ntupleVersion = sys.argv[2]
		T2Name = sys.argv[3]
		multicrabConfigName = sys.argv[4]
	else:
		if len(sys.argv[1:]) == 3:                              # if only T2Name is specified, use the default multicrab.cfg file
			datamc = sys.argv[1]
			ntupleVersion = sys.argv[2]
			T2Name = sys.argv[3]
			multicrabConfigName = "multicrab.cfg"
		else:                                                   # otherwise, use the default T2_CH_CSCS and multicrab.cfg file
			datamc = sys.argv[1]
			ntupleVersion = sys.argv[2]
			T2Name = "T2_CH_CSCS"
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
dirName = 'filelists_'+datamc+'_'+ntupleVersion
os.system('rm -rfv '+dirName)
os.system('mkdir -v -p '+dirName)
for jobName in jobNames:
	command = "../scripts/pyT3 ../scripts/crabJobs_retrieveT2toT3.py \""+jobName+" "+datamc+" "+ntupleVersion+" "+T2Name+"\""
	return_value = os.system(command)
	if return_value != 0:
		print 'Problem in processing task '+jobName+'. Continuing...'
