#!/usr/bin/python
import sys, os, commands, shlex, time

# check the command line arguments
if (len(sys.argv) > 1):        # if no number or filename is specified,
        print 'Usage:' 
	print ' python multicrabJobs_retrieveT2toT3.py [multicrabConfig]'
	sys.exit()
else:
        if len(sys.argv[1:]) == 1:                              # if both T2Name and multicrab file are specified,
                multicrabConfigName = sys.argv[1]
	else:
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
	command = "../scripts/pyT3 ../scripts/crabJobs_retrieveT2toT3.py \""+jobName+"\""
	return_value = os.system(command)
	if return_value != 0:
		print 'Problem in processing task '+jobName+'. Continuing...'
        time.sleep(1) # Stagger process to be safe
