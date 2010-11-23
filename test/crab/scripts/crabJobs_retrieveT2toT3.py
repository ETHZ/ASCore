#!/usr/bin/python
import sys, os, pwd, commands, shlex, glob

print sys.argv[1:]
# check the command line arguments
if (len(sys.argv[1:]) < 3) or (len(sys.argv[1:]) > 4):      # if no argument is specified,
	print 'Usage:' 
	print ' python crabJobs_retrieveT2toT3.py jobName DATAMC ntupleVersion [T2NAME]'
	print ' e.g.'
	print ' 	python crabJobs_retrieveT2toT3.py NTupleProducer MC V0X-0Y-0Z'
	print ' 	(the default T2 is T2_CH_CSCS)'
	print ' 	python crabJobs_retrieveT2toT3.py NTupleProducer MC V0X-0Y-0Z T2_DE_RWTH'
	sys.exit()
else:
	if len(sys.argv[1:]) == 4:                              # if T2NAME is specified,
		jobName = sys.argv[1]
		datamc = sys.argv[2]
		ntupleVersion = sys.argv[3]
		T2NAME = sys.argv[4]
	else:                                                   # otherwise, use the T2_CH_CSCS
		jobName = sys.argv[1]
		datamc = sys.argv[2]
		ntupleVersion = sys.argv[3]
		T2NAME = "T2_CH_CSCS"

# initialize SRM variables
USERNAME=pwd.getpwuid(os.getuid())[0]
jobDir="ntuples/"+datamc+"/"+ntupleVersion+"/"+jobName
T3SRMDIR="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/susy/"+jobDir
T3SRM="srm://t3se01.psi.ch:8443/srm/managerv2?SFN="
SRMOPTIONS="-count=1000"
if ("T2_DE_RWTH" in T2NAME):
	T2SRM="srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="
	T2SRMDIR="srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/"+USERNAME+"/"+jobDir
	nameSplitString="cms"
else:
	T2SRM="srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN="
	T2SRMDIR="srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/"+USERNAME+"/"+jobDir
	nameSplitString="trivcat"

# create the list of files that need to be retrieved from T2 to T3
filelistName = "filelists_"+datamc+"_"+ntupleVersion+"/filelist_"+jobName+".txt"
command="srmls "+T2SRMDIR+" | awk '{split($2,nameParts,\""+nameSplitString+"\"); if (NR>3) print nameParts[2]}' > " + filelistName
print "command: " + command
return_value = os.system(command)
if return_value != 0:
	print 'Problem to determine which jobs need to be retrieved. Exiting...'
	sys.exit()

# run data_replica
command="data_replica.py --debug "+filelistName+" --from-site "+T2NAME+" --to-site T3_CH_PSI /store/user/susy/"+jobDir
#print "command: " + command
return_value = os.system(command)
if return_value == 0:
	print 'Data replicated succesfully. Exiting...'
	sys.exit()
else:
	print 'Problem in replicating data with "data_replica". Trying directly with "lcg-cp"...'

# try to do simple lcg-cp from T2 to T3 in case there was problem with data_replica
numJobs=0
command_copy = 'lcg-cp -v -n 1 -T srmv2 -U srmv2'           # for copying from remote T2 site
while ((numJobs == 0) or (numJobs == 1000) or (numJobs == 2000)): 
	command="srmls -offset="+str(numJobs)+" "+SRMOPTIONS+" "+T2SRMDIR
	print 'Getting the list of '+jobName+' jobs at SE (from job '+str(numJobs+1)+' to max.'+str(numJobs+1000)+')...'
	status, output = commands.getstatusoutput(command)
	filenames = shlex.split(output)                         # get list of files in the jobName directory at the SE path (and their sizes)
	print 'There are ' + str((len(filenames)-2)/2) + ' jobs to be retrieved.'
	if len(filenames) == 2:                                 # get out of the loop if there are no more jobs
		break
	for filename in filenames:
		if (jobName in filename ) and ('.root' in filename) and not ((jobName+'.root') == filename):
			nameParts=filename.split(jobName+'/')
			command = command_copy + ' ' + T2SRMDIR+'/'+nameParts[1] + ' ' + T3SRMDIR+'/'+nameParts[1]
#			print 'with command: ' + command
			print 'copying from SE job ' + str(numJobs+1) + ':'
			return_value = os.system(command)
			numJobs = numJobs + 1
	print "numJobs: " + str(numJobs)
	if numJobs == 0:										# get out of the loop if there are no more jobs
		break

if numJobs > 0:												# print info on num of retrieved jobs
	print str(numJobs) + ' jobs have been processed.'
else:
	print '0 jobs have been retrieved. Exiting...'
	sys.exit()
