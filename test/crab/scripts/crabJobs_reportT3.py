#!/usr/bin/python
import sys, os, shlex, pwd, commands

# check the command line arguments
if (len(sys.argv[1:]) < 3 or len(sys.argv[1:]) > 4):                                # if not all three arguments are specified,
	print 'Usage:' 
	print ' python crabJobs_reportT3.py jobName datamc ntupleVersion [multicrabConfig]'
	print ' e.g.'
	print ' 	python crabJobs_reportT3.py NTupleProducer mc V0X-0Y-0Z'
	print ' 	(the default multicrab config file is multicrab.cfg)'
	print ' 	python crabJobs_reportT3.py NTupleProducer mc V0X-0Y-0Z multicrab.cfg'
	sys.exit()
else:
	if len(sys.argv[1:]) == 3:                              # use the default multicrab.cfg file
		jobName = sys.argv[1]
		datamc = sys.argv[2]
		ntupleVersion = sys.argv[3]
		multicrabConfigName = "multicrab.cfg"
	if len(sys.argv[1:]) == 4:                              # otherwise, use the specifed multicrab cfg file
		jobName = sys.argv[1]
		datamc = sys.argv[2]
		ntupleVersion = sys.argv[3]
		multicrabConfigName = sys.argv[4]

# initialize SRM variables
USERNAME=pwd.getpwuid(os.getuid())[0]
if (USERNAME == "pnef"):
	userNickName = "Pascal"
elif (USERNAME == "pablom"):
	userNickName = "Pablo"
elif (USERNAME == "thea"):
	userNickName = "Alessandro"
elif (USERNAME == "predragm"):
	userNickName = "Pedja"
else:
	userNickName = "Great guy!"
jobDir="ntuples/"+datamc+"/"+ntupleVersion+"/"+jobName
T3SRMDIR="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/susy/"+jobDir
T3SRM="srm://t3se01.psi.ch:8443/srm/managerv2?SFN="
# initialize basic data/mc commands
command_status = 'crab -c '+jobName+' -status'
command_getoutput = 'crab -c '+jobName+' -getoutput'
command_report = 'crab -c '+jobName+' -report'
command_lumiCalc = 'lumiCalc.py -i '+jobName+'/res/lumiSummary.json --nowarning overview > '+jobName+'-'+ntupleVersion+'_lumi.txt'
command_lumiSumm = 'mv '+jobName+'/res/lumiSummary.json '+jobName+'-'+ntupleVersion+'_lumiSummary.json.txt'
command_getSize = "srmls "+T3SRMDIR+" | awk 'BEGIN{size=0}{size+=$1}END{printf( \"Size(GB): %5.1f , Entries: %d\",size/1e9,NR-2)}'"
command_getEvents = "python ../scripts/GetNEvents_CrabReport.py "+jobDir

# check if multicrab config file exists and get the dataset name
datasetFound = 0
if (not os.path.exists(multicrabConfigName)):
	print 'Multicrab config file <'+multicrabConfigName+'> does not exist...'
	datasetName = "\MULTICRAB\CONFIG\NOT\FOUND"
else:
	for line in open(multicrabConfigName):
		if jobName in line:
			datasetFound = 1
		if "CMSSW.datasetpath" in line and datasetFound:
			datasetName = line.split("=")[1].splitlines()[0]
			break
	else:
		datasetName = "\MULTICRAB\CONFIG\NOT\FOUND"

# get the dataset size on SE
status, output = commands.getstatusoutput(command_getSize)
if status != 0:
	print 'Problem with getting size of the job '+jobName+'. Exiting...'
	datasetSize = "XX.X"
else:
	datasetSize = shlex.split(output)[1]

# get number of events ib the dataset
status, output = commands.getstatusoutput(command_getEvents)
if status != 0:
	print 'Problem with getting number of events for the job '+jobName+'. Exiting...'
	datasetEvents = "XX.X"
else:
	datasetEvents = shlex.split(output)[0]

# perform crab tasks for data jobs
if (datamc=="data"):
	return_value = os.system(command_status)
	if return_value != 0:
		print 'Problem with crab -status of the job '+jobName+'. Exiting...'
		sys.exit()

	return_value = os.system(command_getoutput)
	if return_value != 0:
		print 'Problem with crab -getoutput of the job '+jobName+'. Exiting...'
		#sys.exit()

	return_value = os.system(command_report)
	if return_value != 0:
		print 'Problem with crab -report of the job '+jobName+'. Exiting...'
		sys.exit()

	return_value = os.system(command_lumiCalc)
	if return_value != 0:
		print 'Problem with running lumiCalc.py of the job '+jobName+'. Exiting...'
		sys.exit()

	return_value = os.system(command_lumiSumm)
	if return_value != 0:
		print 'Problem with moving/creating json report of the job '+jobName+'. Exiting...'
		sys.exit()
	format_Twiki = "| [[%ATTACHURL%/"+jobName+"-"+ntupleVersion+"_lumiSummary.json.txt][...-...]] | "+datasetName+" | "+datasetSize+" GB | ...M | [[%ATTACHURL%/"+jobName+"-"+ntupleVersion+"_lumi.txt][... /pb]] | %TWISTY{showlink=\"Show...\" hidelink=\"Hide\"}%<br>/store/user/susy/"+jobDir+"/ <br>%ENDTWISTY% | 3_8_6 | "+ntupleVersion+" | "+userNickName+" | |"
else:
	format_Twiki = "| "+jobName+" | "+datasetName+" | "+datasetSize+" GB | "+datasetEvents+" | ... pb | %TWISTY{showlink=\"Show...\" hidelink=\"Hide\"}%<br>/store/user/susy/"+jobDir+"/ <br>%ENDTWISTY% | 3_8_6 | "+ntupleVersion+" | "+userNickName+" | |"

print format_Twiki
