#!/usr/bin/python
import sys, os, shlex, pwd, commands
import re
from optparse import OptionParser


def get_info(lumiFile,jsonFile):
   """Retrieve miscellaneous task information to be put on the twiki"""
   # Run range from json file
   command_print = 'printJSON.py --range '+jsonFile
   status, range = commands.getstatusoutput(command_print)
   if status != 0:
      print 'Problem with getting range of JSON file',jsonFile,':',status
      range = '... - ...'
   else:
      # Strip printJSON.py's output
      range = range.lstrip('runs ')
      range = range.replace(' ','')
        
   # Luminosity from lumi summary
   tbegin = 0 # tag beginning of totals
   pat = re.compile('(\|\s*[\d\.]+\s*)+\|')
   lumiTot = 0
   for line in open(lumiFile):
      if line.find('Total') != -1: tbegin = 1
      elif tbegin:
         m = re.match(pat,line) 
         if m:
            lumiTot = m.group(1)  # Strangely only stores last match...
            lumiTot = lumiTot.strip('| ')
            break

   return [range,lumiTot]

  

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
   elif len(sys.argv[1:]) == 4:                              # otherwise, use the specifed multicrab cfg file
      jobName = sys.argv[1]
      datamc = sys.argv[2]
      ntupleVersion = sys.argv[3]
      multicrabConfigName = sys.argv[4]

# initialize SRM variables
USERNAME=os.getenv('USER')
if (USERNAME == "pnef"):
   userNickName = "Pascal"
elif (USERNAME == "pablom"):
   userNickName = "Pablo"
elif (USERNAME == "thea"):
   userNickName = "Alessandro"
elif (USERNAME == "predragm"):
   userNickName = "Pedja"
elif (USERNAME == "buchmann"):
   userNickName = "Marco-Andrea"
elif (USERNAME == "fronga"):
   userNickName = "Frederic"
else:
   userNickName = "Great guy!"
jobDir="ntuples/"+datamc+"/"+ntupleVersion+"/"+jobName
T3SRMDIR="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/susy/"+jobDir
T3SRM="srm://t3se01.psi.ch:8443/srm/managerv2?SFN="

# Filenames
fJsonSummaryCrab = jobName+'/res/lumiSummary.json'
fJsonSummary = jobName+'-'+ntupleVersion+'_lumiSummary.json.txt'
fLumiSummary = jobName+'-'+ntupleVersion+'_lumi.txt'

# initialize basic data/mc commands
command_status = 'crab -c '+jobName+' -status'
command_getoutput = 'crab -c '+jobName+' -getoutput'
command_report = 'crab -c '+jobName+' -report'
command_lumiCalc = 'lumiCalc2.py -i '+fJsonSummaryCrab+' --nowarning overview > '+fLumiSummary
command_lumiSumm = 'mv '+fJsonSummaryCrab+' '+fJsonSummary
command_getSize = 'awk -F\\" \'$0 ~ "Timing-file-write-totalMegabytes" {s+=$4} END {print s}\' '+jobName+'/res/crab_fjr_*.xml'
#"srmls "+T3SRMDIR+" | awk 'BEGIN{size=0}{size+=$1}END{printf( \"Size(GB): %5.1f , Entries: %d\",size/1e9,NR-2)}'"

command_getEvents = 'awk \'$0 ~ "EventsRead" { getline; s+=$1; getline} END {print s}\' '+jobName+'/res/crab_fjr_*.xml'
#"python ../scripts/GetNEvents_CrabReport.py "+jobDir

# check if multicrab config file exists and get the dataset name
datasetFound = 0
if (not os.path.exists(multicrabConfigName)):
   print 'Multicrab config file <'+multicrabConfigName+'> does not exist...'
   datasetName = "\MULTICRAB\CONFIG\NOT\FOUND"
else:
   datasetName = "\NOT\FOUND\IN\MULTICRAB\CONFIG"
   for line in open(multicrabConfigName):
      if jobName in line:
         datasetFound = 1
      if "CMSSW.datasetpath" in line and datasetFound:
         datasetName = line.split("=")[1].splitlines()[0]
         break


# # get the dataset size on SE
print '--> Getting size of dataset...\n   ',command_getSize
status, datasetSize = commands.getstatusoutput(command_getSize)
if status != 0:
   print 'Problem with getting size of the job '+jobName+'. Exiting...'
   datasetSize = "XX.X"
datasetSize =  "%3.1fGB" % (float(datasetSize)/1024)
print 'datasetSize =',datasetSize

# get number of events ib the dataset
print '--> Getting number of events in dataset...\n   ',command_getEvents
status, datasetEvents = commands.getstatusoutput(command_getEvents)
if status != 0:
   print 'Problem with getting number of events for the job '+jobName+'. Exiting...'
   datasetEvents = "XX.X"
datasetEvents = "%3.1fM" % (float(datasetEvents)/1e6)
print 'datasetEvents =',datasetEvents

# Get CMSSW version
cmsswVersion = os.getenv('CMSSW_VERSION')[6:] # Remove leading "CMSSW_"
print 'cmsswVersion =',cmsswVersion

# perform crab tasks for data jobs
if (datamc=="data"):
   print '--> Getting status of all jobs\n   ',command_status
   return_value,output = commands.getstatusoutput(command_status)
   if return_value != 0:
      print 'Problem with crab -status of the job '+jobName+'. Exiting...'
      print output
      sys.exit()

   print '--> Getting output of all jobs...\n   ',command_getoutput
   return_value,output = commands.getstatusoutput(command_getoutput)
   if return_value != 0:
      print 'Problem with crab -getoutput of the job '+jobName+'. Continuing anyway...'
      print output
      #sys.exit()

   print '--> Getting crab report...\n   ',command_report
   return_value,output = commands.getstatusoutput(command_report)
   if return_value != 0:
      print 'Problem with crab -report of the job '+jobName+'. Exiting...'
      print output
      sys.exit()

   print '--> Running luminosity calculation...\n   ',command_lumiCalc
   return_value,output = commands.getstatusoutput(command_lumiCalc)
   if return_value != 0:
      print 'Problem with running lumiCalc2.py of the job '+jobName+'. Exiting...'
      print output
      sys.exit()

   return_value = os.system(command_lumiSumm)
   if return_value != 0:
      print 'Problem with moving/creating json report of the job '+jobName+'. Exiting...'
      sys.exit()

   # get misc. information to be put on the twiki
   [range,lumi] = get_info(fLumiSummary,fJsonSummary)

   format_Twiki = '|[[%ATTACHURL%/'+fJsonSummary+']['+range+']] |'+datasetName+' |'+datasetSize+' |'+datasetEvents+' |[[%ATTACHURL%/'+fLumiSummary+']['+lumi+'/pb]] | %TWISTY{showlink="Show..." hidelink="Hide"}%<br>/store/user/susy/'+jobDir+'/ <br>%ENDTWISTY% |'+cmsswVersion+' |'+ntupleVersion+' |'+userNickName+' | |'
else:
   format_Twiki = "| "+jobName+" | "+datasetName+" | "+datasetSize+" GB | "+datasetEvents+' | ... pb | %TWISTY{showlink="Show..." hidelink="Hide"}%<br>/store/user/susy/'+jobDir+'/ <br>%ENDTWISTY% | '+cmsswVersion+' | '+ntupleVersion+" | "+userNickName+" | |"

print '---------------------------------------\nAdd the following line to the Twiki:\n\n',format_Twiki
