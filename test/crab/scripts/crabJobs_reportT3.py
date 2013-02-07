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
   patU = re.compile('.*Recorded\((.*)\).*')
   lumiTot = 0
   units = ''
   for line in open(lumiFile):
      # Scroll down to totals
      if line.find('Total') != -1: tbegin = 1
      elif tbegin:
         m = re.match(patU,line)
         if m: 
            units = m.group(1)
            continue
         m = re.match(pat,line) 
         if m:
            lumiTot = m.group(1)  # Strangely only stores last match...
            lumiTot = lumiTot.strip('| ')
            break

   lumiString =  "%3.1f%s" % (float(lumiTot), units)
   print 'lumiString:',lumiString,lumiFile
   return [range,lumiString]

  

# check the command line arguments
if (len(sys.argv[1:]) != 3):                                # if not all three arguments are specified,
   print 'Usage:' 
   print ' python crabJobs_reportT3.py taskName datamc ntupleVersion'
   print ' e.g.'
   print ' 	python crabJobs_reportT3.py DYJetsToLL  mc V0X-0Y-0Z'
   sys.exit()
else:
   jobName = sys.argv[1].rstrip('/')
   datamc = sys.argv[2]
   ntupleVersion = sys.argv[3]

# Intialization: we'll get essential information from crab.log!
crabLogFile = jobName+'/log/crab.log'

# Usernames we know
USERNAME=os.getenv('USER')
if (USERNAME == "pnef"):
   userNickName = "Pascal"
elif (USERNAME == "pablom"):
   userNickName = "Pablo"
elif (USERNAME == "thea"):
   userNickName = "Alessandro"
elif (USERNAME == "buchmann"):
   userNickName = "Marco-Andrea"
elif (USERNAME == "fronga"):
   userNickName = "Frederic"
elif (USERNAME == "mdunser"):
   userNickName = "Marc"
elif (USERNAME == "mmasciov"):
   userNickName = "Super Mario"
elif (USERNAME == "peruzzi"):
   userNickName = "Marco"
else:
   userNickName = "Great guy!"

# Filenames
fJsonSummaryCrab = jobName+'/res/lumiSummary.json'
fJsonSummary     = jobName+'-'+ntupleVersion+'_lumiSummary.json.txt'
fLumiSummary     = jobName+'-'+ntupleVersion+'_lumi.txt'

# initialize basic data/mc commands
command_report = 'crab -c '+jobName+' -report'
command_lumiCalc = 'lumiCalc2.py -i '+fJsonSummaryCrab+' --nowarning overview > '+fLumiSummary
command_lumiSumm = 'mv '+fJsonSummaryCrab+' '+fJsonSummary


# get the input and output dataset names from the crab.cfg file
iDatasetPath = ''
iPrimaryDataset = ''
oDatasetPath = ''

if (not os.path.exists(crabLogFile)):
   print '*** ERROR: crab log file <'+crabLogFile+'> not found! Stopping here.'
   sys.exit(-1)
else:
   for line in open(crabLogFile).readlines():
      if "CMSSW.datasetpath :" in line:
         iDatasetPath = line.split(':')[1].splitlines()[0]
         iPrimaryDataset = iDatasetPath.split('/')[1]
      if "User Dataset Name" in line:
         oDatasetPath = line.split('=')[1].splitlines()[0]
      if ( len(iDatasetPath)>0 and len(oDatasetPath)>0 ):
         break

# We NEED this information!
if ( not (len(iDatasetPath)>0 and len(oDatasetPath)>0) ):
   print 'Couldn\'t find dataset names. Stopping here.'
   sys.exit(-1)

# Remove leading white spaces
iDatasetPath = iDatasetPath.lstrip()
iPrimaryDataset = iPrimaryDataset.lstrip()
oDatasetPath = oDatasetPath.lstrip()

# Get CMSSW version
cmsswVersion = os.getenv('CMSSW_VERSION')[6:] # Remove leading "CMSSW_"
dasURL = 'https://cmsweb.cern.ch/das/request?view=list&limit=10&instance=cms_dbs_ph_analysis_02&input=dataset+dataset%3D'+oDatasetPath

# perform crab tasks for data jobs
if (datamc=="data"):
   print '--> Getting crab report... (this can take a while)\n   ',command_report
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

   format_Twiki = '| [[%ATTACHURL%/'+fJsonSummary+']['+range+']]  | '+iDatasetPath+' | [[%ATTACHURL%/'+fLumiSummary+']['+lumi+']]  | [['+dasURL+'][DAS]]  | '+cmsswVersion+'  | '+ntupleVersion+'  | '+userNickName+'  ||'
else:
   format_Twiki = '| '+iPrimaryDataset+'  | '+iDatasetPath+'  | [['+dasURL+'][DAS]]  | ...pb | '+cmsswVersion+'  | '+ntupleVersion+'  | '+userNickName+'  ||'
   format_Twiki += '\n\nAND DO NOT FORGET TO FILL IN THE CROSS-SECTION INFORMATION!\n'

print '---------------------------------------\nTask '+jobName+':\nAdd the following line to the Twiki:\n\n',format_Twiki
