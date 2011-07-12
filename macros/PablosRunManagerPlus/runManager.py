#!/usr/bin/python

### WATCH OUT: NEEDS TO BE LAUNCHED AS python runManager.py configuration.cfg and NOT ./runManager configuration.cfg

###############################################################
#                       runManager.py                         #
###############################################################
# + Python script which allows to manage multiple jobs in the #
# context of a lsfbatch system.                               #
#                                                             #
# + type runManager -h to see a description of the usage      #   
#                                                             #
# + Feedback to: pablom@cern.ch                               #
# + Feedback about this version: contact marco-andrea         #
#                                                             #
############################################################### 
   

domultiprocessing=True

# System libraries to interact with the shell
import sys
import os
from os import popen
import commands
import time
try:
	from multiprocessing import Pool
except ImportError:
	print "Importing multiprocessing failed. Please use the script like this: python runManager.py template.cfg and not ./runManager.py template.cfg"
	domultiprocessing=False

jobnumbers=[]
totaljobnumber=0

#-------------------------------------------------------------#
# ensure_dir method: Make sure a path exists                  #
#-------------------------------------------------------------#

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

#-------------------------------------------------------------#
# checklist method: Check whether the jobs are still running  #
# to permit post-processing                                   #
#-------------------------------------------------------------#

def checklist(joblist,currlist):
  jobfinished=[]
  for i in joblist:
    if i in currlist: ## if the job is still active
      currlist.remove(int(i))## we don't need this anymore on the list of currently running jobs (we already know it's running)
    else:#job is finished
      jobfinished.append(int(i)) 
  for job in jobfinished:
    joblist.remove(int(job))
  print "Jobs left:"+str(len(joblist))+" / "+str(totaljobnumber)

#-------------------------------------------------------------#
# getFiles method: Obtains a list with the name of the files  #
# that will be used as input for a given job "step"           #        
# returns a string with the list                              # 
#-------------------------------------------------------------#
def getFiles(step, FilesPerJob, NumberOfJobs, OverloadedJobs, listOfFiles):
  
  if(step <= NumberOfJobs-OverloadedJobs):
    first = step*FilesPerJob
  else:
    first = FilesPerJob*(NumberOfJobs-OverloadedJobs)+(FilesPerJob+1)*(step-(NumberOfJobs-OverloadedJobs))-1

  counter = 0
  files = []
  for line in listOfFiles:
    if(step <= NumberOfJobs-OverloadedJobs):
      if(counter >= first and first+FilesPerJob > counter):
        files.append(line)
    else:
      if(counter >= first and first+FilesPerJob+1 > counter):
        files.append(line) 
    counter = counter + 1    
  
  return files
  #Check the ending in case it has a ',' 
  #if(files[-1][-2] == ','):
  #  files[-1] = files[-1][0:-2]
  #stringCad = ""
  #for fileLine in files:
  #  stringCad = stringCad + fileLine
  #return stringCad 
      

#-------------------------------------------------------------#
# createDirectory method: Creates a directory                 #
# returns 0 if everything was all right                       #
#-------------------------------------------------------------# 
def createDirectory(nameOfDirectory):
  result = os.system("mkdir " + nameOfDirectory)
  result2 = os.system("cd " + nameOfDirectory)
  return result 

#-------------------------------------------------------------#
# fixPath method: Just makes sure the path is ending          #
# with slash. Otherwise appends it                            #
#-------------------------------------------------------------# 
def fixPath(path):
  if(path[-1] == '/'):
    return path
  else:
    return path + "/"


#-------------------------------------------------------------#
# createCMSConf method: creates CMS configuration file        #
#                                                             #
#-------------------------------------------------------------# 
def createCMSConf(step, nameOfDirectory, releasePath, nameOfConf, inputString, executablePath, nameOfSRM, hpname, task):

  CFGFile = open(nameOfConf)
  cfgText = CFGFile.read()

  cadFile = "\""
  for file in inputString:
    cadFile = cadFile + file + " " 
  cadFile = cadFile + "\"" 

  nameOfFolder = task[0]
  if(nameOfFolder[len(nameOfFolder)-1] == "/"):
    nameOfFolder = nameOfFolder[0:len(nameOfFolder)-1]
  value = nameOfFolder.rfind("/")
  taskName = nameOfFolder[value+1:len(nameOfFolder)] 


  newcfgText1 = cfgText.replace("sourcefile" , cadFile)
  newcfgText2 = newcfgText1.replace("jobdir" , taskName + "/" + nameOfSRM + "_" + str(step))
  newcfgText3 = newcfgText2.replace("exefile", executablePath)
  newcfgText4 = newcfgText3.replace("srmdir", nameOfSRM + "/" + taskName)
  newcfgText5 = newcfgText4.replace("hpname", hpname)
  newcfgText6 = newcfgText5.replace("argumentsname", "\"" + task[1] + "\"")
  newcfgText7 = newcfgText6.replace("therelease", releasePath)

  nameOfConf2 = nameOfConf.replace(".", "_"+str(step)+ ".")
 
  outputCFGFile = open(nameOfDirectory + taskName + "/" + nameOfConf2, "w")
  outputCFGFile.write(newcfgText7)
  CFGFile.close()
  outputCFGFile.close()
  
  outputName = "output_" + str(step) + ".root"
  thisjobnumber=0
  pipe=popen("qsub -e /tmp/ -o /tmp/ -N " + "RMJ"+str(step)+taskName + " " + nameOfDirectory + taskName + "/" + nameOfConf2 + " " + str(step))
  for l in pipe.readlines():
	  if l.find("Your job")>-1:
		  thisjobnumber=int(l[l.index('job ')+4:l.index(' (')])
		  print str(taskName)+": Submitted job "+str(step)+" with job number: "+str(thisjobnumber)
  return thisjobnumber
###############################################################
# createJob method: Prepares everything for a single job      #
#                                                             #
############################################################### 
def createJob(step, FilesPerJob, NumberOfJobs, OverloadedJobs, stringInfo, listOfFiles, task):
 
  #Definition of the different substrings needed
  nameOfSourceFile = stringInfo[0]
  releasePath = fixPath(stringInfo[1])
  executablePath = stringInfo[2]
  nameOfConfigurationFile = stringInfo[3]
  nameOfSRMPath = stringInfo[4]
  nameOfCurrentDirectory = stringInfo[5]
  hpname = stringInfo[6]
  
  inputFilesForCMS = getFiles(step, FilesPerJob, NumberOfJobs, OverloadedJobs, listOfFiles)
  if(inputFilesForCMS == ""):
    return "Error: No input files available"
  

  result = createCMSConf(step, nameOfCurrentDirectory, releasePath, 
           nameOfConfigurationFile, inputFilesForCMS, executablePath, nameOfSRMPath, hpname, task)
  return result

###############################################################
# parseInputFile method: Parses the input file                #
# More options could be added if needed                       #
############################################################### 

def parseInputFile(name):   
  showMessage("Parsing " + name + " file ...")
  file = open(name)
  stringInfo = [[],[],[],[],[], [], []]
  for line in file.readlines():
    splitLine = line.split()
    if(splitLine[0] == "ReleasePath"):
      stringInfo[1] = splitLine[1]
    if(splitLine[0] == "Executable"):
      stringInfo[2] = splitLine[1]
    if(splitLine[0] == "NameOfConf"):
      stringInfo[3] = splitLine[1]
    if(splitLine[0] == "NameOfSource"):
      stringInfo[0] = splitLine[1]
    if(splitLine[0] == "SrmPath"):
      stringInfo[4] = splitLine[1]
    if(splitLine[0] == "WorkPath"):
      stringInfo[5] = splitLine[1]
    if(splitLine[0] == "HPName"):
      stringInfo[6] = splitLine[1]

  file.close()

  for i in range(0,len(stringInfo)):
    if(stringInfo[i] == []):
      return "Error"   
  showMessage("Configuration Parameters:")
  showMessage("ReleasePath = " + stringInfo[1] + sanitycheck("ReleasePath",stringInfo[1]))
  showMessage("Executable = " + stringInfo[2] + sanitycheck("Executable",stringInfo[2]))
  showMessage("NameOfConf = " + stringInfo[3] + sanitycheck("NameOfConf",stringInfo[3]))
  showMessage("NameOfSource = " + stringInfo[0] + sanitycheck("NameOfSource",stringInfo[0]))
  showMessage("SrmPath = " + stringInfo[4])
  showMessage("WorkPath = " + stringInfo[5] + sanitycheck("WorkPath",stringInfo[5]))
  showMessage("HPName = " + stringInfo[6])
  return stringInfo

  
###############################################################
# showMessage method: Adds label and print message            #
#                                                             #
############################################################### 
def showMessage(message):
  print "[runManager] " + message


###############################################################
# printUsage method: Prints how to use the program            #
#                                                             #
###############################################################
def printUsage():
  print "./runMan configuration.txt"
  print "configuration.txt should look like this:"
  print "---------------------------------------"
  print "ReleasePath path"
  print "Executable name"
  print "NameOfConf path"
  print "NameOfSource path"
  print "SrmPath path"
  print "WorkPath path"
  print "HPName path"
  print "---------------------------------------"
 
###############################################################
# List of tasks: from the source file                         #
#                                                             #
###############################################################
def getListOfTasks(nameOfFile):

  file = open(nameOfFile)
  
  tasks = []

  for line in file.readlines():
    if(line[0] == "#"):
      continue
    splitedLine = line.split()
    argumentBegin = line.find("'")
    argumentEnd = line.rfind("'")
	    
    if(argumentBegin == -1 or argumentEnd == -1 or argumentBegin == argumentEnd):
      showMessage("Ignoring the line : " + line)
    else:  
      arguments = line[argumentBegin+1:argumentEnd]
      endLine = line[argumentEnd+1:len(line)]
      splitedEndLine = endLine.split()
      folder = splitedLine[0] 
      numberOfJobs = splitedEndLine[0]
      numberOfFilesPerJob = splitedEndLine[1]
      task = [folder, arguments, numberOfJobs, numberOfFilesPerJob]
      if(numberOfJobs == '-1' or numberOfFilesPerJob == '-1'):
        tasks.append(task)

  return tasks  


###############################################################
# process: process a single task                              #
#                                                             #
###############################################################
def process(task, conf):
 
  dcapPath = "dcap://t3se01.psi.ch:22125/"
  srmPath = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/"
  folderToList = srmPath + task[0]
  showMessage("Going to fetch list of files pertaining to "+str(task[0]))
  theList = os.popen("lcg-ls " + folderToList)
  list = theList.readlines()
  for li in theList:
    if(li.find("not found") != -1):
      showMessage("Folder not found")
      return "Error"
  
  numberOfFiles = len(list)
  if(numberOfFiles == 0):
    showMessage("No files found")
    return "Error"


  correctList = [];
  for fileName in list:
    auxiliar = dcapPath + fileName
    auxiliar = auxiliar[0:len(auxiliar)-1]
    correctList.append(auxiliar)

  if(task[2] != "-1"):
    FilesPerJob = int(float(numberOfFiles)/int(task[2]))
    OverloadedJobs = numberOfFiles + 1 - FilesPerJob * int(task[2])
    NumberOfJobs = int(task[2])
  else:
    FilesPerJob = int(task[3])
    NumberOfJobs = int(numberOfFiles/FilesPerJob) + 1
    OverloadedJobs = numberOfFiles + 1 - FilesPerJob * NumberOfJobs  
    
  if(numberOfFiles<int(task[2])) : 
    FilesPerJob = 1;
    NumberOfJobs = numberOfFiles
    OverloadedJobs = 0

  nameOfFolder = task[0]
  if(nameOfFolder[len(nameOfFolder)-1] == "/"):
    nameOfFolder = nameOfFolder[0:len(nameOfFolder)-1]
  value = nameOfFolder.rfind("/")
  taskName = nameOfFolder[value+1:len(nameOfFolder)] 
  os.system("mkdir -p " + taskName)
  os.system("mkdir -p ../" + taskName)

  showMessage(str(NumberOfJobs) + " with " + str(FilesPerJob) + " each will be created")
  jobnumbercollection=[]
  if(FilesPerJob > 0):
    for step in range(0, NumberOfJobs):
      result=createJob(step, FilesPerJob, NumberOfJobs, OverloadedJobs, conf, correctList, task)
      jobnumbercollection.append(result)
  return jobnumbercollection

##################################################################################
#                              POST PROCESSING                                   #
##################################################################################

def cb(r): #optional: callback function
    for listentry in r:
      jobnumbers.append(listentry)
    pass

def ensure_dir(f) :
	if not os.path.exists(f):
		os.makedirs(f)

def join_directory(path,filelist,username) :
	localpath="/scratch/"+username+"/ntuples/"+path
	ensure_dir(localpath)
	cleanpath=path;
	if (cleanpath[len(cleanpath)-1]=="/") : # remove trailing slash
		cleanpath=cleanpath[0:len(cleanpath)-2]
	fusecommand="hadd -f /scratch/"+username+"/ntuples/"+cleanpath+".root > /dev/null "
	for item in filelist:
		copycommand="dccp dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/"+username+"/"+item+" /scratch/"+username+"/ntuples/"+item
		commands.getstatusoutput(copycommand)
		fusecommand=fusecommand+" /scratch/"+username+"/ntuples/"+item
	print commands.getoutput(fusecommand)
	deletecommand="rm -r /scratch/"+username+"/ntuples/"+path+"/" 
	print commands.getoutput(deletecommand)
		

def check_directory(path,username) :
	complete_path="/pnfs/psi.ch/cms/trivcat/store/user/"+username+"/"+path
	print "\033[1;34m Going to checkout the subdirectory "+complete_path+" \033[0m "
	listoffiles=[]
	supposedtobejoined=False;
	commandline="lcg-ls -l srm://t3se01.psi.ch:8443/srm/managerv2?SFN="+complete_path
	pipe=popen(commandline)
	for l in pipe.readlines():
		currentline=l.strip("\n")
		if(currentline[0]=="d") :
			check_directory(currentline[currentline.find(path):],username)
		else :
			if(currentline.count(path) > 0) :
				supposedtobejoined=True
				listoffiles.append(currentline[currentline.find(path):])
	if supposedtobejoined==True:
		join_directory(path,listoffiles,username)


##################################################################################
#                              POST PROCESSING                                   #
##################################################################################

def sanitycheck(name,path):
	if os.path.exists(path) == False:
		showMessage("Problem detected with "+name+" which is set to "+path+" -- please fix this. Aborting.")
		sys.exit(-1)
	else :
		return " ... ok !"

##################################################################################
#                                  MAIN                                          #
##################################################################################
if domultiprocessing==True:
	po = Pool()

fusepath=""
uname=""
isdata=0
if(len(sys.argv) != 2 or sys.argv[1] == "-h"):
  printUsage()
else:
  result = parseInputFile(sys.argv[1])
  if(result == "Error"):
    showMessage("Error parsing input file")
  else:
    listOfTasks = getListOfTasks(result[0])
    fusepath=result[4]
    uname=result[6]
    for l in listOfTasks:
      if(l[0].find("/data/")>-1) :
	isdata=1
      showMessage("Processing " + l[0])
      if domultiprocessing==True:
	      po.apply_async(process,(l, result),callback=cb)
      else :
	      print "At this point you could be saving a lot of time with multiprocessing ... "
	      process(l, result)
      		
  if domultiprocessing==True :
	  po.close()
	  po.join()
  totaljobnumber=len(jobnumbers)
  counter=0
  while(len(jobnumbers)>0 and counter<300) :
	  counter+=1
	  currlist=[]
	  pipe=popen("qstat | grep `whoami` | awk '{print $1}'")
	  for line in pipe.readlines():
		currlist.append(int(line))
	  checklist(jobnumbers,currlist)
	  time.sleep(60)
  print "All jobs have finished - need to merge everything and place it in your scratch directory!"
  check_directory(fusepath,uname)
  if isdata==1 and result[2].find("RunJZBAnalyzer")>-1:
    print "We're dealing with data - we still need to merge data files and check for duplicates!"
    pipe=popen("/shome/buchmann/material/flash_remove_duplicates.exec -d /scratch/"+uname+"/ntuples/"+str(fusepath)+"/ -o /scratch/"+uname+"/"+fusepath+".root")
    for l in pipe.readlines():
	  print l
