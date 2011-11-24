#!/usr/bin/python
import sys, os, pwd, commands, shlex, glob

def updateFileList(filelistName,srmDir):
   """Update list of files to process"""
   command='srmls -count=1000 -offset=0 '+srmDir+' | grep .root | perl -an -F"/" -e \'print "/".join("/",@F[5 .. $#F])\' > ' + filelistName
   return_value = os.system(command)
   if return_value != 0:
      print 'Problem to determine which jobs need to be retrieved. Exiting...'
      sys.exit()   
   attempts = 1
   while ( (sum(1 for line in open(filelistName)))==1000*attempts ):
      off = 1000*attempts
      command='srmls -count=1000 -offset=%s'%off+' '+srmDir+' | grep .root | perl -an -F"/" -e \'print "/".join("/",@F[5 .. $#F])\' >> '+filelistName
      return_value = os.system(command)
      if return_value != 0:
         print 'Problem to determine which jobs need to be retrieved. Exiting...'
         sys.exit()   
      attempts+=1



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

# Now run a loop
all_ok = True
logfilename = jobName+'_replica'
doneFiles = []
while (all_ok):
   updateFileList(filelistName,T2SRMDIR)

   
   # run data_replica
   command='/swshare/psit3/bin/data_replica.py --copy-tool=srmcp --delete --debug '+filelistName+' --from-site '+T2NAME+' --to-site T3_CH_PSI /store/user/susy/'+jobDir
   print "Running",command
   return_value = os.system(command)
   if return_value != 0:
      print 'Problem in replicating data with "data_replica".'
      all_ok = False
      sys.exit()

   # In principle one could loop here...
   all_ok = False
   
