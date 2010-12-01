#!/usr/bin/python

#--------------------------------------#
# complaints to pascal.nef@cern.ch     #
#--------------------------------------#

import sys, os, commands

# check the command line arguments
if len(sys.argv[1:]) != 1:        # if no dir-name is given
	print 'Usage:' 
	print ' python GetNEvents_CrabReport.py dir-name'
	print ' example: ./GetNEvents_CrabReport.py /ntuples/data/V01-11-02/MultiJet-Run2010B-Nov4ReReco_v1_RECO '
	sys.exit()
else:
	jobDir = sys.argv[1]
		
T3SRMDIR="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/susy/"

cmd1   = "srmls srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/susy/"+jobDir + \
	 " | grep .root | awk '{print $2}' | awk '{ print \"dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/\" substr( $0, 25, length($0) ) }'"
stats, files = commands.getstatusoutput(cmd1);
f = open('listoffiles.txt', 'w')
f.writelines(files)
f.close()

# run analyzer on it to get the number of events in the chain of root files
runanalyzer = "/shome/pnef/SUSY/SUSY_macros/macros/RunJZBAnalyzer -o tmp.root -v 1 -n 1 -t mc -c -l listoffiles.txt | grep \"total events in ntuples\" | awk '{print $6}'"
stats, nevents = commands.getstatusoutput(runanalyzer)
print  nevents

os.system("rm tmp.root")
os.system("rm listoffiles.txt")

