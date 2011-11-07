#!/usr/bin/env python
#______________________________________________________________
# makeTreeClassBase.py
# 
# This should be called within the macros/ directory, giving
# a rootfile containing the desired version of the Ntuples as
# an argument
#______________________________________________________________


import sys, subprocess, os
from subprocess import call
from ROOT import TTree, TFile, gDirectory

usage = "Usage: makeTreeClassBase.py filename.root"

if len(sys.argv) < 2:
    print usage
    exit(1)

FILENAME = sys.argv[1]
CLASSNAME = 'TreeClassBase'
HEADERNAME = CLASSNAME + '.h'
SOURCENAME = CLASSNAME + '.C'
MAXPFLEPT     = 20
MAXNJETS      = 100
MAXNPHOS      = 50
MAXNGENPHOS   = 100
MAXNSC        = 100
rules = { # Reco
          'NMus'    : 30,
          'NEles'   : 20,
          'NPhotons': 50,
          'NJets'   : MAXNJETS,
          # PF
          'PfMuNObjs'  : MAXPFLEPT,
          'PfMu2NObjs' : MAXPFLEPT,
          'PfMu3NObjs' : MAXPFLEPT,
          'PfElNObjs'  : MAXPFLEPT,
          'PfEl2NObjs' : MAXPFLEPT,
          'PfEl3NObjs' : MAXPFLEPT,
          'PfTauNObjs' : MAXPFLEPT,
          'PfTau2NObjs': MAXPFLEPT,
          'PfTau3NObjs': MAXPFLEPT,
          'PF2PATNJets' : MAXNJETS,
          'PF2PAT2NJets': MAXNJETS,
          'PF2PAT3NJets': MAXNJETS,
          # Anti-iso
          'PfMuAntiIsoNObjs'  : MAXPFLEPT,
          'PfElAntiIsoNObjs'  : MAXPFLEPT,
          'PfTauAntiIsoNObjs' : MAXPFLEPT,
          'PF2PATAntiIsoNJets': MAXNJETS,
          # Other jets
          'CANJets' : MAXNJETS,
          'JPTNJets': MAXNJETS,
          # Generator
          'NGenLeptons': 100,
          'NGenJets'   : 100,
          'NGenPhotons': 100,
	  # Photons
	  'MaxNPhotons':MAXNGENPHOS,
	  'NSuperClusters':MAXNSC,
          # Others
          'NTracks' : 500,
          'NPaths'  : 10,
          'NVrtx'   : 25,
          'NEBhits' : 20,
          'PUnumInteractions': 50,
          'PUnumFilled'      : 50
          }

#______________________________________________________________
def makeClass(filename, classname, treename):
	f = TFile.Open(filename)
	tree = gDirectory.Get(treename)
	tree.MakeClass(classname)
	f.Close()

#______________________________________________________________
if __name__=='__main__':
	print 'makeTreeClassBase >> Creating MakeClass from ' + FILENAME
	makeClass(FILENAME, CLASSNAME, 'analyze/Analysis')

	headerFile = open(HEADERNAME, 'r')
	headerLines = headerFile.readlines()
	headerFile = open(HEADERNAME, 'w')
	
	for line in headerLines:
                for key,val in rules.iteritems():
                    pos = line.find('//['+key+']')
                    if ( pos != -1):
                       line = line[0:line.find('[')+1] + str(val) + line[pos-5:len(line)]
		headerFile.write(line)
	
	headerFile.close()

	sourceFile = open(SOURCENAME, 'r')
	sourceLines = sourceFile.readlines()
	sourceFile = open(SOURCENAME, 'w')
	
	sourceLines[1] = '#include "base/' + HEADERNAME + '"\n'
	for line in sourceLines:
		sourceFile.write(line)
		
	sourceFile.close()
	
	call(['mv', '-f', HEADERNAME, 'include/base/'])
	call(['mv', '-f', SOURCENAME, 'src/base/'])
