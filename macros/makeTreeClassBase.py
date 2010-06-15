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
MAXNLEPT = 20
MAXNGENLEPT = 100
MAXNJETS = 100
MAXNPHOS = 50
MAXNTRKS = 500

#______________________________________________________________
def makeClass(filename, classname, treename):
	f = TFile(filename)
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
		pos = line.find('//[NGenLeptons]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNGENLEPT) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NMus]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNLEPT) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NEles]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNLEPT) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NPhotons]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNPHOS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NJets]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNJETS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[PFNJets]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNJETS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[CANJets]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNJETS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[JPTNJets]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNJETS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NTracks]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNTRKS) + line[pos-5:len(line)]
			line = newline
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
