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
MAXNMUS       = 30
MAXNELES      = 20
MAXNTAUS      = 20
MAXNGENLEPT   = 100
MAXNJETS      = 100
MAXNGENJETS   = 100
MAXNPHOS      = 50
MAXNTRKS      = 500
MAXHLTOBJ     = 10
MAXNVTCS      = 25
MAXNPILEUP    = 50
MAXNEBHITS    = 20

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
			newline = line[0:line.find('[')+1] + str(MAXNMUS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NPfMus]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNMUS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NEles]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNELES) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NPfEls]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNELES) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NPfTaus]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNTAUS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NPhotons]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNPHOS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NJets]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNJETS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NGenJets]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNGENJETS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[PF2PATNJets]')
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
		pos = line.find('//[NPaths]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXHLTOBJ) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NVrtx]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNVTCS) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[PUnumInteractions]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNPILEUP) + line[pos-5:len(line)]
			line = newline
		pos = line.find('//[NEBhits]')
		if(pos != -1):
			newline = line[0:line.find('[')+1] + str(MAXNEBHITS) + line[pos-5:len(line)]
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
