#!/usr/bin/env python
import os

# Small script to make a datacard for the Higgs Limit calculation
# instert the wanted name (.dat will be appended) of the output datacard and the numbers from your
# counting experiment and then execute with ./makeCard.py and you should be fine
# for any more complicated issues (i.e. more channels, more systematic uncertainties etc.)
# one would have to expand this script

cardName = 'card_test'

nSignal = 1
nBkg = 9.24
nBkgErr = 3.91
nObs = 8

bkgUncert = nBkgErr/nBkg
for sigUncert in [0.05, 0.1, 0.15, 0.20, 0.25, 0.35, 0.4, 0.45, 0.5]:
	z= open(cardName+'_sigUncert'+str(int(100*sigUncert))+'.dat', 'w')
	z.write('''# Simple counting experiment, with one signal and one background process
	imax 1  number of channels
	jmax 1  number of backgrounds 
	kmax 2  number of nuisance parameters (sources of systematical uncertainties)
	------------
	# we have just one channel, in which we observe 0 events
	bin         1
	observation {a}
	------------
	# now we list the expected events for signal and all backgrounds in that bin
	# the second \'process\' line must have a positive number for backgrounds, and 0
	# for signal
	# then we list the independent sources of uncertainties, and give their effect
	# (syst. error)
	# on each process and bin
	bin            1      1
	process      signal   Bckg 
	process        0      1
	rate           1      {b}
	------------
	deltaS  lnN    {c}    -    {c1}% uncertainty on signal
	deltaB  lnN      -   {d}   {e}% uncertainty on background, i.e. the error on the {f0} is {f}
	'''.format(a=str(nObs), b=str(nBkg), c=str((1+sigUncert)), c1=str(int(100*sigUncert)), d=str(1+bkgUncert), e=str(int(100*bkgUncert)), f0=str(nBkg), f=nBkgErr ) )
	z.close()
