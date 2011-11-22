#! /usr/bin/python
import commands, os

def print_status(t):
    n_min = t/60
    print '[status] still running after', n_min, 'minutes...'

for m0 in range(20, 2000, 20):
	os.system('qsub -q short.q -N job_limit -e $PWD -o $PWD wrapper.sh ' + str(m0))

time_elapsed = 0
while (int(commands.getoutput('qstat | wc -l'))-2)>0:
	os.system('sleep 1')
	time_elapsed+=1
	if not time_elapsed%60:
		print_status(time_elapsed)
		print '[status] still '+ str(int(commands.getoutput('qstat | wc -l'))-2) +' jobs running'

print '[status] done with running on the files, it took', time_elapsed, 'seconds! you do the math how many minutes this is. hints are just above.'

print '[status] cleaning up...'

os.system('rm job_limit.*')

print '[status] i hope you have sourced a version of ROOT, otherwise this won\'t work...'

if os.path.isfile('limit.root'):
	print '[status] limit.root already exists, i will not merge anything in this case... hrmpf'
	print '[status] find instead all the limits_*.root files in this directory and add them together yourself'
	print '[status] done.'
	sys.exit()

os.system('hadd limit.root limits*.root')
os.system('rm limits*.root')

print '[status] find the exclusion curve the file limit.root'
print '[status] done.'
