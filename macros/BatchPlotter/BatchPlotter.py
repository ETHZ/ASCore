#! /usr/bin/python
import os, sys, commands, subprocess, time

def usage():
	print ' Usage: ./BatchPlotter -c config.cfg'

def check_commands():
	print '[status] checking if `root\' and `hadd\' commands exist...'
	if not commands.getstatusoutput('which root')[0]:
		print '[status] `root\' exists...'
		if not commands.getstatusoutput('which hadd')[0]:
			print '[status] `hadd\' exists...'
			return True
	else:
		return False
			
def read_config(config_name):
	if not os.path.isfile(config_name):
		print '[ERROR] config file doesn\'s exist. exiting'
		sys.exit()
	cfg = open(config_name, 'r')
	info = {}
	for line in cfg.readlines():
		if '#' in line or len(line) == 0:
			continue
		opt = line.split()[0]
		arg = line.split()[1]
		if opt == 'PlotterLocation': info['plotter_location'] = arg
		elif opt == 'DumperConfig' : info['dumper_config']    = arg
		elif opt == 'RegionFile'   : info['region_file']      = arg
		elif opt == 'BatchScript'  : info['batch_script']     = arg
		elif opt == 'ReleaseDir'   : info['release_dir']      = arg
		elif opt == 'DataCard'     : info['datacard_location']= arg
		elif opt == 'SELocation'   : info['storage_location'] = arg
		elif opt == 'ModelName'    : info['model_name'] = arg
		elif opt == 'Date':
			info['output_node'] = '/scratch/mdunser/'+arg+'/'
			info['date']        = arg
		elif opt == 'OutputLocation':
			if not arg[-1]== '/':
				arg+='/'
			info['output_location'] = arg
		## elif opt == 'User':
		## 	info['user'] = arg
	cfg.close()
	return info

def nRunningJobs(jobnames):
	tmp = commands.getoutput('qstat')
	jobs = tmp.split('\n')[2:]
	running_jobs = []
	for job in jobs:
		running_jobs.append(job.split()[2])
	cnt=0
	for job in jobnames:
		if job[:10] in running_jobs:
			cnt+=1
	return cnt

def check_on_jobs(jobnames, time_elapsed):
	tmp = commands.getoutput('qstat')
	jobs = tmp.split('\n')[2:]
	running_jobs = []
	for job in jobs:
		running_jobs.append(job.split()[2])
	cnt=0
	for job in jobnames:
		if job[:10] in running_jobs:
			cnt+=1
	if not time_elapsed%60:
		n_min = time_elapsed/60
		## print '[status] '+str(int(commands.getoutput('qstat | wc -l'))-2),' jobs still running after', n_min, 'minutes...'
		print '[status]', cnt, ' jobs still running after', n_min, 'minutes...'
	return cnt

def nRunningJobs():
	tmp = commands.getoutput('qstat')
	jobs = tmp.split('\n')[2:]
	running_jobs = []
	for job in jobs:
		running_jobs.append(job.split()[2])
	return len(running_jobs)
	
def slowSubmit():
	## print jobstrings
	## sys.exit(0)
	njobs = len(jobstrings)
	nsub = 0
	for job in jobstrings:
		## if nsub < 100:
		## 	print 'submitting', job
		## 	os.system(job)
		## 	nsub += 1
		## 	continue
		## else:                            # if more than 100 jobs have already been submitted
		mustend = time.time() + 1000 # wait max 1000 seconds for each new job
		while time.time() < mustend: # max timout
			if nRunningJobs() < 100: # if less than a hundred jobs are running
				## print 'submitting', job
				os.system(job)       # go ahead and submit one
				nsub += 1            # increment number of submitted jobs
				break                # and move to the next one. continue acts on the while, just fyi
			time.sleep(0.1)          # if more than a hundred jobs are running, wait half a second
			

def getRegions(region_file):
	f=open(region_file, 'r')
	lines = f.readlines()
	regs = []
	for line in lines:
		if line[0] == '#': continue
		if line[0] == 'v': continue
		if len(line) == 0: continue
		if line == '\n': continue
		regs.append(line.split()[0])
	f.close()
	return regs

def getFileList(directory):
	# concatenating the strings for input directories
	files = []
	print '[status] searching for all files to run over'
	raw_files = commands.getoutput('srmls '+srm_path+directory)
	## print raw_files
	raw_file_list = raw_files.split('\n')
	for file in raw_file_list:
		if len(file) == 0 or file.split()[0] == '0': continue
		if not '.root' in file: continue
		files.append(dcap_path+file.split()[1])
	return files
	
def merge_and_clean():
	print '[status] now merging and cleaning up...'
	for region in regions:
		hadd_region = ''
		print '[status] at region:', region
		for ls in os.listdir(output_location):
			if not 'output' in ls: continue
			## if os.path.isdir(output_location+'/'+ls):
			for file in os.listdir(output_location+'/'+ls+'/'): 
				if region in file:
					model = file.split('_')[0]
					if not hadd_region.startswith(model):
						hadd_region += model+'_'+region+'.root '
					absfile = output_location+'/'+ls+'/'+file
					hadd_region += ' '+os.path.abspath(absfile)
		hstring = 'hadd '+output_location+'/'+hadd_region
		os.system(hstring+' >& /dev/null')
	os.system('rm -rf '+output_location+'/output_*')
					
	## os.system('rm -r tmp/ ; rm sgejob-* -rf')
	os.system('rm -r tmp/ ; rm plot_* ; rm sgejob-* -rf')
	#os.system('rm sgejob-* -rf') # no cleaning up, for debugging puposes
	## for ls in os.listdir(output_location):
	## 	if os.path.isdir(output_location+ls) and 'output' in ls:

def do_stuff(config_name):
	print '[status] starting script...'

	global srm_path, dcap_path, plotter_location, dumper_config, output_location, user, noj, release_dir, date, output_node, storage_location
	global regions, model_name
	
	print '[status] reading config...'
	info_dict = read_config(config_name)
	
	plotter_location  = info_dict['plotter_location']
	datacard_location = info_dict['datacard_location']
	dumper_config     = info_dict['dumper_config']
	output_location   = info_dict['output_location']
	output_node       = info_dict['output_node']
	date              = info_dict['date']
	batch_script      = info_dict['batch_script']
	storage_location  = info_dict['storage_location']
	release_dir       = info_dict['release_dir']
	regions           = getRegions(info_dict['region_file'])
	model_name        = info_dict['model_name']
	output_location  += date+'_'+model_name+'/'

	srm_path   = 'srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/'
	dcap_path = 'dcap://t3se01.psi.ch:22125/'

	commit_strings = []
	
	files = getFileList(storage_location)
	# checking if output location exists, otherwise create it
	if not os.path.isdir(output_location):
		print '[status] creating output directory:', output_location
		os.system('mkdir '+output_location)
		## copy the dumper_config into the output directory for future reference and the plotter
		## os.system('cp '+dumper_config+' '+output_location+'/dumperconfig.cfg')
	else:
		print '[WARNING] output directory', output_location, 'already exists, this might lead to problems.'

	# checking if a tmp/ subdir exists in the working dir, if so remove it and create an empty one
	if os.path.isdir('tmp'):
		os.system('rm -r tmp')
	os.mkdir('tmp')
	
	print '[status] reading info and creating cards...'
	print '[status] searching for all the files on the SE'
	
	for f in files:
		astring = f.split('/')[-1].rstrip('.root')
		os.system('mkdir '+output_location+'/'+astring)
		for region in regions:
			outloc = output_node + astring
			commit_strings.append(plotter_location + ' -v 1 -c '+datacard_location+' -p '+dumper_config+' -s ' + region + ' -i '+ f + ' -m '+model_name+' -d '+outloc)
			## print plotter_location + ' -v 1 -c '+datacard_location+' -p '+dumper_config+' -s ' + region + ' -i '+ f + ' -d '+outloc
		
	print '[status] starting to submit jobs...'
	global jobnames, jobstrings
	jobnames = []
	jobstrings = []
	

	for commit in commit_strings:
		ind = commit_strings.index(commit)
		copyline = 'mv ' + commit.split()[-1]+'/* ' +output_location+'/'+commit.split()[-1].split('/')[-1]

		shellScript = open(batch_script, 'r')
		tmpScript_name = 'tmp/tmp_script_'+str(ind)+'.sh'
		tmpScript = open(tmpScript_name, 'w')
		for line in shellScript.readlines():
			if 'USERNAMELINE' in line:
				tmpScript.writelines('HN_NAME='+commands.getoutput('echo $USER')+'\n')
			elif 'RELEASEDIRLINE' in line:
				tmpScript.writelines('cd '+release_dir+'\n')
			elif 'EXECLINE' in line:
				tmpScript.writelines(commit+'\n')
			elif 'COPYLINE' in line:
				tmpScript.writelines(copyline+'\n')
			else:
				tmpScript.writelines(line)
		tmpScript.close()
		#os.system('qsub -q all.q  -N job_'+str(ind)+' '+tmpScript_name)
		## print commit ## this would print the ./RunSSDLDumper-command. good for debugging
		if not dryrun:
			#os.system('qsub -q all.q  -N ssdl_'+str(ind)+' '+tmpScript_name)
			##os.system('qsub -q short.q  -N ssdl_'+str(ind)+' '+tmpScript_name)
			jn = 'plot_'+str(ind)
			jobnames.append(jn)
			jobstrings.append('qsub -q short.q  -N '+jn+' '+tmpScript_name+' > /dev/null')
			## os.system('qsub -q short.q  -N '+jn+' '+tmpScript_name+' > /dev/null')
			## os.system('sleep 1')
	if not dryrun:
		print '[status] i\'m going to slowly submit', len(commit_strings), 'jobs ...'
		slowSubmit()
	print '[status] submitted', len(commit_strings), 'jobs'

	print '[status] done submitting jobs...'

	time_elapsed = 0
	while (check_on_jobs(jobnames, time_elapsed) > 0):
		os.system('sleep 1')
		time_elapsed+=1
	print '[status] done with running on the files, it took', time_elapsed, 'seconds!'
	
	if not dryrun:
		if check_commands():
			pass
			merge_and_clean()
		else:
			clean()

def main(args):
	global dryrun
	dryrun = False
	if len(args)==1:
		print 'You must specify at least something... Let me help you:\n'
		usage()
		sys.exit()
	if not check_commands():
		print '[ERROR] you don\'t have \'root\' (and thus \'hadd\') available. This would only lead to trouble.'
		print '        Obtain the commands first and then re-run the same command\n'
		print '[tip]   You can get root e.g. by doing: \'source /swshare/ROOT/thisroot.sh\'\n'
		print 'EXITING...'
		sys.exit()
	if ('-n' in args):
		print 'Dry running only....'
		dryrun = True
	if ('--help' in args) or ('-h' in args):
		usage()
		sys.exit()
	if ('-c' in args) or ('--config' in args):
		do_stuff(str(args[args.index('-c')+1]))

print '[status] starting...'
main(sys.argv)
print '[status] ...done'
