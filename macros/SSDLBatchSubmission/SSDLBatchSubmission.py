#! /usr/bin/python
import os, sys, commands, subprocess

def usage():
	print ' Usage: ./SSDLBatchSubmission -c BatchConfig.cfg'
	print '--------------------------------------------------------------------'
	print 'Here are your options:'
	print '        -c <your_config_file>'
	print '              Does all the work. Takes the input of the datacard and creates a lot of'
	print '              jobs that go to the batch system. Once the jobs are finished it merges'
	print '              them together into a file called allYields.root in your output directory'
	print '              The directories and some other things must be specified in the config'
	print '              file. IMPORTANT: Read below for more details'
	print '              After finishing everything, the script also cleans up all the mess it made with'
	print '              temporary folders.'
	print '        -h'
	print '              Shows this message.\n'
	print '        -n'
	print '              Starts a dry run, i.e. it doesn\'t submit anything to the batch system.\n'
	print '              This will not clean up anything, the tmp/ directory with all the\n'
	print '              non-submitted scripts will still be there. Good for debugging.\n'
	
	print '--------------------------------------------------------------------'
	print 'The config file must look as follows:'
	print '        DumperLocation     <path to the RunSSDLDumper executable>'
	print '        OutputLocation     <the desired output location. make sure the folder exists and that there is NO allYields.root in it!>'
	print '        OutdirNode         <temporary folder on the workernodes /scratch/ directory'
	print '                           this dir will be created and stuff will be moved out of it, so it doesn\'t really matter what it\'s called>'
	print '        DataCard           <datacard as usual, but you must include the entire path after stiegerb/ (i.e. don\'t forget the 2011(B)-part'
	print '                           omit the \'.root\' in the filename if you want to split up the job, the resulting directory MUST exist on the SE'
	print '                           it is crucial that this datacard has only existing paths!>'
	print '        BatchScript        <There are a few things you must consider:'
	print '                                   replace \'HN_NAME=your-hn-name\' with \'USERNAMELINE\''
	print '                                   The \'YOUR FUNCTIONALITY CODE GOES HERE\' section MUST look as follows:'
	print '                                   export SCRAM_ARCH=slc5_amd64_gcc434'
	print '                                   source $VO_CMS_SW_DIR/cmsset_default.sh'
	print '                                   RELEASEDIRLINE'
	print '                                   eval `scramv1 runtime -sh`'
	print '                                   export LD_LIBRARY_PATH=/swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH'
	print '                                   EXECLINE'
	print '                                   >'
	print '        ReleaseDir         <directory with a working cmssw release>'
	print '        User               <well, your username>\n'
	print '--------------------------------------------------------------------'
	print 'Make sure you have a valid voms-proxy thing running, otherwise the srmls command'
	print '(which is needed in the script) does not work'

def mk_dir_string(line):
	# concatenating the strings for input directories
	files = []
	path = line.split()[1]
	if not path[-1] == '/':
		path += '/'
	print '[status] searching for all files for split up jobs, at', line.split()[0]
	raw_files = commands.getoutput('srmls '+srm_path+path)
	raw_file_list = raw_files.split('\n')
	for file in raw_file_list:
		if len(file) == 0 or file.split()[0] == '0':
			continue
		files.append(dcap_path+file.split()[1])
	strings = []
	for file in files:
		indiv_dir = line.split()[0]+'_'+file.split('/')[-1].replace('.root','')
		strings.append(mk_single_string(line, file, indiv_dir))
	print '[status] found', len(strings), 'files for directory', line.split()[0]
	return strings
		#ind = files.index(file)
		#if not ind%nof:
		#	if not ind == 0:
		#		new_cardFile.close()
		#	new_cardFile = open('tmp/tmp_datacard_'+line.split()[0]+str(ind/nof)+'.dat', 'w')
		#new_cardFile.writelines(line.split()[0]+'\t'+file+'\t'+line.split()[2]+'\t'+line.split()[3]+'\n')
	

def mk_single_string(line, full_location=False, indiv_dir=False):
	string = dumper_location + ' -v 1'
	if not full_location:
		string+= ' -i '+ dcap_path+'/pnfs/psi.ch/cms/trivcat/store/user/'+ line.split()[1]
		#string+= ' -i '+ dcap_path+'/pnfs/psi.ch/cms/trivcat/store/user/stiegerb/'+ line.split()[1]
	else:
		string+=' -i '+ full_location
	string+= ' -n '+line.split()[0]
	string+= ' -m '+line.split()[2]
	string+= ' -c '+line.split()[3]
	if len(line.split()) > 4:
		string+= ' -g '+line.split()[5]
		string+= ' -x '+line.split()[4]
	if not indiv_dir:
		string+= ' -o '+output_node
	else:
		string+= ' -o '+output_node+'/'+indiv_dir
	return string

def mk_card_line(line):
	name = line.split()[0]+'\t'
	loc  = line.split()[1]+'\t'
	xsec = line.split()[4]+'\t'
	ngen = line.split()[1]+'\t'
	card_line = name+dcap_path+'/pnfs/psi.ch/cms/trivcat/store/user/'+loc+line.split()[2]+'\t'+line.split()[3]+'\t'+xsec+'\t'+ngen+'\n'
	#card_line = name+dcap_path+'/pnfs/psi.ch/cms/trivcat/store/user/stiegerb/'+loc+line.split()[2]+'\t'+line.split()[3]+'\n'
	return card_line

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
		if opt == 'DumperLocation':
			info['dumper_location'] = arg
		elif opt == 'OutdirNode':
			if not arg[-1]== '/':
				arg+='/'
			info['output_node'] = arg
		elif opt == 'OutputLocation':
			if not arg[-1]== '/':
				arg+='/'
			info['output_location'] = arg
		elif opt == 'DataCard':
			info['cardfile_name'] = arg
		elif opt == 'User':
			info['user'] = arg
		elif opt == 'BatchScript':
			info['batch_script'] = arg
		elif opt == 'ReleaseDir':
			info['release_dir'] = arg
	cfg.close()
	return info

def print_status(t):
	n_min = t/60
	print '[status] '+str(int(commands.getoutput('qstat | wc -l'))-2),' jobs still running after', n_min, 'minutes...'

def check_commands():
	print '[status] checking if `root\' and `hadd\' commands exist...'
	if not commands.getstatusoutput('which root')[0]:
		print '[status] `root\' exists...'
		if not commands.getstatusoutput('which hadd')[0]:
			print '[status] `hadd\' exists...'
			return True
	else:
		return False
		#print '[status] `root\' and `hadd\' do not exist!'
		#print '[status] trying to obtain them...'
		#os.system('cd '+release_dir+' ; eval `scramv1 runtime -sh`')
		#os.system('eval `scramv1 runtime -sh`')
		#print '[status] checking again...'
		#if commands.getstatusoutput('which root')[0]==0:
		#	print '[status] `root\' exists now...'
		#	if not commands.getstatusoutput('which hadd')[0]:
		#		print '[status] and so does `hadd\'...'
		#elif commands.getstatusoutput('which root')[0] == 256:
		#		print '[ERROR] could not obtain `root\' and `hadd\''
		#		print '[ERROR] you should have the output files in your output directory, you will have to merge them yourself.'
		#		print '[status] cleaning up the mess...'
		#		return False
			
def clean():
	print '[status] cleaning up a bit...'
	ls = os.listdir(output_location)
	for obj in ls:
		#ind = ls.index(obj)
		if not os.path.isdir(output_location+obj):
			continue
		elif os.path.isdir(output_location+obj):
			os.system('cd '+output_location+' ; mv '+output_location+obj+'/* '+obj+'.root' )
			os.system('rm -r '+output_location+obj+'/' )
	os.system('rm -r tmp/ ; rm job_* ; rm -r sgejob*')
	
		
def merge_and_clean():
	print '[status] now merging and cleaning up...'
	print '[status] starting with the special directories...'
	for dir in special_dirs:
		print '[status] at special dir:', dir
		dir_cat = 'cat '
		isdata=False
		for ls in os.listdir(output_location):
			if os.path.isdir(output_location+ls) and ls.startswith(dir+'_output'):
				if os.path.isfile(output_location+ls+'/'+dir+'_SignalEvents.txt'):
					dir_cat+=output_location+ls+'/'+dir+'_SignalEvents.txt '
					isdata=True
		dir_hadd = 'hadd '+output_location+dir+'_Yields.root '+output_location+dir+'_output*/*.root'
		dir_cat+=' >& '+output_location+dir+'_SignalEvents.txt '
		os.system(dir_hadd)
		if isdata:
			os.system(dir_cat)

	# print '[status] done with the special dirs, now merging all together.'
	# hadd_string = 'hadd '+output_location+'SSDLYields.root'
	# for ls in os.listdir(output_location):
	# 	if ls.endswith('.root'):
	# 		hadd_string+=' '+output_location+ls
	# os.system(hadd_string)
	os.system('rm -r tmp/ ; rm job_* ; rm sgejob-* -rf')
	for ls in os.listdir(output_location):
		if os.path.isdir(output_location+ls) and 'output' in ls:
			os.system('rm -rf '+output_location+ls)
		

def do_stuff(config_name):
	print '[status] starting script...'

	global srm_path, dcap_path, dumper_location, output_location, user, noj, release_dir, output_node
	global special_dirs
	special_dirs = []
	
	print '[status] reading config...'
	info_dict = read_config(config_name)
	
	cardFile        =  open(info_dict['cardfile_name'], 'r')
	dumper_location = info_dict['dumper_location']
	output_location = info_dict['output_location']
	output_node     = info_dict['output_node']
	batch_script    = info_dict['batch_script']
	user            = info_dict['user']
	release_dir     = info_dict['release_dir']

	srm_path   = 'srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/'
	dcap_path = 'dcap://t3se01.psi.ch:22125/'

	noj = 10
	
	commit_strings = []
	
	# checking if output location exists, otherwise create it
	if not os.path.isdir(output_location):
		print '[status] creating output directory:', output_location
		os.system('mkdir '+output_location)
	else:
		print '[WARNING] output directory', output_location, 'already exists, this might lead to problems.'

	# checking if a tmp/ subdir exists in the working dir, if so remove it and create an empty one
	if os.path.isdir('tmp'):
		os.system('rm -r tmp')
	os.mkdir('tmp')
	
	#lm_cardFile  = open('tmp/tmp_datacard_lm.dat'  , 'w')
	#qcd_cardFile = open('tmp/tmp_datacard_qcd.dat' , 'w')
	#ww_cardFile  = open('tmp/tmp_datacard_ww.dat'  , 'w')

	
	print '[status] reading info and creating cards...'
	print '[status] searching for all the files on the SE'
	for line in cardFile.readlines():
		if '#' in line or len(line) == 0:
			continue
		if line.split()[1].endswith('.root'):
		#	if 'LM' in line.split()[0]:
		#		lm_cardFile.writelines(mk_card_line(line))
		#	elif 'QCD' in line.split()[0]:
		#		qcd_cardFile.writelines(mk_card_line(line))
		#	elif 'WW' in line.split()[0]:
		#		ww_cardFile.writelines(mk_card_line(line))
		#	else:
			commit_strings.append(mk_single_string(line))
		else:
			print '[status] looking up single files for', line.split()[0]
			special_dirs.append(line.split()[0])
			commit_strings.extend(mk_dir_string(line))
		
	#lm_cardFile.close()
	#qcd_cardFile.close()
	#ww_cardFile.close()
	
	cardFile.close()
	
	for card in os.listdir('tmp'):
		if not 'datacard' in card:
			continue
		if os.path.getsize(os.path.abspath('.')+'/tmp/'+card) == 0:
			continue
		commit_strings.append(dumper_location + ' -v 1 -l ' + os.path.abspath('.')+'/tmp/'+card+' -o '+output_node)
	
	print '[status] starting to submit jobs...'

	for commit in commit_strings:
		ind = commit_strings.index(commit)
		my_line = commit.split()
		name = 'xkcd'
		if '-n' in my_line:
			name = my_line[my_line.index('-n')+1]
		#print 'At job number:', ind
		# if not ind == 1: continue # to test with only the first job
		# this is not to copy everything back directly to the /shome
		# instead it creates a directory /jobnumber/ on the workernode /scratch/ and afterwards moves everything from there to /shome
		commit+=str(ind)+'/'
		if name in special_dirs:
			copyline = 'mv ' + commit.split()[-1]+' '+output_location
		else:
			copyline = 'mv ' + commit.split()[-1]+'* ' +output_location

		shellScript = open(batch_script, 'r')
		tmpScript_name = 'tmp/tmp_script_'+str(ind)+'.sh'
		tmpScript = open(tmpScript_name, 'w')
		for line in shellScript.readlines():
			if 'USERNAMELINE' in line:
				tmpScript.writelines('HN_NAME='+user+'\n')
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
		if not dryrun:
			os.system('qsub -q short.q  -N job_'+str(ind)+' '+tmpScript_name)
	print '[status] submitted', len(commit_strings), 'jobs'

	print '[status] done submitting jobs...'

	time_elapsed = 0
	while (int(commands.getoutput('qstat | wc -l'))-2)>0:
		os.system('sleep 1')
		time_elapsed+=1
		if not time_elapsed%60:
			print_status(time_elapsed)
	print '[status] done with running on the files, it took', time_elapsed, 'seconds!'
	
	if not dryrun:
		#if check_commands():
		merge_and_clean()
		#else:
		#	clean()
	
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
