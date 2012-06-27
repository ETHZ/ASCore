#!/usr/bin/env python

import optparse
import sys, re
import os, subprocess


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
    It must be "yes" (the default), "no" or None (meaning
    an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                "(or 'y' or 'n').\n")

def refreshFileList( srmPath, tmpFile, options ):
#     print 'Reloading list of files in '+tmpFile+'...',
#     sys.stdout.flush()

    cmd = 'srmls -count=1000 '+srmPath+' | grep .root > '+tmpFile
    if options.debug: print '\nRunning',cmd
    srmls = subprocess.Popen(cmd, shell=True)
    srmls.wait()
    attempts = 1
    while ( (sum(1 for line in open(tmpFile)))==1000*attempts ):
        off = 1000*attempts
        srmls = subprocess.Popen('srmls -count=1000 -offset=%s'%off+' '+srmPath+' | grep .root >> '+tmpFile, shell=True)
        srmls.wait()
        attempts+=1
#     print 'Done'


def performDelete( files, srmSite, logFile ):


    toDelete = sorted(files, key = lambda file: file[2])
    print 'The following files will be deleted:'
    print 'Id \t Retry \t Size \t Path #'
    for file in toDelete:
        print file[2],'\t',file[3],'\t',file[1],'\t',file[0]
    
    if query_yes_no('Do you want to proceed?','no'):
        log = open(logFile,'w')
        for file in toDelete:
            lcgdel = subprocess.Popen(['lcg-del','-l','-v',srmSite+file[0]], stdout=log,stderr=log)
            lcgdel.wait()
        log.close()
        print 'Log saved to',logFile
    else:
        print 'No file deleted'

def findDuplicates():

    usage = 'usage: %prog [options] path'
    parser = optparse.OptionParser(usage)
    parser.add_option('--site', dest='site', help='Site where files are located. Can be [t2cscs,t3psi]')
    parser.add_option('--moveTo', dest='moveTo', help='Destination for the duplicates')
    parser.add_option('--refresh', dest='refresh', help='Refresh the list of files', action='store_true')
    parser.add_option('--tryDelete', dest='tryDelete', help='Attempts to delete the duplicates: Keeps the bigger file for each job', action='store_true')
    parser.add_option('--deleteAll', dest='deleteAll', help='Delete all the duplicates', action='store_true')
    parser.add_option('--dcap', dest='dcap', help='Print the list of files in dcap format', action='store_true')
    parser.add_option('--auto', dest='auto', help='Automatically guess remote location (argument should be crab task directory)', action='store_true', default=False)
    parser.add_option('--debug', dest='debug', help='Turn on debugging information', action='store_true', default=False)
    parser.add_option('--dbs', dest='DBS', help='Check also DBS information', action='store_true', default=False)
#     parser.add_option('--tag', dest='tag', default='MC', help='tag to match the files [MC,data]')


    (opt, args) = parser.parse_args()

    if not opt.site:
        parser.error('No site selected')
    if opt.site == 't3psi' or opt.site == 'T3_CH_PSI':
        srmSite = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN="
        pfnPrefix = "/pnfs/psi.ch/cms/trivcat"
        rootpath = srmSite+pfnPrefix
        dcapPrefix= 'dcap://t3se01.psi.ch:22125'
        pfnPrefix = "/pnfs/psi.ch/cms/trivcat"
    elif opt.site == 't2cscs' or opt.site == 'T2_CH_CSCS':
        srmSite = "srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN="
        pfnPrefix = "/pnfs/lcg.cscs.ch/cms/trivcat"
        rootpath = srmSite+pfnPrefix
        dcapPrefix= ''
        pfnPrefix = "/pnfs/lcg.cscs.ch/cms/trivcat"
    elif opt.site == 't2rwth' or opt.site == 'T2_DE_RWTH':
        srmSite = "srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="
        pfnPrefix = "/pnfs/physik.rwth-aachen.de/cms"
        rootpath = srmSite+pfnPrefix
        dcapPrefix= ''
    else:
        parser.error('site can be t3psi, t2cscs, or t2rwth')

    if len(args)!=1:
        parser.error('Wrong number of arguments')
    
    path = args[0]
    if path[0] is not '/' and not opt.auto:
        parser.error('Requires an absolute path: It must start with \'/\'')

    PWD=os.getenv("PWD")
    dbsFileList=os.path.basename(path)+'.'+opt.site+'.dbs.list'
    SiteTranslate = {"t2cscs":"T2_CH_CSCS", "t3psi":"T3_CH_PSI","t2rwth":"T2_DE_RWTH"}
    if opt.DBS:
        ###TODO: check for DBSSQL
        datasetName = os.popen("grep '<User Dataset Name>' "+path+"/log/crab.log | cut -d= -f2 | uniq").readline().strip(" ").strip("\n")
        print "Considering dataset: "+datasetName+", getting files"
        os.system("dbssql --limit=10000 --dbsInst=cms_dbs_ph_analysis_02 --input='find file,site where dataset="+datasetName+" and site="+SiteTranslate[opt.site]+"' | grep -v '#' | awk '{print $1}'>"+dbsFileList)
        

#     fileTag=opt.tag

    if path[-1] == '/':
        path = path[:-1]
    tmpFile = os.path.basename(path)+'.'+opt.site+'.lst'
    logFile = os.path.basename(path)+'.'+opt.site+'.log'
    dcapFile = os.path.basename(path)+'.'+opt.site+'.dcap'

    srmPath = rootpath+path
    # Automatic recognition of srm path: in that case, path is a crab task directory
    if opt.auto:
        query = 'select lfn from bl_runningjob where lfn not null limit 1'
        if opt.debug: print query
        cmd = subprocess.Popen(['sqlite3',path+'/share/crabDB',query],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        cmd.wait()
        if cmd.returncode != 0:
            print "Unable to retrieve storage path automatically:",cmd.stdout.readlines()
            sys.exit(-1)
        else:
            # cmd.stdout.readlines() returns an array of lines: we take the first
            result = eval(cmd.stdout.readlines()[0])
            # the line is actually a python array of files: we take the first
            srmPath = rootpath + os.path.dirname(result[0])
    
    
    hline = '-'*80
    print hline
    print 'Processing:',os.path.basename(path)
    print 'Path:',srmPath
    print hline
    
    if not os.path.exists(tmpFile) or opt.refresh:
        print 'Fetching the list of files to '+tmpFile+'...',
        sys.stdout.flush()

        refreshFileList(srmPath,tmpFile,opt)
        print 'Done'
    else:
        print 'Using',tmpFile,'as source'

    out = open(tmpFile).read()
    if len(out) is 0:
        print 'The file',tmpFile,'exists but is probably empty. Trying to fetch the list once more...',
        sys.stdout.flush()
        refreshFileList(srmPath,tmpFile,opt)
        print 'done'
        out = open(tmpFile).read()
 
#     print 'Processing file list using file tag',fileTag
    lines = out.splitlines()
    # remove the directory name 
    #lines.pop(0)    # not needed with new list method

    #initialize some variables
    sizeOnDisk = 0.
    nFiles = 0
    nDuplicates= 0
    missingFiles=[]

    fileMap = {}
    duplicatesMap = {}
    lastId=0
    #fill the file register
    for line in lines:
        if len(line) is 0:
            continue
        tokens=line.split()
        res = re.search('_([0-9]{1,})_([0-9]{1,})_([0-9a-zA-Z]{3}\.root)',tokens[1])
        if res is None:
            continue

        
        path = tokens[1]
        id = int(res.group(1))
        retry = int(res.group(2))
        size = int(tokens[0])
        lfn= path[len(pfnPrefix):].replace("//","/")

        # make a tuple, path, size, id, retry
        fileTuple=(path,size,id,retry,lfn)
        # count the total size
        sizeOnDisk += size
        # fill the map, where the key is the job id
        if not fileMap.has_key(id):
            fileMap[id] = []
        fileMap[id].append(fileTuple)

    if opt.dcap:
        dcap = open(dcapFile,'w')
        for id,fileArray in fileMap.items():
            for fileTuple in fileArray:
                dcap.write(dcapPrefix+fileTuple[0]+'\n')
        dcap.close()
        print 'dcap file list dumped to',dcapFile

    # sort the keys(job ids) to get the highest
    lastId = sorted(fileMap.keys())[-1]
    print 'Highest Id found:',lastId
    missingFiles = set(range(1,lastId+1))-set(fileMap.keys())
    missingFiles = sorted(missingFiles)

    # re-loop to identify the duplicates
    for id,fileArray in fileMap.items():
        length = len(fileArray)
        nFiles += length
        if length >1:
            nDuplicates += length
            duplicatesMap[id]=fileArray

    sizeGb = sizeOnDisk / 1024**3
    print 'Files found:',nFiles,'Duplicates:',nDuplicates
    print 'Total file size: %.2fGb' % sizeGb

    dbsDuplicates = {}
    if len(duplicatesMap) is not 0:
        print 'Id \t Retry \t Size \t Path #'
        for id in sorted(duplicatesMap.iterkeys()):
            
            fileArray = duplicatesMap[id]
            print '---',id,'-',len(fileArray),'duplicates'
            nInDbs=0

            for fileTuple in fileArray:
                res = os.popen("grep "+fileTuple[4].split("/")[-1]+" "+dbsFileList).readlines()
                dbsString = ""
                if res!=[]:
                    dbsString="DBS"
                    nInDbs+=1
                print fileTuple[2],'\t',fileTuple[3],'\t',fileTuple[1],'\t',fileTuple[0],'\t',dbsString
                dbsDuplicates[str(fileTuple[2])] = nInDbs

    if len(missingFiles) == 0:
        missStr = 'None'
    else:
        missStr='['
        for i in missingFiles:
            missStr += str(i)+','
        missStr = missStr[:-1]+']'

    print 'Missing files','('+str(len(missingFiles))+'):',missStr
    print hline

    nPub=0
    if len(dbsDuplicates)!=0:
        dRange = ""
        for d in dbsDuplicates.keys():
            if dbsDuplicates[d]>1: dRange+=d +" ("+str(dbsDuplicates[d])+"),"
            if dbsDuplicates[d]>0: nPub+=1
        print "Duplicate files in DBS (legend: [ID (nCopies)]):  "+"["+dRange[:-1]+"]"
    print "Published files (unique): ", nPub

    if opt.tryDelete:
        filesToDelete = []
        nUnsafeDuplicates = 0
        for files in duplicatesMap.itervalues():
            canDelete = True
#             # take the size of the first file as a reference
#             refSize = files[0][1]
#             for file in files:
#                 # check if all the sizes are the same
#                 canDelete &= file[1] == refSize

#             if not canDelete:
#                 print '--- ',file[2],'The following files have mismatching size: no action taken'
#                 for file in files:
#                     print file[1],'\t',file[0],'\t',file[3]
#                     nUnsafeDuplicates += 1
            
            #if we are happy with the filesize
            dups = sorted(files, key = lambda file: file[1])
            # don't the one with the highest retry
            dups.pop(-1)
            filesToDelete.extend(dups)

        print 'N of duplicates safe to delete:',len(filesToDelete),' N of duplicates with mismatching size:',nUnsafeDuplicates
        if len(filesToDelete)>0: performDelete( filesToDelete, srmSite, logFile )

    elif opt.deleteAll:
        # make one list of duplicates
        filesToDelete = []
        for files in duplicatesMap.itervalues():
            filesToDelete.extend(files)

        print 'Duplicate files to delete: ',len(filesToDelete)
        if len(filesToDelete)>0: performDelete( filesToDelete, srmSite, logFile )

    return


    print 'Lost files:',lostStr
    if opt.moveTo:
        newPath = rootpath+opt.moveTo
        if newPath[-1] != '/':
            newPath +='/'
        newPath += os.path.basename(path)+'/'
        print 'Destination:',newPath
        log = open(logFile,'w')
        for file in duplicates:
            oldFile=srmSite+file
            newFile=newPath+os.path.basename(file) 
            lcgcp = subprocess.Popen(['lcg-cp','-v',oldFile,newFile],stdout=log,stderr=log)
            lcgcp.wait()
            lcgdel = subprocess.Popen(['lcg-del','-l','-v',oldFile],stdout=log,stderr=log)
            lcgdel.wait()

        log.close()
        print 'Transfer logfile saved to',logFile
    elif opt.delete:
        log = open(logFile,'w')
        for file in duplicates:
            lcgdel = subprocess.Popen(['lcg-del','-l','-v',srmSite+file], stdout=log,stderr=log)
            lcgdel.wait()

        log.close()

if __name__ == '__main__':
        findDuplicates()
