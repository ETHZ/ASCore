#!/usr/bin/python

#
# Reset a runmanager task
#
# Takes a list of task directories as argument

import sys,re
from os import path
import os,subprocess
from optparse import OptionParser

def getSeDir(dir,opt):
   # open the script and get the remote directory
   
   sehome = ''
   sesubdir = ''
   hnname = ''
   sePat  = re.compile('^USER_SRM_HOME\s*=\s*"?(srm:[^"]*)')
   subPat = re.compile('^SEUSERSUBDIR\s*=\s*"?([^"]*)')
   hnPat  = re.compile('^HN_NAME\s*=\s*"?([^"]*)')
   f = open(path.join(dir,"myScript_0.sh"))
   try:
      for line in f:
         line = line.rstrip() # Remove trailing endline
         m = re.search(sePat,line)
         if m: sehome = m.group(1)
         m = re.search(subPat,line)
         if m: sesubdir = m.group(1)
         m = re.search(hnPat,line)
         if m: hnname = m.group(1)
         if sehome!='' and sesubdir != '' and hnname != '': break
   finally: f.close()
   sepath = path.join(sehome,hnname,sesubdir)
   if opt.verb>2: 
      print '  Found remote home',sehome,', HN name',hnname,'and subdir',sesubdir
      print '  Returning',sepath
   return sepath

def getListOfFiles(dir):
   theList = os.popen('lcg-ls --vo cms --srm-timeout=6000 --connect-timeout=6000 -T srmv2 --nobdii ' + dir + ' 2>&1')
   list = []
   baseurl = dir[0:dir.find('SFN=')]+'SFN='
   # Check files
   for li in theList.readlines():
       if(li.find("root") == -1):
           print '*** Problem with folder:\n',li
           return []
       else:
         list.append(baseurl+li.rstrip())
   return list


def process(dir,opt):
   print "Getting remote directory name",
   sedir = getSeDir(dir,opt)             
   print sedir
   
   print "Getting list of files to remove"
   list = getListOfFiles(sedir)
     
   print "Files to remove:",list
   if opt.dryrun: return

   # 1. Remove remote files
   answ = 'y'
   if len(list)>0:
      cmd = ['srmrm']+list
      if opt.verb>1: print "CMD =",' '.join(cmd)
      if not opt.force: 
         answ = raw_input('Remove files? [y/n] ')
         while answ != 'y' and answ != 'n': 
            answ = raw_input('Please answer y or n: ').rstrip()
      if answ == 'y':
         err = subprocess.call(cmd)
         if err:
            print '*** Error:',err
            sys.exit(-1)
         print "--> Removed remote files"
  
   # 2. Remove remote directory
   if len(list)>0:
      answ = 'y'
      cmd = ['srmrmdir',sedir]
      if opt.verb>1: print "CMD =",' '.join(cmd)
      if not opt.force: 
         answ = raw_input('Remove remote directory? [y/n] ')
         while answ != 'y' and answ != 'n':  
            answ = raw_input('Please answer y or n: ').rstrip()
      if answ == 'y':
         err = subprocess.call(cmd)
         if err:
            print '*** Error:',err
            sys.exit(-1)
         print "--> Removed remote directory"

   # 3. Remove local directory
   answ = 'y'
   cmd = ['rm','-r',directory]
   if opt.verb>1: print "CMD =",' '.join(cmd)
   if not opt.force: 
      answ = raw_input('Remove local directory? [y/n] ')
      while answ != 'y' and answ != 'n': 
         answ = raw_input('Please answer y or n: ').rstrip()
   if answ == 'y':
      err = subprocess.call(cmd)
      if err:
         print '*** Error:',err
         sys.exit(-1)
      print "--> Removed local directory"
  

if __name__ == '__main__' :
   usage = """%prog <dir1> [... <dirn>]"""
   parser = OptionParser(usage=usage)
   parser.add_option("-v","--verbose",dest="verb",metavar="LEVEL",default=0,
                     help="Turn on verbosity")
   parser.add_option("-n","--dry-run",dest="dryrun",action="store_true",default=False,
                     help="Do not actually remove the files: just list them")
   parser.add_option("-f","--force",dest="force",action="store_true",default=False,
                     help="Force removal (don't ask questions)")
   (options,args) = parser.parse_args(args=None, values=None)

   if len(args)<1:
       parser.print_usage()
       sys.exit(-1)

   for directory in args:
      print "+++ Processing directory",directory
      process(directory,options)
