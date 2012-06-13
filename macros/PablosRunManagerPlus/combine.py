#!/usr/bin/python
#
#  Combine statistics from several files
# 

import sys,os,re
from optparse import OptionParser
from operator import itemgetter

def usage():
   return "Usage: "+sys.argv[0]+" <file>"
  

def addToCounter(counter,key,count):
    """Increment key entries"""
    # Check if key already exists
    if not key in counter:
        counter[key] = count
    else:
        counter[key] = counter[key] + count

class Counter:
   """
   Counter

   Container for statistics counter
   """
   def __init__(self,title):
      self.title = title
      self.counts = dict()
      self.sum = 0
   def add(self,name,counts):
      #print 'Adding',name,'with counts',counts
      addToCounter(self.counts,name,counts)
   def dump(self):
      print 'COUNTER>',self.title
      format = '%(string)-60s %(#)8d %(##)3.0f%% %(###)6.2f%% '
      sortedCounts = sorted( self.counts.iteritems(), key=itemgetter(1),
                             reverse=True)
      self.sum = sortedCounts[0][1]
      cSum = self.sum
      if cSum>0:
         for k,v in sortedCounts:
            if cSum==0: cSum = 1 # Stupid protection
            print format % {'string':k,'#':v,'##':(v/float(cSum))*100,
                            '###':(v/float(self.sum))*100 }
            cSum = v
         print '-'*80
      
def parseFile( file, counters, indices ):
   """Parse file and increment statistics"""
   try:
       f = open(file)
   except IOError:
       print "File",file,"does not seem to exist"
       print usage()
       sys.exit(3)

   # Patterns to look for
   startPat = re.compile("^COUNTER>\s+(.*)") # FRAGILE: start of new counter
   strPat   = re.compile("^(.*?)\s+(\d+)(\s+\d+\%)(\s+\d+\%)")
   endPat   = re.compile("^-+$")
   
   incounter = False
   cname = ''
   for line in f.readlines():
      # Find start of new counter
      m = re.search(startPat,line)
      if m:
         cname = m.group(1)
         #print 'Found start of counter: ',cname
         if not cname in counters:
            indices.append(cname)
            counters[cname] = Counter(cname)
         incounter = True
      elif incounter: # Reading counter lines
         m = re.search(endPat,line)
         if m: # End of counter
            incounter = False
            continue
         m = re.search(strPat,line)
         if m: # In counter
            # Increment counter's counts (or add to list if new)
            title = m.group(1)
            count = int(m.group(2))
            counters[cname].add(title,count)



#_______________________________________________________________
#### run fill function if invoked on CLI
if __name__ == '__main__':
    
    parser = OptionParser(usage=usage())
    parser.add_option('-i','--include',dest='toInclude',type='string', metavar='LIST', help='Comma-separated list of counters to include')
    parser.add_option('-e','--exclude',dest='toExclude',type='string', metavar='LIST', help='Comma-separated list of counters to exclude')
    (options, args) = parser.parse_args()

    if len(args)<1: parser.error('Need at least one argument')

    # Parse list of counters
    includeList = []
    excludeList = []
    if options.toInclude: includeList = options.toInclude.split(',')
    if options.toExclude: excludeList = options.toExclude.split(',')

    counters = dict()
    indices = []

    for file in args:
      parseFile(file,counters,indices)

    print '-'*80
    for name in indices:
       if (len(includeList)==0 or includeList.count(name)>0) and excludeList.count(name)==0: 
          counters[name].dump()
#    for counter in counters:
#       counters[counter].dump()


