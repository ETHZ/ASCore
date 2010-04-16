#!/bin/env python
# Copies histo files
#   from the NTupleProducer/macros directory 
#   to the PhysQC directory

from optparse import OptionParser
import sys,os,shutil

def copyHistos(inputDir,mainPageFile,options):
# inputDir is the directory where the histos currently reside

    # check that input directory and file exist
    if not os.path.exists(inputDir):
        print " ===> couldn't find directory",inputDir
	sys.exit(-1)
#    if not os.path.exists(mainPageFile):
#        print " ===> couldn't find main page",mainPageFile
#	sys.exit(-1)

    # Form the output directory (default: same name as input)
    if inputDir.endswith('/'): inputDir = inputDir[:-1] # remove trailing /
    outputDir = os.path.basename(inputDir)
    if options.output_dir: outputDir = options.output_dir

    print " input  directory =",inputDir
    print " output directory =",outputDir

    # Check and create output directory
    if outputDir != inputDir:
        if os.path.exists(outputDir):
            # check that the directory does not exist already
            print " ===> directory ", outputDir, "aready exists"
            answ = raw_input("   do you want to overwrite it? (Y/N) ")
            if answ in ("Y", "y"):
                print "                ", outputDir, " will be overwritten"
                shutil.rmtree(outputDir, 1)
            else:
                sys.exit(-2)
                # create the new directory
        os.mkdir(outputDir)
    else:
        print " ===> directory",outputDir,"already copied, will be completed"
        
    
    # copy the files
    if not os.path.exists(outputDir + "/checkList.txt"):
        dest = open(outputDir + "/checkList.txt", "w")
        shutil.copyfileobj(open(inputDir + "/plots_uncleaned_checklist.txt"),dest)
        shutil.copyfileobj(open(inputDir + "/plots_Cleaning_checklist.txt"),dest)
        shutil.copyfileobj(open(inputDir + "/plots_cleaned_checklist.txt"),dest)
    if os.path.exists("refList.txt"):
        shutil.copy("refList.txt", outputDir + "/refList.txt")

    if outputDir != inputDir:
        if os.path.exists(inputDir + "/cleanerStats.txt"):
            shutil.copy(inputDir + "/cleanerStats.txt", outputDir + "/cleanerStats.txt")
        shutil.copy(inputDir + "/info.txt", outputDir + "/info.txt")
        shutil.copytree(inputDir + "/unclean", outputDir + "/unclean")
        shutil.copytree(inputDir + "/Cleaning", outputDir + "/Cleaning")
        shutil.copytree(inputDir + "/clean", outputDir + "/clean")
#        shutil.copytree(inputDir + "/eps", outputDir + "/eps")
        shutil.copytree(inputDir + "/TriggerStats", outputDir + "/TriggerStats")
        shutil.copytree(inputDir + "/MultiplicityPlots", outputDir + "/MultiplicityPlots")
    
    # remove the info file to reinitialize the form
    if os.path.exists("info_files/" + outputDir):
        os.remove("info_files/" + outputDir)
    
    # update the main html page
    import createPage
    createPage.createPage(outputDir, "", mainPageFile)
     

# This is to make copyHistos directly callable

if __name__ == "__main__":
    usage = """%prog [options] <input directory> <HTML file>
    where:
       <input directory>  is the directory to process
       <HTML file>        is the file to add the run to
       (run --help to see options)"""

    # Command-line options
    parser = OptionParser(usage=usage)
    parser.add_option("-o","--output-dir",dest="output_dir",metavar="DIR",
                      default=None,
                      help="Output directory (default is same as input)")
    (options, args) = parser.parse_args()

    # Check that we have an input directory
    if len(args)<2:
        parser.error('Need at least two arguments')

    # Process arguments
    idir     = args[0]
    mainPage = args[1]

    # Run it!
    copyHistos(idir,mainPage,options)
