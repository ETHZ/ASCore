# Copies histo files
#   from the NTupleProducer/macros directory 
#   to the PhysQC directory


def copyHistos(myDir = None, mainPageFile = None):
# myDir is the directory where the histos currently reside

    import sys
    import os
    import shutil
    
    # check that a directory name was given
    if myDir == None:
        print " ===> please give the directory name "
	return
    if mainPageFile == None:
        print " ===> please give the file name of the main page "
	return

    # extract the Data Set name
    print " start directory = ", myDir
    ifnd = -1
    for i in range(len(myDir)):
        if myDir[i] == "/":
	    ifnd = i
    newDir = myDir[ifnd+1:len(myDir)]
    print " new directory = ", newDir
    if newDir == myDir:
        print " ===> directory ", newDir, " already copied, will be completed"
    
    # check that the directory does not exist already
    if os.path.exists(newDir) and newDir != myDir:
        print " ===> directory ", newDir, "aready exists"
	answ = raw_input("   do you want to overwrite it? (Y/N) ")
	if answ in ("Y", "y"):
            print "                ", newDir, " will be overwritten"
	    shutil.rmtree(newDir, 1)
	else:
	    return

    # create the new directory
    try:
        os.mkdir(newDir)
    except OSError:
	pass
    
    # copy the files
    if not os.path.exists(newDir + "/checkList.txt"):
        shutil.copy(myDir + "/plots_uncleaned_checklist.txt", newDir + "/checkList.txt")
        f1 = open(newDir + "/checkList.txt", "a")
        f2 = open(myDir + "/plots_Cleaning_checklist.txt", "r")
        f1.write(f2.read())
        f3 = open(myDir + "/plots_cleaned_checklist.txt", "r")
        f1.write(f3.read())
        f1.close()
        f2.close()
        f3.close()
    if os.path.exists("refList.txt"):
        shutil.copy("refList.txt", newDir + "/refList.txt")
    
    if newDir != myDir:
        if os.path.exists(myDir + "/cleanerStats.txt"):
            shutil.copy(myDir + "/cleanerStats.txt", newDir + "/cleanerStats.txt")
        shutil.copy(myDir + "/info.txt", newDir + "/info.txt")
        shutil.copytree(myDir + "/unclean", newDir + "/unclean")
        shutil.copytree(myDir + "/Cleaning", newDir + "/Cleaning")
        shutil.copytree(myDir + "/clean", newDir + "/clean")
        shutil.copytree(myDir + "/eps", newDir + "/eps")
        shutil.copytree(myDir + "/MultiplicityPlots/eps", newDir + "/MultiplicityPlots")
    
    # remove the info file to reinitialize the form
    if os.path.exists("info_files/" + newDir):
        os.remove("info_files/" + newDir)
    
    # update the main html page
    import createPage
    createPage.createPage(newDir, "", mainPageFile)
     

# This is to make copyHistos directly callable

if __name__ == "__main__":
    import sys
    copyHistos(sys.argv[1], sys.argv[2])
