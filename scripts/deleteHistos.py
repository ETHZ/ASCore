# Deletes a histo files directory from the main ETH QC web page


def deleteHistos(myDir = None, mainPageFile = None):
# myDir is the directory where the histos currently reside

    import sys
    import os
    import shutil
    
    # check that a directory name was given
    if myDir == None:
        print " ===> please give the Data Set directory name "
	return
    if mainPageFile == None:
        print " ===> please give the file name of the main page "
	return

    # check for existence of the Data Set name
    print " directory to be deleted = ", myDir, " from ", mainPageFile
    for i in range(len(mainPageFile)):
#        print i, len(mainPageFile), mainPageFile[i]
        if mainPageFile[i] == ".":
	    if mainPageFile[i+1:len(mainPageFile)] == "html":
	        mainPageFile = mainPageFile[:i]
		break
    
    # check that the main page and the directory exist
    if os.path.exists(mainPageFile + ".html"):
        pass
    else:
        print " ===> file ", mainPageFile, "does not exist, try again"
	return
    if os.path.exists(myDir):
        pass
    else:
        print " ===> directory ", myDir, "does not exist, try again"
	return

    # delete the directory
    try:
        shutil.rmtree(myDir)
    except OSError:
	pass
    
    # remove the info file
    if os.path.exists("info_files/" + myDir):
        os.remove("info_files/" + myDir)
    
    # update the main html page
    import mainPageUpdate
    mainPageUpdate.mainPageUpdate(myDir, mainPageFile, "-")
     

# This is to make deleteHistos directly callable

if __name__ == "__main__":
    import sys
    deleteHistos(sys.argv[1], sys.argv[2])
