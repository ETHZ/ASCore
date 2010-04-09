# Deletes an html main page with the ETHZ checking

def deleteMainPage(mainPageFile):
# deletes the main Physics Quality Checking page for a given type of data
#     deletes the link from the index.html page to the mainPageFile.html 

    import os
    import deleteHistos

    # remove the extension of mainPageFile if there is one
    myMain = mainPageFile
    for i in range(len(myMain)):
        if myMain[i] == ".":
	    if myMain[i+1:len(myMain)] == "html":
	        myMain = myMain[:i]
		break
    
    # open the old main page
    try:
        f = open(myMain + ".html", "r+")
    except IOError:
        print "===> main page = ", myMain + ".html does not exist, try again"
	return

    # read the old main page line by line
    icopy = 0
    save = []
    n = -1
    strt = 24
    for line in f:
        # collect all existing Data Sets
#	print line[:40]
	if line[:4] == "<ul>":
	    icopy = 1
	elif line[:5] == "</ul>":
            break
        elif icopy > 0:
#	    print line[:40]
	    n = n + 1
	    ind = line[strt:].index(",")
	    save[n:n] = [line[strt:strt+ind]]
#    print save
    f.close()
    
    # delete all histos from this main page
    for i in range(len(save)):
        deleteHistos.deleteHistos(save[i], myMain)
    
    # delete the old main page
    os.remove(myMain + ".html")
    
    #delete this main page from the index.html
    import indexPageUpdate
    indexPageUpdate.indexPageUpdate(myMain, "-")
    
    print "===> Mainpage deletion = " + myMain + " completed"
     

# This is to make deleteMainPage directly callable

if __name__ == "__main__":
    import sys
    deleteMainPage(sys.argv[1])

