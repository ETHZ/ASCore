# Creates an html page with the ETHZ checking histograms

def createPage(myDir, histoPrefix = "", mainPageFile = "index"):
# writes the html text for histos in directory myDir
# myDir = directory where the histos reside
# histoPrefix = (optional) prefix for the histo name
# mainPageFile = file name of the main page (default = physQC)

# It is assumed that the histos are residing in the directory myDir.
# Also that the web page will be made of "blocks" of histograms,
# each block corresponding to a given subject (e.g. kinematics, ID variables, ...)
# this structure is described in the file "histoDefs.txt"

# at the end, the main page is updated by linking to this histogram file

# histogram file extension and definition of histo blocks
    ext = "eps"
#FR    extIcon = "jpg"
    extIcon = "png"
    histoDef = "histoDefs.txt"
    nhPerTr = 3

    import os
    import sys
    import glob
    import datetime
    import checkListAnom
    import mainPageUpdate

#   write the page title
    if os.path.exists(myDir) == 0:
        print " CreatePage: directory ", myDir, "does not exist"
	return
    try:
        fhtml = myDir + "/" + myDir + ".html"
        f = open(fhtml, "w")
    except IOError:
        print "I/O error when opening " + fhtml + " for write,\n    EXITING"
        return
    f.write( "<html> \n")
    f.write( "<head> \n")
    f.write( "<title> ETHZ checking histograms </title> \n")
    f.write( "</head> \n")
    f.write( "<body text=\"#111111\" bgcolor=\"#f0f0f0\" " +\
    "link=\"#aa0000\" vlink=\"#aa0000\"> \n")

    f.write("<h1> <font color=\"#aa0000\"><font size=\"+4\">" +\
    "ETHZ checking histograms</font></font></h1>\n")

    f.write( "  \n")
    f.write( "<p> </p> \n")

#   show the directory name
    f.write( "<hr noshade=\"noshade\" width=\"100%\">  \n")
    f.write( "<p><b><font color=\"#0000aa\"><font size=\"+2\"> Data set: " +\
    "<font color =\"#aa0000\">" + myDir + " </font></font></b> \n")
    f.write( "<p> </p> \n")
    f.write( "<hr noshade=\"noshade\" width=\"100%\"> \n")
    f.write("<font size=\"+1\"> ")
    f.write("<p><font color=\"#008000\" + \
    (<i>Click on the image to get a magnified view.</i></font></p> \n")
    f.write( "<p> </p> \n")
    f.write("<p><font color=\"#008000\" + \
    <i>Histograms are ordered according to increasing tightness of cuts </i></font></p> \n")

#   strip off the directory name
#    adir = myDir + "/"
#    nrem = len(adir)
#    for i in range(len(flist)):
#        dum = flist[i]
#        flist[i] = dum[nrem:]
#    print flist

#   open the file with definitions of histogram blocks
    fDef = open(histoDef, "r")
    histo = " "
#    hname = [" ", " ", " ", " ", " "]
    hname = list(nhPerTr*" ")

    #read the file line by line
    for linefull in fDef:
        n = linefull.find("\n")
        if n <= 0: continue
        line = linefull[:n]
#        print line
        # break the loop when the end is reached
	if line == "* end": break
	
	# if a new block title is found,
        elif line[0:2] == "* ":
            if histo != " ": 
                f.write( "    </tr> \n")
		
		#write the histo names of previous block if some are left over
		if nnames >= 0:
		    f.write("    <tr align=\"center\"> \n")
		    for i in range(nnames+1):
		        f.write("    <td><B><font color=\"#00aa00\">" +\
			hname[i] + "</font></B></td> \n")
                f.write( " <BR> \n")
                f.write( "    </tbody> \n")
                f.write( "</table> \n")
            
	    # write the block title as a new header
	    title = line[2:]
            f.write( "  \n")
            f.write( "<p> </p> \n")
            f.write( "<hr noshade=\"noshade\" width=\"100%\">  \n")
            f.write( "<p><b><font color=\"#0000aa\"><font size=\"+2\">" +\
            title +"</font></font></b> \n")
            f.write( " <p>  </p> \n")
            f.write( "<table border=\"3\" cellspacing=\"3\" cellpadding=\"10\">\n")
            f.write( "    <tbody> \n")
            icount = 0
	    nnames = -1
	
	# else if a new subdirectory is found,
        elif line[0:2] == "**":
            subdir = "/" + line[3:]
	    print subdir
            if os.path.exists(myDir) == 0:
                print " wrong subdirectory", subdir, "in histoDefs"
    
#           create the list of available histograms in the subdirectory
            flist = glob.glob(myDir + subdir + "/*." + ext)
#            print flist
        
	# if a new histogram name is found,
	else:
            icount = icount + 1
	    # allow only a limited number of histos on the same table row
            if icount % nhPerTr == 1:
                if icount > 1:
                    f.write( "    </tr> \n")
		    
		    # write the histo names of previous row
		    f.write("    <tr align=\"center\"> \n")
		    for i in range(nnames+1):
		        f.write("    <td><B><font color=\"#00aa00\">" +\
			hname[i] + "</font></B></td> \n")
                    f.write("    </tr> \n")
		    nnames = -1
		f.write( "    <tr align=\"center\"> \n")
            histo = histoPrefix + line + "." + ext
	    icon  = histoPrefix + line + "." + extIcon
            
	    # check that the histo exists in the subdirectory
            ffound = ""
	    try:
                ind = flist.index(myDir + subdir + "/" + histo)
            #if not, make an empty entry
	    except ValueError:
                ffound = ", " + histo + " not found"
                f.write("    <td> <font color=\"#aa0000\">" + histo +\
                " </font> </td> \n")
            # put the histo and its small icon in the table
	    else:
                if os.path.isfile(myDir + subdir + "/" + icon):
                    pass
                else:
#                    print "file " + icon + " does not exist"
                    fIcon = myDir + subdir + "/" + icon
                    fHisto = myDir + subdir + "/" + histo
                    print "convert " + fHisto + " to " + fIcon
                    os.system("convert " + fHisto + " " + fIcon)
                    print "created file ",fIcon
                f.write( "    <td><a href= " + subdir[1:] + "/" + histo + " size=\"2\" " +\
                "<img src= " + subdir[1:] + "/" + icon +\
		" height=\"300\" width=\"280\" align=\"left\" </a></td> \n")
            nnames = nnames + 1
            hname[nnames] = line
#            print icount, nnames, hname[nnames], ffound

    f.write( "    </tr> \n")
    
    # when the table is completed, write the histo names of the last row
    f.write("    <tr align=\"center\"> \n")
    for i in range(nnames+1):
        f.write("    <td><B><font color=\"#00aa00\">" +\
	hname[i] + "</font></B></td> \n")
    f.write("    </tr> \n")
    
    f.write( "  \n")
    f.write( "    </tbody> \n")
    f.write( "</table> \n")
    fDef.close()

    f.write( "  \n")
    
    f.write("<center><BR>")
    f.write("<hr noshade=\"noshade\" width=\"100%\"></center>")
    now = datetime.datetime.now()
    dat = now.strftime("%d %b %Y")
    tim = now.strftime("%H:%M:%S")
    f.write("<font color=\"#AA0000\">")
    f.write("<center> Last update: " + dat + ",   " + tim + "</center>")
    f.write("<BR>")
    f.write("<center> ETHZ SUSY Group </center>")
    f.write( "<BR><BR>")

    f.write( "</body> \n")
    f.write( "</html> \n")

#   create the page with anomalies
    checkListAnom.checkListAnom(myDir)

#   update the main page with this new run histo page
    mainPageUpdate.mainPageUpdate(myDir, mainPageFile, "+")

# This is to make createPage directly callable

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) == 4 :
    	createPage(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
    	createPage(sys.argv[1] , sys.argv[2])
###################################################################################


            
     

