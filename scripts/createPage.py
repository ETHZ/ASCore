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
    f.write( '<link type="text/css" href="../main.css" REL=stylesheet>')
    f.write( "<title> ETHZ checking histograms </title> \n")
    f.write( "</head> \n")
    f.write( "<body text=\"#111111\" bgcolor=\"#f0f0f0\" " +\
    "link=\"#aa0000\" vlink=\"#aa0000\"> \n")

    f.write("<h1> <font color=\"#aa0000\"><font size=\"+4\">" +\
    "ETHZ checking histograms</font></font></h1>\n")

    f.write( "  \n")
    f.write( "<p> </p> \n")

    # Global navigation
    f.write( "<hr noshade=\"noshade\" width=\"100%\">  \n")
    f.write( '<DIV class="global_nav"><A HREF="../index.html">PhysQC</A>')
    pageName = os.path.splitext(os.path.split(mainPageFile)[1])[0] # Extract name of file, without extension
    f.write( ' &gt; <A HREF="../'+mainPageFile+'">'+pageName+'</A>');
    f.write( ' &gt; '+myDir+'</DIV>');
    
    # Show the directory name
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

    # Local navigation menu
    localMenu = '<div id="menu">'

    #   open the file with definitions of histogram blocks
    fDef = open(histoDef, "r")
    histo = " "
    hname = list(nhPerTr*" ")

    # read the file line by line
    iblock = 0
    for linefull in fDef:
        n = linefull.find("\n")
        if n <= 0: continue
        line = linefull[:n]
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
			hname[i] + " &uarr;</font></B></td> \n")
                f.write( " <BR> \n")
                f.write( "    </tbody> \n")
                f.write( "</table> \n")
            
	    # write the block title as a new header
	    title = line[2:]
            iblock = iblock+1
            f.write( "  \n")
            f.write( "<p> </p> \n")
            f.write( '<hr noshade="noshade" width="100%">'+"\n")
            f.write( '<a name="block'+str(iblock)+'"></a>' ) # Local navigation
            f.write( '<p><b><font color="#0000aa"><font size="+2">')
            f.write(  title +'</font></font></b></p>'+"\n")
            f.write( " <p>  </p> \n")
            f.write( '<table border="3" cellspacing="3" cellpadding="10">'+"\n")
            f.write( "    <tbody> \n")
            icount = 0
	    nnames = -1
            # Add to local navigation menu
            localMenu += '<a href="#block'+str(iblock)+'">'+title+'</a><br><br>'+"\n"
            
	
	# else if a new subdirectory is found,
        elif line[0:2] == "**":
            subdir = "/" + line.lstrip('* ')
	    print subdir
            if os.path.exists(myDir) == 0:
                print " wrong subdirectory", subdir, "in histoDefs"
            
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
			hname[i] + " &uarr;</font></B></td> \n")
                    f.write("    </tr> \n")
		    nnames = -1
		f.write( "    <tr align=\"center\"> \n")

            # Form file name and absolute path of plots
	    icon  = histoPrefix + line + "." + extIcon
            histo = histoPrefix + line + "." + ext
            fIcon  = myDir + subdir + "/" + icon
            fHisto = myDir + subdir + '/' + histo
            
	    # check that the histos exist in the subdirectory
            hasIcon  = os.path.isfile(fIcon)
            hasHisto = os.path.isfile(fHisto)
            # put the histo and its small icon in the table
            relsubdir = subdir.lstrip('/')
            f.write('    <td>')
            if hasHisto:
                f.write('<a href="'+relsubdir+'/'+histo+'" size="2">')
            else:
                print "WARNING:",fHisto,"not found"
            if hasIcon:
                f.write('<img src= "'+relsubdir+'/'+icon+'" height="300" width="280" align="left" >')
            else:
                print "WARNING:",fIcon,"not found"                
                f.write('<font color="#aa0000">' + histo + '</font>')
            if hasHisto: f.write('</a>')
            f.write('</td>'+"\n")
                
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

    localMenu += '</div>'

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
    f.write(localMenu)
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


            
     

