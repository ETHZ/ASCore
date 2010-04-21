# Updates an html main page with the ETHZ checking

def mainPageUpdate(myDir, mainPageFile, addel):
# updates the html text for the main Physics Quality Checking page

# if addel = "+"
#     adds a link in the main page mainPageFile.html
#     to the run histos file myDir.html
#     if the main page does not exist, it creates one
# if addel = "-"
#     deletes the link from the main page mainPageFile.html 

    import os
    import datetime
    import mainPageCreate
    
    # remove the extension of mainPageFile if there is one
    myMain = mainPageFile
    for i in range(len(myMain)):
        if myMain[i] == ".":
	    if myMain[i+1:len(myMain)] == "html":
	        myMain = myMain[:i]
		break
    
    # create the header of the main page if no main page existed
    try:
        f = open(myMain + ".html", "r+")
    except IOError:
        if addel == "-":
	    return
	else:
            print "===> also create the main page"
            mainPageCreate.mainPageCreate(myMain)
            f = open(myMain + ".html", "r+")

    # the updated main page is first written on a tmp file
    tmp = open("temp.html", "w")
    
    # read the old main page line by line
    for line in f:
        # test if the Data Set already existed
	try:
            ind = line.index(myDir + ",")
        except ValueError:
            pass
        else:
	    if addel == "-":
	        print "===> Run " + myDir + " deleted "
	    else:
                print "===> Run " + myDir + " existed already, " +\
	        "old Data Set removed from main page"
	    continue
        # write the new Data Set in front of all others
	if line[:4] == "<ul>":
	    tmp.write("<ul> \n")
	    if addel != "-":
#	        print "write ", myDir
	        newRun = "  <LI> <B> Data set:    " + myDir +\
                ", <A HREF=\"FormHandler.py/Load?dsname=" +\
                 myDir + "\"><img src=FormHandler.py/GetImage?dsname=" + myDir + " width=20 height=20 />  info " + " </A>,  " +\
                " <A HREF=\"" + myDir + "/" +\
                myDir + ".html\"> histos " + "</A>, " +\
	        " <A HREF=\"" + myDir + "/" +\
	        "checkList.txt\"> checkList </A>, " +\
	        " <A HREF=\"" + myDir + "/" +\
	        "anomalies.txt\"> anomalies </A>, " +\
	        " <A HREF=\"" + myDir + "/" +\
	        "cleanerStats.txt\"> cleaner stats </A>" +\
	        "</B> </LI> \n"
                tmp.write(newRun)
#                tmp.write("<BR> \n")
        # at the end of the list, write the new footer of the main page
	elif line[:5] == "</ul>":
            tmp.write("</ul> \n")
            tmp.write("<center><BR>")
            tmp.write("<hr noshade=\"noshade\" width=\"100%\"></center> \n")
            now = datetime.datetime.now()
            dat = now.strftime("%d %b %Y")
            tim = now.strftime("%H:%M:%S")
            tmp.write("<font color=\"#AA0000\">")
            tmp.write("<center> Last update: " + dat + ",   " + tim + "</center>\n")
            tmp.write("<BR> \n")
            tmp.write("<center> ETHZ SUSY Group </center> \n")
            tmp.write( "<BR><BR> \n")

            tmp.write( "</body> \n")
            tmp.write( "</html> \n")
            break
        # for all other lines, copy the old to the new one
	else:
            tmp.write(line)

    f.close()
    tmp.close()
    os.chmod(tmp.name,stat.S_IWGRP|stat.S_IRGRP) # Make sure it is group writable
    # replace the old main page with the newly written one
    os.rename("temp.html", myMain + ".html")
    print "===> Mainpage updated with " + myDir
     

# This is to make mainPageUpdate directly callable

if __name__ == "__main__":
    import sys
    mainPageUpdate(sys.argv[1], sys.argv[2], sys.argv[3])

