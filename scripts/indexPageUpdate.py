# Updates an html main page with the ETHZ checking

def indexPageUpdate(mainPageFile, addel):
# updates the index.html file by adding/deleting a main Physics Quality Checking page
#     adds or deletes a link in the index page to mainPageFile.html

    # open the index page and stop if no index page exists
    try:
        f = open("index.html", "r+")
    except IOError:
        print "===> index page does not exist, first create one"
	return

    # the updated index page is first written on a tmp file
    tmp = open("temp.html", "w")
    
    myPage = mainPageFile
    if addel == "-":
        print " Main page to be deleted = ", myPage, " from index.html"
    else:
        print " Main page to be added = ", myPage, " to index.html"
    
    # remove the extension if there is one
    for i in range(len(myPage)):
        if myPage[i] == ".":
	    if myPage[i+1:len(myPage)] == "html":
	        myPage = myPage[:i]
		break
    
    # read the old index page line by line
    for line in f:
        # test if the Data Set already existed
	try:
            ind = line.index(myPage + ",")
        except ValueError:
            pass
        else:
	    if addel == "-":
	        print "===> mainPage " + myPage + " deleted "
	    else:
                print "===> mainPage " + myPage + " existed already, " +\
	        "old main page replaced in index page"
	    continue
        # write the new main page in front of all others
	if line[:4] == "<ul>":
	    tmp.write("<ul> \n")
	    if addel != "-":
#	        print "write ", myPage
	        newRun = "  <LI> <B> Data type:    " + myPage + "," +\
                " <A HREF=\"" + myPage + ".html\"> here " + "</A> " +\
	        "</B> </LI> \n"
                tmp.write(newRun)
#                tmp.write("<BR> \n")
        # at the end of the list, write the new footer of the main page
	elif line[:5] == "</ul>":
            tmp.write("</ul> \n")
            tmp.write("<center><BR>")
            tmp.write("<hr noshade=\"noshade\" width=\"100%\"></center> \n")
            import datetime
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
    # replace the old main page with the newly written one
    import os
    os.rename("temp.html", "index.html")
    print "===> index.html updated with " + myPage
     

# This is to make mainPageUpdate directly callable

if __name__ == "__main__":
    import sys
    indexPageUpdate(sys.argv[1], sys.argv[2])

