# Creates an html main page with the ETHZ checking


def mainPageCreate(fileName):
# writes the html text for the main Physics Quality Checking page

    # remove the html extension if there is one
    for i in range(len(fileName)):
        if fileName[i] == ".":
	    if fileName[i+1:len(fileName)] == "html":
	        fileName = fileName[:i]
		break
    
    # open the file
    f = open(fileName + ".html", "w")

    # write the page title
    f.write( "<html>\n")
    f.write( "<head>\n")
    f.write( "<title> ETHZ Physics Quality Checking </title>\n")
    f.write( "</head>\n")
    f.write( "<body text=\"#111111\" bgcolor=\"#f0f0f0\" " +\
    "link=\"#aa0000\" vlink=\"#aa0000\">\n")

    f.write("<h1> <font color=\"#aa0000\"><font size=\"+4\">" +\
    "ETHZ Physics Quality Checking</font></font></h1>\n")

    f.write( " \n")
    f.write( "<p> </p>\n")

    f.write( "<hr noshade=\"noshade\" width=\"100%\"> \n")
    f.write( "<p><b><font color=\"#0000aa\"><font size=\"+2\"> " +\
    "Histograms by Data Set: </font></b> \n")
    f.write( "<p> </p>\n")

    f.write( "<ul>\n")

    f.write( "</ul>\n")

    f.write( " \n")
    f.write( "</body>\n")
    f.write( "</html>\n")

    f.close()

    print "===> New main page created"
    
    #add this main page to the index.html page
    import indexPageUpdate
    indexPageUpdate.indexPageUpdate(fileName, "+")

# This is to make mainPageCreate directly callable

if __name__ == "__main__":
    import sys
    mainPageCreate(sys.argv[1])
 
###################################################################################



            
     

