# Looks for anomalies in the checkList results

def checkListAnom(myDir):
# creates a text file with anomalous values in the checkList

# adds a link in the main page mainPageFile.html
# to the run histos file myDir.html

    import os
    import shutil

# define the number of s.d. for calling a difference anomalous
    anom = 3.

    # verify whether the checklists have been collected into checkList.txt
    # if not, the pieces of checkLists are supposed to exist in this subdirectory
    if os.path.exists(myDir + "/checkList.txt"):
        pass
    else:
        if os.path.exists(myDir + "/plots_uncleaned_checklist.txt"):
	    pass
	else:
	    print "===> No checklists available, no anomalies file produced"
	    return
	shutil.copy(myDir + "/plots_uncleaned_checkList.txt", myDir + "/checkList.txt")
	f1 = open(myDir + "/checkList.txt", "a")
	f2 = open(myDir + "/plots_Cleaning_checkList.txt", "r")
	f1.write(f2.read())
	f3 = open(myDir + "/plots_cleaned_checkList.txt", "r")
	f1.write(f3.read())
	f1.close()
	f2.close()
	f3.close()
    
    # open the checkList.txt file
    try:
        f = open(myDir + "/checkList.txt", "r")
    except IOError:
        print "===> ", myDir + "/checkList.txt not found"
        return

    # open the anomalies.txt file
    tmp = open(myDir + "/anomalies.txt", "w")
    tmp.write("\n")
    tmp.write(" Anomalous values in file: " + myDir + "/checkList.txt \n")
    tmp.write("\n")
    
    # check whether the refList exists and, if not, copy it over
    try:
        r = open(myDir + "/refList.txt", "r")
    except IOError:
        try:
	    rsv = open("refList.txt", "r")
	except IOError:
            tmp.write(" No refList found, no anomaly check is made \n")
            tmp.close()
	    print "===> refList not found, no anomaly check is made"
            return
	shutil.copy("refList.txt", myDir + "/refList.txt")
	r = open(myDir + "/refList.txt", "r")
    
#   loop over the lines in the checkList
    rposb = 0
    var = ""
    anfnd = 0
    varprnt = ""
    for line in f:
        if line == "":
            print " EOF reached on checkList"
            break
#       skip blank lines
        elif line[0:2] == "\n":
            continue
#       look for a new variable name on checkList (line starts with a *)
        elif line[0:1] == "*":
            var = line[2:]
#            print " checkList found " + var
#       find the number to be checked, after a variable is found
        elif var != "":
            if line[0:2] == "\n":
                continue
            try:
                ista = line.index("=")
            except ValueError:
#                print "*** Going crazy on checkList, no = found"
                continue
            try:
                ichk = line.index("+-")
            except ValueError:
#                print "*** Going crazy on checkList, no +- found"
                continue
#           now search for the corresponding value on the refList
            varFound = 0
            valueFound = 0
            r.seek(rposb)
            for ref in r:
                if line[0:2] == "\n":
                    print " blank line on refList skipped"
                    continue
                elif line == "":
                    print " EOF reached on refList"
                    break
                elif ref[0:1] == "*" and ref[2:] == var:
#                    print " found " + var[:len(var)-2] + " on refList"
                    varFound = 1
                elif varFound:
                    if line[:ista] != ref[:ista]:
                        continue
                    try:
                        istaref = ref.index("=")
                    except ValueError:
                        print "*** Going crazy on refList, no = found"
                        continue
                    try:
                        iref = ref.index("+")
                    except ValueError:
                        print "*** Going crazy on refList, no +- found"
                        continue
                    valueFound = 1
                    break
            if valueFound == False:
                tmp.write("*** No value found for " + var[:len(var)-1] +\
                          ", " + line[:ista] + " on refList, skipped \n")
                continue
            cvalchk = line[ista+1:ichk-1]
            cerrchk = line[ichk+2:]
#            print cvalchk + " +- " + cerrchk
            valchk = eval(cvalchk)
            errchk = eval(cerrchk)
            cvalref = ref[ista+1:iref-1]
#            print cvalref
            valref = eval(cvalref)
	    signif = 999.
	    if errchk > 0.:
                signif = abs(valref - valchk) / errchk
#            print signif
            if signif > anom:
                if varprnt != var:
		    tmp.write("* " + var )
		    varprnt = var
		nchar = len(line)
                tmp.write(line[:nchar-2] + " is anomalous, refList = " +\
                          cvalref + "\n")
		anfnd = 1
#                return
        else:
            print " Could not find a meaningful line on checkList"
        
    if anfnd == 0:
        tmp.write(" No anomaly found")
    f.close()
    r.close()
    tmp.close()
    print "===> Created " + myDir + "/anomalies.txt from " + myDir + "/checkList.txt" +\
          " and " + myDir + "/refList.txt"
     

# This is to make checkListAnom directly callable

if __name__ == "__main__":
    import sys
    checkListAnom(sys.argv[1])

