from string import Template
from mod_python import util
import os

import datetime
def Save(req, dsName,  name, comment ):

    status = req.form.getfirst('status')

    # make sure the user provided all the parameters
    if not (dsName and status and name):
        return "A required parameter is missing, \
               please go back and correct the error"

    directory = os.path.dirname(req.filename)
    outfile = os.path.join(directory,'info_files/' + dsName)

    f=open(outfile, 'w')
    f.write(status + "\n")
    now = datetime.datetime.now()
    f.write( now.strftime("%d %b %Y, %H:%M:%S") + '\n' )
    f.write(name+ "\n")
    f.write(comment)
    f.close()

    util.redirect(req, 'Load?dsname=' + dsName)

    return outfile


def Load(req, dsname):
    req.content_type = "text/html"
    req.send_http_header()

    directory = os.path.dirname(req.filename)


    htmlfile = os.path.join(directory,'info_files/TheForm.html')
    formFile = open(htmlfile , 'r')

    rettem = Template(formFile.read())

    formFile.close()

    timeOfTheRun = '-1'
    nEvents = 'N/A'

    infofile = os.path.join(directory,dsname + '/info.txt')
    try:
        info = open(infofile , 'r')

        timeOfTheRun = info.readline()
        info.readline()
        nEvents = info.readline()

        info.close()
    except:
        nEvents = 'N/A'

    infofile = os.path.join(directory,'info_files/' + dsname )
    try:
        info = open(infofile , 'r')

        status = info.readline()

        GOODCHECKED = ''
        DUBIOUSCHECKED = ''
        BADCHECKED = ''
        if  str(status)[0:-1]==str('GOOD'):
            GOODCHECKED = 'checked=\"checked\"'
        elif  str(status)[0:-1]=='BAD':
            BADCHECKED = 'checked=\"checked\"'
        elif  str(status)[0:-1]=='DUBIOUS':
            DUBIOUSCHECKED = 'checked=\"checked\"'

        modificatedAt = info.readline()    
        name = info.readline()
        comment = info.readline()
        line = info.readline()
        while line:
            comment += line
            line = info.readline()

        info.close()

        allVals = dict(GOODCHECKED=GOODCHECKED,BADCHECKED=BADCHECKED,DUBIOUSCHECKED=DUBIOUSCHECKED,comment=comment,name=name,DSNAME=dsname,TIME=timeOfTheRun,NEvents=nEvents,LastUpdate=modificatedAt)
        
        return str(rettem.safe_substitute( allVals ))

    except IOError:
        allVals = dict(GOODCHECKED='',BADCHECKED='',DUBIOUSCHECKED='',comment='',name='',DSNAME=dsname,TIME=timeOfTheRun,NEvents=nEvents,LastUpdate='Never')

        return str(rettem.safe_substitute( allVals ))

    return str('')


from mod_python import apache
def GetImage(req, dsname):
    req.content_type = "image/gif"
    req.send_http_header()

    directory = os.path.dirname(req.filename)
    infofile = os.path.join(directory,'info_files/' + dsname )

    try:
        info = open(infofile , 'r')

        status = str(info.readline())[0:-1]
        info.close()

        if status=='GOOD':
            util.redirect(req, '../img/green.gif')
        elif status=='BAD':
            util.redirect(req, '../img/red.gif')
        elif status=='DUBIOUS':
            util.redirect(req, '../img/yellow.gif')
        else:
            util.redirect(req, '../img/blue.gif')
    except IOError:
        util.redirect(req, '../img/white.gif')

    return apache.OK
