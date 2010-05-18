#!/usr/bin/env python
#____________________________________________________________
#
#  BeamSpotWorkflow
#
# A very complicate way to automatize the beam spot workflow
#
# Francisco Yumiceva, Lorenzo Uplegger
# yumiceva@fnal.gov, uplegger@fnal.gov
#
# Fermilab, 2010
#
#____________________________________________________________

"""
   BeamSpotWorkflow.py

   A very complicate script to upload the results into the DB

   usage: %prog -d <data file/directory> -t <tag name>
   -o, --overwrite : Overwrite results files when copying.
   -T, --Test      : Upload files to Test dropbox for data validation.   
   -u, --upload    : Upload files to offline drop box via scp.
   -z, --zlarge    : Enlarge sigmaZ to 10 +/- 0.005 cm.

   Francisco Yumiceva (yumiceva@fnal.gov)
   Lorenzo Uplegger   (send an email to Francisco)
   Fermilab 2010
   
"""


import sys,os
import commands, re, time
import datetime
import ConfigParser
import xmlrpclib
from BeamSpotObj import BeamSpot
from IOVObj import IOV
from CommonMethods import *

try: # FUTURE: Python 2.6, prior to 2.6 requires simplejson
    import json
except:
    try:
        import simplejson as json
    except:
        print "Please set a crab environment in order to get the proper JSON lib"
        sys.exit(1)

#####################################################################################
# General functions
#####################################################################################
def getLastUploadedIOV(tagName,destDB="oracle://cms_orcoff_prod/CMS_COND_31X_BEAMSPOT"):
    listIOVCommand = "cmscond_list_iov -c " + destDB + " -P /afs/cern.ch/cms/DB/conddb -t " + tagName 
    aCommand       = listIOVCommand + " | grep DB= | tail -1 | awk \'{print $1}\'"
#    print " >> " + aCommand
    output = commands.getstatusoutput( aCommand )
    if output[0] != 0 :
        dbError = commands.getstatusoutput( listIOVCommand )
        if dbError[1].find("metadata entry \"" + tagName + "\" does not exist"):
            print "Creating a new tag because I got the following error contacting the DB"
            print output[1]
            return 1
            #return 133928
        else:
            exit("ERROR: Can\'t connect to db because:\n" + dbError[1])
    #WARNING when we pass to lumi IOV this should be long long
    return long(output[1])

########################################################################
def getListOfFilesToProcess(dataSet,lastRun=-1):
    queryCommand = "dbs --search --query \"find file where dataset=" + dataSet
    if lastRun != -1:
        queryCommand = queryCommand + " and run > " + str(lastRun)
    queryCommand = queryCommand + "\" | grep .root"    
#    print " >> " + queryCommand
    output = commands.getstatusoutput( queryCommand )
    return output[1].split('\n')

########################################################################
def getListOfRunsAndLumiFromDBS(dataSet,lastRun=-1):
    queryCommand = "dbs --search --query \"find run,lumi where dataset=" + dataSet
    if lastRun != -1:
        queryCommand = queryCommand + " and run > " + str(lastRun)
    queryCommand += "\""
    print " >> " + queryCommand
    output = commands.getstatusoutput( queryCommand )
    outputList = output[1].split('\n')
    runsAndLumis = {}
    for out in outputList:
        regExp = re.search('(\d+)\s+(\d+)',out)
        if regExp:
            run  = long(regExp.group(1))
            lumi = long(regExp.group(2))
            if not run in runsAndLumis:
                runsAndLumis[run] = []
            runsAndLumis[run].append(lumi)

#    print runsAndLumis
#    exit("ok")
    return runsAndLumis

########################################################################
def getListOfRunsAndLumiFromRR(dataSet,lastRun=-1):
    RunReg  ="http://pccmsdqm04.cern.ch/runregistry"
    #RunReg  = "http://localhost:40010/runregistry"
    #Dataset=%Online%
    Group   = "Collisions10"

    # get handler to RR XML-RPC server
    FULLADDRESS=RunReg + "/xmlrpc"
    #print "RunRegistry from: ",FULLADDRESS
    server = xmlrpclib.ServerProxy(FULLADDRESS)
    #sel_runtable="{groupName} ='" + Group + "' and {runNumber} > " + str(lastRun) + " and {datasetName} LIKE '" + Dataset + "'"
    sel_runtable="{groupName} ='" + Group + "' and {runNumber} > " + str(lastRun) 
    sel_dcstable="{groupName} ='" + Group + "' and {runNumber} > " + str(lastRun) + " and {parDcsBpix} = 1 and {parDcsFpix} = 1 and {parDcsTibtid} = 1 and {parDcsTecM} = 1 and {parDcsTecP} = 1 and {parDcsTob} = 1 and {parDcsEbminus} = 1 and {parDcsEbplus} = 1 and {parDcsEeMinus} = 1 and {parDcsEePlus} = 1 and {parDcsEsMinus} = 1 and {parDcsEsPlus} = 1 and {parDcsHbheA} = 1 and {parDcsHbheB} = 1 and {parDcsHbheC} = 1 and {parDcsH0} = 1 and {parDcsHf} = 1"

    tries = 0;
    while tries<10:
        try:
            run_data = server.DataExporter.export('RUN'           , 'GLOBAL', 'csv_runs', sel_runtable)
            dcs_data = server.DataExporter.export('RUNLUMISECTION', 'GLOBAL', 'json'    , sel_dcstable)
            break
        except:
            print "Something wrong in accessing runregistry, retrying in 5s...."
            tries += 1
            time.sleep(5)
        if tries==10:
            error = "Run registry unaccessible.....exiting now"
            sys.exit(error)
    

    listOfRuns=[]
    for line in run_data.split("\n"):
        run=line.split(',')[0]
        if run.isdigit():
            listOfRuns.append(run)


    selected_dcs={}
    jsonList=json.loads(dcs_data)

    for element in jsonList:
        if element in listOfRuns:
            selected_dcs[long(element)]=jsonList[element]

    return selected_dcs

########################################################################
def getLastClosedRun(DBSListOfFiles):
    runs = []
    for file in DBSListOfFiles:
        runNumber = getRunNumberFromDBSName(file)
        if runs.count(runNumber) == 0: 
            runs.append(runNumber)

    if len(runs) <= 1: #No closed run
        return -1
    else:
        runs.sort()
        return long(runs[len(runs)-2])
    
########################################################################
def getRunNumberFromFileName(fileName):
    regExp = re.search('(\D+)_(\d+)_(\d+)_(\d+)',fileName)
    if not regExp:
        return -1
    return long(regExp.group(3))

########################################################################
def getRunNumberFromDBSName(fileName):
    regExp = re.search('(\D+)/(\d+)/(\d+)/(\d+)/(\D+)',fileName)
    if not regExp:
        return -1
    return long(regExp.group(3)+regExp.group(4))
    
########################################################################
def getNewRunList(fromDir,lastUploadedIOV):
    newRunList = []
    listOfFiles = ls(fromDir)
    runFileMap = {}
    for fileName in listOfFiles:
        runNumber = getRunNumberFromFileName(fileName) 
        if runNumber > lastUploadedIOV:
            newRunList.append(fileName)
    return newRunList        

########################################################################
def selectFilesToProcess(listOfRunsAndLumiFromDBS,listOfRunsAndLumiFromRR,newRunList,runListDir,tolerance):
    runsAndLumisProcessed = {}
    runsAndFiles = {}
    for fileName in newRunList:
        file = open(runListDir+fileName)
        for line in file:
            if line.find("Runnumber") != -1:
                run = long(line.replace('\n','').split(' ')[1])
            elif line.find("LumiRange") != -1:
                lumiLine = line.replace('\n','').split(' ')
                begLumi = long(lumiLine[1])
                endLumi = long(lumiLine[3])
                if begLumi != endLumi:
                    error = "The lumi range is greater than 1 for run " + str(run) + " " + line + " in file: " + runListDir + fileName
                    exit(error)
                else:
                    if not run in runsAndLumisProcessed:
                        runsAndLumisProcessed[run] = []
                    if begLumi in runsAndLumisProcessed[run]:
                        print "Lumi " + str(begLumi) + " in event " + str(run) + " already exist. This MUST not happen but right now I will ignore this lumi!"
                    else:    
                        runsAndLumisProcessed[run].append(begLumi)
        if not run in runsAndFiles:
            runsAndFiles[run] = []
        runsAndFiles[run].append(fileName)    
        file.close()
    rrKeys = listOfRunsAndLumiFromRR.keys()
    rrKeys.sort()
    dbsKeys = listOfRunsAndLumiFromDBS.keys()
    dbsKeys.sort()
    procKeys = runsAndLumisProcessed.keys()
    procKeys.sort()
    #print "Run Registry:"    
    #print rrKeys
    #print "DBS:"    
    #print dbsKeys
    #print "List:"    
    #print procKeys
    filesToProcess = []
    #print "WARNING WARNING YOU MUST CHANGE THIS LINE BEFORE YOU CAN REALLY RUN THE SCRIPT!!!!!!"
    for run in rrKeys:
    #for run in dbsKeys:
        if run in runsAndLumisProcessed:
            if not run in listOfRunsAndLumiFromDBS:
                error = "Impossible but run " + str(run) + " has been processed and it is also in the run registry but it is not in DBS!" 
                exit(error)
            print "Working on run " + str(run)    
            errors          = []
            badProcessed    = []
            badDBSProcessed = []
            badDBS          = []
            badRRProcessed  = []
            badRR           = []
            #It is important for runsAndLumisProcessed[run] to be the first because the comparision is not ==
            badDBSProcessed,badDBS = compareLumiLists(runsAndLumisProcessed[run],listOfRunsAndLumiFromDBS[run],errors)
            for i in range(0,len(errors)):
                errors[i] = errors[i].replace("listA","the processed lumis")
                errors[i] = errors[i].replace("listB","DBS")
            #print errors
            #print badProcessed
            #print badDBS
            #exit("ciao")
            if len(badDBS) != 0:
                print "This is weird because I processed more lumis than the ones that are in DBS!"
            if len(badDBSProcessed) != 0 and run in listOfRunsAndLumiFromRR:
                lastError = len(errors)
                RRList = []
                for lumiRange in listOfRunsAndLumiFromRR[run]:
                    #print lumiRange
                    for l in range(lumiRange[0],lumiRange[1]+1):
                        RRList.append(long(l))
                #It is important for runsAndLumisProcessed[run] to be the first because the comparision is not ==
                badRRProcessed,badRR = compareLumiLists(runsAndLumisProcessed[run],RRList,errors)
                for i in range(0,len(errors)):
                    errors[i] = errors[i].replace("listA","the processed lumis")
                    errors[i] = errors[i].replace("listB","Run Registry")
                #print errors
                #print badProcessed
                #print badRunRegistry
                    
                if len(badRRProcessed) != 0:    
                    print "I have not processed some of the lumis that are in the run registry for run: " + str(run)
                    for lumi in badDBSProcessed:
                        if lumi in badRRProcessed:
                            badProcessed.append(lumi)
                    lenA = len(badProcessed)
                    lenB = len(RRList)
                    if 100.*lenA/lenB < tolerance:
                        print "I didn't process " + str(100.*lenA/lenB) + "% of the lumis but I am within the " + str(tolerance) + "% set in the configuration."
                        badProcessed = []
                    else:
                        print "I didn't process " + str(100.*lenA/lenB) + "% of the lumis and I am not within the " + str(tolerance) + "% set in the configuration."  
                        exit(errors)
                    #return filesToProcess
                elif len(errors) != 0:
                    print "The number of lumi sections processed didn't match the one in DBS but they cover all the ones in the Run Registry, so it is ok!"
                    #print errors

            #If I get here it means that I passed or the DBS or the RR test            
            if len(badProcessed) == 0:
                for file in runsAndFiles[run]:
                    filesToProcess.append(file)
            else:
                #print errors
                print "This should never happen because if I have errors I return or exit! Run: " + str(run)
        else:
            error = "Run " + str(run) + " is in the run registry bit it has not been processed!Lore"
            print error
            #exit(error)
                            
    return filesToProcess
########################################################################
def compareLumiLists(listA,listB,errors=[],tolerance=0):
    lenA = len(listA)
    lenB = len(listB)
    if lenA < lenB-(lenB*float(tolerance)/100):
        errors.append("ERROR: The number of lumi sections is different: listA(" + str(lenA) + ")!=(" + str(lenB) + ")listB")
    #else:
        #errors.append("Lumi check ok!listA(" + str(lenA) + ")-(" + str(lenB) + ")listB")
    #print errors
    listA.sort()
    listB.sort()
    #shorter = lenA
    #if lenB < shorter:
    #    shorter = lenB
    #a = 0
    #b = 0
    badA = []
    badB = []
    #print listB
    #print listA
    #print len(listA)
    #print len(listB)
    #counter = 1
    for lumi in listA:
        #print str(counter) + "->" + str(lumi)
        #counter += 1
        if not lumi in listB:
            errors.append("Lumi (" + str(lumi) + ") is in listA but not in listB")
            badB.append(lumi)
            #print "Bad B: " + str(lumi)
    #exit("hola")
    for lumi in listB:
        if not lumi in listA:
            errors.append("Lumi (" + str(lumi) + ") is in listB but not in listA")
            badA.append(lumi)
            #print "Bad A: " + str(lumi)
            
    return badA,badB

########################################################################
def aselectFilesToProcess(listOfFilesToProcess,newRunList):
    selectedFiles = []
    runsToProcess = {}
    processedRuns = {}
    for file in listOfFilesToProcess:
        run = getRunNumberFromDBSName(file)
#        print "To process: " + str(run) 
        if run not in runsToProcess:
            runsToProcess[run] = 1
        else:
            runsToProcess[run] = runsToProcess[run] + 1 

    for file in newRunList:
        run = getRunNumberFromFileName(file)
#        print "Processed: " + str(run)
        if run not in processedRuns:
            processedRuns[run] = 1
        else:
            processedRuns[run] = processedRuns[run] + 1 

    #WARNING: getLastClosedRun MUST also have a timeout otherwise the last run will not be considered
    lastClosedRun = getLastClosedRun(listOfFilesToProcess)
#    print "LastClosedRun:-" + str(lastClosedRun) + "-"

    processedRunsKeys = processedRuns.keys()
    processedRunsKeys.sort()

    for run in processedRunsKeys:
        if run <= lastClosedRun :
            print "For run " + str(run) + " I have processed " + str(processedRuns[run]) + " files and in DBS there are " + str(runsToProcess[run]) + " files!"
            if not run in runsToProcess:
                exit("ERROR: I have a result file for run " + str(run) + " but it doesn't exist in DBS. Impossible but it happened!")
            lumiList = getDBSLumiListForRun(run)
            if processedRuns[run] == runsToProcess[run]:
                for file in newRunList:
                    if run == getRunNumberFromFileName(file):
                        selectedFiles.append(file)
            else:
                exit("ERROR: For run " + str(run) + " I have processed " + str(processedRuns[run]) + " files but in DBS there are " + str(runsToProcess[run]) + " files!")
    return selectedFiles            

########################################################################
def sendEmail(mailList,error):
    print "Sending email to " + mailList + " with body: " + error


########################################################################
def main():
    ######### COMMAND LINE OPTIONS ##############
    option,args = parse(__doc__)

    #Right now always in the test DB
    destDB = 'oracle://cms_orcon_prod/CMS_COND_31X_BEAMSPOT'
    #destDB = 'oracle://cms_orcoff_prep/CMS_COND_BEAMSPOT'
    if option.Test:
        destDB = 'oracle://cms_orcoff_prep/CMS_COND_BEAMSPOT'

    ######### CONFIGURATION FILE ################
    configurationFile = os.getenv("CMSSW_BASE") + "/src/RecoVertex/BeamSpotProducer/scripts/BeamSpotWorkflow.cfg"
    configuration     = ConfigParser.ConfigParser()
    print 'Reading configuration from ', configurationFile
    configuration.read(configurationFile)

    sourceDir   = configuration.get('Common','SOURCE_DIR')
    archiveDir  = configuration.get('Common','ARCHIVE_DIR')
    workingDir  = configuration.get('Common','WORKING_DIR')
    databaseTag = configuration.get('Common','DBTAG')
    dataSet     = configuration.get('Common','DATASET')
    fileIOVBase = configuration.get('Common','FILE_IOV_BASE')
    dbIOVBase   = configuration.get('Common','DB_IOV_BASE')
    dbsTolerance= float(configuration.get('Common','DBS_TOLERANCE'))
    mailList    = configuration.get('Common','EMAIL')

    ######### DIRECTORIES SETUP #################
    if sourceDir[len(sourceDir)-1] != '/':
        sourceDir = sourceDir + '/'
    if not dirExists(sourceDir):
        error = "ERROR: The source directory " + sourceDir + " doesn't exist!"
        sendEmail(mailList,error)
        exit(error)

    if archiveDir[len(archiveDir)-1] != '/':
        archiveDir = archiveDir + '/'
    if not os.path.isdir(archiveDir):
	os.mkdir(archiveDir)

    if workingDir[len(workingDir)-1] != '/':
        workingDir = workingDir + '/'
    if not os.path.isdir(workingDir):
	os.mkdir(workingDir)
    else:
        os.system("rm -f "+ workingDir + "*") 


    print "Getting last IOV for tag: " + databaseTag
    lastUploadedIOV = 1
    if destDB == "oracle://cms_orcon_prod/CMS_COND_31X_BEAMSPOT": 
        lastUploadedIOV = getLastUploadedIOV(databaseTag)
    else:
        lastUploadedIOV = getLastUploadedIOV(databaseTag,destDB)
        
    #lastUploadedIOV = 133885
    #lastUploadedIOV = 575216380019329
    if dbIOVBase == "lumiid":
        lastUploadedIOV = unpackLumiid(lastUploadedIOV)["run"]

    ######### Get list of files processed after the last IOV  
    print "Getting list of files processed after IOV " + str(lastUploadedIOV)
    newProcessedRunList      = getNewRunList(sourceDir,lastUploadedIOV)
    if len(newProcessedRunList) == 0:
        exit("There are no new runs after " + str(lastUploadedIOV))

    ######### Copy files to archive directory
    print "Copying files to archive directory"
    copiedFiles = []
    for i in range(3):
        copiedFiles = cp(sourceDir,archiveDir,newProcessedRunList)    
        if len(copiedFiles) == len(newProcessedRunList):
            break;
    if len(copiedFiles) != len(newProcessedRunList):
        error = "ERROR: I can't copy more than " + str(len(copiedFiles)) + " files out of " + str(len(newProcessedRunList)) 
        sendEmail(mailList,error)
        exit(error)


    ######### Get from DBS the list of files after last IOV    
    #listOfFilesToProcess = getListOfFilesToProcess(dataSet,lastUploadedIOV) 
    print "Getting list of files from DBS"
    listOfRunsAndLumiFromDBS = getListOfRunsAndLumiFromDBS(dataSet,lastUploadedIOV) 
    print "Getting list of files from RR"
    listOfRunsAndLumiFromRR  = getListOfRunsAndLumiFromRR(dataSet,lastUploadedIOV) 

    ######### Get list of files to process for DB
#    selectedFilesToProcess = selectFilesToProcess(listOfFilesToProcess,copiedFiles)
    print "Getting list of files to process"
    selectedFilesToProcess = selectFilesToProcess(listOfRunsAndLumiFromDBS,listOfRunsAndLumiFromRR,copiedFiles,archiveDir,dbsTolerance)
    if len(selectedFilesToProcess) == 0:
       exit("There are no files to process")
        
    #print selectedFilesToProcess
    ######### Copy files to working directory
    print "Copying files from archive to working directory"
    copiedFiles = []
    for i in range(3):
        copiedFiles = cp(archiveDir,workingDir,selectedFilesToProcess)    
        if len(copiedFiles) == len(selectedFilesToProcess):
            break;
    if len(copiedFiles) != len(selectedFilesToProcess):
        error = "ERROR: I can't copy more than " + str(len(copiedFiles)) + " files out of " + str(len(selectedFilesToProcess)) 
        sendEmail(mailList,error)
        exit(error)

    print "Sorting and cleaning beamlist"
    beamSpotObjList = []
    for fileName in copiedFiles:
        readBeamSpotFile(workingDir+fileName,beamSpotObjList,fileIOVBase)

    sortAndCleanBeamList(beamSpotObjList,fileIOVBase)

    if len(beamSpotObjList) == 0:
        error = "WARNING: None of the processed and copied payloads has a valid fit so there are no results. This shouldn't happen since we are filtering using the run register, so there should be at least one good run."
        exit(error)

    payloadFileName = "PayloadFile.txt"

    runBased = False
    if dbIOVBase == "runnumber":
        runBased = True
        
    payloadList = createWeightedPayloads(workingDir+payloadFileName,beamSpotObjList,runBased)
    if len(payloadList) == 0:
        error = "WARNING: I wasn't able to create any payload even if I have some BeamSpot objects."
        exit(error)
       

    tmpPayloadFileName = workingDir + "SingleTmpPayloadFile.txt"
    tmpSqliteFileName  = workingDir + "SingleTmpSqliteFile.db"

    writeDBTemplate = os.getenv("CMSSW_BASE") + "/src/RecoVertex/BeamSpotProducer/test/write2DB_template.py"
    readDBTemplate  = os.getenv("CMSSW_BASE") + "/src/RecoVertex/BeamSpotProducer/test/readDB_template.py"
    payloadNumber = -1
    iovSinceFirst = 0;
    iovTillLast   = 0;

    #Creating the final name for the combined sqlite file
    uuid = commands.getstatusoutput('uuidgen -t')[1]
    final_sqlite_file_name = databaseTag + '@' + uuid
    sqlite_file     = workingDir + final_sqlite_file_name + ".db"
    metadata_file   = workingDir + final_sqlite_file_name + ".txt"

    for payload in payloadList:
        payloadNumber += 1
        if option.zlarge:
            payload.sigmaZ = 10
            payload.sigmaZerr = 2.5e-05
        tmpFile = file(tmpPayloadFileName,'w')
        dumpValues(payload,tmpFile)
        tmpFile.close()
        if not writeSqliteFile(tmpSqliteFileName,databaseTag,dbIOVBase,tmpPayloadFileName,writeDBTemplate,workingDir):
            error = "An error occurred while writing the sqlite file: " + tmpSqliteFileName
            exit(error)
        readSqliteFile(tmpSqliteFileName,databaseTag,readDBTemplate,workingDir)
        
        ##############################################################
        #WARNING I am not sure if I am packing the right values
        if dbIOVBase == "runnumber":
            iov_since = payload.Run
            iov_till  = iov_since
        elif dbIOVBase == "lumiid":
	    iov_since = str( pack(int(payload.Run), int(payload.IOVfirst)) )
            iov_till  = str( pack(int(payload.Run), int(payload.IOVlast)) )
        elif dbIOVBase == "timestamp":
            error = "ERROR: IOV " + dbIOVBase + " still not implemented."
            exit(error)
        else:
            error = "ERROR: IOV " + dbIOVBase + " unrecognized!"
            exit(error)

        if payloadNumber == 0:
            iovSinceFirst = iov_since
        elif payloadNumber == len(payloadList)-1:
            iovTillLast   = iov_till
            
        appendSqliteFile(final_sqlite_file_name + ".db", tmpSqliteFileName, databaseTag, iov_since, iov_till ,workingDir)
        os.system("rm -f " + tmpPayloadFileName + " " + tmpSqliteFileName)

        
    #### CREATE payload for merged output

    print " create MERGED payload card for dropbox ..."

    dfile = open(metadata_file,'w')

    dfile.write('destDB '  + destDB        +'\n')
    dfile.write('tag '     + databaseTag   +'\n')
    dfile.write('inputtag'                 +'\n')
    dfile.write('since '   + iovSinceFirst +'\n')
    #dfile.write('till '    + iov_till      +'\n')
    dfile.write('Timetype '+ dbIOVBase     +'\n')

    ###################################################
    # WARNING tagType forced to offline
    print "WARNING TAG TYPE forced to be just offline"
    tagType = "offline"
    checkType = tagType
    if tagType == "express":
        checkType = "hlt"
    dfile.write('IOVCheck ' + checkType + '\n')
    dfile.write('usertext Beam spot position\n')
            
    dfile.close()

                                                                                                

    if option.upload:
        print " scp files to offline Drop Box"
        dropbox = "/DropBox"
        if option.Test:
            dropbox = "/DropBox_test"
        print "UPLOADING TO TEST DB"
        uploadSqliteFile(workingDir, final_sqlite_file_name, dropbox)
                   
    archive_sqlite_file_name = "Payloads_" + iovSinceFirst + "_" + iovTillLast + "_" + final_sqlite_file_name
    if not os.path.isdir(archiveDir + 'payloads'):
        os.mkdir(archiveDir + 'payloads')
    commands.getstatusoutput('mv ' + sqlite_file   + ' ' + archiveDir + 'payloads/' + archive_sqlite_file_name + '.db')
    commands.getstatusoutput('mv ' + metadata_file + ' ' + archiveDir + 'payloads/' + archive_sqlite_file_name + '.txt')
  
    print archiveDir + "payloads/" + archive_sqlite_file_name + '.db'
    print archiveDir + "payloads/" + archive_sqlite_file_name + '.txt'
   

if __name__ == '__main__':
    main()
