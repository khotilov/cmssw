#!/usr/bin/env python

"""
Create a list of files and sort them by IOV.
Loop on the list and fill the dictionary of values.
"""

import os
import os.path
import time
import calendar


# This function returns a list of the directory content sorted by time (older first)
def listDirSortedByTime(dir):
    # Build up a dictionary of timestamps
    content = {}
    for item in os.listdir(dir):
        content[item] = os.path.getmtime(dir+"/"+item)
    # Sort keys, based on time stamps
    items = content.keys()
    items.sort(lambda x,y: cmp(content[x],content[y]))
    # Report objects in order
    for i in range(0,len(items)):
        items[i] = dir+"/"+items[i]
    return items

def startFile(outputFile, label):
    outputFile.write("{\n")
    outputFile.write("    label: '"+label+"',\n")
    outputFile.write("    data: [")

def fillFile(outputFile, fileBlock, num):
    outputFile.write("["+str(fileBlock[0]+1)+", "+str(num)+"], ")
    # # Add 1 to force it to appear after the previous. The precision is much smaller, so the last three digits always 0 in any case.
    outputFile.write("["+str(fileBlock[1])+", "+str(num)+"]")


totalTrackerModules = 15148

outputFileHV = open("full_hv.js", "w")
startFile(outputFileHV, 'High Voltage ON')

outputFileLV = open("full_lv.js", "w")
startFile(outputFileLV, 'Low Voltage ON')

first = True
for subDir in os.listdir("./"):
    if os.path.isdir(subDir):
        if subDir.startswith("20"):
            if first == True:
                dirList = listDirSortedByTime(subDir)
                first = False
            else:
                dirList.extend(listDirSortedByTime(subDir))

fileList = []

for inputFileName in dirList:
    if inputFileName.find("DetVOffReaderSummary") != -1 and inputFileName.endswith(".log"):
        firstDateString = inputFileName.split("_TO_")[0].split("_FROM_")[1]
        dateArray = firstDateString.replace("__", "_").split("_")
        firstTimeValue = int(time.mktime(time.strptime(dateArray[0]+" "+dateArray[1]+" "+dateArray[2]+" "+dateArray[3]+":"+dateArray[4]+":"+dateArray[5]+" "+dateArray[6])))
        firstTimeValue = firstTimeValue*1000
        
        dateString = inputFileName.split("_TO_")[1].split(".")[0]
        dateArray = dateString.replace("__", "_").split("_")
        lastTimeValue = int(time.mktime(time.strptime(dateArray[0]+" "+dateArray[1]+" "+dateArray[2]+" "+dateArray[3]+":"+dateArray[4]+":"+dateArray[5]+" "+dateArray[6])))
        lastTimeValue = lastTimeValue*1000
        
        fileList.append((firstTimeValue, lastTimeValue, inputFileName))

# Loop on the file list sorted with the time of the start of the IOV
first = True

# Note: sorted is guaranteed to be stable since python 2.2
for fileBlock in sorted(fileList, key=lambda fileTuple: fileTuple[0]):
    inputFile = open(fileBlock[2], "r")
    totHVoff = 0
    totLVoff = 0
    counter = 0
    checkLine = False
    for line in inputFile:
        if "TIB" in line:
            checkLine = True
        if "Summary" in line and checkLine:
            totHVoff = counter
            counter = 0
            checkLine = False
        if "%MSG" in line and checkLine:
            totLVoff = counter
            counter = 0
            checkLine = False
        if checkLine and len(line) != 1:
            counter += int(line.rsplit(" ", 1)[1])

    # print "Total modules with HV off =", totHVoff
    # print "Total modules with LV off =", totLVoff
    if first:
        first = False
    else:
        outputFileHV.write(", ")
        outputFileLV.write(", ")

    fillFile(outputFileHV, fileBlock, totalTrackerModules - totHVoff)
    fillFile(outputFileLV, fileBlock, totalTrackerModules - totLVoff)
    
outputFileHV.write("]\n}")
outputFileLV.write("]\n}")
