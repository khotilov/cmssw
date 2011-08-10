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
        items[i] = dir+"/"+item
        # print "items[",i,"] = ", items[i]
    return items


totalTrackerModules = 15148

outputFile = open("full.js", "w")

outputFile.write("{\n")
outputFile.write("    label: 'High Voltage ON',\n")
outputFile.write("    data: [")

# dirList = os.listdir("./")
first = True
for subDir in os.listdir("./"):
    if os.path.isdir(subDir):
        if subDir.startswith("2"):
            if first == True:
                # print "First in", subDir
                dirList = listDirSortedByTime(subDir)
                first = False
            else:
                # print "Second in", subDir
                dirList.extend(listDirSortedByTime(subDir))

# print dirList

fileList = []

for inputFileName in dirList:
    if inputFileName.find("DetVOffReaderSummary") != -1 and inputFileName.endswith(".log"):
        firstDateString = inputFileName.split("_TO_")[0].split("_FROM_")[1]
        dateArray = firstDateString.replace("__", "_").split("_")
        firstTimeValue = calendar.timegm(time.strptime(dateArray[0]+" "+dateArray[1]+" "+dateArray[2]+" "+dateArray[3]+":"+dateArray[4]+":"+dateArray[5]+" "+dateArray[6]))
        firstTimeValue = firstTimeValue*1000
        
        dateString = inputFileName.split("_TO_")[1].split(".")[0]
        dateArray = dateString.replace("__", "_").split("_")
        lastTimeValue = calendar.timegm(time.strptime(dateArray[0]+" "+dateArray[1]+" "+dateArray[2]+" "+dateArray[3]+":"+dateArray[4]+":"+dateArray[5]+" "+dateArray[6]))
        lastTimeValue = lastTimeValue*1000
        
        fileList.append((firstTimeValue, lastTimeValue, inputFileName))

# Loop on the file list sorted with the time of the start of the IOV
first = True
#firstTimeValue = 0
#for inputFileName in dirList:
#    if inputFileName.endswith(".log"):

# Note: sorted is guaranteed to be stable since python 2.2
for fileBlock in sorted(fileList, key=lambda fileTuple: fileTuple[0]):
    inputFile = open(fileBlock[2], "r")
    totHVoff = 0
    totLVoff = 0
    counter = 0
    checkLine = False
    for line in inputFile:
        # print "line", line
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
        outputFile.write(", ")
    outputFile.write("["+str(fileBlock[0]+1)+", "+str(totalTrackerModules - totHVoff)+"], ")
    # # Add 1 to force it to appear after the previous. The precision is much smaller, so the last three digits always 0 in any case.
    outputFile.write("["+str(fileBlock[1])+", "+str(totalTrackerModules - totHVoff)+"]")

#    if first:
#        print "firstTimeValue =", firstTimeValue
#        valuesDict[firstTimeValue] = totHVoff
#        first = False
#    print "timeValue =", timeValue
#    valuesDict[timeValue] = totHVoff
    # outputFile.write(str(totLVoff))

#    data: [[0, 30], [1, 39], [2, 20], [3, 12], [4, 5], [6, 6], [7, 20], [8, 31], [9, 29], [10, 9]]
# outputFile.write("]\n}")

# print sorted(valuesDict)

#first = True
#for key in sorted(valuesDict):
#    if not first:
#        outputFile.write(" ,")
#    else:
#        # outputFile.write("["+str(firstTimeValue)+", "+str(totHVoff)+"], ")
#        # valuesDict[firstTimeValue] = totHVoff
#        first = False
#    # print "key =", key, "value =", valuesDict[key]
#    outputFile.write("["+str(key)+", "+str(valuesDict[key])+"]")
    
outputFile.write("]\n}")
