#!/bin/bash

# Author: Seth Cooper (UMN), August 4 2011
# Dump the XMLs for each IOV given for a specific given tag from the database
#   and diff them to local copies of the same files.

usage(){
  echo
  echo "Usage: $0 tagName iovsToCheck.txt XMLFiles.txt"
  echo "    where tagName is the name of the tag in the DB"
  echo "          iovsToCheck is a file listing the (since) IOVs, one per line"
  echo "          XMLFiles is a file listing the XML files to be compared in the same order as the IOVs, one per line"
  echo
  echo "Example: $0 EcalTimeCalibConstants_v06_offline iovsToCheck.txt xmlsToCheck.txt"
  exit 1
}

myTagName=$1
iovFile=$2
xmlFile=$3

# validate user input
if [ -z "$myTagName" ] && [ "${myTagName+xxx}" = "xxx" ]; then
  echo "ERROR: no tag name given"
  usage
fi

if [ -z "$iovFile" ] && [ "${iovFile+xxx}" = "xxx" ]; then
  echo "ERROR: no file of IOVs given"
  usage
fi

if [ -z "$xmlFile" ] && [ "${xmlFile+xxx}" = "xxx" ]; then
  echo "ERROR: no file of XMLs given"
  usage
fi

myCondDBTool=$CMSSW_RELEASE_BASE/python/CondTools/Ecal/EcalCondDB.py
myIOVs=( $( cat $iovFile ) )
myXMLs=( $( cat $xmlFile ) )


# was for testing
#echo "a[*] = \"${a[*]}\""
#for i in $(seq 0 $((${#a[*]} - 1))); do
#    echo "a[$i] = \"${a[$i]}\""
#  done

echo "Using tag: $myTagName"

# Loop over IOVs for this tag and dump to XMLs

for i in $(seq 0 $((${#myIOVs[@]} - 1))); do
  myFileName=$myTagName${myIOVs[$i]}.xml
  echo
  echo "Dumping IOV: ${myIOVs[$i]} from tag: $myTagName into file: $myFileName"
  $myCondDBTool --dump=$myFileName -t $myTagName -s ${myIOVs[$i]}
  # Diff to local XML
  echo "Diffing $myFileName to ${myXMLs[$i]}"
  diff $myFileName ${myXMLs[$i]}
done

exit 0

