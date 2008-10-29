#!/bin/bash

Host=`hostname -i`
Port="54125"
Label=DEBUG

source setup.sh

while getopts ":h:p:l:i:" Option
do
case $Option in
h ) Host=$OPTARG ;;
p ) Port=$OPTARG ;;
l ) Label=$OPTARG ;;
esac
done

echo "Host :" $Host
echo "Port :" $Port
echo "Label :" $Label

#Y="$(ps -A | grep -ic 'xdaq.exe')"

#if [ $Y -gt 0 ]; then
#    echo "$Y copies of xdaq.exe running"
#else
    while :
	do
	  CONFIG_FILE="HLXReadoutXDAQConfig.xml" #`hostname`
	  echo "Run xdaq with options -l ${Label} -h ${Host} -p ${Port} -c ${XDAQ_ROOT}/RecoLuminosity/HLXReadOut/XDAQ/xml/HLXReadoutXDAQConfig.xml"
	  /opt/xdaq/bin/xdaq.exe -l $Label -h $Host -p $Port -c /opt/TriDAS/RecoLuminosity/HLXReadOut/XDAQ/xml/${CONFIG_FILE}
	  if [ "$?" -ne "10" ]; then
	  	  echo Exit "$?"
		  exit "$?"
	  fi
    done
#fi
echo "XDAQ launch_supervisor closed"



