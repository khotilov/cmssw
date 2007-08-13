#!/bin/bash
for input in $*
  do
    if [ -z `echo $input | grep -o '^/store'` ]
    then full_file='file:'$input
    else full_file=$input
    fi
    if [ -n "$input_file" ] 
    then input_file=$input_file", '"$full_file"'"
    else input_file="'"$full_file"'"
    fi
done
echo input_file: $input_file

for flavor in midPointCone5CaloJets midPointCone7CaloJets iterativeCone5CaloJets iterativeCone7CaloJets Fastjet10CaloJets Fastjet6CaloJets 
  do
    echo Processing flavor $flavor ...
    cfg_file=$flavor.cfg
    output_file=$flavor.root
    log_file=$flavor.log
    rm -f $cfg_file $log_file
    cat template.cfg | sed s^INPUT_FILES^"$input_file"^g | sed s^OUTPUT_FILE^$output_file^g | sed s^SOURCE^$flavor^g > $cfg_file
      
    eval `scramv1 ru -sh`
    cmsRun $cfg_file > $log_file 2>&1 
done
