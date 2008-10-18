#!/bin/bash

#if [ "X$CMSSW_BASE" == "X" ]    
#    then    
#    echo "Setup CMSSW environment using eval \`scramv1 runtime -(c)sh\`"    
#    exit    
#fi
#usage='Usage: -r <run number> -d <primary dataset> -s <num events per job> -o <output_dir>'
usage='Usage: -r <run number> -d <primary dataset> -s <num events per job> -at <analysis type>'
#
#args=`getopt rdso: -- "$@"`
args=`getopt rds: -- "$@"`
if test $? != 0
     then
         echo $usage
         exit 1
fi

eval set -- "$args"
for i 
  do
  case "$i" in
      -r) shift; run_num=$2;shift;;
      -d) shift; data_set=$2;shift;;      
      -s) shift; nsplit=$2;shift;;
      -at) shift; analy_type=$2;shift;;
#      -o) shift; output_dir=$2;shift;;
  esac      
done

if [ "X"${run_num} == "X" ]
    then
    echo "INVALID RUN NUMBER! Please give a valid run number!"
    echo $usage
    exit 
fi



if [ "X"${data_set} == "X" ]
    then
    echo "INVALID DATASET! Please give a dataset, eg: "
#	/Cosmics/Commissioning08-CRUZET4_v1/RAW"
    echo " /Cosmics/Commissioning08-MW36_v1/RAW "
    echo " /BeamHalo/Commissioning08-MW36_v1/RAW "
    echo $usage
    exit 
fi
echo 'Submitting CRAB cosmics analysis jobs for' ${data_file} 'dataset' ${data_set}

#if [ "X"${output_dir} == "X" ]
#    then
#    echo " no output dir"
#    output_dir=/castor/cern.ch/user/c/ccecal/CRUZET3/CosmicsAnalysis
#    echo " using default output dir" ${output_dir}
#else
#    echo " using output dir "${output_dir}
#    rfmkdir ${output_dir}
#    rfchmod 775 ${output_dir}
#fi


if [ "X"${nsplit} == "X" ]
    then
    nsplit=50000
    echo " using default split of 50000 events per job"
else
    echo " splitting into "${nsplit} "events per job"
fi



if [ "X"${analy_type} == "X" ]
    then
    analy_type="Laser"
    echo " using default analysis type of Laser"
else
    echo " doing analysis on ${analy_type} of events  "
fi



this_dir=`pwd`;

# set up directory
mkdir -p ${analy_type}_${run_num};

# extract time info
#cd /afs/cern.ch/cms/CAF/CMSCOMM/COMM_ECAL/ccecal/usefulScripts/eventTiming;
#pwd
#./getRunTimeAbbrev.csh ${run_num} >& temp.sh
#source temp.sh
#echo " startdate and runlength: "$startdate $runlength
#cd -;

#if [ "X"${startdate} == "X" ]
#then
    startdate=1215107133
    echo "Warning: did not get startdate.. using default: " $startdate
#fi
#
#if [ "X"${runlength} == "X" ]
#then
    runlength=3
    echo "Warning: did not get runlength.. using default: " $runlength
#fi

# set up Laser.cfg

cat ${analy_type}TEMPLATE.cfg | /bin/sed "s@RUNNUMBER@${run_num}@g" > ${analy_type}_${run_num}/${analy_type}Config_${run_num}.cfg

# set up crab.cfg file

cat CRABTEMPLATE.cfg | /bin/sed "s@RUNNUMBER@${run_num}@g; s@EVENTSPERJOB@${nsplit}@g; s@DATASETPLACE@${data_set}@g; s@ANALYTYPE@${analy_type}@g  " > ${analy_type}_${run_num}/crab.cfg




#cat ${run_num}/crab.cfg

cd ${analy_type}_${run_num};
#pwd

echo " setting up crab env"
# setup crab environment
source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh;
eval `scramv1 runtime -sh`;
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh;

echo " launching crab jobs"
crab -create;
sleep 3;
crab -submit;
sleep 5;
crab -status;
cd -; 

#unset startdate
#unset runlength
#rm -f temp
#rm -f temp.sh

exit







