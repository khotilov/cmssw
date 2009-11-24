#!/bin/csh
setenv sample ${1}

echo '===> sample.'  $sample 

if ( $sample == SingleGammaPt10 ) then
setenv outFileName SingleGammaPt10
else if (  $sample == SingleGammaPt35 ) then
setenv outFileName SingleGammaPt35
else if (  $sample ==  H130GGgluonfusion ) then
setenv outFileName H130GGgluonfusion
else if (  $sample == QCD_Pt_80_120 ) then
setenv outFileName  QCD_Pt_80_120
endif

setenv confName  PhotonValidator

setenv MYWORKDIR /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/slc5_ia32_gcc434/CMSSW_3_4_0_pre6/src/Validation/RecoEgamma/test

echo ${MYWORKDIR}

setenv MYOUT ${MYWORKDIR}
#----------------
cd ${MYWORKDIR}
eval `scramv1 runtime -csh`
cp ${MYWORKDIR}/${confName}_${sample}.py    ${WORKDIR}/conf.py


#
cd ${WORKDIR}
echo ${WORKDIR}

cmsRun  conf.py > & ${outFileName}.log
#---------------------------------------------------------------
 rfcp   ${outFileName}.log             ${MYOUT}/.
 rfcp   PhotonValidationRelVal340pre6_${outFileName}.root            ${MYOUT}/.
