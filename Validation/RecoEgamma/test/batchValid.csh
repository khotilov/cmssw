#!/bin/csh
setenv sample ${1}

echo '===> sample.'  $sample 

if ( $sample == SingleGammaPt10 ) then
setenv outFileName SingleGammaPt10
else if (  $sample == SingleGammaPt35 ) then
setenv outFileName SingleGammaPt35
else if (  $sample == SingleGammaFlatPt10To100 ) then
setenv outFileName SingleGammaFlatPt10To100
else if (  $sample ==  H130GGgluonfusion ) then
setenv outFileName H130GGgluonfusion
else if (  $sample == PhotonJets_Pt_10 ) then
setenv outFileName  PhotonJets_Pt_10
else if (  $sample == QCD_Pt_20_30 ) then
setenv outFileName  QCD_Pt_20_30
else if (  $sample == QCD_Pt_80_120 ) then
setenv outFileName  QCD_Pt_80_120
endif

setenv confName  PhotonValidator
setenv MYWORKDIR /afs/cern.ch/user/n/nancy/scratch0/CMSSW/test/CMSSW_3_9_0_pre2/src/Validation/RecoEgamma/test

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
 rfcp   PhotonValidationRelVal390pre2_${outFileName}.root            ${MYOUT}/.


