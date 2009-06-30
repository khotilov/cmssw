#! /bin/bash

current_area=`pwd`
echo $current_area
# Define the directory that will hold the histogram root files for Full Simulation
# Note: Both Full Sim and Fast Sim will produce histogram root files with the same name, e.g METTester_data_QCD_30-50.root, so they need to be output to different directories!!!

FullSimRootFileDirectory=${current_area}/FullSim/
mkdir $FullSimRootFileDirectory -p

#======= Define list of samples that you will be validating ========#
dirlist="QCD_Pt_80_120 QCD_Pt_3000_3500 Wjet_Pt_80_120 LM1_sfts TTbar QCD_FlatPt_15_3000"

#======= Define list of modules that will be run for each sample ========#
RunPath="fileSaver, calotoweroptmaker, analyzeRecHits, analyzecaloTowers, analyzeGenMET, analyzeGenMETFromGenJets, analyzeHTMET, analyzeCaloMET, analyzeTCMET,OB analyzePFMET"


echo "Run path = {" $RunPath "}"
cmssw_version="3_1_0_pre11"
condition="MC_31X_V1-v1"

#==========================================#
cd $current_area

for i in $dirlist; do


#==========================================#
cd $current_area

#======Make RunAnalyzers_cfg.py=================#
echo "import FWCore.ParameterSet.Config as cms

process = cms.Process(\"TEST\")
process.load(\"RecoMET.Configuration.CaloTowersOptForMET_cff\")

process.load(\"RecoMET.Configuration.RecoMET_cff\")

process.load(\"RecoMET.Configuration.RecoHTMET_cff\")

process.load(\"RecoMET.Configuration.RecoGenMET_cff\")

process.load(\"RecoMET.Configuration.GenMETParticles_cff\")

process.load(\"RecoMET.Configuration.RecoPFMET_cff\")

process.load(\"RecoJets.Configuration.CaloTowersRec_cff\")

process.load(\"Validation.RecoMET.CaloMET_cff\")

process.load(\"Validation.RecoMET.GenMET_cff\")

process.load(\"Validation.RecoMET.HTMET_cff\")

process.load(\"Validation.RecoMET.GenMETFromGenJets_cff\")

process.load(\"DQMOffline.JetMET.caloTowers_cff\")

process.load(\"DQMOffline.JetMET.RecHits_cff\")

process.load(\"Validation.RecoMET.PFMET_cff\")

process.load(\"Validation.RecoMET.TCMET_cff\")

process.load(\"Validation.RecoMET.MuonCorrectedCaloMET_cff\")

process.load(\"Configuration.StandardSequences.Geometry_cff\")

process.load(\"Configuration.StandardSequences.MagneticField_cff\")

process.load(\"Configuration.StandardSequences.FrontierConditions_GlobalTag_cff\")

process.GlobalTag.globaltag = cms.string(\"MC_31X_V1::All\")

process.load(\"RecoLocalCalo.Configuration.hcalLocalReco_cff\")

process.DQMStore = cms.Service(\"DQMStore\")

process.source = cms.Source(\"PoolSource\",
    debugFlag = cms.untracked.bool(True),
    debugVebosity = cms.untracked.uint32(10),
    fileNames = cms.untracked.vstring(
"> ${FullSimRootFileDirectory}/RunAnalyzers-${i}_cfg.py


ds=/RelVal$i/CMSSW_$cmssw_version-$condition/GEN-SIM-RECO
./DDSearchCLI.py --verbose=0 --limit=-1 --input="find file where  dataset=$ds" | grep "root" | sed "s/^/'/g" | sed "s/$/',/g">dd
cat dd >> ${FullSimRootFileDirectory}/RunAnalyzers-${i}_cfg.py
rm -f dd
echo "
    )


)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(8000) )


process.fileSaver = cms.EDFilter(\"METFileSaver\",
    OutputFile = cms.untracked.string('METTester_data_${i}.root')
) 
process.p = cms.Path(process.fileSaver*
                     process.calotoweroptmaker*
                     process.analyzeRecHits*
                     process.analyzecaloTowers*
                     process.analyzeGenMET*
                     process.analyzeGenMETFromGenJets*
                     process.analyzeHTMET*
                     process.analyzeCaloMET*
                     process.analyzePFMET*
                     process.analyzeTCMET*
                     process.analyzeMuonCorrectedCaloMET
)
process.schedule = cms.Schedule(process.p)

" >> ${FullSimRootFileDirectory}/RunAnalyzers-${i}_cfg.py
#============================================#
cd $current_area
done
