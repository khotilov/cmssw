import FWCore.ParameterSet.Config as cms

from SUSYBSMAnalysis.CSA07Skims.lepSUSYHLTPaths_cfi import *
lepSUSY_2Muon_0Elec_2Jets_MET_Filter = cms.EDFilter("LepSUSYSkim",
    CaloMETsrc = cms.InputTag("met"),
    Elecsrc = cms.InputTag("pixelMatchGsfElectrons"),
    CaloJet1Ptmin = cms.double(80.0),
    MuonPtmin = cms.double(3.0),
    CaloJetsrc = cms.InputTag("iterativeCone5CaloJets"),
    NminElec = cms.int32(0),
    CaloMetmin = cms.double(50.0),
    ElecPtmin = cms.double(5.0),
    NminCaloJet = cms.int32(2),
    Muonsrc = cms.InputTag("muons"),
    NminMuon = cms.int32(2),
    CaloJet2Ptmin = cms.double(30.0)
)

lepSUSY_2Muon_0Elec_2Jets_MET_Seq = cms.Sequence(lepSUSY_2Muon_0Elec_2Jets_MET_HLTPath+lepSUSY_2Muon_0Elec_2Jets_MET_Filter)

