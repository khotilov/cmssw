import FWCore.ParameterSet.Config as cms

# import config for event selection, event print-out and analysis sequence
from TauAnalysis.Configuration.analyzeZtoElecTau_cfi import *

analyzeZtoElecTau = cms.EDAnalyzer("GenericAnalyzer",
  
  name = cms.string('zElecTauAnalyzer'), 
                            
  eventSelection = cms.VPSet(
    # generator level phase-space selection
    # (NOTE: (1) to be used in case of Monte Carlo samples
    #            overlapping in simulated phase-space only !!
    #        (2) genPhaseSpaceCut needs to be **always** the first entry in the list of cuts
    #           - otherwise the script submitToBatch.csh for submission of cmsRun jobs
    #            to the CERN batch system will not work !!)
    genPhaseSpaceCut,
    
    # generator level selection of Z --> e + tau-jet events
    # passing basic acceptance and kinematic cuts
    # (NOTE: to be used for efficiency studies only !!)
    #genElectronCut,
    #genTauCut,
    
    # trigger selection
    #Trigger,

    # primary event vertex selection
    primaryEventVertex,
    primaryEventVertexQuality,
    primaryEventVertexPosition,
            
    # electron candidate selection
    tightElectronIdCut,
    electronAntiCrackCut,
    electronEtaCut,
    electronPtCut,
    electronTrkIsoCut,
    electronEcalIsoCut,
    electronTrkCut,
    electronTrkIPcut,

    # tau candidate selection
    tauAntiOverlapWithElectronsVeto,
    tauEtaCut,
    tauPtCut,
    tauLeadTrkCut,
    tauLeadTrkPtCut,
    tauTrkIsoCut,
    tauEcalIsoCut,
    tauProngCut,
    tauElectronVeto,

    # di-tau candidate selection
    diTauCandidateForElecTauAntiOverlapVeto,
    diTauCandidateForElecTauZeroChargeCut,
    diTauCandidateForElecTauMt1METCut,

    # veto events containing additional central jets with Et > 20 GeV
    #centralJetVeto
  ),
  
  histManagers = cms.VPSet(
    genPhaseSpaceEventInfoHistManager,
    electronHistManager,
    tauHistManager,
    diTauCandidateHistManagerForElecTau,
    metHistManager,
    vertexHistManager,
    triggerHistManager
  ),

  eventDumps = cms.VPSet(
    elecTauEventDump
  ),
   
  analysisSequence = elecTauAnalysisSequence
)
