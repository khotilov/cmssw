import FWCore.ParameterSet.Config as cms

topHLTDiMuonDQM = cms.EDAnalyzer("TopHLTDiMuonDQM",

    ### Select a level: TEV, L1 or L3
    ### 
    Level          = cms.untracked.string('TEV'),
    monitorName    = cms.untracked.string('HLT/Top/HLTDiMuon'),
    ### 
    TriggerResults = cms.InputTag('TriggerResults',        '','HLT'),
    TriggerEvent   = cms.InputTag('hltTriggerSummaryAOD',  '','HLT'),
    TriggerFilter  = cms.InputTag('hltL1MuOpenL1Filtered0','','HLT'),
    ### 
    hltPaths_L1    = cms.vstring('HLT_L1MuOpen','HLT_L1Mu','HLT_L1Mu20','HLT_L1DoubleMuOpen'),
    hltPaths_L3    = cms.vstring('HLT_Mu3','HLT_IsoMu3','HLT_Mu5','HLT_Mu9','HLT_IsoMu9','HLT_Mu15','HLT_DoubleMu0','HLT_DoubleMu3'),
    hltPaths_sig   = cms.vstring('HLT_IsoMu9', 'HLT_Mu15', 'HLT_DoubleMu3', 'HLT_Mu9'),
    hltPaths_trig  = cms.vstring('HLT_Mu9',    'HLT_Mu9',  'HLT_Mu9',       'HLT_Mu5'),
    ### 
    L1_Collection  = cms.untracked.InputTag('hltL1extraParticles'),
    L3_Collection  = cms.untracked.InputTag('hltL3MuonCandidates'),
    L3_Isolation   = cms.untracked.InputTag('hltL3MuonIsolations'),
    ### 
    vertex_X_cut     = cms.double(  1.0 ),
    vertex_Y_cut     = cms.double(  1.0 ),
    vertex_Z_cut     = cms.double( 20.0 ),
    ### 
    muon_pT_cut    = cms.double( 1.0 ),
    muon_eta_cut   = cms.double( 5.0 ),
    ### 
    MassWindow_up   = cms.double( 1000. ),
    MassWindow_down = cms.double(    1. )

)

topHLTDiMuonAnalyzer = cms.Sequence(topHLTDiMuonDQM)
