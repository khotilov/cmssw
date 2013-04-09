import FWCore.ParameterSet.Config as cms

from MuonAnalysis.PFMuonTools.particleFlowCleaned_cfi import particleFlowCleaned
from MuonAnalysis.PFMuonTools.muonsCleaned_cfi import muonsCleaned



#Add a PF Muon selector
pfMuons = cms.EDFilter("PdgIdPFCandidateSelector",
                       src = cms.InputTag("particleFlowCleaned"),
                       pdgId = cms.vint32( -13, 13)
)


muonCleaning = cms.Sequence(particleFlowCleaned+pfMuons+muonsCleaned)


