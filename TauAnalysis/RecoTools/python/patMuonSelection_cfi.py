import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------  
# produce collections of pat::Muons passing selection criteria
#
# NOTE: the final cut values are (re)defined in
#
#         TauAnalysis/RecoTools/python/patLeptonSelection_cff.py
#
#--------------------------------------------------------------------------------

# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string parser

# require muon candidate to be a global muon
# (track in muon system linked to track in Pixel + SiTracker detectors)
selectedPatMuonsGlobal = cms.EDFilter("PATMuonSelector",
    cut = cms.string('isGlobalMuon()'),
    filter = cms.bool(False)
)

# require muon candidate to be within geometric acceptance of muon trigger
selectedPatMuonsEta21 = cms.EDFilter("PATMuonSelector",
    cut = cms.string('abs(eta) < 2.1'),
    filter = cms.bool(False)
)

# require muon candidate to have transverse momentum above threshold
selectedPatMuonsPt10 = cms.EDFilter("PATMuonSelector",
    cut = cms.string('pt > 10.'),
    filter = cms.bool(False)
)

# require muon candidate to pass VBTF selection
# (selection criteria defined by Vector Boson Task Force
#  and documented in CMS AN-10-264)
selectedPatMuonsVbTfId = cms.EDFilter("PATMuonVbTfSelector",
    beamSpotSource = cms.InputTag("offlineBeamSpot")
)                                      

# require muon candidate to be isolated
# with respect to tracks/charged hadrons
selectedPatMuonsTrkIso = cms.EDFilter("PATMuonIsoDepositSelector",
    type = cms.string('pfChargedHadrons'),
    vetos = cms.vstring("0.01"),                          
    dRisoCone = cms.double(0.4),
    sumPtMax = cms.double(0.10),
    sumPtMethod = cms.string("relative"),                                 
    filter = cms.bool(False)
)

# require muon candidate to be isolated
# with respect to ECAL energy deposits/photons
selectedPatMuonsEcalIso = cms.EDFilter("PATMuonIsoDepositSelector",
    type = cms.string('pfPhotons'),
    vetos = cms.vstring("0.01"),                          
    dRisoCone = cms.double(0.4),
    sumPtMax = cms.double(0.10),
    sumPtMethod = cms.string("relative"),                                 
    filter = cms.bool(False)
)

# require muon candidate to be isolated
# with respect to sum of tracks + ECAL and HCAL energy deposits/
# charged and neutral hadrons + photons
selectedPatMuonsCombIso = cms.EDFilter("PATMuonIsoDepositSelector",
    type = cms.string('pfAllParticles'),
    vetos = cms.vstring("0.01"),                          
    dRisoCone = cms.double(0.4),
    sumPtMax = cms.double(0.10),
    sumPtMethod = cms.string("relative"),                                 
    filter = cms.bool(False)
)

# require muon candidate to pass pion veto
selectedPatMuonsPionVeto = cms.EDFilter("PATMuonAntiPionSelector",
    CaloCompCoefficient = cms.double(0.8),
    SegmCompCoefficient = cms.double(1.2),
    AntiPionCut = cms.double(1.0),
    filter = cms.bool(False)
)

# require muon candidate to be linked to track in silicon strip + pixel detectors
# (all global muons should be linked to tracks in the "inner" tracking detectors;
#  in case the muon is not linked to an "inner" track,
#  the track impact parameter selection will cause processing of the entire event to be skipped !!)
selectedPatMuonsTrk = cms.EDFilter("PATMuonSelector",
    cut = cms.string('innerTrack.isNonnull'),
    filter = cms.bool(False)
)

# require track of muon candidate to have small transverse impact parameter
# (in order to veto muons resulting from b-quark decays)
selectedPatMuonsTrkIP = cms.EDFilter("PATMuonIpSelector",
    vertexSource = cms.InputTag("selectedPrimaryVertexPosition"),
    IpMax = cms.double(0.05),
    filter = cms.bool(False)                                               
)

#--------------------------------------------------------------------------------
# define additional collections of muon candidates
# with loose track and ECAL isolation applied
#
# (NOTE: to be used for the purpose of factorizing efficiencies
#        of muon isolation from other event selection criteria,
#        in order to avoid problems with limited Monte Carlo statistics)
#--------------------------------------------------------------------------------

selectedPatMuonsTrkIsoLooseIsolation = selectedPatMuonsTrkIso.clone()
selectedPatMuonsTrkIsoLooseIsolation.sumPtMax = cms.double(0.25)

selectedPatMuonsEcalIsoLooseIsolation = selectedPatMuonsEcalIso.clone()
selectedPatMuonsEcalIsoLooseIsolation.sumPtMax = cms.double(0.25)

selectedPatMuonsCombIsoLooseIsolation = selectedPatMuonsCombIso.clone()
selectedPatMuonsCombIsoLooseIsolation.sumPtMax = cms.double(0.25)

selectedPatMuonsPionVetoLooseIsolation = selectedPatMuonsPionVeto.clone()

selectedPatMuonsTrkLooseIsolation = selectedPatMuonsTrk.clone()

selectedPatMuonsTrkIPlooseIsolation = selectedPatMuonsTrkIP.clone()
