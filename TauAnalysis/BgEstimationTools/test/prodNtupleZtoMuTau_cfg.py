import FWCore.ParameterSet.Config as cms

process = cms.Process("prodNtupleZtoMuTau")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('MC_31X_V2::All')

#--------------------------------------------------------------------------------
# import sequences for PAT-tuple production
process.load("TauAnalysis.Configuration.producePatTuple_cff")
process.load("TauAnalysis.Configuration.producePatTupleZtoMuTauSpecific_cff")
process.load("TauAnalysis.BgEstimationTools.producePatTupleZtoMuTauSpecificForBgEst_cff")

# import sequence for event selection
process.load("TauAnalysis.Configuration.selectZtoMuTau_cff")
process.load("TauAnalysis.BgEstimationTools.bgEstPreselZtoMuTau_cff")

# import configuration parameters for submission of jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_cfi import *
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_10TeV_cfi import *
#from TauAnalysis.BgEstimationTools.bgEstSampleDefinitionsZtoMuTau_7TeV_cfi import *
from TauAnalysis.BgEstimationTools.bgEstSampleDefinitionsZtoMuTau_10TeV_cfi import *

# import event-content definition of products to be stored in patTuple
from TauAnalysis.Configuration.patTupleEventContent_cff import *
from TauAnalysis.Skimming.EventContent_cff import *
#--------------------------------------------------------------------------------

# print event content 
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/A4DD1FAE-B178-DE11-B608-001D09F24EAC.root',
        #'/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/9408B54D-CB78-DE11-9AEB-001D09F2503C.root'
        'rfio:/castor/cern.ch/user/l/lusito/SkimOctober09/ZtautauSkimMT314_3/muTauSkim_1.root',
        'rfio:/castor/cern.ch/user/l/lusito/SkimOctober09/ZtautauSkimMT314_3/muTauSkim_2.root'
    ),
    skipEvents = cms.untracked.uint32(0)            
)

# produce ntuple
kineReweight_fileName = cms.string('rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_1_2/bgEstKineReweights/bgEstKineEventReweightsZtoMuTau.root')
kineReweight_dqmDirectory = "DQMData/bgEstKineEventReweights"
kineReweight_meName = "diTauMvis"
kineVarExtractor_config = cms.PSet(
    pluginType = cms.string("PATMuTauPairValExtractor"),
    src = cms.InputTag('muTauPairsBgEstPreselection'),
    value = cms.string("p4Vis.mass"),
    indices = cms.vuint32(0)
)

process.ntupleProducer = cms.EDAnalyzer("ObjValNtupleProducer",
    treeName = cms.string("bgEstEvents"),
    branches = cms.PSet(
        # variables specific to selection of Z --> mu+ mu- background enriched sample
        muonPtZmumu = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstZmumuEnriched'),
            value = cms.string("leg1.pt"),
            indices = cms.vuint32(0,1)
        ),
        muonAbsEtaZmumu = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstZmumuEnriched'),
            value = cms.string("abs(leg1.eta)"),
            indices = cms.vuint32(0,1)
        ),
        muonTrackIsoZmumu = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstZmumuEnriched'),
            value = cms.string("leg1.trackIso"),
            indices = cms.vuint32(0,1)
        ),
        muonEcalIsoZmumu = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstZmumuEnriched'),
            value = cms.string("leg1.ecalIso"),
            indices = cms.vuint32(0,1)
        ),
        tauDiscrAgainstMuonsZmumu = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstZmumuEnriched'),
            value = cms.string("leg2.tauID('againstMuon')"),
            indices = cms.vuint32(0,1)
        ),
        diTauAbsChargeZmumu = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstZmumuEnriched'),
            value = cms.string("abs(charge)"),
            indices = cms.vuint32(0,1)
        ),
        numDiTausZmumu = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('muTauPairsForBgEstZmumuEnriched')
        ),
        numJetsAlpha0point1Zmumu = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('jetsAlpha0point1ForMuTauBgEstZmumuEnriched')
        ),

        # variables specific to selection of W + jets background enriched sample
        muonPtWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("leg1.pt"),
            indices = cms.vuint32(0,1)
        ),
        muonAbsEtaWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("abs(leg1.eta)"),
            indices = cms.vuint32(0,1)
        ),
        muonTrackIsoWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("leg1.trackIso"),
            indices = cms.vuint32(0,1)
        ),
        muonEcalIsoWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("leg1.ecalIso"),
            indices = cms.vuint32(0,1)
        ),
        tauTrackIsoWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("leg2.chargedHadronIso"),
            indices = cms.vuint32(0,1)
        ),
        tauTrackIsoDiscrWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("leg2.tauID('trackIsolation')"),
            indices = cms.vuint32(0,1)
        ),
        tauEcalIsoWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("leg2.photonIso"),
            indices = cms.vuint32(0,1)
        ),
        tauEcalIsoDiscrWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("leg2.tauID('ecalIsolation')"),
            indices = cms.vuint32(0,1)
        ),
        tauDiscrAgainstMuonsWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("leg2.tauID('againstMuon')"),
            indices = cms.vuint32(0,1)
        ),
        diTauDphi12WplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("dPhi12"),
            indices = cms.vuint32(0,1)
        ),
        diTauMt1MEtWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("mt1MET"),
            indices = cms.vuint32(0,1)
        ),
        diTauPzetaDiffWplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("pZeta - 1.5*pZetaVis"),
            indices = cms.vuint32(0,1)
        ),
        numDiTausWplusJets = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched')
        ), 
        numJetsAlpha0point1WplusJets = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('jetsAlpha0point1ForMuTauBgEstWplusJetsEnriched')
        ),
        
        # variables specific to selection of ttbar + jets background enriched sample
        muonPtTTplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstTTplusJetsEnriched'),
            value = cms.string("leg1.pt"),
            indices = cms.vuint32(0,1)
        ),
        muonAbsEtaTTplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstTTplusJetsEnriched'),
            value = cms.string("abs(leg1.eta)"),
            indices = cms.vuint32(0,1)
        ),
        muonTrackIsoTTplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstTTplusJetsEnriched'),
            value = cms.string("leg1.trackIso"),
            indices = cms.vuint32(0,1)
        ),
        muonEcalIsoTTplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstTTplusJetsEnriched'),
            value = cms.string("leg1.ecalIso"),
            indices = cms.vuint32(0,1)
        ),        
        tauDiscrAgainstMuonsTTplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstTTplusJetsEnriched'),
            value = cms.string("leg2.tauID('againstMuon')"),
            indices = cms.vuint32(0,1)
        ),
        diTauAbsChargeTTplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstTTplusJetsEnriched'),
            value = cms.string("abs(charge)"),
            indices = cms.vuint32(0,1)
        ),
        numDiTausTTplusJets = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('muTauPairsForBgEstTTplusJetsEnriched')
        ), 
        jetEt40bTaggingDiscrSimpleSecondaryVertexTTplusJets = cms.PSet(
            pluginType = cms.string("PATJetValExtractor"),
            src = cms.InputTag('jetsEt40ForMuTauBgEstTTplusJetsEnriched'),
            value = cms.string("bDiscriminator('simpleSecondaryVertex')"),
            indices = cms.vuint32(0,1,2)
        ),
        jetEt40bTaggingDiscrCombinedSecondaryVertexTTplusJets = cms.PSet(
            pluginType = cms.string("PATJetValExtractor"),
            src = cms.InputTag('jetsEt40ForMuTauBgEstTTplusJetsEnriched'),
            value = cms.string("bDiscriminator('combinedSecondaryVertex')"),
            indices = cms.vuint32(0,1,2)
        ),
        jetEt40bTaggingDiscrTrackCountingHighEffTTplusJets = cms.PSet(
            pluginType = cms.string("PATJetValExtractor"),
            src = cms.InputTag('jetsEt40ForMuTauBgEstTTplusJetsEnriched'),
            value = cms.string("bDiscriminator('trackCountingHighEffBJetTags')"),
            indices = cms.vuint32(0,1,2)
        ),
        jetEt40bTaggingDiscrSoftElectronTTplusJets = cms.PSet(
            pluginType = cms.string("PATJetValExtractor"),
            src = cms.InputTag('jetsEt40ForMuTauBgEstTTplusJetsEnriched'),
            value = cms.string("bDiscriminator('softElectronBJetTags')"),
            indices = cms.vuint32(0,1,2)
        ),
        jetEt40bTaggingDiscrSoftMuonTTplusJets = cms.PSet(
            pluginType = cms.string("PATJetValExtractor"),
            src = cms.InputTag('jetsEt40ForMuTauBgEstTTplusJetsEnriched'),
            value = cms.string("bDiscriminator('softMuonBJetTags')"),
            indices = cms.vuint32(0,1,2)
        ),
        numJetsEt40TTplusJets = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('jetsEt40ForMuTauBgEstTTplusJetsEnriched')
        ),
        numJetsEt60TTplusJets = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('jetsEt60ForMuTauBgEstTTplusJetsEnriched')
        ),

        # variables for selection of QCD background enriched sample
        # (loose muon track and ECAL isolation applied)
        muonPtQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("leg1.pt"),
            indices = cms.vuint32(0,1)
        ),
        muonAbsEtaQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("abs(leg1.eta)"),
            indices = cms.vuint32(0,1)
        ),
        muonTrackIsoQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("leg1.trackIso"),
            indices = cms.vuint32(0,1)
        ),
        muonEcalIsoQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("leg1.ecalIso"),
            indices = cms.vuint32(0,1)
        ),
        tauTrackIsoQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("leg2.chargedHadronIso"),
            indices = cms.vuint32(0,1)
        ),
        tauEcalIsoQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("leg2.photonIso"),
            indices = cms.vuint32(0,1)
        ),
        tauDiscrAgainstMuonsQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("leg2.tauID('againstMuon')"),
            indices = cms.vuint32(0,1)
        ),
        diTauAbsChargeQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("abs(charge)"),
            indices = cms.vuint32(0,1)
        ),
        diTauDphi12QCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("dPhi12"),
            indices = cms.vuint32(0,1)
        ),
        diTauMt1MEtQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("mt1MET"),
            indices = cms.vuint32(0,1)
        ),
        diTauPzetaDiffQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("pZeta - 1.5*pZetaVis"),
            indices = cms.vuint32(0,1)
        ),
        numDiTausQCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation')
        ),
        # (no muon track and ECAL isolation applied)
        muonPtQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("leg1.pt"),
            indices = cms.vuint32(0,1)
        ),
        muonAbsEtaQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("abs(leg1.eta)"),
            indices = cms.vuint32(0,1)
        ),
        muonTrackIsoQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("leg1.trackIso"),
            indices = cms.vuint32(0,1)
        ),
        muonEcalIsoQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("leg1.ecalIso"),
            indices = cms.vuint32(0,1)
        ),
        tauTrackIsoQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("leg2.chargedHadronIso"),
            indices = cms.vuint32(0,1)
        ),
        tauEcalIsoQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("leg2.photonIso"),
            indices = cms.vuint32(0,1)
        ),
        tauDiscrAgainstMuonsQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("leg2.tauID('againstMuon')"),
            indices = cms.vuint32(0,1)
        ),
        diTauAbsChargeQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("abs(charge)"),
            indices = cms.vuint32(0,1)
        ),
        diTauDphi12QCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("dPhi12"),
            indices = cms.vuint32(0,1)
        ),
        diTauMt1MEtQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("mt1MET"),
            indices = cms.vuint32(0,1)
        ),
        diTauPzetaDiffQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation'),
            value = cms.string("pZeta - 1.5*pZetaVis"),
            indices = cms.vuint32(0,1)
        ),
        numDiTausQCDnoMuonIso = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoMuonIsolation')
        ),
        # (no muon track and ECAL isolation applied;
        #  no tau leading track requirements, track and ECAL isolation applied)
        muonPtQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("leg1.pt"),
            indices = cms.vuint32(0,1)
        ),
        muonAbsEtaQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("abs(leg1.eta)"),
            indices = cms.vuint32(0,1)
        ),
        muonTrackIsoQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("leg1.trackIso"),
            indices = cms.vuint32(0,1)
        ),
        muonEcalIsoQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("leg1.ecalIso"),
            indices = cms.vuint32(0,1)
        ),
        tauTrackIsoQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("leg2.chargedHadronIso"),
            indices = cms.vuint32(0,1)
        ),
        tauEcalIsoQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("leg2.photonIso"),
            indices = cms.vuint32(0,1)
        ),
        tauDiscrAgainstMuonsQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("leg2.tauID('againstMuon')"),
            indices = cms.vuint32(0,1)
        ),
        diTauAbsChargeQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("abs(charge)"),
            indices = cms.vuint32(0,1)
        ),
        diTauDphi12QCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("dPhi12"),
            indices = cms.vuint32(0,1)
        ),
        diTauMt1MEtQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("mt1MET"),
            indices = cms.vuint32(0,1)
        ),
        diTauPzetaDiffQCDnoIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
            value = cms.string("pZeta - 1.5*pZetaVis"),
            indices = cms.vuint32(0,1)
        ),
        numDiTausQCDnoIso = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation')
        ),

        # "global" variables not specific to any particular background
        metPt = cms.PSet(
            pluginType = cms.string("PATMetValExtractor"),
            src = cms.InputTag('layer1METs'),
            value = cms.string("pt"),
            indices = cms.vuint32(0)
        ),
        numGlobalMuons = cms.PSet(
            pluginType = cms.string("NumCandidateExtractor"),
            src = cms.InputTag('selectedPatMuonsGlobalIndividual')
        ),

        # variables used for template shape extraction
        diTauMvis12Zmumu = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstZmumuEnriched'),
            value = cms.string("p4Vis.mass"),
            indices = cms.vuint32(0,1)
        ),
        diTauMvis12combinedHypothesesZmumu = cms.PSet(
            pluginType = cms.string("ZtautauVisMassHypothesisMuTauValExtractor"),
            src = cms.InputTag('muTauPairVisMassHypothesesForBgEstZmumuEnriched'),
            srcDiTauCandidates = cms.InputTag('muTauPairsForBgEstZmumuEnriched'),
            indices = cms.vuint32(0,1)
        ),
        diTauMvis12WplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched'),
            value = cms.string("p4Vis.mass"),
            indices = cms.vuint32(0,1)
        ),
        diTauMvis12TTplusJets = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstTTplusJetsEnriched'),
            value = cms.string("p4Vis.mass"),
            indices = cms.vuint32(0,1)
        ),
        diTauMvis12QCDlooseMuonIso = cms.PSet(
            pluginType = cms.string("PATMuTauPairValExtractor"),
            src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation'),
            value = cms.string("p4Vis.mass"),
            indices = cms.vuint32(0,1)
        ),
       
        # event weight variables
        # (accounting for different simulated luminosities of background Monte Carlo samples)
        eventWeight = cms.PSet(
            pluginType = cms.string("ConstObjValExtractor"),
            value = cms.double(1.)
        ),

        # additional event weight variables for correcting "bias"
        # of visible invariant muon + tau-jet mass distribution
        # caused by Mt(muon + tau-jet) transverse mass cut
        # and cut on CDF (Pzeta - 1.5*PzetaVis) variable
        kineEventReweightWplusJets = cms.PSet(
            pluginType = cms.string("KineEventReweightExtractor"),
            weightLookupTable = cms.PSet(
                fileName = kineReweight_fileName,
                meName = cms.string(kineReweight_dqmDirectory + "/" + "WplusJets" + "/" + kineReweight_meName)
            ),
            variables = kineVarExtractor_config.clone(
                src = cms.InputTag('muTauPairsForBgEstWplusJetsEnriched')
            )
        ),
        kineEventReweightTTplusJets = cms.PSet(
            pluginType = cms.string("KineEventReweightExtractor"),
            weightLookupTable = cms.PSet(
                fileName = kineReweight_fileName,
                meName = cms.string(kineReweight_dqmDirectory + "/" + "TTplusJets" + "/" + kineReweight_meName)
            ),
            variables = kineVarExtractor_config.clone(
                src = cms.InputTag('muTauPairsForBgEstTTplusJetsEnriched')
            )
        ),
        kineEventReweightQCD = cms.PSet(
            pluginType = cms.string("KineEventReweightExtractor"),
            weightLookupTable = cms.PSet(
                fileName = kineReweight_fileName,
                meName = cms.string(kineReweight_dqmDirectory + "/" + "QCD" + "/" + kineReweight_meName)
            ),
            variables = kineVarExtractor_config.clone(
                src = cms.InputTag('muTauPairsForBgEstQCDenrichedLooseMuonIsolation')
            )
        ),

        # additional event weight variable for correcting "bias"
        # of muon |eta| distribution caused by cuts on muon track and ECAL isolation variables
        # in QCD background events
        #
        # NOTE: this event weight variable is specific to the tau iso. && id. efficiency analysis
        #       and **not** used for determining the background contributions to the Z --> mu + tau-jet
        #       cross-section analysis via the template technique
        #
        kineEventReweightTauIdEffQCDnoIso = cms.PSet(
            pluginType = cms.string("KineEventReweightExtractor"),
            weightLookupTable = cms.PSet(
                fileName = cms.string('rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_1_2/bgEstKineReweights/muonKineReweightsTauIdEffZtoMuTau.root'),
                meName = cms.string("DQMData/muonKineReweights/QCDbgEnrichedCombIso_data/muonPtVsAbsEta")
            ),
            variables = cms.PSet(
                x = cms.PSet(
                    pluginType = cms.string("PATMuTauPairValExtractor"),
                    src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
                    value = cms.string("abs(leg1.eta)"),
                    indices = cms.vuint32(0)
                ),
                y = cms.PSet(
                    pluginType = cms.string("PATMuTauPairValExtractor"),
                    src = cms.InputTag('muTauPairsForBgEstQCDenrichedNoIsolation'),
                    value = cms.string("leg1.pt"),
                    indices = cms.vuint32(0)
                )
            )
        )
    )
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string('bgEstNtuple.root')
)

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system
#
#__process.source.fileNames = #bgEstSampleFileNames#
#__process.source.skipEvents = cms.untracked.uint32(#skipEvents#)
#__process.maxEvents.input = cms.untracked.int32(#maxEvents#)
#__process.genPhaseSpaceFilter.selector = copy.deepcopy(#genPhaseSpaceCut#)
#__process.TFileService.fileName = #bgEstNtupleOutputFileName#
#__process.ntupleProducer.branches.eventWeight.value = cms.double(#corrFactor#*#intLumiData#/#intLumi#)
#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for switching pat::Tau input
# to different reco::Tau collection stored on AOD
from PhysicsTools.PatAlgos.tools.tauTools import * 

# comment-out to take reco::CaloTaus instead of reco::PFTaus
# as input for pat::Tau production
#switchToCaloTau(process)

# comment-out to take shrinking dR = 5.0/Et(PFTau) signal cone
# instead of fixed dR = 0.07 signal cone reco::PFTaus
# as input for pat::Tau production
switchToPFTauShrinkingCone(process)
#switchToPFTauFixedCone(process)
#--------------------------------------------------------------------------------

process.p = cms.Path(
    process.producePatTuple
   + process.producePatTupleZtoMuTauSpecific
   + process.producePatTupleZtoMuTauSpecificForBgEst
#   + process.printEventContent   # uncomment to enable dump of event content after PAT-tuple production
   + process.selectZtoMuTauEvents
   + process.genPhaseSpaceFilter
   + process.produceBoolEventSelFlags
   + process.selectEventsByBoolEventSelFlags
   + process.ntupleProducer
)

#--------------------------------------------------------------------------------
# disable estimation of systematic uncertainties
from TauAnalysis.Configuration.tools.sysUncertaintyTools import disableSysUncertainties_runZtoMuTau
#
disableSysUncertainties_runZtoMuTau(process)
#--------------------------------------------------------------------------------

# print-out all python configuration parameter information
#print process.dumpPython()
