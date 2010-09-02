import FWCore.ParameterSet.Config as cms

gamgammumuanalysis = cms.EDAnalyzer("GammaGammaMuMu",
    ElectronCollectionLabel = cms.InputTag("selectedPatElectrons"),
    outfilename = cms.untracked.string('mumu.pat.root'),
    JetCollectionLabel = cms.InputTag("selectedPatJets"),
    PhotonCollectionLabel = cms.InputTag("selectedPatPhotons"),
    CaloTowerLabel = cms.InputTag("towerMaker"),
    GlobalMuonCollectionLabel = cms.InputTag("selectedPatMuons"),
    RecoTrackLabel = cms.InputTag("generalTracks"),
    RecoVertexLabel = cms.InputTag("offlinePrimaryVertices"),
    CastorTowerLabel = cms.InputTag("CastorTowerReco"),
    ZDCRecHitsLabel = cms.InputTag("zdcreco"),                             
    CastorRecHitsLabel = cms.InputTag("castorreco"),
    CaloTowerdR = cms.double(0.3),
    DimuonMindphi = cms.double(0.0),
    MetLabel = cms.InputTag("met"),
                                    
    DimuonMaxdpt = cms.double(2000.0),
    MinMuMuVertexSeparation = cms.double(0.1),
    KeepSameSignDimuons = cms.bool(False),
    ReadMCEffCorrections = cms.bool(False),
    AlgoNames = cms.vstring('TrackerMuonLSAT_Data_CaloMuonProbe_JPsi',
                            'HLT_L1DoubleMuOpen_Data_CaloMuonProbe_JPsi',
                            'TrackerMuonLSAT_MC_CaloMuonProbe_JPsi',
                            'HLT_L1DoubleMuOpen_MC_CaloMuonProbe_JPsi'),
    HLTMenuLabel = cms.string("HLT")                                  
)



