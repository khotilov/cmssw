import FWCore.ParameterSet.Config as cms

hltanalysis = cms.EDAnalyzer("HLTAnalyzer",
### Generator objects
##    mctruth = cms.InputTag("genParticleCandidates"),
    mctruth = cms.InputTag("genParticles"),
    genEventScale = cms.InputTag("genEventScale"),
### Reconstructed objects
    genjets = cms.InputTag("iterativeCone5GenJets"),
    genmet = cms.InputTag("genMet"),
    recjets = cms.InputTag("hltIterativeCone5CaloJets"),
    recmet = cms.InputTag("hltMet"),
    ht = cms.InputTag("hltHtMet"),
    calotowers = cms.InputTag("towerMaker"),
    muon = cms.InputTag("muons"),
    Electron = cms.InputTag("pixelMatchGsfElectrons"),
##    Photon = cms.InputTag("correctedPhotons"),
    Photon = cms.InputTag("photons"),
### Trigger objects
##    l1GctCounts = cms.InputTag("l1GctEmulDigis"),
    l1GctCounts = cms.InputTag("hltGctDigis"),
##    l1GtObjectMapRecord = cms.InputTag("l1GtEmulDigis"),
    l1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap"),
##    l1GtReadoutRecord = cms.InputTag("l1GmtEmulDigis"),
    l1GtReadoutRecord = cms.InputTag("hltGtDigis"),
    l1extramc = cms.string('hltL1extraParticles'),
##    hltresults = cms.InputTag("TriggerResults"),
    hltresults = cms.InputTag("TriggerResults::HLT"),
### Muon OpenHLT objects                             
    MuCandTag2 = cms.InputTag("hltL2MuonCandidates"),
    MuCandTag3 = cms.InputTag("hltL3MuonCandidates"),
    MuIsolTag3 = cms.InputTag("hltL3MuonIsolations"),
    MuIsolTag2 = cms.InputTag("hltL2MuonIsolations"),
    MuLinkTag = cms.InputTag("hltL3Muons"),
### Egamma OpenHLT objects                             
    CandIso = cms.InputTag("hltL1IsoRecoEcalCandidate"),
    CandNonIso = cms.InputTag("hltL1NonIsoRecoEcalCandidate"),
    EcalIso =cms.InputTag("hltL1IsolatedPhotonEcalIsol"),
    EcalNonIso =cms.InputTag("hltL1NonIsolatedPhotonEcalIsol"),
    HcalIsoPho =cms.InputTag("hltL1IsolatedPhotonHcalIsol"),
    HcalNonIsoPho =cms.InputTag("hltL1NonIsolatedPhotonHcalIsol"),
    IsoPhoTrackIsol =cms.InputTag("hltL1IsoPhotonTrackIsol"),
    NonIsoPhoTrackIsol =cms.InputTag("hltL1NonIsoPhotonTrackIsol"),
    HcalIsoEle =cms.InputTag("hltL1IsolatedElectronHcalIsol"),
    HcalNonIsoEle  =cms.InputTag("hltL1NonIsolatedElectronHcalIsol"),
    #### Standard or Startup windows                         
##    IsoElectrons =cms.InputTag("hltPixelMatchElectronsL1Iso"),
##    NonIsoElectrons =cms.InputTag("hltPixelMatchElectronsL1NonIso"),
##    PixelSeedL1Iso =cms.InputTag("hltL1IsoElectronPixelSeeds"),
##    PixelSeedL1NonIso =cms.InputTag("hltL1NonIsoElectronPixelSeeds"),
##    IsoEleTrackIsol =cms.InputTag("hltL1IsoElectronTrackIsol"),
##    NonIsoEleTrackIsol =cms.InputTag("hltL1NonIsoElectronTrackIsol"),
    IsoElectrons =cms.InputTag("hltPixelMatchStartUpElectronsL1Iso"),
    NonIsoElectrons =cms.InputTag("hltPixelMatchStartUpElectronsL1NonIso"),
    PixelSeedL1Iso =cms.InputTag("hltL1IsoStartUpElectronPixelSeeds"),
    PixelSeedL1NonIso =cms.InputTag("hltL1NonIsoStartUpElectronPixelSeeds"),
    IsoEleTrackIsol =cms.InputTag("hltL1IsoStartUpElectronTrackIsol"),
    NonIsoEleTrackIsol =cms.InputTag("hltL1NonIsoStartupElectronTrackIsol"),
    ### Large windows
    IsoElectronsLargeWindows =cms.InputTag("hltPixelMatchElectronsL1IsoLargeWindow"),
    NonIsoElectronsLargeWindows =cms.InputTag("hltPixelMatchElectronsL1NonIsoLargeWindow"),
    PixelSeedL1IsoLargeWindows =cms.InputTag("hltL1IsoLargeWindowElectronPixelSeeds"),
    PixelSeedL1NonIsoLargeWindows =cms.InputTag("hltL1NonIsoLargeWindowElectronPixelSeeds"),
    IsoEleTrackIsolLargeWindows =cms.InputTag("hltL1IsoLargeWindowElectronTrackIsol"),
    NonIsoEleTrackIsolLargeWindows =cms.InputTag("hltL1NonIsoLargeWindowElectronTrackIsol"),
### Tau HLT related objects
    HLTTau =cms.InputTag("TauOpenHLT"),
    ########
    RunParameters = cms.PSet(
        GenJetMin = cms.double(0.0),
        Monte = cms.bool(True),
        CalJetMin = cms.double(0.0),
        HistogramFile = cms.string('TEST.root'),
        EtaMin = cms.double(-5.2),
        Debug = cms.bool(False),
        EtaMax = cms.double(5.2)
    )

)



