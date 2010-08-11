import FWCore.ParameterSet.Config as cms

def onia2MuMuPAT(process, GlobalTag, MC=False, HLT='HLT', Filter=True):
    # Setup the process
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True),
        # fileMode = cms.untracked.string('MERGE'),
    )
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = 100
    process.load('Configuration.StandardSequences.GeometryExtended_cff')
    process.load("Configuration.StandardSequences.Reconstruction_cff")
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = GlobalTag

    # Drop the DQM stuff on input
    process.source = cms.Source("PoolSource",
        inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*"),
        fileNames = cms.untracked.vstring()
    )

    # Scraping filter
    process.scrapingFilter = cms.EDFilter("FilterOutScraping",
        applyfilter = cms.untracked.bool(True),
        debugOn = cms.untracked.bool(False),
        numtrack = cms.untracked.uint32(10),
        thresh = cms.untracked.double(0.25)
    )

    # Merge CaloMuons and General Tracks into the collection of reco::Muons
    process.load("RecoMuon.MuonIdentification.mergedMuons_cfi")
    process.mergedMuons.mergeTracks = True
    process.mergedMuons.tracksCut = '(abs(eta) <= 1.3 && pt > 3.3) || (1.3 < abs(eta) <= 2.2 && p > 2.9) || (2.2 < abs(eta) <= 2.4  && pt > 0.8)'
    process.mergedMuons.caloMuonsCut = process.mergedMuons.tracksCut

    # Prune generated particles, needed fot the Tag and Probe
    process.genOniaDecay = cms.EDProducer("GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
            "keep++ pdgId = 443", # J/psi and its daughters
            "keep++ pdgId = 100443", # psi' and its daughters
            "keep++ pdgId = 553", # Upsilon(1S) and its daughters
            "keep++ pdgId = 100553", # Upsilon (2S,3S) and its daughters
        )
    )

    # Make PAT Muons
    process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
    from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
    # with some customization
    if MC:
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it one-to-one
        process.muonMatch.maxDeltaR = 0.05
        process.muonMatch.resolveByMatchQuality = True
    changeRecoMuonInput(process, "mergedMuons")
    useL1MatchingWindowForSinglets(process)
    changeTriggerProcessName(process, HLT)
    switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
    process.muonMatchHLTL3.maxDeltaR = 0.1
    process.muonMatchHLTL3.maxDPtRel = 10.0
    process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
    process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
    process.muonMatchHLTTrackMu.maxDeltaR = 0.1
    process.muonMatchHLTTrackMu.maxDPtRel = 10.0

    # Make a sequence
    process.patMuonSequence = cms.Sequence(
        process.scrapingFilter *
        process.mergedMuons *
        process.genOniaDecay *
        process.patMuonsWithTriggerSequence
    )
    if not MC:
        process.patMuonSequence.remove(process.genOniaDecay)
      
    # Make dimuon candidates
    process.onia2MuMuPatTrkTrk = cms.EDProducer('Onia2MuMuPAT',
        muons = cms.InputTag("patMuonsWithTrigger"),
        beamSpotTag = cms.InputTag("offlineBeamSpot"),
        primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
        higherPuritySelection = cms.string("isGlobalMuon || isTrackerMuon"), ## At least one muon must pass this selection
        lowerPuritySelection  = cms.string("isGlobalMuon || isTrackerMuon"), ## BOTH muons must pass this selection
        dimuonSelection  = cms.string("2 < mass && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25"), ## The dimuon must pass this selection before vertexing
        addCommonVertex = cms.bool(True), ## Embed the full reco::Vertex out of the common vertex fit
        addMuonlessPrimaryVertex = cms.bool(True), ## Embed the primary vertex re-made from all the tracks except the two muons
        addMCTruth = cms.bool(MC),      ## Add the common MC mother of the two muons, if any
        resolvePileUpAmbiguity = cms.bool(False)   ## Order PVs by their vicinity to the J/psi vertex, not by sumPt                            
    )

    # check if there is at least one (inclusive) tracker+tracker di-muon
    process.onia2MuMuPatTrkTrkFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('onia2MuMuPatTrkTrk'),
        minNumber = cms.uint32(1),
    )

    # the onia2MuMu path
    process.Onia2MuMuPAT = cms.Path(
        process.patMuonSequence *
        process.onia2MuMuPatTrkTrk *
        process.onia2MuMuPatTrkTrkFilter
    )

    # Make Tag and Probe pairs for efficiency measurements
    TAG_TRIGGER = "(!triggerObjectMatchesByFilter('hltSingleMu3L3Filtered3').empty())"
    TAG_TRIGGER += "|| (!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3TrackJpsiTrackMassFiltered'))"
    TAG_TRIGGER += "|| (!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5TrackJpsiTrackMassFiltered'))"
    TAG_TRIGGER += "|| (!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu3L2Mu0L3Filtered3'))"
    TAG_TRIGGER += "|| (!triggerObjectMatchesByCollection('hltL3MuonCandidates').empty() && triggerObjectMatchesByCollection('hltL3MuonCandidates').at(0).hasFilterLabel('hltMu5L2Mu0L3Filtered5'))"

    process.tagMuons = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("patMuonsWithTrigger"),
        cut = cms.string("isGlobalMuon && "+TAG_TRIGGER),
    )

    process.probeMuons = cms.EDFilter("PATMuonSelector",
        src = cms.InputTag("patMuonsWithTrigger"),
        cut = cms.string("track.isNonnull"),
    )

    process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
        cut = cms.string('2.6 < mass < 3.5'),
        decay = cms.string('tagMuons@+ probeMuons@-')
    )

    # check if there is at least one Tag and Probe pair
    process.tpPairsFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('tpPairs'),
        minNumber = cms.uint32(1),
    )
    
    # the Tag and Probe path
    process.TagAndProbe = cms.Path(
        process.patMuonSequence *
        process.tagMuons *
        process.probeMuons *
        process.tpPairs *
        process.tpPairsFilter
    )

    if MC:
        process.tagMuonsMCMatch = process.muonMatch.clone(
            src = cms.InputTag("tagMuons"),
            matched = cms.InputTag("genOniaDecay"),
        )
        process.probeMuonsMCMatch = process.tagMuonsMCMatch.clone(src = "probeMuons")
        process.TagAndProbe.replace(process.tpPairs, process.tagMuonsMCMatch * process.probeMuonsMCMatch * process.tpPairs)

    # output
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('onia2MuMuPAT.root'),
        outputCommands = cms.untracked.vstring('drop *',
            'keep *_genOniaDecay_*_Onia2MuMuPAT',                  # generated decay of the Onia
            'keep patMuons_patMuonsWithTrigger_*_Onia2MuMuPAT',    # All PAT muos including general tracks and matches to triggers
            'keep patCompositeCandidates_*__Onia2MuMuPAT',         # PAT di-muons
            'keep patMuons_tagMuons__Onia2MuMuPAT',                # tagMuons for efficiency
            'keep patMuons_probeMuons__Onia2MuMuPAT',              # probeMuons for efficiency
            'keep *_tagMuonsMCMatch__Onia2MuMuPAT',                # tagMuons MC matches for efficiency
            'keep *_probeMuonsMCMatch__Onia2MuMuPAT',              # probeMuons MC matches for efficiency
            'keep recoCompositeCandidates_*__Onia2MuMuPAT',        # RECO di-muons, tpPairs for efficiency
            'keep *_offlinePrimaryVertices_*_*',                   # Primary vertices: you want these to compute impact parameters
            'keep *_offlineBeamSpot_*_*',                          # Beam spot: you want this for the same reason                                   
            'keep edmTriggerResults_TriggerResults_*_*',           # HLT info, per path (cheap)
            'keep l1extraL1MuonParticles_l1extraParticles_*_*',    # L1 info (cheap)
        ),
        SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('Onia2MuMuPAT', 'TagAndProbe') ) if Filter else cms.untracked.PSet()
    )
    process.e = cms.EndPath(process.out)

