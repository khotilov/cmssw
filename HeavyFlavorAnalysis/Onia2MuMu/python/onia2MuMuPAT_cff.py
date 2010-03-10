import FWCore.ParameterSet.Config as cms

##    __  __                        ____      _       __  __                       
##   |  \/  | ___ _ __ __ _  ___   / ___|__ _| | ___ |  \/  |_   _  ___  _ __  ___ 
##   | |\/| |/ _ \ '__/ _` |/ _ \ | |   / _` | |/ _ \| |\/| | | | |/ _ \| '_ \/ __|
##   | |  | |  __/ | | (_| |  __/ | |__| (_| | | (_) | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\___|_|  \__, |\___|  \____\__,_|_|\___/|_|  |_|\__,_|\___/|_| |_|___/
##                    |___/                                                        
##   
## ==== Merge CaloMuons into the collection of reco::Muons  ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
mergedMuons = cms.EDProducer("CaloMuonMerger",
    muons     = cms.InputTag("muons"), 
    caloMuons = cms.InputTag("calomuons"),
    minCaloCompatibility = calomuons.minCaloCompatibility
)

##    __  __  ____   _____           _   _     
##   |  \/  |/ ___| |_   _| __ _   _| |_| |__  
##   | |\/| | |       | || '__| | | | __| '_ \ 
##   | |  | | |___    | || |  | |_| | |_| | | |
##   |_|  |_|\____|   |_||_|   \__,_|\__|_| |_|
##                                             
##   
from PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi import muonMatch
muonMatch.src = 'mergedMuons'

##    __  __       _          ____   _  _____   __  __                       
##   |  \/  | __ _| | _____  |  _ \ / \|_   _| |  \/  |_   _  ___  _ __  ___ 
##   | |\/| |/ _` | |/ / _ \ | |_) / _ \ | |   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | (_| |   <  __/ |  __/ ___ \| |   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|_|\_\___| |_| /_/   \_\_|   |_|  |_|\__,_|\___/|_| |_|___/
##                                                                           
##   
### ==== Make PAT Muons ====
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
patMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone(
    muonSource = 'mergedMuons',
    # embed the tracks, so we don't have to carry them around
    embedTrack          = True,
    embedCombinedMuon   = True,
    embedStandAloneMuon = True,
    # then switch off some features we don't need
    #addTeVRefits = False, ## <<--- this doesn't work. PAT bug ??
    embedPickyMuon = False,
    embedTpfmsMuon = False, 
    userIsolation = cms.PSet(),   # no extra isolation beyond what's in reco::Muon itself
    isoDeposits = cms.PSet(), # no heavy isodeposits
    addGenMatch = True,       # no mc: T&P doesn't take it from here anyway.
    embedGenMatch = True,
    genParticleMatch = 'muonMatch'
)

##    _____     _                         __  __       _       _     
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__  
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \ 
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|
##                |___/ |___/                                        
##   
### ==== Unpack trigger, and match ====
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
patTrigger.onlyStandAlone = True

### ==== Then perform a match for all HLT triggers of interest
from PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi import muonTriggerMatchHLTMu3
muonTriggerMatchHLTMu3.src = 'patMuonsWithoutTrigger'
muonTriggerMatchHLTMu3.andOr = False # i.e. 'AND'

muonMatchHLTMu3       = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_Mu3" ])
muonMatchHLTMu5       = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_Mu5" ])
muonMatchHLTDoubleMu0 = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_DoubleMu0" ])
muonMatchHLTDoubleMu3  = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_DoubleMu3" ])
muonMatchHLTOniaMu = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_Onia" ], collectionTags = ['hltL3MuonCandidates::HLT'   ] )
muonMatchHLTOniaPx = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_Onia" ], collectionTags = ['hltOniaPixelTrackCands::HLT'] )
muonMatchHLTOniaTk = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_Onia" ], collectionTags = ['hltOniaCtfTrackCands::HLT'  ] )
muonMatchHLTOniaMu8E29 = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_Onia_8E29" ], collectionTags = ['hltL3MuonCandidates::HLT'       ] )
muonMatchHLTOniaPx8E29 = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_Onia_8E29" ], collectionTags = ['hltOniaPixelTrackCands8E29::HLT'] )
muonMatchHLTOniaTk8E29 = muonTriggerMatchHLTMu3.clone(pathNames = [ "HLT_Onia_8E29" ], collectionTags = ['hltOniaCtfTrackCands8E29::HLT'  ] )

### == For HLT triggers which are just L1s, we need a different matcher
from MuonAnalysis.MuonAssociators.muonHLTL1Match_cfi import muonHLTL1Match

muonMatchHLTL1MuOpen = muonHLTL1Match.clone(
    # Input
    src     = muonTriggerMatchHLTMu3.src,
    matched = muonTriggerMatchHLTMu3.matched,
    # Matching criteria
    maxDeltaR   = cms.double(0.3),
    # Which trigger to use
    pathName = cms.string('HLT_L1MuOpen'),
)
muonMatchHLTL1DoubleMuOpen = muonMatchHLTL1MuOpen.clone(pathName = 'HLT_L1DoubleMuOpen')


## ==== Embed ====
patMuonsWithTrigger = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
    src     = cms.InputTag(  "patMuonsWithoutTrigger" ),
    matches = cms.VInputTag( 
        # HLT Matches
        cms.InputTag('muonMatchHLTMu3'),
        cms.InputTag('muonMatchHLTMu5'),
        cms.InputTag('muonMatchHLTDoubleMu0'),
        cms.InputTag('muonMatchHLTDoubleMu3'),
        cms.InputTag('muonMatchHLTOniaMu'),
        cms.InputTag('muonMatchHLTOniaPx'),
        cms.InputTag('muonMatchHLTOniaTk'),
        cms.InputTag('muonMatchHLTOniaMu8E29'),
        cms.InputTag('muonMatchHLTOniaPx8E29'),
        cms.InputTag('muonMatchHLTOniaTk8E29'),
        # L1 Matches
        cms.InputTag('muonMatchHLTL1MuOpen','propagatedReco'), # will match if the muon did propagate to station 2
        cms.InputTag('muonMatchHLTL1MuOpen'),
        cms.InputTag('muonMatchHLTL1DoubleMuOpen'),
    )

)

## ==== Trigger Sequence ====
patTriggerMatching = cms.Sequence(
    patTrigger * 
    ( muonMatchHLTL1MuOpen   +
      muonMatchHLTMu3        +
      muonMatchHLTMu5        +
      muonMatchHLTL1DoubleMuOpen  +
      muonMatchHLTDoubleMu0  +
      muonMatchHLTDoubleMu3  +
      muonMatchHLTOniaMu     +
      muonMatchHLTOniaPx     +
      muonMatchHLTOniaTk     +
      muonMatchHLTOniaMu8E29 +
      muonMatchHLTOniaPx8E29 +
      muonMatchHLTOniaTk8E29 
    ) *
    patMuonsWithTrigger
)

### ==== Apply some final selection (none by default) ====
patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(""), 
)

### ==== Sequence ====
patMuonSequence = cms.Sequence( 
    mergedMuons *
    muonMatch *
    patMuonsWithoutTrigger *
    patTriggerMatching *
    patMuons  
)

##    ____  _       __  __                       
##   |  _ \(_)     |  \/  |_   _  ___  _ __  ___ 
##   | | | | |_____| |\/| | | | |/ _ \| '_ \/ __|
##   | |_| | |_____| |  | | |_| | (_) | | | \__ \
##   |____/|_|     |_|  |_|\__,_|\___/|_| |_|___/
##                                               
##   
import HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi 

onia2MuMuPatGlbGlb = HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi.onia2MuMuPAT.clone()

onia2MuMuPatGlbTrk = HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi.onia2MuMuPAT.clone()
onia2MuMuPatGlbTrk.lowerPuritySelection  = cms.string("isGlobalMuon || isTrackerMuon")

onia2MuMuPatGlbCal = HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi.onia2MuMuPAT.clone()
onia2MuMuPatGlbCal.lowerPuritySelection  = cms.string("isGlobalMuon || isTrackerMuon || (track.isNonnull && isCaloMuon)")

onia2MuMuPatTrkTrk = HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi.onia2MuMuPAT.clone()
onia2MuMuPatTrkTrk.higherPuritySelection  = cms.string("isGlobalMuon || isTrackerMuon")
onia2MuMuPatTrkTrk.lowerPuritySelection   = cms.string("isGlobalMuon || isTrackerMuon")

onia2MuMuPatTrkCal = HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi.onia2MuMuPAT.clone()
onia2MuMuPatTrkCal.higherPuritySelection  = cms.string("isGlobalMuon || isTrackerMuon")
onia2MuMuPatTrkCal.lowerPuritySelection   = cms.string("isGlobalMuon || isTrackerMuon || (track.isNonnull && isCaloMuon)")

##    _____           _     
##   |_   _|__   ___ | |___ 
##     | |/ _ \ / _ \| / __|
##     | | (_) | (_) | \__ \
##     |_|\___/ \___/|_|___/
##                          
##   
def onia2MuMu_isNotMC(process):
    process.patMuonSequence.remove(muonMatch)
    process.patMuonsWithoutTrigger.addGenMatch = False
    process.onia2MuMuPatGlbGlb.addMCTruth = False
    process.onia2MuMuPatGlbCal.addMCTruth = False
    process.onia2MuMuPatGlbTrk.addMCTruth = False
    process.onia2MuMuPatTrkTrk.addMCTruth = False
    process.onia2MuMuPatTrkCal.addMCTruth = False
