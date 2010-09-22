import FWCore.ParameterSet.Config as cms

'''

Configuration for combinatoric PFTau producer plugins.

Note that this plugin produces many taus for each PFJet!
To be useful the output from this module must be cleaned
using an implementation of the RecoTauCleaner module.

Author: Evan K. Friis, UC Davis


'''

# N.B. for combinatoric taus that worst-case scaling
# is (maxTracks choose dmTracks) * (maxPiZeros choose dmPiZeros)
#
# So for decay mode 11 (3 tracks, 1 pizero), with 10 for both 
#
# (10 choose 3) * (10 choose 1) = 1200!

_combinatoricTauConfig = cms.PSet(
    name = cms.string("combinatoric"),
    plugin = cms.string("RecoTauBuilderCombinatoricPlugin"),
    pfCandSrc = cms.InputTag("particleFlow"),
    decayModes = cms.VPSet(
        cms.PSet(
            # One prong no pizero mode
            nCharged = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            maxTracks = cms.uint32(10),
            maxPiZeros = cms.uint32(0),
        ),
        cms.PSet(
            #One prong one pizero mode
            nCharged = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            maxTracks = cms.uint32(10),
            maxPiZeros = cms.uint32(10),
        ),
        cms.PSet(
            #One prong two pizero mode
            nCharged = cms.uint32(1),
            nPiZeros = cms.uint32(2),
            maxTracks = cms.uint32(10),
            maxPiZeros = cms.uint32(5),
        ),
        cms.PSet(
            # Three prong no pizero mode
            nCharged = cms.uint32(3),
            nPiZeros = cms.uint32(0),
            maxTracks = cms.uint32(10),
            maxPiZeros = cms.uint32(0),
        ),
        cms.PSet(
            # Three prong one pizero mode
            nCharged = cms.uint32(3),
            nPiZeros = cms.uint32(1),
            maxTracks = cms.uint32(10),
            maxPiZeros = cms.uint32(2),
        )
    )
)

combinatoricRecoTaus = cms.EDProducer(
    "RecoTauProducer",
    jetSrc = cms.InputTag("ak5PFJets"),
    piZeroSrc = cms.InputTag("ak5PFJetsRecoTauPiZeros"),
    builders = cms.VPSet(
        _combinatoricTauConfig,
    ),
    modifiers = cms.VPSet(),
)
