import FWCore.ParameterSet.Config as cms
import copy

'''

Sequences for HPS taus

'''

# Define the discriminators for this tau
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByIsolation_cfi                      import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByMVAIsolation_cfi                   import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByLeadingTrackFinding_cfi            import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstElectron_cfi                  import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstElectronMVA_cfi               import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstElectronMVA2_cfi              import *
from RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstMuon_cfi                      import *

# Load helper functions to change the source of the discriminants
from RecoTauTag.RecoTau.TauDiscriminatorTools import *

# Select those taus that pass the HPS selections
#  - pt > 15, mass cuts, tauCone cut
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByHPSSelection_cfi import hpsSelectionDiscriminator
hpsPFTauDiscriminationByDecayModeFinding = hpsSelectionDiscriminator.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer')
    )

# Define decay mode prediscriminant
requireDecayMode = cms.PSet(
    BooleanOperator = cms.string("and"),
    decayMode = cms.PSet(
        Producer = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
        cut = cms.double(0.5)
    )
)

# First apply only charged isolation
hpsPFTauDiscriminationByLooseChargedIsolation = pfRecoTauDiscriminationByIsolation.clone(
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = requireDecayMode.clone(),
    ApplyDiscriminationByECALIsolation = False
)

#Building the prototype for  the Discriminator by Isolation
hpsPFTauDiscriminationByLooseIsolation = pfRecoTauDiscriminationByIsolation.clone(
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = requireDecayMode.clone(),
    ApplyDiscriminationByTrackerIsolation = False,
    ApplyDiscriminationByECALIsolation = True,
    applyOccupancyCut = True
)
hpsPFTauDiscriminationByLooseIsolation.Prediscriminants.preIso = cms.PSet(
    Producer = cms.InputTag("hpsPFTauDiscriminationByLooseChargedIsolation"),
    cut = cms.double(0.5))


# Make an even looser discriminator
hpsPFTauDiscriminationByVLooseChargedIsolation = hpsPFTauDiscriminationByLooseChargedIsolation.clone(
    customOuterCone = cms.double(0.3),
    isoConeSizeForDeltaBeta = cms.double(0.3)
    )
hpsPFTauDiscriminationByVLooseChargedIsolation.qualityCuts.isolationQualityCuts.minTrackPt = 1.5
hpsPFTauDiscriminationByVLooseChargedIsolation.qualityCuts.isolationQualityCuts.minGammaEt = 2.0

hpsPFTauDiscriminationByVLooseIsolation = hpsPFTauDiscriminationByLooseIsolation.clone(
    customOuterCone = cms.double(0.3),
    isoConeSizeForDeltaBeta = cms.double(0.3),
    )
hpsPFTauDiscriminationByVLooseIsolation.qualityCuts.isolationQualityCuts.minTrackPt = 1.5
hpsPFTauDiscriminationByVLooseIsolation.qualityCuts.isolationQualityCuts.minGammaEt = 2.0
hpsPFTauDiscriminationByVLooseIsolation.Prediscriminants.preIso.Producer =  cms.InputTag("hpsPFTauDiscriminationByVLooseChargedIsolation")


hpsPFTauDiscriminationByMediumChargedIsolation = hpsPFTauDiscriminationByLooseChargedIsolation.clone()
hpsPFTauDiscriminationByMediumChargedIsolation.qualityCuts.isolationQualityCuts.minTrackPt = 0.8
hpsPFTauDiscriminationByMediumChargedIsolation.qualityCuts.isolationQualityCuts.minGammaEt = 0.8
hpsPFTauDiscriminationByMediumChargedIsolation.Prediscriminants.preIso = cms.PSet(
    Producer = cms.InputTag("hpsPFTauDiscriminationByLooseChargedIsolation"),
    cut = cms.double(0.5))

hpsPFTauDiscriminationByMediumIsolation = hpsPFTauDiscriminationByLooseIsolation.clone()
hpsPFTauDiscriminationByMediumIsolation.qualityCuts.isolationQualityCuts.minTrackPt = 0.8
hpsPFTauDiscriminationByMediumIsolation.qualityCuts.isolationQualityCuts.minGammaEt = 0.8
hpsPFTauDiscriminationByMediumIsolation.Prediscriminants.preIso.Producer = cms.InputTag("hpsPFTauDiscriminationByMediumChargedIsolation")


hpsPFTauDiscriminationByTightChargedIsolation = hpsPFTauDiscriminationByLooseChargedIsolation.clone()
hpsPFTauDiscriminationByTightChargedIsolation.qualityCuts.isolationQualityCuts.minTrackPt = 0.5
hpsPFTauDiscriminationByTightChargedIsolation.qualityCuts.isolationQualityCuts.minGammaEt = 0.5
hpsPFTauDiscriminationByTightChargedIsolation.Prediscriminants.preIso = cms.PSet(
    Producer = cms.InputTag("hpsPFTauDiscriminationByMediumChargedIsolation"),
    cut = cms.double(0.5))

hpsPFTauDiscriminationByTightIsolation = hpsPFTauDiscriminationByLooseIsolation.clone()
hpsPFTauDiscriminationByTightIsolation.qualityCuts.isolationQualityCuts.minTrackPt = 0.5
hpsPFTauDiscriminationByTightIsolation.qualityCuts.isolationQualityCuts.minGammaEt = 0.5
hpsPFTauDiscriminationByTightIsolation.Prediscriminants.preIso.Producer = cms.InputTag("hpsPFTauDiscriminationByTightChargedIsolation")

hpsPFTauDiscriminationByChargedIsolationSeq = cms.Sequence(
    hpsPFTauDiscriminationByVLooseChargedIsolation*
    hpsPFTauDiscriminationByLooseChargedIsolation*
    hpsPFTauDiscriminationByMediumChargedIsolation*
    hpsPFTauDiscriminationByTightChargedIsolation
)

hpsPFTauDiscriminationByIsolationSeq = cms.Sequence(
    hpsPFTauDiscriminationByVLooseIsolation*
    hpsPFTauDiscriminationByLooseIsolation*
    hpsPFTauDiscriminationByMediumIsolation*
    hpsPFTauDiscriminationByTightIsolation
)

_isolation_types = ['VLoose', 'Loose', 'Medium', 'Tight']
# Now build the sequences that apply PU corrections

# Make Delta Beta corrections (on SumPt quantity)
hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr = hpsPFTauDiscriminationByLooseIsolation.clone(
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    applyDeltaBetaCorrection = True,
    isoConeSizeForDeltaBeta = 0.8,
    deltaBetaFactor = "%0.4f"%(0.0123/0.1687),
    applyOccupancyCut = False,
    applySumPtCut = True,
    customOuterCone = cms.double(0.5)
    )

hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr.maximumSumPtCut=hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minGammaEt

hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr = hpsPFTauDiscriminationByLooseIsolation.clone(
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    applyDeltaBetaCorrection = True,
    isoConeSizeForDeltaBeta = 0.8,
    deltaBetaFactor = "%0.4f"%(0.0123/0.1687),
    applyOccupancyCut = False,
    applySumPtCut = True,
    )
hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr.maximumSumPtCut=hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minGammaEt


hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr = hpsPFTauDiscriminationByMediumIsolation.clone(
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    applyDeltaBetaCorrection = True,
    isoConeSizeForDeltaBeta = 0.8,
    deltaBetaFactor = "%0.4f"%(0.0462/0.1687),
    applyOccupancyCut = False,
    applySumPtCut = True,
    )
hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr.maximumSumPtCut=hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minGammaEt


hpsPFTauDiscriminationByTightIsolationDBSumPtCorr = hpsPFTauDiscriminationByTightIsolation.clone(
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    applyDeltaBetaCorrection = True,
    isoConeSizeForDeltaBeta = 0.8,
    deltaBetaFactor = "%0.4f"%(0.0772/0.1687),
    applyOccupancyCut = False,
    applySumPtCut = True,
    )
hpsPFTauDiscriminationByTightIsolationDBSumPtCorr.maximumSumPtCut=hpsPFTauDiscriminationByTightIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minGammaEt

hpsPFTauDiscriminationByIsolationSeqDBSumPtCorr = cms.Sequence(
    hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr*
    hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr*
    hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr*
    hpsPFTauDiscriminationByTightIsolationDBSumPtCorr
    )


hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr = hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr.clone(
    ApplyDiscriminationByTrackerIsolation = True,
    ApplyDiscriminationByECALIsolation = True,
    deltaBetaFactor = "%0.4f"%(0.0772/0.1687),
    applyOccupancyCut = False,
    applySumPtCut = True,
    maximumSumPtCut = 3.0,
    Prediscriminants = requireDecayMode.clone()
    )
hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minTrackPt = 0.5
hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minGammaEt = 0.5

hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr = hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr.clone(
    applySumPtCut = False,
    storeRawSumPt = cms.bool(True)
)

hpsPFTauDiscriminationByRawChargedIsolationDBSumPtCorr = hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr.clone(
    applySumPtCut = False,
    ApplyDiscriminationByECALIsolation = False,
    storeRawSumPt = cms.bool(True)
)

hpsPFTauDiscriminationByRawGammaIsolationDBSumPtCorr = hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr.clone(
    applySumPtCut = False,
    ApplyDiscriminationByTrackerIsolation = False,
    storeRawSumPt = cms.bool(True)
)

hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr = hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr.clone(
    ApplyDiscriminationByTrackerIsolation = True,
    ApplyDiscriminationByECALIsolation = True,
    deltaBetaFactor = "%0.4f"%(0.0772/0.1687),
    applyOccupancyCut = False,
    applySumPtCut = True,
    maximumSumPtCut = 2.0,
    Prediscriminants = requireDecayMode.clone()
    )
hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minTrackPt = 0.5
hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minGammaEt = 0.5

hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr = hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr.clone(
    ApplyDiscriminationByTrackerIsolation = True,
    ApplyDiscriminationByECALIsolation = True,
    deltaBetaFactor = "%0.4f"%(0.0772/0.1687),
    applyOccupancyCut = False,
    applySumPtCut = True,
    maximumSumPtCut = 1.0,
    Prediscriminants = requireDecayMode.clone()
    )
hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minTrackPt = 0.5
hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minGammaEt = 0.5

hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr = hpsPFTauDiscriminationByTightIsolationDBSumPtCorr.clone(
    ApplyDiscriminationByTrackerIsolation = True,
    ApplyDiscriminationByECALIsolation = True,
    deltaBetaFactor = "%0.4f"%(0.0772/0.1687),
    applyOccupancyCut = False,
    applySumPtCut = True,
    maximumSumPtCut = 0.8,
    Prediscriminants = requireDecayMode.clone()
    )
hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minTrackPt = 0.5
hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr.qualityCuts.isolationQualityCuts.minGammaEt = 0.5

hpsPFTauDiscriminationByCombinedIsolationSeqDBSumPtCorr = cms.Sequence(
    hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr*
    hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr*
    hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr*
    hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr
    )

# Define MVA based isolation discrimators
hpsPFTauDiscriminationByIsolationMVAraw = pfRecoTauDiscriminationByMVAIsolation.clone(
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = requireDecayMode.clone(),
    returnMVA = cms.bool(True),
    )

hpsPFTauDiscriminationByLooseIsolationMVA = hpsPFTauDiscriminationByDecayModeFinding.clone(
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string("and"),
        mva = cms.PSet(
            Producer = cms.InputTag('hpsPFTauDiscriminationByIsolationMVAraw'),
            cut = cms.double(0.795)
        )
    ))
hpsPFTauDiscriminationByMediumIsolationMVA = copy.deepcopy(hpsPFTauDiscriminationByLooseIsolationMVA)
hpsPFTauDiscriminationByMediumIsolationMVA.Prediscriminants.mva.cut = cms.double(0.884)
hpsPFTauDiscriminationByTightIsolationMVA = copy.deepcopy(hpsPFTauDiscriminationByLooseIsolationMVA)
hpsPFTauDiscriminationByTightIsolationMVA.Prediscriminants.mva.cut = cms.double(0.921)

from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets as dummy
kt6PFJetsForRhoComputationVoronoi = dummy.clone(
    doRhoFastjet = True,
    voronoiRfact = 0.9
)

hpsPFTauDiscriminationByMVAIsolationSeq = cms.Sequence(
    kt6PFJetsForRhoComputationVoronoi*
    hpsPFTauDiscriminationByIsolationMVAraw*
    hpsPFTauDiscriminationByLooseIsolationMVA*
    hpsPFTauDiscriminationByMediumIsolationMVA*
    hpsPFTauDiscriminationByTightIsolationMVA
    )

#copying discriminator against electrons and muons
hpsPFTauDiscriminationByLooseElectronRejection = pfRecoTauDiscriminationAgainstElectron.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = noPrediscriminants,
    PFElectronMVA_maxValue = cms.double(0.6)
    )

hpsPFTauDiscriminationByMediumElectronRejection = pfRecoTauDiscriminationAgainstElectron.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = noPrediscriminants,
    ApplyCut_EcalCrackCut = cms.bool(True)
    )

hpsPFTauDiscriminationByTightElectronRejection = pfRecoTauDiscriminationAgainstElectron.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = noPrediscriminants,
    ApplyCut_EcalCrackCut = cms.bool(True),
    ApplyCut_BremCombined = cms.bool(True)
    )

hpsPFTauDiscriminationByLooseMuonRejection = pfRecoTauDiscriminationAgainstMuon.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = noPrediscriminants
    )

hpsPFTauDiscriminationByMediumMuonRejection = pfRecoTauDiscriminationAgainstMuon.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = noPrediscriminants,
    discriminatorOption = cms.string('noAllArbitrated')
    )

hpsPFTauDiscriminationByTightMuonRejection = pfRecoTauDiscriminationAgainstMuon.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = noPrediscriminants,
    discriminatorOption = cms.string('noAllArbitratedWithHOP')
    )

hpsPFTauDiscriminationByMVAElectronRejection = pfRecoTauDiscriminationAgainstElectronMVA.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = requireDecayMode.clone(),
)

# Additionally require that the MVA electrons pass electron medium
# (this discriminator was used on the training sample)
hpsPFTauDiscriminationByMVAElectronRejection.Prediscriminants.electronMedium = \
        cms.PSet(
            Producer = cms.InputTag('hpsPFTauDiscriminationByMediumElectronRejection'),
            cut = cms.double(0.5)
        )

hpsPFTauDiscriminationByMVA2rawElectronRejection = pfRecoTauDiscriminationAgainstElectronMVA2.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = requireDecayMode.clone()
)

from RecoTauTag.RecoTau.RecoTauDiscriminantCutMultiplexer_cfi import recoTauDiscriminantCutMultiplexer
hpsPFTauDiscriminationByMVA2VLooseElectronRejection = recoTauDiscriminantCutMultiplexer.clone(
    PFTauProducer = cms.InputTag('hpsPFTauProducer'),
    Prediscriminants = requireDecayMode.clone(),
    toMultiplex = cms.InputTag('hpsPFTauDiscriminationByMVA2rawElectronRejection'),
    key = cms.InputTag('hpsPFTauDiscriminationByMVA2rawElectronRejection:category'),
    mapping = cms.VPSet(
        cms.PSet(
            category = cms.uint32(0), # minMVA1prongNoEleMatchBL
            cut = cms.double(-0.141383)
        ),
        cms.PSet(
            category = cms.uint32(1), # minMVA1prongBL
            cut = cms.double(-0.122568)
        ),
        cms.PSet(
            category = cms.uint32(2), # minMVA1prongStripsWOgsfBL
            cut = cms.double(-0.138286)
        ),
        cms.PSet(
            category = cms.uint32(3), # minMVA1prongStripsWgsfWOpfEleMvaBL
            cut = cms.double(-0.100279)
        ),
        cms.PSet(
            category = cms.uint32(4), # minMVA1prongStripsWgsfWpfEleMvaBL
            cut = cms.double(-0.116113)
        ),
        cms.PSet(
            category = cms.uint32(5), # minMVA1prongNoEleMatchEC
            cut = cms.double(-0.191557)
        ),
        cms.PSet(
            category = cms.uint32(6), # minMVA1prongEC
            cut = cms.double(-0.0921955)
        ),
        cms.PSet(
            category = cms.uint32(7), # minMVA1prongStripsWOgsfEC
            cut = cms.double(-0.0936173)
        ),
        cms.PSet(
            category = cms.uint32(8), # minMVA1prongStripsWgsfWOpfEleMvaEC
            cut = cms.double(-0.119732)
        ),
        cms.PSet(
            category = cms.uint32(9), # minMVA1prongStripsWgsfWpfEleMvaEC
            cut = cms.double(-0.1042)
        )
    )
)

hpsPFTauDiscriminationByMVA2LooseElectronRejection = copy.deepcopy(hpsPFTauDiscriminationByMVA2VLooseElectronRejection)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[0].cut = cms.double(-0.0639254)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[1].cut = cms.double(-0.0220708)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[2].cut = cms.double(-0.102071)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[3].cut = cms.double(-0.0233814)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[4].cut = cms.double(-0.0391565)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[5].cut = cms.double(-0.142564)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[6].cut = cms.double(+0.00982555)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[7].cut = cms.double(-0.0596019)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[8].cut = cms.double(-0.0381238)
hpsPFTauDiscriminationByMVA2LooseElectronRejection.mapping[9].cut = cms.double(-0.100381)

hpsPFTauDiscriminationByMVA2MediumElectronRejection = copy.deepcopy(hpsPFTauDiscriminationByMVA2VLooseElectronRejection)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[0].cut = cms.double(+0.011729)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[1].cut = cms.double(+0.0203646)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[2].cut = cms.double(+0.177502)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[3].cut = cms.double(+0.0103449)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[4].cut = cms.double(+0.257798)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[5].cut = cms.double(-0.0966083)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[6].cut = cms.double(-0.0466023)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[7].cut = cms.double(+0.0467638)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[8].cut = cms.double(+0.0863876)
hpsPFTauDiscriminationByMVA2MediumElectronRejection.mapping[9].cut = cms.double(+0.233436)

hpsPFTauDiscriminationByMVA2TightElectronRejection = copy.deepcopy(hpsPFTauDiscriminationByMVA2VLooseElectronRejection)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[0].cut = cms.double(+0.0306715)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[1].cut = cms.double(+0.992195)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[2].cut = cms.double(+0.308324)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[3].cut = cms.double(-0.0370998)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[4].cut = cms.double(+0.864643)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[5].cut = cms.double(+0.0832094)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[6].cut = cms.double(+0.791665)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[7].cut = cms.double(+0.675537)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[8].cut = cms.double(+0.87047)
hpsPFTauDiscriminationByMVA2TightElectronRejection.mapping[9].cut = cms.double(+0.233711)

# Define the HPS selection discriminator used in cleaning
hpsSelectionDiscriminator.PFTauProducer = cms.InputTag("combinatoricRecoTaus")

import RecoTauTag.RecoTau.RecoTauCleanerPlugins as cleaners

hpsPFTauProducerSansRefs = cms.EDProducer(
    "RecoTauCleaner",
    src = cms.InputTag("combinatoricRecoTaus"),
    cleaners = cms.VPSet(
        # Prefer taus that dont' have charge == 3
        cleaners.unitCharge,
        # Prefer taus that are within DR<0.1 of the jet axis
        cleaners.matchingConeCut,
        # Prefer taus that pass HPS selections
        cms.PSet(
            name = cms.string("HPS_Select"),
            plugin = cms.string("RecoTauDiscriminantCleanerPlugin"),
            src = cms.InputTag("hpsSelectionDiscriminator"),
        ),
        cleaners.combinedIsolation
    )
)

hpsPFTauProducer = cms.EDProducer(
    "RecoTauPiZeroUnembedder",
    src = cms.InputTag("hpsPFTauProducerSansRefs")
)

produceHPSPFTaus = cms.Sequence(
    hpsSelectionDiscriminator
    #*hpsTightIsolationCleaner
    #*hpsMediumIsolationCleaner
    #*hpsLooseIsolationCleaner
    #*hpsVLooseIsolationCleaner
    *hpsPFTauProducerSansRefs
    *hpsPFTauProducer
)

produceAndDiscriminateHPSPFTaus = cms.Sequence(
    produceHPSPFTaus*
    hpsPFTauDiscriminationByDecayModeFinding*
    hpsPFTauDiscriminationByChargedIsolationSeq*
    hpsPFTauDiscriminationByIsolationSeq*
    #hpsPFTauDiscriminationByIsolationSeqRhoCorr*
    #hpsPFTauDiscriminationByIsolationSeqCustomRhoCorr*
    hpsPFTauDiscriminationByIsolationSeqDBSumPtCorr*
    hpsPFTauDiscriminationByMVAIsolationSeq*

    hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr*
    hpsPFTauDiscriminationByRawChargedIsolationDBSumPtCorr*
    hpsPFTauDiscriminationByRawGammaIsolationDBSumPtCorr*

    hpsPFTauDiscriminationByCombinedIsolationSeqDBSumPtCorr*
    hpsPFTauDiscriminationByLooseElectronRejection*
    hpsPFTauDiscriminationByMediumElectronRejection*
    hpsPFTauDiscriminationByTightElectronRejection*
    hpsPFTauDiscriminationByMVAElectronRejection*
    hpsPFTauDiscriminationByMVA2rawElectronRejection*
    hpsPFTauDiscriminationByMVA2VLooseElectronRejection*
    hpsPFTauDiscriminationByMVA2LooseElectronRejection*
    hpsPFTauDiscriminationByMVA2MediumElectronRejection*
    hpsPFTauDiscriminationByMVA2TightElectronRejection*
    hpsPFTauDiscriminationByLooseMuonRejection*
    hpsPFTauDiscriminationByMediumMuonRejection*
    hpsPFTauDiscriminationByTightMuonRejection
)
