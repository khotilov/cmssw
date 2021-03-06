import FWCore.ParameterSet.Config as cms
import copy

#from SimGeneral.HepPDTESSource.pdt_cfi import *

#--------------------------------------------------------------------------------
# produce data-formats providing information 
# about compatibility of an electron + muon pair
# with the hypothesis of being a pair of muons,
# resulting from a Z --> mu+ mu- decay
#--------------------------------------------------------------------------------

elecMuPairZmumuHypotheses = cms.EDProducer("ZllHypothesisElecMuProducer",
    diCandidatePairSource = cms.InputTag('selectedElecMuPairsPzetaDiffCumulative'),

    genLeptonsFromZsSource = cms.InputTag('genMuonsFromZs'),

    # try all possible combinations
    caloJetSource = cms.InputTag('ak5CaloJets'),
    pfJetSource = cms.InputTag('ak5PFJets'),                                       
    trackSource = cms.InputTag('generalTracks'),
    gsfElectronSource = cms.InputTag('gsfElectrons'),
    gsfTrackSource = cms.InputTag('electronGsfTracks'),

    # combined (inner) tracks only                                      
    #trackSource = cms.InputTag('generalTracks'),
                                          
    tkminPixelHits = cms.int32(1),
    tkminTrackerHits = cms.int32(8),	
    tkmaxChi2 = cms.double(100.),

    dRmatch = cms.double(0.5),

    verbosity = cms.untracked.int32(0)                                      
)

elecMuPairVisMassHypotheses = cms.EDProducer("ZtautauVisMassHypothesisElecMuProducer",
    diCandidatePairSource = elecMuPairZmumuHypotheses.diCandidatePairSource,

    ZllHypotheses = cms.VPSet(
        cms.PSet(
            src = cms.InputTag('elecMuPairZmumuHypotheses'),
            minZllMass = cms.double(85.),
            maxZllMass = cms.double(100.)
        )
    )
)

elecMuPairZmumuHypothesesLooseElectronIsolation = copy.deepcopy(elecMuPairZmumuHypotheses)
elecMuPairZmumuHypothesesLooseElectronIsolation.diCandidatePairSource = cms.InputTag('selectedElecMuPairsPzetaDiffLooseElectronIsolationCumulative')

elecMuPairVisMassHypothesesLooseElectronIsolation = copy.deepcopy(elecMuPairVisMassHypotheses)
elecMuPairVisMassHypothesesLooseElectronIsolation.diCandidatePairSource = elecMuPairZmumuHypothesesLooseElectronIsolation.diCandidatePairSource
elecMuPairVisMassHypothesesLooseElectronIsolation.ZllHypotheses[0].src = cms.InputTag('elecMuPairZmumuHypothesesLooseElectronIsolation')

produceElecMuPairZmumuHypotheses = cms.Sequence(
    elecMuPairZmumuHypotheses
   * elecMuPairVisMassHypotheses
   * elecMuPairZmumuHypothesesLooseElectronIsolation
   * elecMuPairVisMassHypothesesLooseElectronIsolation
)
                                                  
