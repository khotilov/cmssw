import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# preselection of events considered in data-driven background estimation methods
#--------------------------------------------------------------------------------

electronsLooseIdForBgEst = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("cleanPatElectrons"),                     
    cut = cms.string('electronID("loose") > 0'),
    filter = cms.bool(False)
)

electronsRobustIdForBgEst = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("cleanPatElectrons"),                     
    cut = cms.string('electronID("robust") > 0'),
    filter = cms.bool(False)
)

electronsTrkLooseIsolationForBgEst = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag('selectedLayer1ElectronsEcalIsoLooseIsolationCumulative'),                                        
    cut = cms.string('gsfTrack.isNonnull'),
    filter = cms.bool(False)
)

electronsTrkIsoForBgEst = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag('selectedLayer1ElectronsPt15Cumulative'),
    cut = cms.string('trackIso < 2.'),
    filter = cms.bool(False)
)

electronsEcalIsoForBgEst = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag('electronsTrkIsoForBgEst'),                          
    cut = cms.string('ecalIso < 2.'),
    filter = cms.bool(False)
)

electronsTrkTightIsolationForBgEst = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag('electronsEcalIsoForBgEst'),                                        
    cut = cms.string('gsfTrack.isNonnull'),
    filter = cms.bool(False)
)

electronsConversionVetoForBgEst = cms.EDFilter("PATElectronConversionFinder",
    src = cms.InputTag('electronsTrkTightIsolationForBgEst'),
    trackSource = cms.InputTag('generalTracks'),
    conversionSource = cms.InputTag('conversions'),
    photonSource = cms.InputTag('photons'),
    cotThetaCut = cms.double(0.045),
    docaElecTrack = cms.double(0),
    dRElecTrack = cms.double(0.1),
    doPixCut = cms.bool(True),
    useInnerParsForElec = cms.bool(True),
    useInnerParsForTrks = cms.bool(True),
    useConversionColl = cms.bool(False),                                                           
    nTrkMax = cms.double(1),
    doHists = cms.bool(False),
    filter = cms.bool(False)
)

selectElectronsForBgEst = cms.Sequence(
    electronsLooseIdForBgEst * electronsRobustIdForBgEst
   * electronsTrkLooseIsolationForBgEst
   * electronsTrkIsoForBgEst * electronsEcalIsoForBgEst * electronsTrkTightIsolationForBgEst * electronsConversionVetoForBgEst
)

