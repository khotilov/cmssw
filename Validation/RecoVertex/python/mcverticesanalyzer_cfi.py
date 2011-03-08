import FWCore.ParameterSet.Config as cms

mcverticesanalyzer = cms.EDAnalyzer("MCVerticesAnalyzer",
                                    pileupSummaryCollection = cms.InputTag("addPileupInfo"),
                                    mcTruthCollection = cms.InputTag("generator"),
                                    useWeight = cms.bool(False),
                                    weightProduct = cms.InputTag("mcvertexweight")
                           )

