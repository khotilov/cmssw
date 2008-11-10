import FWCore.ParameterSet.Config as cms

EcalEndcapGeometryEP = cms.ESProducer("EcalEndcapGeometryEP",
                                          applyAlignment = cms.untracked.bool(False)
                                      )

EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
                                          applyAlignment = cms.untracked.bool(False)
                                      )
