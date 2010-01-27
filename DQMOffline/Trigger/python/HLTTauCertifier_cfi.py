import FWCore.ParameterSet.Config as cms

hltTauOfflineCertification = cms.EDFilter("HLTTauCertifier",
                                   targetDir = cms.string("HLT/EventInfo/reportSummaryContents"),
                                   targetME  = cms.string("HLT_Tau"),
                                   inputMEs = cms.vstring(
                                      "HLT/TauOffline/Inclusive/DoubleLooseIsoTauTau/TriggerBits",
                                      "HLT/TauOffline/Inclusive/DoubleLooseIsoTauTau/TriggerBits",

                                   ),
                                   setBadRunOnWarnings = cms.bool(False),
                                   setBadRunOnErrors   = cms.bool(True)
)



