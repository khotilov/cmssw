import FWCore.ParameterSet.Config as cms

# matches to HLT_Mu9
#muonTriggerMatchHLTMu9 = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
muonTriggerMatchHLTMu9 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
				      matchedCuts = cms.string( 'path( "HLT_Mu9")' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu11
muonTriggerMatchHLTMu11 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
				      matchedCuts = cms.string( 'path( "HLT_Mu11")' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu13
muonTriggerMatchHLTMu13 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
				      matchedCuts = cms.string( 'path( "HLT_Mu13")' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu13_v1
#muonTriggerMatchHLTMu13v1 = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
muonTriggerMatchHLTMu13v1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
				      matchedCuts = cms.string( 'path( "HLT_Mu13_v1")' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu15
muonTriggerMatchHLTMu15 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
				      matchedCuts = cms.string( 'path( "HLT_Mu15")' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu15_v1
muonTriggerMatchHLTMu15v1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
				      matchedCuts = cms.string( 'path( "HLT_Mu15_v1")' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu12_v1
muonTriggerMatchHLTMu12v1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_Mu12_v1" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu15_v2
muonTriggerMatchHLTMu15v2 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_Mu15_v2" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu20_v1
muonTriggerMatchHLTMu20v1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_Mu20_v1" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu24_v1
muonTriggerMatchHLTMu24v1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_Mu24_v1" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Mu30_v1
muonTriggerMatchHLTMu30v1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_Mu30_v1" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_IsoMu12_v1
muonTriggerMatchHLTIsoMu12v1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_IsoMu12_v1" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_IsoMu15_v5
muonTriggerMatchHLTIsoMu15v5 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_IsoMu15_v5" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_IsoMu17_v5
muonTriggerMatchHLTIsoMu17v5 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_IsoMu17_v5" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_IsoMu24_v1
muonTriggerMatchHLTIsoMu24v1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_IsoMu24_v1" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_IsoMu30_v1
muonTriggerMatchHLTIsoMu30v1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                      src     = cms.InputTag( "cleanPatMuons" ),
                                      matched = cms.InputTag( "patTrigger" ),
                                      matchedCuts = cms.string( 'path( "HLT_IsoMu30_v1" )' ),
                                      maxDPtRel = cms.double( 0.5 ),
                                      maxDeltaR = cms.double( 0.3 ),
                                      resolveAmbiguities    = cms.bool( True ),
                                      resolveByMatchQuality = cms.bool( True )
                                      )

# matches to HLT_Ele15_LW_L1R
electronTriggerMatchHLTEle15LWL1R = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                  src     = cms.InputTag( "cleanPatElectrons" ),
                                                  matched = cms.InputTag( "patTrigger" ),
						  matchedCuts = cms.string( 'path( "HLT_Ele15_LW_L1R")' ),
                                                  maxDPtRel = cms.double( 0.5 ),
                                                  maxDeltaR = cms.double( 0.3 ),
                                                  resolveAmbiguities    = cms.bool( True ),
                                                  resolveByMatchQuality = cms.bool( True )
                                                  )

# matches to HLT_Ele15_SW_L1R
electronTriggerMatchHLTEle15SWL1R = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                  src     = cms.InputTag( "cleanPatElectrons" ),
                                                  matched = cms.InputTag( "patTrigger" ),
						  matchedCuts = cms.string( 'path( "HLT_Ele15_SW_L1R")' ),
                                                  maxDPtRel = cms.double( 0.5 ),
                                                  maxDeltaR = cms.double( 0.3 ),
                                                  resolveAmbiguities    = cms.bool( True ),
                                                  resolveByMatchQuality = cms.bool( True )
                                                  )

# matches to HLT_Ele20_SW_L1R
electronTriggerMatchHLTEle20SWL1R = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                  src     = cms.InputTag( "cleanPatElectrons" ),
                                                  matched = cms.InputTag( "patTrigger" ),
						  matchedCuts = cms.string( 'path( "HLT_Ele20_SW_L1R")' ),
                                                  maxDPtRel = cms.double( 0.5 ),
                                                  maxDeltaR = cms.double( 0.3 ),
                                                  resolveAmbiguities    = cms.bool( True ),
                                                  resolveByMatchQuality = cms.bool( True )
                                                  )

# matches to HLT_Ele15_SW_EleId_L1R
electronTriggerMatchHLTEle15SWEleIdL1R = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                  src     = cms.InputTag( "cleanPatElectrons" ),
                                                  matched = cms.InputTag( "patTrigger" ),
						  matchedCuts = cms.string( 'path( "HLT_Ele15_SW_EleId_L1R")' ),
                                                  maxDPtRel = cms.double( 0.5 ),
                                                  maxDeltaR = cms.double( 0.3 ),
                                                  resolveAmbiguities    = cms.bool( True ),
                                                  resolveByMatchQuality = cms.bool( True )
                                                  )

# matches to HLT_Ele15_SW_CaloEleId_L1R
electronTriggerMatchHLTEle15SWCaloEleIdL1R = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Ele15_SW_CaloEleId_L1R")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele17_SW_CaloEleId_L1R 
electronTriggerMatchHLTEle17SWCaloEleIdL1R  = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Ele17_SW_CaloEleId_L1R")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele17_SW_TightEleId_L1R 
electronTriggerMatchHLTEle17SWTightEleIdL1R   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Ele17_SW_TightEleId_L1R")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1
electronTriggerMatchHLTEle17SWTightCaloEleIdEle8HEL1Rv1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2
electronTriggerMatchHLTEle17SWTightCaloEleIdEle8HEL1Rv2   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele17_SW_TighterEleIdIsol_L1R 
electronTriggerMatchHLTEle17SWTighterEleIdIsolL1R    = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Ele17_SW_TighterEleIdIsol_L1R")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele17_SW_TighterEleIdIsol_L1R_v2
electronTriggerMatchHLTEle17SWTighterEleIdIsolL1Rv2    = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Ele17_SW_TighterEleIdIsol_L1R_v2")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele17_SW_TighterEleIdIsol_L1R_v3
electronTriggerMatchHLTEle17SWTighterEleIdIsolL1Rv3    = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Ele17_SW_TighterEleIdIsol_L1R_v3")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele22_SW_TighterEleId_L1R_v2
electronTriggerMatchHLTEle22SWTighterEleIdL1Rv2    = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Ele22_SW_TighterEleId_L1R_v2")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1
electronTriggerMatchHLTEle27CaloIdVTCaloIsoTTrkIdTTrkIsoTv1    = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele45_CaloIdVT_TrkIdT_v1
electronTriggerMatchHLTEle45CaloIdVTTrkIdTv1    = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Ele45_CaloIdVT_TrkIdT_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Ele90_NoSpikeFilter_v1
electronTriggerMatchHLTEle90NoSpikeFilterv1    = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatElectrons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Ele90_NoSpikeFilter_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Photon10_Cleaned_L1R
photonTriggerMatchHLTPhoton10CleanedL1R   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Photon10_Cleaned_L1R")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Photon15_Cleaned_L1R
photonTriggerMatchHLTPhoton15CleanedL1R   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Photon15_Cleaned_L1R")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Photon20_Cleaned_L1R
photonTriggerMatchHLTPhoton20CleanedL1R   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Photon20_Cleaned_L1R")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Photon30_Cleaned_L1R
photonTriggerMatchHLTPhoton30CleanedL1R   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Photon30_Cleaned_L1R")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Photon50_Cleaned_L1R_v1
photonTriggerMatchHLTPhoton50CleanedL1Rv1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Photon50_Cleaned_L1R_v1")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Photon70_Cleaned_L1R_v1
photonTriggerMatchHLTPhoton70CleanedL1Rv1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Photon70_Cleaned_L1R_v1")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_DoublePhoton17_L1R
photonTriggerMatchHLTDoublePhoton17L1R = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                         src     = cms.InputTag( "cleanPatPhotons" ),
                                                         matched = cms.InputTag( "patTrigger" ),
							 matchedCuts = cms.string( 'path( "HLT_DoublePhoton17_L1R")' ),
                                                         maxDPtRel = cms.double( 0.5 ),
                                                         maxDeltaR = cms.double( 0.3 ),
                                                         resolveAmbiguities    = cms.bool( True ),
                                                         resolveByMatchQuality = cms.bool( True )
                                                         )

# matches to HLT_Photon10_L1R
photonTriggerMatchHLTPhoton10L1R   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Photon10_L1R")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Photon75_CaloIdVL_IsoL_v1
photonTriggerMatchHLTPhoton75CaloIdVLIsoLv1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Photon75_CaloIdVL_IsoL_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Photon75_CaloIdVL_v1
photonTriggerMatchHLTPhoton75CaloIdVLv1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Photon75_CaloIdVL_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet15U
jetTriggerMatchHLTJet15U   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet15U")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet30U
jetTriggerMatchHLTJet30U   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet30U")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet50U
jetTriggerMatchHLTJet50U   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet50U")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet70U
jetTriggerMatchHLTJet70U   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet70U")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet70U_v2
jetTriggerMatchHLTJet70Uv2   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet70U_v2")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet100U
jetTriggerMatchHLTJet100U   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet100U")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet100U_v2
jetTriggerMatchHLTJet100Uv2   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet100U_v2")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet140U_v1
jetTriggerMatchHLTJet140Uv1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet140U_v1")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet15U_v3
jetTriggerMatchHLTJet15Uv3   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet15U_v3")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet30U_v3
jetTriggerMatchHLTJet30Uv3   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet30U_v3")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet50U_v3
jetTriggerMatchHLTJet50Uv3   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet50U_v3")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet70U_v3
jetTriggerMatchHLTJet70Uv3   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet70U_v3")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet100U_v3
jetTriggerMatchHLTJet100Uv3   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet100U_v3")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet140U_v3
jetTriggerMatchHLTJet140Uv3   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
							   matchedCuts = cms.string( 'path( "HLT_Jet140U_v3")' ),
                                                           #andOr          = cms.bool( False ),
                                                           #filterIdsEnum  = cms.vstring( '*' ),
                                                           #filterIds      = cms.vint32( 0 ),
                                                           #filterLabels   = cms.vstring( '*' ),
                                                           #pathNames      = cms.vstring( 'HLT_Jet140U_v3' ),
                                                           #collectionTags = cms.vstring( '*' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet30_v1
jetTriggerMatchHLTJet30v1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Jet30_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet60_v1
jetTriggerMatchHLTJet60v1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Jet60_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet80_v1
jetTriggerMatchHLTJet80v1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Jet80_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet110_v1
jetTriggerMatchHLTJet110v1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Jet110_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet150_v1
jetTriggerMatchHLTJet150v1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Jet150_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet190_v1
jetTriggerMatchHLTJet190v1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Jet190_v1")' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet240_v1
jetTriggerMatchHLTJet240v1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Jet240_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

# matches to HLT_Jet370_v1
jetTriggerMatchHLTJet370v1   = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "cleanPatJets" ),
                                                           matched = cms.InputTag( "patTrigger" ),
                                                           matchedCuts = cms.string( 'path( "HLT_Jet370_v1" )' ),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

vgTriggerMatcherElectron = cms.Sequence(
        electronTriggerMatchHLTEle15LWL1R +
        electronTriggerMatchHLTEle15SWL1R +
        electronTriggerMatchHLTEle20SWL1R +
        electronTriggerMatchHLTEle15SWEleIdL1R +
        electronTriggerMatchHLTEle15SWCaloEleIdL1R +
        electronTriggerMatchHLTEle17SWCaloEleIdL1R +
        electronTriggerMatchHLTEle17SWTightEleIdL1R +
	electronTriggerMatchHLTEle17SWTighterEleIdIsolL1R +
	electronTriggerMatchHLTEle17SWTighterEleIdIsolL1Rv2 +
	electronTriggerMatchHLTEle22SWTighterEleIdL1Rv2 +
	electronTriggerMatchHLTEle17SWTightCaloEleIdEle8HEL1Rv1 +
	electronTriggerMatchHLTEle17SWTightCaloEleIdEle8HEL1Rv2 +
	electronTriggerMatchHLTEle17SWTighterEleIdIsolL1Rv3 +
	electronTriggerMatchHLTEle27CaloIdVTCaloIsoTTrkIdTTrkIsoTv1 + 
	electronTriggerMatchHLTEle45CaloIdVTTrkIdTv1 +
	electronTriggerMatchHLTEle90NoSpikeFilterv1
        )

vgTriggerMatcherPhoton = cms.Sequence(
	photonTriggerMatchHLTPhoton10CleanedL1R +
	photonTriggerMatchHLTPhoton15CleanedL1R +
	photonTriggerMatchHLTPhoton20CleanedL1R +
	photonTriggerMatchHLTPhoton30CleanedL1R +
	photonTriggerMatchHLTPhoton50CleanedL1Rv1 +
	photonTriggerMatchHLTPhoton70CleanedL1Rv1 +
        photonTriggerMatchHLTDoublePhoton17L1R +
	photonTriggerMatchHLTPhoton10L1R +
	photonTriggerMatchHLTPhoton75CaloIdVLIsoLv1 +
	photonTriggerMatchHLTPhoton75CaloIdVLv1
	)

vgTriggerMatcherMuon = cms.Sequence(
    muonTriggerMatchHLTMu9 +
    muonTriggerMatchHLTMu11 +
    muonTriggerMatchHLTMu13 +
    muonTriggerMatchHLTMu13v1 +
    muonTriggerMatchHLTMu15 +
    muonTriggerMatchHLTMu15v1 +
    muonTriggerMatchHLTMu12v1 +
    muonTriggerMatchHLTMu15v2 +
    muonTriggerMatchHLTMu20v1 +
    muonTriggerMatchHLTMu24v1 +
    muonTriggerMatchHLTMu30v1 +
    muonTriggerMatchHLTIsoMu12v1 +
    muonTriggerMatchHLTIsoMu15v5 +
    muonTriggerMatchHLTIsoMu17v5 +
    muonTriggerMatchHLTIsoMu24v1 +
    muonTriggerMatchHLTIsoMu30v1 
    )

vgTriggerMatcherJet = cms.Sequence(
	jetTriggerMatchHLTJet15U +
	jetTriggerMatchHLTJet30U +
	jetTriggerMatchHLTJet50U +
	jetTriggerMatchHLTJet70U +
	jetTriggerMatchHLTJet70Uv2 +
	jetTriggerMatchHLTJet100U +
	jetTriggerMatchHLTJet100Uv2 +
	jetTriggerMatchHLTJet140Uv1 +
	jetTriggerMatchHLTJet15Uv3 +
	jetTriggerMatchHLTJet30Uv3 +
	jetTriggerMatchHLTJet50Uv3 +
	jetTriggerMatchHLTJet70Uv3 +
	jetTriggerMatchHLTJet100Uv3 +
	jetTriggerMatchHLTJet140Uv3 +
	jetTriggerMatchHLTJet30v1 +
	jetTriggerMatchHLTJet60v1 +
	jetTriggerMatchHLTJet80v1 +
	jetTriggerMatchHLTJet110v1 +
	jetTriggerMatchHLTJet150v1 +
	jetTriggerMatchHLTJet190v1 +
	jetTriggerMatchHLTJet240v1 +
	jetTriggerMatchHLTJet370v1
	)

vgTriggerMatcher = cms.Sequence(
    vgTriggerMatcherElectron + vgTriggerMatcherMuon + vgTriggerMatcherPhoton + vgTriggerMatcherJet
    )
