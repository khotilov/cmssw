import FWCore.ParameterSet.Config as cms

hiCentrality = cms.EDFilter("reco::CentralityProducer",

                            doFilter = cms.bool(False),
                            
                            produceHFhits = cms.bool(True),
                            produceHFtowers = cms.bool(True),
                            produceEcalhits = cms.bool(False),
                            produceBasicClusters = cms.bool(True),
                            produceZDChits = cms.bool(False),
                            produceETmidRapidity = cms.bool(True),
                            producePixelhits = cms.bool(True),
                            produceTracks = cms.bool(True),
                            
                            midRapidityRange = cms.double(1),
                            
                            srcHFhits = cms.InputTag("hfreco"),
                            srcTowers = cms.InputTag("towerMaker"),
                            srcEBhits = cms.InputTag("EcalRecHitsEB"),
                            srcEEhits = cms.InputTag("EcalRecHitsEE"),
                            srcBasicClustersEB = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
                            srcBasicClustersEE = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters"),
                            srcZDChits = cms.InputTag(""),
                            srcPixelhits = cms.InputTag("siPixelRecHits"),
                            srcTracks = cms.InputTag("hiSelectedTracks"),                            
                            srcReUse = cms.InputTag("hiCentrality")
                            
                              )


