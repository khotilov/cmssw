import FWCore.ParameterSet.Config as cms

TreesEle = cms.EDAnalyzer('SingleTopSystematicsTreesDumper',                              
systematics = cms.untracked.vstring("BTagUp","BTagDown","MisTagUp","MisTagDown","JESUp","JESDown","UnclusteredMETUp","UnclusteredMETDown"),
rateSystematics = cms.untracked.vstring("WLightRate","TTBarRate","WqqRate","WqRate"),
leptonsID = cms.InputTag("nTupleElectrons","tightElectronsSimpleEleId95cIso"),  

channelInfo = cms.PSet(
    crossSection = cms.untracked.double(20.93),
    channel = cms.untracked.string("TChannel"),
#    originalEvents = cms.untracked.double(14800000),
    originalEvents = cms.untracked.double(480000),
    finalLumi = cms.untracked.double(14.5),
    MTWCut = cms.untracked.double(50.0),#Default 50.0 GeV
    loosePtCut = cms.untracked.double(30.0),#Default 30.0 GeV
    RelIsoCut = cms.untracked.double(0.1)
    ),


#Part of the kin quantities:
leptonsEta = cms.InputTag("nTupleElectrons","tightElectronsEta"),  
leptonsPt = cms.InputTag("nTupleElectrons","tightElectronsPt"),  
leptonsPhi = cms.InputTag("nTupleElectrons","tightElectronsPhi"),  
leptonsEnergy = cms.InputTag("nTupleElectrons","tightElectronsE"),  
leptonsCharge = cms.InputTag("nTupleElectrons","tightElectronsCharge"),  
leptonsRelIso = cms.InputTag("nTupleElectrons","tightElectronsRelIso"),  

looseMuonsRelIso = cms.InputTag("nTupleLooseElectrons","looseElectronsRelIso"),  
looseElectronsRelIso = cms.InputTag("nTupleLooseMuons","looseMuonsRelIso"),  

leptonsFlavour = cms.untracked.string("electron"),

jetsPt = cms.InputTag("nTupleTopJetsPF","topJetsPFPt"),  
jetsPhi = cms.InputTag("nTupleTopJetsPF","topJetsPFPhi"),  
jetsEta = cms.InputTag("nTupleTopJetsPF","topJetsPFEta"),  
jetsEnergy = cms.InputTag("nTupleTopJetsPF","topJetsPFE"),  

jetsBTagAlgo = cms.InputTag("nTupleTopJetsPF","topJetsPFTrackCountingHighPur"),  
jetsAntiBTagAlgo =  cms.InputTag("nTupleTopJetsPF","topJetsPFTrackCountingHighEff"),  
jetsFlavour = cms.InputTag("nTupleTopJetsPF","topJetsPFFlavour"),   

jetsCorrTotal = cms.InputTag("nTupleTopJetsPF","topJetsPFJetCorrTotal"),   

METPhi = cms.InputTag("nTuplePatMETsPF","patMETsPFPhi"),
METPt = cms.InputTag("nTuplePatMETsPF","patMETsPFPt"),

UnclusteredMETPx = cms.InputTag("UnclusteredMETPF","UnclusteredMETPx"),
UnclusteredMETPy = cms.InputTag("UnclusteredMETPF","UnclusteredMETPy"),

)


TreesMu = TreesEle.clone(
    channelInfo = cms.PSet(
        crossSection = cms.untracked.double(20.93),
            channel = cms.untracked.string("TChannel"),
        #    originalEvents = cms.untracked.double(14800000),
            originalEvents = cms.untracked.double(480000),
            finalLumi = cms.untracked.double(14.5),
            MTWCut = cms.untracked.double(40.0),#Default 50.0 GeV
        RelIsoCut = cms.untracked.double(0.05),
            ),

    leptonsEta = cms.InputTag("nTupleMuons","tightMuonsEta"),  
    leptonsPt = cms.InputTag("nTupleMuons","tightMuonsPt"),  
    leptonsPhi = cms.InputTag("nTupleMuons","tightMuonsPhi"),  
    leptonsEnergy = cms.InputTag("nTupleMuons","tightMuonsE"),  
    leptonsCharge = cms.InputTag("nTupleMuons","tightMuonsCharge"),  
    leptonsRelIso = cms.InputTag("nTupleMuons","tightMuonsRelIso"),  
    leptonsID = cms.InputTag("nTupleElectrons","tightElectronsSimpleEleId95cIso"),  
    leptonsFlavour = cms.untracked.string("muon"),

    
    )

