import FWCore.ParameterSet.Config as cms

# $Id: FourVectorHLTriggerOffline_cfi.py,v 1.2 2009/01/26 05:46:15 berryhil Exp $
hltResults = cms.EDFilter("FourVectorHLTriggerOffline",
    dirname = cms.string("HLT/FourVector/"),
    plotAll = cms.untracked.bool(False),
    ptMax = cms.untracked.double(100.0),
    ptMin = cms.untracked.double(0.0),
    paths = cms.VPSet(
# single jet triggers
             cms.PSet(
              pathname = cms.string("HLT_Jet50"),
              denompathname = cms.string("HLT_Jet30"),  
             ),
             cms.PSet(
              pathname = cms.string("HLT_Jet110"),
              denompathname = cms.string("HLT_Jet50"),  
             )
            ),
                          
     # this is I think MC and CRUZET4
    triggerSummaryLabel = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),
    gtObjectMapRecordLabel = cms.InputTag("hltL1GtObjectMap","","HLT"),
    l1GTRRLabel = cms.InputTag("hltGtDigis","","HLT"),
    processname = cms.string("HLT")

    # this is data (CRUZET I or II best guess)
    #triggerSummaryLabel = cms.InputTag("hltTriggerSummaryAOD","","FU"),
    #triggerResultsLabel = cms.InputTag("TriggerResults","","FU"),
    #processname = cms.string("FU")

 )
