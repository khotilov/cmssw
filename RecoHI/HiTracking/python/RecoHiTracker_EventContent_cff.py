import FWCore.ParameterSet.Config as cms

#Full Event content 
RecoHiTrackerFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_hiGlobalPrimTracks_*_*', 
		'keep *_hiSelectedTracks_*_*', 
		'keep *_hiPixel3PrimTracks_*_*', 
		'keep *_hiPixel3ProtoTracks_*_*',	
		'keep *_hiSelectedProtoTracks_*_*',	
		'keep recoVertexs_hiPixelMedianVertex_*_*',
		'keep recoVertexs_hiPixelAdaptiveVertex_*_*',
		'keep recoVertexs_hiSelectedVertex_*_*'		
    )
)

RecoHiTrackerLocalFEVT = cms.PSet(
   outputCommands = cms.untracked.vstring(
   'keep *_*_APVCM_*'
   )
)

#RECO content
RecoHiTrackerRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_hiGlobalPrimTracks_*_*', 
		'keep *_hiSelectedTracks_*_*', 
		#'keep *_hiPixel3PrimTracks_*_*', 		
		'keep recoVertexs_hiPixelMedianVertex_*_*', # do we need to keep these anymore?
		'keep recoVertexs_hiPixelAdaptiveVertex_*_*', # 
		'keep recoVertexs_hiSelectedVertex_*_*'		
    )
)

RecoHiTrackerLocalRECO = cms.PSet(
   outputCommands = cms.untracked.vstring(
   'keep *_*_APVCM_*'
   )
)

#AOD content
RecoHiTrackerAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoTracks_hiSelectedTracks_*_*',
		'keep recoVertexs_hiSelectedVertex_*_*'		
    )
)
