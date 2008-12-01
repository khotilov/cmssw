import FWCore.ParameterSet.Config as cms

TrackerSystematicMisalignments = cms.EDAnalyzer("TrackerSystematicMisalignments",
								  # grab an existing geometry
								  fromDBGeom = cms.untracked.bool(True),
								  #epsilons
								  radialEpsilon = cms.untracked.double(-999.0), # default 5e-4 ~ 600 um
								  telescopeEpsilon = cms.untracked.double(-999.0), # default 5e-4 ~ 600 um
								  layerRotEpsilon = cms.untracked.double(-999.0), # 9.43e-6
								  bowingEpsilon = cms.untracked.double(-999.0), #6.77e-9
								  zExpEpsilon = cms.untracked.double(-999.0), # 2.02e-4
								  twistEpsilon = cms.untracked.double(-999.0),	# 2.04e-6
								  ellipticalEpsilon = cms.untracked.double(-999.0), # 5e-4
								  skewEpsilon = cms.untracked.double(-999.0), # 5.5e-2
								  saggitaEpsilon = cms.untracked.double(-999.0) #5.0e-4
)

