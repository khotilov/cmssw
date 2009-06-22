import FWCore.ParameterSet.Config as cms

# First register all the hit matching algorithms, then specify preferred
# ones at end.

# Hit matching algorithm a

HitMatchingAlgorithm_a_PSimHit_ = cms.ESProducer(
    "HitMatchingAlgorithm_a_PSimHit_"
    )

HitMatchingAlgorithm_a_PixelDigi_ = cms.ESProducer(
    "HitMatchingAlgorithm_a_PixelDigi_"
    )

HitMatchingAlgorithm_a_TTHit_ = cms.ESProducer(
	"HitMatchingAlgorithm_a_TTHit_"
	)


# Global geometry hit matching algorithm

HitMatchingAlgorithm_globalgeometry_PSimHit_ = cms.ESProducer(
    "HitMatchingAlgorithm_globalgeometry_PSimHit_",
    ipWidth = cms.double(15.0),
    minPtThreshold = cms.double(2.0)
    )

HitMatchingAlgorithm_globalgeometry_PixelDigi_ = cms.ESProducer(
    "HitMatchingAlgorithm_globalgeometry_PixelDigi_",
    ipWidth = cms.double(200.0),
    minPtThreshold = cms.double(2.0)
    )

HitMatchingAlgorithm_globalgeometry_TTHit_ = cms.ESProducer(
    "HitMatchingAlgorithm_globalgeometry_TTHit_",
    ipWidth = cms.double(200.0),
    minPtThreshold = cms.double(2.0)
    )


# Pixel ray hit matching algorithm

HitMatchingAlgorithm_pixelray_PixelDigi_ = cms.ESProducer(
    "HitMatchingAlgorithm_pixelray_PixelDigi_",
    ipWidth = cms.double(40.0),
    minPtThreshold = cms.double(2.0)
    )

HitMatchingAlgorithm_pixelray_TTHit_ = cms.ESProducer(
    "HitMatchingAlgorithm_pixelray_TTHit_",
    ipWidth = cms.double(40.0),
    minPtThreshold = cms.double(2.0)
    )


# Window hit matching algorithm

HitMatchingAlgorithm_window_PixelDigi_ = cms.ESProducer(
    "HitMatchingAlgorithm_window_PixelDigi_",
    ipWidth = cms.double(15.0),
    minPtThreshold = cms.double(5.0),
    RowResolution = cms.double(0.5),
    ColResolution = cms.double(0.5)
    )

HitMatchingAlgorithm_window_TTHit_ = cms.ESProducer(
    "HitMatchingAlgorithm_window_TTHit_",
    ipWidth = cms.double(15.0),
    minPtThreshold = cms.double(5.0),
    RowResolution = cms.double(0.5),
    ColResolution = cms.double(0.5)
    )


# Thresholds hit matching algorithm

HitMatchingAlgorithm_thresholds_PixelDigi_ = cms.ESProducer(
    "HitMatchingAlgorithm_thresholds_PixelDigi_",
    Thresholds = cms.VPSet(
        cms.PSet(
            # thresholds for layer 0
            Layer = cms.uint32(0),
            RowCuts = cms.vuint32(35,51,145,162,600),
            RowOffsets = cms.vuint32(7,7,8,8,9),
            RowWindows = cms.vuint32(2,2,2,2,2),
            ColumnCutMin = cms.uint32(0),
            ColumnCutMax = cms.uint32(2)
            ),
        cms.PSet(
            # thresholds for layer 1
            Layer = cms.uint32(1),
            RowCuts = cms.vuint32(24,140,170,600),
            RowOffsets = cms.vuint32(7,8,8,9),
            RowWindows = cms.vuint32(3,2,3,2),                          
            ColumnCutMin = cms.uint32(0),
            ColumnCutMax = cms.uint32(2)
            )
        )
    )

HitMatchingAlgorithm_thresholds_TTHit_ = cms.ESProducer(
    "HitMatchingAlgorithm_thresholds_TTHit_",
    Thresholds = cms.VPSet(
        cms.PSet(
            # thresholds for layer 0
            Layer = cms.uint32(0),
            RowCuts = cms.vuint32(35,51,145,162,600),
            RowOffsets = cms.vuint32(7,7,8,8,9),
            RowWindows = cms.vuint32(2,2,2,2,2),
            ColumnCutMin = cms.uint32(0),
            ColumnCutMax = cms.uint32(2)
            ),
        cms.PSet(
            # thresholds for layer 1
            Layer = cms.uint32(1),
            RowCuts = cms.vuint32(24,140,170,600),
            RowOffsets = cms.vuint32(7,8,8,9),
            RowWindows = cms.vuint32(3,2,3,2),                          
            ColumnCutMin = cms.uint32(0),
            ColumnCutMax = cms.uint32(2)
            )
        )
    )


# Set the preferred hit matching algorithms.
# We prefer the global geometry algorithm for now in order not to break
# anything. Override with process.myPreferPSimHit = ..., etc. in your
# configuration.

myPreferPSimHit = cms.ESPrefer("HitMatchingAlgorithm_globalgeometry_PSimHit_")

myPreferPixelDigi = cms.ESPrefer("HitMatchingAlgorithm_globalgeometry_PixelDigi_")

myPreferTTHit = cms.ESPrefer("HitMatchingAlgorithm_globalgeometry_TTHit_")
