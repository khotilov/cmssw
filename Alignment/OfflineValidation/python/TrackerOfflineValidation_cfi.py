import FWCore.ParameterSet.Config as cms

# Put here the modules you want the cfg file to use,
# then include this file in your cfg file.
# i.e. in Validator.cfg replace 'module demo = Validator {} '
# with 'include "anlyzerDir/Validator/data/Validator.cfi" '.
# (Remember that filenames are case sensitive.)
TrackerOfflineValidation = cms.EDFilter("TrackerOfflineValidation",
    Tracks = cms.InputTag("TrackRefitter"),
    trajectoryInput           = cms.string('TrackRefitter'),
    localCoorHistosOn         = cms.bool(False),
    moduleLevelHistsTransient = cms.bool(False),
    stripYResiduals           = cms.bool(False),                                        
    overlappOn                = cms.bool(False),                                      
    useFwhm                   = cms.bool(True),
    useFit                    = cms.bool(False),
    useCombinedTrajectory     = cms.bool(False),
    useOverflowForRMS         = cms.bool(False),
    # Normalized X Residuals, normal local coordinates (Strip)
    TH1NormXResStripModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-5.0), xmax = cms.double(5.0)
    ),

    # X Residuals, normal local coordinates (Strip)                      
    TH1XResStripModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    ),

    # Normalized X Residuals, native coordinates (Strip)
    TH1NormXprimeResStripModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-5.0), xmax = cms.double(5.0)
    ),

    # X Residuals, native coordinates (Strip)
    TH1XprimeResStripModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    ),

    # Normalized Y Residuals, native coordinates (Strip -> hardly defined)
    TH1NormYResStripModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-5.0), xmax = cms.double(5.0)
    ),
    # -> very broad distributions expected                                         
    TH1YResStripModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-10.0), xmax = cms.double(10.0)
    ),

    # Normalized X residuals normal local coordinates (Pixel)                                        
    TH1NormXResPixelModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-5.0), xmax = cms.double(5.0)
    ),
    # X residuals normal local coordinates (Pixel)                                        
    TH1XResPixelModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    ),
    # Normalized X residuals native coordinates (Pixel)                                        
    TH1NormXprimeResPixelModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-5.0), xmax = cms.double(5.0)
    ),
    # X residuals native coordinates (Pixel)                                        
    TH1XprimeResPixelModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    ),                                        
    # Normalized Y residuals native coordinates (Pixel)                                         
    TH1NormYResPixelModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-5.0), xmax = cms.double(5.0)
    ),
    # Y residuals native coordinates (Pixel)                                         
    TH1YResPixelModules = cms.PSet(
        Nbinx = cms.int32(100), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    )
)


