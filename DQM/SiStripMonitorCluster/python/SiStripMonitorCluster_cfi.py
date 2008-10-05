import FWCore.ParameterSet.Config as cms

# SiStripMonitorCluster
SiStripMonitorCluster = cms.EDFilter("SiStripMonitorCluster",
    # by default do not write out any file with histograms
    # can overwrite this in .cfg file with: replace SiStripMonitorCluster.OutputMEsInRootFile = true
    ClusterProducer = cms.string('siStripClusters'),
                                     
    OutputMEsInRootFile = cms.bool(False),
    OutputFileName = cms.string('SiStripMonitorCluster.root'),
                                     
    CreateTrendMEs = cms.bool(False),
    ResetMEsEachRun = cms.bool(False),

    StripQualityLabel = cms.string(''),

    SelectAllDetectors = cms.bool(False),
    ShowMechanicalStructureView = cms.bool(True),

    ClusterLabel = cms.string(''),

    Trending = cms.PSet(
        UpdateMode = cms.int32(1),
        Nbins      = cms.int32(10),
        ymax       = cms.double(10000.0),
        Steps      = cms.int32(10),
        xmax       = cms.double(10.0),
        xmin       = cms.double(0.0),
        ymin       = cms.double(0.0)
    ),
    TH1ClusterNoise = cms.PSet(
        Nbinx          = cms.int32(20),
        xmin           = cms.double(-0.5),
        xmax           = cms.double(9.5),
        layerswitchon  = cms.bool(True),
        moduleswitchon = cms.bool(True)
    ),

    TH1NrOfClusterizedStrips = cms.PSet(
        Nbinx          = cms.int32(20),
        xmin           = cms.double(-0.5),
        xmax           = cms.double(99.5),
        layerswitchon  = cms.bool(True),
        moduleswitchon = cms.bool(True)
    ),
    TH1ClusterPos = cms.PSet(
        Nbinx          = cms.int32(768),
        xmin           = cms.double(-0.5),
        xmax           = cms.double(767.5),
        layerswitchon  = cms.bool(False),
        moduleswitchon = cms.bool(True)
    ),
    TH1ModuleLocalOccupancy = cms.PSet(
        Nbinx          = cms.int32(20),
        xmin           = cms.double(-0.5),
        xmax           = cms.double(0.95),
        layerswitchon  = cms.bool(True),
        moduleswitchon = cms.bool(True)
    ),
    TH1nClusters = cms.PSet(
        Nbinx          = cms.int32(11),
        xmin           = cms.double(-0.5),
        xmax           = cms.double(10.5),
        layerswitchon  = cms.bool(False),
        moduleswitchon = cms.bool(True)
    ),
    TH1ClusterStoN = cms.PSet(
        Nbinx          = cms.int32(100),
        xmin           = cms.double(-0.5),
        xmax           = cms.double(299.5),
        layerswitchon  = cms.bool(False),
        moduleswitchon = cms.bool(True)
    ),
    TH1ClusterCharge = cms.PSet(
        Nbinx          = cms.int32(200),
        xmin           = cms.double(-0.5),        
        xmax           = cms.double(799.5),
        layerswitchon  = cms.bool(True),
        moduleswitchon = cms.bool(True)
    ),
    TH1ClusterWidth = cms.PSet(
        Nbinx          = cms.int32(20),
        xmin           = cms.double(-0.5),
        xmax           = cms.double(19.5),
        layerswitchon  = cms.bool(True),        
        moduleswitchon = cms.bool(True)
    ),

    TProfNumberOfCluster = cms.PSet(
        Nbinx            = cms.int32(100),
        xmin             = cms.double(-0.5),
        xmax             = cms.double(499.5),
        layerswitchon    = cms.bool(True),        
        moduleswitchon   = cms.bool(False)        
    ),
      
    TProfClusterWidth    = cms.PSet(
        Nbinx            = cms.int32(100),
        xmin             = cms.double(-0.5),
        xmax             = cms.double(499.5),
        layerswitchon    = cms.bool(True),        
        moduleswitchon   = cms.bool(False)        
    ),
                                     
    ClusterConditions = cms.PSet(
        minWidth   = cms.double(0.0),
        On         = cms.bool(True),
        maxStoN    = cms.double(10000.0),
        minStoN    = cms.double(0.0),
        maxWidth   = cms.double(10000.0)
    ),
                                     
    #select detectors
    detectorson = cms.PSet(
        tidon = cms.bool(True),
        tibon = cms.bool(True),
        tecon = cms.bool(True),
        tobon = cms.bool(True)
    ),

    ShowControlView = cms.bool(False),
    ShowReadoutView = cms.bool(False)
)
