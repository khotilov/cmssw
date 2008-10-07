import FWCore.ParameterSet.Config as cms

from SimGeneral.HepPDTESSource.pythiapdt_cfi import *

# The Geometries
from FastSimulation.Configuration.Geometries_cff import *

#The Magnetic Field ESProducer's
from FastSimulation.ParticlePropagator.MagneticFieldMapESProducer_cfi import *

# The muon tracker trajectory, to be fit without rechit refit
from RecoMuon.GlobalTrackingTools.GlobalTrajectoryBuilderCommon_cff import *
GlobalTrajectoryBuilderCommon.TrackRecHitBuilder = 'WithoutRefit'
GlobalTrajectoryBuilderCommon.TrackTransformer.TrackerRecHitBuilder = 'WithoutRefit'

## The condDB setup (the global tag refers to DevDB, IntDB or ProDB whenever needed)
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *




