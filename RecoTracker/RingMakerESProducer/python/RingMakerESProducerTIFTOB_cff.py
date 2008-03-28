import FWCore.ParameterSet.Config as cms

# geometry
from Geometry.CMSCommonData.cmsIdealGeometryXML_cfi import *
# tracker geometry
from Geometry.TrackerGeometryBuilder.trackerGeometry_cfi import *
# tracker numbering
from Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi import *
import copy
from RecoTracker.RingMakerESProducer.RingMakerESProducer_cfi import *
# rings esproducer
ringsTIFTOB = copy.deepcopy(rings)
ringsTIFTOB.Configuration = 'TIFTOB'
ringsTIFTOB.RingAsciiFileName = 'rings_tiftob.dat'
ringsTIFTOB.ComponentName = 'TIFTOB'

