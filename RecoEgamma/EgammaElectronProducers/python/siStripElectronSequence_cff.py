import FWCore.ParameterSet.Config as cms

#
# $Id: siStripElectronSequence_cff.py,v 1.3 2008/05/19 23:54:05 rpw Exp $
#
# complete sequence to 1) make siStripElectrons
#                      2) do tracking based on these siStripElectrons
#                      3) associate tracks to SiStripElectrons
#
# Created by Shahram Rahatlou, University of Rome & INFN, 4 Aug 2006
# based on the cfg files from Jim Pivarsky, Cornell
#
# tracker geometry
from Geometry.TrackerGeometryBuilder.trackerGeometry_cfi import *
# tracker numbering
from Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi import *
# standard geometry
from Geometry.CMSCommonData.cmsIdealGeometryXML_cfi import *
# magnetic field
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEParmError_cfi import *
# use super clusters and si strip hits to make siStripElectrons
from RecoEgamma.EgammaElectronProducers.siStripElectrons_cfi import *
# do tracking seeded by siStripElectrons
from RecoEgamma.EgammaElectronProducers.egammaCTFFinalFitWithMaterial_cff import *
# asscoiate tracks to siStripElectrons
from RecoEgamma.EgammaElectronProducers.siStripElectronToTrackAssociator_cfi import *
siStripElectronSequence = cms.Sequence(siStripElectrons*egammaCTFFinalFitWithMaterial*siStripElectronToTrackAssociator)

