import FWCore.ParameterSet.Config as cms

import Geometry.TrackerGeometryBuilder.trackerGeometry_cfi
#
# This cff provides a TrackerGeometry with the label 'idealForDigi' that is for sure matching
# the ideal one and thus should be used in the digitisers.
#
idealForDigiTrackerGeometry = Geometry.TrackerGeometryBuilder.trackerGeometry_cfi.TrackerDigiGeometryESModule.clone()
# The es_module providing fake (i.e. empty) alignment constants:
from Alignment.CommonAlignmentProducer.fakeForIdealAlignmentProducer_cfi import *
#replace idealForDigiTrackerGeometry.applyAlignment = true # GF: See below
# Replace although false is default as protection against foreseen removal:
idealForDigiTrackerGeometry.applyAlignment = False
# Label of the produced TrackerGeometry:
idealForDigiTrackerGeometry.appendToDataLabel = 'idealForDigi'
# Alignments are looked for with this label:
idealForDigiTrackerGeometry.alignmentsLabel = 'fakeForIdeal'

