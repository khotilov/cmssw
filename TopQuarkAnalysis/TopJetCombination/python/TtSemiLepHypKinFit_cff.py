import FWCore.ParameterSet.Config as cms

#
# produce kinFit hypothesis with all necessary 
# ingredients
#

## std sequence to perform kinematic fit
from TopQuarkAnalysis.TopKinFitter.TtSemiLepKinFitProducer_Muons_cfi import *

## configure kinFit hypothesis
from TopQuarkAnalysis.TopJetCombination.TtSemiLepHypKinFit_cfi import *

## make hypothesis
makeHypothesis_kinFit = cms.Sequence(kinFitTtSemiEvent *
                                     ttSemiLepHypKinFit)

