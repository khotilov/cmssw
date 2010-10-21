import FWCore.ParameterSet.Config as cms

from HLTriggerOffline.Tau.Validation.HLTTauReferences_cfi import *
#from HLTriggerOffline.Tau.Validation.HLTTauValidation_8E29_cfi import *
#from HLTriggerOffline.Tau.Validation.HLTTauValidation_1E31_cfi import *
from HLTriggerOffline.Tau.Validation.HLTTauValidation_6E31_cfi import *
from HLTriggerOffline.Tau.Validation.HLTTauValidation_2E32_cfi import *

HLTTauVal    = cms.Sequence(hltTauRef+hltTauValIdeal2E32+hltTauValIdeal6E31)

# for FS, hltTauRef Producers go into separate "prevalidation" sequence
# (this fixes the "no EDProducer in EndPath" problem)
HLTTauValFS  = cms.Sequence(hltTauValIdeal2E32+hltTauValIdeal6E31)



