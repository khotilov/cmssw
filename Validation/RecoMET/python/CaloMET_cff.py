import FWCore.ParameterSet.Config as cms

# File: CaloMET.cfi
# Author: B. Scurlock
# Date: 03.04.2008
#
# Fill validation histograms for MET. Assumes met, metNoHF, metOpt, metOptNoHF are in the event
from Validation.RecoMET.CaloMET_cfi import *
analyzeCaloMET = cms.Sequence(metAnalyzer*metNoHFAnalyzer*metOptAnalyzer*metOptNoHFAnalyzer)

