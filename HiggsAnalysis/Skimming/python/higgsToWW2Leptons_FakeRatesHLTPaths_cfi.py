import FWCore.ParameterSet.Config as cms

import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
higgsToWW2LeptonsFakeRatesHLTFilter = copy.deepcopy(hltHighLevel)

#2.2.X
#higgsToWW2LeptonsHLTFilter.HLTPaths = ['HLT_IsoMu11', 'HLT_Mu15_L1Mu7', 'HLT_IsoEle15_L1I', 'HLT_IsoEle18_L1R', 'HLT_DoubleMu3', 
#    'HLT_DoubleIsoEle10_L1I', 'HLT_DoubleIsoEle12_L1R', 'HLT_IsoEle8_IsoMu7', 'HLT_IsoEle10_Mu10_L1R']

#Lum 8E29
#higgsToWW2LeptonsHLTFilter.HLTPaths = ['HLT_Mu3','HLT_Ele10_LW_L1R']

#Lum 1E31
higgsToWW2LeptonsFakeRatesHLTFilter.HLTPaths = ['HLT_Photon15_L1R','HLT_DoublePhoton10_L1R',]
