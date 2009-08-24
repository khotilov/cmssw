import FWCore.ParameterSet.Config as cms
import copy

# define configuration parameters for submission of AH --> e + mu jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)

from TauAnalysis.Configuration.recoSampleDefinitionsZtoElecMu_cfi import *
plotsOutputFileNameZtautau = cms.string('plotsAHtoElecMu_Ztautau.root')
plotsOutputFileNameWplusJets = cms.string('plotsAHtoElecMu_WplusJets_partXX.root')
plotsOutputFileNameZee = cms.string('plotsAHtoElecMu_Zee.root')
plotsOutputFileNameZmumu = cms.string('plotsAHtoElecMu_Zmumu.root')
plotsOutputFileNameZplusJets = cms.string('plotsAHtoElecMu_ZplusJets.root')
plotsOutputFileNameInclusivePPmuX = cms.string('plotsAHtoElecMu_InclusivePPmuX.root')
plotsOutputFileNamePPmuXptGt20 = cms.string('plotsAHtoElecMu_PPmuXptGt20.root')
patTupleOutputFileNameZtautau = cms.untracked.string('patTupleAHtoElecMu_Ztautau.root')
patTupleOutputFileNameWplusJets = cms.untracked.string('patTupleAHtoElecMu_WplusJets_partXX.root')
patTupleOutputFileNameZee = cms.untracked.string('patTupleAHtoElecMu_Zee.root')
patTupleOutputFileNameZmumu = cms.untracked.string('patTupleAHtoElecMu_Zmumu.root')
patTupleOutputFileNameZplusJets = cms.untracked.string('patTupleAHtoElecMu_ZplusJets.root')
patTupleOutputFileNameInclusivePPmuX = cms.untracked.string('patTupleAHtoElecMu_InclusivePPmuX.root')
patTupleOutputFileNamePPmuXptGt20 = cms.untracked.string('patTupleAHtoElecMu_PPmuXptGt20.root')

#--------------------------------------------------------------------------------
# AH115 --> tau+ tau- sample
#
fileNamesAH115tautau = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115v2_tautau_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115v2_tautau_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115v2_tautau_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115v2_tautau_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115v2_tautau_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115v2_tautau_6.root'
)
genPhaseSpaceCutAH115tautau = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameAH115tautau = cms.string('plotsAHtoElecMu_AH115tautau.root')
patTupleOutputFileNameAH115tautau = cms.untracked.string('patTupleAHtoElecMu_AH115tautau.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# AH115bb --> tau+ tau- sample
#
fileNamesAH115bbtautau = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_5.root'
)
genPhaseSpaceCutAH115bbtautau = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameAH115bbtautau = cms.string('plotsAHtoElecMu_AH115bbtautau.root')
patTupleOutputFileNameAH115bbtautau = cms.untracked.string('patTupleAHtoElecMu_AH115bbtautau.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# AH115 --> tau+ tau- --> 2l sample
#
fileNamesAH115tautau2l = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115_tautau_2l_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115_tautau_2l_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115_tautau_2l_3.root'
)
genPhaseSpaceCutAH115tautau2l = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameAH115tautau2l = cms.string('plotsAHtoElecMu_AH115tautau2l.root')
patTupleOutputFileNameAH115tautau2l = cms.untracked.string('patTupleAHtoElecMu_AH115tautau2l.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# AH115bb --> tau+ tau- --> 2l sample
#
fileNamesAH115bbtautau2l_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_2.root'
)
fileNamesAH115bbtautau2l_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_4.root'
)
fileNamesAH115bbtautau2l_part03 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_6.root'
)
fileNamesAH115bbtautau2l_part04 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_8.root'
)
fileNamesAH115bbtautau2l_part05 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_10.root'
)
fileNamesAH115bbtautau2l_part06 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_12.root'
)
fileNamesAH115bbtautau2l_part07 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_13.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_14.root'
)
fileNamesAH115bbtautau2l_part08 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_15.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_16.root'
)
fileNamesAH115bbtautau2l_part09 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_17.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_18.root'
)
fileNamesAH115bbtautau2l_part10 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_19.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_20.root'
)
fileNamesAH115bbtautau2l_part11 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_21.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_22.root'
)
fileNamesAH115bbtautau2l_part12 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_23.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_24.root'
)
fileNamesAH115bbtautau2l_part13 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_25.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_26.root'
)
fileNamesAH115bbtautau2l_part14 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_27.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_28.root'
)
fileNamesAH115bbtautau2l_part15 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_29.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_30.root'
)
fileNamesAH115bbtautau2l_part16 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_31.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_32.root'
)
fileNamesAH115bbtautau2l_part17 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_33.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_34.root'
)
fileNamesAH115bbtautau2l_part18 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_35.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_36.root'
)
fileNamesAH115bbtautau2l_part19 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_37.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_38.root'
)
fileNamesAH115bbtautau2l_part20 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_39.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_40.root'
)
fileNamesAH115bbtautau2l_part21 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_41.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_42.root'
)
fileNamesAH115bbtautau2l_part22 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_43.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH115bb_tautau_2l_44.root'
)
genPhaseSpaceCutAH115bbtautau2l = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameAH115bbtautau2l = cms.string('plotsAHtoElecMu_AH115bbtautau2l_partXX.root')
patTupleOutputFileNameAH115bbtautau2l = cms.untracked.string('patTupleAHtoElecMu_AH115bbtautau2l_partXX.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# AH160 --> tau+ tau- sample
#
fileNamesAH160tautau = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160v2_tautau_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160v2_tautau_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160v2_tautau_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160v2_tautau_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160v2_tautau_5.root'
)
genPhaseSpaceCutAH160tautau = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameAH160tautau = cms.string('plotsAHtoElecMu_AH160tautau.root')
patTupleOutputFileNameAH160tautau = cms.untracked.string('patTupleAHtoElecMu_AH160tautau.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# AH160bb --> tau+ tau- sample
#
fileNamesAH160bbtautau = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_5.root'
)
genPhaseSpaceCutAH160bbtautau = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameAH160bbtautau = cms.string('plotsAHtoElecMu_AH160bbtautau.root')
patTupleOutputFileNameAH160bbtautau = cms.untracked.string('patTupleAHtoElecMu_AH160bbtautau.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# AH160 --> tau+ tau- --> 2l sample
#
fileNamesAH160tautau2l = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160_tautau_2l_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160_tautau_2l_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160_tautau_2l_3.root'
)
genPhaseSpaceCutAH160tautau2l = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameAH160tautau2l = cms.string('plotsAHtoElecMu_AH160tautau2l.root')
patTupleOutputFileNameAH160tautau2l = cms.untracked.string('patTupleAHtoElecMu_AH160tautau2l.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# AH160bb --> tau+ tau- --> 2l sample
#
fileNamesAH160bbtautau2l_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_2.root'
)
fileNamesAH160bbtautau2l_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_4.root'
)
fileNamesAH160bbtautau2l_part03 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_6.root'
)
fileNamesAH160bbtautau2l_part04 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_8.root'
)
fileNamesAH160bbtautau2l_part05 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_10.root'
)
fileNamesAH160bbtautau2l_part06 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_12.root'
)
fileNamesAH160bbtautau2l_part07 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_13.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_14.root'
)
fileNamesAH160bbtautau2l_part08 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_15.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_16.root'
)
fileNamesAH160bbtautau2l_part09 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_17.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_18.root'
)
fileNamesAH160bbtautau2l_part10 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_19.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_20.root'
)
fileNamesAH160bbtautau2l_part11 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_21.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_22.root'
)
fileNamesAH160bbtautau2l_part12 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_23.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_24.root'
)
fileNamesAH160bbtautau2l_part13 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_25.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_26.root'
)
fileNamesAH160bbtautau2l_part14 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_27.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_28.root'
)
fileNamesAH160bbtautau2l_part15 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_29.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_30.root'
)
fileNamesAH160bbtautau2l_part16 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_31.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_32.root'
)
fileNamesAH160bbtautau2l_part17 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_33.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_34.root'
)
fileNamesAH160bbtautau2l_part18 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_35.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_36.root'
)
fileNamesAH160bbtautau2l_part19 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_37.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_38.root'
)
fileNamesAH160bbtautau2l_part20 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_39.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_40.root'
)
fileNamesAH160bbtautau2l_part21 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_41.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_42.root'
)
fileNamesAH160bbtautau2l_part22 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_43.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_44.root'
)
fileNamesAH160bbtautau2l_part23 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_45.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_46.root'
)
fileNamesAH160bbtautau2l_part24 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_47.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_48.root'
)
fileNamesAH160bbtautau2l_part25 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_49.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_50.root'
)
fileNamesAH160bbtautau2l_part26 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_51.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_AH160bb_tautau_2l_52.root'
)
genPhaseSpaceCutAH160bbtautau2l = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameAH160bbtautau2l = cms.string('plotsAHtoElecMu_AH160bbtautau2l_partXX.root')
patTupleOutputFileNameAH160bbtautau2l = cms.untracked.string('patTupleAHtoElecMu_AH160bbtautau2l_partXX.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# V+QQ Fall08 (Full Sim) MadGraph dataset, 1M; V=Z->ll, Mll>50 GeV, W->lnu, l=e, mu, tau, Q=c,b, pTb,c > 10 GeV for at least one b,c;
#
fileNamesVQQ = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_6.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_8.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_10.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_12.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_13.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_14.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_15.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_16.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_17.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_18.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_19.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_20.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_VQQ_21.root'
)
genPhaseSpaceCutVQQ = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameVQQ = cms.string('plotsAHtoElecMu_VQQ.root')
patTupleOutputFileNameVQQ = cms.untracked.string('patTupleAHtoElecMu_VQQ.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# PYTHIA WW inclusive
#
fileNamesWW = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_WW_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_WW_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_WW_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_WW_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_WW_5.root'
)
genPhaseSpaceCutWW = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameWW = cms.string('plotsAHtoElecMu_WW.root')
patTupleOutputFileNameWW = cms.untracked.string('patTupleAHtoElecMu_WW.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Wt Summer08 MadGraph dataset (inclusive, no any selections)
#
fileNamesTW = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TW_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TW_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TW_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TW_4.root'
)
genPhaseSpaceCutTW = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameTW = cms.string('plotsAHtoElecMu_TW.root')
patTupleOutputFileNameTW = cms.untracked.string('patTupleAHtoElecMu_TW.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# tt+jets Summer08 MadGraph dataset, 1M; tt->WbWb->all
#
fileNamesTTJets_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_2.root'
)
fileNamesTTJets_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_4.root'
)
fileNamesTTJets_part03 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_6.root'
)
fileNamesTTJets_part04 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_8.root'
)
fileNamesTTJets_part05 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_10.root'
)
fileNamesTTJets_part06 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_12.root'
)
fileNamesTTJets_part07 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_13.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_14.root'
)
fileNamesTTJets_part08 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_15.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_16.root'
)
fileNamesTTJets_part09 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_17.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_18.root'
)
fileNamesTTJets_part10 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_19.root'
#    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_20.root'
)
fileNamesTTJets_part11 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_21.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_22.root'
)
fileNamesTTJets_part12 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_23.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_24.root'
)
fileNamesTTJets_part13 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_25.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_26.root'
)
fileNamesTTJets_part14 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_27.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_28.root'
)
fileNamesTTJets_part15 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_29.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_30.root'
)
fileNamesTTJets_part16 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_31.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_32.root'
)
fileNamesTTJets_part17 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_33.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_34.root'
)
fileNamesTTJets_part18 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_35.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_36.root'
)
fileNamesTTJets_part19 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_37.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_38.root'
)
fileNamesTTJets_part20 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_39.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_40.root'
)
fileNamesTTJets_part21 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_41.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_42.root'
)
fileNamesTTJets_part22 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_43.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_TTJets_44.root'
)
genPhaseSpaceCutTTJets = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameTTJets = cms.string('plotsAHtoElecMu_TTJets_partXX.root')
patTupleOutputFileNameTTJets = cms.untracked.string('patTupleAHtoElecMu_TTJets_partXX.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Z + jets sample
#
fileNamesZplusJets = cms.untracked.vstring(
#    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_1.root',
#    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_2.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_3.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_4.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_5.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_6.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_7.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_8.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_9.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_10.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_11.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_12.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_13.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_14.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_15.root',
    'rfio:/castor/cern.ch/user/s/sunil/SkimFabruary09/test2/ZJets-madgraph/elecMuSkim_16.root'
)
genPhaseSpaceCutZplusJets = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
)
plotsOutputFileNameZplusJets = cms.string('plotsAHtoElecMu_ZplusJets.root')
patTupleOutputFileNameZplusJets = cms.untracked.string('patTupleAHtoElecMu_ZplusJets.root')
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# gamma+jets, pT > 15 GeV
#
fileNamesPhotonJet15 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_PhotonJet15_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_PhotonJet15_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_PhotonJet15_3.root',
#    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_PhotonJet15_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_PhotonJet15_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_PhotonJet15_6.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_PhotonJet15_7.root'
)
genPhaseSpaceCutPhotonJet15 = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
#  cut = cms.string('ptHat > 15.')
)
plotsOutputFileNamePhotonJet15 = cms.string('plotsAHtoElecMu_PhotonJet15.root')
patTupleOutputFileNamePhotonJet15 = cms.untracked.string('patTupleAHtoElecMu_PhotonJet15.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# gamma+jets, pT > 30 GeV
#
fileNamesPhotonJet30 = cms.untracked.vstring(
)
genPhaseSpaceCutPhotonJet30 = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('')
#  cut = cms.string('ptHat > 30.')
)
plotsOutputFileNamePhotonJet30 = cms.string('plotsAHtoElecMu_PhotonJet30.root')
patTupleOutputFileNamePhotonJet30 = cms.untracked.string('patTupleAHtoElecMu_PhotonJet30.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# e.m. QCD enrihed PYTHIA; EM enriched, pthat 20-30 GeV
#
fileNamesQCDem20to30_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_6.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_8.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_10.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_12.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_13.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_14.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_15.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_16.root'
)    
fileNamesQCDem20to30_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_17.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_18.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_19.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_20.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_21.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_22.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_23.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_24.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_25.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_26.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_27.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_28.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_29.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_30.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_31.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_32.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDem20to30_33.root'
)
genPhaseSpaceCutQCDem20to30 = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('ptHat < 30.')
)
plotsOutputFileNameQCDem20to30 = cms.string('plotsAHtoElecMu_QCDem20to30_partXX.root')
patTupleOutputFileNameQCDem20to30 = cms.untracked.string('patTupleAHtoElecMu_QCDem20to30_partXX.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# e.m. QCD enrihed PYTHIA; EM enriched, pthat 30-80 GeV
#
fileNamesQCDem30to80_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_6.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_8.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_12.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_14.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_15.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_17.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_19.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_20.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_21.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_22.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_24.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_25.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_26.root'
)
fileNamesQCDem30to80_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_27.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_28.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_29.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_30.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_31.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_32.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_33.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_34.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_35.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_36.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_38.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_40.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_41.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_43.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_44.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_46.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_48.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_50.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_52.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_53.root'
)
fileNamesQCDem30to80_part03 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_54.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_56.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_57.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_59.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_60.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_61.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_62.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_64.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_65.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_66.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_67.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_68.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_70.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_72.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_73.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_74.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_75.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_76.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_78.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_79.root'
)
fileNamesQCDem30to80_part04 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_80.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_83.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_84.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_85.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_87.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_88.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_89.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_90.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_91.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_92.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_94.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_95.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_96.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_97.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_98.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_99.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_100.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_102.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_103.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_104.root'
)
fileNamesQCDem30to80_part05 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_105.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_106.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_109.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_111.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_112.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_113.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_116.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_117.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_118.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_119.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_120.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_121.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_123.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_126.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_128.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_129.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_132.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_133.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_134.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_135.root'
)
fileNamesQCDem30to80_part06 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_136.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_138.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_142.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_143.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_145.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_146.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_147.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_148.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_149.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_150.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_151.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_152.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_153.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_156.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_157.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_158.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_160.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_161.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_162.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_163.root'
)
fileNamesQCDem30to80_part07 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_164.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_165.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_169.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_170.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_171.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_175.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_176.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_177.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_178.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_179.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_180.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_181.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_182.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_183.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_184.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_186.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_187.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_190.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_191.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_192.root'
)
fileNamesQCDem30to80_part08 = cms.untracked.vstring(   
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_194.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_196.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_197.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_200.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_201.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_203.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_206.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_207.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_208.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_209.root'
)
fileNamesQCDem30to80_part09 = cms.untracked.vstring(   
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_222.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_224.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_225.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_227.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_228.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_229.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_231.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_233.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_234.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_235.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_238.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_239.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_240.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_241.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_242.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_243.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_245.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_246.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_247.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_248.root'
)
fileNamesQCDem30to80_part10 = cms.untracked.vstring(   
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_249.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_250.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_251.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_252.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_253.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_256.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_257.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_258.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_259.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_260.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_261.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_263.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_264.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_265.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_266.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_269.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_272.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_273.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_274.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_275.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_276.root'
)
fileNamesQCDem30to80_part11 = cms.untracked.vstring(   
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_3.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_10.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_13.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_16.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_18.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_23.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_37.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_39.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_42.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_45.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_47.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_49.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_51.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_55.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_58.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_63.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_69.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_71.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_77.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_81.root'
)
fileNamesQCDem30to80_part12 = cms.untracked.vstring(   
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_82.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_86.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_93.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_101.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_107.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_108.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_110.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_114.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_115.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_122.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_124.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_125.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_127.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_130.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_131.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_137.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_139.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_140.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_141.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_144.root'
)
fileNamesQCDem30to80_part13 = cms.untracked.vstring(   
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_154.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_155.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_159.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_166.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_167.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_168.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_172.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_173.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_174.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_185.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_188.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_189.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_193.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_195.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_198.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_199.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_202.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_204.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_205.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_214.root'
)
fileNamesQCDem30to80_part14 = cms.untracked.vstring(   
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_216.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_223.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_226.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_230.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_232.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_236.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_237.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_244.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_254.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_255.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_262.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_267.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_268.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_270.root',
   'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_271.root'
)
fileNamesQCDem30to80_part15 = cms.untracked.vstring(   
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_210.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_211.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_212.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_213.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_215.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_217.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_218.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_219.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_220.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem30to80_221.root'
)
genPhaseSpaceCutQCDem30to80 = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('ptHat > 30. & ptHat < 80. ')
)
plotsOutputFileNameQCDem30to80 = cms.string('plotsAHtoElecMu_QCDem30to80_partXX.root')
patTupleOutputFileNameQCDem30to80 = cms.untracked.string('patTupleAHtoElecMu_QCDem30to80_partXX.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# e.m. QCD enrihed PYTHIA; EM enriched, pthat 80-170 GeV,
#
fileNamesQCDem80to170_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_6.root'
)
fileNamesQCDem80to170_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_39.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_24.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_25.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_27.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_30.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_31.root'
)
fileNamesQCDem80to170_part03 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_13.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_19.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_26.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_28.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_32.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_17.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_18.root'
)
fileNamesQCDem80to170_part04 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_38.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_40.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_41.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_42.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_43.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_44.root'
    )
fileNamesQCDem80to170_part05 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_29.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_33.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_34.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_35.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_36.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_37.root'
    )
fileNamesQCDem80to170_part06 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_8.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_10.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_12.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_15.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_14.root'
    )
fileNamesQCDem80to170_part07 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_45.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_16.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_20.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_21.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_22.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_newQCDem80to170_23.root'
    )
genPhaseSpaceCutQCDem80to170 = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('ptHat > 80.')
)
plotsOutputFileNameQCDem80to170 = cms.string('plotsAHtoElecMu_QCDem80to170_partXX.root')
patTupleOutputFileNameQCDem80to170 = cms.untracked.string('patTupleAHtoElecMu_QCDem80to170_partXX.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# e.m. QCD enrihed PYTHIA; BC to E, pthat 20-30 GeV
#
fileNamesQCDbc20to30 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_6.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_8.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_10.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_12.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc20to30_13.root'
)    
genPhaseSpaceCutQCDbc20to30 = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('ptHat < 30.')
)
plotsOutputFileNameQCDbc20to30 = cms.string('plotsAHtoElecMu_QCDbc20to30.root')
patTupleOutputFileNameQCDbc20to30 = cms.untracked.string('patTupleAHtoElecMu_QCDbc20to30.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# e.m. QCD enrihed PYTHIA; EM enriched, pthat 30-80 GeV
#
fileNamesQCDbc30to80 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_6.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_8.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_10.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_12.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_13.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc30to80_14.root'
)
genPhaseSpaceCutQCDbc30to80 = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('ptHat > 30. & ptHat < 80. ')
)
plotsOutputFileNameQCDbc30to80 = cms.string('plotsAHtoElecMu_QCDbc30to80.root')
patTupleOutputFileNameQCDbc30to80 = cms.untracked.string('patTupleAHtoElecMu_QCDbc30to80.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# e.m. QCD enrihed PYTHIA; EM enriched, pthat 80-170 GeV,
#
fileNamesQCDbc80to170_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_1.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_2.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_3.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_4.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_5.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_6.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_7.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_8.root'
)
fileNamesQCDbc80to170_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_9.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_10.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_11.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_12.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_13.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_14.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_15.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_16.root',
    'rfio:/castor/cern.ch/user/c/cerati/SkimmingElecMu/elecMuSkim_QCDbc80to170_17.root'
)
genPhaseSpaceCutQCDbc80to170 = cms.PSet(
  pluginName = cms.string('genPhaseSpaceCut'),
  pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
  src = cms.InputTag('genPhaseSpaceEventInfo'),
  cut = cms.string('ptHat > 80.')
)
plotsOutputFileNameQCDbc80to170 = cms.string('plotsAHtoElecMu_QCDbc80to170_partXX.root')
patTupleOutputFileNameQCDbc80to170 = cms.untracked.string('patTupleAHtoElecMu_QCDbc80to170_partXX.root')
#--------------------------------------------------------------------------------
