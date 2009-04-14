import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------  
# produce collections of electron + tau-jet pairs passing selection criteria
#--------------------------------------------------------------------------------

# require electron and tau-jet to be separated in eta-phi,
# in order to ensure that both do not refer to one and the same physical particle
# (NOTE: cut is already applied during skimming,
#        so should not reject any events)
selectedElecTauPairsAntiOverlapVeto = cms.EDFilter("PATElecTauPairSelector",
    cut = cms.string('dR12 > 0.7'),
    filter = cms.bool(False)
)

# require electron and tau to form a zero-charge pair
selectedElecTauPairsZeroCharge = cms.EDFilter("PATElecTauPairSelector",
    cut = cms.string('charge = 0'),
    #cut = cms.string('(leg1.charge + leg2.leadTrack.charge) = 0'), # NOTE: to be used for background studies only !!                    
    filter = cms.bool(False)
)

#require cut transverse mass of electron and MET
selectedElecTauPairsMt1MET = cms.EDFilter("PATElecTauPairSelector",
    cut = cms.string('mt1MET < 60.'),
    filter = cms.bool(False)
)
