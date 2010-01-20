import FWCore.ParameterSet.Config as cms

# module to produce jet correction factors associated in a valuemap
patJetCorrFactors = cms.EDProducer("JetCorrFactorsProducer",
     ## the use of emf in the JEC is not yet implemented
     useEMF     = cms.bool(False),
     ## input collection of jets
     jetSource  = cms.InputTag("ak5CaloJets"),
     ## set of correction factors
     corrSample = cms.string("Summer09"),
     ## correction levels
     corrLevels = cms.PSet(
       ## tags for the individual jet corrections; when
       ## not available the string should be set to 'none'    
       L1Offset   = cms.string('none'),
       L2Relative = cms.string('L2Relative_AK5Calo'),
       L3Absolute = cms.string('L3Absolute_AK5Calo'),
       L4EMF      = cms.string('none'),
       L5Flavor   = cms.string('L5Flavor_IC5'),       # to be changed to L5Flavor   = cms.string('L5Flavor_AK5'),
       L6UE       = cms.string('none'),
       L7Parton   = cms.string('L7Parton_SC5'),       # to be changed to L7Parton   = cms.string('L7Parton_AK5'),
     ),
     ## choose sample type for flavor dependend corrections:
     sampleType = cms.string('dijet')  ##  'dijet': from dijet sample
                                       ##  'ttbar': from ttbar sample

)
