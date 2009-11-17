import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *


# extra l1 producers for kt jets
kt5CaloJetsL1 = kt4CaloJetsL1.clone( src = 'kt5CaloJets' )
kt7CaloJetsL1 = kt6CaloJetsL1.clone( src = 'kt7CaloJets' )
kt5PFJetsL1   = kt4PFJetsL1.clone( src  = 'kt5PFJets' )
kt7PFJetsL1   = kt6PFJetsL1.clone( src  = 'kt7PFJets' )

# extra l1 producers for ca jets
ca4CaloJetsL1 = kt4CaloJetsL1.clone( src = 'ca4CaloJets' )
ca5CaloJetsL1 = ca4CaloJetsL1.clone( src = 'ca5CaloJets' )
ca6CaloJetsL1 = ca4CaloJetsL1.clone( src = 'ca6CaloJets' )
ca7CaloJetsL1 = ca6CaloJetsL1.clone( src = 'ca7CaloJets' )

ca4PFJetsL1   = kt4PFJetsL1.clone( src = 'ca4PFJets' )
ca5PFJetsL1   = ca4PFJetsL1.clone( src = 'ca5PFJets' )
ca6PFJetsL1   = ca4PFJetsL1.clone( src = 'ca6PFJets' )
ca7PFJetsL1   = ca6PFJetsL1.clone( src = 'ca7PFJets' )


# extra l2l3 producers for kt jets
kt5CaloJetsL2L3 = kt4CaloJetsL2L3.clone( src = 'kt5CaloJets' )
kt7CaloJetsL2L3 = kt6CaloJetsL2L3.clone( src = 'kt7CaloJets' )
kt5PFJetsL2L3   = kt4PFJetsL2L3.clone( src  = 'kt5PFJets' )
kt7PFJetsL2L3   = kt6PFJetsL2L3.clone( src  = 'kt7PFJets' )

# extra l2l3 producers for ca jets
ca4CaloJetsL2L3 = kt4CaloJetsL2L3.clone( src = 'ca4CaloJets' )
ca5CaloJetsL2L3 = ca4CaloJetsL2L3.clone( src = 'ca5CaloJets' )
ca6CaloJetsL2L3 = ca4CaloJetsL2L3.clone( src = 'ca6CaloJets' )
ca7CaloJetsL2L3 = ca6CaloJetsL2L3.clone( src = 'ca7CaloJets' )

ca4PFJetsL2L3   = kt4PFJetsL2L3.clone( src = 'ca4PFJets' )
ca5PFJetsL2L3   = ca4PFJetsL2L3.clone( src = 'ca5PFJets' )
ca6PFJetsL2L3   = ca4PFJetsL2L3.clone( src = 'ca6PFJets' )
ca7PFJetsL2L3   = ca6PFJetsL2L3.clone( src = 'ca7PFJets' )


# extra l1l2l3 producers for kt jets
kt5CaloJetsL1L2L3 = kt4CaloJetsL1L2L3.clone( src = 'kt5CaloJets' )
kt7CaloJetsL1L2L3 = kt6CaloJetsL1L2L3.clone( src = 'kt7CaloJets' )
kt5PFJetsL1L2L3   = kt4PFJetsL1L2L3.clone( src  = 'kt5PFJets' )
kt7PFJetsL1L2L3   = kt6PFJetsL1L2L3.clone( src  = 'kt7PFJets' )

# extra l1l2l3 producers for ca jets
ca4CaloJetsL1L2L3 = kt4CaloJetsL1L2L3.clone( src = 'ca4CaloJets' )
ca5CaloJetsL1L2L3 = ca4CaloJetsL1L2L3.clone( src = 'ca5CaloJets' )
ca6CaloJetsL1L2L3 = ca4CaloJetsL1L2L3.clone( src = 'ca6CaloJets' )
ca7CaloJetsL1L2L3 = ca6CaloJetsL1L2L3.clone( src = 'ca7CaloJets' )

ca4PFJetsL1L2L3   = kt4PFJetsL1L2L3.clone( src = 'ca4PFJets' )
ca5PFJetsL1L2L3   = ca4PFJetsL1L2L3.clone( src = 'ca5PFJets' )
ca6PFJetsL1L2L3   = ca4PFJetsL1L2L3.clone( src = 'ca6PFJets' )
ca7PFJetsL1L2L3   = ca6PFJetsL1L2L3.clone( src = 'ca7PFJets' )
