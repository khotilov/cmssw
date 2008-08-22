import FWCore.ParameterSet.Config as cms

# $Id: RecoPFJets_cff.py,v 1.3 2008/04/21 03:27:29 rpw Exp $
#
# ShR 27 Mar 07: move modules producing candidates for Jets into separate cff file due to scheduling problem
#
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.FastjetParameters_cfi import *
from RecoJets.JetProducers.KtJetParameters_cfi import *
from RecoJets.JetProducers.SISConeJetParameters_cfi import *
from RecoJets.JetProducers.IconeJetParameters_cfi import *
kt4PFJets = cms.EDProducer("KtJetProducer",
    FastjetNoPU,
    KtJetParameters,
    PFJetParameters,
    alias = cms.untracked.string('KT4PFJet'),
    FJ_ktRParam = cms.double(0.4)
)

kt6PFJets = cms.EDProducer("KtJetProducer",
    FastjetNoPU,
    KtJetParameters,
    PFJetParameters,
    alias = cms.untracked.string('KT6PFJet'),
    FJ_ktRParam = cms.double(0.6)
)

iterativeCone5PFJets = cms.EDProducer("IterativeConeJetProducer",
    IconeJetParameters,
    PFJetParameters,
    alias = cms.untracked.string('IC5PFJet'),
    coneRadius = cms.double(0.5)
)

sisCone5PFJets = cms.EDProducer("SISConeJetProducer",
    PFJetParameters,
    SISConeJetParameters,
    FastjetNoPU,
    alias = cms.untracked.string('SISC5PFJet'),
    coneRadius = cms.double(0.5)
)

sisCone7PFJets = cms.EDProducer("SISConeJetProducer",
    PFJetParameters,
    SISConeJetParameters,
    FastjetNoPU,
    alias = cms.untracked.string('SISC7PFJet'),
    coneRadius = cms.double(0.7)
)

recoPFJets = cms.Sequence(kt4PFJets+kt6PFJets+iterativeCone5PFJets+sisCone5PFJets+sisCone7PFJets)

