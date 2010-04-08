import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------  
# produce collections of pat::Electrons passing selection criteria
#--------------------------------------------------------------------------------

# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string parser

# require electron candidate to pass the eidRobustTight electron id. criteria
selectedLayer1ElectronsTightId = cms.EDFilter("PATElectronSelector",
    cut = cms.string('(abs(superCluster.eta) < 1.479 & electronID("eidRobustTight") > 0 & eSuperClusterOverP < 1.05 & eSuperClusterOverP > 0.95) | (abs(superCluster.eta) > 1.479 & electronID("eidRobustTight") > 0 & eSuperClusterOverP < 1.12 & eSuperClusterOverP > 0.95)'),
    filter = cms.bool(False)
)

# require electron candidate to pass the eidRobustLoose electron id. criteria
selectedLayer1ElectronsLooseId = cms.EDFilter("PATElectronSelector",
    cut = cms.string('(abs(superCluster.eta) < 1.479 & electronID("eidRobustLoose") > 0 & eSuperClusterOverP < 1.4 & eSuperClusterOverP > 0.8) | (abs(superCluster.eta) > 1.479 & electronID("eidRobustLoose") > 0 & eSuperClusterOverP < 1.6 & eSuperClusterOverP > 0.8)'),
    filter = cms.bool(False)
)

# require electron candidate to not be within eta-crack
# between Barrel and Encap ECAL calorimeter
selectedLayer1ElectronsAntiCrackCut = cms.EDFilter("PATElectronSelector",
    cut = cms.string('abs(superCluster.eta) < 1.442 | abs(superCluster.eta) > 1.560'),
    filter = cms.bool(False)
)

# require electron candidate to be within geometric acceptance of electron trigger
selectedLayer1ElectronsEta21 = cms.EDFilter("PATElectronSelector",
    cut = cms.string('abs(eta) < 2.1'),
    filter = cms.bool(False)
)

# require electron candidate to have transverse momentum above threshold
selectedLayer1ElectronsPt15 = cms.EDFilter("PATElectronSelector",
    cut = cms.string('pt > 15.'),
    filter = cms.bool(False)
)

# require electron candidate to be isolated
# with respect to tracks (of Pt >~ 0.3 GeV)
selectedLayer1ElectronsTrkIso = cms.EDFilter("PATElectronSelector",
    cut = cms.string('userIsolation("pat::TrackIso") < 1.'),
    filter = cms.bool(False)
)                                    

# require electron candidate to be isolated
# with respect to energy deposits in ECAL
# (not associated to electron candidate)
selectedLayer1ElectronsEcalIso = cms.EDFilter("PATElectronSelector",
    cut = cms.string('(abs(superCluster.eta) < 1.479 & userIsolation("pat::EcalIso") < 2.5) | (abs(superCluster.eta) > 1.479 & userIsolation("pat::EcalIso") < 3.5)'),
    filter = cms.bool(False)
)

# require electron candidate to be linked to (GSF) track
selectedLayer1ElectronsTrk = cms.EDFilter("PATElectronSelector",
    cut = cms.string('gsfTrack.isNonnull'),
    filter = cms.bool(False)
)

# require track of electron candidate to have small transverse impact parameter
# (in order to veto electrons resulting from b-quark decays)
selectedLayer1ElectronsTrkIP = cms.EDFilter("PATElectronIpSelector",
    vertexSource = cms.InputTag("selectedPrimaryVertexPosition"),
    IpMax = cms.double(0.05),
    filter = cms.bool(False)                                               
)

#--------------------------------------------------------------------------------
# define additional collections of muon candidates
# with loose track and ECAL isolation applied
#
# (NOTE: to be used for the purpose of factorizing efficiencies
#        of muon isolation from other event selection criteria,
#        in order to avoid problems with limited Monte Carlo statistics)
#--------------------------------------------------------------------------------

selectedLayer1ElectronsTrkIsoLooseIsolation = copy.deepcopy(selectedLayer1ElectronsTrkIso)
selectedLayer1ElectronsTrkIsoLooseIsolation.cut = cms.string('userIsolation("pat::TrackIso") < 8.')

selectedLayer1ElectronsEcalIsoLooseIsolation = copy.deepcopy(selectedLayer1ElectronsEcalIso)
selectedLayer1ElectronsEcalIsoLooseIsolation.cut = cms.string('userIsolation("pat::EcalIso") < 8.')

selectedLayer1ElectronsTrkLooseIsolation = copy.deepcopy(selectedLayer1ElectronsTrk)

selectedLayer1ElectronsTrkIPlooseIsolation = copy.deepcopy(selectedLayer1ElectronsTrkIP)
