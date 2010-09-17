import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# produce all objects specific to Z --> tau-jet + tau-jet channel
# which do **not** get stored in PATTuple
#--------------------------------------------------------------------------------

#
# produce collections of pat::Electrons, pat::Muons and
# pat::(PF)Taus passing different selection criteria
#
from TauAnalysis.RecoTools.patLeptonSelection_cff import *
#
# produce collection of pat::Jets passing Et threshold and
# Eta acceptance cuts and not overlapping with any object
# passing selection criteria for pat::Electron, pat::Muon or pat::(PF)Tau
# (pat::Jet collection to be considered for central-jet veto)
#
from TauAnalysis.RecoTools.patJetSelection_cff import *
#
# produce collections of pat::(Calo)MET objects
# passing different selection criteria
#
from TauAnalysis.RecoTools.patMetSelection_cff import *
#
# produce collections of tau-jet + tau-jet pairs
# passing different selection criteria
#
from TauAnalysis.CandidateTools.diTauPairProduction_cff import *
from TauAnalysis.CandidateTools.diTauPairSelectionAllKinds_cff import *

producePatTupleZtoDiTauSpecific = cms.Sequence(
    selectPatMuons
   + selectPatElectrons
   + selectPatTaus + selectPatTausForDiTau
   + produceDiTauPairs
   + selectDiTauPairs + selectDiTauPairsLoose2ndTau
   + selectPatJets
)
