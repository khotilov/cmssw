import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.svFitAlgorithm_cfi import *

svFitLikelihoodDiTauPairKinematicsPhaseSpace = copy.deepcopy(svFitLikelihoodDiTauKinematicsPhaseSpace)
svFitLikelihoodDiTauPairKinematicsPhaseSpace.pluginType = "SVfitLikelihoodDiTauPairKinematics"
svFitLikelihoodDiTauPairKinematicsPhaseSpace.leg1.pluginType = "SVfitTauLikelihoodPhaseSpace"
svFitLikelihoodDiTauPairKinematicsPhaseSpace.leg2.pluginType = "SVfitTauLikelihoodPhaseSpace"

svFitLikelihoodDiTauPairKinematicsPolarized = copy.deepcopy(svFitLikelihoodDiTauKinematicsPolarized)
svFitLikelihoodDiTauPairKinematicsPolarized.pluginType = "SVfitLikelihoodDiTauPairKinematics"
svFitLikelihoodDiTauPairKinematicsPolarized.leg2.pluginType = "SVfitTauLikelihoodPolarization"
svFitLikelihoodDiTauPairKinematicsPolarized.leg1 = svFitLikelihoodDiTauPairKinematicsPolarized.leg2

svFitLikelihoodDiTauPairMEt = copy.deepcopy(svFitLikelihoodMEt)
svFitLikelihoodDiTauPairMEt.pluginType = cms.string("SVfitLikelihoodMEtDiTau")

svFitLikelihoodDiTauPairPtBalance = copy.deepcopy(svFitLikelihoodDiTauPtBalance)
svFitLikelihoodDiTauPairPtBalance.pluginType = cms.string("SVfitLikelihoodDiTauPairPtBalance")

#--------------------------------------------------------------------------------
# produce combinations of tau-jet + tau-jet pairs
#--------------------------------------------------------------------------------

allDiTauPairs = cms.EDProducer("PATDiTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('cleanPatTaus'),
    srcLeg2 = cms.InputTag('cleanPatTaus'),
    dRmin12 = cms.double(0.3),
    srcMET = cms.InputTag('patMETs'),
    srcPrimaryVertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),
    srcGenParticles = cms.InputTag('genParticles'),                  
    recoMode = cms.string(""),
    doSVreco = cms.bool(False),                          
    SVOptions = cms.PSet(
        usePtBalanceInFit = cms.bool(True),
        useMEtInFit = cms.bool(True),
        useLeg1TrackingInFit = cms.bool(False),
        useLeg2TrackingInFit = cms.bool(False),
        correctPrimaryVertexInFit = cms.bool(False)
    ),                           
    scaleFuncImprovedCollinearApprox = cms.string('1'),                           
    verbosity = cms.untracked.int32(0)
)

#---------------------------------------------------------------------------------------
# produce tau-jet + tau-jet pairs with tau-jet selection criteria applied on first leg
#---------------------------------------------------------------------------------------

# require 1st tau candidate to be within geometric acceptance of Pixel + SiTracker detectors
selectedDiTauPairs1stTauEta21Individual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTauEta21Individual')
)

selectedDiTauPairs1stTauEta21Cumulative = selectedDiTauPairs1stTauEta21Individual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTauEta21Cumulative')
)

# require 1st tau candidate to have transverse energy above threshold
selectedDiTauPairs1stTauPt20Individual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTauPt20Individual')
)

selectedDiTauPairs1stTauPt20Cumulative = selectedDiTauPairs1stTauPt20Individual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTauPt20Cumulative')
)

# require 1st tau candidate to have a leading track
# (track of Pt > 1. GeV within matching cone of size dR = 0.2 around jet-axis)
selectedDiTauPairs1stTauLeadTrkIndividual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTauLeadTrkIndividual')
)

selectedDiTauPairs1stTauLeadTrkCumulative = selectedDiTauPairs1stTauLeadTrkIndividual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTauLeadTrkCumulative')
)

# require leading track of 1st tau candidate to have Pt > 12. GeV
selectedDiTauPairs1stTauLeadTrkPtIndividual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stLeadTrkPtIndividual')
)

selectedDiTauPairs1stTauLeadTrkPtCumulative = selectedDiTauPairs1stTauLeadTrkPtIndividual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stLeadTrkPtCumulative')
)

# require 1st tau candidate to pass TaNC discriminator for 0.25% nominal fake-rate
selectedDiTauPairs1stTauTaNCdiscrIndividual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stTaNCdiscrIndividual')
)

selectedDiTauPairs1stTauTaNCdiscrCumulative = selectedDiTauPairs1stTauTaNCdiscrIndividual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stTaNCdiscrCumulative')
)

# require 1st tau candidate to have no tracks of Pt > 1. GeV
# in isolation cone of size dR = 0.8, surrounding signal cone of size dR = 5./Et
selectedDiTauPairs1stTauTrkIsoIndividual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stTrkIsoIndividual')
)

selectedDiTauPairs1stTauTrkIsoCumulative = selectedDiTauPairs1stTauTrkIsoIndividual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stTrkIsoCumulative')
)

# require 1st tau candidate to be isolated
# with respect to energy deposits in ECAL
selectedDiTauPairs1stTauEcalIsoIndividual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stEcalIsoIndividual')
)

selectedDiTauPairs1stTauEcalIsoCumulative = selectedDiTauPairs1stTauEcalIsoIndividual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stEcalIsoCumulative')
)

# require 1st tau candidate to have either one or three tracks within signal cone
selectedDiTauPairs1stTauProngIndividual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stProngIndividual')
)

selectedDiTauPairs1stTauProngCumulative = selectedDiTauPairs1stTauProngIndividual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stProngCumulative')
)

# require 1st tau candidate to have charge either +1 or -1
# (computed as sum of charges of tracks within signal cone)
selectedDiTauPairs1stTauChargeIndividual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stChargeIndividual')
)

selectedDiTauPairs1stTauChargeCumulative = selectedDiTauPairs1stTauChargeIndividual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stChargeCumulative')
)

# require 1st tau candidate to pass muon veto
selectedDiTauPairs1stTauMuonVetoIndividual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stMuonVetoIndividual')
)

selectedDiTauPairs1stTauMuonVetoCumulative = selectedDiTauPairs1stTauMuonVetoIndividual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stMuonVetoCumulative')
)

# require 1st tau candidate to pass electron veto
selectedDiTauPairs1stTauElectronVetoIndividual = allDiTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stElectronVetoIndividual')
)

selectedDiTauPairs1stTauElectronVetoCumulative = selectedDiTauPairs1stTauElectronVetoIndividual.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau1stElectronVetoCumulative')
)

#----------------------------------------------------------------------------------------
# produce tau-jet + tau-jet pairs with tau-jet selection criteria applied on second leg
#----------------------------------------------------------------------------------------

# require 2nd tau candidate to be within geometric acceptance of Pixel + SiTracker detectors
selectedDiTauPairs2ndTauEta21Individual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTauEta21Individual')
)    

selectedDiTauPairs2ndTauEta21Cumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTauEta21Cumulative')
)    

# require 2nd tau candidate to have transverse energy above threshold
selectedDiTauPairs2ndTauPt20Individual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTauPt20Individual')
)    

selectedDiTauPairs2ndTauPt20Cumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTauPt20Cumulative')
)    

# require 2nd tau candidate to have a leading track
# (track of Pt > 1. GeV within matching cone of size dR = 0.2 around jet-axis)
selectedDiTauPairs2ndTauLeadTrkIndividual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTauLeadTrkIndividual')
)    

selectedDiTauPairs2ndTauLeadTrkCumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTauLeadTrkCumulative')
)    

# require leading track of 2nd tau candidate to have Pt > 8. GeV
selectedDiTauPairs2ndTauLeadTrkPtIndividual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndLeadTrkPtIndividual'),
    
)    

selectedDiTauPairs2ndTauLeadTrkPtCumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndLeadTrkPtCumulative')
)                           

# require 2nd tau candidate to pass TaNC discriminator for 0.25% nominal fake-rate
selectedDiTauPairs2ndTauTaNCdiscrIndividual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
        srcLeg1 = cms.InputTag('selectedPatTausForDiTau2ndTaNCdiscrIndividual')
)

selectedDiTauPairs2ndTauTaNCdiscrCumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau2ndTaNCdiscrCumulative')
)

# require 2nd tau candidate to have no tracks of Pt > 1. GeV
# in isolation cone of size dR = 0.8, surrounding signal cone of size dR = 5./Et
selectedDiTauPairs2ndTauTrkIsoIndividual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndTrkIsoIndividual')
)    

selectedDiTauPairs2ndTauTrkIsoCumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndTrkIsoCumulative')
)    

# require 2nd tau candidate to be isolated
# with respect to energy deposits in ECAL
selectedDiTauPairs2ndTauEcalIsoIndividual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndEcalIsoIndividual')
)    

selectedDiTauPairs2ndTauEcalIsoCumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndEcalIsoCumulative')
)    

# require 2nd tau candidate to have either one or three tracks within signal cone
selectedDiTauPairs2ndTauProngIndividual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndProngIndividual')
)    

selectedDiTauPairs2ndTauProngCumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndProngCumulative')
)    

# require 2nd tau candidate to have charge either +1 or -1
# (computed as sum of charges of tracks within signal cone)
selectedDiTauPairs2ndTauChargeIndividual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndChargeIndividual')
)    

selectedDiTauPairs2ndTauChargeCumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg2 = cms.InputTag('selectedPatTausForDiTau2ndChargeCumulative')
)

# require 2nd tau candidate to pass muon veto
selectedDiTauPairs2ndTauMuonVetoIndividual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau2ndMuonVetoIndividual')
)

selectedDiTauPairs2ndTauMuonVetoCumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau2ndMuonVetoCumulative')
)

# require 2nd tau candidate to pass electron veto
selectedDiTauPairs2ndTauElectronVetoIndividual = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau2ndElectronVetoIndividual')
)

selectedDiTauPairs2ndTauElectronVetoCumulative = selectedDiTauPairs1stTauElectronVetoCumulative.clone(
    srcLeg1 = cms.InputTag('selectedPatTausForDiTau2ndElectronVetoCumulative')
)

selectedDiTauPairs2ndTauElectronVetoCumulative.svFit = cms.PSet(
    psKine = cms.PSet(
        likelihoodFunctions = cms.VPSet(
            svFitLikelihoodDiTauPairKinematicsPhaseSpace         
        ),
        estUncertainties = cms.PSet(
            numSamplings = cms.int32(-1)
        )
    ),
    psKine_MEt = cms.PSet(
        likelihoodFunctions = cms.VPSet(
            svFitLikelihoodDiTauPairKinematicsPhaseSpace,
            svFitLikelihoodDiTauPairMEt
        ),
        estUncertainties = cms.PSet(
            numSamplings = cms.int32(-1)
        )
    ),
    psKine_MEt_ptBalance = cms.PSet(
        likelihoodFunctions = cms.VPSet(
            svFitLikelihoodDiTauPairKinematicsPhaseSpace,
            svFitLikelihoodDiTauPairMEt,
            svFitLikelihoodDiTauPairPtBalance
    ),
        estUncertainties = cms.PSet(
            #numSamplings = cms.int32(1000)
            numSamplings = cms.int32(-1)
        )
    ),
    polKine = cms.PSet(
        likelihoodFunctions = cms.VPSet(
            svFitLikelihoodDiTauPairKinematicsPolarized        
        ),
        estUncertainties = cms.PSet(
            numSamplings = cms.int32(-1)
        )
    ),
    polKine_MEt = cms.PSet(
        likelihoodFunctions = cms.VPSet(
            svFitLikelihoodDiTauPairKinematicsPolarized,
            svFitLikelihoodDiTauPairMEt
        ),
        estUncertainties = cms.PSet(
            numSamplings = cms.int32(-1)
        )
    ),
    polKine_MEt_ptBalance = cms.PSet(
        likelihoodFunctions = cms.VPSet(
            svFitLikelihoodDiTauPairKinematicsPolarized,
            svFitLikelihoodDiTauPairMEt,
            svFitLikelihoodDiTauPairPtBalance
        ),
        estUncertainties = cms.PSet(
            #numSamplings = cms.int32(1000)
            numSamplings = cms.int32(-1)
        )
    )
)

produceDiTauPairs = cms.Sequence(
    allDiTauPairs
   * selectedDiTauPairs1stTauEta21Individual * selectedDiTauPairs1stTauEta21Cumulative
   * selectedDiTauPairs1stTauPt20Individual * selectedDiTauPairs1stTauPt20Cumulative
   * selectedDiTauPairs1stTauLeadTrkIndividual * selectedDiTauPairs1stTauLeadTrkCumulative
   * selectedDiTauPairs1stTauLeadTrkPtIndividual * selectedDiTauPairs1stTauLeadTrkPtCumulative
   * selectedDiTauPairs1stTauTaNCdiscrIndividual * selectedDiTauPairs1stTauTaNCdiscrCumulative 
   * selectedDiTauPairs1stTauTrkIsoIndividual * selectedDiTauPairs1stTauTrkIsoCumulative
   * selectedDiTauPairs1stTauEcalIsoIndividual * selectedDiTauPairs1stTauEcalIsoCumulative
   * selectedDiTauPairs1stTauProngIndividual * selectedDiTauPairs1stTauProngCumulative
   * selectedDiTauPairs1stTauChargeIndividual * selectedDiTauPairs1stTauChargeCumulative
   * selectedDiTauPairs1stTauMuonVetoIndividual * selectedDiTauPairs1stTauMuonVetoCumulative
   * selectedDiTauPairs1stTauElectronVetoIndividual * selectedDiTauPairs1stTauElectronVetoCumulative
   * selectedDiTauPairs2ndTauEta21Individual * selectedDiTauPairs2ndTauEta21Cumulative
   * selectedDiTauPairs2ndTauPt20Individual * selectedDiTauPairs2ndTauPt20Cumulative
   * selectedDiTauPairs2ndTauLeadTrkIndividual * selectedDiTauPairs2ndTauLeadTrkCumulative
   * selectedDiTauPairs2ndTauLeadTrkPtIndividual * selectedDiTauPairs2ndTauLeadTrkPtCumulative
   * selectedDiTauPairs2ndTauTaNCdiscrIndividual * selectedDiTauPairs2ndTauTaNCdiscrCumulative  
   * selectedDiTauPairs2ndTauTrkIsoIndividual * selectedDiTauPairs2ndTauTrkIsoCumulative
   * selectedDiTauPairs2ndTauEcalIsoIndividual * selectedDiTauPairs2ndTauEcalIsoCumulative
   * selectedDiTauPairs2ndTauProngIndividual * selectedDiTauPairs2ndTauProngCumulative
   * selectedDiTauPairs2ndTauChargeIndividual * selectedDiTauPairs2ndTauChargeCumulative
   * selectedDiTauPairs2ndTauMuonVetoIndividual * selectedDiTauPairs2ndTauMuonVetoCumulative
   * selectedDiTauPairs2ndTauElectronVetoIndividual * selectedDiTauPairs2ndTauElectronVetoCumulative
)
