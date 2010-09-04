import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.tools.objProdConfigurator import *
from TauAnalysis.CandidateTools.resolutions_cfi import *
from TauAnalysis.CandidateTools.svFitAlgorithm_cfi import *

#--------------------------------------------------------------------------------
# produce combinations of muon + tau-jet pairs
#--------------------------------------------------------------------------------

allMuTauPairs = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    srcLeg2 = cms.InputTag('selectedPatTausForMuTauElectronVetoCumulative'),                               
    dRmin12 = cms.double(0.3),
    srcMET = cms.InputTag('patMETs'),
    srcPrimaryVertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),
    srcGenParticles = cms.InputTag('genParticles'),                  
    recoMode = cms.string(""),
    doSVreco = cms.bool(True),                          
    SVOptions = cms.PSet(
        usePtBalanceInFit = cms.bool(True),
        useMEtInFit = cms.bool(True),
        useLeg1TrackingInFit = cms.bool(False),
        useLeg2TrackingInFit = cms.bool(False),
        correctPrimaryVertexInFit = cms.bool(False)
    ),
    collinearApproxMassCompatibility = cms.PSet(
        mZ = cms.PSet(
            resonanceMass = cms.double(91.2),
            resonanceWidth = cms.double(2.5),
            metResolutionPx = pfMEtResolutionPx,
            metResolutionPy = pfMEtResolutionPy
        ),
        mAH120 = cms.PSet(
            resonanceMass = cms.double(120),
            resonanceWidth = cms.double(1.),
            metResolutionPx = pfMEtResolutionPx,
            metResolutionPy = pfMEtResolutionPy
        )
    ),
    svFit = cms.PSet(
        psKine = cms.PSet(
            likelihoodFunctions = cms.VPSet(
                svFitLikelihoodDiTauKinematicsPhaseSpace         
            ),
            estUncertainties = cms.PSet(
                numSamplings = cms.int32(-1)
            )
        ),
        psKine_MEt = cms.PSet(
            likelihoodFunctions = cms.VPSet(
                svFitLikelihoodDiTauKinematicsPhaseSpace,
                svFitLikelihoodMEt
            ),
            estUncertainties = cms.PSet(
                numSamplings = cms.int32(-1)
            )
        ),
        psKine_MEt_ptBalance = cms.PSet(
            likelihoodFunctions = cms.VPSet(
                svFitLikelihoodDiTauKinematicsPhaseSpace,
                svFitLikelihoodMEt,
                svFitLikelihoodDiTauPtBalance
            ),
            estUncertainties = cms.PSet(
                #numSamplings = cms.int32(1000)
                numSamplings = cms.int32(-1)
            )
        ##),
        ##psKine_MEt_ptBalance_diTauPt = cms.PSet(
        ##    likelihoodFunctions = cms.VPSet(
        ##        svFitLikelihoodDiTauKinematicsPhaseSpace,
        ##        svFitLikelihoodMEt,
        ##        svFitLikelihoodDiTauPtBalance,
        ##        svFitLikelihoodDiTauPt
        ##    ),
        ##    estUncertainties = cms.PSet(
        ##        numSamplings = cms.int32(-1)
        ##    )
        )
    ),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                           
    verbosity = cms.untracked.int32(0)
)

muTauPairProdConfigurator = objProdConfigurator(
    allMuTauPairs,
    pyModuleName = __name__
)

produceMuTauPairs = muTauPairProdConfigurator.configure(pyNameSpace = locals())

# define additional collections of muon + tau-jet candidates
# with loose track and ECAL isolation applied on muon leg
# (NOTE: to be used for the purpose of factorizing efficiencies
#        of muon isolation from other event selection criteria,
#        in order to avoid problems with limited Monte Carlo statistics)

allMuTauPairsLooseMuonIsolation = allMuTauPairs.clone(
    srcLeg1 = cms.InputTag('selectedPatMuonsTrkIPlooseIsolationCumulative'),
)

muTauPairProdConfiguratorLooseMuonIsolation = objProdConfigurator(
    allMuTauPairsLooseMuonIsolation,
    pyModuleName = __name__
)

produceMuTauPairsLooseMuonIsolation = muTauPairProdConfiguratorLooseMuonIsolation.configure(pyNameSpace = locals())

produceMuTauPairsAll = cms.Sequence(produceMuTauPairs * produceMuTauPairsLooseMuonIsolation)
