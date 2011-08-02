#include "TauAnalysis/Core/plugins/CompositePtrCandidateT1T2MEtHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>
#include <TVector2.h>

const double epsilon = 0.01;

template<typename T1, typename T2>
bool matchesGenCandidatePair(const CompositePtrCandidateT1T2MEt<T1,T2>& compositePtrCandidate)
{
  bool isGenMatched = false;
// not implemented yet...
  return isGenMatched;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

template<typename T1, typename T2>
CompositePtrCandidateT1T2MEtHistManager<T1,T2>::CompositePtrCandidateT1T2MEtHistManager(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg),
    pfMEtCovInverse_(2, 2),
    pfMEtCov_(2, 2),
    pfMEtCovEigenvalues_(2)
{
  //std::cout << "<CompositePtrCandidateT1T2MEtHistManager::CompositePtrCandidateT1T2MEtHistManager>:" << std::endl;

  diTauCandidateSrc_ = cfg.getParameter<edm::InputTag>("diTauCandidateSource");
  //std::cout << " diTauCandidateSrc = " << diTauCandidateSrc_ << std::endl;

  genParticleSrc_ = cfg.getParameter<edm::InputTag>("genParticleSource");
  //std::cout << " genParticleSrc = " << genParticleSrc_ << std::endl;

  pfCandidateSrc_ = cfg.getParameter<edm::InputTag>("pfCandidateSource");
  //std::cout << " pfCandidateSrc = " << pfCandidateSrc_ << std::endl;

  vertexSrc_ = cfg.getParameter<edm::InputTag>("vertexSource");
  //std::cout << " vertexSrc = " << vertexSrc_ << std::endl;

  visMassHypothesisSrc_ = ( cfg.exists("visMassHypothesisSource") ) ?  
    cfg.getParameter<edm::InputTag>("visMassHypothesisSource") : edm::InputTag();
  //std::cout << " visMassHypothesisSrc = " << visMassHypothesisSrc_ << std::endl;

  diTauLeg1WeightExtractors_ = getTauJetWeightExtractors<T1>(cfg, "diTauLeg1WeightSource");
  diTauLeg2WeightExtractors_ = getTauJetWeightExtractors<T2>(cfg, "diTauLeg2WeightSource");

  requireGenMatch_ = cfg.getParameter<bool>("requireGenMatch");
  //std::cout << " requireGenMatch = " << requireGenMatch_ << std::endl;

  pdgIdsElectron_.push_back(+11);
  pdgIdsElectron_.push_back(-11);
  pdgIdsMuon_.push_back(+13);
  pdgIdsMuon_.push_back(-13);
  pdgIdsPhoton_.push_back(22);
  for ( int iQuarkType = 1; iQuarkType <= 6; ++iQuarkType ) {
    pdgIdsJet_.push_back(+iQuarkType);
    pdgIdsJet_.push_back(-iQuarkType);
  }
  pdgIdsJet_.push_back(21);

  std::string normalization_string = cfg.getParameter<std::string>("normalization");
  normMethod_ = getNormMethod(normalization_string, "diTauCandidates");
}

template<typename T1, typename T2>
CompositePtrCandidateT1T2MEtHistManager<T1,T2>::~CompositePtrCandidateT1T2MEtHistManager()
{
  for ( typename std::vector<FakeRateJetWeightExtractor<T1>*>::iterator it = diTauLeg1WeightExtractors_.begin();
	it != diTauLeg1WeightExtractors_.end(); ++it ) {
    delete (*it);
  }

  for ( typename std::vector<FakeRateJetWeightExtractor<T2>*>::iterator it = diTauLeg2WeightExtractors_.begin();
	it != diTauLeg2WeightExtractors_.end(); ++it ) {
    delete (*it);
  }
}

template<typename T1, typename T2>
void CompositePtrCandidateT1T2MEtHistManager<T1,T2>::bookHistogramsImp()
{
  //std::cout << "<CompositePtrCandidateT1T2MEtHistManager::bookHistogramsImp>:" << std::endl;

  hGenDiTauCandidatePt_ = book1D("GenDiTauCandidatePt", "gen. Composite P_{T}", 75, 0., 150.);
  hGenDiTauCandidateVisPt_ = book1D("GenDiTauCandidateVisPt", "gen. Composite visible P_{T}", 75, 0., 150.);
  hGenDiTauCandidateEta_ = book1D("GenDiTauCandidateEta", "gen. Composite #eta", 100, -5., +5.);
  hGenDiTauCandidatePhi_ = book1D("GenDiTauCandidatePhi", "gen. Composite #phi", 36, -TMath::Pi(), +TMath::Pi());
  
  hGenLeg1En_ = book1D("GenLeg1En", "gen. leg_{1} Energy", 125, 0., 250.);
  hGenLeg2En_ = book1D("GenLeg2En", "gen. leg_{2} Energy", 125, 0., 250.);
  hGenLeg1PtVsLeg2Pt_ = book2D("GenLeg1PtVsLeg2Pt", "gen. leg_{1} P_{T} vs. leg_{2} P_{T}", 20, 0., 100., 20, 0., 100.);
  hGenLeg1EtaVsLeg2Eta_ = book2D("GenLeg1EtaVsLeg2Eta", "gen. leg_{1} #eta vs. leg_{2} #eta", 20, -2.5, 2.5, 20, -2.5, 2.5);
  hGenDeltaRleg1VisNu_ = book1D("GenDeltaRleg1VisNu", "gen. leg_{1} dR(pVis, pNu)", 102, -0.01, 1.01);
  hGenDeltaRleg2VisNu_ = book1D("GenDeltaRleg2VisNu", "gen. leg_{2} dR(pVis, pNu)", 102, -0.01, 1.01);
  hGenLeg1DecayTime_ = book1D("GenLeg1DecayTime", "gen. leg_{1} Decay eigentime", 100, 0., 1000.);
  hGenLeg2DecayTime_ = book1D("GenLeg2DecayTime", "gen. leg_{2} Decay eigentime", 100, 0., 1000.);
  hGenMass_ = book1D("GenMass", "gen. leg_{1} + leg_{2} Invariant Mass", 50, 0., 250.);

  hGenSqrtX1X2VsDPhi12_ = bookProfile1D("GenSqrtX1X2VsDPhi12", 
					"gen. sqrt(X_{1} * X_{2}) vs. #Delta#phi_{1,2}", 18, -epsilon, TMath::Pi() + epsilon);
  hGenSqrtX1X2VsVisEnRatio_ = bookProfile1D("GenSqrtX1X2VsVisEnRatio", 
					    "gen. sqrt(X_{1} * X_{2}) vs. vis. E_{1} / E_{2}", 25, -epsilon, 2.51);
  hGenSqrtX1X2VsVisPtRatio_ = bookProfile1D("GenSqrtX1X2VsVisPtRatio", 
					    "gen. sqrt(X_{1} * X_{2}) vs. vis. P_{T}^{1} / P_{T}^{2}", 25, -epsilon, 2.51);

  hGenLeg1TauPlusDecayAngleLepton_ = book1D("GenLeg1TauPlusDecayAngleLepton", 
					    "gen. leg_{1} #theta(#tau, #ell) for lep. Tau+ decays", 36, 0., TMath::Pi());
  hGenLeg1TauPlusDecayAngleOneProng_ = book1D("GenLeg1TauPlusDecayAngleOneProng", 
					      "gen. leg_{1} #theta(#tau, #nu) for had. one-prong Tau+ decays", 36, 0., TMath::Pi());
  hGenLeg1TauPlusDecayAngleThreeProng_ = book1D("GenLeg1TauPlusDecayAngleThreeProng", 
					      "gen. leg_{1} #theta(#tau, #nu) for had. three-prong Tau+ decays", 36, 0., TMath::Pi());
  hGenLeg1TauMinusDecayAngleLepton_ = book1D("GenLeg1TauMinusDecayAngleLepton", 
					     "gen. leg_{1} #theta(#tau, #ell) for lep. Tau- decays", 36, 0., TMath::Pi());
  hGenLeg1TauMinusDecayAngleOneProng_ = book1D("GenLeg1TauMinusDecayAngleOneProng", 
					       "gen. leg_{1} #theta(#tau, #nu) for had. one-prong Tau- decays", 36, 0., TMath::Pi());
  hGenLeg1TauMinusDecayAngleThreeProng_ = book1D("GenLeg1TauMinusDecayAngleThreeProng", 
						 "gen. leg_{1} #theta(#tau, #nu) for had. three-prong Tau- decays", 36, 0., TMath::Pi());
  hGenLeg2TauPlusDecayAngleLepton_ = book1D("GenLeg2TauPlusDecayAngleLepton", 
					    "gen. leg_{2} #theta(#tau, #ell) for lep. Tau+ decays", 36, 0., TMath::Pi());
  hGenLeg2TauPlusDecayAngleOneProng_ = book1D("GenLeg2TauPlusDecayAngleOneProng", 
					      "gen. leg_{2} #theta(#tau, #nu) for had. one-prong Tau+ decays", 36, 0., TMath::Pi());
  hGenLeg2TauPlusDecayAngleThreeProng_ = book1D("GenLeg2TauPlusDecayAngleThreeProng", 
					      "gen. leg_{2} #theta(#tau, #nu) for had. three-prong Tau+ decays", 36, 0., TMath::Pi());
  hGenLeg2TauMinusDecayAngleLepton_ = book1D("GenLeg2TauMinusDecayAngleLepton", 
					     "gen. leg_{2} #theta(#tau, #ell) for lep. Tau- decays", 36, 0., TMath::Pi());
  hGenLeg2TauMinusDecayAngleOneProng_ = book1D("GenLeg2TauMinusDecayAngleOneProng", 
					       "gen. leg_{2} #theta(#tau, #nu) for had. one-prong Tau- decays", 36, 0., TMath::Pi());
  hGenLeg2TauMinusDecayAngleThreeProng_ = book1D("GenLeg2TauMinusDecayAngleThreeProng", 
						 "gen. leg_{2} #theta(#tau, #nu) for had. three-prong Tau- decays", 36, 0., TMath::Pi());
  
  hDiTauCandidatePt_ = book1D("DiTauCandidatePt", "Composite P_{T}", 75, 0., 150.);
  hDiTauCandidateVisPt_ = book1D("DiTauCandidateVisPt", "Composite visible P_{T}", 75, 0., 150.);
  hDiTauCandidateEta_ = book1D("DiTauCandidateEta", "Composite #eta", 100, -5., +5.);
  hDiTauCandidatePhi_ = book1D("DiTauCandidatePhi", "Composite #phi", 36, -TMath::Pi(), +TMath::Pi());
  hDiTauCandidateCharge_ = book1D("DiTauCandidateCharge", "Composite Charge", 11, -5.5, +5.5);
  hDiTauCandidateMass_ = book1D("DiTauCandidateMass", "Composite Invariant Mass", 50, 0., 250.);
  
  hLeg1PtVsLeg2Pt_ = book2D("Leg1PtVsLeg2Pt", "leg_{1} P_{T} vs. leg_{2} P_{T}", 20, 0., 100., 20, 0., 100.);
  hLeg1EtaVsLeg2Eta_ = book2D("Leg1EtaVsLeg2Eta", "leg_{1} #eta vs. leg_{2} #eta", 20, -2.5, 2.5, 20, -2.5, 2.5);
  hLeg1IsoVsLeg2Iso_ = book2D("Leg1IsoVsLeg2Iso", "leg_{1} Iso vs. leg_{2} Iso", 20, 0., 2., 20, 0., 2.);
  
  bookWeightHistograms(*dqmStore_, "DiTauCandidateWeight", "Composite Weight", 
		       hDiTauCandidateWeightPosLog_, hDiTauCandidateWeightNegLog_, hDiTauCandidateWeightZero_, 
		       hDiTauCandidateWeightLinear_);

  hDiTauCandidateImpParSig_ = book1D("DiTauCandidateImpParSig", "#sqrt{#frac{dxy_{1}}{#Delta dxy_{1}}^{2}+#frac{dxy_{2}}{#Delta dxy_{2}}^{2}}", 30, 0., 15.);
  
  hVisPt_ = book1D("VisPt", "Visible P_{T}", 50, 0., 100.);
  hVisPhi_ = book1D("VisPhi", "Visible #phi", 36, -TMath::Pi(), +TMath::Pi());
  hVisMass_ = book1D("VisMass", "Visible Mass", 50, 0., 250.);
  hVisMassL_ = book1D("VisMassL", "Visible Mass", 100, 0., 500.);
  hVisMassXL_ = book1D("VisMassXL", "Visible Mass", 150, 0., 750.);
  hVisMassGenLeg2Electron_ = book1D("VisMassGenLeg2Electron", "Visible Mass (rec. Tau matching gen. Electron)", 50, 0., 250.);
  hVisMassGenLeg2Muon_ = book1D("VisMassGenLeg2Muon", "Visible Mass (rec. Tau matching gen. Muon)", 50, 0., 250.);
  hVisMassGenLeg2Photon_ = book1D("VisMassGenLeg2Photon", "Visible Mass (rec. Tau matching gen. Photon)", 50, 0., 250.);
  hVisMassGenLeg2Jet_ = book1D("VisMassGenLeg2Jet", "Visible Mass (rec. Tau matching gen. quark/gluon Jet)", 50, 0., 250.);
  hVisMassZllCombinedHypothesis_ = ( visMassHypothesisSrc_.label() != "" ) ?
    book1D("VisMassZllCombinedHypothesis", "Visible Mass (combined Value of different Event Hypotheses)", 40, 0., 200.) : 0;

  hVisMassDPhiGt175_ = book1D("VisMassDPhiGt175", "VisMassDPhiGt175", 100, 0., 500.);
  hVisMassDPhi170to175_ = book1D("VisMassDPhi170to175", "VisMassDPhi170to175", 100, 0., 500.);
  hVisMassDPhi160to170_ = book1D("VisMassDPhi160to170", "VisMassDPhi160to170", 100, 0., 500.);
  hVisMassDPhi140to160_ = book1D("VisMassDPhi140to160", "VisMassDPhi140to160", 100, 0., 500.);
  hVisMassDPhi120to140_ = book1D("VisMassDPhi120to140", "VisMassDPhi120to140", 100, 0., 500.);
  hVisMassDPhiLt120_ = book1D("VisMassDPhiLt120", "VisMassDPhiLt120", 100, 0., 500.);
  
  hCollinearApproxEta_ = book1D("CollinearApproxEta", "Collinear Approximation #eta", 100, -5., +5.);
  hCollinearApproxMass_ = book1D("CollinearApproxMass", "Collinear Approximation Mass", 50, 0., 250.);
  hCollinearApproxX1_ = book1D("CollinearApproxX1", "Collinear Approximation X_{1}", 100, -2.5, +2.5);
  hCollinearApproxX2_ = book1D("CollinearApproxX2", "Collinear Approximation X_{2}", 100, -2.5, +2.5);

  hCDFmethodMass_ = book1D("CDFmethodMass", "CDF Method Mass", 50, 0., 250.);

  hTauPairMassMin_ = book1D("TauPairMassMin", "TauPairMassMin", 50, 0., 250.);
  hTauPairMassMinL_ = book1D("TauPairMassMinL", "TauPairMassMinL", 100, 0., 500.);
  hTauPairMassMinXL_ = book1D("TauPairMassMinXL", "TauPairMassMinXL", 150, 0., 750.);

  hTauPairMassMinDPhiGt175_ = book1D("TauPairMassMinDPhiGt175", "TauPairMassMinDPhiGt175", 100, 0., 500.);
  hTauPairMassMinDPhi170to175_ = book1D("TauPairMassMinDPhi170to175", "TauPairMassMinDPhi170to175", 100, 0., 500.);
  hTauPairMassMinDPhi160to170_ = book1D("TauPairMassMinDPhi160to170", "TauPairMassMinDPhi160to170", 100, 0., 500.);
  hTauPairMassMinDPhi140to160_ = book1D("TauPairMassMinDPhi140to160", "TauPairMassMinDPhi140to160", 100, 0., 500.);
  hTauPairMassMinDPhi120to140_ = book1D("TauPairMassMinDPhi120to140", "TauPairMassMinDPhi120to140", 100, 0., 500.);
  hTauPairMassMinDPhiLt120_ = book1D("TauPairMassMinDPhiLt120", "TauPairMassMinDPhiLt120", 100, 0., 500.);

  hEvMETresParGenDiTau_ = book1D("EvMETresParGenDiTau", "MET resolution parallel gen. diTau", 100, -50., 50.);
  hEvMETresPerpGenDiTau_ = book1D("EvMETresPerpGenDiTau", "MET resolution perp. gen. diTau", 50, -0.01, 50.);
  hDtMETresParGenDiTau_ = book1D("DtMETresParGenDiTau", "MET resolution parallel gen. diTau", 100, -50., 50.);
  hDtMETresPerpGenDiTau_ = book1D("DtMETresPerpGenDiTau", "MET resolution perp. gen. diTau", 50, -0.01, 50.);  
  
  hMt12MET_ = book1D("Mt12MET", "Mass_{T 1,2,MET}", 50, 0., 250.);
  
  hMt1MET_ = book1D("Mt1MET", "Mass_{T 1,MET}", 40, 0., 200.);
  hPt1MET_ = book1D("Pt1MET", "P_{T}^{T 1,MET}", 40, 0., 200.);
  hMt1METvisMassRatio_ = book1D("Mt1METvisMassRatio", "Mass_{T 1,MET} / Visible Mass", 100, -5.01, +5.01);
  hMt2MET_ = book1D("Mt2MET", "Mass_{T 2,MET}", 40, 0., 200.);
  hPt2MET_ = book1D("Pt2MET", "P_{T}^{T 2,MET}", 40, 0., 200.);
  hMt2METvisMassRatio_ = book1D("Mt2METvisMassRatio", "Mass_{T 2,MET} / Visible Mass", 100, -5.01, +5.01);
  
  hHt12MET_ = book1D("Ht12MET", "P_{T}^{1} + P_{T}^{2} + MET", 50, 0., 250.);

  hDPhi12_ = book1D("DPhi12", "#Delta#phi_{1,2}", 36, -epsilon, TMath::Pi() + epsilon);
  hDR12_ = book1D("DR12", "#Delta R_{1,2}", 51, -0.1, 10.1);
  
  hVisEtaMin_ = book1D("VisEtaMin", "VisEtaMin", 60, -3., +3.);
  hVisEtaMax_ = book1D("VisEtaMax", "VisEtaMax", 60, -3., +3.);
  
  hDPhi1MET_ = book1D("DPhi1MET", "#Delta#phi_{1,MET}", 36, -epsilon, TMath::Pi() + epsilon);
  hDPhi2MET_ = book1D("DPhi2MET", "#Delta#phi_{2,MET}", 36, -epsilon, TMath::Pi() + epsilon);
  hDPhi1METvsDPhi2MET_ = book2D("DPhi1METvsDPhi2MET", "#Delta#phi_{1,MET} vs. #Delta#phi_{2,MET}", 18, -epsilon, TMath::Pi() + epsilon, 18, -epsilon, TMath::Pi() + epsilon); 

  hPzetaCorr_ = book2D("PzetaCorr", "P_{#zeta} vs. P_{#zeta}^{vis}", 20, 0., 200., 40., -200., +200.);
  hPzetaDiff_ = book1D("PzetaDiff", "P_{#zeta} - 1.5*P_{#zeta}^{vis}", 40, -100., +100.);

  hPzetaDiffMEtSignRatio_ = book1D("PzetaDiffMEtSignRatio", "(P_{#zeta} - 1.5*P_{#zeta}^{vis}) / MET sign.", 100, -5.01, +5.01);
  
  hMEtSignProb_ = book1D("MEtSignProb", "Probability of (gen. MET - rec. MET)/sigmaMET", 102., -0.01, 1.01);
  hMEtSignEigenValue1_ = book1D("MEtSignEigenValue1", "MET Sign. 1st Eigenvalue", 100., -0.01, 50.);
  hMEtSignEigenValue2_ = book1D("MEtSignEigenValue2", "MEt Sign. 2nd Eigenvalue", 100., -0.01, 50.);

  hPzetaDiffVsDPhi12_ = 
    book2D("PzetaDiffVsDPhi12",  "P_{#zeta} - 1.5*P_{#zeta}^{vis} vs. #Delta#phi_{1,2}", 36, -epsilon, TMath::Pi() + epsilon, 20, -100., +100.);
  hPzetaDiffVsMt1MET_  = 
    book2D("PzetaDiffVsMt1MET",  "P_{#zeta} - 1.5*P_{#zeta}^{vis} vs. Mass_{T 1,MET}", 20, 0., 200., 20, -100., +100.);
  hPzetaDiffVsMET_  = 
    book2D("PzetaDiffVsMET",     "P_{#zeta} - 1.5*P_{#zeta}^{vis} vs. MET",            20, 0., 200., 20, -100., +100.);
  hPzetaDiffVsVisMass_ = 
    book2D("PzetaDiffVsVisMass", "P_{#zeta} - 1.5*P_{#zeta}^{vis} vs. Visible Mass",   20, 0., 500., 20, -100., +100.);
  hMt1METvsMET_        = 
    book2D("Mt1METvsMET",        "Mass_{T 1,MET} vs. MET",                             20, 0., 200., 20,    0.,  200.);
  hMt1METvsVisMass_    = 
    book2D("Mt1METvsVisMass",    "Mass_{T 1,MET} vs. Visible Mass",                    20, 0., 500., 20, -100., +100.);
}

template<typename T1, typename T2>
double CompositePtrCandidateT1T2MEtHistManager<T1,T2>::getDiTauCandidateWeight(const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate)
{
  double diTauLeg1Weight = getTauJetWeight<T1>(*diTauCandidate.leg1(), diTauLeg1WeightExtractors_);
  double diTauLeg2Weight = getTauJetWeight<T2>(*diTauCandidate.leg2(), diTauLeg2WeightExtractors_);
  return (diTauLeg1Weight*diTauLeg2Weight);
}

void fillDeltaRvisNuHistogram(MonitorElement* h, 
			      const reco::Candidate::LorentzVector& p4, const reco::Candidate::LorentzVector& p4Vis, double weight)
{
  reco::Candidate::LorentzVector p4Nu = p4 - p4Vis;
  double dR = deltaR(p4Vis, p4Nu);
  h->Fill(TMath::Min(dR, 1.), weight);
}

void fillGenTauHistograms(MonitorElement* hGenTauPlusDecayAngleLepton, 
			  MonitorElement* hGenTauPlusDecayAngleOneProng, 
			  MonitorElement* hGenTauPlusDecayAngleThreeProng,
			  MonitorElement* hGenTauMinusDecayAngleLepton, 
			  MonitorElement* hGenTauMinusDecayAngleOneProng, 
			  MonitorElement* hGenTauMinusDecayAngleThreeProng,
			  const reco::Candidate::LorentzVector& genTauMomentum, 
			  const reco::GenParticleCollection& genParticles,
			  double weight)
{
  //std::cout << "<fillGenTauHistograms>:" << std::endl;

  const reco::GenParticle* genTau = findGenParticle(genTauMomentum, genParticles);
  //std::cout << " genTau = " << genTau << std::endl;

  if ( genTau ) {
    std::vector<const reco::GenParticle*> genTauDecayProducts;
    findDaughters(genTau, genTauDecayProducts);

    reco::Candidate::LorentzVector genTauVisMomentum = getVisMomentum(genTauDecayProducts);
    reco::Candidate::LorentzVector genTauInvisMomentum = getInvisMomentum(genTauDecayProducts);
    
    reco::Candidate::LorentzVector genTauVisMomentum_restframe = boostToRestFrame(genTauVisMomentum, genTauMomentum);
    reco::Candidate::LorentzVector genTauInvisMomentum_restframe = boostToRestFrame(genTauInvisMomentum, genTauMomentum);

    std::string genTauDecayMode = getGenTauDecayMode(genTau);
    //std::cout << " genTauDecayMode = " << genTauDecayMode << std::endl;

    if ( genTauDecayMode == "electron" || 
	 genTauDecayMode == "muon" ) {	  
      double genTauVisDecayAngle = angle(genTauVisMomentum_restframe, genTau->p4());
      if ( genTau->pdgId() == -15 ) hGenTauPlusDecayAngleLepton->Fill(genTauVisDecayAngle, weight);
      if ( genTau->pdgId() == +15 ) hGenTauMinusDecayAngleLepton->Fill(genTauVisDecayAngle, weight);
    }
	
    if ( genTauDecayMode == "oneProng0Pi0" ||
	 genTauDecayMode == "oneProng1Pi0" ||
	 genTauDecayMode == "oneProng2Pi0" ) {
      double genTauInvisDecayAngle = angle(genTauInvisMomentum_restframe, genTau->p4());
      if ( genTau->pdgId() == -15 ) hGenTauPlusDecayAngleOneProng->Fill(genTauInvisDecayAngle, weight);
      if ( genTau->pdgId() == +15 ) hGenTauMinusDecayAngleOneProng->Fill(genTauInvisDecayAngle, weight);
    }
	
    if ( genTauDecayMode == "threeProng0Pi0" ||
	 genTauDecayMode == "threeProng1Pi0" ) {
      double genTauInvisDecayAngle = angle(genTauInvisMomentum_restframe, genTau->p4());
      if ( genTau->pdgId() == -15 ) hGenTauPlusDecayAngleThreeProng->Fill(genTauInvisDecayAngle, weight);
      if ( genTau->pdgId() == +15 ) hGenTauMinusDecayAngleThreeProng->Fill(genTauInvisDecayAngle, weight);
    }
  }
}
     
void fillMEtResolutionHistograms(const reco::Candidate::LorentzVector& p4RecMEt, const reco::Candidate::LorentzVector& p4GenMEt,
				 double projCosPhi, double projSinPhi, 
				 MonitorElement* hResPar, MonitorElement* hResPerp, double weight)
{
  reco::Candidate::LorentzVector residual = p4RecMEt - p4GenMEt;
  	  
  double resPar  = residual.px()*projCosPhi + residual.py()*projSinPhi;
  double resPerp = TMath::Abs(residual.px()*projSinPhi - residual.py()*projCosPhi);
	  
  hResPar->Fill(resPar, weight);
  hResPerp->Fill(resPerp, weight);
}

template<typename T1, typename T2>
void CompositePtrCandidateT1T2MEtHistManager<T1,T2>::fillHistogramsImp(const edm::Event& evt, const edm::EventSetup& es, double evtWeight)
{  
  //std::cout << "<CompositePtrCandidateT1T2MEtHistManager::fillHistogramsImp>:" << std::endl; 

  typedef std::vector<CompositePtrCandidateT1T2MEt<T1,T2> > CompositePtrCandidateCollection;
  edm::Handle<CompositePtrCandidateCollection> diTauCandidates;
  getCollection(evt, diTauCandidateSrc_, diTauCandidates);

  //std::cout << " diTauCandidateSrc = " << diTauCandidateSrc_.label() << ":" 
  //	      << " diTauCandidates.size = " << diTauCandidates->size() << std::endl;

  edm::Handle<reco::GenParticleCollection> genParticles;
  if ( genParticleSrc_.label() != "" ) evt.getByLabel(genParticleSrc_, genParticles);

//--- CV: collection of PFCandidate only needed in case histograms of MET resolution 
//        projected parallel and perpendicular to diTau direction are to be filled
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  if ( makeMEtProjResolutionHistograms_ ) {
    evt.getByLabel(pfCandidateSrc_, pfCandidates);
  }

  edm::Handle<std::vector<reco::Vertex> > recoVertices;
  evt.getByLabel(vertexSrc_, recoVertices);

  double diTauCandidateWeightSum = 0.;
  for ( typename CompositePtrCandidateCollection::const_iterator diTauCandidate = diTauCandidates->begin(); 
	diTauCandidate != diTauCandidates->end(); ++diTauCandidate ) {
    if ( requireGenMatch_ && !matchesGenCandidatePair(*diTauCandidate) ) continue;

    diTauCandidateWeightSum += getDiTauCandidateWeight(*diTauCandidate);
  }

  for ( typename CompositePtrCandidateCollection::const_iterator diTauCandidate = diTauCandidates->begin(); 
	diTauCandidate != diTauCandidates->end(); ++diTauCandidate ) {

    //bool isGenMatched = matchesGenCandidatePair(*diTauCandidate);
    //std::cout << " Pt = " << diTauCandidate->pt() << ", phi = " << diTauCandidate->phi() << ", visMass = " << diTauCandidate->p4Vis().mass() << std::endl;
    //std::cout << " isGenMatched = " << isGenMatched << std::endl;

    if ( requireGenMatch_ && !matchesGenCandidatePair(*diTauCandidate) ) continue;

    double diTauCandidateWeight = getDiTauCandidateWeight(*diTauCandidate);
    double weight = getWeight(evtWeight, diTauCandidateWeight, diTauCandidateWeightSum);

    if ( diTauCandidate->p4Leg1gen().energy() > epsilon && 
	 diTauCandidate->p4Leg2gen().energy() > epsilon ) {
      hGenDiTauCandidatePt_->Fill(diTauCandidate->p4gen().pt(), weight);
      hGenDiTauCandidateVisPt_->Fill(diTauCandidate->p4VisGen().pt(), weight);
      hGenDiTauCandidateEta_->Fill(diTauCandidate->p4gen().eta(), weight);
      hGenDiTauCandidatePhi_->Fill(diTauCandidate->p4gen().phi(), weight);
      
      hGenLeg1En_->Fill(diTauCandidate->p4Leg1gen().energy(), weight);
      hGenLeg2En_->Fill(diTauCandidate->p4Leg2gen().energy(), weight);

      hGenLeg1PtVsLeg2Pt_->Fill(diTauCandidate->p4VisLeg1gen().pt(), diTauCandidate->p4VisLeg2gen().pt(), weight);
      hGenLeg1EtaVsLeg2Eta_->Fill(diTauCandidate->p4VisLeg1gen().eta(), diTauCandidate->p4VisLeg2gen().eta(), weight);

      fillDeltaRvisNuHistogram(hGenDeltaRleg1VisNu_, diTauCandidate->p4Leg1gen(), diTauCandidate->p4VisLeg1gen(), weight);
      fillDeltaRvisNuHistogram(hGenDeltaRleg2VisNu_, diTauCandidate->p4Leg2gen(), diTauCandidate->p4VisLeg2gen(), weight);

      double sqrtX1X2gen = TMath::Sqrt(diTauCandidate->x1gen()*diTauCandidate->x2gen());
      double dPhi12gen = TMath::Abs(normalizedPhi(diTauCandidate->p4VisLeg1gen().phi() - diTauCandidate->p4VisLeg2gen().phi()));
/*
  
  CV: temporary work-around until MonitorElement::Fill(double, double, double) is fixed for TProfiles

      hGenSqrtX1X2VsDPhi12_->Fill(dPhi12gen, TMath::Sqrt(diTauCandidate->x1gen()*diTauCandidate->x2gen()), weight);
 */
      hGenSqrtX1X2VsDPhi12_->getTProfile()->Fill(dPhi12gen, sqrtX1X2gen, weight);
      double visEnRatio = ( diTauCandidate->p4VisLeg2gen().energy() > 0. ) ?
	diTauCandidate->p4VisLeg1gen().energy()/diTauCandidate->p4VisLeg2gen().energy() : -1.;
      hGenSqrtX1X2VsVisEnRatio_->getTProfile()->Fill(TMath::Min(visEnRatio, 2.5), sqrtX1X2gen, weight);
      double visPtRatio = ( diTauCandidate->p4VisLeg2gen().pt() > 0. ) ?
	diTauCandidate->p4VisLeg1gen().pt()/diTauCandidate->p4VisLeg2gen().pt() : -1.;
      hGenSqrtX1X2VsVisPtRatio_->getTProfile()->Fill(TMath::Min(visPtRatio, 2.5), sqrtX1X2gen, weight);

      hGenLeg1DecayTime_->Fill(compDecayEigenTime(diTauCandidate->decayVertexPosLeg1gen(), 
						  diTauCandidate->primaryVertexPosGen(), diTauCandidate->p4Leg1gen().energy()), weight);
      hGenLeg2DecayTime_->Fill(compDecayEigenTime(diTauCandidate->decayVertexPosLeg2gen(), 
						  diTauCandidate->primaryVertexPosGen(), diTauCandidate->p4Leg2gen().energy()), weight);

      hGenMass_->Fill(diTauCandidate->p4gen().mass(), weight);

      if ( genParticles.isValid() ) { 
	fillGenTauHistograms(
          hGenLeg1TauPlusDecayAngleLepton_, hGenLeg1TauPlusDecayAngleOneProng_, hGenLeg1TauPlusDecayAngleThreeProng_,
	  hGenLeg1TauMinusDecayAngleLepton_, hGenLeg1TauMinusDecayAngleOneProng_, hGenLeg1TauMinusDecayAngleThreeProng_,
	  diTauCandidate->p4Leg1gen(), *genParticles, weight);
	fillGenTauHistograms(
	  hGenLeg2TauPlusDecayAngleLepton_, hGenLeg2TauPlusDecayAngleOneProng_, hGenLeg2TauPlusDecayAngleThreeProng_,
	  hGenLeg2TauMinusDecayAngleLepton_, hGenLeg2TauMinusDecayAngleOneProng_, hGenLeg2TauMinusDecayAngleThreeProng_,
	  diTauCandidate->p4Leg2gen(), *genParticles, weight);
      }
    }

    hDiTauCandidatePt_->Fill(diTauCandidate->pt(), weight);
    hDiTauCandidateVisPt_->Fill(diTauCandidate->p4Vis().pt(), weight);
    hDiTauCandidateEta_->Fill(diTauCandidate->eta(), weight);
    hDiTauCandidatePhi_->Fill(diTauCandidate->phi(), weight);
    hDiTauCandidateCharge_->Fill(diTauCandidate->charge(), weight);
    hDiTauCandidateMass_->Fill(diTauCandidate->mass(), weight);

    if ( diTauCandidate->p4gen().pt() > 10. && diTauCandidate->met().isNonnull() ) {
      double projCosPhi = TMath::Cos(diTauCandidate->p4gen().phi());
      double projSinPhi = TMath::Sin(diTauCandidate->p4gen().phi());

      if ( dynamic_cast<const pat::MET*>(diTauCandidate->met().get()) != 0 ) {
	const pat::MET* patMEt = dynamic_cast<const pat::MET*>(diTauCandidate->met().get());
	if ( patMEt->genMET() ) 
	  fillMEtResolutionHistograms(diTauCandidate->met()->p4(), patMEt->genMET()->p4(), projCosPhi, projSinPhi,
				      hEvMETresParGenDiTau_, hEvMETresPerpGenDiTau_, weight);	
      }

      fillMEtResolutionHistograms(diTauCandidate->met()->p4(), diTauCandidate->p4InvisGen(), projCosPhi, projSinPhi,
				  hDtMETresParGenDiTau_, hDtMETresPerpGenDiTau_, weight);
    }

    hLeg1PtVsLeg2Pt_->Fill(diTauCandidate->leg1()->pt(), diTauCandidate->leg2()->pt(), weight);
    hLeg1EtaVsLeg2Eta_->Fill(diTauCandidate->leg1()->eta(), diTauCandidate->leg2()->eta(), weight);
    double isoLeg1 = (diTauCandidate->leg1().get()->userIsolation(pat::TrackIso)
                     + diTauCandidate->leg1().get()->userIsolation(pat::EcalIso)
                     + diTauCandidate->leg1().get()->userIsolation(pat::HcalIso))/diTauCandidate->leg1()->pt();
    double isoLeg2 = (diTauCandidate->leg2().get()->userIsolation(pat::TrackIso)
                     + diTauCandidate->leg2().get()->userIsolation(pat::EcalIso)
		     + diTauCandidate->leg2().get()->userIsolation(pat::HcalIso))/diTauCandidate->leg2()->pt();
    hLeg1IsoVsLeg2Iso_->Fill(isoLeg1, isoLeg2, weight);

    fillWeightHistograms(hDiTauCandidateWeightPosLog_, hDiTauCandidateWeightNegLog_, hDiTauCandidateWeightZero_, 
			 hDiTauCandidateWeightLinear_, diTauCandidateWeight);

    const reco::Track* trackLeg1 = trackExtractorLeg1_(*diTauCandidate->leg1());
    const reco::Track* trackLeg2 = trackExtractorLeg2_(*diTauCandidate->leg2());
    if ( trackLeg1 && trackLeg2 &&
	 recoVertices->size() >= 1 ) {
      const reco::Vertex& thePrimaryEventVertex = (*recoVertices->begin());

      double trackLeg1IpSig = trackLeg1->dxy(thePrimaryEventVertex.position())/trackLeg1->dxyError();
      double trackLeg2IpSig = trackLeg2->dxy(thePrimaryEventVertex.position())/trackLeg2->dxyError();

      hDiTauCandidateImpParSig_->Fill(TMath::Sqrt(trackLeg1IpSig*trackLeg1IpSig + trackLeg2IpSig*trackLeg2IpSig), weight);
    }

    hVisPt_->Fill(diTauCandidate->p4Vis().pt(), weight);
    hVisPhi_->Fill(diTauCandidate->p4Vis().phi(), weight);
    double visMass = diTauCandidate->p4Vis().mass();   
    hVisMass_->Fill(visMass, weight);
    hVisMassL_->Fill(visMass, weight);
    hVisMassXL_->Fill(visMass, weight);
    if ( genParticles.isValid() ) { 
      fillHistogramGenMatch(hVisMassGenLeg2Electron_, visMass, diTauCandidate->leg2()->p4(), *genParticles, pdgIdsElectron_, weight);
      fillHistogramGenMatch(hVisMassGenLeg2Muon_, visMass, diTauCandidate->leg2()->p4(), *genParticles, pdgIdsMuon_, weight);
      fillHistogramGenMatch(hVisMassGenLeg2Photon_, visMass, diTauCandidate->leg2()->p4(), *genParticles, pdgIdsPhoton_, weight);
      fillHistogramGenMatch(hVisMassGenLeg2Jet_, visMass, diTauCandidate->leg2()->p4(), *genParticles, pdgIdsJet_, weight);
    }
    if ( visMassHypothesisSrc_.label() != "" ) {
      typedef edm::RefProd<CompositePtrCandidateCollection> CompositePtrCandidateRefProd;
      typedef std::vector<float> vfloat;
      typedef edm::AssociationVector<CompositePtrCandidateRefProd, vfloat> ZtautauVisMassHypothesisCollection;
      edm::Handle<ZtautauVisMassHypothesisCollection> visMassHypotheses;
      evt.getByLabel(visMassHypothesisSrc_, visMassHypotheses);

      unsigned numDiTauCandidates = diTauCandidates->size();
      for ( unsigned idxDiTauCandidate = 0; idxDiTauCandidate < numDiTauCandidates; ++idxDiTauCandidate ) {
	edm::Ref<CompositePtrCandidateCollection> diTauCandidateRef(diTauCandidates, idxDiTauCandidate);

	if ( diTauCandidateRef.get() != &(*diTauCandidate) ) continue;

	if ( visMassHypotheses->keyProduct().id() == diTauCandidateRef.id() ) {
	  hVisMassZllCombinedHypothesis_->Fill((*visMassHypotheses)[diTauCandidateRef], weight);
	} else {
          // CV: print warning, but for first event only !!
	  static bool isFirstWarning = true;
	  if ( isFirstWarning )
  	    edm::LogWarning("CompositePtrCandidateT1T2MEtHistManager::fillHistogramsImp") 
	      << " Mismatch between Product ids. of ZtautauVisMassHypothesisCollection  = " << visMassHypotheses->keyProduct().id()
	      << " and CompositePtrCandidateCollection = " << diTauCandidateRef.id() << ", "
	      << " presumably CompositePtrCandidateCollection has been filtered without re-keying ZtautauVisMassHypothesisCollection !!";
          isFirstWarning = false;
	}
      }
    }

    double dPhi = diTauCandidate->dPhi12()*180./TMath::Pi();
    if      ( dPhi > 175. ) hVisMassDPhiGt175_->Fill(visMass, weight);
    else if ( dPhi > 170. ) hVisMassDPhi170to175_->Fill(visMass, weight);
    else if ( dPhi > 160. ) hVisMassDPhi160to170_->Fill(visMass, weight);
    else if ( dPhi > 140. ) hVisMassDPhi140to160_->Fill(visMass, weight);
    else if ( dPhi > 120. ) hVisMassDPhi120to140_->Fill(visMass, weight);
    else                    hVisMassDPhiLt120_->Fill(visMass, weight);
    
    if ( diTauCandidate->collinearApproxIsValid() ) {  
      hCollinearApproxEta_->Fill(diTauCandidate->p4CollinearApprox().eta(), weight);
      hCollinearApproxMass_->Fill(diTauCandidate->p4CollinearApprox().mass(), weight);
      hCollinearApproxX1_->Fill(diTauCandidate->x1CollinearApprox(), weight);
      hCollinearApproxX2_->Fill(diTauCandidate->x2CollinearApprox(), weight);
    }

    hCDFmethodMass_->Fill(diTauCandidate->p4CDFmethod().mass(), weight);

    if ( diTauCandidate->mTauTauMin_isValid() ) {
      double mTauTauMin = diTauCandidate->mTauTauMin();      
      hTauPairMassMin_->Fill(mTauTauMin, weight);
      hTauPairMassMinL_->Fill(mTauTauMin, weight);
      hTauPairMassMinXL_->Fill(mTauTauMin, weight);
      
      if      ( dPhi > 175. ) hTauPairMassMinDPhiGt175_->Fill(mTauTauMin, weight);
      else if ( dPhi > 170. ) hTauPairMassMinDPhi170to175_->Fill(mTauTauMin, weight);
      else if ( dPhi > 160. ) hTauPairMassMinDPhi160to170_->Fill(mTauTauMin, weight);
      else if ( dPhi > 140. ) hTauPairMassMinDPhi140to160_->Fill(mTauTauMin, weight);
      else if ( dPhi > 120. ) hTauPairMassMinDPhi120to140_->Fill(mTauTauMin, weight);
      else                    hTauPairMassMinDPhiLt120_->Fill(mTauTauMin, weight);
    }

    hMt12MET_->Fill(diTauCandidate->mt12MET(), weight);

    hMt1MET_->Fill(diTauCandidate->mt1MET(), weight);
    hPt1MET_->Fill(diTauCandidate->pt1MET(), weight);
    hMt2MET_->Fill(diTauCandidate->mt2MET(), weight);
    hPt2MET_->Fill(diTauCandidate->pt2MET(), weight);
    if ( visMass > 0. ) {
      hMt2METvisMassRatio_->Fill(getBoundedValue(diTauCandidate->mt1MET()/visMass, -5.0, +5.0), weight);
      hMt2METvisMassRatio_->Fill(getBoundedValue(diTauCandidate->mt2MET()/visMass, -5.0, +5.0), weight);
    }

    if ( diTauCandidate->leg1().isNonnull() && diTauCandidate->leg2().isNonnull() && diTauCandidate->met().isNonnull() ) {
      hHt12MET_->Fill(diTauCandidate->leg1()->pt() + diTauCandidate->leg2()->pt() + diTauCandidate->met()->pt(), weight);
    }

    hDPhi12_->Fill(diTauCandidate->dPhi12(), weight);
    hDR12_->Fill(diTauCandidate->dR12(), weight);

    hVisEtaMin_->Fill(diTauCandidate->visEtaMin(), weight);
    hVisEtaMax_->Fill(diTauCandidate->visEtaMax(), weight);

    hDPhi1MET_->Fill(diTauCandidate->dPhi1MET(), weight);
    hDPhi2MET_->Fill(diTauCandidate->dPhi2MET(), weight);
    hDPhi1METvsDPhi2MET_->Fill(diTauCandidate->dPhi1MET(), diTauCandidate->dPhi2MET(), weight);

    hPzetaCorr_->Fill(diTauCandidate->pZetaVis(), diTauCandidate->pZeta(), weight);
    double pZetaDiff = diTauCandidate->pZeta() - 1.5*diTauCandidate->pZetaVis();
    hPzetaDiff_->Fill(pZetaDiff, weight);
    if ( diTauCandidate->hasMEtSignMatrix() && TMath::Abs(diTauCandidate->metSignMatrix().Determinant()) > 1.e-4 ) {
      double zetaX = TMath::Cos(diTauCandidate->leg1()->phi()) + TMath::Sin(diTauCandidate->leg2()->phi());
      double zetaY = TMath::Sin(diTauCandidate->leg1()->phi()) + TMath::Sin(diTauCandidate->leg2()->phi());
      double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
      if ( zetaR > 0. ) {
	zetaX /= zetaR;
	zetaY /= zetaR;
      }
      TVectorD zeta(2);
      zeta(0) = zetaX;
      zeta(1) = zetaY;
      double projMEtSign = zeta*(diTauCandidate->metSignMatrix()*zeta);
      if ( projMEtSign != 0. ) hPzetaDiffMEtSignRatio_->Fill(getBoundedValue(pZetaDiff/projMEtSign, -5.0, +5.0), weight);
      if ( dynamic_cast<const pat::MET*>(diTauCandidate->met().get()) != 0 ) {
	const pat::MET* patMEt = dynamic_cast<const pat::MET*>(diTauCandidate->met().get());
	if ( patMEt->genMET() ) {
	  pfMEtCovInverse_ = diTauCandidate->metSignMatrix();
	  pfMEtCovInverse_.Invert();
	  pfMEtCovInverseS_(0, 0) = pfMEtCovInverse_(0, 0);
	  pfMEtCovInverseS_(1, 0) = pfMEtCovInverse_(1, 0);
	  pfMEtCovInverseS_(0, 1) = pfMEtCovInverse_(0, 1);
	  pfMEtCovInverseS_(1, 1) = pfMEtCovInverse_(1, 1);

	  residualS_(0) = patMEt->genMET()->px() - patMEt->px();
	  residualS_(1) = patMEt->genMET()->py() - patMEt->py();

	  // CV: in case pfMEt resolution are accurately modeled,
	  //     residual*cov^-1*residual is expected to be distributed
	  //     like a chi^2 with 2 degrees of freedom
	  //    (cf. CMS AN-10/400, page 1)
	  double chi2 = ROOT::Math::Similarity(residualS_, pfMEtCovInverseS_);
	  double prob = TMath::Prob(chi2, 2);
	  hMEtSignProb_->Fill(getBoundedValue(prob, 0., 1.0), weight);

	  pfMEtCov_ = diTauCandidate->metSignMatrix();
	  pfMEtCov_.EigenVectors(pfMEtCovEigenvalues_);
	  hMEtSignEigenValue1_->Fill(pfMEtCovEigenvalues_(0), weight);
	  hMEtSignEigenValue2_->Fill(pfMEtCovEigenvalues_(1), weight);
	}
      }
    }

    hPzetaDiffVsDPhi12_->Fill(diTauCandidate->dPhi12(), diTauCandidate->pZeta() - 1.5*diTauCandidate->pZetaVis(), weight);
    hPzetaDiffVsMt1MET_->Fill(diTauCandidate->mt1MET(), diTauCandidate->pZeta() - 1.5*diTauCandidate->pZetaVis(), weight);
    if ( diTauCandidate->met().isNonnull() ) {
      hPzetaDiffVsMET_->Fill(diTauCandidate->met()->pt(), diTauCandidate->pZeta() - 1.5*diTauCandidate->pZetaVis(), weight);
      hPzetaDiffVsVisMass_->Fill(visMass, diTauCandidate->pZeta() - 1.5*diTauCandidate->pZetaVis(), weight);
      hMt1METvsMET_->Fill(diTauCandidate->met()->pt(), diTauCandidate->mt1MET(), weight);
      hMt1METvsVisMass_->Fill(visMass, diTauCandidate->mt1MET(), weight);
    }
  }
}

#include "DataFormats/Candidate/interface/Candidate.h" 
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

//typedef CompositePtrCandidateT1T2MEtHistManager<reco::Candidate, reco::Candidate> DiCandidatePairHistManager;
typedef CompositePtrCandidateT1T2MEtHistManager<pat::Electron, pat::Tau> PATElecTauPairHistManager;
typedef CompositePtrCandidateT1T2MEtHistManager<pat::Muon, pat::Tau> PATMuTauPairHistManager;
typedef CompositePtrCandidateT1T2MEtHistManager<pat::Tau, pat::Tau> PATDiTauPairHistManager;
typedef CompositePtrCandidateT1T2MEtHistManager<pat::Electron, pat::Muon> PATElecMuPairHistManager;
typedef CompositePtrCandidateT1T2MEtHistManager<pat::Electron, pat::Electron> PATDiElecPairHistManager;
typedef CompositePtrCandidateT1T2MEtHistManager<pat::Muon, pat::Muon> PATDiMuPairHistManager;

#include "FWCore/Framework/interface/MakerMacros.h"

//DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, DiCandidatePairHistManager, "DiCandidatePairHistManager");
//DEFINE_EDM_PLUGIN(HistManagerPluginFactory, DiCandidatePairHistManager, "DiCandidatePairHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATElecTauPairHistManager, "PATElecTauPairHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATElecTauPairHistManager, "PATElecTauPairHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATMuTauPairHistManager, "PATMuTauPairHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATMuTauPairHistManager, "PATMuTauPairHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATDiTauPairHistManager, "PATDiTauPairHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATDiTauPairHistManager, "PATDiTauPairHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATElecMuPairHistManager, "PATElecMuPairHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATElecMuPairHistManager, "PATElecMuPairHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATDiElecPairHistManager, "PATDiElecPairHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATDiElecPairHistManager, "PATDiElecPairHistManager");
DEFINE_EDM_PLUGIN(AnalyzerPluginFactory, PATDiMuPairHistManager, "PATDiMuPairHistManager");
DEFINE_EDM_PLUGIN(HistManagerPluginFactory, PATDiMuPairHistManager, "PATDiMuPairHistManager");
  
#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

//typedef HistManagerAdapter<DiCandidatePairHistManager> DiCandidatePairAnalyzer;
typedef HistManagerAdapter<PATElecTauPairHistManager> PATElecTauPairAnalyzer;
typedef HistManagerAdapter<PATMuTauPairHistManager> PATMuTauPairAnalyzer;
typedef HistManagerAdapter<PATDiTauPairHistManager> PATDiTauPairAnalyzer;
typedef HistManagerAdapter<PATElecMuPairHistManager> PATElecMuPairAnalyzer;
typedef HistManagerAdapter<PATDiElecPairHistManager> PATDiElecPairAnalyzer;
typedef HistManagerAdapter<PATDiMuPairHistManager> PATDiMuPairAnalyzer;

//DEFINE_FWK_MODULE(DiCandidatePairAnalyzer);
DEFINE_FWK_MODULE(PATElecTauPairAnalyzer);
DEFINE_FWK_MODULE(PATMuTauPairAnalyzer);
DEFINE_FWK_MODULE(PATDiTauPairAnalyzer);
DEFINE_FWK_MODULE(PATElecMuPairAnalyzer);
DEFINE_FWK_MODULE(PATDiElecPairAnalyzer);
DEFINE_FWK_MODULE(PATDiMuPairAnalyzer);

