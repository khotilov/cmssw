#include "TauAnalysis/RecoTools/interface/ZllRecoilCorrectionHistManager.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>

ZllRecoilCorrectionHistManager::ZllRecoilCorrectionHistManager(const edm::ParameterSet& cfg)
{}

ZllRecoilCorrectionHistManager::~ZllRecoilCorrectionHistManager()
{
  for ( std::vector<histogramsUvsQtNumVtxType*>::iterator it = histogramsUvsQtNumVtxBinned_.begin();
	it != histogramsUvsQtNumVtxBinned_.end(); ++it ) {
    delete (*it);
  }
}

void ZllRecoilCorrectionHistManager::bookHistograms(TFileDirectory& dir)
{
  histogramLepPlusPt_        = book1D(dir, "lPlusPt",            "P_{T}^{l+}",                            40,          0. ,         100.);
  histogramLepPlusEta_       = book1D(dir, "lPlusEta",           "#eta_{l+}",                             50,         -2.5,         +2.5);
  histogramLepPlusPhi_       = book1D(dir, "lPlusPhi",           "#phi_{l+}",                             36, -TMath::Pi(), +TMath::Pi());

  histogramLepMinusPt_       = book1D(dir, "lMinusPt",           "P_{T}^{l-}",                            40,          0. ,         100.);
  histogramLepMinusEta_      = book1D(dir, "lMinusEta",          "#eta_{l-}",                             50,         -2.5,         +2.5);
  histogramLepMinusPhi_      = book1D(dir, "lMinusPhi",          "#phi_{l-}",                             36, -TMath::Pi(), +TMath::Pi());

  histogramZllCandPt_        = book1D(dir, "ZllCandPt",          "P_{T}^{Z}",                             40,          0. ,         100.);
  histogramZllCandEta_       = book1D(dir, "ZllCandEta",         "#eta_{Z}",                              50,         -2.5,         +2.5);
  histogramZllCandPhi_       = book1D(dir, "ZllCandPhi",         "#phi_{Z}",                              36, -TMath::Pi(), +TMath::Pi());
  histogramZllCandMass_      = book1D(dir, "ZllCandMass",        "M(l+ l-)",                              60,         60. ,         120.);
  
  histogramMEtS_             = book1D(dir, "metS",               "E_{T}^{miss}",                          30,          0.0,         60.0);
  histogramMEtL_             = book1D(dir, "metL",               "E_{T}^{miss}",                          75,          0.0,        150.0);
  histogramMEtProjParlZ_     = book1D(dir, "metProjParlZ",       "E_{T}^{miss} Proj. parallel Z",         50,        -50.0,        +50.0);
  histogramMEtProjPerpZ_     = book1D(dir, "metProjPerpZ",       "E_{T}^{miss} Proj. perp. Z",            50,        -50.0,        +50.0);

  const int qTnumBins = 22;
  double qTbinning[qTnumBins + 1] = { 
    0., 2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30., 35., 40., 45., 50., 60., 70., 80., 100., 120., 150. 
  };
  histogramUparlDivQtVsQt_ = book2D(dir, "uParlDivQtVsQt", "u_{#parallel}/q_{T} vs q_{T}",                           
				    qTnumBins, qTbinning, 400,  -5.0,   +5.0);
  histogramUparlVsQt_      = book2D(dir, "uParlVsQt",      "u_{#parallel} vs q_{T}",                           
				    qTnumBins, qTbinning, 120, -50.0, +250.0);
  histogramUperpDivQtVsQt_ = book2D(dir, "uPerpDivQtVsQt", "u_{#perp}/q_{T} vs q_{T}",                           
				    qTnumBins, qTbinning, 400,  -5.0,   +5.0);
  histogramUperpVsQt_      = book2D(dir, "uPerpVsQt",      "u_{#perp} vs q_{T}",                           
				    qTnumBins, qTbinning,  40, -50.0,  +50.0);
  
  histogramsUvsQtNumVtxBinned_.push_back(new histogramsUvsQtNumVtxType(this, dir, qTnumBins, qTbinning, -1,  2));
  histogramsUvsQtNumVtxBinned_.push_back(new histogramsUvsQtNumVtxType(this, dir, qTnumBins, qTbinning,  3,  5));
  histogramsUvsQtNumVtxBinned_.push_back(new histogramsUvsQtNumVtxType(this, dir, qTnumBins, qTbinning,  6,  8));
  histogramsUvsQtNumVtxBinned_.push_back(new histogramsUvsQtNumVtxType(this, dir, qTnumBins, qTbinning,  9, 11));
  histogramsUvsQtNumVtxBinned_.push_back(new histogramsUvsQtNumVtxType(this, dir, qTnumBins, qTbinning, 12, -1));

  histogramVtxMultiplicity_  = book1D(dir, "numVertices",        "Num. Vertices",                         20,         -0.5,         19.5);
  histogramRhoNeutral_       = book1D(dir, "rhoNeutral",         "#rho_{neutral}",                        50,          0. ,          12.);
  histogramRhoChargedHadronNoPileUp_ = 
    book1D(dir, "rhoChargedHadronNoPileUp",                      "#rho_{h#pm}",                           50,          0. ,          12.);
  histogramRhoNeutralToChargedHadronNoPileUpRatio_ = 
    book1D(dir, "rhoNeutralToChargedHadronNoPileUpRatio",        "#rho_{neutral}/#rho_{h#pm}",            51,      -0.025 ,        2.525);
  histogramRhoPFNoPileUp_ = 
    book1D(dir, "rhoPFNoPileUp",                                 "#rho_{noPU}",                           50,          0. ,          20.);
}

void ZllRecoilCorrectionHistManager::fillHistograms(
       const reco::CompositeCandidate& ZllCand, const pat::MET& met, size_t vtxMultiplicity, double rhoNeutral, 
       double rhoChargedHadronNoPileUp, double rhoPFNoPileUp, double evtWeight)
{
  assert(ZllCand.numberOfDaughters() == 2);

  reco::Candidate::LorentzVector p4LepPlus;
  if      ( ZllCand.daughter(0)->charge() > +0.5 ) p4LepPlus = ZllCand.daughter(0)->p4();
  else if ( ZllCand.daughter(1)->charge() > +0.5 ) p4LepPlus = ZllCand.daughter(1)->p4();
  else assert(0);

  histogramLepPlusPt_->Fill(p4LepPlus.pt(), evtWeight);
  histogramLepPlusEta_->Fill(p4LepPlus.eta(), evtWeight); 
  histogramLepPlusPhi_->Fill(p4LepPlus.phi(), evtWeight); 

  reco::Candidate::LorentzVector p4LepMinus;
  if      ( ZllCand.daughter(0)->charge() < -0.5 ) p4LepMinus = ZllCand.daughter(0)->p4();
  else if ( ZllCand.daughter(1)->charge() < -0.5 ) p4LepMinus = ZllCand.daughter(1)->p4();
  else assert(0);

  histogramLepMinusPt_->Fill(p4LepMinus.pt(), evtWeight);
  histogramLepMinusEta_->Fill(p4LepMinus.eta(), evtWeight); 
  histogramLepMinusPhi_->Fill(p4LepMinus.phi(), evtWeight); 

  histogramZllCandPt_->Fill(ZllCand.pt(), evtWeight);
  histogramZllCandEta_->Fill(ZllCand.eta(), evtWeight);
  histogramZllCandPhi_->Fill(ZllCand.phi(), evtWeight);
  histogramZllCandMass_->Fill(ZllCand.mass(), evtWeight);
  
  histogramMEtS_->Fill(met.pt(), evtWeight);
  histogramMEtL_->Fill(met.pt(), evtWeight);

  if ( ZllCand.pt() > 0. ) {
    double qT = ZllCand.pt();
    double qX = ZllCand.px();
    double qY = ZllCand.py();
    
    double metPx = met.px();
    double metPy = met.py();

    double metProjParlZ = (qX*metPx + qX*metPy)/qT;
    double metProjPerpZ = (qX*metPy - qY*metPx)/qT;
    
    histogramMEtProjParlZ_->Fill(metProjParlZ, evtWeight);
    histogramMEtProjPerpZ_->Fill(metProjPerpZ, evtWeight);

    int errorFlag = 0;
    std::pair<double, double> uT = compMEtProjU(ZllCand.p4(), met.px(), met.py(), errorFlag);
    if ( !errorFlag ) {
      double uParl = uT.first;
      double uPerp = uT.second;
      if ( qT > 0. ) histogramUparlDivQtVsQt_->Fill(qT, uParl/qT);
      histogramUparlVsQt_->Fill(qT, uParl);
      if ( qT > 0. ) histogramUperpDivQtVsQt_->Fill(qT, uPerp/qT);
      histogramUperpVsQt_->Fill(qT, uPerp);

      for ( std::vector<histogramsUvsQtNumVtxType*>::iterator it = histogramsUvsQtNumVtxBinned_.begin();
	    it != histogramsUvsQtNumVtxBinned_.end(); ++it ) {
	if ( ((*it)->numVtxMin_ == -1 || (int)vtxMultiplicity >= (*it)->numVtxMin_) &&
	     ((*it)->numVtxMax_ == -1 || (int)vtxMultiplicity <= (*it)->numVtxMax_) ) {
	  if ( qT > 0. ) (*it)->histogramUparlDivQtVsQt_->Fill(qT, uParl/qT);
	  (*it)->histogramUparlVsQt_->Fill(qT, uParl);
	  if ( qT > 0. ) (*it)->histogramUperpDivQtVsQt_->Fill(qT, uPerp/qT);
	  (*it)->histogramUperpVsQt_->Fill(qT, uPerp);
	}
      }
    }
  }

  histogramVtxMultiplicity_->Fill(vtxMultiplicity, evtWeight);
  histogramRhoNeutral_->Fill(rhoNeutral, evtWeight);
  histogramRhoChargedHadronNoPileUp_->Fill(rhoChargedHadronNoPileUp, evtWeight);
  double rhoNeutralToChargedHadronNoPileUpRatio = ( rhoChargedHadronNoPileUp >= 0. ) ?
    (rhoNeutral/rhoChargedHadronNoPileUp) : -1.;
  histogramRhoNeutralToChargedHadronNoPileUpRatio_->Fill(rhoNeutralToChargedHadronNoPileUpRatio, evtWeight);
  histogramRhoPFNoPileUp_->Fill(rhoPFNoPileUp, evtWeight);
}

void ZllRecoilCorrectionHistManager::scaleHistograms(double factor)
{
  for ( std::vector<TH1*>::iterator histogram = histograms_.begin();
	histogram != histograms_.end(); ++histogram ) {
    if ( !(*histogram)->GetSumw2N() ) (*histogram)->Sumw2(); // CV: compute "proper" errors before scaling histogram
    (*histogram)->Scale(factor);
  }
}

TH1* ZllRecoilCorrectionHistManager::book1D(
       TFileDirectory& dir, const std::string& distribution, const std::string& title, int numBins, double min, double max)
{
  TH1* retVal = dir.make<TH1D>(distribution.data(), title.data(), numBins, min, max);
  histograms_.push_back(retVal);
  return retVal;
}

TH2* ZllRecoilCorrectionHistManager::book2D(
       TFileDirectory& dir, const std::string& distribution, const std::string& title, 
       int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax)				 
{
  TH2* retVal = dir.make<TH2D>(distribution.data(), title.data(), numBinsX, xMin, xMax, numBinsY, yMin, yMax);
  histograms_.push_back(retVal);
  return retVal;
}

TH2* ZllRecoilCorrectionHistManager::book2D(
       TFileDirectory& dir, const std::string& distribution, const std::string& title, 
       int numBinsX, double* xBinning, int numBinsY, double yMin, double yMax)				 
{
  TH2* retVal = dir.make<TH2D>(distribution.data(), title.data(), numBinsX, xBinning, numBinsY, yMin, yMax);
  histograms_.push_back(retVal);
  return retVal;
}
