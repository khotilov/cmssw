#ifndef TauAnalysis_TauIdEfficiency_ZllRecoilCorrectionHistManager_h
#define TauAnalysis_TauIdEfficiency_ZllRecoilCorrectionHistManager_h

/** \class ZllRecoilCorrectionHistManager
 *
 * Fill Z-recoil correction control plots
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.9 $
 *
 * $Id: ZllRecoilCorrectionHistManager.h,v 1.9 2012/02/02 10:18:27 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

class ZllRecoilCorrectionHistManager
{

 public:
  /// constructor
  ZllRecoilCorrectionHistManager(edm::ParameterSet const&);

  /// destructor
  virtual ~ZllRecoilCorrectionHistManager();

  /// book and fill histograms
  void bookHistograms(TFileDirectory&);
  void fillHistograms(const reco::CompositeCandidate&, const std::vector<pat::Muon>&, 
		      const std::vector<pat::Jet>&, const pat::MET&, 
		      const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, 
		      const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, 
		      const reco::Candidate::LorentzVector&,
		      int, int, int, const reco::VertexCollection&, double, double);
  
  /// scale all bin-contents/bin-errors by factor given as function argument
  /// (to account for events lost, due to aborted skimming/crab or PAT-tuple production/lxbatch jobs)
  void scaleHistograms(double);

 protected:

  TH1* book1D(TFileDirectory&, const std::string&, const std::string&, int, double, double);
  TH2* book2D(TFileDirectory&, const std::string&, const std::string&, int, double, double, int, double, double);
  TH2* book2D(TFileDirectory&, const std::string&, const std::string&, int, double*, int, double, double);

 private:

  TH1* histogramLepPlusPt_;
  TH1* histogramLepPlusEta_;
  TH1* histogramLepPlusPhi_;

  TH1* histogramLepMinusPt_;
  TH1* histogramLepMinusEta_;
  TH1* histogramLepMinusPhi_;

  TH1* histogramZllCandPt_;
  TH1* histogramZllCandEta_;
  TH1* histogramZllCandPhi_;
  TH1* histogramZllCandMass_;

  TH1* histogramNumJetsRawPtGt10_;
  TH1* histogramNumJetsCorrPtGt10_;
  TH1* histogramNumJetsRawPtGt12_;
  TH1* histogramNumJetsCorrPtGt12_;
  TH1* histogramNumJetsRawPtGt15_;
  TH1* histogramNumJetsCorrPtGt15_;
  TH1* histogramNumJetsRawPtGt17_;
  TH1* histogramNumJetsCorrPtGt17_;
  TH1* histogramNumJetsRawPtGt20_;
  TH1* histogramNumJetsCorrPtGt20_;

  TH1* histogramJetPtAbsEtaLt11_;
  TH1* histogramJetResAbsEtaLt11_;
  TH1* histogramJetPtAbsEta11to17_;
  TH1* histogramJetResAbsEta11to17_;
  TH1* histogramJetPtAbsEta17to23_;
  TH1* histogramJetResAbsEta17to23_;
  TH1* histogramJetPtAbsEtaGt23_;
  TH1* histogramJetResAbsEtaGt23_;

  TH1* histogramMEtS_;
  TH1* histogramMEtL_;
  TH1* histogramMEtProjParlZ_;
  TH1* histogramMEtProjPerpZ_;
  TH1* histogramMEtX_;
  TH1* histogramMEtY_;

  TH1* histogramUparl_;
  TH1* histogramUperp_;

  TH2* histogramUparlDivQtVsQt_;
  TH2* histogramUparlVsQt_;
  TH2* histogramUperpDivQtVsQt_;
  TH2* histogramUperpVsQt_;
  TH1* histogramQt_;

  struct histogramsUvsQtNumObjType
  {
    histogramsUvsQtNumObjType(ZllRecoilCorrectionHistManager* histManager, 
			      TFileDirectory& dir, int qTnumBins, double* qTbinning, 
			      const std::string& numObjLabel, int numObjMin, int numObjMax)
      : numObjMin_(numObjMin),
	numObjMax_(numObjMax)
    {
      TString label = numObjLabel.data();
      if      ( numObjMin_ == -1         ) label.Append(Form("Le%i",   numObjMax_));
      else if ( numObjMax_ == -1         ) label.Append(Form("Ge%i",   numObjMin_));
      else if ( numObjMin_ == numObjMax_ ) label.Append(Form("Eq%i",   numObjMin_));
      else                                 label.Append(Form("%ito%i", numObjMin_, numObjMax_));

      histogramUparlDivQtVsQt_ = histManager->book2D(
	dir, TString("uParlDivQtVsQt").Append(label).Data(), "u_{#parallel}/q_{T} vs q_{T}", 
	qTnumBins, qTbinning, 400,  -5.0,   +5.0);
      histogramUparlVsQt_ = histManager->book2D(
	dir, TString("uParlVsQt").Append(label).Data(),      "u_{#parallel} vs q_{T}",   
	qTnumBins, qTbinning, 230, -500.0,  +75.0);
      histogramUperpDivQtVsQt_ = histManager->book2D(
        dir, TString("uPerpDivQtVsQt").Append(label).Data(), "u_{#perp}/q_{T} vs q_{T}", 
	qTnumBins, qTbinning, 400,  -5.0,   +5.0);
      histogramUperpVsQt_ = histManager->book2D(
        dir, TString("uPerpVsQt").Append(label).Data(),      "u_{#perp} vs q_{T}",  
	qTnumBins, qTbinning,  60, -75.0,  +75.0);
      histogramQt_ = histManager->book1D(
        dir, TString("qT").Append(label).Data(),             "q_{T}",  
	600, 0., 300.);
    }
    ~histogramsUvsQtNumObjType() {}
    int numObjMin_;
    int numObjMax_;

    TH2* histogramUparlDivQtVsQt_;
    TH2* histogramUparlVsQt_;
    TH2* histogramUperpDivQtVsQt_;
    TH2* histogramUperpVsQt_;
    TH1* histogramQt_;
  };

  std::vector<histogramsUvsQtNumObjType*> histogramsUvsQtNumVtxBinned_;

  std::vector<histogramsUvsQtNumObjType*> histogramsUvsQtNumJetsBinned_;

  struct histogramsMEtPhiAsymmetryVsQtType
  {
    histogramsMEtPhiAsymmetryVsQtType(ZllRecoilCorrectionHistManager* histManager, 
				      TFileDirectory& dir, 
				      double qTmin, double qTmax)
      : qTmin_(qTmin),
	qTmax_(qTmax)
    {
      TString label;
      if      ( qTmin_ <= 0. ) label.Append(Form("_qTle%1.0f",      qTmax_));
      else if ( qTmax_ <= 0. ) label.Append(Form("_qTge%1.0f",      qTmin_));
      else                     label.Append(Form("_qT%1.0fto%1.0f", qTmin_, qTmax_));
      label = label.ReplaceAll(".", "_");

      histogramMEtVsPhiZ_ = histManager->book2D(
	dir, TString("metVsPhiZ").Append(label).Data(), "MET vs #phi_{Z}", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);
      histogramMEtVsPhiZforEtaHadGt0_ = histManager->book2D(
	dir, TString("metVsPhiZforEtaHadGt0").Append(label).Data(), "MET vs #phi_{Z} (#eta_{had} > 0)", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);
      histogramMEtVsPhiZforEtaHadLt0_ = histManager->book2D(
	dir, TString("metVsPhiZforEtaHadLt0").Append(label).Data(), "MET vs #phi_{Z} (#eta_{had} < 0)", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);

      histogramSumPFChargedHadronsVsPhiZ_ = histManager->book2D(
	dir, TString("sumPFChargedHadronsVsPhiZ").Append(label).Data(), "#Sigma #vec{P}_{T}^{h} vs #phi_{Z}", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);
      histogramSumPFChargedHadronsVsPhiZforEtaHadGt0_ = histManager->book2D(
	dir, TString("sumPFChargedHadronsVsPhiZforEtaHadGt0").Append(label).Data(), "#Sigma #vec{P}_{T}^{h} vs #phi_{Z} (#eta_{had} > 0)", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);
      histogramSumPFChargedHadronsVsPhiZforEtaHadLt0_ = histManager->book2D(
	dir, TString("sumPFChargedHadronsVsPhiZforEtaHadLt0").Append(label).Data(), "#Sigma #vec{P}_{T}^{h} vs #phi_{Z} (#eta_{had} < 0)", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);

      histogramSumPFNeutralHadronsVsPhiZ_ = histManager->book2D(
	dir, TString("sumPFNeutralHadronsVsPhiZ").Append(label).Data(), "#Sigma #vec{P}_{T}^{h0} vs #phi_{Z}", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);
      histogramSumPFNeutralHadronsVsPhiZforEtaHadGt0_ = histManager->book2D(
	dir, TString("sumPFNeutralHadronsVsPhiZforEtaHadGt0").Append(label).Data(), "#Sigma #vec{P}_{T}^{h0} vs #phi_{Z} (#eta_{had} > 0)", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);
      histogramSumPFNeutralHadronsVsPhiZforEtaHadLt0_ = histManager->book2D(
	dir, TString("sumPFNeutralHadronsVsPhiZforEtaHadLt0").Append(label).Data(), "#Sigma #vec{P}_{T}^{h0} vs #phi_{Z} (#eta_{had} < 0)", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);

      histogramSumPFGammasVsPhiZ_ = histManager->book2D(
	dir, TString("sumPFGammasVsPhiZ").Append(label).Data(), "#Sigma #vec{P}_{T}^{#gamma} vs #phi_{Z}", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);
      histogramSumPFGammasVsPhiZforEtaHadGt0_ = histManager->book2D(
	dir, TString("sumPFGammasVsPhiZforEtaHadGt0").Append(label).Data(), "#Sigma #vec{P}_{T}^{#gamma} vs #phi_{Z} (#eta_{had} > 0)", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);
      histogramSumPFGammasVsPhiZforEtaHadLt0_ = histManager->book2D(
	dir, TString("sumPFGammasVsPhiZforEtaHadLt0").Append(label).Data(), "#Sigma #vec{P}_{T}^{#gamma} vs #phi_{Z} (#eta_{had} < 0)", 
	36, -TMath::Pi(), +TMath::Pi(), 60, 0., 60.);
    }
    ~histogramsMEtPhiAsymmetryVsQtType() {}
    double qTmin_;
    double qTmax_;

    TH2* histogramMEtVsPhiZ_;
    TH2* histogramMEtVsPhiZforEtaHadGt0_;
    TH2* histogramMEtVsPhiZforEtaHadLt0_;

    TH2* histogramSumPFChargedHadronsVsPhiZ_;
    TH2* histogramSumPFChargedHadronsVsPhiZforEtaHadGt0_;
    TH2* histogramSumPFChargedHadronsVsPhiZforEtaHadLt0_;
    TH2* histogramSumPFNeutralHadronsVsPhiZ_;
    TH2* histogramSumPFNeutralHadronsVsPhiZforEtaHadGt0_;
    TH2* histogramSumPFNeutralHadronsVsPhiZforEtaHadLt0_;
    TH2* histogramSumPFGammasVsPhiZ_;
    TH2* histogramSumPFGammasVsPhiZforEtaHadGt0_;
    TH2* histogramSumPFGammasVsPhiZforEtaHadLt0_;
  };

  std::vector<histogramsMEtPhiAsymmetryVsQtType*> histogramsMEtPhiAsymmetryVsQtBinned_;

  TH1* histogramNumPU_bxMinus1_;
  TH1* histogramNumPU_bx0_;
  TH1* histogramNumPU_bxPlus1_;

  TH1* histogramVtxMultiplicity_;
  TH1* histogramVtxZ_;
  TH1* histogramVtxRho_;
  TH1* histogramRhoNeutral_;
  
  std::vector<TH1*> histograms_;
};

#endif
