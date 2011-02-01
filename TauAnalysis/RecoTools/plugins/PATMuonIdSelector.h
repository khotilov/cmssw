
/** \class PATMuonIdSelector
 *
 * Selection of "good quality" muons according to criteria
 * defined by Vector Boson Task Force and documented in Analysis Note CMS AN-10/264
 * (including extension to WW cross-section analysis documented in CMS AN-10/344)
 *
 * \author Michail Bachtis,
 *         Christian Veelken
 *
 * \version $Revision: 1.2 $
 *
 * $Id: PATMuonIdSelector.h,v 1.2 2010/09/28 11:23:36 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include <vector>

class PATMuonIdSelectorImp
{
 public:
  typedef pat::MuonCollection collection;

  explicit PATMuonIdSelectorImp(const edm::ParameterSet&);
  ~PATMuonIdSelectorImp();

  std::vector<const pat::Muon*>::const_iterator begin() const { return selected_.begin(); }
  std::vector<const pat::Muon*>::const_iterator end()   const { return selected_.end();   }

  void select(const edm::Handle<collection>&, edm::Event&, const edm::EventSetup&);
    
  size_t size() const { return selected_.size(); }

 private:
  void print();

  std::vector<const pat::Muon*> selected_;

//--- configuration parameters
  edm::InputTag srcBeamSpot_;
  edm::InputTag srcVertex_;

  bool     applyGlobalMuonPromptTight_;
  bool     applyAllArbitrated_;

  double   maxIPxy_;            // max. transverse   impact parameter of muon track
  double   maxIPz_;             // max. longitudinal impact parameter of muon track
  int      IPtrackType_;        // compute impact parameters for inner/global track of muon
  int      IPrefType_;          // compute impact parameters wrt. beam spot/reconstructed event vertex
  double   maxChi2red_;         // max. (normalized) chi^2 of global muon track fit per degree of freedom
  double   maxDptOverPt_;       // max. relative error on muon momentum (computed for inner track)
  unsigned minTrackerHits_;     // min. number of hits in SiStrip + Pixel detectors
  unsigned minPixelHits_;       // min. number of hits in Pixel detector
  unsigned minMuonChamberHits_; // min. number of hits in Muon chambers
  unsigned minMuonSegments_;    // min. number of segments in Muon stations matched to inner track
};
