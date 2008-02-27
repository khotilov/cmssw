#ifndef RecoTauTag_RecoTau_CaloRecoTauAlgorithm_H
#define RecoTauTag_RecoTau_CaloRecoTauAlgorithm_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h" 
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h" 
#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/CaloTauTagInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTowerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "RecoTauTag/TauTagTools/interface/CaloTauElementsOperators.h"
#include "RecoTauTag/TauTagTools/interface/TauTagTools.h"

#include "RecoJets/JetAlgorithms/interface/JetMatchingTools.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

using namespace std;
using namespace reco;
using namespace edm;

class  CaloRecoTauAlgorithm  {
 public:
  CaloRecoTauAlgorithm();  
  CaloRecoTauAlgorithm(const ParameterSet& iConfig);
  ~CaloRecoTauAlgorithm(){}
  void setTransientTrackBuilder(const TransientTrackBuilder*);
  void setMagneticField(const MagneticField*);
  CaloTau buildCaloTau(Event&,const EventSetup&,const CaloTauTagInfoRef&,const Vertex&); 
 private:
  vector<CaloTowerDetId> getCaloTowerneighbourDetIds(const CaloSubdetectorGeometry*,CaloTowerDetId);
  const TransientTrackBuilder* TransientTrackBuilder_;
  const MagneticField* MagneticField_;
  double LeadTrack_minPt_;
  double Track_minPt_;
  bool UseTrackLeadTrackDZconstraint_;
  double TrackLeadTrack_maxDZ_;
  double ECALRecHit_minEt_;
  string MatchingConeMetric_;
  string MatchingConeSizeFormula_;
  double MatchingConeSize_min_;
  double MatchingConeSize_max_;
  string TrackerSignalConeMetric_;
  string TrackerSignalConeSizeFormula_;
  double TrackerSignalConeSize_min_;
  double TrackerSignalConeSize_max_;
  string TrackerIsolConeMetric_;
  string TrackerIsolConeSizeFormula_;
  double TrackerIsolConeSize_min_;
  double TrackerIsolConeSize_max_;
  string ECALSignalConeMetric_;
  string ECALSignalConeSizeFormula_;
  double ECALSignalConeSize_min_;
  double ECALSignalConeSize_max_;
  string ECALIsolConeMetric_;
  string ECALIsolConeSizeFormula_;
  double ECALIsolConeSize_min_;
  double ECALIsolConeSize_max_;
  double AreaMetric_recoElements_maxabsEta_;
  const double chargedpi_mass_; //PDG Particle Physics Booklet, 2004
};
#endif 

