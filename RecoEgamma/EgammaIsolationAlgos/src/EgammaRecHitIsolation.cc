//*****************************************************************************
// File:      EgammaRecHitIsolation.cc
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer, hacked by Sam Harper (ie the ugly stuff is mine)
// Institute: IIHE-VUB, RAL
//=============================================================================
//*****************************************************************************
//C++ includes
#include <vector>
#include <functional>

//ROOT includes
#include <Math/VectorUtil.h>

//CMSSW includes
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"

using namespace std;

EgammaRecHitIsolation::EgammaRecHitIsolation (double extRadius,
					      double intRadius,
					      double etaSlice,
					      double etLow,
					      double eLow,
					      edm::ESHandle<CaloGeometry> theCaloGeom,
					      CaloRecHitMetaCollectionV* caloHits,
					      DetId::Detector detector):  // not used anymore, kept for compatibility
  extRadius_(extRadius),
  intRadius_(intRadius),
  etaSlice_(etaSlice),
  etLow_(etLow),
  eLow_(eLow),
  theCaloGeom_(theCaloGeom) ,  
  caloHits_(caloHits)
{
  //set up the geometry and selector
  const CaloGeometry* caloGeom = theCaloGeom_.product();
  subdet_[0] = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
  subdet_[1] = caloGeom->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);
}

EgammaRecHitIsolation::~EgammaRecHitIsolation ()
{
}

double EgammaRecHitIsolation::getSum_(const reco::Candidate* emObject,bool returnEt) const
{

  double energySum = 0.;
  if (caloHits_){
    //Take the SC position
    reco::SuperClusterRef sc = emObject->get<reco::SuperClusterRef>();
    math::XYZPoint theCaloPosition = sc.get()->position();
    GlobalPoint pclu (theCaloPosition.x () ,
		      theCaloPosition.y () ,
		      theCaloPosition.z () );
    double etaclus = pclu.eta();
    double phiclus = pclu.phi();
    double r2 = intRadius_*intRadius_;

    for(int subdetnr=0; subdetnr<=1 ; subdetnr++){  // look in barrel and endcap
      CaloSubdetectorGeometry::DetIdSet chosen = subdet_[subdetnr]->getCells(pclu,extRadius_);// select cells around cluster
      CaloRecHitMetaCollectionV::const_iterator j=caloHits_->end();
      for (CaloSubdetectorGeometry::DetIdSet::const_iterator  i = chosen.begin ();i!= chosen.end ();++i){//loop selected cells

	j=caloHits_->find(*i); // find selected cell among rechits
	if( j!=caloHits_->end()){ // add rechit only if available 
	  const  GlobalPoint & position = theCaloGeom_.product()->getPosition(*i);
	  double eta = position.eta();
          double phi = position.phi();
          double etaDiff = eta - etaclus;
          double phiDiff= deltaPhi(phi,phiclus);
          double energy = j->energy();
	  
	  if ( fabs(etaDiff) < etaSlice_) continue;  // jurassic strip cut
	  if ( etaDiff*etaDiff + phiDiff*phiDiff < r2) continue; // jurassic exclusion cone cut

	  double et = energy*position.perp()/position.mag();
	  if ( et > etLow_ && fabs(energy) > eLow_){ //Changed energy --> fabs(energy)
	    if(returnEt) energySum+=et;
	    else energySum+=energy;
	  }
	}
	
      } 
    }
  }
  return energySum;
}

