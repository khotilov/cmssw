#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/TrackReco/interface/Track.h" 
#include "DataFormats/EgammaReco/interface/SuperCluster.h" 
#include "DataFormats/EgammaReco/interface/ClusterShape.h" 

using namespace reco;

Conversion::Conversion(  const reco::SuperClusterRef sc, 
                                   const std::vector<reco::TrackRef> tr, 
			           const std::vector<math::XYZPoint> trackPositionAtEcal, 
				   const reco::Vertex  & convVtx,
				   const std::vector<reco::BasicClusterRef> & matchingBC):  

  superCluster_(sc), tracks_(tr), 
  thePositionAtEcal_(trackPositionAtEcal), 
  theConversionVertex_(convVtx), 
  theMatchingBCs_(matchingBC)  {
  
}


Conversion::~Conversion() { }


Conversion * Conversion::clone() const { 
  return new Conversion( * this ); 
}

reco::SuperClusterRef Conversion::superCluster() const {
  return superCluster_;
}



std::vector<reco::TrackRef>  Conversion::tracks() const { 
   return tracks_;
}



bool Conversion::isConverted() const {
  
  if ( this->nTracks() > 0) 
    return true;
  else
    return false;
}


double  Conversion::makePrimaryVertexZ() const  {
  double theZOfPrimaryVertexFromTracks=-9999.;

  float pTrkMag=this->pairMomentum().mag();
  
  if ( pTrkMag>0 && sqrt(this->conversionVertex().position().perp2()) !=0 ) {
    float theta=acos(this->pairMomentum().z() /pTrkMag);
    theZOfPrimaryVertexFromTracks = this->conversionVertex().position().z()  - sqrt(this->conversionVertex().position().perp2())*(1./tan(theta));
    
  }

  return  theZOfPrimaryVertexFromTracks;

}


double Conversion::makePairInvariantMass() const{
  double invMass=-99.;
  const float mElec= 0.000511;
  if ( nTracks()==2 ) {
    double px= tracks()[0]->innerMomentum().x() + tracks()[1]->innerMomentum().x();
    double py= tracks()[0]->innerMomentum().y() + tracks()[1]->innerMomentum().y();
    double pz= tracks()[0]->innerMomentum().z() + tracks()[1]->innerMomentum().z();
    double mom1=tracks()[0]->innerMomentum().Mag2() ;
    double mom2=tracks()[1]->innerMomentum().Mag2() ;
    double e = sqrt( mom1+ mElec*mElec ) + sqrt( mom2 + mElec*mElec );
    invMass= ( e*e - px*px -py*py - pz*pz);
  }
  
  return invMass;
}

double  Conversion::makePairCotThetaSeparation() const  {
  double dCotTheta=-99.;
  
  if ( nTracks()==2 ) {
    double theta1=tracks()[0]->innerMomentum().Theta();
    double theta2=tracks()[1]->innerMomentum().Theta();
    dCotTheta =  1./tan(theta1) - 1./tan(theta2) ;
  }

  return dCotTheta;

}

GlobalVector  Conversion::makePairMomentum() const  {

  double px=0.;
  double py=0.;
  double pz=0.;
  
  if ( nTracks()==2 ) {
    px= tracks()[0]->innerMomentum().x() + tracks()[1]->innerMomentum().x();
    py= tracks()[0]->innerMomentum().y() + tracks()[1]->innerMomentum().y();
    pz= tracks()[0]->innerMomentum().z() + tracks()[1]->innerMomentum().z();

  } else if (  nTracks()==1 ) {
    px= tracks()[0]->innerMomentum().x() ;
    py= tracks()[0]->innerMomentum().y() ;
    pz= tracks()[0]->innerMomentum().z() ;
  }

  GlobalVector momTracks(px,py,pz);
  return momTracks;


}

double  Conversion::makePairMomentumEta() const {

  double etaTracks=-99.;
  
  if ( nTracks() > 0 ) {
    etaTracks= pairMomentum().eta();
  } 

  return etaTracks;
}

double  Conversion::makePairMomentumPhi() const  {

  double phiTracks=-99.;
  
  if ( nTracks() > 0 ) {
    phiTracks= pairMomentum().phi();
  } 

  return phiTracks;

}

double  Conversion::makePairPtOverEtSC() const  {

  double ptOverEtSC=-99.;
  double etaSC=superCluster()->eta();
  double EtSC= superCluster()->energy()/cosh(etaSC);
  
   if ( nTracks() ==2 ) {
    ptOverEtSC= (pairMomentum().perp()) /EtSC; 
   }

   return ptOverEtSC;

}

double  Conversion::makeEoverP() const  {

  double ep=-99.;
  if ( nTracks() > 0  ) 
    ep=  this->superCluster()->energy()/this->pairMomentum().mag();

  return ep;  

}

 
