#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/TrackReco/interface/Track.h" 
#include "DataFormats/CaloRecHit/interface/CaloCluster.h" 
#include "CLHEP/Units/GlobalPhysicalConstants.h"

// Temporary hack workaround for algoByName string "above array bounds"
// error reported by gcc from Conversion::algoByName below. This could
// should probably be demangled a bit in any case, but for now just
// turn off the warning/error in this file. 
// Do not copy this into other files without checking with the release
// coordinators!   PE 20091231
#if defined(__GNUC__) && __GNUC__ == 4 && __GNUC_MINOR__ == 4
# pragma GCC diagnostic ignored "-Warray-bounds"
#endif

using namespace reco;


Conversion::Conversion(  const reco::CaloClusterPtrVector sc, 
			 const std::vector<reco::TrackRef> tr, 
			 const std::vector<math::XYZPoint> trackPositionAtEcal, 
			 const reco::Vertex  & convVtx,
			 const std::vector<reco::CaloClusterPtr> & matchingBC,
                         const float DCA,
			 const std::vector<math::XYZPoint> & innPoint,
			 const std::vector<math::XYZVector> & trackPin,
			 const std::vector<math::XYZVector> & trackPout,
                         const float mva,
			 ConversionAlgorithm algo):  
			 

  caloCluster_(sc), tracks_(tr), 
  thePositionAtEcal_(trackPositionAtEcal), 
  theConversionVertex_(convVtx), 
  theMatchingBCs_(matchingBC), 
  theMinDistOfApproach_(DCA),
  theTrackInnerPosition_(innPoint),
  theTrackPin_(trackPin),
  theTrackPout_(trackPout),
  nSharedHits_(0),  
  theMVAout_(mva),
  algorithm_(algo),
  qualityMask_(0)
 { 
   
 }




Conversion::Conversion(  const reco::CaloClusterPtrVector sc, 
			 const std::vector<edm::RefToBase<reco::Track> > tr, 
			 const std::vector<math::XYZPoint> trackPositionAtEcal, 
			 const reco::Vertex  & convVtx,
			 const std::vector<reco::CaloClusterPtr> & matchingBC,
                         const float DCA,
			 const std::vector<math::XYZPoint> & innPoint,
			 const std::vector<math::XYZVector> & trackPin,
			 const std::vector<math::XYZVector> & trackPout,
                         const std::vector<uint8_t> nHitsBeforeVtx,                  
                         const std::vector<Measurement1DFloat> & dlClosestHitToVtx,
                         uint8_t nSharedHits,
                         const float mva,
			 ConversionAlgorithm algo):  
			 

  caloCluster_(sc), trackToBaseRefs_(tr), 
  thePositionAtEcal_(trackPositionAtEcal), 
  theConversionVertex_(convVtx), 
  theMatchingBCs_(matchingBC), 
  theMinDistOfApproach_(DCA),
  theTrackInnerPosition_(innPoint),
  theTrackPin_(trackPin),
  theTrackPout_(trackPout),
  nHitsBeforeVtx_(nHitsBeforeVtx),
  dlClosestHitToVtx_(dlClosestHitToVtx),
  nSharedHits_(nSharedHits),
  theMVAout_(mva),
  algorithm_(algo),
  qualityMask_(0)
 { 
   
 }




Conversion::Conversion(  const reco::CaloClusterPtrVector sc, 
			 const std::vector<reco::TrackRef> tr, 
			 const reco::Vertex  & convVtx,
			 ConversionAlgorithm algo):  
  caloCluster_(sc), tracks_(tr), 
  theConversionVertex_(convVtx),
  nSharedHits_(0),
  algorithm_(algo),
  qualityMask_(0)
 { 


  theMinDistOfApproach_ = 9999.;
  theMVAout_ = 9999.;
  thePositionAtEcal_.push_back(math::XYZPoint(0.,0.,0.));
  thePositionAtEcal_.push_back(math::XYZPoint(0.,0.,0.));
  theTrackInnerPosition_.push_back(math::XYZPoint(0.,0.,0.));
  theTrackInnerPosition_.push_back(math::XYZPoint(0.,0.,0.));
  theTrackPin_.push_back(math::XYZVector(0.,0.,0.));
  theTrackPin_.push_back(math::XYZVector(0.,0.,0.));
  theTrackPout_.push_back(math::XYZVector(0.,0.,0.));
  theTrackPout_.push_back(math::XYZVector(0.,0.,0.));

   
 }


Conversion::Conversion(  const reco::CaloClusterPtrVector sc, 
			 const std::vector<edm::RefToBase<reco::Track> >  tr, 
			 const reco::Vertex  & convVtx,
			 ConversionAlgorithm algo):  
  caloCluster_(sc), trackToBaseRefs_(tr), 
  theConversionVertex_(convVtx), 
  nSharedHits_(0),
  algorithm_(algo),
  qualityMask_(0)
 { 


  theMinDistOfApproach_ = 9999.;
  theMVAout_ = 9999.;
  thePositionAtEcal_.push_back(math::XYZPoint(0.,0.,0.));
  thePositionAtEcal_.push_back(math::XYZPoint(0.,0.,0.));
  theTrackInnerPosition_.push_back(math::XYZPoint(0.,0.,0.));
  theTrackInnerPosition_.push_back(math::XYZPoint(0.,0.,0.));
  theTrackPin_.push_back(math::XYZVector(0.,0.,0.));
  theTrackPin_.push_back(math::XYZVector(0.,0.,0.));
  theTrackPout_.push_back(math::XYZVector(0.,0.,0.));
  theTrackPout_.push_back(math::XYZVector(0.,0.,0.));

   
 }



Conversion::Conversion() { 

  algorithm_=0;
  qualityMask_=0;
  theMinDistOfApproach_ = 9999.;
  nSharedHits_ = 0;
  theMVAout_ = 9999.;
  thePositionAtEcal_.push_back(math::XYZPoint(0.,0.,0.));
  thePositionAtEcal_.push_back(math::XYZPoint(0.,0.,0.));
  theTrackInnerPosition_.push_back(math::XYZPoint(0.,0.,0.));
  theTrackInnerPosition_.push_back(math::XYZPoint(0.,0.,0.));
  theTrackPin_.push_back(math::XYZVector(0.,0.,0.));
  theTrackPin_.push_back(math::XYZVector(0.,0.,0.));
  theTrackPout_.push_back(math::XYZVector(0.,0.,0.));
  theTrackPout_.push_back(math::XYZVector(0.,0.,0.));
    
}


Conversion::~Conversion() { }


std::string const Conversion::algoNames[] = { "undefined","ecalSeeded","trackerOnly","mixed"};  

Conversion::ConversionAlgorithm Conversion::algoByName(const std::string &name){
  ConversionAlgorithm size = algoSize;
  int index = std::find(algoNames, algoNames+size, name)-algoNames;
  if(index == size) return undefined;

  return ConversionAlgorithm(index);
}

Conversion * Conversion::clone() const { 
  return new Conversion( * this ); 
}

reco::CaloClusterPtrVector Conversion::caloCluster() const {
  return caloCluster_;
}



std::vector<edm::RefToBase<reco::Track> >  Conversion::tracks() const { 
  if (trackToBaseRefs_.size() ==0 ) {
 
    for (std::vector<reco::TrackRef>::const_iterator ref=tracks_.begin(); ref!=tracks_.end(); ref++ ) 
      {
	edm::RefToBase<reco::Track> tt(*ref);
	trackToBaseRefs_.push_back(tt);
	
      }  
  }

  return trackToBaseRefs_;
}



bool Conversion::isConverted() const {
  
  if ( this->nTracks() == 2 ) 
    return true;
  else
    return false;
}


double  Conversion::zOfPrimaryVertexFromTracks()  const  {
  double theZOfPrimaryVertexFromTracks=-9999.;

  float pTrkMag=this->pairMomentum().mag();
  
  if ( pTrkMag>0 && this->conversionVertex().isValid() && sqrt(this->conversionVertex().position().perp2()) !=0 ) {
    float theta=acos(this->pairMomentum().z() /pTrkMag);
    theZOfPrimaryVertexFromTracks = this->conversionVertex().position().z()  - sqrt(this->conversionVertex().position().perp2())*(1./tan(theta));
    
  }

  return  theZOfPrimaryVertexFromTracks;

}


double Conversion::pairInvariantMass() const{
  double invMass=-99.;
  const float mElec= 0.000511;
  if ( nTracks()==2 ) {
    double px= tracksPin()[0].x() +  tracksPin()[1].x();
    double py= tracksPin()[0].y() +  tracksPin()[1].y();
    double pz= tracksPin()[0].z() +  tracksPin()[1].z();
    double mom1= tracksPin()[0].Mag2();
    double mom2= tracksPin()[1].Mag2();
    double e = sqrt( mom1+ mElec*mElec ) + sqrt( mom2 + mElec*mElec );
    invMass= ( e*e - px*px -py*py - pz*pz);
    if ( invMass>0) invMass = sqrt(invMass);
    else 
      invMass = -1;
  }
  
  return invMass;
}

double  Conversion::pairCotThetaSeparation() const  {
  double dCotTheta=-99.;
  
  if ( nTracks()==2 ) {
    double theta1=tracksPin()[0].Theta();
    double theta2=tracksPin()[1].Theta();
    dCotTheta =  1./tan(theta1) - 1./tan(theta2) ;
  }

  return dCotTheta;

}


GlobalVector  Conversion::pairMomentum() const  {

  double px=0.;
  double py=0.;
  double pz=0.;
  
  if ( nTracks()==2 ) {
    px= tracksPin()[0].x() +  tracksPin()[1].x();
    py= tracksPin()[0].y() +  tracksPin()[1].y();
    pz= tracksPin()[0].z() +  tracksPin()[1].z();

  } else if (  nTracks()==1 ) {
    px= tracksPin()[0].x();
    py= tracksPin()[0].y();
    pz= tracksPin()[0].z();
  }

  GlobalVector momTracks(px,py,pz);
  return momTracks;


}


math::XYZTLorentzVectorD Conversion::refittedPair4Momentum() const  {

  math::XYZTLorentzVectorD p4;
  if ( this->conversionVertex().isValid() ) 
    p4 = this->conversionVertex().p4( 0.000511, 0.5);

  return p4;


}


GlobalVector  Conversion::refittedPairMomentum() const  {

  double px=0.;
  double py=0.;
  double pz=0.;
  
  if (  this->conversionVertex().isValid() ) {
    px= this->refittedPair4Momentum().px();
    py= this->refittedPair4Momentum().py();
    pz= this->refittedPair4Momentum().pz();
  
  }

  GlobalVector momTracks(px,py,pz);
  return momTracks;


}



double  Conversion::EoverP() const  {


  double ep=-99.;

  if ( nTracks() > 0  ) {
    unsigned int size= this->caloCluster().size();
    float etot=0.;
    for ( unsigned int i=0; i<size; i++) {
      etot+= caloCluster()[i]->energy();
    }
    if (this->pairMomentum().mag() !=0) ep= etot/this->pairMomentum().mag();
  }



  return ep;  

}



double  Conversion::EoverPrefittedTracks() const  {


  double ep=-99.;
 
  if ( nTracks() > 0  ) {
    unsigned int size= this->caloCluster().size();
    float etot=0.;
    for ( unsigned int i=0; i<size; i++) {
      etot+= caloCluster()[i]->energy();
    }
    if (this->refittedPairMomentum().mag() !=0) ep= etot/this->refittedPairMomentum().mag();
  }



  return ep;  

}

 

std::vector<double>  Conversion::tracksSigned_d0() const  {
  std::vector<double> result;

  for (unsigned int i=0; i< nTracks(); i++)   
    result.push_back(tracks()[i]->d0()* tracks()[i]->charge()) ;

  return result;


}

double  Conversion::dPhiTracksAtVtx() const  {
  double result=-99;
  if  ( nTracks()==2 ) {
    float phiTk1=  tracksPin()[0].phi();
    float phiTk2=  tracksPin()[1].phi();
    result = phiTk1-phiTk2;
    
    if( result   > pi)  { result = result - twopi;}
    if( result  < -pi)  { result = result + twopi;}

  }

  return result;


}

double  Conversion::dPhiTracksAtEcal() const  {
  double result =-99.;
  
  if (  nTracks()==2  && bcMatchingWithTracks()[0].isNonnull() && bcMatchingWithTracks()[1].isNonnull() ) {
    
    float recoPhi1 = ecalImpactPosition()[0].phi();
    if( recoPhi1   > pi)  { recoPhi1 = recoPhi1 - twopi;}
    if( recoPhi1  < -pi)  { recoPhi1 = recoPhi1 + twopi;}

    float recoPhi2 = ecalImpactPosition()[1].phi();
    if( recoPhi2   > pi)  { recoPhi2 = recoPhi2 - twopi;}
    if( recoPhi2  < -pi)  { recoPhi2 = recoPhi2 + twopi;}

    result = recoPhi1 -recoPhi2;

    if( result   > pi)  { result = result - twopi;}
    if( result  < -pi)  { result = result + twopi;}

  }

  return result;


}

double  Conversion::dEtaTracksAtEcal() const  {
  double result=-99.;


  if ( nTracks()==2 && bcMatchingWithTracks()[0].isNonnull() && bcMatchingWithTracks()[1].isNonnull() ) {

    result =ecalImpactPosition()[0].eta() - ecalImpactPosition()[1].eta();

  }



  return result;


}

