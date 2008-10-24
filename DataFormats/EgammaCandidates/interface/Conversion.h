#ifndef EgammaCandidates_Conversion_h
#define EgammaCandidates_Conversion_h
/** \class reco::Conversion Conversion.h DataFormats/EgammaCandidates/interface/Conversion.h
 *
 * 
 *
 * \author N.Marinelli  University of Notre Dame, US
 *
 * \version $Id: Conversion.h,v 1.4 2008/05/08 19:42:45 nancy Exp $
 *
 */

#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"


namespace reco {
    class Conversion  {
  public:

      // Default constructor
      Conversion();
      
      Conversion( const reco::CaloClusterPtrVector clu, 
		  const std::vector<reco::TrackRef> tr,
		  const std::vector<math::XYZPoint> trackPositionAtEcal , 
		  const reco::Vertex  &  convVtx,
		  const std::vector<reco::CaloClusterPtr> & matchingBC,
		  const float DCA,        
		  const std::vector<math::XYZVector> & trackPin ,
		  const std::vector<math::XYZVector> & trackPout);
      
      
      
      
      /// destructor
      virtual ~Conversion();
      /// returns a clone of the candidate
      Conversion * clone() const;
      /// reference to a SuperCluster
      reco::CaloClusterPtrVector caloCluster() const ;
      /// vector of references to  tracks
      std::vector<reco::TrackRef> tracks() const ; 
      /// returns  the reco conversion vertex
      const reco::Vertex & conversionVertex() const  { return theConversionVertex_ ; }
      /// Bool flagging objects having track size >0
      bool isConverted() const;
      /// Number of tracks= 0,1,2
      unsigned int nTracks() const {return  tracks().size(); }
      /// if nTracks=2 returns the pair invariant mass
      double pairInvariantMass() const;
      /// Delta cot(Theta) where Theta is the angle in the (y,z) plane between the two tracks 
      double pairCotThetaSeparation() const;
      /// Conversion tracks momentum 
      GlobalVector  pairMomentum() const;
      /// Super Cluster energy divided by track pair momentum if Standard seeing method. If a pointer to two (or more clusters)
      /// is stored in the conversion, this method returns the energy sum of clusters divided by the  track pair momentum
      double EoverP() const;
      /// set primary event vertex used to define photon direction
      double zOfPrimaryVertexFromTracks() const;
      // Dist of minimum approach between tracks
      double distOfMinimumApproach() const {return  theMinDistOfApproach_;}
      // deltaPhi tracks at innermost point
      double dPhiTracksAtVtx() const;
      // deltaPhi tracks at ECAl
      double dPhiTracksAtEcal() const;
      // deltaEta tracks at ECAl
      double dEtaTracksAtEcal() const;
      
      ///// The following are variables provided per each track
      /// positions of the track extrapolation at the ECAL front face
      const std::vector<math::XYZPoint> & ecalImpactPosition() const  {return thePositionAtEcal_;}
      //  pair of BC matching a posteriori the tracks
      const std::vector<reco::CaloClusterPtr>&  bcMatchingWithTracks() const { return theMatchingBCs_;}
      /// signed transverse impact parameter for each track
      std::vector<double> tracksSigned_d0() const ;
      /// Vector of track momentum measured at the outermost hit
      const std::vector<math::XYZVector>& tracksPout() const {return theTrackPout_;}
      /// Vector of track momentum measured at the innermost hit
      const std::vector<math::XYZVector>& tracksPin() const  {return theTrackPin_;}
      
      
      
    private:
      
      /// vector pointer to a/multiple seed CaloCluster(s)
      reco::CaloClusterPtrVector caloCluster_;
      /// reference to a vector Track references
      std::vector<reco::TrackRef>  tracks_;
      /// position at the ECAl surface of the track extrapolation
      std::vector<math::XYZPoint>  thePositionAtEcal_;
      /// Fitted Kalman conversion vertex
      reco::Vertex theConversionVertex_;
      /// Clusters mathing the tracks (these are not the seeds)
      std::vector<reco::CaloClusterPtr> theMatchingBCs_;
      /// Distance of min approach of the two tracks
      float theMinDistOfApproach_;
      /// P_in of tracks
      std::vector<math::XYZVector> theTrackPin_;    
      /// P_out of tracks
      std::vector<math::XYZVector> theTrackPout_;    
      

  };
  
}

#endif
