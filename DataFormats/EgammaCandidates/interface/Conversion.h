#ifndef EgammaCandidates_Conversion_h
#define EgammaCandidates_Conversion_h
/** \class reco::Conversion Conversion.h DataFormats/EgammaCandidates/interface/Conversion.h
 *
 * 
 *
 * \author N.Marinelli  University of Notre Dame, US
 *
 * \version $Id: Conversion.h,v 1.13 2010/09/23 16:59:09 nancy Exp $
 *
 */

#include <bitset>
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"


namespace reco {
    class Conversion  {
  public:

      enum ConversionAlgorithm {undefined=0, 
				ecalSeeded=1, 
				trackerOnly=2, 
				mixed=3, 
				algoSize=4}; 

      enum ConversionQuality {generalTracksOnly=0, 
			      arbitratedEcalSeeded=1, 
			      arbitratedMerged=2, 
			      highPurity=8, 
			      highEfficiency=9,
			      ecalMatched1Track=10,
			      ecalMatched2Track=11};

      static const std::string algoNames[];      

      // Default constructor
      Conversion();
      
      Conversion( const reco::CaloClusterPtrVector clu, 
		  const std::vector<edm::RefToBase<reco::Track> > tr,
		  const std::vector<math::XYZPoint> trackPositionAtEcal , 
		  const reco::Vertex  &  convVtx,
		  const std::vector<reco::CaloClusterPtr> & matchingBC,
		  const float DCA,        
		  const std::vector<math::XYZPoint> & innPoint,
		  const std::vector<math::XYZVector> & trackPin ,
		  const std::vector<math::XYZVector> & trackPout,
                  const float mva,
		  ConversionAlgorithm=undefined);


      Conversion( const reco::CaloClusterPtrVector clu, 
		  const std::vector<reco::TrackRef> tr,
		  const std::vector<math::XYZPoint> trackPositionAtEcal , 
		  const reco::Vertex  &  convVtx,
		  const std::vector<reco::CaloClusterPtr> & matchingBC,
		  const float DCA,        
		  const std::vector<math::XYZPoint> & innPoint,
		  const std::vector<math::XYZVector> & trackPin ,
		  const std::vector<math::XYZVector> & trackPout,
                  const float mva,
		  ConversionAlgorithm=undefined);
      


      
      Conversion( const reco::CaloClusterPtrVector clu, 
		  const std::vector<reco::TrackRef> tr,
		  const reco::Vertex  &  convVtx,
		  ConversionAlgorithm=undefined);
      
      Conversion( const reco::CaloClusterPtrVector clu, 
		  const std::vector<edm::RefToBase<reco::Track> > tr,
		  const reco::Vertex  &  convVtx,
		  ConversionAlgorithm=undefined);
      
           
      
      /// destructor
      virtual ~Conversion();
      /// returns a clone of the candidate
      Conversion * clone() const;
      /// Pointer to CaloCluster (foe Egamma Conversions it points to  a SuperCluster)
      reco::CaloClusterPtrVector caloCluster() const ;
      /// vector of track to base references 
      std::vector<edm::RefToBase<reco::Track> > tracks() const ; 
      /// returns  the reco conversion vertex
      const reco::Vertex & conversionVertex() const  { return theConversionVertex_ ; }
      /// Bool flagging objects having track size >0
      bool isConverted() const;
      /// Number of tracks= 0,1,2
      unsigned int nTracks() const {return  tracks().size(); }
      /// set the value  of the TMVA output
      void setMVAout(const float& mva) { theMVAout_=mva;}
      /// get the value  of the TMVA output
      double MVAout() const { return theMVAout_;}
      /// if nTracks=2 returns the pair invariant mass. Original tracks are used here
      double pairInvariantMass() const;
      /// Delta cot(Theta) where Theta is the angle in the (y,z) plane between the two tracks. Original tracks are used
      double pairCotThetaSeparation() const;
      /// Conversion tracks momentum from the tracks inner momentum
      GlobalVector  pairMomentum() const;
      /// Conversion track pair 4-momentum from the tracks refitted with vertex constraint
      math::XYZTLorentzVectorD   refittedPair4Momentum() const;
      /// Conversion tracks momentum from the tracks refitted with vertex constraint
      GlobalVector  refittedPairMomentum() const;
      /// Super Cluster energy divided by track pair momentum if Standard seeding method. If a pointer to two (or more clusters)
      /// is stored in the conversion, this method returns the energy sum of clusters divided by the track pair momentum
      /// Track innermost momentum is used here
      double EoverP() const;
      /// Super Cluster energy divided by track pair momentum if Standard seeing method. If a pointer to two (or more clusters)
      /// is stored in the conversion, this method returns the energy sum of clusters divided by the track pair momentum
      ///  Track momentum refitted with vertex constraint is used
      double EoverPrefittedTracks() const;
      /// z coordinate of the photon origin calculated from the track-pair direction and the position of the vertex
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
      /// Vector containing the position of the innermost hit of each track
      const std::vector<math::XYZPoint>& tracksInnerPosition() const {return theTrackInnerPosition_;}
      /// Vector of track momentum measured at the outermost hit
      const std::vector<math::XYZVector>& tracksPout() const {return theTrackPout_;}
      /// Vector of track momentum measured at the innermost hit
      const std::vector<math::XYZVector>& tracksPin() const  {return theTrackPin_;}
      
      /// Conversion Track algorithm/provenance
      void setConversionAlgorithm(const ConversionAlgorithm a, bool set=true) { if (set) algorithm_=a; else algorithm_=undefined;}
      ConversionAlgorithm algo() const ;
      std::string algoName() const;
      static std::string algoName(ConversionAlgorithm );
      static ConversionAlgorithm  algoByName(const std::string &name);      

      bool quality(ConversionQuality q) const { return  (qualityMask_ & (1<<q))>>q; }
      void setQuality(ConversionQuality q, bool b);


      
    private:
      
      /// vector pointer to a/multiple seed CaloCluster(s)
      reco::CaloClusterPtrVector caloCluster_;
      ///  vector of Track references
      std::vector<reco::TrackRef>  tracks_;
      /// vector Track RefToBase
      mutable std::vector<edm::RefToBase<reco::Track> >  trackToBaseRefs_;
      /// position at the ECAl surface of the track extrapolation
      std::vector<math::XYZPoint>  thePositionAtEcal_;
      /// Fitted Kalman conversion vertex
      reco::Vertex theConversionVertex_;
      /// Clusters mathing the tracks (these are not the seeds)
      std::vector<reco::CaloClusterPtr> theMatchingBCs_;
      /// Distance of min approach of the two tracks
      float theMinDistOfApproach_;
      /// P_in of tracks
      std::vector<math::XYZPoint> theTrackInnerPosition_;    
      /// P_in of tracks
      std::vector<math::XYZVector> theTrackPin_;    
      /// P_out of tracks
      std::vector<math::XYZVector> theTrackPout_;    
      /// TMVA output
      float theMVAout_;
      /// conversion algorithm/provenance
      uint8_t algorithm_;
      uint16_t qualityMask_;


  };

    inline Conversion::ConversionAlgorithm Conversion::algo() const {
      return (ConversionAlgorithm) algorithm_;
    }
    
    
    inline std::string Conversion::algoName() const{
            
      switch(algorithm_)
	{
	case undefined: return "undefined";
	case ecalSeeded: return "ecalSeeded";
	case trackerOnly: return "trackerOnly";
	case mixed: return "mixed";

	}
      return "undefined";
    }

    inline std::string Conversion::algoName(ConversionAlgorithm a){
      if(int(a) < int(algoSize) && int(a)>0) return algoNames[int(a)];
      return "undefined";
    }

    inline void Conversion::setQuality(ConversionQuality q, bool b){
      if (b)//regular OR if setting value to true
        qualityMask_ |= (1<<q) ;
      else // doing "half-XOR" if unsetting value
        qualityMask_ &= (~(1<<q));

    }
  
}

#endif
