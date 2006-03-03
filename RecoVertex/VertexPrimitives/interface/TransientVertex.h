#ifndef _Vtx_TransientVertex_H_
#define _Vtx_TransientVertex_H_

// #include "CommonReco/CommonVertex/interface/RecoVertex.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/CommonDetAlgo/interface/GlobalError.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexPrimitives/interface/TTtoTTmap.h"

#include <vector>
#include <map>

class CachingVertex;

/** \class TransientVertex
 *  The concrete RecVertex implementation.
 */

class TransientVertex {

public:

  typedef std::map<reco::TransientTrack, float> TransientTrackToFloatMap;

  /** Empty constructor, produces invalid vertex
   */
  TransientVertex();

  /** Constructor defining the RecVertex by its 3D position
   *  and position uncertainty, its associated tracks
   *  and its chi-squared.
   *  The number of degrees of freedom is equal to
   *  2*nb of tracks - 3.
   */
  TransientVertex(const GlobalPoint & pos, const GlobalError & posError,
	    const std::vector<reco::TransientTrack> & tracks, float chi2);

  /** Constructor defining the RecVertex by its 3D position
   *  and position uncertainty, its associated tracks, its chi-squared
   *  and its number of degrees of freedom.
   *  The ndf can be a float.
   */
  TransientVertex(const GlobalPoint & pos, const GlobalError & posError,
	    const std::vector<reco::TransientTrack> & tracks, float chi2, float ndf);

  /** Constructor defining the RecVertex by the prior,
   *  the vertex 3D position and uncertainty, the associated tracks
   *  and the chi-squared. Since the prior brings information on
   *  3 coordinates, the number of degrees of freedom is equal to
   *  2*nb of tracks.
   */
  TransientVertex(const GlobalPoint & priorPos, const GlobalError & priorErr,
  	    const GlobalPoint & pos, const GlobalError & posError,
	    const std::vector<reco::TransientTrack> & tracks, float chi2);

  /** Constructor defining the RecVertex by the prior,
   *  the vertex 3D position and uncertainty, the associated tracks,
   *  the chi-squared and the number of degrees of freedom.
   *  The ndf can be a float.
   */
  TransientVertex(const GlobalPoint & priorPos, const GlobalError & priorErr,
  	    const GlobalPoint & pos, const GlobalError & posError,
	    const std::vector<reco::TransientTrack> & tracks, float chi2, float ndf);

  /** Constructor defining the RecVertex by its 3D position 
   *  and position uncertainty, its associated tracks 
   *  and its chi-squared. 
   *  The number of degrees of freedom is equal to 
   *  2*nb of tracks - 3.
   */
  TransientVertex(const VertexState & state, 
		    const std::vector<reco::TransientTrack> & tracks, float chi2);

  /** Constructor defining the RecVertex by its 3D position
   *  and position uncertainty, its associated tracks, its chi-squared
   *  and its number of degrees of freedom.
   *  The ndf can be a float.
   */
  TransientVertex(const VertexState & state, 
		    const std::vector<reco::TransientTrack> & tracks, float chi2, float ndf);

  /** Constructor defining the RecVertex by the prior,
   *  the vertex 3D position and uncertainty, the associated tracks
   *  and the chi-squared. Since the prior brings information on
   *  3 coordinates, the number of degrees of freedom is equal to
   *  2*nb of tracks.
   */
  TransientVertex(const VertexState & prior,
		    const VertexState & state,
		    const std::vector<reco::TransientTrack> & tracks, float chi2);

  /** Constructor defining the RecVertex by the prior,
   *  the vertex 3D position and uncertainty, the associated tracks,
   *  the chi-squared and the number of degrees of freedom.
   *  The ndf can be a float.
   */
  TransientVertex(const VertexState & prior,
		    const VertexState & state,
		    const std::vector<reco::TransientTrack> & tracks, float chi2, float ndf);


//   /** Constructor defining the RecVertex by its 3D position
//    *  and position uncertainty, its associated tracks, its chi-squared
//    *  and its number of degrees of freedom, and the track weights. 
//    *  The ndf can be a float.
//    */
//   TransientVertex(const VertexState & state, 
// 		    const std::vector<reco::TransientTrack> & tracks, float chi2, float ndf, 
// 		    const reco::TransientTrackToFloatMap & weightMap);


  /** Access methods
   */
  VertexState vertexState() const { return theVertexState; }
  GlobalPoint position() const { return theVertexState.position(); }
  GlobalError positionError() const { return theVertexState.error(); }
  GlobalPoint priorPosition() const { return thePriorVertexState.position(); }
  GlobalError priorError() const { return thePriorVertexState.error(); }
  bool hasPrior() const { return withPrior; }

//   /** Implements method of abstract Vertex.
//    *  Returns track pointer container by value
//    */
//   Vertex::TrackPtrContainer tracks() const;

  float totalChiSquared() const { return theChi2; }
  float normalisedChiSquared() const {
    return totalChiSquared() / degreesOfFreedom();
  }
  float degreesOfFreedom() const { return theNDF; }

  /** Returns true if vertex is valid.
   *  An invalid RecVertex is created e.g. when vertex fitting fails.
   *  WARNING: obsolete. All RecVertex objects are valid.
   */
  bool isValid() const {
    return vertexValid;
  }

  /** Access to the original tracks used to make the vertex.
   *  Returns track container by value.
   */
  std::vector<reco::TransientTrack> originalTracks() const {
    return theOriginalTracks;
  }


  /**
   * Returns true if at for at least one of the original tracks,
   * the refitted track is available
   */
  bool hasRefittedTracks() const { return withRefittedTracks; }


  /** Access to the refitted tracks used to make the vertex.
   *  Returns track container by value.
   */
  std::vector<reco::TransientTrack> refittedTracks() const {
    return theRefittedTracks;
  }

  /**
   * Returns true if the track-weights are available.
   */
  bool hasTrackWeight() const { return theWeightMapIsAvailable; }


  /**
   *   Returns the weight with which a track has been used in the fit.
   *   If the track is not present in the list, it is obviously not used, and
   *   this method returns a weight of 0.
   *   If this information has not been provided at construction, a weight of
   *   1.0 is assumed for all tracks used and present in the originalTracks() std::vector.
   */
  float trackWeight(const reco::TransientTrack track) const;

  TransientTrackToFloatMap weightMap() const { return theWeightMap; }

  void weightMap(const TransientTrackToFloatMap & theMap);

  /**
   * Returns true if the Track-to-track covariance matrices have been calculated.
   */
  bool tkToTkCovarianceIsAvailable() const { return theCovMapAvailable; }

  /**
   *   Returns the Track-to-track covariance matrix for two specified tracks.
   *   In case these do not exist, or one of the tracks does not belong to the
   *   vertex, an exception is thrown.
   */
  AlgebraicMatrix tkToTkCovariance(const reco::TransientTrack& t1, 
  				const reco::TransientTrack& t2) const;
  void tkToTkCovariance(const TTtoTTmap covMap);

private:


  mutable VertexState thePriorVertexState;
  mutable VertexState theVertexState;

  std::vector<reco::TransientTrack> theOriginalTracks;
  std::vector<reco::TransientTrack> theRefittedTracks;

  float theChi2;
  float theNDF;
  bool vertexValid;
  bool withPrior, theWeightMapIsAvailable, theCovMapAvailable;
  bool withRefittedTracks;
  TTtoTTmap theCovMap;
  TransientTrackToFloatMap theWeightMap;

};

#endif
