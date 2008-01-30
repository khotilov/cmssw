#ifndef CSCSegment_CSCSegAlgoShowering_h
#define CSCSegment_CSCSegAlgoShowering_h

/**
 * class CSCSegAlgoShowering
 *
 *  \author: D. Fortin - UC Riverside
 *
 * Handle case where too many hits are reconstructed in the chamber, even after preclustering
 * for normal segment reconstruction to properly handle these.
 * In this case, determine the average local (x,y) for each layer and find the hit closest to that localpoint
 * for that given layer.  From these hits, reconstruct a segment.
 * The idea is to give the reconstruction (seeding) a valid starting point for the kalman filter.
 */

#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>
#include <DataFormats/CSCRecHit/interface/CSCSegment.h>
#include <Geometry/CSCGeometry/interface/CSCChamber.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <vector>

class CSCSegAlgoShowering {

 public:

  typedef std::vector<const CSCRecHit2D*> ChamberHitContainer;

  /// Constructor
  explicit CSCSegAlgoShowering(const edm::ParameterSet& ps);

  /// Destructor
  virtual ~CSCSegAlgoShowering();

  CSCSegment showerSeg( const CSCChamber* aChamber, ChamberHitContainer rechits );


 private:

  /// Utility functions 	
  bool isHitNearSegment(const CSCRecHit2D* h) const;
  bool addHit(const CSCRecHit2D* hit, int layer);
  void updateParameters(void);
  bool hasHitOnLayer(int layer) const;
  void compareProtoSegment(const CSCRecHit2D* h, int layer);
  AlgebraicSymMatrix calculateError(void) const;
  void pruneFromResidual();

  // Member variables
  const std::string myName; 
  const CSCChamber* theChamber;

  ChamberHitContainer protoSegment;
  float       protoSlope_u;
  float       protoSlope_v;
  LocalPoint  protoIntercept;		
  double      protoChi2;
  LocalVector protoDirection;

  // input from .cfi file
  bool   debug;
  int    minHitsPerSegment;
  double dRPhiFineMax;
  double dPhiFineMax;
  float  tanPhiMax;
  float  tanThetaMax;
  float  chi2Max;
  float  maxRatioResidual;

};
#endif

