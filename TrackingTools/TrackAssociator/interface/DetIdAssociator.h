#ifndef TrackingTools_TrackAssociator_DetIdAssociator_h
#define TrackingTools_TrackAssociator_DetIdAssociator_h 1

// -*- C++ -*-
//
// Package:    TrackingTools/TrackAssociator
// Class:      DetIdAssociator
// 
/**\

 Description: Abstract base class for 3D point -> std::set<DetId>

 Implementation:
     A look up map of active detector elements in eta-phi space is 
     built to speed up access to the detector element geometry as well 
     as associated hits. The map is uniformly binned in eta and phi 
     dimensions. It is expected that the map is used to find a set of
     DetIds close to a given point, but since all methods are virtual 
     implementation may vary for various subdetectors.
**/
//
// Original Author:  Dmytro Kovalskyi
//         Created:  Fri Apr 21 10:59:41 PDT 2006
// $Id: DetIdAssociator.h,v 1.15 2009/09/06 16:34:11 dmytro Exp $
//
//

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "TrackingTools/TrackAssociator/interface/FiducialVolume.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"
#include <set>
#include <vector>

class DetIdAssociator{
 public:
   enum PropagationTarget { Barrel, ForwardEndcap, BackwardEndcap };
   struct MapRange {
      float dThetaPlus;
      float dThetaMinus;
      float dPhiPlus;
      float dPhiMinus;
   };
   typedef std::vector<GlobalPoint>::const_iterator const_iterator;
	
   DetIdAssociator();
   DetIdAssociator(const int nPhi, const int nEta, const double etaBinSize);
   
   virtual ~DetIdAssociator();
   
   /// Preselect DetIds close to a point on the inner surface of the detector. 
   /// "iN" is a number of the adjacent bins of the map to retrieve 
   virtual std::set<DetId> getDetIdsCloseToAPoint(const GlobalPoint&,
						  const int iN = 0) const;
   virtual std::set<DetId> getDetIdsCloseToAPoint(const GlobalPoint& direction,
						  const unsigned int iNEtaPlus,
						  const unsigned int iNEtaMinus,
						  const unsigned int iNPhiPlus,
						  const unsigned int iNPhiMinus) const;
   virtual std::set<DetId> getDetIdsCloseToAPoint(const GlobalPoint& direction,
						  const MapRange& mapRange) const;
   /// Preselect DetIds close to a point on the inner surface of the detector. 
   /// "d" defines the allowed range in theta-phi space:
   /// - theta is in [point.theta()-d, point.theta()+d]
   /// - phi is in [point.phi()-d, point.phi()+d]
   virtual std::set<DetId> getDetIdsCloseToAPoint(const GlobalPoint& point,
						  const double d = 0) const;
   /// - theta is in [point.theta()-dThetaMinus, point.theta()+dThetaPlus]
   /// - phi is in [point.phi()-dPhiMinus, point.phi()+dPhiPlus]
   virtual std::set<DetId> getDetIdsCloseToAPoint(const GlobalPoint& point,
						  const double dThetaPlus,
						  const double dThetaMinus,
						  const double dPhiPlus,
						  const double dPhiMinus) const;
   /// Find DetIds that satisfy given requirements
   /// - inside eta-phi cone of radius dR
   virtual std::set<DetId> getDetIdsInACone(const std::set<DetId>&,
					    const std::vector<GlobalPoint>& trajectory,
					    const double dR) const;
   /// - DetIds crossed by the track
   ///   tolerance is the radius of the trajectory used for matching
   ///   -1 is default and represent the case with no uncertainty
   ///   on the trajectory direction. It's the fastest option
   /// - DetIds crossed by the track, ordered according to the order
   ///   that they were crossed by the track flying outside the detector
   virtual std::vector<DetId> getCrossedDetIds(const std::set<DetId>&,
					       const std::vector<GlobalPoint>& trajectory) const;
   virtual std::vector<DetId> getCrossedDetIds(const std::set<DetId>&,
					       const std::vector<SteppingHelixStateInfo>& trajectory,
					       const double toleranceInSigmas = -1) const;
   /// look-up map eta index
   virtual int iEta (const GlobalPoint&) const;
   /// look-up map phi index
   virtual int iPhi (const GlobalPoint&) const;
   /// set a specific track propagator to be used
   virtual void setPropagator(Propagator* ptr){	ivProp_ = ptr; };
   /// number of bins of the look-up map in phi dimension
   int nPhiBins() const { return nPhi_;}
   /// number of bins of the look-up map in eta dimension
   int nEtaBins() const { return nEta_;}
   /// look-up map bin size in eta dimension
   double etaBinSize() const { return etaBinSize_;};
   /// make the look-up map
   virtual void buildMap();
   /// get active detector volume
   const FiducialVolume& volume() const;

   virtual void setGeometry(const DetIdAssociatorRecord&) = 0;
   virtual const GeomDet* getGeomDet(const DetId&) const = 0;

   virtual void setConditions(const DetIdAssociatorRecord&) {};
   
 protected:
   virtual void check_setup() const;
   
   virtual void dumpMapContent( int, int ) const;
   virtual void dumpMapContent( int, int, int, int ) const;
   
   virtual GlobalPoint getPosition(const DetId&) const = 0;
   virtual std::set<DetId> getASetOfValidDetIds() const = 0;
   virtual std::pair<const_iterator, const_iterator> getDetIdPoints(const DetId&) const = 0;
   
   virtual bool insideElement(const GlobalPoint&, const DetId&) const = 0;
   virtual bool crossedElement(const GlobalPoint&, 
			       const GlobalPoint&, 
			       const DetId&,
			       const double toleranceInSigmas = -1,
			       const SteppingHelixStateInfo* = 0 ) const { return false; }
   virtual bool nearElement(const GlobalPoint& point, 
			    const DetId& id, 
			    const double distance) const;
   
   // map parameters
   const int nPhi_;
   const int nEta_;
   std::set<DetId> **theMap_;
   bool theMapIsValid_;
   const double etaBinSize_;
   double maxEta_;
   double minTheta_;
   
   Propagator *ivProp_;
   
   // Detector fiducial volume 
   // approximated as a closed cylinder with non-zero width.
   // Parameters are extracted from the active detector elements.
   FiducialVolume volume_;
};
#endif
