#ifndef TrackingTools_MuonUpdatorAtVertex_H
#define TrackingTools_MuonUpdatorAtVertex_H

 /**  \class MuonUpdatorAtVertex
  *
  *   Extrapolate a muon trajectory to 
  *   a given vertex and 
  *   apply a vertex constraint
  *
  *   $Date: 2006/07/21 04:02:43 $
  *   $Revision: 1.4 $
  *
  *   \author   N. Neumeister            Purdue University
  *
  */



#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/CommonDetAlgo/interface/GlobalError.h"
#include "RecoMuon/TrackingTools/interface/MuonVertexMeasurement.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

class TrajectoryStateOnSurface;
class TransverseImpactPointExtrapolator;
class KFUpdator;
class MeasurementEstimator;

namespace edm {class ParameterSet; class EventSetup;}

class MuonUpdatorAtVertex {

  public:
 
    /// constructor from parameter set
    MuonUpdatorAtVertex(const edm::ParameterSet&);

    /// default constructor
    MuonUpdatorAtVertex();

    /// constructor from propagator
    MuonUpdatorAtVertex(const Propagator*);

    /// destructor
    virtual ~MuonUpdatorAtVertex();

    /// initialize propagator EventSetup
    void setES(const edm::EventSetup&);

    /// set Propagator directly
    void setPropagator(Propagator*);

    /// set Propagator from 2 propagators, tk & gen
    void setPropagator(Propagator*,Propagator*,const MagneticField*);

    void setPropagator(const Propagator&,const Propagator&, const MagneticField*);


    /// return vertex measurement
    MuonVertexMeasurement update(const TrajectoryStateOnSurface&) const;

    /// only return the state on outer tracker surface
    TrajectoryStateOnSurface stateAtTracker(const TrajectoryStateOnSurface&) const;

    void setVertex(const GlobalPoint, const GlobalError);

    
  private:
 
    GlobalPoint theVertexPos;
    GlobalError theVertexErr;

    Propagator* thePropagator;
    TransverseImpactPointExtrapolator* theExtrapolator;
    KFUpdator* theUpdator;
    MeasurementEstimator* theEstimator;
    std::string theOutPropagatorName;
    std::string theInPropagatorName;

};

#endif

