#ifndef BasicVertexState_H
#define BasicVertexState_H

#include "TrackingTools/TrajectoryState/interface/ProxyBase.h"
#include "Geometry/Surface/interface/ReferenceCounted.h"
#include "TrackingTools/TrajectoryState/interface/CopyUsingClone.h"

#include "Geometry/Vector/interface/GlobalPoint.h"
#include "Geometry/CommonDetAlgo/interface/GlobalError.h"
#include "Geometry/CommonDetAlgo/interface/GlobalWeight.h"
#include "Geometry/CommonDetAlgo/interface/AlgebraicObjects.h"
//#include "CommonReco/CommonVertex/interface/RefCountedVertexSeed.h"

#include <vector>

class VertexState;

/** Class containing a measurement of a vertex.
 */

class BasicVertexState  : private ReferenceCounted {

public:

  typedef ProxyBase< BasicVertexState, CopyUsingClone<BasicVertexState> > Proxy;
  typedef ReferenceCountingPointer<BasicVertexState>    		  RCPtr;

private:
  friend class Proxy;
  friend class RCPtr;

public:

  virtual ~BasicVertexState() {}

  virtual BasicVertexState* clone() const = 0;

  /** Access methods
   */
  virtual GlobalPoint position() const = 0;
  virtual GlobalError error() const = 0;
  virtual GlobalWeight weight() const = 0;
  virtual AlgebraicVector weightTimesPosition() const = 0;
  virtual double weightInMixture() const = 0;
  virtual std::vector<VertexState> components() const;

  /** conversion to VertexSeed
   */
//   virtual RefCountedVertexSeed seedWithoutTracks() const = 0;
};

#endif
