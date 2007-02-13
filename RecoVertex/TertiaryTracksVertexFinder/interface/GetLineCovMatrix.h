#ifndef GETLINECOVMATRIX_H
#define GETLINECOVMATRIX_H

#include "Geometry/CommonDetAlgo/interface/GlobalError.h" 
#include "CLHEP/Matrix/Matrix.h"

// Class that calculates the Covariance Matrix of a Globalpoint
// The GlobalPoint is on a straight line, that is defined by two 
// GlobalPoints( plus their Covariance Matrices)


class GetLineCovMatrix {

 public:

  GetLineCovMatrix(GlobalPoint, GlobalPoint, GlobalError, GlobalError);

  ~GetLineCovMatrix() {}

  GlobalError GetMatrix(GlobalPoint); 
  
 private:  

  GlobalPoint PointOne;
  GlobalPoint PointTwo;
  HepMatrix CombinedErrorMatrix; // CombinedErrorMatrix of the two points that define the straight line   
  HepMatrix B;                   // derivatives of the linear equation

};

#endif
