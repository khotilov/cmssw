#ifndef IOMC_FlatEvtVtxGenerator_H
#define IOMC_FlatEvtVtxGenerator_H

/**
 * Generate event vertices according to a Flat distribution. 
 * Attention: All values are assumed to be cm!
 *
 * $Id: FlatEvtVtxGenerator.h,v 1.2 2006/11/07 19:38:31 wdd Exp $
 */

#include "IOMC/EventVertexGenerators/interface/BaseEvtVtxGenerator.h"

namespace CLHEP {
   class RandFlat;
}

class FlatEvtVtxGenerator : public BaseEvtVtxGenerator 
{
public:
  FlatEvtVtxGenerator(const edm::ParameterSet & p);
  virtual ~FlatEvtVtxGenerator();

  /// return a new event vertex
  //virtual CLHEP::Hep3Vector* newVertex();
  virtual HepMC::FourVector* newVertex() ;

  /// set min in X in cm
  void minX(double m=0.0);
  /// set min in Y in cm
  void minY(double m=0.0);
  /// set min in Z in cm
  void minZ(double m=0.0);

  /// set max in X in cm
  void maxX(double m=0);
  /// set max in Y in cm
  void maxY(double m=0);
  /// set max in Z in cm
  void maxZ(double m=0);
  
private:
  /** Copy constructor */
  FlatEvtVtxGenerator(const FlatEvtVtxGenerator &p);
  /** Copy assignment operator */
  FlatEvtVtxGenerator&  operator = (const FlatEvtVtxGenerator & rhs );
private:
  double fMinX, fMinY, fMinZ;
  double fMaxX, fMaxY, fMaxZ;
  CLHEP::RandFlat*  fRandom ; 
};

#endif
