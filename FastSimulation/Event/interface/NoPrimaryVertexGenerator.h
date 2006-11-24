#ifndef FastSimulation_Event_NoPrimaryVertexGenerator_H
#define FastSimulation_Event_NoPrimaryVertexGenerator_H

// CLHEP Headers
#include "FastSimulation/Event/interface/PrimaryVertexGenerator.h"

/** A class that generates a primary vertex with a Gaussian profile, in cm*/ 

class NoPrimaryVertexGenerator : public PrimaryVertexGenerator {

public:
  /// Default constructor
  NoPrimaryVertexGenerator() : PrimaryVertexGenerator() {;}

  /// Destructor
  ~NoPrimaryVertexGenerator() {;}
  
  /// Generation process (to be implemented)
  virtual void generate() {;}

 private:

};

#endif // NoPrimaryVertexGenerator_H
