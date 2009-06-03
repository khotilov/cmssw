#ifndef RPCDigitizer_RPCSimSimple_h
#define RPCDigitizer_RPCSimSimple_h

/** \class RPCSimSimple
 *   Class for the RPC strip response simulation based
 *   on a very simple model
 *
 *  \author Marcello Maggi -- INFN Bari
 */
#include "SimMuon/RPCDigitizer/src/RPCSim.h"
#include "SimMuon/RPCDigitizer/src/RPCSynchronizer.h"

class RPCGeometry;

namespace CLHEP {
  class HepRandomEngine;
  class RandFlat;
  class RandPoissonQ;
}

class RPCSimSimple : public RPCSim
{
 public:
  RPCSimSimple(const edm::ParameterSet& config);
  ~RPCSimSimple();

  void simulate(const RPCRoll* roll,
		const edm::PSimHitContainer& rpcHits);

  void simulateNoise(const RPCRoll*);

 private:
  void init(){};

  RPCSynchronizer* _rpcSync;
  int N_hits;
  int nbxing;
  double rate;
  double gate;

  //  CLHEP::HepRandomEngine* rndEngine;
  //  CLHEP::RandFlat* flatDistribution_;
  //  CLHEP::RandPoissonQ *poissonDistribution_;
};
#endif
