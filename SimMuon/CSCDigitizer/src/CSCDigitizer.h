#ifndef MU_END_DIGITIZER
#define MU_END_DIGITIZER

// This is CSCDigitizer.h

/** \class CSCDigitizer
 *  Digitizer class for endcap muon CSCs.
 *
 *  \author Rick Wilkinson
 *
 */

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"


class CSCDriftSim;
class CSCWireHitSim;
class CSCStripHitSim;
class CSCWireElectronicsSim;
class CSCStripElectronicsSim;
class CSCLayer;
class CSCNeutronFactory;

class CSCDigitizer 
{
public:
  CSCDigitizer();
  ~CSCDigitizer();

  /**  digitize
   */
  void doAction(MixCollection<PSimHit> & simHits,
                CSCWireDigiCollection & wireDigis,
                CSCStripDigiCollection & stripDigis,
                CSCComparatorDigiCollection & comparators);

  /// sets geometry
  void setGeometry(const CSCGeometry * geom) {theCSCGeometry = geom;}

  /// sets the magnetic field
  void setMagneticField(const MagneticField * field);

  /// finds the layer in the geometry associated with this det ID
  const CSCLayer * findLayer(int detId) const;

private:
  CSCDriftSim            * theDriftSim;
  CSCWireHitSim          * theWireHitSim;
  CSCStripHitSim         * theStripHitSim;
  CSCWireElectronicsSim  * theWireElectronicsSim;
  CSCStripElectronicsSim * theStripElectronicsSim;
  CSCNeutronFactory      * theNeutronFactory;
  const CSCGeometry      * theCSCGeometry;
  bool doNeutrons;
};

#endif

