#ifndef MagGeometry_H
#define MagGeometry_H

/** \class MagGeometry
 *  Entry point to the geometry of magnetic volumes.
 *
 *  $Date: 2005/09/27 15:15:52 $
 *  $Revision: 1.3 $
 *  \author N. Amapane - INFN Torino
 */

#include "Geometry/Surface/interface/BoundPlane.h"
#include "MagneticField/Layers/src/MagBinFinders.h"
/* #include "Utilities/Notification/interface/Singleton.h" */
//#include "DetectorDescription/Base/interface/Singleton.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"

#include <vector>

class MagBLayer;
class MagESector;
class MagVolume;
class MagVolume6Faces;
template <class T> class PeriodicBinFinderInPhi;
namespace edm {class ParameterSet;}

class MagGeometry {
public:

  typedef Surface::GlobalVector   GlobalVector;
  typedef Surface::GlobalPoint    GlobalPoint;

  /// Constructor
  MagGeometry(const edm::ParameterSet& config, std::vector<MagBLayer *> ,
			     std::vector<MagESector *> ,
			     std::vector<MagVolume6Faces*> ,
			     std::vector<MagVolume6Faces*> );

  /// Destructor
  ~MagGeometry();

  /// Return field vector at the specified global point
  GlobalVector fieldInTesla(const GlobalPoint & gp) const;

  /// Find a volume
  MagVolume * findVolume(const GlobalPoint & gp) const;

  // FIXME: only for temporary tests, should be removed.
  const std::vector<MagVolume6Faces*> & barrelVolumes() {return theBVolumes;}
  const std::vector<MagVolume6Faces*> & endcapVolumes() {return theEVolumes;}

private:
//  friend class DDI::Singleton<MagGeometry>;

  friend class MagGeometryExerciser; // for debug purposes

  // Linear search (for debug purposes only)
  MagVolume * findVolume1(const GlobalPoint & gp, double tolerance=0.) const;

  // Use hierarchical structure for the barrel with non-optimized R search
  MagVolume * findVolume2(const GlobalPoint & gp, double tolerance=0.) const;

  bool inBarrel(const GlobalPoint& gp) const;

  mutable MagVolume * lastVolume; // Cache last volume found

  std::vector<MagBLayer *> theBLayers;
  std::vector<MagESector *> theESectors;

  // FIXME: only for temporary tests, should be removed.
  std::vector<MagVolume6Faces*> theBVolumes;
  std::vector<MagVolume6Faces*> theEVolumes;

  MagBinFinders::GeneralBinFinderInR<double>* theBarrelBinFinder;
  PeriodicBinFinderInPhi<float> * theEndcapBinFinder;

  double tolerance;
  bool cacheLastVolume;
  bool timerOn;
};
#endif

