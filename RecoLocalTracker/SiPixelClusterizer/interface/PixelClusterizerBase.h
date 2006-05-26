#ifndef RecoLocalTracker_SiPixelClusterizer_PixelClusterizerBase_H
#define RecoLocalTracker_SiPixelClusterizer_PixelClusterizerBase_H

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include <vector>

class PixelGeomDetUnit;

/**
 * Abstract interface for Pixel Clusterizers
 */
class PixelClusterizerBase {
public:
  typedef edm::DetSet<PixelDigi>::const_iterator    DigiIterator;

  // Virtual destructor, this is a base class.
  virtual ~PixelClusterizerBase() {}

  // Build clusters in a DetUnit. Both digi and cluster stored in a DetSet

  virtual void clusterizeDetUnit( const edm::DetSet<PixelDigi> & input,	
				  const PixelGeomDetUnit * pixDet,
				  const std::vector<short>& badChannels,
				  edm::DetSet<SiPixelCluster>& output) = 0;

  // TO DO: the way we pass bad channels is most likely have to change later.

};

#endif
