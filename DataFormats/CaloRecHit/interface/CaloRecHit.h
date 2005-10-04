#ifndef DATAFORMATS_CALORECHIT_CALORECHIT_H
#define DATAFORMATS_CALORECHIT_CALORECHIT_H 1

#include "DataFormats/DetId/interface/DetId.h"
#include <ostream>


/** \class CaloRecHit
 * 
 * $Date: 2005/09/26 14:10:28 $
 * $Revision: 1.3 $
 *\author J. Mans - Minnesota
 */
class CaloRecHit {
public:
  CaloRecHit(); // for persistence
  explicit CaloRecHit(const DetId& id, float energy, float time);
  virtual ~CaloRecHit();
  float energy() const { return energy_; }
  float time() const { return time_; }
  const DetId& detid() const { return id_; }
private:
  DetId id_;
  float energy_;
  float time_;
};

std::ostream& operator<<(std::ostream& s, const CaloRecHit& hit);
  
#endif
