#ifndef CALORECHIT_H
#define CALORECHIT_H 1

#include "DataFormats/DetId/interface/DetId.h"
#include <ostream>

namespace cms {

  /** \class CaloRecHit
    
  $Date: 2005/07/26 15:51:28 $
  $Revision: 1.1 $
  \author J. Mans - Minnesota
  */
  class CaloRecHit {
  public:
    CaloRecHit(); // for persistence
    explicit CaloRecHit(const DetId& id, float energy, float time);
    virtual ~CaloRecHit();
    float energy() const { return energy_; }
    float time() const { return time_; }
    const DetId& id() const { return id_; }
  private:
    DetId id_;
    float energy_;
    float time_;
  };

  std::ostream& operator<<(std::ostream& s, const CaloRecHit& hit);
}
  
#endif
