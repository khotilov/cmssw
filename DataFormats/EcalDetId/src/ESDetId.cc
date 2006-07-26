#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include <stdexcept>

ESDetId::ESDetId() : DetId() {
}
  
ESDetId::ESDetId(uint32_t rawid) : DetId(rawid) {
}
  
ESDetId::ESDetId(int strip, int ixs, int iys, int plane, int iz) : DetId(Ecal,EcalPreshower) {
  if ((strip<ISTRIP_MIN) || (strip > ISTRIP_MAX) ||
      (ixs<IX_MIN) || (ixs > IX_MAX) ||
      (iys<IY_MIN) || (iys > IY_MAX) ||
      (plane != 1 && plane != 2)) 
    throw(std::runtime_error("ESDetId:  Cannot create object.  Indexes out of bounds."));
  id_ |=
    (strip&0x3F) |
    ((ixs&0x3F)<<6) |
    ((iys&0x3F)<<12) |
    (((plane-1)&0x1)<<18) |
    ((iz>0)?(1<<19):(0));
}
  
ESDetId::ESDetId(const DetId& gen) {
  if (!gen.null() && ( gen.det()!=Ecal || gen.subdetId()!=EcalPreshower )) {
    throw new std::exception();
  }
  id_=gen.rawId();
}
  
ESDetId& ESDetId::operator=(const DetId& gen) {
  if (!gen.null() && ( gen.det()!=Ecal || gen.subdetId()!=EcalPreshower )) {
    throw new std::exception();
  }
  id_=gen.rawId();
  return *this;
}
  
int ESDetId::hashedIndex() const {
  // TODO: more efficient index!
  return id_&0xFFFFFF;
}
  
std::ostream& operator<<(std::ostream& s,const ESDetId& id) {
  return s << "(ES z=" << id.zside() << "  plane " << id.plane() << " " <<
    id.six() << ':' << id.siy() << " " << id.strip() << ')';
}
