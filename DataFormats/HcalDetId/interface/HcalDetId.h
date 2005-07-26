#ifndef CALOID_HCALCELLID_H
#define CALOID_HCALCELLID_H

#include <ostream>
#include <boost/cstdint.hpp>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

namespace cms {

/** \class HcalDetId
 *  Cell identifier class for the HCAL subdetectors
 *
 *  $Date: 2005/07/20 00:10:52 $
 *  $Revision: 1.2 $
 *  \author J. Mans - Minnesota
 */
class HcalDetId : public DetId {
public:
  /** Create a null cellid*/
  HcalDetId();
  /** Create cellid from raw id (0=invalid tower id) */
  HcalDetId(uint32_t rawid);
  /** Constructor from subdetector, signed tower ieta,iphi,and depth */
  HcalDetId(HcalSubdetector subdet, int tower_ieta, int tower_iphi, int depth);
  /** Constructor from a generic cell id */
  HcalDetId(const DetId& id);
  /** Assignment from a generic cell id */
  HcalDetId& operator=(const DetId& id);

  /// get the subdetector
  HcalSubdetector subdet() const { return (HcalSubdetector)(subdetId()); }
  /// get the z-side of the cell (1/-1)
  int zside() const { return (id_&0x2000)?(1):(-1); }
  /// get the absolute value of the cell ieta
  int ietaAbs() const { return (id_>>7)&0x3f; }
  /// get the cell ieta
  int ieta() const { return zside()*ietaAbs(); }
  /// get the cell iphi
  int iphi() const { return id_&0x7F; }
  /// get the tower depth
  int depth() const { return (id_>>14)&0x7; }
  /// get the smallest crystal_ieta of the crystal in front of this tower (HB and HE tower 17 only)
  int crystal_ieta_low() const { return ((ieta()-zside())*5)+zside(); }
  /// get the largest crystal_ieta of the crystal in front of this tower (HB and HE tower 17 only)
  int crystal_ieta_high() const { return ((ieta()-zside())*5)+5*zside(); }
  /// get the smallest crystal_iphi of the crystal in front of this tower (HB and HE tower 17 only)
  int crystal_iphi_low() const { return ((iphi()-1)*5)+1; }
  /// get the largest crystal_iphi of the crystal in front of this tower (HB and HE tower 17 only)
  int crystal_iphi_high() const { return ((iphi()-1)*5)+5; }
  /// compact index for arrays, etc [assumes ieta/iphi/depth sufficient to fully id tower]
  int hashedIndex() const;

};

std::ostream& operator<<(std::ostream&,const HcalDetId& id);

}

#endif
