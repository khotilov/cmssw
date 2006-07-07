#ifndef L1CALOREGIONDETID_H
#define L1CALOREGIONDETID_H

#include "DataFormats/DetId/interface/DetId.h"

/** \class L1CaloRegionDetId
 *  Cell identifier class for L1 Calo Trigger Regions (4x4 trigger tower sums)
 *
 *  $Date: 2006/07/06 12:38:11 $
 *  $Revision: 1.8 $
 *  \author Jim Brooke 
*/

/**
 * Stores eta value in bits 4-0, phi value in bits 9-5
 *
 *
 */

class L1CaloRegionDetId : public DetId {

 public:

  static const unsigned N_PHI;
  static const unsigned N_ETA;

  /// create null id
  L1CaloRegionDetId();
  
  /// create id from raw data (0=invalid code?)
  L1CaloRegionDetId(uint32_t rawid);

  /// create id from global eta, phi indices (eta=0-21, phi=0-17)
  L1CaloRegionDetId(unsigned ieta, unsigned iphi);

  /// create id from RCT crate, RCT card, RCT region (within card)
  /// set isForward to true to create forward regions (ignoring card argument)
  /// or to false to create central regions (including card argument)
  L1CaloRegionDetId(bool isForward, unsigned icrate, unsigned icard, unsigned irgn);

  /// create id from GCT card and input number
  /// NB - isForward has no effect; dummy argument to differentiate from global eta/phi indices!
  L1CaloRegionDetId(bool isForward, unsigned icard, unsigned irgn);

  /// global eta index (0-21)
  unsigned ieta() const { return id_&0x1f; }

  /// global phi index (0-17)
  unsigned iphi() const { return (id_>>5)&0x1f; }

  /// return GCT source card number
  unsigned gctCard() const;

  /// return GCT region index (within source card)
  unsigned gctRegion() const;

  /// return central or forward type
  bool isForward() const { return (ieta()<4 || ieta()>17); }

  /// return RCT crate number (0-17)
  unsigned rctCrate() const;

  /// return RCT card number (0-6)
  unsigned rctCard() const;  /// TODO - check

  /// return RCT region index (0-1 for barrel, 0-7 for HF)
  unsigned rctRegion() const; /// TODO - check

  /// return local RCT eta index (0-10)
  unsigned rctEta() const { return (ieta()<11 ? 10-ieta() : ieta()-11); }  /// TODO - check this is correct

  /// return local RCT phi index (0-1)
  unsigned rctPhi() const { return (iphi()%2); }  /// TODO - check this is correct

};

#endif
