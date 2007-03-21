#ifndef DataFormats_SiStripDetId_TOBDetId_H
#define DataFormats_SiStripDetId_TOBDetId_H

#include <ostream>
#include <vector>
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

/** 
 *  Det identifier class for the TOB
 */
class TOBDetId;

std::ostream& operator<<(std::ostream& os,const TOBDetId& id);

class TOBDetId : public SiStripDetId {
 public:
  /** Constructor of a null id */
  TOBDetId();
  /** Constructor from a raw value */
  TOBDetId(uint32_t rawid);
  /**Construct from generic DetId */
  TOBDetId(const DetId& id); 
  
  TOBDetId(uint32_t layer,
	   uint32_t rod_fw_bw,
	   uint32_t rod,
	   uint32_t module,
	   uint32_t ster) : SiStripDetId(DetId::Tracker,StripSubdetector::TOB){
    id_ |= (layer& layerMask_) << layerStartBit_ |
      (rod_fw_bw& rod_fw_bwMask_) << rod_fw_bwStartBit_ |
      (rod& rodMask_) << rodStartBit_ |
      (module& moduleMask_) << moduleStartBit_ |
      (ster& sterMask_) << sterStartBit_ ;
  }
  
  
  /// layer id
  unsigned int layer() const{
    return int((id_>>layerStartBit_) & layerMask_);
  }
  
  /// rod id
  /**
   * vector[0] = 1 -> bw rod (TOB-)
   * vector[0] = 2 -> fw rod (TOB+)
   * vector[1] -> rod
   */
  std::vector<unsigned int> rod() const
    { std::vector<unsigned int> num;
    num.push_back(((id_>>rod_fw_bwStartBit_) & rod_fw_bwMask_));
    num.push_back(((id_>>rodStartBit_) & rodMask_));
    return num ;}
  
  /// detector id
  unsigned int module() const 
    { return ((id_>>moduleStartBit_)& moduleMask_) ;}
  
 private:
  /// two bits would be enough, but  we could use the number "0" as a wildcard
  static const unsigned int layerStartBit_=     14;
  static const unsigned int rod_fw_bwStartBit_= 12;
  static const unsigned int rodStartBit_=       5;
  static const unsigned int moduleStartBit_=    2;
  static const unsigned int sterStartBit_=      0;
  /// two bits would be enough, but  we could use the number "0" as a wildcard
  
  static const unsigned int layerMask_=       0x7;
  static const unsigned int rod_fw_bwMask_=   0x3;
  static const unsigned int rodMask_=         0x7F;
  static const unsigned int moduleMask_=      0x7;
  static const unsigned int sterMask_=        0x3;
};


#endif
