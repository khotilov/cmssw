#ifndef HcalChannelStatus_h
#define HcalChannelStatus_h

/* 
\class HcalChannelStatus
\author Radek Ofierzynski
contains one channel status and corresponding DetId
*/

#include <boost/cstdint.hpp>


class HcalChannelStatus
{
 public:
  // contains the defined bits for easy access, see https://twiki.cern.ch/twiki/bin/view/CMS/HcalDataValidationWorkflow
  enum StatusBit {       
    HcalCellOff=0,      // 1=Hcal cell is off
    HcalCellL1Mask=1,   // 1=Hcal cell is masked/to be masked by L1 trigger
    HcalCellDead=5,     // 1=Hcal cell is dead (from DQM algo)
    HcalCellHot=6,      // 1=Hcal cell is hot (from DQM algo)
    HcalCellStabErr=7,  // 1=Hcal cell has stability error
    HcalCellTimErr=8    // 1=Hcal cell has timing error
  };

  HcalChannelStatus(): mId(0), mStatus(0) {}
  HcalChannelStatus(unsigned long fid, uint32_t status): mId(fid), mStatus(status) {}

  //  void setDetId(unsigned long fid) {mId = fid;}
  void setValue(uint32_t value) {mStatus = value;}

  // for the following, one can use unsigned int values or HcalChannelStatus::StatusBit values
  //   e.g. 5 or HcalChannelStatus::HcalCellDead
  void setBit(unsigned int bitnumber) 
  {
    uint32_t statadd = 0x1<<(bitnumber);
    mStatus = mStatus|statadd;
  }
  void unsetBit(unsigned int bitnumber) 
  {
    uint32_t statadd = 0x1<<(bitnumber);
    statadd = ~statadd;
    mStatus = mStatus&statadd;
  }
  
  bool isBitSet(unsigned int bitnumber) const
  {
    uint32_t statadd = 0x1<<(bitnumber);
    return (mStatus&statadd)?(true):(false);
  }
  
  uint32_t rawId() const {return mId;}
  
  uint32_t getValue() const {return mStatus;}
  
 private:
  uint32_t mId;
  uint32_t mStatus;

};
#endif
