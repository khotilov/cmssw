#ifndef HcalGain_h
#define HcalGain_h

/** 
\class HcalGain
\author Fedor Ratnikov (UMd)
POOL object to store Gain values 4xCapId
$Author: ratnikov
$Date: 2007/01/09 22:49:20 $
$Revision: 1.4 $
*/

#include <boost/cstdint.hpp>

class HcalGain {
 public:
  /// get value for all capId = 0..3
  const float* getValues () const {return &mValue0;}
  /// get value for capId = 0..3
  float getValue (int fCapId) const {return *(getValues () + fCapId);}

  // functions below are not supposed to be used by consumer applications

  HcalGain () : mId (0), mValue0 (0), mValue1 (0), mValue2 (0), mValue3 (0) {}
  
  HcalGain (unsigned long fId, float fCap0, float fCap1, float fCap2, float fCap3) :
    mId (fId),
    mValue0 (fCap0),
    mValue1 (fCap1),
    mValue2 (fCap2),
    mValue3 (fCap3) {}

    // because of an oracle digestion problem with uint32_t 
    // use unsigned long long
    //  uint32_t rawId () const {return mId;}
  unsigned long long rawId () const {return mId;}

 private:
  //  uint32_t mId;
  unsigned long long mId;
  float mValue0;
  float mValue1;
  float mValue2;
  float mValue3;
};

#endif
