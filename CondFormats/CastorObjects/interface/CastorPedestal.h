#ifndef CastorPedestal_h
#define CastorPedestal_h

/** 
\class CastorPedestal
\author Panos Katsas
POOL object to store Pedestal values 4xCapId
$Author: katsas
*/
#include <boost/cstdint.hpp>

class CastorPedestal {
 public:
  /// get value for all capId = 0..3
  const float* getValues () const {return &mValue0;}
  /// get value for capId = 0..3
  float getValue (int fCapId) const {return *(getValues () + fCapId);}

  // functions below are not supposed to be used by consumer applications

  CastorPedestal () : mId (0), mValue0 (0), mValue1 (0), mValue2 (0), mValue3 (0) {}
  
  CastorPedestal (unsigned long fId, float fCap0, float fCap1, float fCap2, float fCap3) :
    mId (fId),
    mValue0 (fCap0),
    mValue1 (fCap1),
    mValue2 (fCap2),
    mValue3 (fCap3) {}

    uint32_t rawId () const {return mId;}
 private:
  uint32_t mId;
  float mValue0;
  float mValue1;
  float mValue2;
  float mValue3;
};

#endif
