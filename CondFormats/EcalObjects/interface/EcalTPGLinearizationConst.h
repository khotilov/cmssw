#ifndef EcalTPGLinearizationConst_h
#define EcalTPGLinearizationConst_h

#include "CondFormats/EcalObjects/interface/EcalCondObjectContainer.h"

struct EcalTPGLinearizationConstant
{
  struct Item 
  {
    uint32_t mult_x12 ;
    uint32_t mult_x6 ;
    uint32_t mult_x1 ;
    uint32_t shift_x12 ;
    uint32_t shift_x6 ;
    uint32_t shift_x1 ;
  };
};

typedef EcalCondObjectContainer<EcalTPGLinearizationConstant> EcalTPGLinearizationConstMap;
typedef EcalCondObjectContainer<EcalTPGLinearizationConstant>::const_iterator EcalTPGLinearizationConstMapIterator;
typedef EcalTPGLinearizationConstMap EcalTPGLinearizationConst;

//class EcalTPGLinearizationConst 
//{
// public:
//  EcalTPGLinearizationConst() ;
//  ~EcalTPGLinearizationConst() ;
//
//  struct Item 
//  {
//    uint32_t mult_x12 ;
//    uint32_t mult_x6 ;
//    uint32_t mult_x1 ;
//    uint32_t shift_x12 ;
//    uint32_t shift_x6 ;
//    uint32_t shift_x1 ;
//  };
//
//  const std::map<uint32_t, Item> & getMap() const { return map_; }
//  void  setValue(const uint32_t & id, const Item & value) ;
//
// private:
//  std::map<uint32_t, Item> map_ ;
//
//};
//
//typedef std::map<uint32_t, EcalTPGLinearizationConst::Item>                 EcalTPGLinearizationConstMap;
//typedef std::map<uint32_t, EcalTPGLinearizationConst::Item>::const_iterator EcalTPGLinearizationConstMapIterator;

#endif
