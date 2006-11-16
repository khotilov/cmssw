#ifndef L1Scales_L1CaloEtScale_h
#define L1Scales_L1CaloEtScale_h
// -*- C++ -*-
//
// Package:     L1Scales
// Class  :     L1CaloEtScale
// 
/**\class L1CaloEtScale L1CaloEtScale.h L1Trigger/L1Scales/interface/L1CaloEtScale.h

 Description: Class to handle conversion between Et scales in L1 hardware

 Usage:
    <usage>

*/
//
// Author:      Jim Brooke
// Created:     Wed Sep 27 17:18:27 CEST 2006
// $Id: 
//

#include <boost/cstdint.hpp>
#include <vector>
#include <ostream>

class L1CaloEtScale {

 public:

  /// linear scale maximum
  static uint16_t linScaleMax;
  
  /// rank scale maximum
  static uint16_t rankScaleMax;

  /// default constructor (out = in)
  L1CaloEtScale();

  /// constructor takes physically meaningful quantities
  L1CaloEtScale(const double linearLsbInGeV, const std::vector<double> thresholdsInGeV);

  // destructor
  ~L1CaloEtScale();

  /// get LSB of linear input scale
  double linearLsb() const { return m_linearLsb; }

  /// convert from linear Et scale to rank scale
  uint16_t rank(const uint16_t linear) const;

  /// convert from rank to physically meaningful quantity
  double et(const uint16_t rank) const;

  void print(std::ostream& s) const;

 private:

  /// LSB of linear scale in GeV
  double m_linearLsb;

  /// thresholds associated with rank scale in GeV
  std::vector<double> m_thresholds;

};

#endif
