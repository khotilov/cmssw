// -*- C++ -*-
//
// Package:     L1Scales
// Class  :     L1CaloEtScale
// 
// Implementation:
//     <Notes on implementation>
//
// Author:      
// Created:     Wed Sep 27 17:18:27 CEST 2006
// $Id: 

#include "L1Trigger/L1Scales/interface/L1CaloEtScale.h"

using std::vector;
using std::ostream;
using std::endl;

uint16_t L1CaloEtScale::linScaleMax = 0x3ff;
uint16_t L1CaloEtScale::rankScaleMax = 0x3f;

// default constructor (testing only!)
L1CaloEtScale::L1CaloEtScale() :
  m_linearLsb(1.0)
{

  for (unsigned i=0; i<rankScaleMax; i++) {
    m_thresholds[i] = m_linearLsb * i;
  }

}

// real constructor
L1CaloEtScale::L1CaloEtScale(const double linearLsbInGeV, const vector<double> thresholdsInGeV) :
  m_linearLsb(linearLsbInGeV),
  m_thresholds(thresholdsInGeV) {

  // protect against too many thresholds!
  //  while ( m_threshold.size() > (L1GctJetScale::maxRank+1) ) {
  //    m_thresholds.pop_back();
  //  }

}


L1CaloEtScale::~L1CaloEtScale() {

}

// convert from linear Et to rank
uint16_t L1CaloEtScale::rank(const uint16_t linear) const {

  uint16_t out = 0;

  for (unsigned i=0; i<m_thresholds.size() && i<(rankScaleMax+1); i++) {
    if ( ( (linear & linScaleMax) * m_linearLsb) > m_thresholds[i] ) { out = i; }
  }

  return out & rankScaleMax;

}

// convert from rank to Et/GeV
double L1CaloEtScale::et(const uint16_t rank) const {

  if (rank < m_thresholds.size()-1) {
    return (m_thresholds[rank+1]+m_thresholds[rank]) / 2;
  }
  else {
    return m_thresholds.back();
  }

}

void L1CaloEtScale::print(ostream& s) const {
  s << "L1CaloEtScale" << endl;
  s << "L1CaloEtScale : linear LSB = " << m_linearLsb << " GeV" << endl;
  for (unsigned i=0; i<m_thresholds.size(); i++) {
    s << "L1CaloEtScale : threshold " << i << " = " << m_thresholds[i] << " GeV" << endl;
  }
}


