#ifndef HCALSIMPLERECALGO_H
#define HCALSIMPLERECALGO_H 1

#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalDigi/interface/HFDataFrame.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"

/** \class HcalSimpleRecAlgo

   This class reconstructs RecHits from Digis for HBHE, HF, and HO by addition
   of selected time samples, pedestal subtraction, and gain application. The
   time of the hit is reconstructed using a weighted peak bin calculation
   supplemented by precise time lookup table. A consumer of this class also
   has the option of correcting the reconstructed time for energy-dependent
   time slew associated with the QIE.
    
   $Date: 2006/02/07 23:29:10 $
   $Revision: 1.4 $
   \author J. Mans - Minnesota
*/
class HcalSimpleRecAlgo {
public:
  HcalSimpleRecAlgo(int firstSample, int samplesToAdd, bool correctForTimeslew);
  HBHERecHit reconstruct(const HBHEDataFrame& digi, const HcalCoder& coder, const HcalCalibrations& calibs) const;
  HFRecHit reconstruct(const HFDataFrame& digi, const HcalCoder& coder, const HcalCalibrations& calibs) const;
  HORecHit reconstruct(const HODataFrame& digi, const HcalCoder& coder, const HcalCalibrations& calibs) const;
private:
  int firstSample_, samplesToAdd_;
  bool correctForTimeslew_;
};

#endif
