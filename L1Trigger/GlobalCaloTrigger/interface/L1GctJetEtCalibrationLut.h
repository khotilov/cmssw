#ifndef L1GCTJETETCALIBRATIONLUT_H_
#define L1GCTJETETCALIBRATIONLUT_H_

#include "L1Trigger/GlobalCaloTrigger/src/L1GctLut.h"

class L1GctJetEtCalibrationFunction;

/*!
 * \author Robert Frazier & Greg Heath
 * \date May 2006
 */

/*! \class L1GctJetEtCalibrationLut
 * \brief Jet Et calibration LUT
 * 
 * Input is 10 bit Et and 4 bit eta
 * Outputs are 6 bit rank (for jet sorting) and 10 bit Et (for Ht calculation)
 * 
 * Modified March 2007 to remove the actual calculation to a separate class
 *
 */

class L1GctJetEtCalibrationLut : public L1GctLut<15,16>
{
 public:
  static const int NAddress;
  static const int NData;
  static const unsigned JET_ENERGY_BITWIDTH;
  static L1GctJetEtCalibrationLut* setupLut(const L1GctJetEtCalibrationFunction* lutfn);
  virtual ~L1GctJetEtCalibrationLut();

  void setFunction(const L1GctJetEtCalibrationFunction* lutfn);
  const L1GctJetEtCalibrationFunction* getFunction() const { return m_lutFunction; }

  /// Overload << operator
  friend std::ostream& operator << (std::ostream& os, const L1GctJetEtCalibrationLut& lut);
  
 protected:
  
  L1GctJetEtCalibrationLut();

  virtual uint16_t value (const uint16_t lutAddress) const;

 private:

  const L1GctJetEtCalibrationFunction* m_lutFunction;

};

std::ostream& operator << (std::ostream& os, const L1GctJetEtCalibrationLut& lut);

#endif /*L1GCTJETETCALIBRATIONLUT_H_*/
