#ifndef MU_END_WIRE_ELECTRONICS_SIM_H
#define MU_END_WIRE_ELECTRONICS_SIM_H

/** \class CSCWireElectronicsSim
 * Model the readout electronics chain for EMU CSC wires
 *
 * \author Rick Wilkinson
 *
 * Last mod: <BR>
 * 30-Jun-00 ptc Minor printout mods. <BR>
 * 06-Jul-00 ptc Update beamCrossingTag to better match reality. <BR>
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimMuon/CSCDigitizer/src/CSCBaseElectronicsSim.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"

// declarations
class CSCLayer;
class CSCDetectorHit;
class CSCWireDigi;
class CSCAnalogSignal;


class CSCWireElectronicsSim : public CSCBaseElectronicsSim
{
public:
  /// configurable parameters
  CSCWireElectronicsSim(const edm::ParameterSet &p);

  /// hardcoed default parameters
  CSCWireElectronicsSim();

  void setFraction(float newFraction)  {theFraction = newFraction;};

  void fillDigis(CSCWireDigiCollection & digis);

private:
  // helper functions
  void init();
  /// initialization for each layer
  virtual void initParameters();

  // will return wire group, given wire.
  virtual int readoutElement(int element) const;

  float calculateAmpResponse(float t) const;
 
  virtual float signalDelay(int element, float pos) const;
  virtual float timeOfFlightCalibration(int wireGroup) const;

  /// we code strip indices from 1-80, and wire indices start at 100
  virtual int channelIndex(int channel) const {return channel+100;}

  // member data
  // the fractional discriminator returns the time when the signal
  // reaches this fraction of its maximum
  float theFraction;
  float theWireNoise;
  float theWireThreshold;
  float theTimingCalibrationError; // in ns
};

#endif
