#ifndef MU_END_STRIP_ELECTRONICS_SIM_H
#define MU_END_STRIP_ELECTRONICS_SIM_H

/** \class CSCStripElectronicsSim
 * Model the readout electronics chain for EMU CSC strips
 *
 * \author Rick Wilkinson
 *
 */

#include "SimMuon/CSCDigitizer/src/CSCBaseElectronicsSim.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h"
#include "SimMuon/CSCDigitizer/src/CSCStripAmpResponse.h"

class CSCDetectorHit;
class CSCComparatorDigi;
class CSCCrosstalkGenerator;
class CSCScaNoiseGenerator;
class CSCStripConditions;
#include <vector>
#include <string>

class CSCStripElectronicsSim : public CSCBaseElectronicsSim
{
public:
  /// configurable parameters
  explicit CSCStripElectronicsSim(const edm::ParameterSet & p);

  virtual ~CSCStripElectronicsSim();

  void fillDigis(CSCStripDigiCollection & digis,
                 CSCComparatorDigiCollection & comparators);

  void setStripConditions(CSCStripConditions * cond) {theStripConditions = cond;}

private:
  /// initialization for each layer
  void initParameters();

  virtual int readoutElement(int strip) const;

  float calculateAmpResponse(float t) const;
  CSCStripAmpResponse theAmpResponse;

  void runComparator(std::vector<CSCComparatorDigi> & result);

  /// calculates the comparator reading, including saturation and offsets
  float comparatorReading(const CSCAnalogSignal & signal, float time) const;

  CSCAnalogSignal makeNoiseSignal(int element);

  virtual float signalDelay(int element, float pos) const;
  
  // tells which strips to read out around the input strip
  void getReadoutRange(int inputStrip, 
                       int & minStrip, int & maxStrip);

  /// finds the key strips from these comparators
  std::list<int>
  getKeyStrips(const std::vector<CSCComparatorDigi> & comparators) const;

  /// get ths strips that have detector hits
  std::list<int>
  getKeyStripsFromMC() const;
  /// finds what strips to read.  Will either take 5 strips around
  /// the keystrip, or the whole CFEB, based on doSuppression_
  std::list<int>
  channelsToRead(const std::list<int> & keyStrips, int window) const;

  void addCrosstalk();
  void addCrosstalk(const CSCAnalogSignal & signal,
                    int thisStrip, int otherStrip);


  void selfTest() const;

  void createDigi(int istrip, float startTime, std::vector<CSCStripDigi> & result);

  // saturation of the 12-bit ADC.  Max reading is 4095
  void doSaturation(CSCStripDigi & digi);

  // useful constants
  float theComparatorThreshold;      // in fC
  float theComparatorNoise;
  float theComparatorRMSOffset;
  // note that we don't implement the effect of the x3.5 amplifier
  float theComparatorSaturation;
  // all of these times are in nanoseconds
  float theComparatorWait;
  float theComparatorDeadTime;
  float theDaqDeadTime;
  // save the calculation of time-of-flight+drift+shaping
  float theTimingOffset;

  int nScaBins_;
  bool doSuppression_;
  bool doCrosstalk_;
  CSCStripConditions * theStripConditions;
  CSCCrosstalkGenerator * theCrosstalkGenerator;

  int theComparatorClockJump;
  // the length of each SCA time bin, in ns.  50 by default
  float sca_time_bin_size;
  // the SCA bin which holds the peak signal.  4, by default.
  // that's really the 5th, since We start counting at 0
  int   sca_peak_bin;
  // which time bin the trigger crossing goes in
  double theComparatorTimeBinOffset;
  // to center comparator signals
  double theComparatorTimeOffset;

};

#endif

