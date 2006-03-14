#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimMuon/CSCDigitizer/src/CSCWireElectronicsSim.h"
#include "SimMuon/CSCDigitizer/src/CSCDetectorHit.h"
#include "SimMuon/CSCDigitizer/src/CSCAnalogSignal.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigi.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamberSpecs.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussQ.h"
#include <iostream>


CSCWireElectronicsSim::CSCWireElectronicsSim(const edm::ParameterSet & p) 
 : CSCBaseElectronicsSim()
{
  theSignalStartTime = p.getParameter<double>("wireSignalStartTime");
  theSignalStopTime = p.getParameter<double>("wireSignalStopTime");
  theSamplingTime = p.getParameter<double>("wireSamplingTime");
  theTimingCalibrationError = p.getParameter<double>("wireTimingError");
  init();
}


CSCWireElectronicsSim::CSCWireElectronicsSim()
 : CSCBaseElectronicsSim()
{
  theSignalStartTime = -250.;
  theSignalStopTime = 500.;
  theSamplingTime = 2.0;
  theTimingCalibrationError = 0.;
  init();
}


void CSCWireElectronicsSim::init() {
  theFraction = 0.5;
  theShapingTime = 30;
  theAmpGainVariance = 0.;
  thePeakTimeVariance = 0.;
  theTailShaping = RADICAL;

  theBunchTimingOffsets.resize(11);
  theBunchTimingOffsets[1] = 33.8;
  theBunchTimingOffsets[2] = 34.8;
  theBunchTimingOffsets[3] = 42.0;
  theBunchTimingOffsets[4] = 42.0;
  theBunchTimingOffsets[5] = 43.6;
  theBunchTimingOffsets[6] = 42.0;
  theBunchTimingOffsets[7] = 42.6;
  theBunchTimingOffsets[8] = 41.5;
  theBunchTimingOffsets[9] = 43.6;
  theBunchTimingOffsets[10] = 41.5;
  theNumberOfSamples = (int)((theSignalStopTime-theSignalStartTime)/theSamplingTime);
  fillAmpResponse();
}


void CSCWireElectronicsSim::initParameters() {
  theLayerGeometry = theLayer->geometry();
  nElements = theLayerGeometry->numberOfWireGroups();
  theWireNoise = theSpecs->wireNoise(theShapingTime)
                * e_SI * pow(10.0,15);
  theWireThreshold = theWireNoise * 8;
}


int CSCWireElectronicsSim::readoutElement(int element) const {
  return theLayerGeometry->wireGroup(element);
}

void CSCWireElectronicsSim::fillDigis(CSCWireDigiCollection & digis) {

  if(theSignalMap.empty()) {
    return;
  }

  const float fiveElectronCharge = 5.*e_SI*1.e15*theSpecs->gasGain()*theSpecs->fractionQS();

  // Loop over analog signals, run the fractional discriminator on each one,
  // and save the DIGI in the layer.
  for(MESignalMap::iterator mapI = theSignalMap.begin(); 
      mapI != theSignalMap.end(); ++mapI) {
    int wireGroup            = (*mapI).first;
    CSCAnalogSignal signal = (*mapI).second;
    LogDebug("CSCWireElectronicsSim") << "CSCWireElectronicsSim: dump of wire signal follows... " <<  signal;
    // the way we handle noise in this chamber is by randomly varying
    // the threshold
    float threshold = theWireThreshold + RandGaussQ::shoot() * theWireNoise;
    for(int ibin = 0; ibin < signal.getSize(); ++ibin)
      if(signal.getBinValue(ibin) > threshold) {
        // jackpot.  Now define this signal as everything up until
        // the signal goes below zero.
        int lastbin = signal.getSize();
        int i;
        for(i = ibin; i < signal.getSize(); ++i) {
          if(signal.getBinValue(i) < 0.) {
            lastbin = i;
            break;
          }
        }

      int bin_of_5e_arrival = 0;
      float qMax = 0.0;
      // in this loop, find the max charge and the 'fifth' electron arrival
      for ( i = ibin; i < lastbin; ++i)
      {
        float next_charge = signal.getBinValue(i);
        if( bin_of_5e_arrival == 0 
            && next_charge >= fiveElectronCharge)
        {
          bin_of_5e_arrival = i;
        }
        if(next_charge > qMax) {
          qMax = next_charge;
        }
      }
     
      int bin_firing_FD = 0;
      for ( i = ibin; i < lastbin; ++i)
      {
        if( bin_firing_FD == 0 && signal.getBinValue(i) >= qMax * theFraction )
        {
           bin_firing_FD = i;
        }
      } 
      float tofOffset = timeOfFlightCalibration(wireGroup);
      int chamberType = theSpecs->chamberType();
      // fill in the rest of the information, for everything that survives the fraction discrim.
      //int t_5thElectron = 
      //   static_cast<int>(theSignalStartTime + theSamplingTime*bin_of_5e_arrival );

      // Shouldn't one measure from theTimeOffset of the CSCAnalogSignal?
      // Well, yes, but unfortunately it seems CSCAnalogSignal::superimpose
      // fails to reset theTimeOffset properly. In any case, if everything
      // is self-consistent, the overall theTimeOffset should BE
      // theSignalStartTime. There is big trouble if any of the
      // overlapping MEAS's happen to have a timeOffset earlier than
      // theSignalStartTime (which is why it defaults to -10bx = -250ns).
      // That could only happen in the case of signals from pile-up events
      // arising way earlier than 10bx before the signal event crossing
      // and so only if pile-up were simulated over an enormous range of
      // bx earlier than the signal bx.
      // (Comments from Tim, Aug-2005)

      float fdTime = theSignalStartTime + theSamplingTime*bin_firing_FD
                   + RandGaussQ::shoot() * theTimingCalibrationError;
      //int fractionalDiscriminatorTime = 
      //   static_cast<int>(fdTime  - tofOffset);
      int beamCrossingTag = 
         static_cast<int>( (fdTime - tofOffset - theBunchTimingOffsets[chamberType])
                            / theBunchSpacing );
      // no pedestal
      //int adcCounts = 
      //   static_cast<int>( qMax/theSpecs->chargePerCount() );
      CSCWireDigi newDigi(wireGroup, beamCrossingTag);
      LogDebug("CSCWireElectronicsSim") << newDigi;
      digis.insertDigi(layerId(), newDigi);

      // we code the channels so strips start at 0, wire groups at 100
      //@@ I believe above comment is wrong... (Tim, Aug-2005)
      //@@ ... We count from 1 unless there's a _strong_ reason to count from 0
      addLinks(channelIndex(wireGroup));
      // skip over all the time bins used for this digi
      ibin = lastbin;
    } // loop over time bins in signal
  } // loop over wire signals   
}


float CSCWireElectronicsSim::calculateAmpResponse(float t) const {
  static const float fC_by_ns = 1000000;
  static const float resistor = 20000;
  static const float amplifier_pole               = 1/7.5;
  static const float fastest_chamber_exp_risetime = 10.;
  static const float p0=amplifier_pole;
  static const float p1=1/fastest_chamber_exp_risetime;

  static const float dp = p0 - p1;

  // ENABLE DISC:

  static const float norm = -12 * resistor * p1 * pow(p0/dp, 4) / fC_by_ns;

  float enable_disc_volts = norm*(  exp(-p0*t) *(1          +
						 t*dp       +
						 pow(t*dp,2)/2 +
						 pow(t*dp,3)/6  )
				    - exp(-p1*t) );
  static const float collectionFraction = 0.12;
  static const float igain = 1./0.005; // volts per fC
  return enable_disc_volts * igain * collectionFraction;
}                                                                               


float CSCWireElectronicsSim::signalDelay(int element, float pos) const {
  // readout is on right edge of chamber, signal speed is c
  // zero calibrated to chamber center
  // pos is assumed to be in wire coordinates, not local
  float distance = -1. * pos; // in cm
  float speed = c_light / cm;
  float delay = distance / speed;
  return delay;
}

float CSCWireElectronicsSim::timeOfFlightCalibration(int wireGroup) const {
  // calibration is done for groups of 8 wire groups, facetiously
  // called wireGroupGroups
  int middleWireGroup = wireGroup - wireGroup%8 + 4;
  int numberOfWireGroups = theLayerGeometry->numberOfWireGroups();
  if(middleWireGroup > numberOfWireGroups) 
     middleWireGroup = numberOfWireGroups;

//  LocalPoint centerOfGroupGroup = theLayerGeometry->centerOfWireGroup(middleWireGroup);
//  float averageDist = theLayer->surface().toGlobal(centerOfGroupGroup).mag();
  GlobalPoint centerOfGroupGroup = theLayer->centerOfWireGroup(middleWireGroup);
  float averageDist = centerOfGroupGroup.mag();


  float averageTOF  = averageDist * cm / c_light; // Units of c_light: mm/ns

  LogDebug("CSCWireElectronicsSim") << "MEWES TofCalib:  wg = " << wireGroup << 
       " mid wg = " << middleWireGroup << 
       " av dist = " << averageDist << 
      " av tof = " << averageTOF;
  
  return averageTOF;
}
 
