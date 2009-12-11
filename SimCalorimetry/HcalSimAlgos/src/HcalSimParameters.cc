#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameters.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CondFormats/HcalObjects/interface/HcalGain.h"
#include "CondFormats/HcalObjects/interface/HcalGainWidth.h"
#include "CLHEP/Random/RandGaussQ.h"
using namespace std;

HcalSimParameters::HcalSimParameters(double simHitToPhotoelectrons, double photoelectronsToAnalog,
                 double samplingFactor, double timePhase,
                 int readoutFrameSize, int binOfMaximum,
                 bool doPhotostatistics, bool syncPhase,
                 int firstRing, const std::vector<double> & samplingFactors)
: CaloSimParameters(simHitToPhotoelectrons, photoelectronsToAnalog, samplingFactor, timePhase,
                    readoutFrameSize, binOfMaximum, doPhotostatistics, syncPhase),
  theDbService(0),
  theFirstRing(firstRing),
  theSamplingFactors(samplingFactors),
  thePixels(0),
  doTimeSmear_(true)
{
  defaultTimeSmearing();
}

HcalSimParameters::HcalSimParameters(const edm::ParameterSet & p)
:  CaloSimParameters(p),
   theDbService(0),
   theFirstRing( p.getParameter<int>("firstRing") ),
   theSamplingFactors( p.getParameter<std::vector<double> >("samplingFactors") ),
   thePixels(0),
   doTimeSmear_( p.getParameter<bool>("timeSmearing"))
{
  if(p.exists("pixels"))
  {
    thePixels = p.getParameter<int>("pixels");
  }
  defaultTimeSmearing();
}

double HcalSimParameters::simHitToPhotoelectrons(const DetId & detId) const 
{
  // the gain is in units of GeV/fC.  We want a constant with pe/dGeV
  // pe/dGeV = (GeV/dGeV) / (GeV/fC) / (fC/pe) 
  double result = CaloSimParameters::simHitToPhotoelectrons(detId);
  if(HcalGenericDetId(detId).genericSubdet() != HcalGenericDetId::HcalGenForward
     || HcalGenericDetId(detId).genericSubdet() != HcalGenericDetId::HcalGenZDC)
    { 
      result = samplingFactor(detId) / fCtoGeV(detId) / photoelectronsToAnalog(detId);
    }
  return result;
}


double HcalSimParameters::fCtoGeV(const DetId & detId) const
{
  assert(theDbService != 0);
  HcalGenericDetId hcalGenDetId(detId);
  const HcalGain* gains = theDbService->getGain(hcalGenDetId);
  const HcalGainWidth* gwidths = theDbService->getGainWidth(hcalGenDetId);
  if (!gains || !gwidths )
  {
    edm::LogError("HcalAmplifier") << "Could not fetch HCAL conditions for channel " << hcalGenDetId;
  }
  // only one gain will be recorded per channel, so just use capID 0 for now
  double result = gains->getValue(0);
//  if(doNoise_)
///  {
//    result += CLHEP::RandGaussQ::shoot(0.,  gwidths->getValue(0));
//  }
  return result;
}

double HcalSimParameters::samplingFactor(const DetId & detId) const
{
  HcalDetId hcalDetId(detId);
  return theSamplingFactors.at(hcalDetId.ietaAbs()-theFirstRing);
}

//static const double GeV2fC = 1.0/0.145;
static const double GeV2fC = 1.0/0.4;

void HcalSimParameters::defaultTimeSmearing() {
  // GeV->ampl (fC), time (ns)
  theSmearSettings.push_back(std::pair<double,double>(  4.00*GeV2fC, 4.050));
  theSmearSettings.push_back(std::pair<double,double>( 20.00*GeV2fC, 3.300));
  theSmearSettings.push_back(std::pair<double,double>( 25.00*GeV2fC, 2.925));
  theSmearSettings.push_back(std::pair<double,double>( 30.00*GeV2fC, 2.714));
  theSmearSettings.push_back(std::pair<double,double>( 37.00*GeV2fC, 2.496));
  theSmearSettings.push_back(std::pair<double,double>( 44.50*GeV2fC, 2.278));
  theSmearSettings.push_back(std::pair<double,double>( 56.00*GeV2fC, 2.138));
  theSmearSettings.push_back(std::pair<double,double>( 63.50*GeV2fC, 2.022));
  theSmearSettings.push_back(std::pair<double,double>( 81.00*GeV2fC, 1.788));
  theSmearSettings.push_back(std::pair<double,double>( 88.50*GeV2fC, 1.695));
  theSmearSettings.push_back(std::pair<double,double>(114.50*GeV2fC, 1.716));
  theSmearSettings.push_back(std::pair<double,double>(175.50*GeV2fC, 1.070));
  theSmearSettings.push_back(std::pair<double,double>(350.00*GeV2fC, 1.564));
  theSmearSettings.push_back(std::pair<double,double>(99999.00*GeV2fC, 1.564));
}

double HcalSimParameters::timeSmearRMS(double ampl) const {
  HcalTimeSmearSettings::size_type i;
  double smearsigma=0;

  for (i=0; i<theSmearSettings.size(); i++)
    if (theSmearSettings[i].first > ampl)
      break;

  // Smearing occurs only within the envelope definitions.
  if (i!=0 && (i < theSmearSettings.size())) {
    double energy1 = theSmearSettings[i-1].first;
    double sigma1  = theSmearSettings[i-1].second;
    double energy2 = theSmearSettings[i].first;
    double sigma2  = theSmearSettings[i].second;
    
    if (energy2 != energy1)
      smearsigma = sigma1 + ((sigma2-sigma1)*(ampl-energy1)/(energy2-energy1));
    else
      smearsigma = (sigma2+sigma1)/2.;    
  }

  return smearsigma;

}
