#include "SimGeneral/NoiseGenerators/interface/GaussianTailNoiseGenerator.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"

#include <math.h>

#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_result.h>

GaussianTailNoiseGenerator::GaussianTailNoiseGenerator(CLHEP::HepRandomEngine& eng ) :
  gaussDistribution_(eng),
  poissonDistribution_(eng),
  flatDistribution_(eng)
{
}

void GaussianTailNoiseGenerator::generate(int NumberOfchannels, 
					  float threshold, 
					  float noiseRMS, 
					  std::map<int,float, std::less<int> >& theMap ) {

   // Gaussian tail probability
  gsl_sf_result result;
  int status = gsl_sf_erf_Q_e(threshold, &result);

  if (status != 0) std::cerr<<"GaussianTailNoiseGenerator::could not compute gaussian tail probability for the threshold chosen"<<std::endl;

  float probabilityLeft = result.val;  
  float meanNumberOfNoisyChannels = probabilityLeft * NumberOfchannels;
  int numberOfNoisyChannels = poissonDistribution_.fire(meanNumberOfNoisyChannels);

  float lowLimit = threshold * noiseRMS;
  for (int i = 0; i < numberOfNoisyChannels; i++) {

    // Find a random channel number    
    int theChannelNumber = (int)flatDistribution_.fire(NumberOfchannels);

    // Find random noise value
    double noise = generate_gaussian_tail(lowLimit, noiseRMS);
              
    // Fill in map
    theMap[theChannelNumber] = noise;
    
  }
}

void GaussianTailNoiseGenerator::generate(int NumberOfchannels, 
					  float threshold, 
					  float noiseRMS, 
					  std::vector<std::pair<int,float> > &theVector ) {
  // Compute number of channels with noise above threshold
  // Gaussian tail probability
  gsl_sf_result result;
  int status = gsl_sf_erf_Q_e(threshold, &result);
  
  if (status != 0) std::cerr<<"GaussianTailNoiseGenerator::could not compute gaussian tail probability for the threshold chosen"<<std::endl;
  
  double probabilityLeft = result.val;  
  double meanNumberOfNoisyChannels = probabilityLeft * NumberOfchannels;
  int numberOfNoisyChannels = poissonDistribution_.fire(meanNumberOfNoisyChannels);

  theVector.reserve(numberOfNoisyChannels);
  float lowLimit = threshold * noiseRMS;
  for (int i = 0; i < numberOfNoisyChannels; i++) {

    // Find a random channel number    
    int theChannelNumber = (int)flatDistribution_.fire(NumberOfchannels);
    
    // Find random noise value
    double noise = generate_gaussian_tail(lowLimit, noiseRMS);
              
    // Fill in the vector
    theVector.push_back(std::pair<int, float>(theChannelNumber, noise));
  }
}

void GaussianTailNoiseGenerator::generateRaw(int NumberOfchannels, 
					     float noiseRMS, 
					     std::vector<std::pair<int,float> > &theVector ) {
  theVector.reserve(NumberOfchannels);
  for (int i = 0; i < NumberOfchannels; i++) {
    // Find random noise value
    float noise = gaussDistribution_.fire(0.,noiseRMS);
    // Fill in the vector
    theVector.push_back(std::pair<int, float>(i,noise));
  }
}

double
GaussianTailNoiseGenerator::generate_gaussian_tail(const double a, const double sigma){
  /* Returns a gaussian random variable larger than a
   * This implementation does one-sided upper-tailed deviates.
   */
  
  double s = a/sigma;
  
  if (s < 1){
    /*
      For small s, use a direct rejection method. The limit s < 1
      can be adjusted to optimise the overall efficiency 
    */
    double x;
    
    do{
      x = gaussDistribution_.fire(0.,1.0);
    }
    while (x < s);
    return x * sigma;
    
  }else{
    
    /* Use the "supertail" deviates from the last two steps
     * of Marsaglia's rectangle-wedge-tail method, as described
     * in Knuth, v2, 3rd ed, pp 123-128.  (See also exercise 11, p139,
     * and the solution, p586.)
     */
    
    double u, v, x;
    
    do{
      u = flatDistribution_.fire();
      do{
	v = flatDistribution_.fire();
      }while (v == 0.0);
      x = sqrt(s * s - 2 * log(v));
    }
    while (x * u > s);
    return x * sigma;
  }
}
