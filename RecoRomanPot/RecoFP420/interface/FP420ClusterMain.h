#ifndef FP420ClusterMain_h
#define FP420ClusterMain_h
   
#include <string>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimRomanPot/SimFP420/interface/DigiCollectionFP420.h"
#include "RecoRomanPot/RecoFP420/interface/ClusterCollectionFP420.h"
#include "RecoRomanPot/RecoFP420/interface/ClusterNoiseFP420.h"
#include "RecoRomanPot/RecoFP420/interface/ClusterFP420.h"
#include <iostream>
#include <vector>
using namespace std;

class ClusterNoiseFP420;
class ClusterProducerFP420;

class FP420ClusterMain 
{
 public:
  

    FP420ClusterMain(const edm::ParameterSet& conf, int sn, int pn);
  //  FP420ClusterMain();

  ~FP420ClusterMain();

  /// Runs the algorithm
    //void run(const DigiCollectionFP420* input,
  void run(const DigiCollectionFP420 &input,
	   ClusterCollectionFP420 &soutput,
	   const std::vector<ClusterNoiseFP420>& noise 
	   );
  //	   ,unsigned int);

 private:


  ClusterProducerFP420 *threeThreshold_;
  std::string clusterMode_;

  //std::vector<HDigiFP420> collector;
  edm::ParameterSet conf_;


  bool validClusterizer_;
  double ElectronPerADC_;
  double ENC_;
  double BadElectrodeProbability_;
  bool UseNoiseBadElectrodeFlagFromDB_;

  
  double ChannelThreshold;
  double SeedThreshold;
  double ClusterThreshold;
  int MaxVoidsInCluster;	

  double ldriftX;
  double ldriftY;
  double ldrift;
  double pitchX;          // pitchX
  double pitchY;          // pitchY
  double pitch;          // pitch automatic
  float moduleThicknessX; // plate thicknessX 
  float moduleThicknessY; // plate thicknessY 
  float moduleThickness; // plate thickness 
  int numStripsX, numStripsXW;    // number of strips in the moduleX
  int numStripsY, numStripsYW;    // number of strips in the moduleY
  int numStrips;    // number of strips in the module

  float Thick300;

 // Number of Stations:
   int sn0;
 // Number of planes:
   int pn0;

   int verbosity;

 //float sigma1_;
 //float sigma2_;

};

#endif
