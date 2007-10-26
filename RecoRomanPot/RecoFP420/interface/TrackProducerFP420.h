#ifndef TrackProducerFP420_h
#define TrackProducerFP420_h

#include "RecoRomanPot/RecoFP420/interface/TrackFP420.h"
#include "RecoRomanPot/RecoFP420/interface/ClusterFP420.h"
#include "RecoRomanPot/RecoFP420/interface/ClusterCollectionFP420.h"

#include <vector>
#include <algorithm>
#include <cmath>


class TrackProducerFP420 {
public:

  typedef std::vector<ClusterFP420>::const_iterator           ClusterFP420Iter;

  //TrackProducerFP420(int, int, double, double, double, double, double, double, double, double, bool, bool, double, double, float, float);
    TrackProducerFP420(int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, bool, bool, bool, bool, double, double, float, float, double);

    std::vector<TrackFP420> trackFinderMaxAmplitude(ClusterCollectionFP420 input);
    std::vector<TrackFP420> trackFinderMaxAmplitude2(ClusterCollectionFP420 input);

    // std::vector<TrackFP420> trackFinderVar1(ClusterCollectionFP420 input);
    //std::vector<TrackFP420> trackFinderVar2(ClusterCollectionFP420 input);

    std::vector<TrackFP420> trackFinderSophisticated(ClusterCollectionFP420 input);

    std::vector<TrackFP420> trackFinder3D(ClusterCollectionFP420 input);

private:
  ClusterCollectionFP420 soutput;

  std::vector<TrackFP420> rhits; 

 // Number of Stations:
 int sn0;
 // Number of planes:
 int pn0;
 // Number of planes:
 int zn0;

 // shift of planes:
	bool UseHalfPitchShiftInX;
	bool UseHalfPitchShiftInY;
	bool UseHalfPitchShiftInXW;
	bool UseHalfPitchShiftInYW;

	//double zUnit; 
	double z420; 
	double zD2; 
	double zD3; 
	double pitchX;
	double pitchY;
	double pitchXW;
	double pitchYW;
        double ZGapLDet;
	//double ZBoundDet;
	double ZSiStep;
	double ZSiPlane;
	double ZSiDetL;
	double ZSiDetR;

	double dXX;
	double dYY;
	float chiCutX;
	float chiCutY;

	double zinibeg;

};



#endif
