///////////////////////////////////////////////////////////////////////////////
// File: FP420TrackMain.cc
// Date: 12.2006
// Description: FP420TrackMain for FP420
// Modifications: 
///////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <iostream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoRomanPot/RecoFP420/interface/FP420TrackMain.h"
#include "SimRomanPot/SimFP420/interface/ClusterFP420.h"
#include "RecoRomanPot/RecoFP420/interface/TrackFP420.h"
#include "RecoRomanPot/RecoFP420/interface/TrackProducerFP420.h"

#include "CLHEP/Random/RandFlat.h"

using namespace std;

//#define mytrackdebug0

//FP420TrackMain::FP420TrackMain(){ 
FP420TrackMain::FP420TrackMain(const edm::ParameterSet& conf):conf_(conf)  { 
  
  edm::ParameterSet m_Anal = conf.getParameter<edm::ParameterSet>("FP420TrackMain");
  verbosity    = m_Anal.getParameter<int>("Verbosity");

  //trackMode_         = "TrackProducerMaxAmplitudeFP420";
  //trackMode_         = "TrackProducerVar1FP420";
  //trackMode_         = "TrackProducerVar2FP420";
  // trackMode_         = "TrackProducerSophisticatedFP420";
  trackMode_  =  m_Anal.getParameter<std::string>("TrackModeFP420");
  
  //  sn0_ = 4;// related to  number of station: sn0=4 mean 3 Stations
  // pn0_ = 9;// related to number of planes: pn0=5 mean 4 Planes
    sn0_ = m_Anal.getParameter<int>("NumberFP420Stations");
    pn0_ = m_Anal.getParameter<int>("NumberFP420SPlanes");

    z420_           = m_Anal.getParameter<double>("z420");
    zD2_            = m_Anal.getParameter<double>("zD2");
    zD3_            = m_Anal.getParameter<double>("zD3");
    //  dXX_ = 12.7+0.05;//(BoxYshft+dYGap) + (YSi - YSiDet)/2. = 12.7+0.05
    dXX_ = m_Anal.getParameter<double>("dXXFP420");//(BoxYshft+dYGap) + (YSi - YSiDet)/2. = 12.7
    dYY_ = m_Anal.getParameter<double>("dYYFP420");//  XSiDet/2. = 5.0
    chiCutX_ = m_Anal.getParameter<double>("chiCutX420");//  =3
    chiCutY_ = m_Anal.getParameter<double>("chiCutY420");//  =3
    
    if (verbosity > 0) {
      std::cout << "FP420TrackMain constructor::" << std::endl;
      std::cout << "sn0=" << sn0_ << " pn0=" << pn0_ << std::endl;
      std::cout << "trackMode = " << trackMode_ << std::endl;
      std::cout << "dXX=" << dXX_ << " dYY=" << dYY_ << std::endl;
      std::cout << "chiCutX=" << chiCutX_ << " chiCutY=" << chiCutY_ << std::endl;
    }
  ///////////////////////////////////////////////////////////////////
      // zD2_ = 1000.;  // dist between centers of 1st and 2nd stations
      // zD3_ = 8000.;  // dist between centers of 1st and 3rd stations
  
  UseHalfPitchShiftInX_= true;
  UseHalfPitchShiftInY_= true;

  pitchX_= 0.050;
  pitchY_= 0.050;// was 0.040

//
  double zBlade = 5.00;
  double gapBlade = 1.6;
  ZSiPlane_=2*(zBlade+gapBlade);
  
  double ZKapton = 0.1;
  ZSiStep_=ZSiPlane_+ZKapton;
  
  double ZBoundDet = 0.020;
  double ZSiElectr = 0.250;
  double ZCeramDet = 0.500;
//
  ZSiDetL_ = 0.250;
  ZSiDetR_ = 0.250;
//
  ZGapLDet_= zBlade/2-(ZSiDetL_+ZSiElectr+ZBoundDet+ZCeramDet/2);
//
    if (verbosity > 1) {
      std::cout << "FP420TrackMain constructor::" << std::endl;
      std::cout << " zD2=" << zD2_ << " zD3=" << zD3_ << std::endl;
      std::cout << " UseHalfPitchShiftInX=" << UseHalfPitchShiftInX_ << " UseHalfPitchShiftInY=" << UseHalfPitchShiftInY_ << std::endl;
      std::cout << " pitchX=" << pitchX_ << " pitchY=" << pitchY_ << std::endl;
      std::cout << " zBlade=" << zBlade << " gapBlade=" << gapBlade << std::endl;
      std::cout << " ZKapton=" << ZKapton << " ZBoundDet=" << ZBoundDet << std::endl;
      std::cout << " ZSiElectr=" << ZSiElectr << " ZCeramDet=" << ZCeramDet << std::endl;
      std::cout << " ZSiDetL=" << ZSiDetL_ << " ZSiDetR=" << ZSiDetR_ << std::endl;
    }
  ///////////////////////////////////////////////////////////////////



  
  
  //trackMode_ == "TrackProducerVar1FP420" ||
  //trackMode_ == "TrackProducerVar2FP420" ||
      if ( trackMode_ == "TrackProducerMaxAmplitudeFP420" ||
	   trackMode_ == "TrackProducerSophisticatedFP420" )  {
	finderParameters_ = new TrackProducerFP420(sn0_, pn0_, z420_, zD2_, zD3_,
						   pitchX_, pitchY_,
						   ZGapLDet_, ZSiStep_,
						   ZSiPlane_, ZSiDetL_, ZSiDetR_,
						   UseHalfPitchShiftInX_, UseHalfPitchShiftInY_,
						   dXX_,dYY_,chiCutX_,chiCutY_);
	validTrackerizer_ = true;
      } 
      else {
	std::cout << "ERROR:FP420TrackMain: No valid finder selected" << std::endl;
	validTrackerizer_ = false;
      }
}

FP420TrackMain::~FP420TrackMain() {
  if ( finderParameters_ != 0 ) {
    delete finderParameters_;
  }
}




void FP420TrackMain::run(const ClusterCollectionFP420 &input, TrackCollectionFP420 &toutput )
{

  if ( validTrackerizer_ ) {

    int number_detunits          = 0;
    int number_localelectroderechits = 0;
    /*
    for (int sector=1; sector<sn0_; sector++) {
      for (int zmodule=1; zmodule<pn0_; zmodule++) {
	for (int zside=1; zside<3; zside++) {
	  int sScale = 2*(pn0-1);
	  //      int index = FP420NumberingScheme::packFP420Index(det, zside, sector, zmodule);
	  // intindex is a continues numbering of FP420
	  int zScale=2;  unsigned int detID = sScale*(sector - 1)+zScale*(zmodule - 1)+zside;
	  ClusterMap.clear();
	  ClusterCollectionFP420::Range clusterRange;
	  clusterRange = input.get(detID);
	  ClusterCollectionFP420::ContainerIterator clusterRangeIteratorBegin = clusterRange.first;
	  ClusterCollectionFP420::ContainerIterator clusterRangeIteratorEnd   = clusterRange.second;
	  for ( ;sort_begin != sort_end; ++sort_begin ) {
	    ClusterMap.push_back(*sort_begin);
	  } // for
	  
	}//for
      }//for
    }//for
*/
    // get vector of detunit ids
    //    const std::vector<unsigned int> detIDs = input->detIDs();
    
  // to be used in put (besause of 0 in track collection for: 1) 1st track and 2) case of no track)
    // ignore 0, but to save info for 1st track record it second time on place 1   .

      bool first = true;
    // loop over detunits
      //    for (int sector=1; sector<5; sector++) {
	  ++number_detunits;
	  int StID = 1111;
	       std::vector<TrackFP420> collector;
// 	    vector<TrackFP420> collector;
		 collector.clear();

	  if ( trackMode_ == "TrackProducerMaxAmplitudeFP420") {
		 collector = finderParameters_->trackFinderMaxAmplitude(input); //std::vector<TrackFP420> collector;
	  }// if ( trackMode
	  /*
	  else if (trackMode_ == "TrackProducerVar1FP420" ) {
		 collector = finderParameters_->trackFinderVar1(input); //
	  }// if ( trackMode
	  else if (trackMode_ == "TrackProducerVar2FP420" ) {
		 collector = finderParameters_->trackFinderVar2(input); //
	  }// if ( trackMode
*/
	  else if (trackMode_ == "TrackProducerSophisticatedFP420" ) {
		 collector = finderParameters_->trackFinderSophisticated(input); //
	  }// if ( trackMode

	  //   if (collector.size()>0){
		 TrackCollectionFP420::Range inputRange;
		 inputRange.first = collector.begin();
		 inputRange.second = collector.end();
		 
		 if ( first ) {
		   // use it only if TrackCollectionFP420 is the TrackCollection of one event, otherwise, do not use (loose 1st cl. of 1st event only)
		   first = false;
		   unsigned int StID0 = 0;
		   toutput.put(inputRange,StID0); // !!! put into adress 0 for detID which will not be used never
		 } //if ( first ) 

		 // !!! put                                        !!! put
		 toutput.put(inputRange,StID);

		 number_localelectroderechits += collector.size();
		 //  } // if collector.size

#ifdef mytrackdebug0
	  std::cout << "FP420TrackMain: execution in mode " << trackMode_ << " generating " << number_localelectroderechits << " tracks in  " << number_detunits << " detectors" << std::endl; 
#endif



#ifdef mytrackdebug0
  //     check of access to the collector:
	  //	  std::vector<TrackFP420> collector;
	collector.clear();
	TrackCollectionFP420::Range outputRange;
	//	  int StID = 1111;
	outputRange = toutput.get(StID);
  // fill output in collector vector (for may be sorting? or other checks)
  TrackCollectionFP420::ContainerIterator sort_begin = outputRange.first;
  TrackCollectionFP420::ContainerIterator sort_end = outputRange.second;
  for ( ;sort_begin != sort_end; ++sort_begin ) {
    collector.push_back(*sort_begin);
  } // for
 std::cout <<" ===" << std::endl;
 std::cout <<" ===" << std::endl;
 std::cout <<"=======FP420TrackMain:check of re-new collector size = " << collector.size() << std::endl;
 std::cout <<" ===" << std::endl;
 std::cout <<" ===" << std::endl;
   vector<TrackFP420>::const_iterator simHitIter = collector.begin();
   vector<TrackFP420>::const_iterator simHitIterEnd = collector.end();
   // loop in #tracks
   for (;simHitIter != simHitIterEnd; ++simHitIter) {
     const TrackFP420 itrack = *simHitIter;
     
 std::cout << "FP420TrackMain:check: nclusterx = " << itrack.nclusterx() << "  nclustery = " << itrack.nclustery() << std::endl;
 std::cout << "  ax = " << itrack.ax() << "  bx = " << itrack.bx() << std::endl;
 std::cout << "  ay = " << itrack.ay() << "  by = " << itrack.by() << std::endl;
 std::cout << " chi2x= " << itrack.chi2x() << " chi2y= " << itrack.chi2y() << std::endl;
  std::cout <<" ===" << std::endl;
  std::cout <<" ===" << std::endl;
  std::cout <<" =======================" << std::endl;
   }

   //==================================

  //     end of check of access to the strip collection
 std::cout <<"=======            FP420TrackMain:                    end of check     " << std::endl;

#endif





  }// if ( validTrackerizer_



};
