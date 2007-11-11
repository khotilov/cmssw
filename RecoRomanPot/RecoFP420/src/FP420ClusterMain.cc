///////////////////////////////////////////////////////////////////////////////
// File: FP420ClusterMain.cc
// Date: 12.2006
// Description: FP420ClusterMain for FP420
// Modifications: 
///////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <iostream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoRomanPot/RecoFP420/interface/FP420ClusterMain.h"
#include "DataFormats/FP420Digi/interface/HDigiFP420.h"
#include "DataFormats/FP420Cluster/interface/ClusterFP420.h"
#include "RecoRomanPot/RecoFP420/interface/ClusterProducerFP420.h"
#include "RecoRomanPot/RecoFP420/interface/ClusterNoiseFP420.h"

#include "CLHEP/Random/RandFlat.h"

using namespace std;


FP420ClusterMain::FP420ClusterMain(const edm::ParameterSet& conf, int sn, int pn):conf_(conf),sn0(sn),pn0(pn)  { 
  
  verbosity         = conf_.getUntrackedParameter<int>("VerbosityLevel");
  ElectronPerADC_   = conf_.getParameter<double>("ElectronFP420PerAdc");
  clusterMode_      = conf_.getParameter<std::string>("ClusterModeFP420");
  ChannelThreshold  = conf_.getParameter<double>("ChannelFP420Threshold");//6
  SeedThreshold     = conf_.getParameter<double>("SeedFP420Threshold");//7
  ClusterThreshold  = conf_.getParameter<double>("ClusterFP420Threshold");//7
  MaxVoidsInCluster = conf_.getParameter<int>("MaxVoidsFP420InCluster");//1
  
  if (verbosity > 0) {
    std::cout << "FP420ClusterMain constructor: ElectronPerADC = " << ElectronPerADC_ << std::endl;
    std::cout << " clusterMode = " << clusterMode_ << std::endl;
    std::cout << " ChannelThreshold = " << ChannelThreshold << std::endl;
    std::cout << " SeedThreshold = " << SeedThreshold << std::endl;
    std::cout << " ClusterThreshold = " << ClusterThreshold << std::endl;
    std::cout << " MaxVoidsInCluster = " << MaxVoidsInCluster << std::endl;
  }
  xytype=2;// only X types of planes
  ENC_ = 960.;    // 
  Thick300 = 0.300;
  BadElectrodeProbability_ = 0.002;
  //UseNoiseBadElectrodeFlagFromDB_ = true;
  UseNoiseBadElectrodeFlagFromDB_ = false;
  //
  // pitches and ldriftes:
  //
  ldriftX = 0.050; 
  ldriftY = 0.050;// was 0.040
  pitchY= 0.050;// was 0.040
  pitchX= 0.050;
  moduleThicknessY = 0.250; // mm
  moduleThicknessX = 0.250; // mm
  numStripsY = 201;        // Y plate number of strips:200*0.050=10mm (xytype=1)
  numStripsX = 401;        // X plate number of strips:400*0.050=20mm (xytype=2)
  numStripsYW = 51;        // Y plate number of W strips:50 *0.400=20mm (xytype=1) - W have ortogonal projection
  numStripsXW = 26;        // X plate number of W strips:25 *0.400=10mm (xytype=2) - W have ortogonal projection
  
  //  sn0 = 4;
  //  pn0 = 9;
  if (verbosity > 1) {
    std::cout << "FP420ClusterMain constructor: sn0 = " << sn0 << " pn0=" << pn0 << std::endl;
    std::cout << "FP420ClusterMain constructor: ENC = " << ENC_ << std::endl;
    std::cout << " Thick300 = " << Thick300 << std::endl;
    std::cout << " BadElectrodeProbability = " << BadElectrodeProbability_ << std::endl;
    std::cout << " ldriftX = " << ldriftX << " ldriftY = " << ldriftY << std::endl;
    std::cout << " pitchY = " << pitchY << " pitchX = " << pitchX << std::endl;
    std::cout << " numStripsY = " << numStripsY << " numStripsX = " << numStripsX << std::endl;
    std::cout << " moduleThicknessY = " << moduleThicknessY << " moduleThicknessX = " << moduleThicknessX << std::endl;
  }
  
  if (UseNoiseBadElectrodeFlagFromDB_==false){	  
    if (verbosity > 0) {
      std::cout << "using a SingleNoiseValue and good electrode flags" << std::endl;
    }
  } else {
    if (verbosity > 0) {
      std::cout << "using Noise and BadElectrode flags accessed from DB" << std::endl;
    }
  }
  
  if ( clusterMode_ == "ClusterProducerFP420" ) {
    
    
    //    ChannelThreshold    = 6.0;// was 2.6.0 7 18
    //   SeedThreshold       = 7.0;//was 3.7.0  8 20
    // ClusterThreshold    = 7.0;// was 2. 7.0 8 20
    // MaxVoidsInCluster   = 1;	
    threeThreshold_ = new ClusterProducerFP420(ChannelThreshold,
					       SeedThreshold,
					       ClusterThreshold,
					       MaxVoidsInCluster);
    validClusterizer_ = true;
  } else {
    std::cout << "ERROR:FP420ClusterMain: No valid clusterizer selected" << std::endl;
    validClusterizer_ = false;
  }
}

FP420ClusterMain::~FP420ClusterMain() {
  if ( threeThreshold_ != 0 ) {
    delete threeThreshold_;
  }
}


//void FP420ClusterMain::run(const DigiCollectionFP420 *input, ClusterCollectionFP420 &soutput,
//			   const std::vector<ClusterNoiseFP420>& electrodnoise)
void FP420ClusterMain::run(edm::Handle<DigiCollectionFP420> &input, std::auto_ptr<ClusterCollectionFP420> &soutput,
			   std::vector<ClusterNoiseFP420>& electrodnoise)

{
  // unpack from iu:
  //  int  sScale = 20, zScale=2;
  //  int  sector = (iu-1)/sScale + 1 ;
  //  int  zmodule = (iu - (sector - 1)*sScale - 1) /zScale + 1 ;
  //  int  zside = iu - (sector - 1)*sScale - (zmodule - 1)*zScale ;
  
    if (verbosity > 0) {
      std::cout << "FP420ClusterMain: OK1" << std::endl;
    }
  if ( validClusterizer_ ) {
    
    int number_detunits          = 0;
    int number_localelectroderechits = 0;
    
    // get vector of detunit ids
    //    const std::vector<unsigned int> detIDs = input->detIDs();
    
    // to be used in put (besause of 0 in cluster collection for: 1) 1st cluster and 2) case of no cluster)
    // ignore 0, but to save info for 1st cluster record it second time on place 1   .
    
    bool first = true;
    // loop over detunits
    for (int sector=1; sector<sn0; sector++) {
      for (int zmodule=1; zmodule<pn0; zmodule++) {
	for (int zside=1; zside<3; zside++) {
	  int sScale = 2*(pn0-1);
	  //	int det= 1;
	  //      int index = CaloNumberingPacker::packCastorIndex(det,zside, sector, zmodule);
	  //      int index = FP420NumberingScheme::packFP420Index(det, zside, sector, zmodule);
	  // intindex is a continues numbering of FP420
	  int zScale=2;  unsigned int detID = sScale*(sector - 1)+zScale*(zmodule - 1)+zside;
	  // int zScale=10;	unsigned int intindex = sScale*(sector - 1)+zScale*(zside - 1)+zmodule;
	  if (verbosity > 0) {
	    std::cout << " FP420ClusterMain:1 run loop   index no  iu = " << detID  << std::endl;
	  }	  
	  // Y:
	  if (xytype ==1) {
	    numStrips = numStripsY*numStripsYW;  
	    moduleThickness = moduleThicknessY; 
	    pitch= pitchY;
	    ldrift = ldriftX;
	  }
	  // X:
	  if (xytype ==2) {
	    numStrips = numStripsX*numStripsXW;  
	    moduleThickness = moduleThicknessX; 
	    pitch= pitchX;
	    ldrift = ldriftY;
	  }
	  
	  
	  //    for ( std::vector<unsigned int>::const_iterator detunit_iterator = detIDs.begin(); detunit_iterator != detIDs.end(); ++detunit_iterator ) {
	  //      unsigned int detID = *detunit_iterator;
	  ++number_detunits;
	  
	  //   .
	  //   GET DIGI collection  !!!!
	  //   .
	  //	  const DigiCollectionFP420::Range digiRange = input->get(detID);
	  DigiCollectionFP420::Range digiRange;
	  	digiRange = input->get(detID);
		//digiRange = input.get(detID);
	  
	  
	  DigiCollectionFP420::ContainerIterator digiRangeIteratorBegin = digiRange.first;
	  DigiCollectionFP420::ContainerIterator digiRangeIteratorEnd   = digiRange.second;
	  
	  if ( clusterMode_ == "ClusterProducerFP420" ) {
	    
	    std::vector<ClusterFP420> collector;
	    // 	    vector<ClusterFP420> collector;
	    
	    if (UseNoiseBadElectrodeFlagFromDB_==false){	  
	      
	      //Case of SingleValueNoise flags for all electrodes of a Detector
	      
	      
	      //float noise = ENC_*ldrift/Thick300/ElectronPerADC_;//Noise is proportional to charge collection path 
	      float noise = ENC_*moduleThickness/Thick300/ElectronPerADC_;//Noise is proportional to moduleThickness
	      
	      //vector<float> noiseVec(numElectrodes,noise);	    
	      //Construct a ElectrodNoiseVector ( in order to be compliant with the DB access)
	      ElectrodNoiseVector vnoise;
	      ClusterNoiseFP420::ElectrodData theElectrodData;       	   
	      
	      for(int electrode=0; electrode < numStrips; ++electrode){
		//   discard  randomly  bad  electrode with probability BadElectrodeProbability_
		bool badFlag= RandFlat::shoot(1.) < BadElectrodeProbability_ ? true : false;
		theElectrodData.setData(noise,badFlag);
		vnoise.push_back(theElectrodData);// fill vector vnoise
	      } // for
	      
	      //clusterizeDetUnit   or    clusterizeDetUnitPixels      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      collector.clear();
	      //	      std::vector<ClusterFP420> collector;
	      //  collector = threeThreshold_->clusterizeDetUnit(digiRangeIteratorBegin,digiRangeIteratorEnd,detID,vnoise);
	      collector = threeThreshold_->clusterizeDetUnitPixels(digiRangeIteratorBegin,digiRangeIteratorEnd,detID,vnoise,xytype);
	      
	      
	    } else {
	      //Case of Noise and BadElectrode flags access from DB
	      /*
		const ElectrodNoiseVector& vnoise = electrodnoise->getElectrodNoiseVector(detID);
		
		if (vnoise.size() <= 0) {
		std::cout << "WARNING requested Noise Vector for detID " << detID << " that isn't in map " << std::endl; 
		continue;
		}
		collector.clear();
		collector = threeThreshold_->clusterizeDetUnit(digiRangeIteratorBegin,digiRangeIteratorEnd,detID,vnoise);
	      */
	      
	      
	    }// if (UseNoiseBadElectrodeFlagFromDB
	    
	    if (collector.size()>0){
	      ClusterCollectionFP420::Range inputRange;
	      inputRange.first = collector.begin();
	      inputRange.second = collector.end();
	      
	      if ( first ) {
		// use it only if ClusterCollectionFP420 is the ClusterCollection of one event, otherwise, do not use (loose 1st cl. of 1st event only)
	      	first = false;
		unsigned int detID0 = 0;
		soutput->put(inputRange,detID0); // !!! put into adress 0 for detID which will not be used never
	      } //if ( first ) 
	      
	      // !!! put
	      soutput->put(inputRange,detID);
	      
	      number_localelectroderechits += collector.size();
	    } // if
	  }// if ( clusterMode
	  if (verbosity > 0) {
	    std::cout << "[FP420ClusterMain] execution in mode " << clusterMode_ << " generating " << number_localelectroderechits << " ClusterFP420s in " << number_detunits << " DetUnits." << std::endl; 
	  }
	  
	}//for
      }//for
    }//for
    
    if (verbosity > 0) {
      //     check of access to the collector:
      for (int sector=1; sector<sn0; sector++) {
	for (int zmodule=1; zmodule<pn0; zmodule++) {
	  for (int zside=1; zside<3; zside++) {
	    int sScale = 2*(pn0-1);
	    int zScale=2;  unsigned int iu = sScale*(sector - 1)+zScale*(zmodule - 1)+zside;
	    // int zScale=10;	unsigned int intindex = sScale*(sector - 1)+zScale*(zside - 1)+zmodule;
	    std::vector<ClusterFP420> collector;
	    collector.clear();
	    ClusterCollectionFP420::Range outputRange;
	    outputRange = soutput->get(iu);
	    // fill output in collector vector (for may be sorting? or other checks)
	    ClusterCollectionFP420::ContainerIterator sort_begin = outputRange.first;
	    ClusterCollectionFP420::ContainerIterator sort_end = outputRange.second;
	    for ( ;sort_begin != sort_end; ++sort_begin ) {
	      collector.push_back(*sort_begin);
	    } // for
	    std::cout <<" ===" << std::endl;
	    std::cout <<" ===" << std::endl;
	    std::cout <<"=======FP420ClusterMain:check of re-new collector size = " << collector.size() << std::endl;
	    std::cout <<" ===" << std::endl;
	    std::cout <<" ===" << std::endl;
	    vector<ClusterFP420>::const_iterator simHitIter = collector.begin();
	    vector<ClusterFP420>::const_iterator simHitIterEnd = collector.end();
	    // loop in #clusters
	    for (;simHitIter != simHitIterEnd; ++simHitIter) {
	      const ClusterFP420 icluster = *simHitIter;
	      
	      std::cout << "FP420ClusterMain:check: zside = " << zside << std::endl;
	      std::cout << " firstStrip = " << icluster.firstStrip() << "  barycenter = " << icluster.barycenter() << std::endl;
	      std::cout << " size of cluster= " << icluster.amplitudes().size() << std::endl;
	      for(unsigned int i = 0; i < icluster.amplitudes().size(); i++ ) {
		std::cout <<  "i = " << i << "   amplitudes = "  << icluster.amplitudes()[i] << std::endl;
	      }
	      std::cout <<" ===" << std::endl;
	      std::cout <<" ===" << std::endl;
	      std::cout <<" =======================" << std::endl;
	    }
	    /* 
	       for (DigitalMapType::const_iterator i=collector.begin(); i!=collector.end(); i++) {
	       std::cout << "DigitizerFP420:check: HDigiFP420::  " << std::endl;
	       std::cout << " strip number is as (*i).first  = " << (*i).first << "  adc is in (*i).second  = " << (*i).second << std::endl;
	       }
	    */
	    
	    //==================================
	    
	  }   // for
	}   // for
      }   // for
      
      //     end of check of access to the strip collection
      std::cout <<"=======            FP420ClusterMain:                    end of check     " << std::endl;
      
    }
    
    
    
    
  }// if ( validClusterizer_
}
