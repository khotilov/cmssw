

#ifndef _ECALRAWTODIGIDEV_H_ 
#define _ECALRAWTODIGIDEV_H_ 

/*
 *\ Class EcalRawToDigi
 *
 * This class takes care of unpacking ECAL's raw data info
 *
 * \file EcalRawToDigi.h
 *
 * $Date: 2007/03/20 01:28:39 $
 * $Revision: 1.1.2.1 $
 * \author N. Almeida
 * \author G. Franzoni
 *
*/

#include <iostream>                                 

#include "ECALUnpackerException.h"
#include "DCCRawDataDefinitions.h"

#include <DataFormats/FEDRawData/interface/FEDRawData.h>
#include <DataFormats/FEDRawData/interface/FEDNumbering.h>
#include <DataFormats/FEDRawData/interface/FEDRawDataCollection.h>
#include <DataFormats/EcalDigi/interface/EcalDigiCollections.h>
#include <DataFormats/EcalRawData/interface/EcalRawDataCollections.h>

#include <FWCore/Framework/interface/Handle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <sys/time.h>

class EcalElectronicsMapper;
class EcalElectronicsMapping;
class DCCDataUnpacker;

using namespace std;
using namespace edm;

class EcalRawToDigiDev : public EDProducer{

 public:
  /**
   * Class constructor
   */
  explicit EcalRawToDigiDev(const edm::ParameterSet& ps);
  
  /**
   * Functions that are called by framework at each event
   */
  virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
  /**
   * Class destructor
   */
  virtual ~EcalRawToDigiDev();
  
 private:

  //list of FEDs to unpack
  std::vector<int> fedUnpackList_;
  
  uint numbXtalTSamples_;
  uint numbTriggerTSamples_;
  
  bool headerUnpacking_;
  bool srpUnpacking_;
  bool tccUnpacking_;
  bool feUnpacking_;
  bool memUnpacking_;
  bool first_;
  bool put_;

  //an electronics mapper class 
  EcalElectronicsMapper * myMap_;

  EcalElectronicsMapping * mmm_; 
 
  //Ecal unpacker
  DCCDataUnpacker * theUnpacker_;

   
  
  uint nevts_; // NA: for testing
  double  RUNNING_TIME_, SETUP_TIME_;
  
  
};



#endif

