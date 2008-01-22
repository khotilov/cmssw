#ifndef ECALRAWTODIGI_ECALDIGIDISPLAY_h
#define ECALRAWTODIGI_ECALDIGIDISPLAY_h

#include <memory>
#include <vector>
#include <string>

#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CaloOnlineTools/EcalTools/interface/EcalFedMap.h"

// class declaration

class EcalDigiDisplay : public edm::EDAnalyzer {
 public:
  //Constractor
  EcalDigiDisplay(const edm::ParameterSet& ps);
  //Distractor
  ~EcalDigiDisplay();
  
 private:
  virtual void analyze( const edm::Event & e, const  edm::EventSetup& c);
  virtual void beginJob(const edm::EventSetup& c);
  virtual void endJob();
  
 protected:
  void readEBDigis (edm::Handle<EBDigiCollection> digis, int Mode);
  void readEEDigis (edm::Handle<EEDigiCollection> digis, int Mode);
  void readPNDigis (edm::Handle<EcalPnDiodeDigiCollection> PNs, int Mode);
  
  EcalFedMap* fedMap;

  std::string ebDigiCollection_;
  std::string eeDigiCollection_;
  std::string digiProducer_;

  std::vector<int> requestedFeds_;
  std::vector<std::string> requestedEbs_;

  bool inputIsOk;
  bool cryDigi;
  bool ttDigi;
  bool pnDigi;
  
  int mode;
 
  std::vector<int> listChannels;
  std::vector<int> listTowers;
  std::vector<int> listPns;
  

};
#endif
