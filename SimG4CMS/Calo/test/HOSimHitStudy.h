#ifndef SimG4CMS_HOSimHitStudy_H
#define SimG4CMS_HOSimHitStudy_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include <TH1F.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

class HOSimHitStudy: public edm::EDAnalyzer{

public:

  HOSimHitStudy(const edm::ParameterSet& ps);
  ~HOSimHitStudy() {}

protected:

  void beginJob () {}
  void endJob   () {}
  void analyze  (const edm::Event& e, const edm::EventSetup& c);

  void analyzeHits  ();

private:

  std::string           sourceLabel, g4Label, hitLab[2];
  std::vector<PCaloHit> ecalHits, hcalHits;
  double                maxEnergy, scaleEB, scaleHB, scaleHO;
  bool                  scheme_, print_;
  double                tcut_;
  TH1F                  *hit_[3],  *time_[3], *edepTW_[3], *edepTWT_[3];
  TH1F                  *edep_[3], *hitTow_[3], *eneInc_, *etaInc_, *phiInc_;
  TH1F                  *edepT_[3], *eEB_, *eEBHB_, *eEBHBHO_, *eEBHBHOT_;
  TH1F                  *edepZon_[3], *edepZonT_[3], *eEBT_, *eEBHBT_;
  TProfile              *eHO1_, *eHO1T_, *eHO17_, *eHO17T_, *eHO18_, *eHO18T_;
  TProfile              *nHO1_, *nHO1T_;
  TProfile2D            *eHO2_, *eHO2T_, *nHO2_, *nHO2T_;
  double                eInc, etaInc, phiInc;
};

#endif
