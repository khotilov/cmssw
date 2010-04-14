// -*- C++ -*-
//
// Package:    HSCP
// Class:      HSCPValidator
// 
/**\class HSCPValidator HSCPValidator.cc HSCPValidation/HSCPValidator/src/HSCPValidator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seth Cooper,27 1-024,+41227672342,
//         Created:  Wed Apr 14 14:27:52 CEST 2010
// $Id: HSCPValidator.h,v 1.1 2010/04/14 13:31:01 scooper Exp $
//
//
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TH2F.h"
//
// class declaration
//

class HSCPValidator : public edm::EDAnalyzer {
   public:
      explicit HSCPValidator(const edm::ParameterSet&);
      ~HSCPValidator();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      std::string intToString(int num);
      void makeGenPlots(const edm::Event& iEvent);
      void makeSimDigiPlotsECAL(const edm::Event& iEvent);
      void makeSimDigiPlotsRPC(const edm::Event& iEvent);

      // ----------member data ---------------------------
      bool doGenPlots_;
      bool doSimDigiPlots_;
      bool doRecoPlots_;

      // GEN section
      edm::InputTag label_;
      std::vector<int> particleIds_;
      int particleStatus_;
      std::map<int,int> particleIdsFoundMap_;
      TH1F* particleEtaHist_;
      TH1F* particlePhiHist_;
      TH1F* particlePHist_;
      TH1F* particlePtHist_;
      TH1F* particleMassHist_;
      TH1F* particleStatusHist_;
      TH1F* particleBetaHist_;
      TH1F* particleBetaInverseHist_;

      // SIM-DIGI section
      edm::InputTag ebSimHitTag_;
      edm::InputTag eeSimHitTag_;
      edm::InputTag simTrackTag_;
      edm::InputTag EBDigiCollection_;
      edm::InputTag EEDigiCollection_;
      edm::InputTag RPCRecHitTag_;
      edm::ESHandle <RPCGeometry> rpcGeo;
      // ECAL      
      TH1F* simHitsEcalEnergyHistEB_;
      TH1F* simHitsEcalTimeHistEB_;
      TH1F* simHitsEcalNumHistEB_;
      TH2F* simHitsEcalEnergyVsTimeHistEB_;
      TH1F* simHitsEcalEnergyHistEE_;
      TH1F* simHitsEcalTimeHistEE_;
      TH1F* simHitsEcalNumHistEE_;
      TH2F* simHitsEcalEnergyVsTimeHistEE_;
      TH1F* simHitsEcalDigiMatchEnergyHistEB_;
      TH1F* simHitsEcalDigiMatchTimeHistEB_;
      TH2F* simHitsEcalDigiMatchEnergyVsTimeHistEB_;
      TH1F* simHitsEcalDigiMatchEnergyHistEE_;
      TH1F* simHitsEcalDigiMatchTimeHistEE_;
      TH2F* simHitsEcalDigiMatchEnergyVsTimeHistEE_;
      TH1F* simHitsEcalDigiMatchIEtaHist_;
      TH1F* simHitsEcalDigiMatchIPhiHist_;
      TH1F* digisEcalNumHistEB_;
      TH1F* digisEcalNumHistEE_;
      TH2F* digiOccupancyMapEB_;
      TH2F* digiOccupancyMapEEP_;
      TH2F* digiOccupancyMapEEM_;
      // RPC
      TH1F* residualsRPCRecHitSimDigis_;
      TH1F* efficiencyRPCRecHitSimDigis_;
      TH1F* cluSizeDistribution_; 
      TH1F* rpcTimeOfFlightBarrel_[6];       
      TH1F* rpcBXBarrel_[6];       
      TH1F* rpcTimeOfFlightEndCap_[3];       
      TH1F* rpcBXEndCap_[3];       
};
