#ifndef MCSingleParticleFilter_h
#define MCSingleParticleFilter_h
// -*- C++ -*-
//
// Package:    MCSingleParticleFilter
// Class:      MCSingleParticleFilter
// 
/* 

 Description: filter events based on the Pythia particleID and the Pt_hat

 Implementation: inherits from generic EDFilter
     
*/
//
// Original Author:  Filip Moortgat
//         Created:  Mon Sept 11 10:57:54 CET 2006
// $Id: MCSingleParticleFilter.h,v 1.0 2006/08/11 15:40:46 fmoortga Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


using namespace edm;
using namespace std;

//
// class decleration
//

class MCSingleParticleFilter : public edm::EDFilter {
   public:
      explicit MCSingleParticleFilter(const edm::ParameterSet&);
      ~MCSingleParticleFilter();


      virtual bool filter(Event&, const EventSetup&);
   private:
      // ----------member data ---------------------------
      
       std::string label_;
       std::vector<int> particleID;  
       std::vector<double> ptMin;
       std::vector<double> etaMin;  
       std::vector<double> etaMax;
       std::vector<int> status;
};
#endif
