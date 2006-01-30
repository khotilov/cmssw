#ifndef PYTHIAFILTER_h
#define PYTHIAFILTER_h
// -*- C++ -*-
//
// Package:    PythiaFilter
// Class:      PythiaFilter
// 
/**\class PythiaFilter PythiaFilter.cc IOMC/PythiaFilter/src/PythiaFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Filip Moortgat
//         Created:  Mon Jan 23 14:57:54 CET 2006
// $Id$
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

class PythiaFilter : public edm::EDFilter {
   public:
      explicit PythiaFilter(const edm::ParameterSet&);
      ~PythiaFilter();


      virtual bool filter(const Event&, const EventSetup&);
   private:
      // ----------member data ---------------------------
      
       std::string label_;
       double minptcut;
      
};
#endif
