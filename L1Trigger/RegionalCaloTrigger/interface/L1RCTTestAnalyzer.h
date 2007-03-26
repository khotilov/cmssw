// -*- C++ -*-
//
// Package:    L1RCTTestAnalyzer
// Class:      L1RCTTestAnalyzer
//
/**\class L1RCTTestAnalyzer L1RCTTestAnalyzer.cc src/L1RCTTestAnalyzer/src/L1RCTTestAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/47
//         Created:  Thu Jul 13 21:38:08 CEST 2006
// $Id: L1RCTTestAnalyzer.h,v 1.4 2007/03/03 09:23:44 dasu Exp $
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

//
// class declaration
//

class L1RCTTestAnalyzer : public edm::EDAnalyzer {
public:
  explicit L1RCTTestAnalyzer(const edm::ParameterSet&);
  ~L1RCTTestAnalyzer();
  
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
private:
  // ----------member data ---------------------------
  std::string m_HepMCProduct;
  bool showEmCands;
  bool showRegionSums;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
