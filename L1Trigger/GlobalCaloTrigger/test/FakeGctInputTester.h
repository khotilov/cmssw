
// -*- C++ -*-
//
// Package:    FakeGctInputTester
// Class:      FakeGctInputTester
// 
/**\class FakeGctInputTester FakeGctInputTester.h L1Trigger/GlobalCalotrigger/src/FakeGctInputTester.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jim Brooke
//         Created:  Tue Nov 21 14:49:14 CET 2006
// $Id: FakeGctInputTester.h,v 1.4 2006/12/08 09:33:11 jbrooke Exp $
//
//



#ifndef FAKEGCTINPUTTESTER_H
#define FAKEGCTINPUTTESTER_H


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>


class TH1F;
class TFile;

class FakeGctInputTester : public edm::EDAnalyzer {
public:
  explicit FakeGctInputTester(const edm::ParameterSet&);
  ~FakeGctInputTester();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------

  std::string hFileName_;
   TFile* hFile_;

   TH1F* isoEmDEta_;
   TH1F* isoEmDPhi_;

   TH1F* nonIsoEmDEta_;
   TH1F* nonIsoEmDPhi_;

   TH1F* jetDEta_;
   TH1F* jetDPhi_;


};

#endif
