#ifndef TtHadLRSignalSelObservables_h
#define TtHadLRSignalSelObservables_h

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <iostream>
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TVectorD.h"
#include "TMatrix.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDSym.h"
#include "TMatrixTSym.h"

#include "AnalysisDataFormats/TopObjects/interface/TtHadEvtSolution.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

const double PI=3.14159265;

using namespace std;

class TtHadLRSignalSelObservables{
  
 public:
  
  TtHadLRSignalSelObservables();
  ~TtHadLRSignalSelObservables();	
  
  void  operator()(TtHadEvtSolution&);
  
 private:
  
  // compare two jets in ET
  struct CompareET {
    bool operator()( pat::Jet j1, pat::Jet j2 ) const
    {
      return j1.recJet().et() > j2.recJet().et();
    }
  };
  
  CompareET EtComparator;
  
  // compare two jets in bdisc
  struct CompareBdisc {
    bool operator()( pat::Jet j1, pat::Jet j2 ) const
    {
      return j1.bDiscriminator("trackCountingJetTags") > j2.bDiscriminator("trackCountingJetTags");
    }
  };
  
  CompareBdisc BdiscComparator;
  
  // compare two double
  struct CompareDouble {
    bool operator()( double j1, double j2 ) const
    {
      return j1 > j2 ;
    }
  };
  
  CompareDouble dComparator;
  
  vector<pair<unsigned int,double> > evtselectVarVal;
  
};

#endif
