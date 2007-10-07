#ifndef TtHadLRSignalSelObservables_h
#define TtHadLRSignalSelObservables_h


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// General C++ stuff
#include <iostream>
#include <string>
#include <vector>

// ROOT classes
#include "TLorentzVector.h"
#include "TVector.h"
#include "TVectorD.h"

#include "TMatrix.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDSym.h"
#include "TMatrixTSym.h"

//own code
#include "AnalysisDataFormats/TopObjects/interface/TtHadEvtSolution.h"
#include "AnalysisDataFormats/TopObjects/interface/TopJet.h"

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
      bool operator()( TopJet j1, TopJet j2 ) const
      {
	return j1.getRecJet().et() > j2.getRecJet().et();
      }
    };
    
    CompareET EtComparator;
    
    // compare two jets in bdisc
    struct CompareBdisc {
      bool operator()( TopJet j1, TopJet j2 ) const
      {
	return j1.getBDiscriminator("trackCountingJetTags") > j2.getBDiscriminator("trackCountingJetTags");
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
