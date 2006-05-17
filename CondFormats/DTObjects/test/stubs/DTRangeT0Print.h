
/*----------------------------------------------------------------------

Toy EDProducers and EDProducts for testing purposes only.

----------------------------------------------------------------------*/

//#include <stdexcept>
//#include <string>
//#include <iostream>
//#include <map>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "CondFormats/DTObjects/interface/DTRangeT0.h"
//#include "CondFormats/DataRecord/interface/DTRangeT0Rcd.h"

using namespace std;

namespace edmtest {
  class DTRangeT0Print : public edm::EDAnalyzer
  {
  public:
    explicit  DTRangeT0Print(edm::ParameterSet const& p);
    explicit  DTRangeT0Print(int i) ;
    virtual ~ DTRangeT0Print();
    virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  private:
  };
}
