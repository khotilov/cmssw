
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

//#include "CondFormats/DTObjects/interface/DTMtime.h"
//#include "CondFormats/DataRecord/interface/DTMtimeRcd.h"

using namespace std;

namespace edmtest {
  class DTMtimePrint : public edm::EDAnalyzer
  {
  public:
    explicit  DTMtimePrint(edm::ParameterSet const& p);
    explicit  DTMtimePrint(int i) ;
    virtual ~ DTMtimePrint();
    virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  private:
  };
}
