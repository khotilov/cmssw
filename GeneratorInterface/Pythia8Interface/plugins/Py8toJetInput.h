#ifndef gen_Py8toJetInput_h
#define gen_Py8toJetInput_h

#include "FastJet3.h" // Py8 overhead on top of FastJets package
#include "Event.h"

//#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
//#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"
//#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"

namespace lhef {

   class LHEEvent;

}

using namespace Pythia8;

class Py8toJetInput
{

   public:
   
      Py8toJetInput(): fJetEtaMax(10.) {}
      ~Py8toJetInput() {}
      
      virtual const std::vector<fastjet::PseudoJet> fillJetAlgoInput( const Event&, const Event&, 
                                                                      const lhef::LHEEvent* lhee=0,
                                                                      const std::vector<int>* partonList=0 );
      void setJetEtaMax( double max ) { fJetEtaMax=max; return; }
      
   protected:

      enum partonTypes { ID_TOP=6, ID_GLUON=21, ID_PHOTON=22 };
      double fJetEtaMax;
            
      int getAncestor( int, const Event&, const Event& );
      
      std::vector<fastjet::PseudoJet> fJetInput;

};

class Py8toJetInputHEPEVT : public Py8toJetInput
{

   public:
   
      Py8toJetInputHEPEVT() {}
      ~Py8toJetInputHEPEVT() {}
      
      const std::vector<fastjet::PseudoJet> fillJetAlgoInput( const Event&, const Event&, 
                                                              const lhef::LHEEvent*,
                                                              const std::vector<int>* partonList=0 ); 
};

#endif
