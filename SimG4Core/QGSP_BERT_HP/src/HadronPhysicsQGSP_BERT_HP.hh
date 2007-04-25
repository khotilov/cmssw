#ifndef HadronPhysicsQGSP_BERT_HP_h
#define HadronPhysicsQGSP_BERT_HP_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4StoppingHadronBuilder.hh"
#include "G4MiscLHEPBuilder.hh"

#include "G4PiKBuilder.hh"
#include "G4LEPPiKBuilder.hh"
#include "G4QGSPPiKBuilder.hh"
#include "G4BertiniPiKBuilder.hh"

#include "G4ProtonBuilder.hh"
#include "G4LEPProtonBuilder.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4BertiniProtonBuilder.hh"

#include "G4NeutronBuilder.hh"
#include "G4LEPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
#include "G4NeutronHPBuilder.hh"

class HadronPhysicsQGSP_BERT_HP : public G4VPhysicsConstructor
{
  public: 
    HadronPhysicsQGSP_BERT_HP(const G4String& name ="hadron");
    virtual ~HadronPhysicsQGSP_BERT_HP();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4NeutronBuilder theNeutrons;
    G4LEPNeutronBuilder theLEPNeutron;
    G4QGSPNeutronBuilder theQGSPNeutron;
    //
    G4BertiniNeutronBuilder theBertiniNeutron;
    G4NeutronHPBuilder theHPNeutron;
  
    G4PiKBuilder thePiK;
    G4LEPPiKBuilder theLEPPiK;
    G4QGSPPiKBuilder theQGSPPiK;
    //
    G4BertiniPiKBuilder theBertiniPiK;
    
    G4ProtonBuilder thePro;
    G4LEPProtonBuilder theLEPPro;
    G4QGSPProtonBuilder theQGSPPro;  
    //
    G4BertiniProtonBuilder theBertiniPro;
    
    G4MiscLHEPBuilder theMiscLHEP;
    G4StoppingHadronBuilder theStoppingHadron;
    //G4HadronQEDBuilder theHadronQED;
};

// 2002 by J.P. Wellisch

#endif

