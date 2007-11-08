#include "SimG4Core/CustomPhysics/interface/CustomPhysicsList.h"
#include "SimG4Core/CustomPhysics/interface/CustomParticleFactory.h"
#include "SimG4Core/CustomPhysics/interface/DummyChargeFlipProcess.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4Decay.hh"
#include "G4MultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ProcessManager.hh"


#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "SimG4Core/CustomPhysics/interface/FullModelHadronicProcess.hh"
#include "SimG4Core/CustomPhysics/interface/ToyModelHadronicProcess.hh"
 

CustomPhysicsList::CustomPhysicsList(std::string name, const edm::ParameterSet & p)  :  G4VPhysicsConstructor(name)
{
  
  myConfig = p;
  edm::FileInPath fp = p.getParameter<edm::FileInPath>("particlesDef");
  particleDefFilePath = fp.fullPath();
  edm::LogInfo("")<<"Path for custom particle definition file: "<<particleDefFilePath<<std::endl;
  myHelper = 0;
  
 }

CustomPhysicsList::~CustomPhysicsList() {
  delete myHelper;
}
 
void CustomPhysicsList::ConstructParticle(){
  CustomParticleFactory::loadCustomParticles(particleDefFilePath);     
}
 
void CustomPhysicsList::ConstructProcess() {
  addCustomPhysics();
}
 
void CustomPhysicsList::addCustomPhysics(){
    LogDebug("") << " CustomPhysics: adding CustomPhysics processes  " <<std::endl;
    theParticleIterator->reset();
/*    while((*theParticleIterator)())
    {	
	G4ParticleDefinition* particle = theParticleIterator->value();
	if(CustomParticleFactory::isCustomParticle(particle))
	{
	    edm::LogInfo("") << particle->GetParticleName() << " is Custom" << std::endl;
	    G4ProcessManager* pmanager = particle->GetProcessManager();
	    if(pmanager)
	    { 
              if(particle->GetParticleType()=="rhadron") {
		pmanager->AddDiscreteProcess(new ToyModelHadronicProcess);
	      }
              if(particle->GetPDGCharge()/eplus != 0)
		{ 
//		  pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
//		  pmanager->AddProcess(new G4hIonisation,       -1, 2,2);

		    edm::LogInfo("") << "    Processes for charged particle added." << std::endl;
                }
                else edm::LogInfo("") << "   It is neutral!!" << std::endl;
               
            }
	    else
		edm::LogInfo("") << "   No pmanager" << std::endl;
	}*/
	

  while((*theParticleIterator)())    {
      int i = 0;
      G4ParticleDefinition* particle = theParticleIterator->value();
      CustomParticle* cp = dynamic_cast<CustomParticle*>(particle);
      if(CustomParticleFactory::isCustomParticle(particle))
        {
          LogDebug("") << particle->GetParticleName()<<", "<<particle->GetPDGEncoding()
		       << " is Custom. Mass is "
		       <<particle->GetPDGMass()/GeV
		       <<" GeV."<<G4endl;
          if(cp->GetCloud()!=0)
            {
              LogDebug("")<<"Cloud mass is "
			  <<cp->GetCloud()->GetPDGMass()/GeV
			  <<" GeV. Spectator mass is "
			  <<static_cast<CustomParticle*>(particle)->GetSpectator()->GetPDGMass()/GeV
			  <<" GeV."       << G4endl;
            }
          G4ProcessManager* pmanager = particle->GetProcessManager();
          if(pmanager)
            {
              if(cp!=0) {
                //              i++;
                //              edm::LogInfo("")<<"1"<<std::endl;
		if(particle->GetParticleType()=="rhadron" || particle->GetParticleType()=="mesonino" || particle->GetParticleType() == "sbaryon" ){
		  if(!myHelper) myHelper = new G4ProcessHelper(myConfig);
		  pmanager->AddDiscreteProcess(new FullModelHadronicProcess(myHelper)); //GHEISHA
		}
                //              edm::LogInfo("")<<"2"<<std::endl;
                //              pmanager->AddDiscreteProcess(new ToyModelHadronicProcess()); //Toy model
              }
              if(particle->GetPDGCharge()/eplus != 0)
                {
                  //              edm::LogInfo("")<<"3"<<std::endl;
		  pmanager->AddProcess(new G4MultipleScattering,-1, 1,i+1);
                  //              edm::LogInfo("")<<"4"<<std::endl;
		  pmanager->AddProcess(new G4hIonisation,       -1, 2,i+2);
                  //              edm::LogInfo("")<<"5"<<std::endl;
                  //              G4cout << "    Processes for charged particle added." << G4endl;
                  //              G4cout << " i = "<<i<<G4endl;		 
                }
              //              else G4cout << "   It is neutral!!" << G4endl;
              //pmanager->DumpInfo();
            }
          else      LogDebug("") << "   No pmanager" << G4endl;
        }





	
    }
}


void CustomPhysicsList::setupRHadronPhycis(G4ParticleDefinition* particle){

  //    LogDebug("")<<"Configuring rHadron: "
  //	<<cp->

  CustomParticle* cp = dynamic_cast<CustomParticle*>(particle);
  if(cp->GetCloud()!=0) 
    LogDebug("")<<"Cloud mass is "
		<<cp->GetCloud()->GetPDGMass()/GeV
		<<" GeV. Spectator mass is "
		<<static_cast<CustomParticle*>(particle)->GetSpectator()->GetPDGMass()/GeV
		<<" GeV."       << G4endl;
  
  G4ProcessManager* pmanager = particle->GetProcessManager();
  if(pmanager){
    if(!myHelper) myHelper = new G4ProcessHelper(myConfig);
    pmanager->AddDiscreteProcess(new FullModelHadronicProcess(myHelper)); //GHEISHA
    if(particle->GetPDGCharge()/eplus != 0){
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4hIonisation,       -1, 2,2);
    }
  }
  else      LogDebug("") << "   No pmanager" << G4endl;
}
					       

void CustomPhysicsList::setupSUSYPhycis(G4ParticleDefinition* particle){

  CustomParticle* cp = dynamic_cast<CustomParticle*>(particle);
  G4ProcessManager* pmanager = particle->GetProcessManager();
  if(pmanager){
    pmanager->AddProcess(new G4Decay,1, 1,1);
    if(particle->GetPDGCharge()/eplus != 0){
      pmanager->AddProcess(new G4MultipleScattering,-1, 2,2);
      pmanager->AddProcess(new G4hIonisation,       -1, 3,3);
    }
  }
  else      LogDebug("") << "   No pmanager" << G4endl;
}
