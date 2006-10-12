
#include "Configuration/CSA06Skimming/interface/MCZll.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include <iostream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace edm;
using namespace std;


MCZll::MCZll(const edm::ParameterSet& iConfig) :
  label_(iConfig.getUntrackedParameter("moduleLabel",std::string("source"))), nEvents_(0), nAccepted_(0)
{
  leptonFlavour_ = iConfig.getUntrackedParameter<int>("leptonFlavour",11);
  leptonPtMin_ = iConfig.getUntrackedParameter<double>("leptonPtMin",5.);
  leptonPtMax_ = iConfig.getUntrackedParameter<double>("leptonPtMax",99999.);
  leptonEtaMin_ = iConfig.getUntrackedParameter<double>("leptonEtaMin",0.);
  leptonEtaMax_ = iConfig.getUntrackedParameter<double>("leptonEtaMax",2.7);
  zMassRange_.first = iConfig.getUntrackedParameter<double>("zMassMin",60.);
  zMassRange_.second = iConfig.getUntrackedParameter<double>("zMassMax",120.);
  //  filter_ = iConfig.getUntrackedParameter<bool>("filter",true);
  std::ostringstream str;
  str << "=========================================================\n" 
      << "Filter MCZll being constructed with parameters: " 
      << "\nleptonFlavour " << leptonFlavour_ 
      << "\nleptonPtMin " << leptonPtMin_
      << "\nleptonPtMax " << leptonPtMax_
      << "\nleptonEtaMin " << leptonEtaMin_
      << "\nleptonEtaMax " << leptonEtaMax_
      << "\nzMassMin " << zMassRange_.first
      << "\nzMassMax " << zMassRange_.second
      << "\n=========================================================" ;
  edm::LogVerbatim("MCZllInfo") <<  str.str() ;
  if (filter_)
    produces< HepMCProduct >();
}


MCZll::~MCZll()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void MCZll::endJob() 
{
  edm::LogVerbatim("MCZllInfo") << "================MCZll report========================================\n" 
	    << "Events read " << nEvents_ << " Events accepted " << nAccepted_ << "\nEfficiency " << ((double)nAccepted_)/((double)nEvents_) 
	    << "\n====================================================================" << std::endl;
}

// ------------ method called to skim the data  ------------
bool MCZll::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<HepMCProduct> bare_product(new HepMCProduct()); 

  nEvents_++;
  using namespace edm;
  bool accepted = false;
  Handle<HepMCProduct> evt;
  iEvent.getByLabel(label_, evt);
  HepMC::GenEvent * myGenEvent = new  HepMC::GenEvent(*(evt->GetEvent()));
  HepMC::GenEvent * zEvent = new HepMC::GenEvent();

  if (myGenEvent->signal_process_id() != 1) 
    return false;
  
  //found a prompt Z
  
  for ( HepMC::GenEvent::particle_iterator p = myGenEvent->particles_begin();
	p != myGenEvent->particles_end(); ++p ) 
    {
      if ( !accepted && ( (*p)->pdg_id() == 23 ) && (*p)->status() == 3 ) 
	{ 
	  accepted=true;
	  HepMC::GenVertex* zVertex = new HepMC::GenVertex();
	  HepMC::GenParticle* myZ= new HepMC::GenParticle(*(*p));
	  zVertex->add_particle_in(myZ);
	  //	  std::cout << (*p)->momentum().invariantMass() << std::endl;
	  if ((*p)->momentum().invariantMass() < zMassRange_.first || (*p)->momentum().invariantMass() > zMassRange_.second)
	    accepted = false; 
	  std::vector<HepMC::GenParticle*> children = (*p)->listChildren();
	  std::vector<HepMC::GenParticle*>::const_iterator aDaughter;
	  for (aDaughter = children.begin();aDaughter != children.end();aDaughter++)
	    {
	      HepMC::GenParticle* myDa= new HepMC::GenParticle(*(*aDaughter));
	      zVertex->add_particle_out(myDa);
	      if ((*aDaughter)->status() == 2)
		continue;
	      //	      (*aDaughter)->print();	      

	      if (! (abs((*aDaughter)->pdg_id()) == abs(leptonFlavour_)) )
		accepted = false; 
	      //		std::cout << (*aDaughter)->momentum().perp() << " " << (*aDaughter)->momentum().eta() << std::endl;
	      if ((*aDaughter)->momentum().perp() < leptonPtMin_)
		accepted = false; 
	      if ((*aDaughter)->momentum().perp() > leptonPtMax_)
		accepted = false; 
	      if (fabs((*aDaughter)->momentum().eta()) > leptonEtaMax_)
		accepted = false; 
	      if (fabs((*aDaughter)->momentum().eta()) < leptonEtaMin_)
		accepted = false; 
	    }
	  zEvent->add_vertex( zVertex );
	  if (accepted)
	    break;
	}

    } 
  

  if (accepted)
    {
      if(zEvent)  
	bare_product->addHepMCData(zEvent);
      if (filter_)
	iEvent.put(bare_product);
      nAccepted_++;
      //    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"<< std::endl;
      LogDebug("MCZll") << "Event " << iEvent.id().event()  << " accepted" << std::endl; 
      //    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"<< std::endl;
      //       myGenEvent->print(); 
      return true; 
    } 
  else 
    { 
      return false;
    }
  delete myGenEvent;   
  delete zEvent;
}

