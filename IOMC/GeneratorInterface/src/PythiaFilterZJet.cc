#include "IOMC/GeneratorInterface/interface/PythiaFilterZJet.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include <iostream>
#include<list>
#include<vector>
#include<cmath>

PythiaFilterZJet::PythiaFilterZJet(const edm::ParameterSet& iConfig) :
label_(iConfig.getUntrackedParameter("moduleLabel",std::string("source"))),
etaMuMax(iConfig.getUntrackedParameter<double>("MaxMuonEta", 2.5)),
ptZMin(iConfig.getUntrackedParameter<double>("MinZPt")),
ptZMax(iConfig.getUntrackedParameter<double>("MaxZPt")),
maxnumberofeventsinrun(iConfig.getUntrackedParameter<int>("MaxEvents",10)){ 
  
  theNumberOfSelected = 0;
}


PythiaFilterZJet::~PythiaFilterZJet(){}


// ------------ method called to produce the data  ------------
bool PythiaFilterZJet::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){

  if(theNumberOfSelected>=maxnumberofeventsinrun)   {
    throw cms::Exception("endJob")<<"we have reached the maximum number of events ";
  }

  bool accepted = false;
  edm::Handle<edm::HepMCProduct> evt;
  iEvent.getByLabel(label_, evt);

  const HepMC::GenEvent * myGenEvent = evt->GetEvent();
  std::vector<const HepMC::GenParticle *> mu;

  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();   p != myGenEvent->particles_end(); ++p ) {
   
    if ( std::abs((*p)->pdg_id())==13 && (*p)->status()==1 )
      mu.push_back(*p);
    if(mu.size()>1) break;
  }

  if(mu.size()>1){
    double ptZ= (mu[0]->momentum() + mu[1]->momentum()).perp();
    if (ptZ > ptZMin && ptZ < ptZMax  && 
	std::abs(mu[0]->momentum().eta()) < etaMuMax &&
	std::abs(mu[1]->momentum().eta()) < etaMuMax) 
      accepted=true;
  }

  if (accepted) {
    theNumberOfSelected++;
    return true; 
  }
  else return false;

}

