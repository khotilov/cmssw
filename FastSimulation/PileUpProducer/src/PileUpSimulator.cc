#include "FastSimulation/PileUpProducer/interface/PileUpSimulator.h"
#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Particle/interface/RawParticle.h"

#include "HepMC/GenEvent.h"

//#include <iostream>

PileUpSimulator::PileUpSimulator(FSimEvent* aSimEvent) 
  : 
  mySimEvent(aSimEvent)  
{}
	       
PileUpSimulator::~PileUpSimulator() {}

void PileUpSimulator::produce(const HepMC::GenEvent* myGenEvent)
{

  // There might be no pile-up event to process (Poisson and/or flag)
  if ( !myGenEvent ) return;

  // Pile-up event iterator
  HepMC::GenEvent::vertex_const_iterator viter;
  HepMC::GenEvent::vertex_const_iterator vbegin = myGenEvent->vertices_begin();
  HepMC::GenEvent::vertex_const_iterator vend = myGenEvent->vertices_end();
  
  int ievt = 0;
  // Loop on all pile-up events
  for ( viter=vbegin; viter!=vend; ++viter ) { 

    // std::cout << "Vertex n0 " << ievt << std::endl;

    // The origin vertex
    HepMC::GenVertex* v = *viter;    
    XYZTLorentzVector smearedVertex =  
      XYZTLorentzVector(v->position().x(),v->position().y(),
			v->position().z(),v->position().t());

    // std::cout << "Vertex position " << smearedVertex << std::endl;

    // Add it to the FBaseSimEvent
    int mainVertex = mySimEvent->addSimVertex(smearedVertex);

    // Particles iterator
    HepMC::GenVertex::particles_out_const_iterator firstDaughterIt = v->particles_out_const_begin();
    HepMC::GenVertex::particles_out_const_iterator lastDaughterIt = v->particles_out_const_end();

    // Loop on particles
    for ( ; firstDaughterIt != lastDaughterIt ; ++firstDaughterIt ) {

      // A particle
      HepMC::GenParticle* daugh = *firstDaughterIt;
      RawParticle myPart(XYZTLorentzVector(daugh->momentum().px(),
					   daugh->momentum().py(),
					   daugh->momentum().pz(),
					   daugh->momentum().e()),
			 smearedVertex);
      
      // Add it to the FBaseSimEvent
      myPart.setID(daugh->pdg_id());
      // myPart.print();

      // Add the particle to the event (with a genpartIndex 
      // indicating the pileup event index)
      mySimEvent->addSimTrack(&myPart,mainVertex,-ievt-2);

      // End particle loop  
    }

    // Increment the number of pile-up events
    ++ievt;
    // End vertex loop
  }
  // Pile-up events now inserted

}
