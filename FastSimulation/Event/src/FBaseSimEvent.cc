//HepMC Headers
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

// CMSSW Sim headers
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"


//FAMOS Headers
#include "FastSimulation/Event/interface/FBaseSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/Event/interface/FSimVertex.h"
#include "FastSimulation/Event/interface/KineParticleFilter.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FastSimulation/Event/interface/GaussianPrimaryVertexGenerator.h"
#include "FastSimulation/Event/interface/FlatPrimaryVertexGenerator.h"
#include "FastSimulation/Event/interface/NoPrimaryVertexGenerator.h"
//#include "FastSimulation/Utilities/interface/Histos.h"

using namespace std;
//using namespace edm;
using namespace HepMC;
//using namespace CLHEP;
using namespace HepPDT;

// system include
#include <iostream>
#include <iomanip>
#include <map>

FBaseSimEvent::FBaseSimEvent(const edm::ParameterSet& kine) 
  : random(0)
{

  theVertexGenerator = new NoPrimaryVertexGenerator();

  // Initialize the vectors of particles and vertices
  theGenParticles = new vector<GenParticle*>(); 
  theSimTracks = new vector<FSimTrack>;
  theSimVertices = new vector<FSimVertex>;
  theChargedTracks = new vector<unsigned>();

  // Reserve some size to avoid mutiple copies
  theSimTracks->reserve(20000);
  theSimVertices->reserve(20000);
  theGenParticles->reserve(20000);
  theChargedTracks->reserve(20000);

  // Initialize the Particle filter
  myFilter = new KineParticleFilter(kine);

}

FBaseSimEvent::FBaseSimEvent(const edm::ParameterSet& vtx,
			     const edm::ParameterSet& kine,
			     const RandomEngine* engine) 
  : theVertexGenerator(0), random(engine)
{


  // Initialize the vertex generator
  string vtxType = vtx.getParameter<string>("type");
  if ( vtxType == "Gaussian" ) 
    theVertexGenerator = new GaussianPrimaryVertexGenerator(vtx,random);
  else if ( vtxType == "Flat" ) 
    theVertexGenerator = new FlatPrimaryVertexGenerator(vtx,random);
  else
    theVertexGenerator = new NoPrimaryVertexGenerator();

  // Initialize the vectors of particles and vertices
  theGenParticles = new vector<GenParticle*>(); 
  theSimTracks = new vector<FSimTrack>;
  theSimVertices = new vector<FSimVertex>;
  theChargedTracks = new vector<unsigned>();

  // Reserve some size to avoid mutiple copies
  theSimTracks->reserve(20000);
  theSimVertices->reserve(20000);
  theGenParticles->reserve(20000);
  theChargedTracks->reserve(20000);

  // Initialize the Particle filter
  myFilter = new KineParticleFilter(kine);

  // Get the Famos Histos pointer
  //  myHistos = Histos::instance();

  // Initialize a few histograms
  /*
  myHistos->book("hvtx",100,-0.1,0.1);
  myHistos->book("hvty",100,-0.1,0.1);
  myHistos->book("hvtz",100,-500.,500.);
  */
}
 
FBaseSimEvent::~FBaseSimEvent(){

  clear();
  delete theGenParticles;
  delete theSimTracks;
  delete theSimVertices;
  delete theChargedTracks;
  delete myFilter;

  //Write the histograms
  //  myHistos->put("histos.root");
  //  delete myHistos;
}

void 
FBaseSimEvent::initializePdt(const DefaultConfig::ParticleDataTable* aPdt) { 

  pdt = aPdt; 

}

const DefaultConfig::ParticleDataTable*
FBaseSimEvent::theTable() const {
  return pdt;
}

void
FBaseSimEvent::fill(const HepMC::GenEvent& myGenEvent) {
  
  // Clear old vectors
  clear();

  // printMCTruth(myGenEvent);

  // Fill the event with the stable particles of the GenEvent 
  // (and their mother), in the order defined by the original 
  // particle barcodes

  // Add the particles in the FSimEvent
  addParticles(myGenEvent);

  // Check 
  // for ( unsigned i=0; i<nTracks(); ++i ) cout << embdTrack(i) << endl;

  // for ( unsigned i=0; i<nVertices(); ++i ) cout << embdVertex(i) << endl;
}

void
FBaseSimEvent::fill(const std::vector<SimTrack>& simTracks, 
		    const std::vector<SimVertex>& simVertices) {

  // Watch out there ! A SimVertex is in mm (stupid), 
  //            while a FSimVertex is in cm (clever).
  
  clear();

  unsigned nVtx = simVertices.size();
  unsigned nTks = simTracks.size();

  // Empty event, do nothin'
  if ( nVtx == 0 ) return;

  // Two arrays for internal use.
  vector<int> myVertices(nVtx,-1);
  vector<int> myTracks(nTks,-1);

  // create a map associating geant particle id and position in the 
  // event SimTrack vector
  
  map<unsigned, unsigned> geantToIndex;
  for( unsigned it=0; it<simTracks.size(); ++it ) {
    geantToIndex[ simTracks[it].trackId() ] = it;
  }  

  // Set the main vertex for the kine particle filter
  // SimVertices were in mm until 110_pre2
  //  HepLorentzVector primaryVertex = simVertices[0].position()/10.;
  // SImVertices are now in cm
  HepLorentzVector primaryVertex = simVertices[0].position();
  myFilter->setMainVertex(primaryVertex);
  // Add the main vertex to the list.
  addSimVertex(myFilter->vertex());
  myVertices[0] = 0;

  for( unsigned trackId=0; trackId<nTks; ++trackId ) {

    // The track
    const SimTrack& track = simTracks[trackId];

    // The origin vertex
    int vertexId = track.vertIndex();
    const SimVertex& vertex = simVertices[vertexId];

    // The mother track 
    int motherId = -1;
    if( !vertex.noParent() ) { // there is a parent to this vertex
      // geant id of the mother
      unsigned motherGeantId =   vertex.parentIndex(); 
      map<unsigned, unsigned >::iterator association  
	= geantToIndex.find( motherGeantId );
      if(association != geantToIndex.end() )
	motherId = association->second;
    }
    int originId = motherId == - 1 ? -1 : myTracks[motherId];

    // Add the vertex (if it does not already exist!)
    if ( myVertices[vertexId] == -1 )
      myVertices[vertexId] = addSimVertex(vertex.position(),originId); 

    // Add the track (with protection for brem'ing electrons)
    int motherType = motherId == -1 ? 0 : simTracks[motherId].type();

    if ( abs(motherType) != 11 || motherType != track.type() ) {
      // SimVertices were in mm until 110_pre2
      // RawParticle part(track.momentum(), vertex.position()/10.);
      // SImVertices are now in cm
      RawParticle part(track.momentum(), vertex.position());
      part.setID(track.type()); 
      myTracks[trackId] = addSimTrack(&part,myVertices[vertexId]);
    } else {
      myTracks[trackId] = myTracks[motherId];
    }
    
  }

  // Now loop over the remaining end vertices !
  for( unsigned vertexId=0; vertexId<nVtx; ++vertexId ) {

    // if the vertex is already saved, just ignore.
    if ( myVertices[vertexId] != -1 ) continue;

    // The yet unused vertex
    const SimVertex& vertex = simVertices[vertexId];

    // The mother track 
    int motherId = -1;
    if( !vertex.noParent() ) { // there is a parent to this vertex

      // geant id of the mother
      unsigned motherGeantId =   vertex.parentIndex(); 
      map<unsigned, unsigned >::iterator association  
	= geantToIndex.find( motherGeantId );
      if(association != geantToIndex.end() )
	motherId = association->second;
    }
    int originId = motherId == - 1 ? -1 : myTracks[motherId];

    // Add the vertex
    myVertices[vertexId] = 
      addSimVertex(vertex.position(),originId);
  }

  // Finally, propagate all particles to the calorimeters
  BaseParticlePropagator myPart;
  HepLorentzVector mom;
  HepLorentzVector pos;
  double enele, enegam;

  // Loop over the tracks
  for( int fsimi=0; fsimi < (int)nTracks() ; ++fsimi) {

    FSimTrack& myTrack = track(fsimi);
    mom = myTrack.momentum();
    enele = mom.e();
    // Special treatment for electrons to account for bremstrahlung photons
    if ( abs(myTrack.type()) == 11 && myTrack.nDaughters() > 0 ) { 
      for ( int idaugh=0; idaugh<myTrack.nDaughters(); ++idaugh ) {
	// Subtract photon energy
	enegam = myTrack.daughter(idaugh).momentum().e();
	enele -= enegam;
	// Give the proper direction (assuming collinear emission)
	mom = myTrack.daughter(idaugh).momentum()*enele/enegam;
	//	mom -= myTrack.daughter(idaugh).momentum();
	pos =  myTrack.daughter(idaugh).vertex().position();
      }
    } else {
      pos = myTrack.vertex().position();
    }

    // The particle to be propagated
    myPart = BaseParticlePropagator(RawParticle(mom,pos),0.,0.,4.);
    myPart.setCharge(myTrack.charge());

    // Propagate to Preshower layer 1
    myPart.propagateToPreshowerLayer1(false);
    if ( myTrack.notYetToEndVertex(myPart.vertex()) && myPart.getSuccess()>0 )
    myTrack.setLayer1(myPart,myPart.getSuccess());
  
    // Propagate to Preshower Layer 2 
    myPart.propagateToPreshowerLayer2(false);
    if ( myTrack.notYetToEndVertex(myPart.vertex()) && myPart.getSuccess()>0 )
    myTrack.setLayer2(myPart,myPart.getSuccess());

    // Propagate to Ecal Endcap
    myPart.propagateToEcalEntrance(false);
    if ( myTrack.notYetToEndVertex(myPart.vertex()) )
      myTrack.setEcal(myPart,myPart.getSuccess());
    
    // Propagate to HCAL entrance
    myPart.propagateToHcalEntrance(false);
    if ( myTrack.notYetToEndVertex(myPart.vertex()) )
      myTrack.setHcal(myPart,myPart.getSuccess());
    
    // Propagate to VFCAL entrance
    myPart.propagateToVFcalEntrance(false);
    if ( myTrack.notYetToEndVertex(myPart.vertex()) )
      myTrack.setVFcal(myPart,myPart.getSuccess());

  }

}


void
FBaseSimEvent::addParticles(const HepMC::GenEvent& myGenEvent) {

  /// Some internal array to work with.
  map<const GenParticle*,int> myGenVertices;

  // If no particles, no work to be done !
  if ( myGenEvent.particles_empty() ) return;

  // Are there particles in the FSimEvent already ? 
  int offset = nGenParts();

  // Primary vertex (already smeared by the SmearedVtx module)
  GenVertex* primaryVertex = *(myGenEvent.vertices_begin());
  HepLorentzVector primaryVertexPosition(primaryVertex->position().x()/10.,
					 primaryVertex->position().y()/10.,
					 primaryVertex->position().z()/10.,
					 primaryVertex->position().t()/10.);

  // Smear the main vertex if needed
  HepLorentzVector smearedVertex; 
  if ( primaryVertex->point3d().mag() < 1E-10 ) {
    theVertexGenerator->generate();
    smearedVertex = HepLorentzVector(*theVertexGenerator);
  }

  // Fill Histos
  /* 
  myHistos->fill("hvtx",smearedVertex.x());
  myHistos->fill("hvty",smearedVertex.y());
  myHistos->fill("hvtz",smearedVertex.z());
  cout << smearedVertex << endl;
  */

  // Set the main vertex
  myFilter->setMainVertex(primaryVertexPosition+smearedVertex);

  // This is the smeared main vertex
  //  GenVertex* mainVertex = new GenVertex(myFilter.vertex());
  //  addSimVertex(mainVertex);
  int mainVertex = addSimVertex(myFilter->vertex());

  // Loop on the particles of the generated event
  for ( HepMC::GenEvent::particle_const_iterator 
	  piter  = myGenEvent.particles_begin();
	  piter != myGenEvent.particles_end(); 
	++piter ) {

    // This is the generated particle pointer - for the signal event only
    GenParticle* p = *piter;
    if  ( !offset && nGenParts() != 20000 ) theGenParticles->push_back(p);
    // Keep only: 
    // 1) Stable particles
    bool testStable = p->status()==1;

    // 2) or particles with stable daughters
    bool testDaugh = false;
    GenVertex::particles_out_const_iterator firstDaughterIt = 
      p->end_vertex()->particles_out_const_begin();
    GenVertex::particles_out_const_iterator lastDaughterIt = 
      p->end_vertex()->particles_out_const_end();
    for ( ; firstDaughterIt != lastDaughterIt ; ++firstDaughterIt ) {
      GenParticle* daugh = *firstDaughterIt;
      if ( daugh->status()==1 ) {
	testDaugh=true;
	break;
      }
    }

    // 3) or particles that fly more than one micron.
    double dist = 0.;
    if ( p->production_vertex() ) {
      HepLorentzVector 
	productionVertexPosition(p->production_vertex()->position().x()/10.,
				 p->production_vertex()->position().y()/10.,
				 p->production_vertex()->position().z()/10.,
				 p->production_vertex()->position().t()/10.);
      dist = (primaryVertexPosition-productionVertexPosition).vect().mag();
    }
    bool testDecay = ( dist > 0.0001 ) ? true : false; 

    // Save the corresponding particle and vertices
    if ( testStable || testDaugh || testDecay ) {

      // The particle is the copy of the original
      //      GenParticle* part = 
      //	new GenParticle(p->momentum(),
      //			p->pdg_id(),
      //			p->status(),
      //			p->flow(),
      //			p->polarization());

      // The origin vertex is either the primary, 
      // or the end vertex of the mother, if saved
      //      GenVertex* originVertex = 
      //	p->mother() &&  
      //	myGenVertices.find(p->mother()) != myGenVertices.end() ? 
      //      	originVertex = myGenVertices[p->mother()] : mainVertex;
      const GenParticle* mother = p->production_vertex() ?
	*(p->production_vertex()->particles_in_const_begin()) : 0;
      int originVertex = 
	mother &&  
	myGenVertices.find(mother) != myGenVertices.end() ? 
      	myGenVertices[mother] : mainVertex;
      
      HepLorentzVector momentum(p->momentum().px(),
				p->momentum().py(),
				p->momentum().pz(),
				p->momentum().e());
      RawParticle part(momentum, vertex(originVertex).position());
      part.setID(p->pdg_id());

      // Add the particle to the event and to the various lists
      int theTrack = addSimTrack(&part,originVertex, nGenParts()-1-offset);

      // It there an end vertex ?
      if ( !p->end_vertex() ) continue; 

      // If yes, create it
      //      GenVertex* decayVertex = 
      //	new GenVertex(p->end_vertex()->position()
      //		      +mainVertex->position());

      // Add the vertex to the event and to the various lists
      HepLorentzVector decayVertex = 
	HepLorentzVector(p->end_vertex()->position().x()/10.,
			 p->end_vertex()->position().y()/10.,
			 p->end_vertex()->position().z()/10.,
			 p->end_vertex()->position().t()/10.) +
	vertex(mainVertex).position();
      int theVertex = addSimVertex(decayVertex,theTrack);

      // And record it for later use 
      //      if ( theVertex != -1 ) myGenVertices[p] = decayVertex;
      if ( theVertex != -1 ) myGenVertices[p] = theVertex;

      // There we are !

    }
  }

  //  printMCTruth(*this);

}

int 
FBaseSimEvent::addSimTrack(const RawParticle* p, int iv, int ig) { 
  
  // Check that the particle is in the Famos "acceptance"
  if ( !myFilter->accept(p) ) return -1;

  // An increasing barcode, corresponding to the list index
  //  part->suggest_barcode(nTracks()+1);
  
  // Protection on the number of tracks
  int trackId = nTracks();
  if ( trackId == 20000 ) {
    cout << "FastSimulation:FBaseSimEvent - try to store more than 20000 tracks" << endl;
    return -1;
  }

  // Attach the particle to the origin vertex, and to the mother
  vertex(iv).addDaughter(trackId);
  if ( !vertex(iv).noParent() )  
    track(vertex(iv).parent().id()).addDaughter(trackId);

  // Attach the vertex to the event (inoccuous if the vertex exists)
  // add_vertex(originVertex);
  
  // Some persistent information for the users
  //  mySimTracks->push_back(SimTrack(p->pid(),*p,iv,ig)); 

  // Some transient information for FAMOS internal use
  //  theSimTracks->push_back(FSimTrack(mySimTracks->size()-1,this));
  theSimTracks->push_back(FSimTrack(p,iv,ig,trackId,this));

  return trackId;

}

int
FBaseSimEvent::addSimVertex(const HepLorentzVector& v,int im) {
  
  // Check that the vertex is in the Famos "acceptance"
  if ( !myFilter->accept(RawParticle(HepLorentzVector(),v)) ) return -1;

  // An increasing -barcode, corresponding to the list index
  //  decayVertex->suggest_barcode(-nVertices()-1);

  // Attach the vertex to the event (inoccuous if the vertex exists)
  //  add_vertex(decayVertex);

  //  if ( im!=-1 ) decayVertex->add_particle_in(track(im).me());
  int vertexId = nVertices();
  if ( vertexId == 20000 ) {
    cout << "FastSimulation:FBaseSimEvent - try to store more than 20000 vertices" << endl;
    return -1;
  }

  // Attach the end vertex to the particle (if accepted)
  if ( im !=-1 ) track(im).setEndVertex(vertexId);

  // Some persistent information for the users
  //  mySimVertices->push_back(SimVertex(v.vect(),v.e(),im));

  // Some transient information for FAMOS internal use
  theSimVertices->push_back(FSimVertex(v,im,vertexId,this));

  return vertexId;

}

void
FBaseSimEvent::printMCTruth(const HepMC::GenEvent& myGenEvent) {
  
  cout << "Id  Gen Name       eta    phi     pT     E    Vtx1   " 
       << " x      y      z   " 
       << "Moth  Vtx2  eta   phi     R      Z   Da1  Da2 Ecal?" << endl;

  for ( HepMC::GenEvent::particle_const_iterator 
	  piter  = myGenEvent.particles_begin();
	  piter != myGenEvent.particles_end(); 
	++piter ) {
  //  for ( int i=1; i !=myGenEvent.particles_size(); ++i ) { 
    
    HepMC::GenParticle* p = *piter;
     /* */
     //     const std::string name = (*p)->particledata().name();
    int partId = p->pdg_id();
    std::string name;

    if ( pdt->particle(ParticleID(partId)) !=0 ) {
      name = (pdt->particle(ParticleID(partId)))->name();
    } else {
      name = "none";
    }
       
    HepLorentzVector momentum1(p->momentum().px(),
			       p->momentum().py(),
			       p->momentum().pz(),
			       p->momentum().e());
    int vertexId1 = 0;
    if ( !p->production_vertex() ) continue;
    Hep3Vector vertex1 = Hep3Vector(p->production_vertex()->position().x()/10.,
				    p->production_vertex()->position().y()/10.,
				    p->production_vertex()->position().z()/10.);
    vertexId1 = p->production_vertex()->barcode();
    
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::right, ios::adjustfield);
    
    cout << setw(4) << p->barcode()-1 << " " 
	 << name;
    
    for(unsigned int k=0;k<11-name.length() && k<12; k++) cout << " ";  
    
    double eta = momentum1.eta();
    if ( eta > +10. ) eta = +10.;
    if ( eta < -10. ) eta = -10.;
    cout << setw(6) << setprecision(2) << eta << " " 
	 << setw(6) << setprecision(2) << momentum1.phi() << " " 
	 << setw(7) << setprecision(2) << momentum1.perp() << " " 
	 << setw(7) << setprecision(2) << momentum1.e() << " " 
	 << setw(4) << vertexId1 << " " 
	 << setw(6) << setprecision(1) << vertex1.x() << " " 
	 << setw(6) << setprecision(1) << vertex1.y() << " " 
	 << setw(6) << setprecision(1) << vertex1.z() << " ";

    const GenParticle* mother = 
      *(p->production_vertex()->particles_in_const_begin());

    if ( mother )
      cout << setw(4) << mother->barcode() << " ";
    else 
      cout << "     " ;
    
    if ( p->end_vertex() ) {  
      HepLorentzVector vertex2(p->end_vertex()->position().x()/10.,
			       p->end_vertex()->position().y()/10.,
			       p->end_vertex()->position().z()/10.,
			       p->end_vertex()->position().t()/10.);
      int vertexId2 = p->end_vertex()->barcode();
      
      GenParticle* firstDaughter 
	= *(p->end_vertex()->particles_out_const_begin());
      GenParticle* lastDaughter 
	= *(p->end_vertex()->particles_out_const_end()--);
      cout << setw(4) << vertexId2 << " "
	   << setw(6) << setprecision(2) << vertex2.eta() << " " 
	   << setw(6) << setprecision(2) << vertex2.phi() << " " 
	   << setw(5) << setprecision(1) << vertex2.perp() << " " 
	   << setw(6) << setprecision(1) << vertex2.z() << " "
	   << setw(4) << firstDaughter->barcode() << " " 
	   << setw(4) << lastDaughter->barcode() << " ";
    }
    cout << endl;

  }

}

void
FBaseSimEvent::print() const {
  //  for(int i=0; i<(int)genparts()->size(); ++i)
  //    cout << i << " " << embdGenpart(i) << endl << endl;

  cout << "  Id  Gen Name       eta    phi     pT     E    Vtx1   " 
       << " x      y      z   " 
       << "Moth  Vtx2  eta   phi     R      Z   Daughters Ecal?" << endl;

  for( int i=0; i<(int)nTracks(); i++ ) 
    cout << track(i) << endl;
}

void 
FBaseSimEvent::clear() {

  // Clear the vectors
  theGenParticles->clear();
  theSimTracks->clear();
  theSimVertices->clear();
  theChargedTracks->clear();

}

void 
FBaseSimEvent::addChargedTrack(int id) { 
  theChargedTracks->push_back(id);
}

static FSimTrack oTrack;
FSimTrack&
FBaseSimEvent::track(int id) const { 
  return  id>=0 && id<(int)theSimTracks->size() ? 
    (*theSimTracks)[id] : oTrack; }

static FSimVertex oVertex;
FSimVertex&
FBaseSimEvent::vertex(int id) const { 
  return   id>=0 && id<(int)theSimVertices->size() ? 
    (*theSimVertices)[id] : oVertex; }

int
FBaseSimEvent::chargedTrack(int id) const {
  if (id>=0 && id<(int)theChargedTracks->size()) return (*theChargedTracks)[id]; 
  else return -1;
}

unsigned int 
FBaseSimEvent::nTracks() const {
  return theSimTracks->size();
}

unsigned int 
FBaseSimEvent::nVertices() const { 
  return theSimVertices->size();
}

unsigned int 
FBaseSimEvent::nGenParts() const {
  return theGenParticles->size();
}

unsigned int 
FBaseSimEvent::nChargedTracks() const {
  return theChargedTracks->size();
}


static  const SimVertex zeroVertex;
const SimVertex & 
FBaseSimEvent::embdVertex(int i) const { 
  if (i>=0 && i<(int)theSimVertices->size()) 
    return (*theSimVertices)[i]; 
  else 
    return zeroVertex;
}

static  const SimTrack zeroTrack;
const SimTrack & 
FBaseSimEvent::embdTrack(int i) const { 
  if (i>=0 && i<(int)theSimTracks->size()) 
    return (*theSimTracks)[i]; 
  else 
    return zeroTrack;
}

const HepMC::GenParticle* 
FBaseSimEvent::embdGenpart(int i) const { 
  if (i>=0 && i<(int)theGenParticles->size()) 
    return (*theGenParticles)[i]; 
  else 
    return 0;
}

std::vector<FSimTrack>* 
FBaseSimEvent::tracks() const { return theSimTracks; }

std::vector<FSimVertex>*
FBaseSimEvent::vertices() const { return theSimVertices; }

std::vector<HepMC::GenParticle*>* 
FBaseSimEvent::genparts() const { return theGenParticles; }


