// -*- C++ -*-
//
// Package:    JetTracksAssociator
// Class:      JetTracksAssociator
// 
/**\class JetTracksAssociator JetTracksAssociator.cc RecoBTag/JetTracksAssociator/src/JetTracksAssociator.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andrea Rizzi
//         Created:  Wed Apr 12 11:12:49 CEST 2006
// $Id: JetTracksAssociator.cc,v 1.3 2006/05/19 15:24:56 arizzi Exp $
//
//


// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BTauReco/interface/JetTracksAssociation.h"



#include "DataFormats/Math/interface/Vector3D.h"


//Math
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"

//
// class decleration
//

#include <iostream>
using namespace std;


using namespace reco;
class JetTracksAssociator : public edm::EDProducer {
   public:
      explicit JetTracksAssociator(const edm::ParameterSet&);
      ~JetTracksAssociator();


      virtual void produce(edm::Event&, const edm::EventSetup&);
   private:
     JetTracksAssociationCollection * associate(const edm::Handle<CaloJetCollection> & jets,
                                const edm::Handle<TrackCollection> & tracks ) const;
     bool trackIsInJetCone ( const Jet & jet , const Track & track ) const;
      // ----------member data ---------------------------
     double m_deltaRCut;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JetTracksAssociator::JetTracksAssociator(const edm::ParameterSet& iConfig)
{
   produces<reco::JetTracksAssociationCollection>();
//produces<float>();
m_deltaRCut = 0.5;
}


JetTracksAssociator::~JetTracksAssociator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
JetTracksAssociator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   Handle<CaloJetCollection> jets;
   iEvent.getByLabel("mcone5",jets);
   Handle<TrackCollection> tracks;
   iEvent.getByType(tracks);
   
   std::auto_ptr<JetTracksAssociationCollection> jetTracks(associate(jets,tracks));
   iEvent.put(jetTracks);
}

JetTracksAssociationCollection * JetTracksAssociator::associate(const edm::Handle<CaloJetCollection> & jets,
                                const edm::Handle<TrackCollection> & tracks ) const
{
JetTracksAssociationCollection * outputCollection = new JetTracksAssociationCollection();
 //loop on jets and associate
 for(size_t j=0; j < jets->size() ; j++)
 {
    cout << "Jet #" << j << " px= " << (*jets)[j].px() << " py= " << (*jets)[j].py() << " pz= " << (*jets)[j].pz() << endl;
  for(size_t t=0; t < tracks->size() ; t++) {
      bool inside= trackIsInJetCone((*jets)[j],(*tracks)[t]);
      cout << "Track #" << t << " " << (*tracks)[t].momentum() << " is inside ? " << inside << endl;
     if(inside)outputCollection->insert(edm::Ref<CaloJetCollection>(jets,j),edm::Ref<TrackCollection>(tracks,t));
   }
 } 
//dummy insert

return outputCollection;
}

bool JetTracksAssociator::trackIsInJetCone ( const Jet & jet , const Track & track ) const {

  double deltaR ;

  // get direction info from jet/track and compute deltaR
  math::XYZVector jet3Vec   (jet.px(),jet.py(),jet.pz()) ;
  math::XYZVector trackMomentum = track.momentum() ;
  deltaR = ROOT::Math::VectorUtil::DeltaR(jet3Vec,trackMomentum ) ;
   cout << "deltaR: " << deltaR << " " ;  
  return ( deltaR < m_deltaRCut ) ;
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetTracksAssociator)
