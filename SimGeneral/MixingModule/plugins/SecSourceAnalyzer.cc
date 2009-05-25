// -*- C++ -*-
//
// Package:    SecSourceAnalyzer
// Class:      SecSourceAnalyzer
// 
/**\class SecSourceAnalyzer SecSourceAnalyzer.cc SecSource/SecSourceAnalyzer/src/SecSourceAnalyzer.cc

 Description: Get the data from the secondary source file using the getProductByTag method

 Implementation:
   
*/
//
// Original Author:  Emilia Lubenova Becheva
//         Created:  Wed Apr 22 16:54:31 CEST 2009
// $Id: SecSourceAnalyzer.cc,v 1.1 2009/04/24 13:35:34 ebecheva Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/CrossingFrame/interface/PCrossingFrame.h"

#include "FWCore/Framework/interface/InputSourceDescription.h"
#include "FWCore/Sources/interface/VectorInputSourceFactory.h"

#include "FWCore/Sources/interface/VectorInputSource.h"

#include "Mixing/Base/interface/PileUp.h"

#include "SecSourceAnalyzer.h"

//
// constructors and destructor
//
namespace edm
{
SecSourceAnalyzer::SecSourceAnalyzer(const edm::ParameterSet& iConfig)
:minBunch_(0),
 maxBunch_(0),
 tag_(InputTag())
{
   int minb = minBunch_;
   int maxb = maxBunch_;
   int averageNumber = 1;
   bool playback = false;
   
   input_.reset(new edm::PileUp(iConfig.getParameter<edm::ParameterSet>("input"),minb,maxb,averageNumber,playback));
      
   dataStep2_ = iConfig.getParameter<bool>("dataStep2");
   
   
   if (dataStep2_)
     // The data file contain the PCrossingFrame<SimTrack>
     label_   = iConfig.getParameter<edm::InputTag>("collPCF");
   else
     // The data file contain the SimTrack
     label_   = iConfig.getParameter<edm::InputTag>("collSimTrack");

}


SecSourceAnalyzer::~SecSourceAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SecSourceAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
   fileSeqNrs_.resize(maxBunch_-minBunch_+1);
   eventIDs_.resize(maxBunch_-minBunch_+1);
   nrEvents_.resize(maxBunch_-minBunch_+1);

   input_->readPileUp(pileup_[0],eventIDs_, fileSeqNrs_, nrEvents_);
   
   std::cout << "-> The std::vector<EventPrincipalVector> of the secondary source 'input' has been filled with " 
   	     << pileup_[0].size() << " element corresponding to " << maxBunch_-minBunch_+1 
	     << " bunch." << std::endl;
  
   // Run for one source (input) and the bunch 0
   Loop(pileup_[0][0]);
}


void SecSourceAnalyzer::Loop(const EventPrincipalVector& vec) {
    //
    // loop over events
    //
    std::cout <<"-> Loop over EventPrincipalVector wich size is "<< vec.size() << "." << std::endl;
    
    for (EventPrincipalVector::const_iterator it = vec.begin(); it != vec.end(); ++it) {
      std::cout <<"-> Get the event:  id " << (*it)->id() << std::endl;

      getBranches(&(**it));
    }
}


void  SecSourceAnalyzer::getBranches(EventPrincipal *ep)
  { 
    //std::cout << " ### SecSourceAnalyzer::getBranches - begin " << std::endl;
    std::cout << "-> dataStep2_ = " << dataStep2_ << std::endl;
    tag_ = InputTag(label_);
    
    std::cout << "-> Will try to get the branch with the tag : " << tag_ << std::endl;    
    std::cout << " and the EventPrincipal ep with a size = " << ep->size() << std::endl;
    
    if (!dataStep2_){
        // Get the SimTrack collection
        
	// default version changed to transmit vertexoffset
        boost::shared_ptr<Wrapper<std::vector<SimTrack> > const> shPtr =
        getProductByTag<std::vector<SimTrack> >(*ep, tag_);
    
        if (shPtr) 
		std::cout << "-> Could get SimTrack !" << std::endl;
	else 
		std::cout << "-> Could not get SimTrack !" << std::endl;

    }
    else{
        // Get the PCrossingFrame collection given as signal
    
        // default version changed to transmit vertexoffset
	tag_ = InputTag("CFwriter","g4SimHits");
        boost::shared_ptr<Wrapper<PCrossingFrame<SimTrack> > const> shPtr =
        getProductByTag<PCrossingFrame<SimTrack> >(*ep, tag_);
        
	if (shPtr) 
		std::cout << "-> Could get PCrossingFrame<SimTrack> !" << std::endl;
	else 
		std::cout << "-> Could not get PCrossingFrame<SimTrack> !" << std::endl;
	
    }
 
    
   //std::cout << " ### SecSourceAnalyzer::getBranches - end " << std::endl;
  }


// ------------ method called once each job just before starting event loop  ------------
void 
SecSourceAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SecSourceAnalyzer::endJob() {
  if (input_) input_->endJob();
}

}//edm
