// File: BMixingModule.cc
// Description:  see BMixingModule.h
// Author:  Ursula Berthon, LLR Palaiseau, Bill Tanenbaum
//
//--------------------------------------------

#include "Mixing/Base/interface/BMixingModule.h"
#include "FWCore/Utilities/interface/GetPassID.h"
#include "FWCore/Utilities/interface/GetReleaseVersion.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Provenance/interface/ModuleDescription.h"

using namespace std;

int edm::BMixingModule::vertexoffset = 0;
const unsigned int edm::BMixingModule::maxNbSources =4;
namespace
{
  boost::shared_ptr<edm::PileUp>
  maybeMakePileUp(edm::ParameterSet const& ps, std::string SourceName)
  {
    boost::shared_ptr<edm::PileUp> pileup; // value to be returned
    // Make sure we have a parameter named 'input'.
    vector<string> names = ps.getParameterNames();
    if (find(names.begin(), names.end(), std::string(SourceName))
	!= names.end())
      {
	pileup.reset(new edm::PileUp(ps.getParameter<edm::ParameterSet>(SourceName)));
      }

    return pileup;
  }
}

namespace edm {

  // Constructor 
  BMixingModule::BMixingModule(const edm::ParameterSet& pset) :
    bunchSpace_(pset.getParameter<int>("bunchspace")),
    checktof_(pset.getUntrackedParameter<bool>("checktof",true)),
    input_(maybeMakePileUp(pset,"input")),
    cosmics_(maybeMakePileUp(pset,"cosmics")),
    beamHalo_p_(maybeMakePileUp(pset,"beamhalo_plus")),
    beamHalo_m_(maybeMakePileUp(pset,"beamhalo_minus")),
    md_()
  {
    md_.parameterSetID_ = pset.id();
    md_.moduleName_ = pset.getParameter<std::string>("@module_type");
    md_.moduleLabel_ = pset.getParameter<std::string>("@module_label");
    //#warning process name is hard coded, for now.  Fix this.
    //#warning the parameter set ID passed should be the one for the full process.  Fix this.
    md_.processConfiguration_ = ProcessConfiguration("PILEUP", pset.id(), getReleaseVersion(), getPassID());
  }

  // Virtual destructor needed.
  BMixingModule::~BMixingModule() { }  

  // Functions that get called by framework every event
  void BMixingModule::produce(edm::Event& e, const edm::EventSetup&) { 

    // Create EDProduct
    createnewEDProduct();
    // Add signals 
    addSignals(e);
    



    // Read the PileUp 
    std::vector<EventPrincipalVector> pileup[maxNbSources];
    bool doit[]={false,false,false,false};
    int bunchCrossing[maxNbSources];
    if ( input_)  {
      LogDebug("MixingModule") <<"\n\n==============================>Adding pileup to signal event "<<e.id(); 
      input_->readPileUp(pileup[0]); 
      bunchCrossing[0]=input_->minBunch();
      if ( input_->doPileup()) doit[0]=true;
    }
    if (cosmics_) {
      cosmics_->readPileUp(pileup[1]);
      bunchCrossing[1]=cosmics_->minBunch();
      LogDebug("MixingModule") <<"\n\n==============================>Adding cosmics to signal event "<<e.id(); 
      if (cosmics_->doPileup()) doit[1]=true;
    }
    if (beamHalo_p_) {
      beamHalo_p_->readPileUp(pileup[2]);
      bunchCrossing[2]=beamHalo_p_->minBunch();
      LogDebug("MixingModule") <<"\n\n==============================>Adding beam halo+ to signal event "<<e.id();
      if (beamHalo_p_->doPileup()) doit[2]=true;
    }
    if (beamHalo_m_) {
      beamHalo_m_->readPileUp(pileup[3]);
      bunchCrossing[3]=beamHalo_m_->minBunch();
      LogDebug("MixingModule") <<"\n\n==============================>Adding beam halo- to signal event "<<e.id();
      if (beamHalo_m_->doPileup()) doit[3]=true;
    }


    // and merge it
    // we have to loop over bunchcrossings first since added objects are all stored in one vector, 
    // ordered by bunchcrossing
      for (unsigned int isource=0;isource<maxNbSources;++isource) {
	if (doit[isource])   {
	  int TheBunch = bunchCrossing[isource];
	for (std::vector<EventPrincipalVector>::const_iterator it = pileup[isource].begin();
	     it != pileup[isource].end(); ++it, ++TheBunch) {
	  merge(TheBunch, *it);
	}}
      }
	    



    // Put output into event
    put(e);

  }

  void BMixingModule::merge(const int bcr, const EventPrincipalVector& vec) {
    //
    // main loop: loop over events and merge 
    //
    eventId_=0;
    LogDebug("merge") <<"For bunchcrossing "<<bcr<<", "<<vec.size()<< " events will be merged";
    vertexoffset=0;
    for (EventPrincipalVector::const_iterator it = vec.begin(); it != vec.end(); ++it) {
      Event e(**it, md_);
      LogDebug("merge") <<" merging Event:  id " << e.id();
      addPileups(bcr, &e, ++eventId_);
    }// end main loop
  }



} //edm
