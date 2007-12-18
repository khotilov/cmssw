/** \file 
 *
 *  $Date: 2007/12/11 00:32:37 $
 *  $Revision: 1.13 $
 *  \author N. Amapane - S. Argiro'
 */

#include "DaqSource.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"

#include "IORawData/DaqSource/interface/DaqBaseReader.h"
#include "IORawData/DaqSource/interface/DaqReaderPluginFactory.h"

#include "DataFormats/Provenance/interface/Timestamp.h" 
#include "FWCore/Framework/interface/EventPrincipal.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/LuminosityBlockPrincipal.h"
#include "FWCore/Framework/interface/RunPrincipal.h"

#include <string>
#include <sys/time.h>


////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

namespace edm {

  //______________________________________________________________________________
  DaqSource::DaqSource(const ParameterSet& pset, 
		     const InputSourceDescription& desc) 
    : InputSource(pset,desc)
    , reader_(0)
    , lumiSegmentSizeInEvents_(pset.getUntrackedParameter<unsigned int>("evtsPerLS",0))
    , fakeLSid_(lumiSegmentSizeInEvents_ != 0)
    , runNumber_(RunID::firstValidRun().run())
    , luminosityBlockNumber_(LuminosityBlockID::firstValidLuminosityBlock().luminosityBlock())
    , noMoreEvents_(false)
    , newRun_(true)
    , newLumi_(true)
    , ep_() {
    produces<FEDRawDataCollection>();
    setTimestamp(Timestamp::beginOfTime());
    
    // Instantiate the requested data source
    std::string reader = pset.getUntrackedParameter<std::string>("readerPluginName");
    
    try{
      reader_=
        DaqReaderPluginFactory::get()->create(reader,
  					    pset.getUntrackedParameter<ParameterSet>("readerPset"));
    }
    catch(edm::Exception &e) {
      if(e.category() == "Configuration" && reader_ == 0) {
  	reader_ = DaqReaderPluginFactoryU::get()->create(reader);
  	if(reader_ == 0) throw;
      } else {
        throw;
      }
    }
  }
  
  //______________________________________________________________________________
  DaqSource::~DaqSource() {
    delete reader_;
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////
  // implementation of member functions
  ////////////////////////////////////////////////////////////////////////////////
  
  //______________________________________________________________________________
  InputSource::ItemType 
  DaqSource::getNextItemType() {
    if (noMoreEvents_) {
      return IsStop;
    }
    if (newRun_) {
      return IsRun;
    }
    if (newLumi_ && luminosityBlockPrincipal()) {
      return IsLumi;
    }
    if (ep_.get() != 0) {
      return IsEvent;
    }
    if(reader_ == 0) {
      throw cms::Exception("LogicError")
        << "DaqSource is used without a reader. Check your configuration !";
    }
    EventID eventId;
    TimeValue_t time = 0LL;
    gettimeofday((timeval *)(&time),0);
    Timestamp tstamp(time);
    
    // pass a 0 pointer to fillRawData()!
    FEDRawDataCollection* fedCollection(0);
  
    // let reader_ fill the fedCollection 
    if(!reader_->fillRawData(eventId, tstamp, fedCollection)) {
      // fillRawData() failed, clean up the fedCollection in case it was allocated!
      if (0 != fedCollection) delete fedCollection;
      noMoreEvents_ = true;
      return IsStop;
    }
    setTimestamp(tstamp);
    if(fakeLSid_ && luminosityBlockNumber_ != (eventId.event()/lumiSegmentSizeInEvents_ + 1)) {
	luminosityBlockNumber_ = eventId.event()/lumiSegmentSizeInEvents_ + 1;
        newLumi_ = true;
	resetLuminosityBlockPrincipal();
    }

    // Framework event numbers start at 1, not at zero.
    eventId = EventID(runNumber_, eventId.event() + 1);
    
    // If there is no luminosity block principal, make one.
    if (luminosityBlockPrincipal().get() == 0 || luminosityBlockPrincipal()->luminosityBlock() != luminosityBlockNumber_) {
      newLumi_ = true;
      setLuminosityBlockPrincipal(boost::shared_ptr<LuminosityBlockPrincipal>(
        new LuminosityBlockPrincipal(luminosityBlockNumber_,
                                     timestamp(),
                                     Timestamp::invalidTimestamp(),
                                     productRegistry(),
                                     runPrincipal(),
                                     processConfiguration())));

    }
    // make a brand new event
    ep_ = std::auto_ptr<EventPrincipal>(
	new EventPrincipal(eventId, timestamp(),
	productRegistry(), luminosityBlockPrincipal(), processConfiguration(), true, EventAuxiliary::Data));
    
    // have fedCollection managed by a std::auto_ptr<>
    std::auto_ptr<FEDRawDataCollection> bare_product(fedCollection);

    std::auto_ptr<Event> e(new Event(*ep_, moduleDescription()));
    // put the fed collection into the transient event store
    e->put(bare_product);
    // The commit is needed to complete the "put" transaction.
    e->commit_();
    if (newLumi_) {
      return IsLumi;
    }
    return IsEvent;
  }

  void
  DaqSource::setRun(RunNumber_t r) {
    assert(ep_.get() == 0);
    newRun_ = newLumi_ = true;
    runNumber_ = r;
    noMoreEvents_ = false;
    resetLuminosityBlockPrincipal();
    resetRunPrincipal();
  }

  boost::shared_ptr<RunPrincipal>
  DaqSource::readRun_() {
    assert(newRun_);
    assert(!noMoreEvents_);
    newRun_ = false;
    return boost::shared_ptr<RunPrincipal>(
	new RunPrincipal(runNumber_,
			 timestamp(),
			 Timestamp::invalidTimestamp(),
			 productRegistry(),
			 processConfiguration()));
  }

  boost::shared_ptr<LuminosityBlockPrincipal>
  DaqSource::readLuminosityBlock_() {
    assert(!newRun_);
    assert(newLumi_);
    assert(!noMoreEvents_);
    assert(luminosityBlockPrincipal());
    assert(ep_.get() != 0);
    newLumi_ = false;
    return luminosityBlockPrincipal();
  }

  std::auto_ptr<EventPrincipal>
  DaqSource::readEvent_(boost::shared_ptr<LuminosityBlockPrincipal> lbp) {
    assert(!newRun_);
    assert(!newLumi_);
    assert(!noMoreEvents_);
    assert(lbp);
    assert(ep_.get() != 0);
    std::auto_ptr<EventPrincipal> result = ep_;
    return result;
  }

  void
  DaqSource::setLumi(LuminosityBlockNumber_t) {
      throw cms::Exception("LogicError","DaqSource::setLumi(LuminosityBlockNumber_t lumiNumber)")
        << "The luminosity block number cannot be set externally for DaqSource.\n"
        << "Contact a Framework developer.\n";
  }

  std::auto_ptr<EventPrincipal>
  DaqSource::readIt(EventID const&) {
      throw cms::Exception("LogicError","DaqSource::readIt(EventID const& eventID)")
        << "Random access read cannot be used for DaqSource.\n"
        << "Contact a Framework developer.\n";
  }

  void
  DaqSource::skip(int) {
      throw cms::Exception("LogicError","DaqSource::skip(int offset)")
        << "Random access skip cannot be used for DaqSource\n"
        << "Contact a Framework developer.\n";
  }

}
