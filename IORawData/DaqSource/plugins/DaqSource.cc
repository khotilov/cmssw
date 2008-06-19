/** \file 
 *
 *  $Date: 2008/06/03 15:20:35 $
 *  $Revision: 1.22 $
 *  \author N. Amapane - S. Argiro'
 */

#include "DaqSource.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "EventFilter/Utilities/interface/GlobalEventNumber.h"

#include "IORawData/DaqSource/interface/DaqBaseReader.h"
#include "IORawData/DaqSource/interface/DaqReaderPluginFactory.h"

#include "DataFormats/Provenance/interface/Timestamp.h" 
#include "FWCore/Framework/interface/EventPrincipal.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Provenance/interface/LuminosityBlockAuxiliary.h"
#include "DataFormats/Provenance/interface/RunAuxiliary.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/LuminosityBlockPrincipal.h"
#include "FWCore/Framework/interface/RunPrincipal.h"

#include <string>
#include <iostream>
#include <sys/time.h>


////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////


namespace daqsource{
  static unsigned int gtpEvmId_ =  FEDNumbering::getTriggerGTPFEDIds().first;
}

namespace edm {
 namespace daqsource{
  static unsigned int gtpEvmId_ =  FEDNumbering::getTriggerGTPFEDIds().first;
 }

  //______________________________________________________________________________
  DaqSource::DaqSource(const ParameterSet& pset, 
		     const InputSourceDescription& desc) 
    : InputSource(pset,desc)
    , evf::ModuleWeb("DaqSource")
    , reader_(0)
    , lumiSegmentSizeInEvents_(pset.getUntrackedParameter<unsigned int>("evtsPerLS",0))
    , fakeLSid_(lumiSegmentSizeInEvents_ != 0)
    , runNumber_(RunID::firstValidRun().run())
    , luminosityBlockNumber_(LuminosityBlockID::firstValidLuminosityBlock().luminosityBlock())
    , noMoreEvents_(false)
    , newRun_(true)
    , newLumi_(true)
    , ep_()
    , lumiSectionIndex_(1)
    , prescaleSetIndex_(0)
    , is_(0)
    , mis_(0)
  {
    count = 0;
    pthread_mutex_init(&mutex_,0);
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
    if(is_)
      {
	is_->fireItemRevoked("lumiSectionIndex");
	is_->fireItemRevoked("prescaleSetIndex");
      }
    if(mis_)
      {
	mis_->fireItemRevoked("lumiSectionIndex");
	mis_->fireItemRevoked("prescaleSetIndex");
      }
    delete reader_;
    pthread_mutex_unlock(&mutex_);
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////
  // implementation of member functions
  ////////////////////////////////////////////////////////////////////////////////
  
  //______________________________________________________________________________
  InputSource::ItemType 
  DaqSource::getNextItemType() {
    if (noMoreEvents_) {
      return IsStop;
      pthread_mutex_unlock(&mutex_);
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
    timeval stv;
    gettimeofday(&stv,0);
    time = stv.tv_sec;
    time = (time << 32) + stv.tv_usec;
    Timestamp tstamp(time);

    int bunchCrossing = EventAuxiliary::invalidBunchXing;
    int orbitNumber   = EventAuxiliary::invalidBunchXing;
    
    // pass a 0 pointer to fillRawData()!
    FEDRawDataCollection* fedCollection(0);
  
    // let reader_ fill the fedCollection 
    if(!reader_->fillRawData(eventId, tstamp, fedCollection)) {
      // fillRawData() failed, clean up the fedCollection in case it was allocated!
      if (0 != fedCollection) delete fedCollection;
      noMoreEvents_ = true;
      pthread_mutex_unlock(&mutex_);
      return IsStop;
    }
    if (eventId.event() == 0) {
      throw cms::Exception("LogicError")
        << "The reader used with DaqSource has returned an invalid (zero) event number!\n"
        << "Event numbers must begin at 1, not 0.";
    }
    setTimestamp(tstamp);
    
    unsigned char *fedAddr = fedCollection->FEDData(daqsource::gtpEvmId_).data();

    if(fakeLSid_ && luminosityBlockNumber_ != ((eventId.event() - 1)/lumiSegmentSizeInEvents_ + 1)) {
	luminosityBlockNumber_ = (eventId.event() - 1)/lumiSegmentSizeInEvents_ + 1;
	pthread_mutex_unlock(&mutex_);
	pthread_mutex_lock(&mutex_);
        newLumi_ = true;
	lumiSectionIndex_.value_ = luminosityBlockNumber_;
	resetLuminosityBlockPrincipal();
    }
    else if(!fakeLSid_){ 

      if(evf::evtn::evm_board_sense(fedAddr)){
	unsigned int thisEventLSid = evf::evtn::getlbn(fedAddr);
       if(luminosityBlockNumber_ != (thisEventLSid + 1)){
         luminosityBlockNumber_ = thisEventLSid + 1;
	pthread_mutex_unlock(&mutex_);
	pthread_mutex_lock(&mutex_);
         newLumi_ = true;
	lumiSectionIndex_.value_ = luminosityBlockNumber_;
         resetLuminosityBlockPrincipal();
       }
      }
    }
    if(evf::evtn::evm_board_sense(fedAddr)){
      bunchCrossing =  int(evf::evtn::getfdlbx(fedAddr));
      orbitNumber =  int(evf::evtn::getorbit(fedAddr));
    }
    eventId = EventID(runNumber_, eventId.event());
    
    // If there is no luminosity block principal, make one.
    if (luminosityBlockPrincipal().get() == 0 || luminosityBlockPrincipal()->luminosityBlock() != luminosityBlockNumber_) {
      newLumi_ = true;
      LuminosityBlockAuxiliary lumiAux(runPrincipal()->run(),
        luminosityBlockNumber_, timestamp(), Timestamp::invalidTimestamp());
      setLuminosityBlockPrincipal(boost::shared_ptr<LuminosityBlockPrincipal>(
        new LuminosityBlockPrincipal(lumiAux,
                                     productRegistry(),
                                     runPrincipal(),
                                     processConfiguration())));

    }
    // make a brand new event
    EventAuxiliary eventAux(eventId,
      processGUID(), timestamp(), luminosityBlockPrincipal()->luminosityBlock(), true, EventAuxiliary::Data, bunchCrossing, EventAuxiliary::invalidStoreNumber,
			    orbitNumber);
    ep_ = std::auto_ptr<EventPrincipal>(
	new EventPrincipal(eventAux, productRegistry(), luminosityBlockPrincipal(), processConfiguration()));
    
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
    reset();
    newRun_ = newLumi_ = true;
    runNumber_ = r;
    noMoreEvents_ = false;
    resetLuminosityBlockPrincipal();
    resetRunPrincipal();
  }

  boost::shared_ptr<RunPrincipal>
  DaqSource::readRun_() {
    pthread_mutex_lock(&mutex_);
    assert(newRun_);
    assert(!noMoreEvents_);
    newRun_ = false;
    RunAuxiliary runAux(runNumber_, timestamp(), Timestamp::invalidTimestamp());
    return boost::shared_ptr<RunPrincipal>(
	new RunPrincipal(runAux,
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

  void DaqSource::publish(xdata::InfoSpace *is)
  {
    is_ = is;
    is->fireItemAvailable("lumiSectionIndex", &lumiSectionIndex_);
    is->fireItemAvailable("prescaleSetIndex", &prescaleSetIndex_);
  }
  void DaqSource::publishToXmas(xdata::InfoSpace *is)
  {
    mis_ = is;
    is->fireItemAvailable("lumiSectionIndex", &lumiSectionIndex_);
    is->fireItemAvailable("prescaleSetIndex", &prescaleSetIndex_);
  }

  void DaqSource::openBackDoor()
  {
    count++;
    if(count==2) throw;
    pthread_mutex_lock(&mutex_);
  }
  
  void DaqSource::closeBackDoor()
  {
    count--;
    pthread_mutex_unlock(&mutex_);
  }
}
