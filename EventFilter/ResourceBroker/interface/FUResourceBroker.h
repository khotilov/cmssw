#ifndef FURESOURCEBROKER_H
#define FURESOURCEBROKER_H 1


#include "EventFilter/ResourceBroker/interface/FUTypes.h"
#include "EventFilter/ResourceBroker/interface/FUResourceTable.h"

#include "EventFilter/Utilities/interface/StateMachine.h"
#include "EventFilter/Utilities/interface/WebGUI.h"

#include "xdaq/include/xdaq/Application.h"
#include "xdaq/NamespaceURI.h"

#include "xdata/include/xdata/InfoSpace.h"
#include "xdata/include/xdata/String.h"
#include "xdata/include/xdata/Boolean.h"
#include "xdata/include/xdata/UnsignedInteger32.h"
#include "xdata/include/xdata/Double.h"

#include "toolbox/include/toolbox/task/TimerFactory.h"
#include "toolbox/include/Task.h"
#include "toolbox/include/toolbox/mem/Reference.h"
#include "toolbox/include/toolbox/fsm/exception/Exception.h"
#include "toolbox/include/BSem.h"

#include "interface/shared/include/frl_header.h"
#include "interface/shared/include/fed_header.h"
#include "interface/shared/include/fed_trailer.h"

#include <vector>
#include <string>
#include <semaphore.h>


namespace evf {

  class BUProxy;
  
  class FUResourceBroker : public xdaq::Application,
			   public toolbox::task::TimerListener,
			   public xdata::ActionListener
  {
  public:
    //
    // typedefs
    //
    typedef std::vector<BUProxy*> BUVec_t;
    
    
    //
    // xdaq instantiator macro
    //
    XDAQ_INSTANTIATOR();
    
    
    //
    // construction/destruction
    //
    FUResourceBroker(xdaq::ApplicationStub *s);
    virtual ~FUResourceBroker();
    
    
    //
    // public member functions
    //
    
    // work loop functions to be executed during transitional states (async)
    bool configuring(toolbox::task::WorkLoop* wl);
    bool enabling(toolbox::task::WorkLoop* wl);
    bool stopping(toolbox::task::WorkLoop* wl);
    bool halting(toolbox::task::WorkLoop* wl);
    
    // fsm soap command callback
    xoap::MessageReference fsmCallback(xoap::MessageReference msg)
      throw (xoap::exception::Exception);
    
    // toolbox::task::TimerListener callback, and init/start/stop the corresp. timer
    void timeExpired(toolbox::task::TimerEvent& e);
    void initTimer();
    void startTimer();
    void stopTimer();

    // xdata::ActionListener callback(s)
    void actionPerformed(xdata::Event& e);
    
    //  connection to builder unit bu_
    void connectToBUs();
    
    // Hyper DAQ web page(s) [see Utilities/WebGUI]
    void webPageRequest(xgi::Input *in,xgi::Output *out)
      throw (xgi::exception::Exception);
    
    // receive data from a builder unit
    void I2O_FU_TAKE_Callback(toolbox::mem::Reference *bufRef);
    
    
  private:
    //
    // private member functions
    //
    void exportParameters();
    void reset();
    
    
  private:    
    //
    // member data
    //
    
    // finite state machine
    evf::StateMachine        fsm_;
    
    // application identifier
    std::string              sourceId_;
    
    // binary semaphore
    BSem                     lock_;
    
    // Hyper DAQ web GUI
    WebGUI*                  gui_;
    
    // application logger
    Logger                   log_;
    
    // vector of connected builder units (BUs)
    BUVec_t                  bu_;
    
    // memory management for bu <-> fu comunication
    toolbox::mem::Pool*      i2oPool_;
    
    // manageme resources (events=fed collections) to be build
    FUResourceTable*         resourceTable_;
    
    // monitored parameters 
    xdata::String            url_;
    xdata::String            class_;
    xdata::UnsignedInteger32 instance_;
    xdata::UnsignedInteger32 runNumber_;
    xdata::UnsignedInteger32 nbShmClients_;
    
    xdata::Double            nbMBTot_;
    xdata::Double            nbMBPerSec_;
    xdata::Double            nbMBPerSecMin_;
    xdata::Double            nbMBPerSecMax_;
    xdata::Double            nbMBPerSecAvg_;
    
    // monitored counters
    xdata::UnsignedInteger32 nbEvents_;
    xdata::UnsignedInteger32 nbEventsPerSec_;
    xdata::UnsignedInteger32 nbEventsPerSecMin_;
    xdata::UnsignedInteger32 nbEventsPerSecMax_;
    xdata::UnsignedInteger32 nbEventsPerSecAvg_;
    xdata::UnsignedInteger32 nbAllocatedEvents_;
    xdata::UnsignedInteger32 nbPendingRequests_;
    xdata::UnsignedInteger32 nbReceivedEvents_;
    xdata::UnsignedInteger32 nbDiscardedEvents_;
    xdata::UnsignedInteger32 nbProcessedEvents_;
    xdata::UnsignedInteger32 nbLostEvents_;
    xdata::UnsignedInteger32 nbDataErrors_;
    xdata::UnsignedInteger32 nbCrcErrors_;

    // standard parameters
    xdata::Boolean           shmMode_;
    xdata::UnsignedInteger32 eventBufferSize_;
    //xdata::Boolean           doDumpFragments_;
    xdata::Boolean           doDropEvents_;
    xdata::Boolean           doFedIdCheck_;
    xdata::UnsignedInteger32 doCrcCheck_;

    xdata::String            buClassName_;
    xdata::UnsignedInteger32 buInstance_;
    xdata::UnsignedInteger32 queueSize_;
    
    // debug parameters
    xdata::UnsignedInteger32 nbAllocateSent_;
    xdata::UnsignedInteger32 nbTakeReceived_;
    
    // internal parameters, not exported
    unsigned int             nbMeasurements_;
    unsigned int             nbEventsLast_;
    
  };

} // namespace evf


#endif
