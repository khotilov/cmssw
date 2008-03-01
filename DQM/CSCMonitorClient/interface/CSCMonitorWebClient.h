#ifndef _CSCMonitorWebClient_h_
#define _CSCMonitorWebClient_h_

/**
  This class is an example DQM client with a web interface.
  Such clients inherit the state machine from the DQMBaseClient class
  of Components. From the UpdateObserver class of Components, they
  inherit the ability to define a function that gets automatically
  called when an update is received.
*/

#include "DQMServices/XdaqCollector/interface/DQMBaseClient.h"
#include "DQMServices/XdaqCollector/interface/Updater.h"
#include "DQMServices/XdaqCollector/interface/UpdateObserver.h"

#include "DQMServices/Core/interface/DQMOldReceiver.h"

#include "DQM/CSCMonitorClient/interface/CSCMonitorWebClientInterface.h"

#include <vector>
#include <string>
#include <iostream>


class CSCMonitorWebClient : public DQMBaseClient,
			       public dqm::UpdateObserver
{
public:

  /// You always need to have this line! Do not remove:
  XDAQ_INSTANTIATOR();

  /// The class constructor:
  CSCMonitorWebClient(xdaq::ApplicationStub *s);

  void Default(xgi::Input * in, xgi::Output * out ) throw (xgi::exception::Exception);

  /// implement the method that outputs the page with the widgets (declared in DQMBaseClient):
  void general(xgi::Input * in, xgi::Output * out ) throw (xgi::exception::Exception);

  /// the method which answers all HTTP requests of the form ".../Request?RequestID=..."
  void handleWebRequest(xgi::Input * in, xgi::Output * out);

  /// this obligatory method is called whenever the client enters the "Configured" state:
  void configure();

  /// this obligatory method is called whenever the client enters the "Enabled" state:
  void newRun();

  /// this obligatory method is called whenever the client enters the "Halted" state:
  void endRun();

  /// this obligatory method is called by the Updater component, whenever there is an update 
  void onUpdate() const;


public:

  /// this client has a web interface:  
  CSCMonitorWebClientInterface * webInterface_p;
};

/// You always need to have this line! Do not remove:
XDAQ_INSTANTIATOR_IMPL(CSCMonitorWebClient)

#endif
