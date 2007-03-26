////////////////////////////////////////////////////////////////////////////////
//
// FUEventProcessor
// ----------------
//
////////////////////////////////////////////////////////////////////////////////

#include "EventFilter/Processor/interface/FUEventProcessor.h"
#include "EventFilter/Utilities/interface/ModuleWebRegistry.h"
#include "EventFilter/Utilities/interface/Exception.h"
#include "EventFilter/Utilities/interface/ParameterSetRetriever.h"
#include "EventFilter/Utilities/interface/MicroStateService.h"
#include "EventFilter/Message2log4cplus/interface/MLlog4cplus.h"

#include "FWCore/Framework/interface/EventProcessor.h"
#include "FWCore/Framework/interface/RawInputSource.h"
#include "FWCore/PluginManager/interface/PresenceFactory.h"
#include "FWCore/PluginManager/interface/ProblemTracker.h"
#include "FWCore/ParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Presence.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Provenance/interface/ModuleDescription.h"

#include "DQMServices/Daemon/interface/MonitorDaemon.h"

#include "xcept/include/xcept/tools.h"
#include "xgi/include/xgi/Method.h"

#include "extern/cgicc/linuxx86/include/cgicc/CgiDefs.h"
#include "extern/cgicc/linuxx86/include/cgicc/Cgicc.h"
#include "extern/cgicc/linuxx86/include/cgicc/FormEntry.h"

#include "xoap/MessageReference.h"
#include "xoap/MessageFactory.h"
#include "xoap/include/xoap/SOAPEnvelope.h"
#include "xoap/include/xoap/SOAPBody.h"
#include "xoap/include/xoap/domutils.h"
#include "xoap/Method.h"

#include <typeinfo>
#include <stdlib.h>


namespace evf {

  namespace internal {
  
    void addService(vector<edm::ParameterSet>& adjust,string const& service)
    {
      edm::ParameterSet newpset;
      newpset.addParameter<string>("@service_type",service);
      adjust.push_back(newpset);
    }

    // Add a service to the services list if it is not already there
    void addServiceMaybe(vector<edm::ParameterSet>& adjust,string const& service)
    {
      std::vector<edm::ParameterSet>::const_iterator it;
      for(it=adjust.begin();it!=adjust.end();++it) {
	string name = it->getParameter<std::string>("@service_type");
	if (name == service) return;
      }
      addService(adjust, service);
    }
    
  } // namespace internal
  
} // namespace evf


using namespace std;
using namespace evf;
using namespace cgicc;


////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
FUEventProcessor::FUEventProcessor(xdaq::ApplicationStub *s) 
  : xdaq::Application(s)
  , fsm_(this)
  , evtProcessor_(0)
  , serviceToken_()
  , servicesDone_(false)
  , prescaleSvc_(0)
  , runNumber_(0)
  , epInitialized_(false)
  , outPut_(true)
  , inputPrescale_(1)
  , outputPrescale_(1)
  , timeoutOnStop_(10)
  , outprev_(true)
{
  // bind relevant callbacks to finite state machine
  fsm_.initialize<evf::FUEventProcessor>(this);
  
  //set sourceId_
  string       xmlClass = getApplicationDescriptor()->getClassName();
  unsigned int instance = getApplicationDescriptor()->getInstance();
  ostringstream sourcename;
  sourcename<<xmlClass<<"_"<<instance;
  LOG4CPLUS_INFO(getApplicationLogger(),sourcename.str()<<" constructor");
  cout<<"FUEventProcessor constructor"<<endl;
  LOG4CPLUS_INFO(getApplicationLogger(),"plugin path:"<<getenv("SEAL_PLUGINS"));
  
  ostringstream ns; ns<<"EP"<<instance;

  dqmCollectorAddr_       = "localhost";
  dqmCollectorPort_       = 9090;
  dqmCollectorDelay_      = 5000;
  dqmCollectorReconDelay_ = 5;
  dqmCollectorSourceName_ = ns.str();
  
  xdata::InfoSpace *ispace = getApplicationInfoSpace();

  // default configuration
  ispace->fireItemAvailable("parameterSet",         &configString_);
  ispace->fireItemAvailable("pluginPath",           &sealPluginPath_);
  ispace->fireItemAvailable("epInitialized",        &epInitialized_);
  ispace->fireItemAvailable("stateName",             fsm_.stateName());
  ispace->fireItemAvailable("runNumber",            &runNumber_);
  ispace->fireItemAvailable("outputEnabled",        &outPut_);
  ispace->fireItemAvailable("globalInputPrescale",  &inputPrescale_);
  ispace->fireItemAvailable("globalOutputPrescale", &outputPrescale_);
  ispace->fireItemAvailable("timeoutOnStop",        &timeoutOnStop_);
  
  ispace->fireItemAvailable("foundRcmsStateListener",fsm_.foundRcmsStateListener());
  
  ispace->fireItemAvailable("collectorAddr",        &dqmCollectorAddr_);
  ispace->fireItemAvailable("collectorPort",        &dqmCollectorPort_);
  ispace->fireItemAvailable("collSendUs",           &dqmCollectorDelay_);
  ispace->fireItemAvailable("collReconnSec",        &dqmCollectorReconDelay_);
  ispace->fireItemAvailable("monSourceName",        &dqmCollectorSourceName_);
  
  ispace->fireItemAvailable("prescalerAsString",    &prescalerAsString_);
  ispace->fireItemAvailable("triggerReportAsString",&triggerReportAsString_);
  
  // Add infospace listeners for exporting data values
  getApplicationInfoSpace()->addItemChangedListener("parameterSet",        this);
  getApplicationInfoSpace()->addItemChangedListener("outputEnabled",       this);
  getApplicationInfoSpace()->addItemChangedListener("globalInputPrescale", this);
  getApplicationInfoSpace()->addItemChangedListener("globalOutputPrescale",this);
 
  // bind prescale related soap callbacks
  xoap::bind(this,&FUEventProcessor::getPsReport ,"GetPsReport",XDAQ_NS_URI);
  xoap::bind(this,&FUEventProcessor::setPsUpdate ,"SetPsUpdate",XDAQ_NS_URI);
  xoap::bind(this,&FUEventProcessor::putPrescaler,"PutPrescaler",XDAQ_NS_URI);
  
  // Bind web interface
  xgi::bind(this, &FUEventProcessor::css           , "styles.css");
  xgi::bind(this, &FUEventProcessor::defaultWebPage, "Default"   );
  xgi::bind(this, &FUEventProcessor::moduleWeb     , "moduleWeb" );
  xgi::bind(this, &FUEventProcessor::microState    , "microState");

  // instantiate the plugin manager, not referenced here after!
  edm::AssertHandler ah;
}


//______________________________________________________________________________
FUEventProcessor::~FUEventProcessor()
{
  if (evtProcessor_) delete evtProcessor_;
}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void FUEventProcessor::getTriggerReport(toolbox::Event::Reference e)
  throw (toolbox::fsm::exception::Exception)
{
  // Calling this method results in calling 
  // evtProcessor_->getTriggerReport, the value returned is encoded as
  // a string. This value is used to set the exported SOAP param :
  // 'triggerReportAsString_'. The FM then picks up this value use getParam...
  LOG4CPLUS_DEBUG(getApplicationLogger(),"getTriggerReport action invoked");
  
  //Get the trigger report.
  edm::TriggerReport tr; 
  evtProcessor_->getTriggerReport(tr);
  
  triggerReportAsString_ = triggerReportToString(tr);
  
  //Print the trigger report message in debug format.
  printTriggerReport(tr);
}


//______________________________________________________________________________
xoap::MessageReference FUEventProcessor::getPsReport(xoap::MessageReference msg)
  throw (xoap::exception::Exception)
{
  // callback to return the trigger statistics as a string
  // cout <<"getPsReport from cout " <<endl;
  LOG4CPLUS_DEBUG(getApplicationLogger(),"getPsReport from log4");

  // print request
  //msg->writeTo(std::cout);
  
  //Get the trigger report.
  edm::TriggerReport tr; 
  evtProcessor_->getTriggerReport(tr);
  
  // xdata::String ReportAsString = triggerReportToString(tr);
  string s = triggerReportToString(tr);
  
  // reply message
  try {
      xoap::MessageReference reply = xoap::createMessage();
      xoap::SOAPEnvelope envelope = reply->getSOAPPart().getEnvelope();
      xoap::SOAPBody body = envelope.getBody();
      xoap::SOAPName responseName = envelope.createName("getPsReportResponse", "xdaq", "XDAQ_NS_URI");
      xoap::SOAPBodyElement responseElement = body.addBodyElement(responseName);
      xoap::SOAPName attributeName = envelope.createName("state", "xdaq", "XDAQ_NS_URI");
      xoap::SOAPElement keyElement = responseElement.addChildElement(attributeName);
      keyElement.addTextNode(s);
      return reply;

  }
  catch (xcept::Exception &e) {
    XCEPT_RETHROW(xoap::exception::Exception,
		  "Failed to create getPsReport response message",e);
  }
}


//______________________________________________________________________________
xoap::MessageReference FUEventProcessor::setPsUpdate(xoap::MessageReference msg)
  throw (xoap::exception::Exception)
{
  // callback to return the trigger statistics as a string
  // LOG4CPLUS_DEBUG(this->getApplicationLogger(), "setPsUpdate from log4");
  // cout <<"setPsUpdate from cout " <<endl;
  // msg->writeTo(std::cout);
  
  // decode 
  xoap::SOAPPart part = msg->getSOAPPart();
  xoap::SOAPEnvelope env = part.getEnvelope();
  xoap::SOAPBody msgbody = env.getBody();
  DOMNode* node = msgbody.getDOMNode();
  
  string requestString;
  DOMNodeList* bodyList = node->getChildNodes();
  for (unsigned int i = 0; i < bodyList->getLength(); i++) {
    DOMNode* command = bodyList->item(i);
    if (command->getNodeType() == DOMNode::ELEMENT_NODE) {
      std::string commandName = xoap::XMLCh2String (command->getLocalName());
      if ( commandName == "state" ) {
	if ( command->hasAttributes() ) {
	  DOMNamedNodeMap * map = command->getAttributes();
	  for (int l=0 ; l< (int)map->getLength() ; l++) {
	    // loop over attributes of node
	    DOMNode * anode = map->item(l);
	    string attributeName = XMLString::transcode(anode->getNodeName());
	    if (attributeName == "xdaq:stateName")
	      requestString = xoap::XMLCh2String(anode->getNodeValue());
	  }
	}
      }
    }
  }
  
  // reply message
  try {
      xoap::MessageReference reply = xoap::createMessage();
      xoap::SOAPEnvelope envelope = reply->getSOAPPart().getEnvelope();
      xoap::SOAPBody body = envelope.getBody();
      xoap::SOAPName responseName = envelope.createName("getPsReportResponse", "xdaq", "XDAQ_NS_URI");
      xoap::SOAPBodyElement responseElement = body.addBodyElement(responseName);
      xoap::SOAPName attributeName = envelope.createName("state", "xdaq", "XDAQ_NS_URI");
      xoap::SOAPElement keyElement = responseElement.addChildElement(attributeName);
      keyElement.addTextNode(requestString);
      return reply;
    
  }
  catch (xcept::Exception &e) {
    XCEPT_RETHROW(xoap::exception::Exception,
		  "Failed to create setPsUpdate response message", e);
  }
}


//______________________________________________________________________________
xoap::MessageReference FUEventProcessor::putPrescaler(xoap::MessageReference msg)
  throw (xoap::exception::Exception)
{
  //The EPSM has an exported SOAP param 'nextPrescalerAsString_' this is always
  //set first from the FM with the new prescaler value encoded as a string.
  //Next this function is called to pick up the new string value and fill the 
  //appropriate prescaler structure for addition to the prescaler cache...
  
  LOG4CPLUS_INFO(getApplicationLogger(),"putPrescaler action invoked");
  
  //  msg->writeTo(std::cout);
  //  cout << endl;
  
  string prescalerAsString = "INITIAL_VALUE";
  string replyPs = "";
  
  // decode 
  xoap::SOAPPart part = msg->getSOAPPart();
  xoap::SOAPEnvelope env = part.getEnvelope();
  xoap::SOAPBody msgbody = env.getBody();
  DOMNode* node = msgbody.getDOMNode();
  
  DOMNodeList* bodyList = node->getChildNodes();
  for (unsigned int i = 0; i < bodyList->getLength(); i++) {
    DOMNode* command = bodyList->item(i);
    if (command->getNodeType() == DOMNode::ELEMENT_NODE) {
      std::string commandName = xoap::XMLCh2String (command->getLocalName());
      if ( commandName == "PutPrescaler" ) {
	if ( command->hasAttributes() ) {
	  DOMNamedNodeMap * map = command->getAttributes();
	  for (int l=0 ; l< (int)map->getLength() ; l++) {
	    // loop over attributes of node
	    DOMNode * anode = map->item(l);
	    string attributeName = XMLString::transcode(anode->getNodeName());
	    if (attributeName == "prescalerAsString")
	      prescalerAsString =  xoap::XMLCh2String(anode->getNodeValue());
	  }
	}
      }
    }
  }
  
  //Get the prescaler string value. (Which was set by the FM)
  LOG4CPLUS_INFO(getApplicationLogger(),
		 "Using new prescaler string setting: "<<prescalerAsString);


  if ( prescalerAsString == "INITIAL_VALUE" ) {
    // cout << "prescalerAsString not updated, is " << prescalerAsString << endl;
  }
  else {
    if(prescaleSvc_ != 0) {
      //The prescale value associated with the LS# and module name.
      int prescaleValue = prescaleSvc_->getPrescale(1,"prescale2");
      LOG4CPLUS_DEBUG(getApplicationLogger(),
		      "prescaleSvc_->getPrescale(1,prescale2): "<<prescaleValue);
      
      //The number of LS# to prescale module set associations in the prescale
      //cache.
      int storeSize = prescaleSvc_->putPrescale(prescalerAsString);
      LOG4CPLUS_DEBUG(getApplicationLogger(),
		      "prescaleSvc_->putPrescale(s): " << storeSize);
      replyPs = prescaleSvc_->getTriggerCounters();
    }
    else {
      LOG4CPLUS_DEBUG(getApplicationLogger(),"PrescaleService pointer == 0"); 
    }
  }

    xoap::MessageReference reply = xoap::createMessage();
    xoap::SOAPEnvelope envelope = reply->getSOAPPart().getEnvelope();
    xoap::SOAPBody body = envelope.getBody();
    xoap::SOAPName responseName = envelope.createName("PutPrescalerResponse", "xdaq", "XDAQ_NS_URI");
    xoap::SOAPBodyElement responseElement = body.addBodyElement(responseName);
    xoap::SOAPName attributeName = envelope.createName("prescalerAsString", "xdaq", "XDAQ_NS_URI");
    xoap::SOAPElement keyElement = responseElement.addChildElement(attributeName);
    keyElement.addTextNode(replyPs);

  return reply;
}


//______________________________________________________________________________
bool FUEventProcessor::configuring(toolbox::task::WorkLoop* wl)
{
  try {
    LOG4CPLUS_INFO(getApplicationLogger(),"Start configuring ...");
    initEventProcessor();
    LOG4CPLUS_INFO(getApplicationLogger(),"Finished configuring!");
    
    fsm_.fireEvent("ConfigureDone",this);
  }
  catch (xcept::Exception &e) {
    string msg = "configuring FAILED: " + (string)e.what();
    fsm_.fireFailed(msg,this);
  }

  return false;
}


//______________________________________________________________________________
bool FUEventProcessor::enabling(toolbox::task::WorkLoop* wl)
{
  try {
    LOG4CPLUS_INFO(getApplicationLogger(),"Start enabling ...");
    
    // if the ep is intialized already, the initialization will be skipped
    initEventProcessor();
    
    int sc = 0;
    evtProcessor_->setRunNumber(runNumber_.value_);
    try {
      evtProcessor_->runAsync();
      sc = evtProcessor_->statusAsync();
    }
    catch(seal::Error& e) {
      fsm_.fireFailed(e.explainSelf(),this);
      return false;
    }
    catch(cms::Exception &e) {
      fsm_.fireFailed(e.explainSelf(),this);
      return false;
    }    
    catch(std::exception &e) {
      fsm_.fireFailed(e.what(),this);
      return false;
    }
    catch(...) {
      fsm_.fireFailed("Unknown Exception",this);
      return false;
    }
    
    if(sc != 0) {
      ostringstream oss;
      oss<<"EventProcessor::runAsync returned status code " << sc;
      fsm_.fireFailed(oss.str(),this);
      return false;
    }
    
    LOG4CPLUS_INFO(getApplicationLogger(),"Finished enabling!");
    
    fsm_.fireEvent("EnableDone",this);
  }
  catch (xcept::Exception &e) {
    string msg = "enabling FAILED: " + (string)e.what();
    fsm_.fireFailed(msg,this);
  }
  
  return false;
}


//______________________________________________________________________________
bool FUEventProcessor::stopping(toolbox::task::WorkLoop* wl)
{
  try {
    LOG4CPLUS_INFO(getApplicationLogger(),"Start stopping :) ...");
    edm::EventProcessor::StatusCode rc = stopEventProcessor();
    LOG4CPLUS_INFO(getApplicationLogger(),"Finished stopping!");
    if(rc != edm::EventProcessor::epTimedOut) 
      fsm_.fireEvent("StopDone",this);
    else
      fsm_.fireFailed("EventProcessor stop timed out",this);
  }
  catch (xcept::Exception &e) {
    string msg = "stopping FAILED: " + (string)e.what();
    fsm_.fireFailed(msg,this);
  }
  
  return false;
}


//______________________________________________________________________________
bool FUEventProcessor::halting(toolbox::task::WorkLoop* wl)
{
  try {
    LOG4CPLUS_INFO(getApplicationLogger(),"Start halting ...");
    edm::EventProcessor::StatusCode rc = stopEventProcessor();
    if(rc != edm::EventProcessor::epTimedOut)
      {
	evtProcessor_->endJob();
	delete evtProcessor_;
	evtProcessor_ = 0;
	epInitialized_ = false;
	LOG4CPLUS_INFO(getApplicationLogger(),"Finished halting!");
  
	fsm_.fireEvent("HaltDone",this);
      }
    else
      fsm_.fireFailed("EventProcessor stop timed out",this);
  }
  catch (xcept::Exception &e) {
    string msg = "halting FAILED: " + (string)e.what();
    fsm_.fireFailed(msg,this);
  }
  
  return false;
}


//______________________________________________________________________________
xoap::MessageReference FUEventProcessor::fsmCallback(xoap::MessageReference msg)
  throw (xoap::exception::Exception)
{
  return fsm_.commandCallback(msg);
}


//______________________________________________________________________________
void FUEventProcessor::initEventProcessor()
{
  if (epInitialized_) {
    LOG4CPLUS_INFO(getApplicationLogger(),
		   "CMSSW EventProcessor already initialized: skip!");
    return;
  }
  
  LOG4CPLUS_INFO(getApplicationLogger(),"Initialize CMSSW EventProcessor.");
  
  if (0!=setenv("SEAL_PLUGINS",sealPluginPath_.value_.c_str(),0)) {
    LOG4CPLUS_ERROR(getApplicationLogger(),"Failed to set SEAL_PLUGINS search path");
  }
  else {
    LOG4CPLUS_INFO(getApplicationLogger(),"plugin path: "<<getenv("SEAL_PLUGINS"));
  }
  
  
  // job configuration string
  ParameterSetRetriever pr(configString_.value_);
  string configString = pr.getAsString();
  
  
  boost::shared_ptr<edm::ParameterSet> params; // change this name!
  boost::shared_ptr<vector<edm::ParameterSet> > pServiceSets;
  makeParameterSets(configString, params, pServiceSets);
  
  // add default set of services
  if(!servicesDone_) {
    internal::addServiceMaybe(*pServiceSets,"DaqMonitorROOTBackEnd");
    internal::addServiceMaybe(*pServiceSets,"MonitorDaemon");
    internal::addServiceMaybe(*pServiceSets,"MLlog4cplus");
    internal::addServiceMaybe(*pServiceSets,"MicroStateService");
    internal::addServiceMaybe(*pServiceSets,"PrescaleService");
    
    try{
      serviceToken_ = edm::ServiceRegistry::createSet(*pServiceSets);
    }
    catch(seal::Error& e) {
      LOG4CPLUS_ERROR(getApplicationLogger(),e.explainSelf());
    }
    catch(cms::Exception &e) {
      LOG4CPLUS_ERROR(getApplicationLogger(),e.explainSelf());
    }    
    catch(std::exception &e) {
      LOG4CPLUS_ERROR(getApplicationLogger(),e.what());
    }
    catch(...) {
      LOG4CPLUS_ERROR(getApplicationLogger(),"Unknown Exception");
    }
  }

  
  edm::ServiceRegistry::Operate operate(serviceToken_);
  try{
    edm::Service<MonitorDaemon>()->rmt(dqmCollectorAddr_,
				       dqmCollectorPort_,
				       dqmCollectorDelay_,
				       dqmCollectorSourceName_,
				       dqmCollectorReconDelay_);
    edm::Service<ML::MLlog4cplus>()->setAppl(this);
  }
  catch(...) { 
    LOG4CPLUS_INFO(getApplicationLogger(),
		   "exception when trying to get service MonitorDaemon");
  }

  
  if(!servicesDone_) {
    try{
      LOG4CPLUS_DEBUG(getApplicationLogger(),
		      "Trying to create message service presence ");
      edm::PresenceFactory *pf = edm::PresenceFactory::get();
      if(pf != 0) {
	pf->makePresence("MessageServicePresence").release();
      }
      else {
	LOG4CPLUS_ERROR(getApplicationLogger(),
			"Unable to create message service presence ");
      }
      
      servicesDone_ = true;
      
    } 
    catch(seal::Error& e) {
      LOG4CPLUS_ERROR(getApplicationLogger(),e.explainSelf());
    }
    catch(cms::Exception &e) {
      LOG4CPLUS_ERROR(getApplicationLogger(),e.explainSelf());
    }    
    catch(std::exception &e) {
      LOG4CPLUS_ERROR(getApplicationLogger(),e.what());
    }
    catch(...) {
      LOG4CPLUS_ERROR(getApplicationLogger(),"Unknown Exception");
    }
    
  }
  
  //test rerouting of fwk logging to log4cplus
  edm::LogInfo("FUEventProcessor")<<"started MessageLogger Service.";
  edm::LogInfo("FUEventProcessor")<<"Using config string \n"<<configString;


  // instantiate the event processor
  try{
    vector<string> defaultServices;
    defaultServices.push_back("MessageLogger");
    defaultServices.push_back("InitRootHandlers");
    defaultServices.push_back("JobReportService");
    
    if (0!=evtProcessor_) delete evtProcessor_;
    
    evtProcessor_ = new edm::EventProcessor(configString,
					    serviceToken_,
					    edm::serviceregistry::kTokenOverrides,
					    defaultServices);
    //    evtProcessor_->setRunNumber(runNumber_.value_);

    if(!outPut_)
      //evtProcessor_->toggleOutput();
      //evtProcessor_->prescaleInput(inputPrescale_);
      //evtProcessor_->prescaleOutput(outputPrescale_);
      evtProcessor_->enableEndPaths(outPut_);
    
    outprev_=outPut_;
    
    // to publish all module names to XDAQ infospace
    ModuleWebRegistry *mwr = 0;
    try{
      if(edm::Service<ModuleWebRegistry>().isAvailable())
	mwr = edm::Service<ModuleWebRegistry>().operator->();
    }
    catch(...) { 
      LOG4CPLUS_INFO(getApplicationLogger(),
		     "exception when trying to get service ModuleWebRegistry");
    }

    if(mwr) mwr->publish(getApplicationInfoSpace());

    // get the prescale service
    LOG4CPLUS_INFO(getApplicationLogger(),
		   "Checking for edm::service::PrescaleService!");
    try {
      if(edm::Service<edm::service::PrescaleService>().isAvailable())
	{
	  LOG4CPLUS_INFO(getApplicationLogger(),
			 "edm::service::PrescaleService is available!");
	  prescaleSvc_ = edm::Service<edm::service::PrescaleService>().operator->();
	  LOG4CPLUS_INFO(getApplicationLogger(),
			 "Obtained pointer to PrescaleService");
	  prescaleSvc_->putHandle(evtProcessor_);
	  LOG4CPLUS_INFO(getApplicationLogger(),
			 "PrescaleService::putHandle called");
	}
    }
    catch(...) {
      LOG4CPLUS_INFO(getApplicationLogger(),
		     "exception when trying to get service "
		     <<"edm::service::PrescaleService");
    }
  }
  catch(seal::Error& e) {
    fsm_.fireFailed(e.explainSelf(),this);
    return;
  }
  catch(cms::Exception &e) {
    fsm_.fireFailed(e.explainSelf(),this);
    return;
  }    
  catch(std::exception &e) {
    fsm_.fireEvent(e.what(),this);
    return;
  }
  catch(...) {
    fsm_.fireFailed("Unknown Exception",this);
    return;
  }
  
  LOG4CPLUS_INFO(getApplicationLogger(),"FUEventProcessor configuration finished.");
  
  epInitialized_ = true;

  return;
}


//______________________________________________________________________________
edm::EventProcessor::StatusCode FUEventProcessor::stopEventProcessor()
{
  edm::EventProcessor::StatusCode rc = edm::EventProcessor::epSuccess;
  try  {
    rc = evtProcessor_->waitTillDoneAsync(timeoutOnStop_.value_);
  }
  catch(seal::Error& e) {
    XCEPT_RAISE(evf::Exception,e.explainSelf());
  }
  catch(cms::Exception &e) {
    XCEPT_RAISE(evf::Exception,e.explainSelf());
  }    
  catch(std::exception &e) {
    XCEPT_RAISE(evf::Exception,e.what());
  }
  catch(...) {
    XCEPT_RAISE(evf::Exception,"Unknown Exception");
  }
  return rc;

}


//______________________________________________________________________________
void FUEventProcessor::actionPerformed(xdata::Event& e)
{
  if (e.type()=="ItemChangedEvent" && !(fsm_.stateName()->toString()=="Halted")) {
    string item = dynamic_cast<xdata::ItemChangedEvent&>(e).itemName();
    
    if ( item == "parameterSet") {
      epInitialized_ = false;
    }
    
    if ( item == "outputEnabled") {
      if(outprev_ != outPut_) {
	LOG4CPLUS_WARN(getApplicationLogger(),
		       (outprev_ ? "Disabling " : "Enabling ")<<"global output");
	evtProcessor_->enableEndPaths(outPut_);
	outprev_ = outPut_;
      }
    }
    
    if (item == "globalInputPrescale") {
      //evtProcessor_->prescaleInput(inputPrescale_);
      //LOG4CPLUS_WARN(this->getApplicationLogger(),
      //"Setting global input prescale factor to" << inputPrescale_);
      //
      LOG4CPLUS_WARN(getApplicationLogger(),
		     "Setting global input prescale has no effect "
		     <<"in this version of the code");
    }
    if ( item == "globalOutputPrescale") {
      //evtProcessor_->prescaleOutput(outputPrescale_);
      //LOG4CPLUS_WARN(this->getApplicationLogger(),
      //"Setting global output prescale factor to" << outputPrescale_);
      LOG4CPLUS_WARN(getApplicationLogger(),
		     "Setting global output prescale has no effect "
		     <<"in this version of the code");
    }
  }

}


//______________________________________________________________________________
string FUEventProcessor::triggerReportToString(const edm::TriggerReport& tr)
{
  // Add an array length indicator so that the resulting string will have a 
  // little more readability.
  string ARRAY_LEN = "_";
  string SEPARATOR = " ";
  
  ostringstream oss;
  
  //TriggerReport::eventSummary
  oss<<tr.eventSummary.totalEvents<<SEPARATOR 
     <<tr.eventSummary.totalEventsPassed<<SEPARATOR
     <<tr.eventSummary.totalEventsFailed<<SEPARATOR;
  
  //TriggerReport::trigPathSummaries
  oss<<ARRAY_LEN<<tr.trigPathSummaries.size()<<SEPARATOR;
  for(unsigned int i=0; i<tr.trigPathSummaries.size(); i++) {
    oss<<tr.trigPathSummaries[i].bitPosition<<SEPARATOR 
       <<tr.trigPathSummaries[i].timesRun<<SEPARATOR
       <<tr.trigPathSummaries[i].timesPassed<<SEPARATOR
       <<tr.trigPathSummaries[i].timesFailed<<SEPARATOR
       <<tr.trigPathSummaries[i].timesExcept<<SEPARATOR
       <<tr.trigPathSummaries[i].name<<SEPARATOR;
    
    //TriggerReport::trigPathSummaries::moduleInPathSummaries
    oss<<ARRAY_LEN<<tr.trigPathSummaries[i].moduleInPathSummaries.size()<<SEPARATOR;
    for(unsigned int j=0;j<tr.trigPathSummaries[i].moduleInPathSummaries.size();j++) {
      oss<<tr.trigPathSummaries[i].moduleInPathSummaries[j].timesVisited<<SEPARATOR
	 <<tr.trigPathSummaries[i].moduleInPathSummaries[j].timesPassed <<SEPARATOR
	 <<tr.trigPathSummaries[i].moduleInPathSummaries[j].timesFailed <<SEPARATOR
	 <<tr.trigPathSummaries[i].moduleInPathSummaries[j].timesExcept <<SEPARATOR
	 <<tr.trigPathSummaries[i].moduleInPathSummaries[j].moduleLabel <<SEPARATOR;
    }
  }
  
  //TriggerReport::endPathSummaries
  oss<<ARRAY_LEN<<tr.endPathSummaries.size()<<SEPARATOR;
  for(unsigned int i=0; i<tr.endPathSummaries.size(); i++) {
    oss<<tr.endPathSummaries[i].bitPosition<<SEPARATOR 
       <<tr.endPathSummaries[i].timesRun<<SEPARATOR
       <<tr.endPathSummaries[i].timesPassed<<SEPARATOR
       <<tr.endPathSummaries[i].timesFailed<<SEPARATOR
       <<tr.endPathSummaries[i].timesExcept<<SEPARATOR
       <<tr.endPathSummaries[i].name<<SEPARATOR;
    
    //TriggerReport::endPathSummaries::moduleInPathSummaries
    oss<<ARRAY_LEN<<tr.endPathSummaries[i].moduleInPathSummaries.size()<<SEPARATOR;
    for(unsigned int j=0;j<tr.endPathSummaries[i].moduleInPathSummaries.size();j++) {
      oss<<tr.endPathSummaries[i].moduleInPathSummaries[j].timesVisited<<SEPARATOR
	 <<tr.endPathSummaries[i].moduleInPathSummaries[j].timesPassed <<SEPARATOR
	 <<tr.endPathSummaries[i].moduleInPathSummaries[j].timesFailed <<SEPARATOR
	 <<tr.endPathSummaries[i].moduleInPathSummaries[j].timesExcept <<SEPARATOR
	 <<tr.endPathSummaries[i].moduleInPathSummaries[j].moduleLabel <<SEPARATOR;
    }
  }
  
  //TriggerReport::workerSummaries
  oss<<ARRAY_LEN<<tr.workerSummaries.size()<<SEPARATOR;
  for(unsigned int i=0; i<tr.workerSummaries.size(); i++) {
    oss<<tr.workerSummaries[i].timesVisited<<SEPARATOR 
       <<tr.workerSummaries[i].timesRun    <<SEPARATOR
       <<tr.workerSummaries[i].timesPassed <<SEPARATOR
       <<tr.workerSummaries[i].timesFailed <<SEPARATOR
       <<tr.workerSummaries[i].timesExcept <<SEPARATOR
       <<tr.workerSummaries[i].moduleLabel <<SEPARATOR;
  }
  
  return oss.str();
}


//______________________________________________________________________________
void FUEventProcessor::printTriggerReport(const edm::TriggerReport& tr)
{
  ostringstream oss;
  
  oss<<"================================="<<"\n";
  oss<<"== BEGINNING OF TRIGGER REPORT =="<<"\n";
  oss<<"================================="<<"\n";
  oss<<"tr.eventSummary.totalEvents= "<<tr.eventSummary.totalEvents<<"\n" 
     <<"tr.eventSummary.totalEventsPassed= "<<tr.eventSummary.totalEventsPassed<<"\n"
     <<"tr.eventSummary.totalEventsFailed= "<<tr.eventSummary.totalEventsFailed<<"\n";
  
  oss<<"TriggerReport::trigPathSummaries"<<"\n";
  for(unsigned int i=0; i<tr.trigPathSummaries.size(); i++) {
    oss<<"tr.trigPathSummaries["<<i<<"].bitPosition = "
       <<tr.trigPathSummaries[i].bitPosition <<"\n" 
       <<"tr.trigPathSummaries["<<i<<"].timesRun = "
       <<tr.trigPathSummaries[i].timesRun <<"\n"
       <<"tr.trigPathSummaries["<<i<<"].timesPassed = "
       <<tr.trigPathSummaries[i].timesPassed <<"\n"
       <<"tr.trigPathSummaries["<<i<<"].timesFailed = "
       <<tr.trigPathSummaries[i].timesFailed <<"\n"
       <<"tr.trigPathSummaries["<<i<<"].timesExcept = "
       <<tr.trigPathSummaries[i].timesExcept <<"\n"
       <<"tr.trigPathSummaries["<<i<<"].name = "
       <<tr.trigPathSummaries[i].name <<"\n";
    
    //TriggerReport::trigPathSummaries::moduleInPathSummaries
    for(unsigned int j=0;j<tr.trigPathSummaries[i].moduleInPathSummaries.size();j++) {
      oss<<"tr.trigPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].timesVisited = "
	 <<tr.trigPathSummaries[i].moduleInPathSummaries[j].timesVisited<<"\n"
	 <<"tr.trigPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].timesPassed = "
	 <<tr.trigPathSummaries[i].moduleInPathSummaries[j].timesPassed<<"\n"
	 <<"tr.trigPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].timesFailed = "
	 <<tr.trigPathSummaries[i].moduleInPathSummaries[j].timesFailed<<"\n"
	 <<"tr.trigPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].timesExcept = "
	 <<tr.trigPathSummaries[i].moduleInPathSummaries[j].timesExcept<<"\n"
	 <<"tr.trigPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].moduleLabel = "
	 <<tr.trigPathSummaries[i].moduleInPathSummaries[j].moduleLabel<<"\n";
    }
  }
  
  //TriggerReport::endPathSummaries
  for(unsigned int i=0;i<tr.endPathSummaries.size();i++) {
    oss<<"tr.endPathSummaries["<<i<<"].bitPosition = "
       <<tr.endPathSummaries[i].bitPosition <<"\n" 
       <<"tr.endPathSummaries["<<i<<"].timesRun = "
       <<tr.endPathSummaries[i].timesRun <<"\n"
       <<"tr.endPathSummaries["<<i<<"].timesPassed = "
       <<tr.endPathSummaries[i].timesPassed <<"\n"
       <<"tr.endPathSummaries["<<i<<"].timesFailed = "
       <<tr.endPathSummaries[i].timesFailed <<"\n"
       <<"tr.endPathSummaries["<<i<<"].timesExcept = "
       <<tr.endPathSummaries[i].timesExcept <<"\n"
       <<"tr.endPathSummaries["<<i<<"].name = "
       <<tr.endPathSummaries[i].name <<"\n";
    
    //TriggerReport::endPathSummaries::moduleInPathSummaries
    for(unsigned int j=0;j<tr.endPathSummaries[i].moduleInPathSummaries.size();j++) {
      oss<<"tr.endPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].timesVisited = "
	 <<tr.endPathSummaries[i].moduleInPathSummaries[j].timesVisited <<"\n"
	 <<"tr.endPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].timesPassed = "
	 <<tr.endPathSummaries[i].moduleInPathSummaries[j].timesPassed <<"\n"
	 <<"tr.endPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].timesFailed = "
	 <<tr.endPathSummaries[i].moduleInPathSummaries[j].timesFailed <<"\n"
	 <<"tr.endPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].timesExcept = "
	 <<tr.endPathSummaries[i].moduleInPathSummaries[j].timesExcept <<"\n"
	 <<"tr.endPathSummaries["<<i
	 <<"].moduleInPathSummaries["<<j<<"].moduleLabel = "
	 <<tr.endPathSummaries[i].moduleInPathSummaries[j].moduleLabel <<"\n";
    }
  }
  
  //TriggerReport::workerSummaries
  for(unsigned int i=0; i<tr.workerSummaries.size(); i++) {
    oss<<"tr.workerSummaries["<<i<<"].timesVisited = "
       <<tr.workerSummaries[i].timesVisited<<"\n" 
       <<"tr.workerSummaries["<<i<<"].timesRun = "
       <<tr.workerSummaries[i].timesRun<<"\n"
       <<"tr.workerSummaries["<<i<<"].timesPassed = "
       <<tr.workerSummaries[i].timesPassed <<"\n"
       <<"tr.workerSummaries["<<i<<"].timesFailed = "
       <<tr.workerSummaries[i].timesFailed <<"\n"
       <<"tr.workerSummaries["<<i<<"].timesExcept = "
       <<tr.workerSummaries[i].timesExcept <<"\n"
       <<"tr.workerSummaries["<<i<<"].moduleLabel = "
       <<tr.workerSummaries[i].moduleLabel <<"\n";
  }
  
  oss<<"==========================="<<"\n";
  oss<<"== END OF TRIGGER REPORT =="<<"\n";
  oss<<"==========================="<<"\n";
  
  LOG4CPLUS_DEBUG(getApplicationLogger(),oss.str());
}


//______________________________________________________________________________
void FUEventProcessor::defaultWebPage(xgi::Input  *in, xgi::Output *out)
  throw (xgi::exception::Exception)
{
  string urn = getApplicationDescriptor()->getURN();
  ostringstream ourl;
  ourl << "'/" <<  urn << "/microState'";
  *out << "<!-- base href=\"/" <<  urn
       << "\"> -->" << endl;
  *out << "<html>"                                                   << endl;
  *out << "<head>"                                                   << endl;
  //insert javascript code
  jsGen(in,out,ourl.str());
  *out << "<STYLE type=\"text/css\"> #T1 {border-width: 2px; border: solid blue; text-align: center} </STYLE> "                                      << endl; 
  *out << "<link type=\"text/css\" rel=\"stylesheet\"";
  *out << " href=\"/" <<  urn
       << "/styles.css\"/>"                   << endl;
  *out << "<title>" << getApplicationDescriptor()->getClassName() 
       << getApplicationDescriptor()->getInstance() 
       << " MAIN</title>"     << endl;
  *out << "</head>"                                                  << endl;
  *out << "<body onload=\"loadXMLDoc()\">"                           << endl;
  *out << "<table border=\"0\" width=\"100%\">"                      << endl;
  *out << "<tr>"                                                     << endl;
  *out << "  <td align=\"left\">"                                    << endl;
  *out << "    <img"                                                 << endl;
  *out << "     align=\"middle\""                                    << endl;
  *out << "     src=\"/daq/evb/examples/fu/images/fu64x64.gif\""     << endl;
  *out << "     alt=\"main\""                                        << endl;
  *out << "     width=\"64\""                                        << endl;
  *out << "     height=\"64\""                                       << endl;
  *out << "     border=\"\"/>"                                       << endl;
  *out << "    <b>"                                                  << endl;
  *out << getApplicationDescriptor()->getClassName() 
       << getApplicationDescriptor()->getInstance()                  << endl;
  *out << "      " << fsm_.stateName()->toString()                   << endl;
  *out << "    </b>"                                                 << endl;
  *out << "  </td>"                                                  << endl;
  *out << "  <td width=\"32\">"                                      << endl;
  *out << "    <a href=\"/urn:xdaq-application:lid=3\">"             << endl;
  *out << "      <img"                                               << endl;
  *out << "       align=\"middle\""                                  << endl;
  *out << "       src=\"/daq/xdaq/hyperdaq/images/HyperDAQ.jpg\""    << endl;
  *out << "       alt=\"HyperDAQ\""                                  << endl;
  *out << "       width=\"32\""                                      << endl;
  *out << "       height=\"32\""                                     << endl;
  *out << "       border=\"\"/>"                                     << endl;
  *out << "    </a>"                                                 << endl;
  *out << "  </td>"                                                  << endl;
  *out << "  <td width=\"32\">"                                      << endl;
  *out << "  </td>"                                                  << endl;
  *out << "  <td width=\"32\">"                                      << endl;
  *out << "    <a href=\"/" << urn 
       << "/debug\">"                                                << endl;
  *out << "      <img"                                               << endl;
  *out << "       align=\"middle\""                                  << endl;
  *out << "       src=\"/daq/evb/bu/images/debug32x32.gif\""         << endl;
  *out << "       alt=\"debug\""                                     << endl;
  *out << "       width=\"32\""                                      << endl;
  *out << "       height=\"32\""                                     << endl;
  *out << "       border=\"\"/>"                                     << endl;
  *out << "    </a>"                                                 << endl;
  *out << "  </td>"                                                  << endl;
  *out << "</tr>"                                                    << endl;
  *out << "</table>"                                                 << endl;
  
  *out << "<hr/>"                                                    << endl;
  *out << "<table>"                                                  << endl;
  *out << "<tr valign=\"top\">"                                      << endl;
  *out << "  <td>"                                                   << endl;
  *out << "<div id=\"T1\" style=\"border:2px solid blue;height:80;width:150\">microState</div><br /> " << endl;
  *out << "  </td>"                                                  << endl;

  *out << "  <td>"                                                   << endl;
  if(evtProcessor_)
    taskWebPage(in,out,urn);
  else
    *out << "Unconfigured" << endl;
  *out << "  </td>"                                                  << endl;
  *out << "</table>"                                                 << endl;
  
  *out << "<textarea rows=" << 10 << " cols=80 scroll=yes>"          << endl;
  *out << configString_.value_                                       << endl;
  *out << "</textarea><P>"                                           << endl;
  
  *out << "</body>"                                                  << endl;
  *out << "</html>"                                                  << endl;

}


//______________________________________________________________________________
void FUEventProcessor::taskWebPage(xgi::Input *in, xgi::Output *out,const string &urn)
{
  evf::filter *filt = 0;
  ModuleWebRegistry *mwr = 0;
  edm::ServiceRegistry::Operate operate(evtProcessor_->getToken());
  vector<edm::ModuleDescription const*> descs_ =
    evtProcessor_->getAllModuleDescriptions();

  try{
    if(edm::Service<ModuleWebRegistry>().isAvailable())
      mwr = edm::Service<ModuleWebRegistry>().operator->();
  }
  catch(...) {
    cout<<"exception when trying to get the service registry " << endl;
  }

  *out << "<table frame=\"void\" rules=\"groups\" class=\"states\">" << endl;
  *out << "<colgroup> <colgroup align=\"rigth\">"                    << endl;
  *out << "  <tr>"                                                   << endl;
  *out << "    <th colspan=2>"                                       << endl;
  *out << "      " << "Configuration"                                << endl;
  *out << "    </th>"                                                << endl;
  *out << "  </tr>"                                                  << endl;
  
  *out << "<tr>" << endl;
  *out << "<th >" << endl;
  *out << "Parameter" << endl;
  *out << "</th>" << endl;
  *out << "<th>" << endl;
  *out << "Value" << endl;
  *out << "</th>" << endl;
  *out << "</tr>" << endl;
  *out << "<tr>" << endl;
  *out << "<td >" << endl;
  *out << "EP state" << endl;
  *out << "</td>" << endl;
  *out << "<td>" << endl;
  *out << evtProcessor_->getState() << endl;
  *out << "</td>" << endl;
  *out << "</tr>"                                            << endl;
  *out << "<tr>" << endl;
  *out << "<td>" << endl;
  *out << "edm::EP initialized" << endl;
  *out << "</td>" << endl;
  *out << "<td>" << endl;
  *out << epInitialized_ << endl;
  *out << "</td>" << endl;
  *out << "  </tr>"                                            << endl;
  *out << "<tr>" << endl;
  *out << "<td >" << endl;
  *out << "Processed Events/Accepted Events" << endl;
  *out << "</td>" << endl;
  *out << "<td>" << endl;
  *out << evtProcessor_->totalEvents() << "/" 
       << evtProcessor_->totalEventsPassed() << endl;
  *out << "</td>" << endl;
  *out << "  </tr>"                                            << endl;
  *out << "<tr>" << endl;
  *out << "<td >" << endl;
  *out << "Endpaths State" << endl;
  *out << "</td>" << endl;
  *out << "<td";
  *out << (evtProcessor_->endPathsEnabled() ?  "> enabled" : 
	   " bgcolor=\"red\"> disabled" ) << endl;
  //*out << "> N/A this version" << endl;
  *out << "</td>" << endl;
  *out << "  </tr>"                                            << endl;
  *out << "<tr>" << endl;
  *out << "<td >" << endl;
  *out << "Global Input Prescale" << endl;
  *out << "</td>" << endl;
  *out << "<td";
  //*out << (sched_->global_input_prescale_!=1 ? " bgcolor=\"red\">" : ">") << endl;
  //  *out <<  sched_->global_input_prescale_ << endl;
  *out << "> N/A this version" << endl;
  *out << "</td>" << endl;
  *out << "  </tr>"                                            << endl;
  *out << "<tr>" << endl;
  *out << "<td >" << endl;
  *out << "Global Output Prescale" << endl;
  *out << "</td>" << endl;
  *out << "<td";
  //*out  << (sched_->global_output_prescale_!=1 ? " bgcolor=\"red\">" : ">") << endl;
  //  *out <<  sched_->global_output_prescale_ << endl;
  *out << ">N/A this version" << endl;
  *out << "</td>" << endl;
  *out << "  </tr>"                                            << endl;
  
  
  
  *out << "</table>" << endl;
  
  *out << "<table frame=\"void\" rules=\"rows\" class=\"modules\">" << endl;
  *out << "  <tr>"                                                   << endl;
  *out << "    <th colspan=3>"                                       << endl;
  *out << "      " << "Application"                                  << endl;
  
  if(descs_.size()>0)
    *out << " (Process name=" << descs_[0]->processName() << ")"       << endl;
  
  
  
  *out << "    </th>"                                                << endl;
  *out << "  </tr>"                                                  << endl;
  
  *out << "<tr >" << endl;
  *out << "<th >" << endl;
  *out << "Module" << endl;
  *out << "</th>" << endl;
  *out << "<th >" << endl;
  *out << "Label" << endl;
  *out << "</th>" << endl;
  *out << "<th >" << endl;
  *out << "Version" << endl;
  *out << "</th>" << endl;
  *out << "</tr>" << endl;
  
  for(unsigned int idesc = 0; idesc < descs_.size(); idesc++)
    {
      *out << "<tr>" << endl;
      *out << "<td >" << endl;
      if(mwr && mwr->checkWeb(descs_[idesc]->moduleName()))
	*out << "<a href=\"/" << urn << "/moduleWeb?module=" << descs_[idesc]->moduleName() << "\">" 
	     << descs_[idesc]->moduleName() << "</a>" << endl;
      else
	*out << descs_[idesc]->moduleName() << endl;
      *out << "</td>" << endl;
      *out << "<td >" << endl;
      *out << descs_[idesc]->moduleLabel() << endl;
      *out << "</td>" << endl;
      *out << "<td >" << endl;
      *out << descs_[idesc]->releaseVersion() << endl;
      *out << "</td>" << endl;
      *out << "</tr>" << endl;
    }
  *out << "</table>" << endl;
  *out << "<table border=1 bgcolor=\"#CFCFCF\">" << endl;
  *out << "<tr>" << endl;

  if(filt) {
    //HLT summary status goes here
  }
  else {      
    *out << "<td >" << endl;
    *out << "No Filter Module" << endl;
    *out << "</td>" << endl;
  }
  *out << "</tr>" << endl;
  *out << "</table>" << endl;
}


//______________________________________________________________________________
void FUEventProcessor::moduleWeb(xgi::Input  *in, xgi::Output *out)
  throw (xgi::exception::Exception)
{
  Cgicc cgi(in);
  vector<FormEntry> el1;
  cgi.getElement("module",el1);
  if(evtProcessor_)  {
    if(el1.size()!=0) {
      string mod = el1[0].getValue();
      edm::ServiceRegistry::Operate operate(evtProcessor_->getToken());
      ModuleWebRegistry *mwr = 0;
      try{
	if(edm::Service<ModuleWebRegistry>().isAvailable())
	  mwr = edm::Service<ModuleWebRegistry>().operator->();
      }
      catch(...) { 
	cout<<"exception when trying to get the service registry "<<endl;
      }
      mwr->invoke(in,out,mod);
    }
  }
  else {
    *out<<"EventProcessor just disappeared "<<endl;
  }
}


//______________________________________________________________________________
void FUEventProcessor::jsGen(xgi::Input *in, xgi::Output *out, string url)
  throw (xgi::exception::Exception)
{
  *out << "<script type=\"text/javascript\"> \n";
  *out << "var xmlhttp \n";
  *out << " \n";
  *out << "function loadXMLDoc() \n";
  *out << "{ \n";
  *out << "xmlhttp=null \n";
  *out << " \n";
  *out << "if (window.XMLHttpRequest) \n";
  *out << "  { \n";
  *out << "  xmlhttp=new XMLHttpRequest() \n";
  *out << "  } \n";
  *out << " \n";
  *out << "else if (window.ActiveXObject) \n";
  *out << "  { \n";
  *out << "  xmlhttp=new ActiveXObject(\"Microsoft.XMLHTTP\") \n";
  *out << "  } \n";
  *out << "if (xmlhttp!=null) \n";
  *out << "  { \n";
  *out << "  xmlhttp.onreadystatechange=state_Change \n";
  *out << "  xmlhttp.open(\"GET\"," << url << ",true) \n";
  *out << "  xmlhttp.send(null) \n";
  *out << "  setTimeout('loadXMLDoc()',500) \n";
  *out << "  } \n";
  *out << "else \n";
  *out << "  { \n";
  *out << "  alert(\"Your browser does not support XMLHTTP.\") \n";
  *out << "  } \n";
  *out << "} \n";
  *out << " \n";
  *out << "function state_Change() \n";
  *out << "{ \n";
  // if xmlhttp shows "loaded"
  *out << "if (xmlhttp.readyState==4) \n";
  *out << "  { \n";
  // if "OK" 
  *out << " if (xmlhttp.status==200) \n";
  *out << "  { \n";
  *out << "  document.getElementById('T1').innerHTML=xmlhttp.responseText \n";
  *out << "  } \n";
  *out << "  else \n";
  *out << "  { \n";
  *out << "  document.getElementById('T1').innerHTML=xmlhttp.statusText \n";
  *out << "  } \n";
  *out << "  } \n";
  *out << "} \n";
  *out << " \n";
  *out << "</script> \n";
}


//______________________________________________________________________________
void FUEventProcessor::microState(xgi::Input  *in, xgi::Output *out)
  throw (xgi::exception::Exception)
{
  edm::ServiceRegistry::Operate operate(serviceToken_);
  MicroStateService *mss = 0;
  string micro1 = "unavailable";
  string micro2 = "unavailable";
  try{
    mss = edm::Service<MicroStateService>().operator->();
  }
  catch(...) { 
    LOG4CPLUS_INFO(getApplicationLogger(),
		   "exception when trying to get service MicroStateService");
  }
  if(mss) {
    micro1 = mss->getMicroState1();
    micro2 = mss->getMicroState2();
  }
  
  *out << "<br>  " << micro1 << endl;
  *out << "<br>  " << micro2 << endl;
}


XDAQ_INSTANTIATOR_IMPL(evf::FUEventProcessor)
