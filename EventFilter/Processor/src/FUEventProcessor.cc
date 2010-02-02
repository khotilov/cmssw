////////////////////////////////////////////////////////////////////////////////
//
// FUEventProcessor
// ----------------
//
////////////////////////////////////////////////////////////////////////////////

#include "FUEventProcessor.h"
#include "procUtils.h"


#include "EventFilter/Utilities/interface/Exception.h"

#include "EventFilter/Message2log4cplus/interface/MLlog4cplus.h"
#include "EventFilter/Modules/interface/FUShmDQMOutputService.h"
#include "EventFilter/Utilities/interface/ServiceWebRegistry.h"
#include "EventFilter/Utilities/interface/ServiceWeb.h"


#include "FWCore/PluginManager/interface/ProblemTracker.h"
#include "FWCore/PluginManager/interface/PresenceFactory.h"
#include "FWCore/Utilities/interface/Presence.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Version/interface/GetReleaseVersion.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <boost/tokenizer.hpp>

// to handle pt file descriptors left open at fork
#include "pt/PeerTransportReceiver.h"
#include "pt/PeerTransportAgent.h"

#include "xcept/tools.h"
#include "xgi/Method.h"

#include "cgicc/CgiDefs.h"
#include "cgicc/Cgicc.h"
#include "cgicc/FormEntry.h"


#include <sys/wait.h>
#include <sys/utsname.h>

#include <typeinfo>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

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
  , log_(getApplicationLogger())
  , evtProcessor_(log_, getApplicationDescriptor()->getInstance())
  , runNumber_(0)
  , epInitialized_(false)
  , outPut_(true)
  , autoRestartSlaves_(false)
  , slaveRestartDelaySecs_(10)
  , hasShMem_(true)
  , hasPrescaleService_(true)
  , hasModuleWebRegistry_(true)
  , hasServiceWebRegistry_(true)
  , isRunNumberSetter_(true)
  , outprev_(true)
  , reasonForFailedState_()
  , squidnet_(3128,"http://localhost:8000/RELEASE-NOTES.txt")
  , logRing_(logRingSize_)
  , logRingIndex_(logRingSize_)
  , logWrap_(false)
  , nbSubProcesses_(0)
  , sq_(0)
  , nblive_(0)
  , nbdead_(0)
  , wlReceiving_(0)
  , asReceiveMsgAndExecute_(0)
  , receiving_(false) 
  , wlReceivingMonitor_(0)
  , asReceiveMsgAndRead_(0)
  , receivingM_(false)
  , isChildProcess_(false)
  , wlSupervising_(0)
  , asSupervisor_(0)
  , supervising_(false)
  , monitorInfoSpace_(0)
  , applicationInfoSpace_(0)
  , nbProcessed(0)
  , nbAccepted(0)
  , scalersInfoSpace_(0)
  , wlScalers_(0)
  , asScalers_(0)
  , wlScalersActive_(false)
  , superSleepSec_(1)
{
  names_.push_back("nbProcessed"    );
  names_.push_back("nbAccepted"     );
  names_.push_back("epMacroStateInt");
  names_.push_back("epMicroStateInt");
  // create pipe for web communication
  int retpipe = pipe(anonymousPipe_);
  if(retpipe != 0)
        LOG4CPLUS_ERROR(getApplicationLogger(),"Failed to create pipe");
  // check squid presence
  squidPresent_ = squidnet_.check();
  //pass application parameters to FWEPWrapper
  evtProcessor_.setAppDesc(getApplicationDescriptor());
  evtProcessor_.setAppCtxt(getApplicationContext());
  // bind relevant callbacks to finite state machine
  fsm_.initialize<evf::FUEventProcessor>(this);
  
  //set sourceId_
  url_ =
    getApplicationDescriptor()->getContextDescriptor()->getURL()+"/"+
    getApplicationDescriptor()->getURN();
  class_   =getApplicationDescriptor()->getClassName();
  instance_=getApplicationDescriptor()->getInstance();
  sourceId_=class_.toString()+"_"+instance_.toString();
  LOG4CPLUS_INFO(getApplicationLogger(),sourceId_ <<" constructor"         );
  LOG4CPLUS_INFO(getApplicationLogger(),"CMSSW_BASE:"<<getenv("CMSSW_BASE"));
  
  getApplicationDescriptor()->setAttribute("icon", "/evf/images/epicon.jpg");
  
  xdata::InfoSpace *ispace = getApplicationInfoSpace();
  applicationInfoSpace_ = ispace;

  // default configuration
  ispace->fireItemAvailable("parameterSet",         &configString_                );
  ispace->fireItemAvailable("epInitialized",        &epInitialized_               );
  ispace->fireItemAvailable("stateName",             fsm_.stateName()             );
  ispace->fireItemAvailable("runNumber",            &runNumber_                   );
  ispace->fireItemAvailable("outputEnabled",        &outPut_                      );

  ispace->fireItemAvailable("hasSharedMemory",      &hasShMem_);
  ispace->fireItemAvailable("hasPrescaleService",   &hasPrescaleService_          );
  ispace->fireItemAvailable("hasModuleWebRegistry", &hasModuleWebRegistry_        );
  ispace->fireItemAvailable("hasServiceWebRegistry", &hasServiceWebRegistry_      );
  ispace->fireItemAvailable("isRunNumberSetter",    &isRunNumberSetter_           );
  ispace->fireItemAvailable("rcmsStateListener",     fsm_.rcmsStateListener()     );
  ispace->fireItemAvailable("foundRcmsStateListener",fsm_.foundRcmsStateListener());
  ispace->fireItemAvailable("nbSubProcesses",       &nbSubProcesses_              );
  ispace->fireItemAvailable("superSleepSec",        &superSleepSec_               );
  ispace->fireItemAvailable("autoRestartSlaves",    &autoRestartSlaves_           );
  ispace->fireItemAvailable("slaveRestartDelaySecs",&slaveRestartDelaySecs_       );
  
  // Add infospace listeners for exporting data values
  getApplicationInfoSpace()->addItemChangedListener("parameterSet",        this);
  getApplicationInfoSpace()->addItemChangedListener("outputEnabled",       this);

  // findRcmsStateListener
  fsm_.findRcmsStateListener();
  
  // initialize monitoring infospace

  std::stringstream oss2;
  oss2<<"urn:xdaq-monitorable-"<<class_.toString();
  string monInfoSpaceName=oss2.str();
  toolbox::net::URN urn = this->createQualifiedInfoSpace(monInfoSpaceName);
  monitorInfoSpace_ = xdata::getInfoSpaceFactory()->get(urn.toString());

  
  monitorInfoSpace_->fireItemAvailable("url",                      &url_            );
  monitorInfoSpace_->fireItemAvailable("class",                    &class_          );
  monitorInfoSpace_->fireItemAvailable("instance",                 &instance_       );
  monitorInfoSpace_->fireItemAvailable("runNumber",                &runNumber_      );
  monitorInfoSpace_->fireItemAvailable("stateName",                 fsm_.stateName()); 

  monitorInfoSpace_->fireItemAvailable("squidPresent",             &squidPresent_   );

  std::stringstream oss3;
  oss3<<"urn:xdaq-scalers-"<<class_.toString();
  string monInfoSpaceName2=oss3.str();
  toolbox::net::URN urn2 = this->createQualifiedInfoSpace(monInfoSpaceName2);

  xdata::InfoSpace *scalersInfoSpace_ = xdata::getInfoSpaceFactory()->get(urn2.toString());
  evtProcessor_.setScalersInfoSpace(scalersInfoSpace_);
  scalersInfoSpace_->fireItemAvailable("instance", &instance_);

  evtProcessor_.setApplicationInfoSpace(ispace);
  evtProcessor_.setMonitorInfoSpace(monitorInfoSpace_);
  evtProcessor_.publishConfigAndMonitorItems(nbSubProcesses_.value_);

  //subprocess state vectors for MP
  monitorInfoSpace_->fireItemAvailable("epMacroStateInt",             &spMStates_); 
  monitorInfoSpace_->fireItemAvailable("epMicroStateInt",             &spmStates_); 
  
  // Bind web interface
  xgi::bind(this, &FUEventProcessor::css,              "styles.css");
  xgi::bind(this, &FUEventProcessor::defaultWebPage,   "Default"   );
  xgi::bind(this, &FUEventProcessor::spotlightWebPage, "Spotlight" );
  xgi::bind(this, &FUEventProcessor::scalersWeb,       "scalersWeb");
  xgi::bind(this, &FUEventProcessor::pathNames,        "pathNames" );
  xgi::bind(this, &FUEventProcessor::subWeb,           "SubWeb"    );
  xgi::bind(this, &FUEventProcessor::moduleWeb,        "moduleWeb" );
  xgi::bind(this, &FUEventProcessor::serviceWeb,       "serviceWeb");
  xgi::bind(this, &FUEventProcessor::microState,       "microState");
  xgi::bind(this, &FUEventProcessor::updater,          "updater"   );
  xgi::bind(this, &FUEventProcessor::procStat,         "procStat"  );

  // instantiate the plugin manager, not referenced here after!

  edm::AssertHandler ah;

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
  } 
  catch(...) {
    LOG4CPLUS_ERROR(getApplicationLogger(),"Unknown Exception");
  }
  ML::MLlog4cplus::setAppl(this);      

  typedef set<xdaq::ApplicationDescriptor*> AppDescSet_t;
  typedef AppDescSet_t::iterator            AppDescIter_t;
  
  AppDescSet_t rcms=
    getApplicationContext()->getDefaultZone()->
    getApplicationDescriptors("RCMSStateListener");
  if(rcms.size()==0) 
    {
      LOG4CPLUS_WARN(getApplicationLogger(),
		       "MonitorReceiver not found, perhaphs it has not been defined ? Scalers updater wl will bail out!");
      //	localLog("-W- MonitorReceiver not found, perhaphs it has not been defined ? Scalers updater wl will bail out!");
    }
  else
    {
      AppDescIter_t it = rcms.begin();
      evtProcessor_.setRcms(*it);
    }
  pthread_mutex_init(&start_lock_,0);
  pthread_mutex_init(&stop_lock_,0);
  pthread_mutex_init(&pickup_lock_,0);

  std::ostringstream ost;
  ost  << "<div id=\"ve\">2.0.10 (" << edm::getReleaseVersion() <<")</div>"
       << "<div id=\"ou\">" << outPut_.toString() << "</div>"
       << "<div id=\"sh\">" << hasShMem_.toString() << "</div>"
       << "<div id=\"mw\">" << hasModuleWebRegistry_.toString() << "</div>"
       << "<div id=\"sw\">" << hasServiceWebRegistry_.toString() << "</div>";
  
  xdata::Serializable *monsleep = 0;
  xdata::Serializable *lstimeout = 0;
  try{
    monsleep = ispace->find("monSleepSec");
    lstimeout = ispace->find("lsTimeOut");
  }
  catch(xdata::exception::Exception e){
  }
  
  if(monsleep!=0)
    ost << "<div id=\"ms\">" << monsleep->toString() << "</div>";
  if(lstimeout!=0)
    ost << "<div id=\"lst\">" << lstimeout->toString() << "</div>";
  char cbuf[sizeof(struct utsname)];
  struct utsname* buf = (struct utsname*)cbuf;
  uname(buf);
  ost << "<div id=\"sysinfo\">" << buf->sysname << " " << buf->nodename 
      << " " << buf->release << " " << buf->version << " " << buf->machine << "</div>";
  updaterStatic_ = ost.str();
  startSupervisorLoop();  
}



//______________________________________________________________________________
FUEventProcessor::~FUEventProcessor()
{
  // no longer needed since the wrapper is a member of the class and one can rely on 
  // implicit destructor - to be revised
  //  if (evtProcessor_) delete evtProcessor_;
}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
bool FUEventProcessor::configuring(toolbox::task::WorkLoop* wl)
{
  std::cout << "values " << ((nbSubProcesses_.value_!=0) ? 0x10 : 0) << " "
	    << ((instance_.value_==0) ? 0x8 : 0) << " "
	    << (hasServiceWebRegistry_.value_ ? 0x4 : 0) << " "
	    << (hasModuleWebRegistry_.value_ ? 0x2 : 0) << " "
	    << (hasPrescaleService_.value_ ? 0x1 : 0) <<std::endl;
  unsigned short smap 
    = ((nbSubProcesses_.value_!=0) ? 0x10 : 0)
    + ((instance_.value_==0) ? 0x8 : 0)
    + (hasServiceWebRegistry_.value_ ? 0x4 : 0) 
    + (hasModuleWebRegistry_.value_ ? 0x2 : 0) 
    + (hasPrescaleService_.value_ ? 0x1 : 0);

  if(nbSubProcesses_.value_==0) 
    {
      spMStates_.setSize(1); 
      spmStates_.setSize(1); 
    }
  else
    {
      spMStates_.setSize(nbSubProcesses_.value_);
      spmStates_.setSize(nbSubProcesses_.value_);
      for(unsigned int i = 0; i < spMStates_.size(); i++)
	{
	  spMStates_[i] = edm::event_processor::sInit; 
	  spmStates_[i] = 0; 
	}
    }
  try {
    LOG4CPLUS_INFO(getApplicationLogger(),"Start configuring ...");
    std::string cfg = configString_.toString(); evtProcessor_.init(smap,cfg);
    epInitialized_=true;
    if(evtProcessor_)
      {
	// moved to wrapper class
	configuration_ = evtProcessor_.configuration();
	if(nbSubProcesses_.value_==0) evtProcessor_.startMonitoringWorkLoop(); 
	evtProcessor_->beginJob(); 
	fsm_.fireEvent("ConfigureDone",this);
	LOG4CPLUS_INFO(getApplicationLogger(),"Finished configuring!");
	localLog("-I- Configuration completed");
      }
  }
  catch (xcept::Exception &e) {
    reasonForFailedState_ = "configuring FAILED: " + (string)e.what();
    fsm_.fireFailed(reasonForFailedState_,this);
    localLog(reasonForFailedState_);
  }
  catch(cms::Exception &e) {
    reasonForFailedState_ = e.explainSelf();
    fsm_.fireFailed(reasonForFailedState_,this);
    localLog(reasonForFailedState_);
  }    
  catch(std::exception &e) {
    reasonForFailedState_ = e.what();
    fsm_.fireFailed(reasonForFailedState_,this);
    localLog(reasonForFailedState_);
  }
  catch(...) {
    fsm_.fireFailed("Unknown Exception",this);
  }
  return false;
}


//______________________________________________________________________________
bool FUEventProcessor::enabling(toolbox::task::WorkLoop* wl)
{
  std::cout << "values " << ((nbSubProcesses_.value_!=0) ? 0x10 : 0) << " "
	    << ((instance_.value_==0) ? 0x8 : 0) << " "
	    << (hasServiceWebRegistry_.value_ ? 0x4 : 0) << " "
	    << (hasModuleWebRegistry_.value_ ? 0x2 : 0) << " "
	    << (hasPrescaleService_.value_ ? 0x1 : 0) <<std::endl;
  unsigned short smap 
    = ((nbSubProcesses_.value_!=0) ? 0x10 : 0)
    + ((instance_.value_==0) ? 0x8 : 0)
    + (hasServiceWebRegistry_.value_ ? 0x4 : 0) 
    + (hasModuleWebRegistry_.value_ ? 0x2 : 0) 
    + (hasPrescaleService_.value_ ? 0x1 : 0);

  LOG4CPLUS_INFO(getApplicationLogger(),"Start enabling ...");
  if(!epInitialized_) evtProcessor_.forceInitEventProcessorMaybe();

  std::string cfg = configString_.toString(); evtProcessor_.init(smap,cfg);
  configuration_ = evtProcessor_.configuration(); // get it again after init has been carried out...
  //classic appl will return here 
  if(nbSubProcesses_.value_==0) return enableClassic();

  //protect manipulation of subprocess array
  pthread_mutex_lock(&start_lock_);
  subs_.clear();

  pid_t retval = -1;
  subs_.resize(nbSubProcesses_.value_);
  for(unsigned int i=0; i<nbSubProcesses_.value_; i++)
    {
      subs_[i]=SubProcess(i,retval); //this will replace all the scattered variables
      retval = subs_[i].forkNew();
      if(retval>0)
	{
	  pthread_mutex_unlock(&start_lock_);
	}
      if(retval==0)
	{
	  isChildProcess_=true;
	  try{
	    pt::PeerTransport * ptr =
	      pt::getPeerTransportAgent()->getPeerTransport("http","soap",pt::Receiver);
	    delete ptr;
	  }
	  catch (pt::exception::PeerTransportNotFound & e ){
	    LOG4CPLUS_WARN(getApplicationLogger()," ***Slave Failed to shutdown ptHTTP");
	  }
	  fsm_.disableRcmsStateNotification();
	  ostringstream ost1;
	  ost1 << "-I- Slave Process " << retval << " forked for slot " << i; 
	  localLog(ost1.str());
	  return enableMPEPSlave(i);
	  // the loop is broken in the child 
	}
    }
  
  startScalersWorkLoop();
  LOG4CPLUS_INFO(getApplicationLogger(),"Finished enabling!");
  fsm_.fireEvent("EnableDone",this);
  localLog("-I- Start completed");
  return false;
}


//______________________________________________________________________________
bool FUEventProcessor::stopping(toolbox::task::WorkLoop* wl)
{
  if(nbSubProcesses_.value_!=0) 
    stopSlavesAndAcknowledge();
  return stopClassic();
}


//______________________________________________________________________________
bool FUEventProcessor::halting(toolbox::task::WorkLoop* wl)
{
  LOG4CPLUS_INFO(getApplicationLogger(),"Start halting ...");
  if(nbSubProcesses_.value_!=0) 
    stopSlavesAndAcknowledge();
  try{
    evtProcessor_.stopAndHalt();
  }
  catch (evf::Exception &e) {
    reasonForFailedState_ = "halting FAILED: " + (string)e.what();
    localLog(reasonForFailedState_);
    fsm_.fireFailed(reasonForFailedState_,this);
  }
  //  if(hasShMem_) detachDqmFromShm();

  LOG4CPLUS_INFO(getApplicationLogger(),"Finished halting!");
  fsm_.fireEvent("HaltDone",this);

  localLog("-I- Halt completed");
  return false;
}


//______________________________________________________________________________
xoap::MessageReference FUEventProcessor::fsmCallback(xoap::MessageReference msg)
  throw (xoap::exception::Exception)
{
  return fsm_.commandCallback(msg);
}


//______________________________________________________________________________
void FUEventProcessor::actionPerformed(xdata::Event& e)
{

  if (e.type()=="ItemChangedEvent" && fsm_.stateName()->toString()!="Halted") {
    string item = dynamic_cast<xdata::ItemChangedEvent&>(e).itemName();
    
    if ( item == "parameterSet") {
      LOG4CPLUS_WARN(getApplicationLogger(),
		     "HLT Menu changed, force reinitialization of EventProcessor");
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
      LOG4CPLUS_WARN(getApplicationLogger(),
		     "Setting global input prescale has no effect "
		     <<"in this version of the code");
    }
    if ( item == "globalOutputPrescale") {
      LOG4CPLUS_WARN(getApplicationLogger(),
		     "Setting global output prescale has no effect "
		     <<"in this version of the code");
    }
  }

}

//______________________________________________________________________________
void FUEventProcessor::subWeb(xgi::Input  *in, xgi::Output *out)
{
  using namespace cgicc;
  pid_t pid = 0;
  std::ostringstream ost;
  ost << "&";

  Cgicc cgi(in);
  internal::MyCgi *mycgi = (internal::MyCgi*)in;
  for(std::map<std::string, std::string, std::less<std::string> >::iterator mit = 
	mycgi->getEnvironment().begin();
      mit != mycgi->getEnvironment().end(); mit++)
    ost << mit->first << "%" << mit->second << ";";
  std::vector<FormEntry> els = cgi.getElements() ;
  std::vector<FormEntry> el1;
  cgi.getElement("method",el1);
  std::vector<FormEntry> el2;
  cgi.getElement("process",el2);
  if(el1.size()!=0) {
    std::string meth = el1[0].getValue();
    if(el2.size()!=0) {
      unsigned int i = 0;
      std::string mod = el2[0].getValue();
      pid = atoi(mod.c_str()); // get the process id to be polled
      for(; i < subs_.size(); i++)
	if(subs_[i].pid()==pid) break;
      if(i>=subs_.size()){ //process was not found, let the browser know
	*out << "ERROR 404 : Process " << pid << " Not Found !" << std::endl;
	return;
      } 
      MsgBuf msg1(meth.length()+ost.str().length(),MSQM_MESSAGE_TYPE_WEB);
      strcpy(msg1->mtext,meth.c_str());
      strcpy(msg1->mtext+meth.length(),ost.str().c_str());
      subs_[i].post(msg1,true);
      MsgBuf msg(MAX_MSG_SIZE,MSQS_MESSAGE_TYPE_WEB);
      bool done = false;
      std::vector<char *>pieces;
      while(!done){
	subs_[i].rcv(msg,true);
	unsigned int nbytes = atoi(msg->mtext);
	if(nbytes < MAX_PIPE_BUFFER_SIZE) done = true; // this will break the while loop
	char *buf= new char[nbytes];
	read(anonymousPipe_[PIPE_READ],buf,nbytes);
	pieces.push_back(buf);
      }
      for(unsigned int j = 0; j < pieces.size(); j++){
	*out<<pieces[j];    // chain the buffers into the output strstream
	delete[] pieces[j]; //make sure to release all buffers used for reading the pipe
      }
    }
  }
}


//______________________________________________________________________________
void FUEventProcessor::defaultWebPage(xgi::Input  *in, xgi::Output *out)
  throw (xgi::exception::Exception)
{


  *out << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">" 
       << "<html><head><title>" << getApplicationDescriptor()->getClassName() << (nbSubProcesses_.value_ > 0 ? "MP " : " ") 
       << getApplicationDescriptor()->getInstance() << "</title>"
       << "<meta http-equiv=\"REFRESH\" content=\"0;url=/evf/html/defaultBasePage.html\">"
       << "</head></html>";
}


//______________________________________________________________________________


void FUEventProcessor::spotlightWebPage(xgi::Input  *in, xgi::Output *out)
  throw (xgi::exception::Exception)
{

  string urn = getApplicationDescriptor()->getURN();

  *out << "<!-- base href=\"/" <<  urn
       << "\"> -->" << endl;
  *out << "<html>"                                                   << endl;
  *out << "<head>"                                                   << endl;
  *out << "<link type=\"text/css\" rel=\"stylesheet\"";
  *out << " href=\"/evf/html/styles.css\"/>"                   << endl;
  *out << "<title>" << getApplicationDescriptor()->getClassName() 
       << getApplicationDescriptor()->getInstance() 
       << " MAIN</title>"     << endl;
  *out << "</head>"                                                  << endl;
  *out << "<body>"                                                   << endl;
  *out << "<table border=\"0\" width=\"100%\">"                      << endl;
  *out << "<tr>"                                                     << endl;
  *out << "  <td align=\"left\">"                                    << endl;
  *out << "    <img"                                                 << endl;
  *out << "     align=\"middle\""                                    << endl;
  *out << "     src=\"/evf/images/spoticon.jpg\""			     << endl;
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
  *out << "       src=\"/hyperdaq/images/HyperDAQ.jpg\""             << endl;
  *out << "       alt=\"HyperDAQ\""                                  << endl;
  *out << "       width=\"32\""                                      << endl;
  *out << "       height=\"32\""                                     << endl;
  *out << "       border=\"\"/>"                                     << endl;
  *out << "    </a>"                                                 << endl;
  *out << "  </td>"                                                  << endl;
  *out << "  <td width=\"32\">"                                      << endl;
  *out << "  </td>"                                                  << endl;
  *out << "  <td width=\"32\">"                                      << endl;
  *out << "    <a href=\"/" << urn << "/\">"                         << endl;
  *out << "      <img"                                               << endl;
  *out << "       align=\"middle\""                                  << endl;
  *out << "       src=\"/evf/images/epicon.jpg\""		     << endl;
  *out << "       alt=\"main\""                                      << endl;
  *out << "       width=\"32\""                                      << endl;
  *out << "       height=\"32\""                                     << endl;
  *out << "       border=\"\"/>"                                     << endl;
  *out << "    </a>"                                                 << endl;
  *out << "  </td>"                                                  << endl;
  *out << "</tr>"                                                    << endl;
  *out << "</table>"                                                 << endl;

  *out << "<hr/>"                                                    << endl;
  
  std::ostringstream ost;
  if(isChildProcess_) 
    ost << "/SubWeb?process=" << getpid() << "&method=moduleWeb&";
  else
    ost << "/moduleWeb?";
  urn += ost.str();
  if(evtProcessor_ && (isChildProcess_ || nbSubProcesses_.value_==0))
    evtProcessor_.taskWebPage(in,out,urn);
  else if(evtProcessor_)
    evtProcessor_.summaryWebPage(in,out,urn);
  else
    *out << "<td>HLT Unconfigured</td>" << endl;
  *out << "</table>"                                                 << endl;
  
  *out << "<br><textarea rows=" << 10 << " cols=80 scroll=yes>"      << endl;
  *out << configuration_                                             << endl;
  *out << "</textarea><P>"                                           << endl;
  
  *out << "</body>"                                                  << endl;
  *out << "</html>"                                                  << endl;


}
void FUEventProcessor::scalersWeb(xgi::Input  *in, xgi::Output *out)
  throw (xgi::exception::Exception)
{

  out->getHTTPResponseHeader().addHeader( "Content-Type",
					  "application/octet-stream" );
  out->getHTTPResponseHeader().addHeader( "Content-Transfer-Encoding",
					  "binary" );
  if(evtProcessor_ != 0){
    out->write( (char*)(evtProcessor_.getPackedTriggerReport()->mtext), sizeof(TriggerReportStatic) );
  }
}

void FUEventProcessor::pathNames(xgi::Input  *in, xgi::Output *out)
  throw (xgi::exception::Exception)
{

  if(evtProcessor_ != 0){
    xdata::Serializable *legenda = scalersInfoSpace_->find("scalersLegenda");
    if(legenda !=0){
      std::string slegenda = ((xdata::String*)legenda)->value_;
      *out << slegenda << std::endl;
    }
  }
}

void FUEventProcessor::attachDqmToShm() throw (evf::Exception)  
{
  string errmsg;
  bool success = false;
  try {
    edm::ServiceRegistry::Operate operate(evtProcessor_->getToken());
    if(edm::Service<FUShmDQMOutputService>().isAvailable())
      success = edm::Service<FUShmDQMOutputService>()->attachToShm();
    if (!success) errmsg = "Failed to attach DQM service to shared memory";
  }
  catch (cms::Exception& e) {
    errmsg = "Failed to attach DQM service to shared memory: " + (string)e.what();
  }
  if (!errmsg.empty()) XCEPT_RAISE(evf::Exception,errmsg);
}



void FUEventProcessor::detachDqmFromShm() throw (evf::Exception)
{
  string errmsg;
  bool success = false;
  try {
    edm::ServiceRegistry::Operate operate(evtProcessor_->getToken());
    if(edm::Service<FUShmDQMOutputService>().isAvailable())
      success = edm::Service<FUShmDQMOutputService>()->detachFromShm();
    if (!success) errmsg = "Failed to detach DQM service from shared memory";
  }
  catch (cms::Exception& e) {
    errmsg = "Failed to detach DQM service from shared memory: " + (string)e.what();
  }
  if (!errmsg.empty()) XCEPT_RAISE(evf::Exception,errmsg);
}


std::string FUEventProcessor::logsAsString()
{
  ostringstream oss;
  if(logWrap_)
    {
      for(unsigned int i = logRingIndex_; i < logRing_.size(); i++)
	oss << logRing_[i] << std::endl;
      for(unsigned int i = 0; i <  logRingIndex_; i++)
	oss << logRing_[i] << std::endl;
    }
  else
      for(unsigned int i = logRingIndex_; i < logRing_.size(); i++)
	oss << logRing_[i] << std::endl;
    
  return oss.str();
}
  
void FUEventProcessor::localLog(string m)
{
  timeval tv;

  gettimeofday(&tv,0);
  tm *uptm = localtime(&tv.tv_sec);
  char datestring[256];
  strftime(datestring, sizeof(datestring),"%c", uptm);

  if(logRingIndex_ == 0){logWrap_ = true; logRingIndex_ = logRingSize_;}
  logRingIndex_--;
  ostringstream timestamp;
  timestamp << " at " << datestring;
  m += timestamp.str();
  logRing_[logRingIndex_] = m;
}

void FUEventProcessor::startSupervisorLoop()
{
  try {
    wlSupervising_=
      toolbox::task::getWorkLoopFactory()->getWorkLoop("Supervisor",
						       "waiting");
    if (!wlSupervising_->isActive()) wlSupervising_->activate();
    asSupervisor_ = toolbox::task::bind(this,&FUEventProcessor::supervisor,
					"Supervisor");
    wlSupervising_->submit(asSupervisor_);
    supervising_ = true;
  }
  catch (xcept::Exception& e) {
    std::string msg = "Failed to start workloop 'Supervisor'.";
    XCEPT_RETHROW(evf::Exception,msg,e);
  }
}

void FUEventProcessor::startReceivingLoop()
{
  try {
    wlReceiving_=
      toolbox::task::getWorkLoopFactory()->getWorkLoop("Receiving",
						       "waiting");
    if (!wlReceiving_->isActive()) wlReceiving_->activate();
    asReceiveMsgAndExecute_ = toolbox::task::bind(this,&FUEventProcessor::receiving,
					"Receiving");
    wlReceiving_->submit(asReceiveMsgAndExecute_);
    receiving_ = true;
  }
  catch (xcept::Exception& e) {
    std::string msg = "Failed to start workloop 'Receiving'.";
    XCEPT_RETHROW(evf::Exception,msg,e);
  }
}
void FUEventProcessor::startReceivingMonitorLoop()
{
  try {
    wlReceivingMonitor_=
      toolbox::task::getWorkLoopFactory()->getWorkLoop("ReceivingM",
						       "waiting");
    if (!wlReceivingMonitor_->isActive()) 
      wlReceivingMonitor_->activate();
    asReceiveMsgAndRead_ = 
      toolbox::task::bind(this,&FUEventProcessor::receivingAndMonitor,
			  "ReceivingM");
    wlReceivingMonitor_->submit(asReceiveMsgAndRead_);
    receivingM_ = true;
  }
  catch (xcept::Exception& e) {
    std::string msg = "Failed to start workloop 'ReceivingM'.";
    XCEPT_RETHROW(evf::Exception,msg,e);
  }
}

bool FUEventProcessor::receiving(toolbox::task::WorkLoop *)
{
  MsgBuf msg;
  try{
    sq_->rcv(msg); //will receive only messages from Master
    if(msg->mtype==MSQM_MESSAGE_TYPE_STOP)
      {
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
	} 
	catch(...) {
	  LOG4CPLUS_ERROR(getApplicationLogger(),"Unknown Exception");
	}
	
	fsm_.fireEvent("Stop",this); // need to set state in fsm first to allow stopDone transition
	stopClassic(); // call the normal sequence of stopping - as this is allowed to fail provisions must be made ...@@@EM
	// tried to cure the coredump on exit(0) for 227 by orderly halting first but it does not work...
	//	evtProcessor_.stopAndHalt();
	MsgBuf msg1(0,MSQS_MESSAGE_TYPE_STOP);
	sq_->post(msg1);
	fclose(stdout);
	fclose(stderr);
	exit(EXIT_SUCCESS);
      }
  }
  catch(evf::Exception &e){}
  return true;
}

bool FUEventProcessor::supervisor(toolbox::task::WorkLoop *)
{
  pthread_mutex_lock(&stop_lock_);
  if(subs_.size()!=nbSubProcesses_.value_)
    {
      subs_.resize(nbSubProcesses_.value_);
      spMStates_.resize(nbSubProcesses_.value_);
      spmStates_.resize(nbSubProcesses_.value_);
      for(unsigned int i = 0; i < spMStates_.size(); i++)
	{
	  spMStates_[i] = edm::event_processor::sInit; 
	  spmStates_[i] = 0; 
	}
    }
  bool running = fsm_.stateName()->toString()=="Enabled";
  bool stopping = fsm_.stateName()->toString()=="stopping";
  for(unsigned int i = 0; i < subs_.size(); i++)
    {
      if(subs_[i].alive()==-1000) continue;
      int sl;

      pid_t killedOrNot = waitpid(subs_[i].pid(),&sl,WNOHANG);

      if(killedOrNot==subs_[i].pid()) subs_[i].setStatus((WIFEXITED(sl) != 0 ? 0 : -1));
      else continue;
      pthread_mutex_lock(&pickup_lock_);
      ostringstream ost;
      if(subs_[i].alive()==0) ost << " process exited with status " << WEXITSTATUS(sl);
      else if(WIFSIGNALED(sl)!=0) ost << " process terminated with signal " << WTERMSIG(sl);
      else ost << " process stopped ";
      subs_[i].countdown()=slaveRestartDelaySecs_.value_;
      subs_[i].setReasonForFailed(ost.str());
      spMStates_[i] = edm::event_processor::sInvalid;
      spmStates_[i] = 0;
      ostringstream ost1;
      ost1 << "-E- Slave " << subs_[i].pid() << ost.str();
      localLog(ost1.str());
      if(!autoRestartSlaves_.value_) subs_[i].disconnect();
      pthread_mutex_unlock(&pickup_lock_);
    }
  pthread_mutex_unlock(&stop_lock_);	
  if(stopping) return true; // if in stopping we are done

  if(running)
    {
      // if enabled, this loop will periodically check if dead slaves countdown has expired and restart them
      // this is only active while running, hence, the stop lock is acquired and only released at end of loop
      if(autoRestartSlaves_.value_){
	pthread_mutex_lock(&stop_lock_); //lockout slave killing at stop while you check for restarts
	for(unsigned int i = 0; i < subs_.size(); i++)
	  {
	    if(subs_[i].alive() != 1){
	      if(subs_[i].countdown()-- == 0)
		{
		  pid_t rr = subs_[i].forkNew();
		  if(rr==0)
		    {
		      isChildProcess_=true;
		      try{
			pt::PeerTransport * ptr =
			  pt::getPeerTransportAgent()->getPeerTransport("http","soap",pt::Receiver);
			delete ptr;
		      }
		      catch (pt::exception::PeerTransportNotFound & e ){
			LOG4CPLUS_WARN(getApplicationLogger()," ***Slave Failed to shutdown ptHTTP");
		      }
		      fsm_.disableRcmsStateNotification();
		      fsm_.fireEvent("Stop",this); // need to set state in fsm first to allow stopDone transition
		      fsm_.fireEvent("StopDone",this); // need to set state in fsm first to allow stopDone transition
		      fsm_.fireEvent("Enable",this); // need to set state in fsm first to allow stopDone transition
		      enableMPEPSlave(i);
		      return false; // exit the supervisor loop immediately in the child !!!
		    }
		  else
		    {
		      ostringstream ost1;
		      ost1 << "-I- New Process " << rr << " forked for slot " << i; 
		      localLog(ost1.str());
		    }
		}
	    }
	  }
	pthread_mutex_unlock(&stop_lock_);
      } // finished handling replacement of dead slaves once they've been reaped
    }
  xdata::Serializable *lsid = 0; 
  xdata::Serializable *psid = 0;
  xdata::Serializable *epMAltState = 0; 
  xdata::Serializable *epmAltState = 0;
  
  MsgBuf msg1(0,MSQM_MESSAGE_TYPE_PRG);
  MsgBuf msg2(MAX_MSG_SIZE,MSQS_MESSAGE_TYPE_PRR);
  if(running){  
    try{
      lsid = applicationInfoSpace_->find("lumiSectionIndex");
      psid = applicationInfoSpace_->find("prescaleSetIndex");
      nbProcessed = monitorInfoSpace_->find("nbProcessed");
      nbAccepted  = monitorInfoSpace_->find("nbAccepted");
      epMAltState = monitorInfoSpace_->find("epSPMacroStateInt");
      epmAltState = monitorInfoSpace_->find("epSPMicroStateInt");
    
    }
    catch(xdata::exception::Exception e){
      LOG4CPLUS_INFO(getApplicationLogger(),"could not retrieve some data - " << e.what());    
    }

    try{
      if(nbProcessed !=0 && nbAccepted !=0)
	{
	  xdata::UnsignedInteger32*nbp = ((xdata::UnsignedInteger32*)nbProcessed);
	  xdata::UnsignedInteger32*nba = ((xdata::UnsignedInteger32*)nbAccepted);
	  xdata::UnsignedInteger32*ls  = ((xdata::UnsignedInteger32*)lsid);
	  xdata::UnsignedInteger32*ps  = ((xdata::UnsignedInteger32*)psid);
	  nbp->value_ = 0;
	  nba->value_ = 0;
	  nblive_ = 0;
	  nbdead_ = 0;
	  
	  for(unsigned int i = 0; i < subs_.size(); i++)
	    {
	      if(subs_[i].alive()>0)
		{
		  nblive_++;
		  subs_[i].post(msg1,true);
		  
		    subs_[i].rcvNonBlocking(msg2,true);

		    prg* p = (struct prg*)(msg2->mtext);
		    subs_[i].setParams(p);
		    spMStates_[i] = p->Ms;
		    spmStates_[i] = p->ms;
		    ((xdata::UnsignedInteger32*)nbProcessed)->value_ += p->nbp;
		    ((xdata::UnsignedInteger32*)nbAccepted)->value_  += p->nba;
		    if(p->ls > ls->value_) ls->value_ = p->ls;
		    if(p->ps != ps->value_) ps->value_ = p->ps;
		  }
		else
		  nbdead_++;
	      }
	  }

      }
      catch(evf::Exception &e){
	LOG4CPLUS_INFO(getApplicationLogger(),"could not send/receive msg - " << e.what());    
      }
      catch(std::exception &e){
	LOG4CPLUS_INFO(getApplicationLogger(),"std exception - " << e.what());    
      }
      catch(...){
	LOG4CPLUS_INFO(getApplicationLogger(),"unknown exception ");    
      }
  }
  else{
    for(unsigned int i = 0; i < subs_.size(); i++)
      {
	if(subs_[i].alive()==-1000)
	  {
	    spMStates_[i] = edm::event_processor::sInit;
	    spmStates_[i] = 0;
	  }
      }
  }
  try{
    monitorInfoSpace_->lock();
    monitorInfoSpace_->fireItemGroupChanged(names_,0);
    monitorInfoSpace_->unlock();
  }
  catch(xdata::exception::Exception &e)
    {
      LOG4CPLUS_ERROR(log_, "Exception from fireItemGroupChanged: " << e.what());
      //	localLog(e.what());
    }
  ::sleep(superSleepSec_.value_);	
  return true;
}

void FUEventProcessor::startScalersWorkLoop() throw (evf::Exception)
{
  try {
    if(isChildProcess_ || nbSubProcesses_.value_==0)
      {
	wlScalers_=
	  toolbox::task::getWorkLoopFactory()->getWorkLoop("Scalers",
							   "waiting");
	if (!wlScalers_->isActive()) wlScalers_->activate();
	asScalers_ = toolbox::task::bind(this,&FUEventProcessor::scalers,
					 "Scalers");
      }
    else
      {
	wlScalers_=
	  toolbox::task::getWorkLoopFactory()->getWorkLoop("Summary",
							   "waiting");
	if (!wlScalers_->isActive()) wlScalers_->activate();

	asScalers_ = toolbox::task::bind(this,&FUEventProcessor::summarize,
				       "Summary");
      }
    wlScalers_->submit(asScalers_);
    wlScalersActive_ = true;
  }
  catch (xcept::Exception& e) {
    std::string msg = "Failed to start workloop 'Scalers'.";
    XCEPT_RETHROW(evf::Exception,msg,e);
  }
}

//______________________________________________________________________________
bool FUEventProcessor::scalers(toolbox::task::WorkLoop* wl)
{
  ::sleep(1); //avoid synchronization issues at the start of the event loop

  monitorInfoSpace_->lock();

  if(evtProcessor_)
    {

      edm::event_processor::State st = evtProcessor_->getState();
      monitorInfoSpace_->unlock();

      if(st == edm::event_processor::sRunning)
	{

	  if(!evtProcessor_.getTriggerReport(true)) {

	    wlScalersActive_ = false;

	    return false;
	  }
	}
      else 
	{
	  std::cout << getpid()<< " Scalers workloop, bailing out, not running " << std::endl;
	  wlScalersActive_ = false;
	  return false;
	}
    }
  else
    {
      monitorInfoSpace_->unlock();
      std::cout << getpid()<< " Scalers workloop, bailing out, no evtProcessor " << std::endl;
      wlScalersActive_ = false;
      return false;
    }
  if(isChildProcess_) 
    {
      std::cout << "going to post on control queue from scalers" << std::endl;
      int ret = sq_->post(evtProcessor_.getPackedTriggerReport());
      if(ret!=0)      std::cout << "scalers workloop, error posting to sq_ " << errno << std::endl;
    }
  else
    evtProcessor_.fireScalersUpdate();
  return true;
}

//______________________________________________________________________________
bool FUEventProcessor::summarize(toolbox::task::WorkLoop* wl)
{
  MsgBuf msg(MAX_MSG_SIZE,MSQS_MESSAGE_TYPE_TRR);
  evtProcessor_.resetPackedTriggerReport();
  for (unsigned int i = 0; i < subs_.size(); i++)
    {
      if(subs_[i].alive()>0)
	{

	  int ret = 0;
	  try{
	    ret = subs_[i].rcv(msg,false);
	  }
	  catch(evf::Exception &e)
	    {
	      std::cout << "exception in msgrcv on " << i 
			<< " " << subs_[i].alive() << " " << errno << std::endl;
	      continue;
	      //do nothing special
	    }
	  

	  if(ret==MSQS_MESSAGE_TYPE_TRR) evtProcessor_.sumAndPackTriggerReport(msg);
	  else std::cout << "msgrcv returned error " << errno << std::endl;
	}
    }
  evtProcessor_.updateRollingReport();
  evtProcessor_.fireScalersUpdate();
  if(fsm_.stateName()->toString()!="Enabled"){
    wlScalersActive_ = false;
    return false;
  }
  return true;
}



bool FUEventProcessor::receivingAndMonitor(toolbox::task::WorkLoop *)
{
  MsgBuf msg;
  try{
    sqm_->rcv(msg); //will receive only messages from Master
    switch(msg->mtype)
      {
      case MSQM_MESSAGE_TYPE_MCS:
	{
	  xgi::Input *in = 0;
	  xgi::Output out;
	  evtProcessor_.microState(in,&out);
	  MsgBuf msg1(out.str().size(),MSQS_MESSAGE_TYPE_MCR);
	  strncpy(msg1->mtext,out.str().c_str(),out.str().size());
	  sqm_->post(msg1);
	  break;
	}
      
      case MSQM_MESSAGE_TYPE_PRG:
	{
	  MsgBuf msg1(sizeof(prg),MSQS_MESSAGE_TYPE_PRR);
	  // 	  monitorInfoSpace_->lock();  
	  prg * data = (prg*)msg1->mtext;
	  data->ls=evtProcessor_.lsid_;
	  data->ps=evtProcessor_.psid_;
	  data->nbp=evtProcessor_->totalEvents();
	  data->nba=evtProcessor_->totalEventsPassed();
	  data->Ms=evtProcessor_.epMAltState_.value_;
	  data->ms=evtProcessor_.epmAltState_.value_;

	  //	  monitorInfoSpace_->unlock();  
	  sqm_->post(msg1);
	  break;
	}
      case MSQM_MESSAGE_TYPE_WEB:
	{
	  xgi::Input  *in = 0;
	  xgi::Output out;
	  unsigned int bytesToSend = 0;
	  MsgBuf msg1(NUMERIC_MESSAGE_SIZE,MSQS_MESSAGE_TYPE_WEB);
	  std::string query = msg->mtext;
	  size_t pos = query.find_first_of("&");
	  std::string method;
	  std::string args;
	  if(pos!=std::string::npos)  
	    {
	      method = query.substr(0,pos);
	      args = query.substr(pos+1,query.length()-pos-1);
	    }
	  else
	    method=query;

	  if(method=="Spotlight")
	    {
	      spotlightWebPage(in,&out);
	    }
	  else if(method=="procStat")
	    {
	      procStat(in,&out);
	    }
	  else if(method=="moduleWeb")
	    {
	      internal::MyCgi mycgi;
	      boost::char_separator<char> sep(";");
	      boost::tokenizer<boost::char_separator<char> > tokens(args, sep);
	      for (boost::tokenizer<boost::char_separator<char> >::iterator tok_iter = tokens.begin();
		   tok_iter != tokens.end(); ++tok_iter){
		size_t pos = (*tok_iter).find_first_of("%");
		if(pos != std::string::npos){
		  std::string first  = (*tok_iter).substr(0    ,                        pos);
		  std::string second = (*tok_iter).substr(pos+1, (*tok_iter).length()-pos-1);
		  mycgi.getEnvironment()[first]=second;
		}
	      }
	      moduleWeb(&mycgi,&out);
	    }
	  else if(method=="Default")
	    {
	      defaultWebPage(in,&out);
	    }
	  else 
	    {
	      out << "Error 404!!!!!!!!" << std::endl;
	    }


	  bytesToSend = out.str().size();
	  unsigned int cycle = 0;
	  if(bytesToSend==0)
	    {
	      snprintf(msg1->mtext, NUMERIC_MESSAGE_SIZE, "%d", bytesToSend);
	      sqm_->post(msg1);
	    }
	  while(bytesToSend !=0){
	    unsigned int msgSize = bytesToSend>MAX_PIPE_BUFFER_SIZE ? MAX_PIPE_BUFFER_SIZE : bytesToSend;
	    snprintf(msg1->mtext, NUMERIC_MESSAGE_SIZE, "%d", msgSize);
	    write(anonymousPipe_[PIPE_WRITE],
		  out.str().c_str()+MAX_PIPE_BUFFER_SIZE*cycle,
		  msgSize);
	    sqm_->post(msg1);
	    bytesToSend -= msgSize;
	    cycle++;
	  }
	  break;
	}
      case MSQM_MESSAGE_TYPE_TRP:
	{
	  break;
	}
      }
  }
  catch(evf::Exception &e){std::cout << "exception caught in recevingM: " << e.what() << std::endl;}
  return true;
}

bool FUEventProcessor::enableCommon()
{
  try {    
    if(hasShMem_) attachDqmToShm();
    int sc = 0;
    evtProcessor_->clearCounters();
    if(isRunNumberSetter_)
      evtProcessor_->setRunNumber(runNumber_.value_);
    else
      evtProcessor_->declareRunNumber(runNumber_.value_);
    try{
      ::sleep(1);
      evtProcessor_->runAsync();
      sc = evtProcessor_->statusAsync();
    }
    catch(cms::Exception &e) {
      reasonForFailedState_ = e.explainSelf();
      fsm_.fireFailed(reasonForFailedState_,this);
      localLog(reasonForFailedState_);
      return false;
    }    
    catch(std::exception &e) {
      reasonForFailedState_  = e.what();
      fsm_.fireFailed(reasonForFailedState_,this);
      localLog(reasonForFailedState_);
      return false;
    }
    catch(...) {
      reasonForFailedState_ = "Unknown Exception";
      fsm_.fireFailed(reasonForFailedState_,this);
      localLog(reasonForFailedState_);
      return false;
    }
    if(sc != 0) {
      ostringstream oss;
      oss<<"EventProcessor::runAsync returned status code " << sc;
      reasonForFailedState_ = oss.str();
      fsm_.fireFailed(reasonForFailedState_,this);
      localLog(reasonForFailedState_);
      return false;
    }
  }
  catch (xcept::Exception &e) {
    reasonForFailedState_ = "enabling FAILED: " + (string)e.what();
    fsm_.fireFailed(reasonForFailedState_,this);
    localLog(reasonForFailedState_);
    return false;
  }
  try{
    fsm_.fireEvent("EnableDone",this);
  }
  catch (xcept::Exception &e) {
    std::cout << "exception " << (string)e.what() << std::endl;
    throw;
  }

  return false;
}
  
bool FUEventProcessor::enableClassic()
{
  bool retval = enableCommon();
  while(evtProcessor_->getState()!= edm::event_processor::sRunning){
    LOG4CPLUS_INFO(getApplicationLogger(),"waiting for edm::EventProcessor to start before enabling watchdog");
    ::sleep(1);
  }
  
  //  implementation moved to EPWrapper
  startScalersWorkLoop();
  localLog("-I- Start completed");
  return retval;
}
bool FUEventProcessor::enableMPEPSlave(int ind)
{
  //all this happens only in the child process
  sq_ = new SlaveQueue(ind);
  sqm_ = new SlaveQueue(200+ind);
  startReceivingLoop();
  startReceivingMonitorLoop();
  ::sleep(1);
  evtProcessor_.startMonitoringWorkLoop();
  try{
    //    evtProcessor_.makeServicesOnly();
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
    } 
    catch(...) {
      LOG4CPLUS_ERROR(getApplicationLogger(),"Unknown Exception");
    }
  ML::MLlog4cplus::setAppl(this);      
  }	  
  catch (xcept::Exception &e) {
    reasonForFailedState_ = "enabling FAILED: " + (string)e.what();
    fsm_.fireFailed(reasonForFailedState_,this);
    localLog(reasonForFailedState_);
  }
  bool retval =  enableCommon();
  while(evtProcessor_->getState()!= edm::event_processor::sRunning){
    LOG4CPLUS_INFO(getApplicationLogger(),"waiting for edm::EventProcessor to start before enabling watchdog");
    ::sleep(1);
  }
  startScalersWorkLoop();
  return retval;
}

bool FUEventProcessor::stopClassic()
{
  try {
    LOG4CPLUS_INFO(getApplicationLogger(),"Start stopping :) ...");
    edm::EventProcessor::StatusCode rc = evtProcessor_.stop();
    if(rc != edm::EventProcessor::epTimedOut) 
      fsm_.fireEvent("StopDone",this);
    else
      {
	//	epMState_ = evtProcessor_->currentStateName();
	reasonForFailedState_ = "EventProcessor stop timed out";
	fsm_.fireFailed(reasonForFailedState_,this);
	localLog(reasonForFailedState_);
      }
    if(hasShMem_) detachDqmFromShm();
  }
  catch (xcept::Exception &e) {
    reasonForFailedState_ = "stopping FAILED: " + (string)e.what();
    localLog(reasonForFailedState_);
    fsm_.fireFailed(reasonForFailedState_,this);
  }
  LOG4CPLUS_INFO(getApplicationLogger(),"Finished stopping!");
  localLog("-I- Stop completed");
  return false;
}

void FUEventProcessor::stopSlavesAndAcknowledge()
{
  MsgBuf msg(0,MSQM_MESSAGE_TYPE_STOP);
  MsgBuf msg1(MAX_MSG_SIZE,MSQS_MESSAGE_TYPE_STOP);

  for(unsigned int i = 0; i < subs_.size(); i++)
    {
      pthread_mutex_lock(&stop_lock_);
      if(subs_[i].alive()>0)subs_[i].post(msg,false);

      if(subs_[i].alive()<=0) 
	{
	  pthread_mutex_unlock(&stop_lock_);
	  continue;
	}
      try{
	subs_[i].rcv(msg1,false);
      }
      catch(evf::Exception &e){
	ostringstream ost;
	ost << "failed to get STOP - errno ->" << errno << " " << e.what(); 
	reasonForFailedState_ = ost.str();
	LOG4CPLUS_WARN(getApplicationLogger(),reasonForFailedState_);
	fsm_.fireFailed(reasonForFailedState_,this);
	localLog(reasonForFailedState_);
	break;
      }
      pthread_mutex_unlock(&stop_lock_);
      if(msg1->mtype==MSQS_MESSAGE_TYPE_STOP)
	while(subs_[i].alive()>0) ::usleep(10000);
      subs_[i].disconnect();
    }
  //  subs_.clear();

}

void FUEventProcessor::microState(xgi::Input *in,xgi::Output *out)
{
  string urn = getApplicationDescriptor()->getURN();
  try{
    evtProcessor_.stateNameFromIndex(0);
    evtProcessor_.moduleNameFromIndex(0);
  if(isChildProcess_) {std::cout << "microstate called for child! bail out" << std::endl; return;}
  *out << "<tr><td>" << fsm_.stateName()->toString() 
       << "</td><td>"<< (isChildProcess_ ? "S" : "M") <<"</td><td>" << nblive_ << "</td><td>"
       << nbdead_ << "</td><td><a href=\"/" << urn << "/procStat\">" << getpid() <<"</a></td>";
  evtProcessor_.microState(in,out);
  *out << "</tr>";
  if(nbSubProcesses_.value_!=0 && !isChildProcess_) 
    {
      pthread_mutex_lock(&start_lock_);
      for(unsigned int i = 0; i < subs_.size(); i++)
	{
	  try{
	    if(subs_[i].alive()>0)
	      {
		*out << "<tr><td  bgcolor=\"#00FF00\" id=\"a"
		     << i << "\">""Alive</td><td>S</td><td>"
		     << subs_[i].queueId() << "<td>" 
		     << subs_[i].queueStatus()<< "/"
		     << subs_[i].queueOccupancy() << "/"
		     << subs_[i].queuePidOfLastSend() << "/"
		     << subs_[i].queuePidOfLastReceive() 
		     << "</td><td><a id=\"p"<< i << "\" href=\"SubWeb?process=" 
		     << subs_[i].pid() << "&method=procStat\">" 
		     << subs_[i].pid()<<"</a></td>" //<< msg->mtext;
		     << "<td>" << evtProcessor_.stateNameFromIndex(subs_[i].params().Ms) << "</td><td>" 
		     << evtProcessor_.moduleNameFromIndex(subs_[i].params().ms) << "</td><td>" 
		     << subs_[i].params().nba << "/" << subs_[i].params().nbp 
		     << " (" << float(subs_[i].params().nba)/float(subs_[i].params().nbp)*100. <<"%)" 
		     << "</td><td>" << subs_[i].params().ls << "/" << subs_[i].params().ls 
		     << "</td><td>" << subs_[i].params().ps << "</td>";
	      }
	    else 
	      {
		pthread_mutex_lock(&pickup_lock_);
		*out << "<tr><td id=\"a"<< i << "\" ";
		if(subs_[i].alive()==-1000)
		  *out << " bgcolor=\"#bbaabb\">NotInitialized";
		else
		  *out << (subs_[i].alive()==0 ? ">Done" : " bgcolor=\"#FF0000\">Dead");
		*out << "</td><td>S</td><td>"<< subs_[i].queueId() << "<td>" 
		     << subs_[i].queueStatus() << "/"
		     << subs_[i].queueOccupancy() << "/"
		     << subs_[i].queuePidOfLastSend() << "/"
		     << subs_[i].queuePidOfLastReceive() 
		     << "</td><td id=\"p"<< i << "\">"
		     <<subs_[i].pid()<<"</td><td colspan=\"5\">" << subs_[i].reasonForFailed();
		if(subs_[i].alive()!=0 && subs_[i].alive()!=-1000) 
		  {
		    if(autoRestartSlaves_) *out << " will restart in " << subs_[i].countdown() << " s";
		    else *out << " autoRestart is disabled ";
		  }
		*out << "</td>";
		pthread_mutex_unlock(&pickup_lock_);
	      }
	    *out << "</tr>";
	  }
	  catch(evf::Exception &e){
	    *out << "<tr><td id=\"a"<< i << "\" " 
		 <<"bgcolor=\"#FFFF00\">NotResp</td><td>S</td><td>"<< subs_[i].queueId() << "<td>" 
		 << subs_[i].queueStatus() << "/"
		 << subs_[i].queueOccupancy() << "/"
		 << subs_[i].queuePidOfLastSend() << "/"
		 << subs_[i].queuePidOfLastReceive() 
		 << "</td><td id=\"p"<< i << "\">"
		 <<subs_[i].pid()<<"</td>";
	  }
	}
      pthread_mutex_unlock(&start_lock_); 
    }
  }
      catch(evf::Exception &e)
      {
	LOG4CPLUS_INFO(getApplicationLogger(),"evf::Exception caught in microstate - " << e.what());    
      }
    catch(cms::Exception &e)
      {
	LOG4CPLUS_INFO(getApplicationLogger(),"cms::Exception caught in microstate - " << e.what());    
      }
    catch(std::exception &e)
      {
	LOG4CPLUS_INFO(getApplicationLogger(),"std::Exception caught in microstate - " << e.what());    
      }
    catch(...)
      {
	LOG4CPLUS_INFO(getApplicationLogger(),"unknown exception caught in microstate - ");    
      }

}


void FUEventProcessor::updater(xgi::Input *in,xgi::Output *out)
{
  *out << updaterStatic_;
  *out << "<div id=\"loads\">"; 
  utils::uptime(out);
  *out << "</div>";
  *out << "<div id=\"st\">" << fsm_.stateName()->toString() << "</div>"
       << "<div id=\"ru\">" << runNumber_.toString() << "</div>"
       << "<div id=\"nsl\">" << nbSubProcesses_.value_ << "</div>"
       << "<div id=\"cl\">" << getApplicationDescriptor()->getClassName() 
       << (nbSubProcesses_.value_ > 0 ? "MP " : " ") << "</div>"
       << "<div id=\"in\">" << getApplicationDescriptor()->getInstance() << "</div>";
  if(fsm_.stateName()->toString() != "Halted" && fsm_.stateName()->toString() != "halting")
    *out << "<div id=\"hlt\"><a href=\"" << configString_.toString() << "\">HLT Config</a></div>" << endl;
  else
    *out << "<div id=\"hlt\">Not yet...</div>" << endl;  
  *out << "<div id=\"sq\">" << squidPresent_.toString() << "</div>"
       << "<div id=\"vwl\">" << (supervising_ ? "Active" : "not initialized") << "</div>"
       << "<div id=\"mwl\">" << evtProcessor_.wlMonitoring() << "</div>";
  if(nbProcessed != 0 && nbAccepted != 0)
    {
      *out << "<div id=\"tt\">" << ((xdata::UnsignedInteger32*)nbProcessed)->value_ << "</div>"
	   << "<div id=\"ac\">" << ((xdata::UnsignedInteger32*)nbAccepted)->value_ << "</div>";
    }
  else
    {
      *out << "<div id=\"tt\">" << 0 << "</div>"
	   << "<div id=\"ac\">" << 0 << "</div>";
    }

  *out<< "<div id=\"swl\">" << (wlScalersActive_ ? "active" : "inactive") << "</div>";
  if(evtProcessor_)
    *out << "<div id=\"ll\">" << evtProcessor_.lastLumi().ls
	 << "," << evtProcessor_.lastLumi().proc << "," << evtProcessor_.lastLumi().acc << "</div>";
    
  *out << "<div id=\"lg\">";
  for(unsigned int i = logRingIndex_; i<logRingSize_; i++)
    *out << logRing_[i] << std::endl;
  if(logWrap_)
    for(unsigned int i = 0; i<logRingIndex_; i++)
      *out << logRing_[i] << std::endl;
  *out << "</div>" << std::endl;
}
void FUEventProcessor::procStat(xgi::Input *in, xgi::Output *out)
{
  evf::utils::procStat(out);
}
void FUEventProcessor::sendMessageOverMonitorQueue(MsgBuf &buf)
{
  sqm_->post(buf);
}


XDAQ_INSTANTIATOR_IMPL(evf::FUEventProcessor)
