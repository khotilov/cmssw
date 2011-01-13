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
  , exitOnError_(true)
  , reasonForFailedState_()
  , squidnet_(3128,"http://localhost:8000/RELEASE-NOTES.txt")
  , logRing_(logRingSize_)
  , logRingIndex_(logRingSize_)
  , logWrap_(false)
  , nbSubProcesses_(0)
  , nblive_(0)
  , nbdead_(0)
  , nbTotalDQM_(0)
  , wlReceiving_(0)
  , asReceiveMsgAndExecute_(0)
  , receiving_(false) 
  , wlReceivingMonitor_(0)
  , asReceiveMsgAndRead_(0)
  , receivingM_(false)
  , myProcess_(0)
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
  , scalersUpdates_(0)
  , wlSummarize_(0)
  , asSummarize_(0)
  , wlSummarizeActive_(false)
  , superSleepSec_(1)
  , iDieUrl_("none")
  , vulture_(0)
  , vp_(0)
{
  using namespace utils;

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
  ispace->fireItemAvailable("iDieUrl",              &iDieUrl_                     );
  
  // Add infospace listeners for exporting data values
  getApplicationInfoSpace()->addItemChangedListener("parameterSet",        this);
  getApplicationInfoSpace()->addItemChangedListener("outputEnabled",       this);

  // findRcmsStateListener
  fsm_.findRcmsStateListener();
  
  // initialize monitoring infospace

  std::stringstream oss2;
  oss2<<"urn:xdaq-monitorable-"<<class_.toString();
  std::string monInfoSpaceName=oss2.str();
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
  std::string monInfoSpaceName2=oss3.str();
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

  typedef std::set<xdaq::ApplicationDescriptor*> AppDescSet_t;
  typedef AppDescSet_t::iterator                 AppDescIter_t;
  
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

  makeStaticInfo();
  startSupervisorLoop();  

  if(vulture_==0) vulture_ = new Vulture(true);

  ////////////////////////////////  

  AppDescSet_t setOfiDie=
    getApplicationContext()->getDefaultZone()->
    getApplicationDescriptors("evf::iDie");
  
  for (AppDescIter_t it=setOfiDie.begin();it!=setOfiDie.end();++it)
    if ((*it)->getInstance()==0) // there has to be only one instance of iDie
      iDieUrl_ = (*it)->getContextDescriptor()->getURL() + "/" + (*it)->getURN();
}
//___________here ends the *huge* constructor___________________________________


//______________________________________________________________________________
FUEventProcessor::~FUEventProcessor()
{
  // no longer needed since the wrapper is a member of the class and one can rely on 
  // implicit destructor - to be revised - at any rate the most common exit path is via "kill"...
  //  if (evtProcessor_) delete evtProcessor_;
  if(vulture_ != 0) delete vulture_;
}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
bool FUEventProcessor::configuring(toolbox::task::WorkLoop* wl)
{
//   std::cout << "values " << ((nbSubProcesses_.value_!=0) ? 0x10 : 0) << " "
// 	    << ((instance_.value_==0) ? 0x8 : 0) << " "
// 	    << (hasServiceWebRegistry_.value_ ? 0x4 : 0) << " "
// 	    << (hasModuleWebRegistry_.value_ ? 0x2 : 0) << " "
// 	    << (hasPrescaleService_.value_ ? 0x1 : 0) <<std::endl;
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
    reasonForFailedState_ = "configuring FAILED: " + (std::string)e.what();
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

  if(vulture_!=0 && vp_ == 0) vp_ = vulture_->makeProcess();
  return false;
}




//______________________________________________________________________________
bool FUEventProcessor::enabling(toolbox::task::WorkLoop* wl)
{
  nbTotalDQM_ = 0;
  scalersUpdates_ = 0;
//   std::cout << "values " << ((nbSubProcesses_.value_!=0) ? 0x10 : 0) << " "
// 	    << ((instance_.value_==0) ? 0x8 : 0) << " "
// 	    << (hasServiceWebRegistry_.value_ ? 0x4 : 0) << " "
// 	    << (hasModuleWebRegistry_.value_ ? 0x2 : 0) << " "
// 	    << (hasPrescaleService_.value_ ? 0x1 : 0) <<std::endl;
  unsigned short smap 
    = ((nbSubProcesses_.value_!=0) ? 0x10 : 0)
    + ((instance_.value_==0) ? 0x8 : 0)
    + (hasServiceWebRegistry_.value_ ? 0x4 : 0) 
    + (hasModuleWebRegistry_.value_ ? 0x2 : 0) 
    + (hasPrescaleService_.value_ ? 0x1 : 0);

  LOG4CPLUS_INFO(getApplicationLogger(),"Start enabling...");
  if(!epInitialized_) evtProcessor_.forceInitEventProcessorMaybe();

  std::string cfg = configString_.toString(); evtProcessor_.init(smap,cfg);
  configuration_ = evtProcessor_.configuration(); // get it again after init has been carried out...
  evtProcessor_.resetLumiSectionReferenceIndex();
  //classic appl will return here 
  if(nbSubProcesses_.value_==0) return enableClassic();
  //protect manipulation of subprocess array
  pthread_mutex_lock(&start_lock_);
  subs_.clear();
  subs_.resize(nbSubProcesses_.value_);
  pthread_mutex_unlock(&start_lock_);
  pid_t retval = -1;
  for(unsigned int i=0; i<nbSubProcesses_.value_; i++)
    {
      subs_[i]=SubProcess(i,retval); //this will replace all the scattered variables
      retval = subs_[i].forkNew();
      if(retval==0)
	{
	  myProcess_ = &subs_[i];
	  int retval = pthread_mutex_destroy(&stop_lock_);
	  if(retval != 0) perror("error");
	  retval = pthread_mutex_init(&stop_lock_,0);
	  if(retval != 0) perror("error");
	  try{
	    pt::PeerTransport * ptr =
	      pt::getPeerTransportAgent()->getPeerTransport("http","soap",pt::Receiver);
	    delete ptr;
	  }
	  catch (pt::exception::PeerTransportNotFound & e ){
	    LOG4CPLUS_WARN(getApplicationLogger()," ***Slave Failed to shutdown ptHTTP");
	  }
	  fsm_.disableRcmsStateNotification();
	  return enableMPEPSlave();
	  // the loop is broken in the child 
	}
    }
  
  startSummarizeWorkLoop();
  vp_ = vulture_->start(iDieUrl_.value_,runNumber_.value_);

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
  vulture_->stop();
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
    reasonForFailedState_ = "halting FAILED: " + (std::string)e.what();
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
    std::string item = dynamic_cast<xdata::ItemChangedEvent&>(e).itemName();
    
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
      if(subs_[i].alive() != 1){
	*out << "ERROR 405 : Process " << pid << " Not Alive !" << std::endl;
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

  std::string urn = getApplicationDescriptor()->getURN();

  *out << "<!-- base href=\"/" <<  urn
       << "\"> -->" << std::endl;
  *out << "<html>"                                                   << std::endl;
  *out << "<head>"                                                   << std::endl;
  *out << "<link type=\"text/css\" rel=\"stylesheet\"";
  *out << " href=\"/evf/html/styles.css\"/>"                   << std::endl;
  *out << "<title>" << getApplicationDescriptor()->getClassName() 
       << getApplicationDescriptor()->getInstance() 
       << " MAIN</title>"     << std::endl;
  *out << "</head>"                                                  << std::endl;
  *out << "<body>"                                                   << std::endl;
  *out << "<table border=\"0\" width=\"100%\">"                      << std::endl;
  *out << "<tr>"                                                     << std::endl;
  *out << "  <td align=\"left\">"                                    << std::endl;
  *out << "    <img"                                                 << std::endl;
  *out << "     align=\"middle\""                                    << std::endl;
  *out << "     src=\"/evf/images/spoticon.jpg\""			     << std::endl;
  *out << "     alt=\"main\""                                        << std::endl;
  *out << "     width=\"64\""                                        << std::endl;
  *out << "     height=\"64\""                                       << std::endl;
  *out << "     border=\"\"/>"                                       << std::endl;
  *out << "    <b>"                                                  << std::endl;
  *out << getApplicationDescriptor()->getClassName() 
       << getApplicationDescriptor()->getInstance()                  << std::endl;
  *out << "      " << fsm_.stateName()->toString()                   << std::endl;
  *out << "    </b>"                                                 << std::endl;
  *out << "  </td>"                                                  << std::endl;
  *out << "  <td width=\"32\">"                                      << std::endl;
  *out << "    <a href=\"/urn:xdaq-application:lid=3\">"             << std::endl;
  *out << "      <img"                                               << std::endl;
  *out << "       align=\"middle\""                                  << std::endl;
  *out << "       src=\"/hyperdaq/images/HyperDAQ.jpg\""             << std::endl;
  *out << "       alt=\"HyperDAQ\""                                  << std::endl;
  *out << "       width=\"32\""                                      << std::endl;
  *out << "       height=\"32\""                                     << std::endl;
  *out << "       border=\"\"/>"                                     << std::endl;
  *out << "    </a>"                                                 << std::endl;
  *out << "  </td>"                                                  << std::endl;
  *out << "  <td width=\"32\">"                                      << std::endl;
  *out << "  </td>"                                                  << std::endl;
  *out << "  <td width=\"32\">"                                      << std::endl;
  *out << "    <a href=\"/" << urn << "/\">"                         << std::endl;
  *out << "      <img"                                               << std::endl;
  *out << "       align=\"middle\""                                  << std::endl;
  *out << "       src=\"/evf/images/epicon.jpg\""		     << std::endl;
  *out << "       alt=\"main\""                                      << std::endl;
  *out << "       width=\"32\""                                      << std::endl;
  *out << "       height=\"32\""                                     << std::endl;
  *out << "       border=\"\"/>"                                     << std::endl;
  *out << "    </a>"                                                 << std::endl;
  *out << "  </td>"                                                  << std::endl;
  *out << "</tr>"                                                    << std::endl;
  *out << "</table>"                                                 << std::endl;

  *out << "<hr/>"                                                    << std::endl;
  
  std::ostringstream ost;
  if(myProcess_) 
    ost << "/SubWeb?process=" << getpid() << "&method=moduleWeb&";
  else
    ost << "/moduleWeb?";
  urn += ost.str();
  if(evtProcessor_ && (myProcess_ || nbSubProcesses_.value_==0))
    evtProcessor_.taskWebPage(in,out,urn);
  else if(evtProcessor_)
    evtProcessor_.summaryWebPage(in,out,urn);
  else
    *out << "<td>HLT Unconfigured</td>" << std::endl;
  *out << "</table>"                                                 << std::endl;
  
  *out << "<br><textarea rows=" << 10 << " cols=80 scroll=yes>"      << std::endl;
  *out << configuration_                                             << std::endl;
  *out << "</textarea><P>"                                           << std::endl;
  
  *out << "</body>"                                                  << std::endl;
  *out << "</html>"                                                  << std::endl;


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
  std::string errmsg;
  bool success = false;
  try {
    edm::ServiceRegistry::Operate operate(evtProcessor_->getToken());
    if(edm::Service<FUShmDQMOutputService>().isAvailable())
      success = edm::Service<FUShmDQMOutputService>()->attachToShm();
    if (!success) errmsg = "Failed to attach DQM service to shared memory";
  }
  catch (cms::Exception& e) {
    errmsg = "Failed to attach DQM service to shared memory: " + (std::string)e.what();
  }
  if (!errmsg.empty()) XCEPT_RAISE(evf::Exception,errmsg);
}



void FUEventProcessor::detachDqmFromShm() throw (evf::Exception)
{
  std::string errmsg;
  bool success = false;
  try {
    edm::ServiceRegistry::Operate operate(evtProcessor_->getToken());
    if(edm::Service<FUShmDQMOutputService>().isAvailable())
      success = edm::Service<FUShmDQMOutputService>()->detachFromShm();
    if (!success) errmsg = "Failed to detach DQM service from shared memory";
  }
  catch (cms::Exception& e) {
    errmsg = "Failed to detach DQM service from shared memory: " + (std::string)e.what();
  }
  if (!errmsg.empty()) XCEPT_RAISE(evf::Exception,errmsg);
}


std::string FUEventProcessor::logsAsString()
{
  std::ostringstream oss;
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
  
void FUEventProcessor::localLog(std::string m)
{
  timeval tv;

  gettimeofday(&tv,0);
  tm *uptm = localtime(&tv.tv_sec);
  char datestring[256];
  strftime(datestring, sizeof(datestring),"%c", uptm);

  if(logRingIndex_ == 0){logWrap_ = true; logRingIndex_ = logRingSize_;}
  logRingIndex_--;
  std::ostringstream timestamp;
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
    myProcess_->rcvSlave(msg,false); //will receive only messages from Master
    if(msg->mtype==MSQM_MESSAGE_TYPE_STOP)
      {
	pthread_mutex_lock(&stop_lock_);
	fsm_.fireEvent("Stop",this); // need to set state in fsm first to allow stopDone transition
	pthread_mutex_unlock(&stop_lock_);
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
	stopClassic(); // call the normal sequence of stopping - as this is allowed to fail provisions must be made ...@@@EM
	MsgBuf msg1(0,MSQS_MESSAGE_TYPE_STOP);
	myProcess_->postSlave(msg1,false);
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
      std::ostringstream ost;
      if(subs_[i].alive()==0) ost << " process exited with status " << WEXITSTATUS(sl);
      else if(WIFSIGNALED(sl)!=0) ost << " process terminated with signal " << WTERMSIG(sl);
      else ost << " process stopped ";
      subs_[i].countdown()=slaveRestartDelaySecs_.value_;
      subs_[i].setReasonForFailed(ost.str());
      spMStates_[i] = evtProcessor_.notstarted_state_code();
      spmStates_[i] = 0;
      std::ostringstream ost1;
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
		      myProcess_=&subs_[i];
		      scalersUpdates_ = 0;
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
		      try{
			xdata::Serializable *lsid = applicationInfoSpace_->find("lumiSectionIndex");
			if(lsid) {
			  ((xdata::UnsignedInteger32*)(lsid))->value_--; // need to reset to value before end of ls in which process died
			}
		      }
		      catch(...){
			std::cout << "trouble with lsindex during restart" << std::endl;
		      }
		      try{
			xdata::Serializable *lstb = applicationInfoSpace_->find("lsToBeRecovered");
			if(lstb) {
			  ((xdata::Boolean*)(lstb))->value_ = false; // do not issue eol/bol for all Ls when restarting
			}
		      }
		      catch(...){
			std::cout << "trouble with resetting flag for eol recovery " << std::endl;
		      }

		      evtProcessor_.adjustLsIndexForRestart();
		      evtProcessor_.resetPackedTriggerReport();
		      enableMPEPSlave();
		      return false; // exit the supervisor loop immediately in the child !!!
		    }
		  else
		    {
		      std::ostringstream ost1;
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
  xdata::Serializable *dqmp = 0;
  xdata::UnsignedInteger32 *dqm = 0;


  
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
      dqmp = applicationInfoSpace_-> find("nbDqmUpdates");      
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
	  if(dqmp!=0)
	    dqm = (xdata::UnsignedInteger32*)dqmp;
	  if(dqm) dqm->value_ = 0;
	  nbTotalDQM_ = 0;
	  nbp->value_ = 0;
	  nba->value_ = 0;
	  nblive_ = 0;
	  nbdead_ = 0;
	  scalersUpdates_ = 0;

	  for(unsigned int i = 0; i < subs_.size(); i++)
	    {
	      if(subs_[i].alive()>0)
		{
		  nblive_++;
		  try{
		    subs_[i].post(msg1,true);
		    
		    unsigned long retval = subs_[i].rcvNonBlocking(msg2,true);
		    if(retval == (unsigned long) msg2->mtype){
		      prg* p = (struct prg*)(msg2->mtext);
		      subs_[i].setParams(p);
		      spMStates_[i] = p->Ms;
		      spmStates_[i] = p->ms;
		      if(!subs_[i].inInconsistentState() && 
			 (p->Ms == edm::event_processor::sError 
			  || p->Ms == edm::event_processor::sInvalid
			  || p->Ms == edm::event_processor::sStopping))
			{
			  std::ostringstream ost;
			  ost << "edm::eventprocessor slot " << i << " process id " 
			      << subs_[i].pid() << " not in Running state : Mstate=" 
			      << evtProcessor_.stateNameFromIndex(p->Ms) << " mstate="
			      << evtProcessor_.moduleNameFromIndex(p->ms) 
			      << " - Look into possible error messages from HLT process";
			  XCEPT_DECLARE(evf::Exception,
					sentinelException, ost.str());
			  notifyQualified("error",sentinelException);
			  subs_[i].setReportedInconsistent();
			}
		      ((xdata::UnsignedInteger32*)nbProcessed)->value_ += p->nbp;
		      ((xdata::UnsignedInteger32*)nbAccepted)->value_  += p->nba;
		      if(dqm)dqm->value_ += p->dqm;
		      nbTotalDQM_ +=  p->dqm;
		      scalersUpdates_ += p->trp;
		      if(p->ls > ls->value_) ls->value_ = p->ls;
		      if(p->ps != ps->value_) ps->value_ = p->ps;
		    }
		  } 
		  catch(evf::Exception &e){
		    LOG4CPLUS_INFO(getApplicationLogger(),
				   "could not send/receive msg on slot " 
				   << i << " - " << e.what());    
		  }
		    
		}
	      else
		nbdead_++;
	    }
	}
      
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
    wlScalers_=
      toolbox::task::getWorkLoopFactory()->getWorkLoop("Scalers",
						       "waiting");
    if (!wlScalers_->isActive()) wlScalers_->activate();
    asScalers_ = toolbox::task::bind(this,&FUEventProcessor::scalers,
				     "Scalers");
    
  wlScalers_->submit(asScalers_);
  wlScalersActive_ = true;
  }
  catch (xcept::Exception& e) {
    std::string msg = "Failed to start workloop 'Scalers'.";
    XCEPT_RETHROW(evf::Exception,msg,e);
  }
}

//______________________________________________________________________________

void FUEventProcessor::startSummarizeWorkLoop() throw (evf::Exception)
{
  try {
    wlSummarize_=
      toolbox::task::getWorkLoopFactory()->getWorkLoop("Summary",
						       "waiting");
    if (!wlSummarize_->isActive()) wlSummarize_->activate();
    
    asSummarize_ = toolbox::task::bind(this,&FUEventProcessor::summarize,
				       "Summary");

    wlSummarize_->submit(asSummarize_);
    wlSummarizeActive_ = true;
  }
  catch (xcept::Exception& e) {
    std::string msg = "Failed to start workloop 'Summarize'.";
    XCEPT_RETHROW(evf::Exception,msg,e);
  }
}

//______________________________________________________________________________

bool FUEventProcessor::scalers(toolbox::task::WorkLoop* wl)
{
  if(evtProcessor_)
    {
      if(!evtProcessor_.getTriggerReport(true)) {
	wlScalersActive_ = false;
	return false;
      }
    }
  else
    {
      std::cout << getpid()<< " Scalers workloop, bailing out, no evtProcessor " << std::endl;
      wlScalersActive_ = false;
      return false;
    }
  if(myProcess_) 
    {
      //      std::cout << getpid() << "going to post on control queue from scalers" << std::endl;
      int ret = myProcess_->postSlave(evtProcessor_.getPackedTriggerReport(),false);
      if(ret!=0)      std::cout << "scalers workloop, error posting to sqs_ " << errno << std::endl;
      scalersUpdates_++;
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
  bool atLeastOneProcessUpdatedSuccessfully = false;
  int msgCount = 0;
  for (unsigned int i = 0; i < subs_.size(); i++)
    {
      if(subs_[i].alive()>0)
	{

	  int ret = 0;
	  try{
	    ret = subs_[i].rcv(msg,false);
	    msgCount++;
	  }
	  catch(evf::Exception &e)
	    {
	      std::cout << "exception in msgrcv on " << i 
			<< " " << subs_[i].alive() << " " << strerror(errno) << std::endl;
	      continue;
	      //do nothing special
	    }
	  

	  if(ret==MSQS_MESSAGE_TYPE_TRR) {
	    TriggerReportStatic *trp = (TriggerReportStatic *)msg->mtext;
	    if(trp->lumiSection > evtProcessor_.getLumiSectionReferenceIndex()){
	      std::cout << "postpone handling of msg from slot " << i << " with Ls " <<  trp->lumiSection
			<< " should be " << evtProcessor_.getLumiSectionReferenceIndex() << std::endl;
	      subs_[i].post(msg,false);
	    }else{
	      atLeastOneProcessUpdatedSuccessfully = true;
	      evtProcessor_.sumAndPackTriggerReport(msg);
	    }
	  }
	  else std::cout << "msgrcv returned error " << errno << std::endl;
	}
    }
  if(atLeastOneProcessUpdatedSuccessfully){
    evtProcessor_.updateRollingReport();
    evtProcessor_.fireScalersUpdate();
  }
  else{
    LOG4CPLUS_WARN(getApplicationLogger(),"Summarize loop: no process updated successfully ");          
    if(msgCount==0) evtProcessor_.withdrawLumiSectionIncrement();
  }
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
    myProcess_->rcvSlave(msg,true); //will receive only messages from Master
    switch(msg->mtype)
      {
      case MSQM_MESSAGE_TYPE_MCS:
	{
	  xgi::Input *in = 0;
	  xgi::Output out;
	  evtProcessor_.microState(in,&out);
	  MsgBuf msg1(out.str().size(),MSQS_MESSAGE_TYPE_MCR);
	  strncpy(msg1->mtext,out.str().c_str(),out.str().size());
	  myProcess_->postSlave(msg1,true);
	  break;
	}
      
      case MSQM_MESSAGE_TYPE_PRG:
	{
	  xdata::Serializable *dqmp = 0;
	  xdata::UnsignedInteger32 *dqm = 0;
	  try{
	    dqmp = applicationInfoSpace_-> find("nbDqmUpdates");
	  }  catch(xdata::exception::Exception e){}
	  if(dqmp!=0)
	    dqm = (xdata::UnsignedInteger32*)dqmp;
	  MsgBuf msg1(sizeof(prg),MSQS_MESSAGE_TYPE_PRR);
	  // 	  monitorInfoSpace_->lock();  
	  prg * data           = (prg*)msg1->mtext;
	  data->ls             = evtProcessor_.lsid_;
	  data->ps             = evtProcessor_.psid_;
	  data->nbp            = evtProcessor_->totalEvents();
	  data->nba            = evtProcessor_->totalEventsPassed();
	  data->Ms             = evtProcessor_.epMAltState_.value_;
	  data->ms             = evtProcessor_.epmAltState_.value_;
	  if(dqm) data->dqm    = dqm->value_; else data->dqm = 0;
	  data->trp            = scalersUpdates_;
	  //	  monitorInfoSpace_->unlock();  
	  myProcess_->postSlave(msg1,true);
	  if(exitOnError_.value_)
	  { 
	    // after each monitoring cycle check if we are in inconsistent state and exit if configured to do so  
	    //	    std::cout << getpid() << "receivingAndMonitor: trying to acquire stop lock " << std::endl;
	    int retval = pthread_mutex_lock(&stop_lock_);
	    if(retval != 0) perror("error");
	    //	    std::cout << getpid() << " stop lock acquired" << std::endl;
	    bool running = fsm_.stateName()->toString()=="Enabled";
	    if(!running) pthread_mutex_unlock(&stop_lock_);
	    else if(data->Ms == edm::event_processor::sStopping || data->Ms == edm::event_processor::sError) 
	      {::sleep(5); exit(-1); /* no need to unlock mutex after exit :-)*/}
	    pthread_mutex_unlock(&stop_lock_);
	    
	  }
	  //	  scalersUpdates_++;
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
	      myProcess_->postSlave(msg1,true);
	    }
	  while(bytesToSend !=0){
	    unsigned int msgSize = bytesToSend>MAX_PIPE_BUFFER_SIZE ? MAX_PIPE_BUFFER_SIZE : bytesToSend;
	    snprintf(msg1->mtext, NUMERIC_MESSAGE_SIZE, "%d", msgSize);
	    write(anonymousPipe_[PIPE_WRITE],
		  out.str().c_str()+MAX_PIPE_BUFFER_SIZE*cycle,
		  msgSize);
	    myProcess_->postSlave(msg1,true);
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
      std::ostringstream oss;
      oss<<"EventProcessor::runAsync returned status code " << sc;
      reasonForFailedState_ = oss.str();
      fsm_.fireFailed(reasonForFailedState_,this);
      localLog(reasonForFailedState_);
      return false;
    }
  }
  catch (xcept::Exception &e) {
    reasonForFailedState_ = "enabling FAILED: " + (std::string)e.what();
    fsm_.fireFailed(reasonForFailedState_,this);
    localLog(reasonForFailedState_);
    return false;
  }
  try{
    fsm_.fireEvent("EnableDone",this);
  }
  catch (xcept::Exception &e) {
    std::cout << "exception " << (std::string)e.what() << std::endl;
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
bool FUEventProcessor::enableMPEPSlave()
{
  //all this happens only in the child process
  startReceivingLoop();
  startReceivingMonitorLoop();
  evtProcessor_.resetWaiting();
  startScalersWorkLoop();
  while(!evtProcessor_.isWaitingForLs())
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
    reasonForFailedState_ = "enabling FAILED: " + (std::string)e.what();
    fsm_.fireFailed(reasonForFailedState_,this);
    localLog(reasonForFailedState_);
  }
  bool retval =  enableCommon();
  //  while(evtProcessor_->getState()!= edm::event_processor::sRunning){
  //    LOG4CPLUS_INFO(getApplicationLogger(),"waiting for edm::EventProcessor to start before enabling watchdog");
  //    ::sleep(1);
  //  }
  return retval;
}

bool FUEventProcessor::stopClassic()
{
  try {
    LOG4CPLUS_INFO(getApplicationLogger(),"Start stopping :) ...");
    edm::EventProcessor::StatusCode rc = evtProcessor_.stop();
    if(rc == edm::EventProcessor::epSuccess) 
      fsm_.fireEvent("StopDone",this);
    else
      {
	//	epMState_ = evtProcessor_->currentStateName();
	if(rc == edm::EventProcessor::epTimedOut)
	  reasonForFailedState_ = "EventProcessor stop timed out";
	else
	  reasonForFailedState_ = "EventProcessor did not receive STOP event";
	fsm_.fireFailed(reasonForFailedState_,this);
	localLog(reasonForFailedState_);
      }
    if(hasShMem_) detachDqmFromShm();
  }
  catch (xcept::Exception &e) {
    reasonForFailedState_ = "stopping FAILED: " + (std::string)e.what();
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
	std::ostringstream ost;
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
  std::string urn = getApplicationDescriptor()->getURN();
  try{
    evtProcessor_.stateNameFromIndex(0);
    evtProcessor_.moduleNameFromIndex(0);
  if(myProcess_) {std::cout << "microstate called for child! bail out" << std::endl; return;}
  *out << "<tr><td>" << fsm_.stateName()->toString() 
       << "</td><td>"<< (myProcess_ ? "S" : "M") <<"</td><td>" << nblive_ << "</td><td>"
       << nbdead_ << "</td><td><a href=\"/" << urn << "/procStat\">" << getpid() <<"</a></td>";
  evtProcessor_.microState(in,out);
  *out << "<td>" << nbTotalDQM_ 
       << "</td><td>" << evtProcessor_.getScalersUpdates() << "</td></tr>";
  if(nbSubProcesses_.value_!=0 && !myProcess_) 
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
		     << "</td><td>" << subs_[i].params().ls  << "/" << subs_[i].params().ls 
		     << "</td><td>" << subs_[i].params().ps 
		     << "</td><td>" << subs_[i].params().dqm 
		     << "</td><td>" << subs_[i].params().trp << "</td>";
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
		*out << "</td><td>" << subs_[i].params().dqm 
		     << "</td><td>" << subs_[i].params().trp << "</td>";
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
  using namespace utils;

  *out << updaterStatic_;
  mDiv(out,"loads");
  uptime(out);
  cDiv(out);
  mDiv(out,"st",fsm_.stateName()->toString());
  mDiv(out,"ru",runNumber_.toString());
  mDiv(out,"nsl",nbSubProcesses_.value_);
  mDiv(out,"cl");
  *out << getApplicationDescriptor()->getClassName() 
       << (nbSubProcesses_.value_ > 0 ? "MP " : " ");
  cDiv(out);
  mDiv(out,"in",getApplicationDescriptor()->getInstance());
  if(fsm_.stateName()->toString() != "Halted" && fsm_.stateName()->toString() != "halting"){
    mDiv(out,"hlt");
    *out << "<a href=\"" << configString_.toString() << "\">HLT Config</a>";
    cDiv(out);
    *out << std::endl;
  }
  else
    mDiv(out,"hlt","Not yet...");

  mDiv(out,"sq",squidPresent_.toString());
  mDiv(out,"vwl",(supervising_ ? "Active" : "Not Initialized"));
  mDiv(out,"mwl",evtProcessor_.wlMonitoring());
  if(nbProcessed != 0 && nbAccepted != 0)
    {
      mDiv(out,"tt",((xdata::UnsignedInteger32*)nbProcessed)->value_);
      mDiv(out,"ac",((xdata::UnsignedInteger32*)nbAccepted)->value_);
    }
  else
    {
      mDiv(out,"tt",0);
      mDiv(out,"ac",0);
    }
  if(!myProcess_)
    mDiv(out,"swl",(wlSummarizeActive_ ? "Active" : "Inactive"));
  else
    mDiv(out,"swl",(wlScalersActive_ ? "Active" : "Inactive"));

  mDiv(out,"idi",iDieUrl_.value_);
  if(vp_!=0){
    mDiv(out,"vpi",(unsigned int) vp_);
    if(vulture_->hasStarted()>=0)
      mDiv(out,"vul","Prowling");
    else
      mDiv(out,"vul","Dead");
  }
  else{
    mDiv(out,"vul",(vulture_==0 ? "Nope" : "Hatching"));
  }    
  if(evtProcessor_){
    mDiv(out,"ll");
    *out << evtProcessor_.lastLumi().ls
	 << "," << evtProcessor_.lastLumi().proc << "," << evtProcessor_.lastLumi().acc;
    cDiv(out);
  }
  mDiv(out,"lg");
  for(unsigned int i = logRingIndex_; i<logRingSize_; i++)
    *out << logRing_[i] << std::endl;
  if(logWrap_)
    for(unsigned int i = 0; i<logRingIndex_; i++)
      *out << logRing_[i] << std::endl;
  cDiv(out);
  
}

void FUEventProcessor::procStat(xgi::Input *in, xgi::Output *out)
{
  evf::utils::procStat(out);
}

void FUEventProcessor::sendMessageOverMonitorQueue(MsgBuf &buf)
{
  if(myProcess_) myProcess_->postSlave(buf,true);
}

void FUEventProcessor::makeStaticInfo()
{
  using namespace utils;
  std::ostringstream ost;
  mDiv(&ost,"ve");
  ost<< "$Revision: 1.112 $ (" << edm::getReleaseVersion() <<")";
  cDiv(&ost);
  mDiv(&ost,"ou",outPut_.toString());
  mDiv(&ost,"sh",hasShMem_.toString());
  mDiv(&ost,"mw",hasModuleWebRegistry_.toString());
  mDiv(&ost,"sw",hasServiceWebRegistry_.toString());
  
  xdata::Serializable *monsleep = 0;
  xdata::Serializable *lstimeout = 0;
  try{
    monsleep  = applicationInfoSpace_->find("monSleepSec");
    lstimeout = applicationInfoSpace_->find("lsTimeOut");
  }
  catch(xdata::exception::Exception e){
  }
  
  if(monsleep!=0)
    mDiv(&ost,"ms",monsleep->toString());
  if(lstimeout!=0)
    mDiv(&ost,"lst",lstimeout->toString());
  char cbuf[sizeof(struct utsname)];
  struct utsname* buf = (struct utsname*)cbuf;
  uname(buf);
  mDiv(&ost,"sysinfo");
  ost << buf->sysname << " " << buf->nodename 
      << " " << buf->release << " " << buf->version << " " << buf->machine;
  cDiv(&ost);
  updaterStatic_ = ost.str();
}

XDAQ_INSTANTIATOR_IMPL(evf::FUEventProcessor)
