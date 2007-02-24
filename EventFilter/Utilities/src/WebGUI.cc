////////////////////////////////////////////////////////////////////////////////
//
// WebGUI
// ------
//
//            10/19/2006 Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
////////////////////////////////////////////////////////////////////////////////


#include "EventFilter/Utilities/interface/WebGUI.h"

#include "xcept/include/xcept/Exception.h"
#include "xcept/include/xcept/tools.h"

#include "extern/cgicc/linuxx86/include/cgicc/CgiDefs.h"
#include "extern/cgicc/linuxx86/include/cgicc/Cgicc.h"
#include "extern/cgicc/linuxx86/include/cgicc/HTMLClasses.h"

#include <sstream>


using namespace std;
using namespace evf;
using namespace cgicc;


////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
WebGUI::WebGUI(xdaq::Application* app,toolbox::fsm::FiniteStateMachine* fsm)
  : app_(app)
  , fsm_(fsm)
  , log_(app->getApplicationContext()->getLogger())
  , appInfoSpace_(0)
  , monInfoSpace_(0)
  , parametersExported_(false)
  , countersAddedToParams_(false)
  , largeAppIcon_("/daq/evb/examples/fu/images/fu64x64.gif")
  , smallAppIcon_("/daq/evb/examples/fu/images/fu32x32.gif")
  , smallDbgIcon_("/daq/evb/bu/images/debug32x32.gif")
  , hyperDAQIcon_("/daq/xdaq/hyperdaq/images/HyperDAQ.jpg")
{
  // initialize application information
  string       xmlClass=app_->getApplicationDescriptor()->getClassName();
  unsigned int instance=app_->getApplicationDescriptor()->getInstance();
  stringstream oss;
  oss<<xmlClass<<instance;
  
  sourceId_=oss.str();
  urn_     ="/"+app_->getApplicationDescriptor()->getURN();
  
  std::stringstream oss2;
  oss2<<"urn:xdaq-monitorable:"<<xmlClass<<":"<<instance;
  string monInfoSpaceName=oss2.str();

  appInfoSpace_=app_->getApplicationInfoSpace();
  monInfoSpace_=xdata::InfoSpace::get(monInfoSpaceName);
  

  // bind xgi callbacks
  xgi::bind(this,&WebGUI::defaultWebPage,"defaultWebPage");
  xgi::bind(this,&WebGUI::debugWebPage,  "debugWebPage");
  xgi::bind(this,&WebGUI::css,           "styles.css");
}


//______________________________________________________________________________
WebGUI::~WebGUI()
{
  
}


////////////////////////////////////////////////////////////////////////////////
// implementation of public member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void WebGUI::defaultWebPage(Input_t *in,Output_t *out) throw (WebGUI::XgiException_t)
{
  updateParams();
  
  *out<<html()<<endl;
  htmlHead(in,out,sourceId_);
  *out<<body()<<endl;
  htmlHeadline(in,out,urn_+"/debugWebPage",smallDbgIcon_);
  *out<<"<table cellpadding=\"25\"><tr valign=\"top\"><td>"<<endl;
  htmlTable(in,out,"Standard Parameters",standardParams_);
  *out<<"</td><td>"<<endl;
  htmlTable(in,out,"Monitored Parameters",monitorParams_);
  *out<<"</td></tr></table>"<<endl;
  *out<<body()<<endl<<html()<<endl;
  return;
}


//______________________________________________________________________________
void WebGUI::debugWebPage(Input_t *in,Output_t *out) throw (WebGUI::XgiException_t)
{
  updateParams();

  *out<<html()<<endl;
  htmlHead(in,out,sourceId_+" [DEBUG]");
  *out<<body()<<endl;
  htmlHeadline(in,out,urn_+"/defaultWebPage",smallAppIcon_);
  *out<<"<table cellpadding=\"25\"><tr valign=\"top\"><td>"<<endl;
  htmlTable(in,out,"Debug Parameters",debugParams_);
  *out<<"</td></tr></table>"<<endl;
  *out<<body()<<endl<<html()<<endl;
  return;
}


//______________________________________________________________________________
void WebGUI::css(Input_t *in,Output_t *out) throw (WebGUI::XgiException_t)
{
  css_.css(in,out);
}


//______________________________________________________________________________
void WebGUI::addStandardParam(CString_t& name,Param_t* param)
{
  if (parametersExported_) {
    LOG4CPLUS_ERROR(log_,"Failed to add standard parameter '"<<name<<"'.");
    return;
  }
  standardParams_.push_back(make_pair(name,param));
}


//______________________________________________________________________________
void WebGUI::addMonitorParam(CString_t& name,Param_t* param)
{
  if (parametersExported_) {
    LOG4CPLUS_ERROR(log_,"Failed to add monitor parameter '"<<name<<"'.");
    return;
  }
  monitorParams_.push_back(make_pair(name,param));
}


//______________________________________________________________________________
void WebGUI::addDebugParam(CString_t& name,Param_t* param)
{
  if (parametersExported_) {
    LOG4CPLUS_ERROR(log_,"Failed to add debug parameter '"<<name<<"'.");
    return;
  }
  debugParams_.push_back(make_pair(name,param));
}


//______________________________________________________________________________
void WebGUI::addStandardCounter(CString_t& name,Counter_t* counter)
{
  if (countersAddedToParams_) {
    LOG4CPLUS_ERROR(log_,"can't add standard counter '"<<name
		    <<"' to WebGUI of "<<sourceId_);
  }
  standardCounters_.push_back(make_pair(name,counter));
}


//______________________________________________________________________________
void WebGUI::addMonitorCounter(CString_t& name,Counter_t* counter)
{
  if (countersAddedToParams_) {
    LOG4CPLUS_ERROR(log_,"can't add monitor counter '"<<name
		    <<"' to WebGUI of "<<sourceId_);
  }
  monitorCounters_.push_back(make_pair(name,counter));
}


//______________________________________________________________________________
void WebGUI::addDebugCounter(CString_t& name,Counter_t* counter)
{
  if (countersAddedToParams_) {
    LOG4CPLUS_ERROR(log_,"can't add debug counter '"<<name
		    <<"' to WebGUI of "<<sourceId_);
  }
  debugCounters_.push_back(make_pair(name,counter));
}


//______________________________________________________________________________
void WebGUI::exportParameters()
{
  if (parametersExported_) return;

  if (!countersAddedToParams_) addCountersToParams();

  addParamsToInfoSpace(standardParams_,appInfoSpace());
  addParamsToInfoSpace(monitorParams_, appInfoSpace());
  addParamsToInfoSpace(debugParams_,   appInfoSpace());

  addParamsToInfoSpace(monitorParams_,monInfoSpace());
  
  parametersExported_=true;
}


//______________________________________________________________________________
void WebGUI::resetCounters()
{
  // standard counters
  for (unsigned int i=0;i<standardCounters_.size();i++) {
    Counter_t* counter=standardCounters_[i].second;
    *counter=0;
  }
  // monitor counters
  for (unsigned int i=0;i<monitorCounters_.size();i++) {
    Counter_t* counter=monitorCounters_[i].second;
    *counter=0;
  }
  // debug counters
  for (unsigned int i=0;i<debugCounters_.size();i++) {
    Counter_t* counter=debugCounters_[i].second;
    *counter=0;
  }
}


//______________________________________________________________________________
void WebGUI::addItemChangedListener(CString_t& name,xdata::ActionListener* l)
{
  if (!parametersExported_) {
    LOG4CPLUS_ERROR(log_,"Can't add ItemChangedListener for parameter '"<<name
		    <<"' before WebGUI::exportParameters() is called.");
    return;
  }
  
  try {
    appInfoSpace()->addItemChangedListener(name,l);
  }
  catch (xcept::Exception) {
    LOG4CPLUS_ERROR(log_,"failed to add ItemChangedListener to "
		    <<"application infospace for parameter '"<<name<<"'.");
  }
  
  if (isMonitorParam(name)) {
    try {
      monInfoSpace()->addItemChangedListener(name,l);
    }
    catch (xcept::Exception) {
      LOG4CPLUS_ERROR(log_,"failed to add ItemChangedListener to "
		      <<"monitor infospace for parameter '"<<name<<"'.");
    }
  }
  
}


//______________________________________________________________________________
void WebGUI::addItemRetrieveListener(CString_t& name,xdata::ActionListener* l)
{
  if (!parametersExported_) {
    LOG4CPLUS_ERROR(log_,"Can't add ItemRetrieveListener for parameter '"<<name
		    <<"' before WebGUI::exportParameters() is called.");
    return;
  }
  
  try {
    appInfoSpace()->addItemRetrieveListener(name,l);
    updateParams_.push_back(make_pair(name,l));
  }
  catch (xcept::Exception) {
    LOG4CPLUS_ERROR(log_,"failed to add ItemRetrieveListener to "
		    <<"application infospace for parameter '"<<name<<"'.");
  }
  
  if (isMonitorParam(name)) {
    try {
      monInfoSpace()->addItemRetrieveListener(name,l);
    }
    catch (xcept::Exception) {
      LOG4CPLUS_ERROR(log_,"failed to add ItemRetrieveListener to "
		      <<"monitor infospace for parameter '"<<name<<"'.");
    }
  }
}


//______________________________________________________________________________
void WebGUI::lockInfoSpaces()
{
  appInfoSpace()->lock();
  monInfoSpace()->lock();
}


//______________________________________________________________________________
void WebGUI::unlockInfoSpaces()
{
  appInfoSpace()->unlock();
  monInfoSpace()->unlock();
}



////////////////////////////////////////////////////////////////////////////////
// implementation of private member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void WebGUI::addParamsToInfoSpace(const ParamVec_t& params,
				  xdata::InfoSpace* infoSpace)
{
  for (unsigned int i=0;i<params.size();i++) {
    string   name =params[i].first;
    Param_t* value=params[i].second;
    try {
      infoSpace->fireItemAvailable(name,value);
    }
    catch (xcept::Exception &e) {
      LOG4CPLUS_ERROR(log_,"Can't add parameter '"<<name<<"' to info space '"
		      <<infoSpace->name()<<"': "
		      <<xcept::stdformat_exception_history(e));
    }
  }
}

    
//______________________________________________________________________________
void WebGUI::addCountersToParams()
{
  if (countersAddedToParams_) return;

  // standard counters 
  for (unsigned int i=0;i<standardCounters_.size();i++) {
    standardParams_.push_back(make_pair(standardCounters_[i].first,
					standardCounters_[i].second));
  }
  // monitor counters 
  for (unsigned int i=0;i<monitorCounters_.size();i++) {
    monitorParams_.push_back(make_pair(monitorCounters_[i].first,
				       monitorCounters_[i].second));
  }
  // debug counters 
  for (unsigned int i=0;i<debugCounters_.size();i++) {
    debugParams_.push_back(make_pair(debugCounters_[i].first,
				     debugCounters_[i].second));
  }
  countersAddedToParams_=true;
}


//______________________________________________________________________________
bool WebGUI::isMonitorParam(CString_t& name)
{
  ParamVec_t::const_iterator it;
  for (it=monitorParams_.begin();it!=monitorParams_.end();++it)
    if (it->first==name) return true;
  return false;
}


//______________________________________________________________________________
void WebGUI::updateParams()
{
  UpdateVec_t::iterator it;
  for (it=updateParams_.begin();it!=updateParams_.end();++it)
    appInfoSpace()->fireItemValueRetrieve(it->first,it->second);
}


//______________________________________________________________________________
void WebGUI::htmlTable(Input_t*in,Output_t*out,
		       CString_t& title,const ParamVec_t& params)
{
  *out<<table().set("frame","void").set("rules","rows")
               .set("class","modules").set("width","300")<<endl
      <<tr()<<th(title).set("colspan","2")<<tr()<<endl
      <<tr()
      <<th("Parameter").set("align","left")
      <<th("Value").set("align","right")
      <<tr()
      <<endl;

  for (unsigned int i=0;i<params.size();i++) {
    string valueAsString;
    try {
      valueAsString = params[i].second->toString();
    }
    catch (xcept::Exception& e) {
      valueAsString = e.what();
    }
    *out<<tr()
	<<td(params[i].first).set("align","left")
	<<td(valueAsString).set("align","right")
	<<tr()<<endl;
  }
  *out<<table()<<endl;
}


//______________________________________________________________________________
void WebGUI::htmlHead(Input_t *in,Output_t* out,CString_t& pageTitle)
{
  *out<<head()<<endl<<cgicc::link().set("type","text/css")
                                   .set("rel","stylesheet")
                                   .set("href",urn_+"/styles.css")
      <<endl<<title(pageTitle.c_str())<<endl<<head()<<endl;
}


//______________________________________________________________________________
void WebGUI::htmlHeadline(Input_t *in,Output_t *out,CString_t& link,CString_t& icon)
{
  int    state    =fsm_->getCurrentState();
  string stateName=fsm_->getStateName(state);
  
  *out<<table().set("border","0").set("width","100%")<<endl
      <<tr()<<td().set("align","left")<<endl
      <<img().set("align","middle").set("src",largeAppIcon_)
             .set("alt","main")    .set("width","64")
             .set("height","64")   .set("border","")
      <<endl
      <<b()<<sourceId_<<"    "<<stateName<<b()<<endl
      <<td()<<endl
      <<td().set("width","32")<<endl
      <<a().set("href","/urn:xdaq-application:lid=3")
      <<img().set("align","middle").set("src",hyperDAQIcon_)
             .set("alt","HyperDAQ").set("width","32")
             .set("height","32")   .set("border","")
      <<a()
      <<td()<<endl
      <<td().set("width","32")<<td() 
      <<td().set("width","32")
      <<a().set("href",link)
      <<img().set("align","middle").set("src",icon)
             .set("alt","Debug")   .set("width","32")
             .set("height","32")   .set("border","")
      <<a()
      <<td()<<tr()<<table()<<endl;
  *out<<hr()<<endl;
}
