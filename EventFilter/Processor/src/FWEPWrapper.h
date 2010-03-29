#ifndef EVENTFILTER_PROCESSOR_FWEPWRAPPER_H
#define EVENTFILTER_PROCESSOR_FWEPWRAPPER_H

#include "FWCore/Framework/interface/EventProcessor.h"

#include "FWCore/PrescaleService/interface/PrescaleService.h"

#include "EventFilter/Utilities/interface/Exception.h"
#include "TriggerReportHelpers.h"

#include "toolbox/lang/Class.h"
#include "toolbox/task/WorkLoop.h"
#include "toolbox/net/URL.h"


#include "xgi/Input.h"
#include "xgi/Output.h"
#include "xgi/exception/Exception.h"

#include "xoap/MessageReference.h"

#include "xdata/InfoSpace.h"
#include "xdata/Integer.h"
#include "log4cplus/logger.h"

#include <string>
#include "MsgBuf.h"

namespace xdaq{
  class ApplicationDescriptor;
  class ApplicationContext;
}

namespace evf{

  struct lsTriplet{
    lsTriplet(): ls(0), proc(0), acc(0){}
    unsigned int ls;
    unsigned int proc;
    unsigned int acc;
  };
  class FWEPWrapper : public toolbox::lang::Class{

  public:

    FWEPWrapper(log4cplus::Logger &, unsigned int instance);
    virtual ~FWEPWrapper();



    void taskWebPage(xgi::Input *,xgi::Output *,const std::string &);
    void summaryWebPage(xgi::Input *,xgi::Output *,const std::string &);
    void moduleWeb(xgi::Input *in,xgi::Output *out);
    void serviceWeb(xgi::Input *in,xgi::Output *out);
    void microState(xgi::Input *in,xgi::Output *out);

    edm::EventProcessor *operator->(){return evtProcessor_;}
    operator bool(){return epInitialized_;}

    // initialize the cmssw event processor
    void init(unsigned short, std::string &);
    void makeServicesOnly();
    void forceInitEventProcessorMaybe(){epInitialized_ = false;}

    // stop the cmssw event processor
    edm::EventProcessor::StatusCode stop();

    // stop and halt (check if stop done)
    void stopAndHalt();

    // infospace pointers transfer
    void setScalersInfoSpace(xdata::InfoSpace *is){scalersInfoSpace_ = is;}
    void setApplicationInfoSpace(xdata::InfoSpace *is){applicationInfoSpace_ = is;}
    void setMonitorInfoSpace(xdata::InfoSpace *is){monitorInfoSpace_ = is;}
    void setMonitorInfoSpaceAlt(xdata::InfoSpace *is){monitorInfoSpaceAlt_ = is;}
    void setMonitorInfoSpaceLegend(xdata::InfoSpace *is){monitorInfoSpaceLegend_ = is;}
    void setRcms(xdaq::ApplicationDescriptor* rcms){rcms_ = rcms;}
    void setAppDesc(xdaq::ApplicationDescriptor *ad){xappDesc_ = ad;}
    void setAppCtxt(xdaq::ApplicationContext *ctx){xappCtxt_ = ctx;}


    void publishConfigAndMonitorItems(unsigned int);
    void publishConfigAndMonitorItemsSP();

    std::string wlMonitoring(){
      if(wlMonitoring_!=0 && wlMonitoring_->isActive()) return (wlMonitoringActive_ ? "active" : "inactive");
      else return "not initialized"; 
    }
    std::string const &configuration() const {return configuration_;}
    // calculate monitoring information in separate thread
    void startMonitoringWorkLoop() throw (evf::Exception);
    bool monitoring(toolbox::task::WorkLoop* wl);
    // trigger report callback
    bool getTriggerReport(bool useLock);
    void updateRollingReport();

    void sumAndPackTriggerReport(MsgBuf &buf);
    void resetPackedTriggerReport(){trh_.resetPackedTriggerReport();}
    void adjustLsIndexForRestart(){trh_.adjustLsIndexForRestart();}
    void resetTriggerReport(){trh_.resetTriggerReport();}
    MsgBuf &getPackedTriggerReport(){return trh_.getPackedTriggerReport();}
    bool fireScalersUpdate(); 
    void createAndSendScalersMessage();

    void lumiSumTable(xgi::Output *out);

    // some accessors for FUEventProcessor
    std::string const &moduleNameFromIndex(unsigned int i) const
      {
	if(i<mapmod_.size()) return mapmod_[i];
	else return unknown;
      }
    std::string const &stateNameFromIndex(unsigned int i) const
      {
	if(i<statmod_.size()) return statmod_[i];
	else return unknown;
      }

    lsTriplet &lastLumi(){return lumiSectionsCtr_[rollingLsIndex_];}
    void resetWaiting(){waitingForLs_ = false;}
    bool isWaitingForLs(){return waitingForLs_;}
    unsigned int getScalersUpdates(){return scalersUpdateCounter_;}
    void resetLumiSectionReferenceIndex(){trh_.resetLumiSectionReferenceIndex();}
    void withdrawLumiSectionIncrement(){trh_.withdrawLumiSectionIncrement();}
    unsigned int getLumiSectionReferenceIndex(){return trh_.getLumiSectionReferenceIndex();}
  private:
    static const std::string        unknown;
    edm::EventProcessor             *evtProcessor_;

    edm::ServiceToken                serviceToken_;    
    edm::ServiceToken                slaveServiceToken_;    
    bool                             servicesDone_;
    bool                             epInitialized_;
    std::string                      configString_;

    // prescale (cmssw framework-) service
    edm::service::PrescaleService*  prescaleSvc_;
    
    
    fuep::TriggerReportHelpers       trh_;
    log4cplus::Logger                log_;
    bool                             isPython_;
    bool                             hasPrescaleService_;
    bool                             hasModuleWebRegistry_;
    bool                             hasServiceWebRegistry_;
						 
    // infospace variables, general
    xdata::InfoSpace                *applicationInfoSpace_;

    // infospace variables, monitor
    xdata::InfoSpace                *monitorInfoSpace_;
    xdata::InfoSpace                *monitorInfoSpaceAlt_;
    xdata::InfoSpace                *monitorInfoSpaceLegend_;

    xdata::UnsignedInteger32         timeoutOnStop_; // in seconds

    std::vector<edm::ModuleDescription const*> descs_; //module description array
    std::map<std::string,int>        modmap_;
    std::vector<std::string>         mapmod_;
    std::vector<std::string>         statmod_;

    // monitoring workloop params
    xdata::UnsignedInteger32         monSleepSec_;
    struct timeval                   monStartTime_;

    // flahslist variables
    xdata::String                    epMState_;
    xdata::String                    epmState_;
    xdata::UnsignedInteger32         nbProcessed_;
    xdata::UnsignedInteger32         nbAccepted_;

    // workloop / action signature for monitoring
    toolbox::task::WorkLoop         *wlMonitoring_;      
    toolbox::task::ActionSignature  *asMonitoring_;
    bool                             wlMonitoringActive_;
    bool                             watching_;

    // flahslist variables, alt
    xdata::Integer                   epMAltState_;
    xdata::Integer                   epmAltState_;


    // flahslist variables, legend
    xdata::String                    macro_state_legend_;
    xdata::String                    micro_state_legend_;

    // LS stuff
    unsigned int                     allPastLumiProcessed_;
    std::list<std::string>           names_;
    std::string                      lsidTimedOutAsString_;
    unsigned int                     lsid_;
    unsigned int                     psid_;

    xdata::InfoSpace                *scalersInfoSpace_;
    xdata::UnsignedInteger32         localLsIncludingTimeOuts_;
    xdata::UnsignedInteger32         lsTimeOut_;
    xdata::Table                     scalersComplete_;
    unsigned int                     scalersUpdateAttempted_;    
    unsigned int                     scalersUpdateCounter_;
    std::vector<lsTriplet>           lumiSectionsCtr_;
    std::vector<bool>                lumiSectionsTo_;
    unsigned int                     rollingLsIndex_;
    bool                             rollingLsWrap_;
    static const unsigned int        lsRollSize_ = 20;
    xdaq::ApplicationDescriptor*     rcms_;
    xdaq::ApplicationDescriptor*     xappDesc_;
    xdaq::ApplicationContext*        xappCtxt_;
    std::string                      configuration_;
    xdata::UnsignedInteger32         instance_;
    bool                             hasSubProcesses;
    bool                             waitingForLs_;
    friend class FUEventProcessor;
  };
}
#endif
