

#include "EventFilter/Utilities/interface/MicroStateService.h"


namespace evf{

  MicroStateService::MicroStateService(const edm::ParameterSet& iPS, 
				       edm::ActivityRegistry& reg)
  {
  
    reg.watchPostBeginJob(this,&MicroStateService::postBeginJob);
    reg.watchPostEndJob(this,&MicroStateService::postEndJob);
  
    reg.watchPreProcessEvent(this,&MicroStateService::preEventProcessing);
    reg.watchPostProcessEvent(this,&MicroStateService::postEventProcessing);
    reg.watchPreSource(this,&MicroStateService::preSource);
    reg.watchPostSource(this,&MicroStateService::postSource);
  
    reg.watchPreModule(this,&MicroStateService::preModule);
    reg.watchPostModule(this,&MicroStateService::postModule);
    microstate1_ = "BJ";
    microstate2_ = "INIT";

  }


  MicroStateService::~MicroStateService()
  {
  }

  void MicroStateService::postBeginJob()
  {
    boost::mutex::scoped_lock sl(lock_);
    microstate1_ = "BJD";
  }

  void MicroStateService::postEndJob()
  {
    boost::mutex::scoped_lock sl(lock_);
    microstate1_ = "EJ";
    microstate2_ = "done";
  }

  void MicroStateService::preEventProcessing(const edm::EventID& iID,
					     const edm::Timestamp& iTime)
  {
    boost::mutex::scoped_lock sl(lock_);
    microstate1_ = "PRO";
  }

  void MicroStateService::postEventProcessing(const edm::Event& e, const edm::EventSetup&)
  {
    boost::mutex::scoped_lock sl(lock_);
    microstate2_ = "IN";
  }
  void MicroStateService::preSource()
  {
    boost::mutex::scoped_lock sl(lock_);
    microstate2_ = "IN";
  }

  void MicroStateService::postSource()
  {
    boost::mutex::scoped_lock sl(lock_);
    microstate2_ = "IND";
  }

  void MicroStateService::preModule(const edm::ModuleDescription& desc)
  {
    boost::mutex::scoped_lock sl(lock_);
    microstate2_ = desc.moduleLabel_;
  }

  void MicroStateService::postModule(const edm::ModuleDescription& desc)
  {
  }
  
  std::string MicroStateService::getMicroState1()
  { 
	boost::mutex::scoped_lock sl(lock_);
	return microstate1_;
  }

  std::string MicroStateService::getMicroState2()
  { 
	boost::mutex::scoped_lock sl(lock_);
	return microstate2_;
  }

  void MicroStateService::setMicroState(std::string &in)
  {
    	boost::mutex::scoped_lock sl(lock_);
	microstate2_ = in;
  }

} //end namespace evf

