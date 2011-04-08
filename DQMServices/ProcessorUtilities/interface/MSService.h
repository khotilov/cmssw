#ifndef DQMEVF_MSSERVICE_H
#define DQMEVF_MSSERVICE_H


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Provenance/interface/ModuleDescription.h"

#include "boost/thread/thread.hpp"

#include <string>
#include <vector>

namespace dqmevf {

    class MSService
    {
    public:
      MSService(const edm::ParameterSet&,edm::ActivityRegistry&);
      ~MSService();
      
      std::string getMicroState1();
      
      std::string getMicroState2();

      void postBeginJob();
      void postEndJob();
      
      void preEventProcessing(const edm::EventID&, const edm::Timestamp&);
      void postEventProcessing(const edm::Event&, const edm::EventSetup&);
      
      void preSource();
      void postSource();
      
      
      void preModule(const edm::ModuleDescription&);
      void postModule(const edm::ModuleDescription&);

      void setMicroState(std::string &);
      
    private:

      std::string microstate1_;
      std::string microstate2_;

      boost::mutex lock_;

    };

}

#endif
