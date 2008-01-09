#include "DataFormats/HLTReco/interface/HLTResult.h"
#include "DataFormats/HLTReco/interface/ModuleTiming.h"
#include "DataFormats/HLTReco/interface/HLTPerformanceInfo.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/Common/interface/Wrapper.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Ref.h"

namespace {
  namespace {
    boost::transform_iterator<HLTPerformanceInfo::Path::Adapter,__gnu_cxx::__normal_iterator<const unsigned int*,std::vector<unsigned int> >,boost::use_default,boost::use_default> hltfubar1;
    boost::transform_iterator<HLTPerformanceInfo::Path::Adapter,__gnu_cxx::__normal_iterator<const unsigned long*,std::vector<unsigned long> >,boost::use_default,boost::use_default> hltfubar2;
  }
}

namespace {
  namespace {

    reco::HLTResult< 8> h1;
    reco::HLTResult<16> h2;
    reco::HLTResult<24> h3;

    edm::Wrapper<reco::HLTResult< 8> > w1;
    edm::Wrapper<reco::HLTResult<16> > w2;
    edm::Wrapper<reco::HLTResult<24> > w3;

    edm::EventTime                                et0;

    edm::Wrapper<edm::EventTime>                wet10;

    // Performance Information
    HLTPerformanceInfo pw0;
    edm::Wrapper<HLTPerformanceInfo> pw1;
    HLTPerformanceInfoCollection pw2; 
    edm::Wrapper<HLTPerformanceInfoCollection> pw3; 

    HLTPerformanceInfo::Module pw4;
    HLTPerformanceInfo::Path pw6;
    std::vector<HLTPerformanceInfo::Module> pw8;
    std::vector<HLTPerformanceInfo::Module>::const_iterator pw9;
    std::vector<HLTPerformanceInfo::Path> pw10;
    std::vector<HLTPerformanceInfo::Path>::const_iterator pw11;
    HLTPerformanceInfo::Path::Adapter pw12;
    HLTPerformanceInfo::Path::const_iterator pw13;

    edm::Ref<reco::CompositeCandidateCollection> rxx;

    trigger::TriggerObjectCollection toc;
    trigger::TriggerRefsCollections trc;
    trigger::TriggerFilterObjectWithRefs tfowr;
    trigger::TriggerEvent te;
    trigger::TriggerEventWithRefs tewr;

    edm::Wrapper<trigger::TriggerObjectCollection> wtoc;
    edm::Wrapper<trigger::TriggerFilterObjectWithRefs> wtfowr;
    edm::Wrapper<trigger::TriggerEvent> wte;
    edm::Wrapper<trigger::TriggerEventWithRefs> wtewr;

  }
}
