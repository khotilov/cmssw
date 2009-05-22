#ifndef DQMSERVICES_CORE_STANDALONE_H
# define DQMSERVICES_CORE_STANDALONE_H
# if !WITHOUT_CMS_FRAMEWORK
#  include "FWCore/ParameterSet/interface/ParameterSet.h"
#  include "FWCore/ServiceRegistry/interface/ServiceRegistry.h"
#  include "FWCore/ServiceRegistry/interface/Service.h"
#  include "FWCore/MessageLogger/interface/JobReport.h"
#  include "FWCore/Utilities/interface/GetReleaseVersion.h"
# else
#  include <string>
#  include <vector>
#  include <map>

namespace edm
{
  std::string getReleaseVersion(void)
  { return "CMSSW_STANDALONE"; }

  struct ParameterSet
  {
    template <class T> static const T &
    getUntrackedParameter(const char * /* key */, const T &value)
    { return value; }
  };

  struct ServiceToken
  {
    ServiceToken(int) {}
  };

  struct ServiceRegistry
  {
    struct Operate
    {
      Operate(const ServiceToken &) {}
    };

    static int createSet(const std::vector<ParameterSet> &) { return 0; }
  };

  template <class T>
  struct Service
  {
    bool isAvailable(void) { return false; }
    T *operator->(void)
    {
      static char buf[sizeof(T)]; static T *x;
      if (! x) x = new (buf) T(ParameterSet());
      return x;
    }
    T &operator*(void) { return * operator->(); }
  };

  struct JobReport
  {
    JobReport(const edm::ParameterSet &) {}
    void reportAnalysisFile(const std::string &, const std::map<std::string, std::string> &) {}
  };
}
# endif // WITHOUT_CMS_FRAMEWORK
#endif // DQMSERVICES_CORE_STANDALONE_H
