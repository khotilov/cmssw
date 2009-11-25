#ifndef CondTools_Luminosity_MixedSource_h
#define CondTools_Luminosity_MixedSource_h
#include <vector>
#include <string>
#include "CondTools/Luminosity/interface/LumiRetrieverBase.h"
class TFile;
namespace coral{
  class ConnectionService;
}
namespace edm{
  class ParameterSet;
}
namespace lumi{
  class LumiSectionData;
  /**This source takes lumi measurement data from Lumi Root file and L1 data from L1 database 
   **/
  class MixedSource:public LumiRetrieverBase{
  public:
    explicit MixedSource(const edm::ParameterSet& pset);
    virtual ~MixedSource(){}
    virtual const std::string 
      fill(std::vector< std::pair<lumi::LumiSectionData*,cond::Time_t> >& result, bool allowForceFirstSince=false);
  private:
    //per run information
    typedef std::vector<std::string> TriggerNameResult_Algo;
    typedef std::vector<std::string> TriggerNameResult_Tech;
    typedef std::vector<unsigned int> PrescaleResult_Algo;
    typedef std::vector<unsigned int> PrescaleResult_Tech;
    //per lumisection information
    typedef unsigned int DEADCOUNT;
    typedef std::vector<DEADCOUNT> TriggerDeadCountResult;
    typedef std::vector<unsigned int> BITCOUNT;
    typedef std::vector<BITCOUNT> TriggerCountResult_Algo;
    typedef std::vector<BITCOUNT> TriggerCountResult_Tech;
    void initDB();
    void getTrgData(unsigned int runnumber,
		    coral::ISessionProxy* session,
		    TriggerNameResult_Algo& algonames,
		    TriggerNameResult_Tech& technames,
		    PrescaleResult_Algo& algoprescale,
		    PrescaleResult_Tech& techprescale,
		    TriggerCountResult_Algo& algocount,
		    TriggerCountResult_Tech& techcount,
		    TriggerDeadCountResult& deadtime
		    );
    std::string int2str(int t);
    void printCountResult(const TriggerCountResult_Algo& algo,const TriggerCountResult_Tech& tech);
    void printDeadTimeResult(const TriggerDeadCountResult& result);
    void printTriggerNameResult(const TriggerNameResult_Algo& algonames,const TriggerNameResult_Tech& technames);
    void printPrescaleResult(const PrescaleResult_Algo& algo,const PrescaleResult_Tech& tech);
    
  private:
    std::string m_filename;
    TFile* m_source;
    std::string m_authpath;
    std::string m_lumiversion;
    std::string m_trgdb;
    coral::ConnectionService* m_dbservice;
  };
}//ns lumi
#endif
