
/** \class RPCTriggerPatternBuilder
 *
 * Description:
 *      Class to parse xml files and fill Offline DB
 *
 * $Date: 2008/03/05 12:25:52 $
 * $Revision: 1.4 $
 * \author Michal Szleper -- INS Warsaw
 *
 */
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <time.h>
#include "OnlineDB/Oracle/interface/Oracle.h"

#include <FWCore/ParameterSet/interface/FileInPath.h>
#include "CondFormats/RPCObjects/interface/L1RPCConfig.h"
#include "L1Trigger/RPCTrigger/interface/RPCPatternsParser.h"



using namespace std;
using namespace oracle::occi;

class RPCTriggerPatternBuilder : public edm::EDAnalyzer {
  
public:
  explicit RPCTriggerPatternBuilder(const edm::ParameterSet& iConfig);
  ~RPCTriggerPatternBuilder();
  virtual void beginJob(const edm::EventSetup& iSetup);
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){}
  void WritePoolDB();
  void parsePatterns();

private:

  L1RPCConfig *pL1RPCConfig;
    
  string tag_;
  int m_ppt;
  string m_dataDir;
  string m_patternsDir;

};
                                                                                
RPCTriggerPatternBuilder::RPCTriggerPatternBuilder(const edm::ParameterSet& iConfig) 
  : tag_(iConfig.getUntrackedParameter<std::string>("tag",(string)"RPCTriggerPatterns_v1")),
    m_ppt(iConfig.getUntrackedParameter<int>("PACsPerTower")),
    m_dataDir(iConfig.getUntrackedParameter<std::string>("filedir"))
{
  cout<<"CTor called"<<endl;

  edm::FileInPath fp(m_dataDir+"pacPat_t0sc0sg0.xml");
  std::string patternsDirNameUnstriped = fp.fullPath();
  m_patternsDir = patternsDirNameUnstriped.substr(0,patternsDirNameUnstriped.find_last_of("/")+1);
  
}

RPCTriggerPatternBuilder::~RPCTriggerPatternBuilder()
{ 
  cout<<"DTor called"<<endl;
}

void RPCTriggerPatternBuilder::endJob()
{ 
  try{
    WritePoolDB();
  }
  catch (SQLException &e) {
    cerr << "SQLException:  " << e.getMessage() << endl;
  } catch (std::exception &e) {
    cerr << "std::exception:  " << e.what() << endl;
  } catch (...) {
    cerr << "Unknown exception caught" << endl;
  }
}

void RPCTriggerPatternBuilder::beginJob(const edm::EventSetup& iSetup)
{ 
  cout<<"beginJob called"<<endl;
  try{
    parsePatterns();
  }
  catch (std::exception &e) { 
    cout <<"std::exception:  "<< e.what(); 
  } catch (...) { 
    cout << "Unknown error caught "<<endl; 
  }
}

void RPCTriggerPatternBuilder::WritePoolDB()
{
  cout << endl << "Start writing to PoolDB" << flush << endl;

  edm::Service<cond::service::PoolDBOutputService> mydbservice;
  if( !mydbservice.isAvailable() ){
    cout<<"DB service unavailable"<<endl;
    return;
  } else { 
    cout<<"DB service OK"<<endl; 
  }
  
  int iov=mydbservice->currentTime();
  std::cout<<"Current IOV is "<<iov<<std::endl;
  mydbservice->writeOne<L1RPCConfig>(pL1RPCConfig, iov, "L1RPCConfigRcd");
  cout<<"... all done, end"<<endl;
}

// Method called to parse xml files and fill L1RPCConfig
void RPCTriggerPatternBuilder::parsePatterns()
{

  pL1RPCConfig =  new L1RPCConfig();
  pL1RPCConfig->setPPT(m_ppt);
// parse and insert patterns
  int scCnt = 0, sgCnt = 0;
  if(m_ppt == 1) {
    scCnt = 1;
    sgCnt = 1;
  }
  else if(m_ppt == 12) {
    scCnt = 1;
    sgCnt = 12;
  }
  else if(m_ppt == 144) {
    scCnt = 12;
    sgCnt = 12;
  }
  else {
    throw cms::Exception("BadConfig") << "Bad number of ppt requested: " << m_ppt << "\n";
  }

  for (int tower = 0; tower < RPCConst::m_TOWER_COUNT; ++tower) {
    for (int logSector = 0; logSector < scCnt; ++logSector) {
      for (int logSegment = 0; logSegment < sgCnt; ++logSegment) {

        std::stringstream fname;
        fname << m_patternsDir
              << "pacPat_t" << tower 
              << "sc" << logSector 
	      << "sg" <<logSegment 
              << ".xml";

        cout<<"Parsing "<<fname.str()<<flush;
        RPCPatternsParser parser;
        parser.parse(fname.str());
        cout<<" - done "<<endl;

        RPCPattern::RPCPatVec npats = parser.getPatternsVec(tower, logSector, logSegment);
           
        for (int ip=0; ip<npats.size(); ip++) {
          npats[ip].setCoords(tower,logSector,logSegment);
          pL1RPCConfig->m_pats.push_back(npats[ip]);
        }

        RPCPattern::TQualityVec nquals = parser.getQualityVec(); 
        for (int iq=0; iq<nquals.size(); iq++) {
          nquals[iq].m_tower=tower;
          nquals[iq].m_logsector=logSector;
          nquals[iq].m_logsegment=logSegment;
          pL1RPCConfig->m_quals.push_back(nquals[iq]);
        }
	    
      } // segments
    } // sectors
  } // towers

}


//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCTriggerPatternBuilder);

