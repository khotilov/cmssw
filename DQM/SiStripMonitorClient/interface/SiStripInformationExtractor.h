#ifndef _SiStripInformationExtractor_h_
#define _SiStripInformationExtractor_h_

#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "xgi/Utils.h"
#include "xgi/Method.h"


#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <map>

class SiStripLayoutParser;
class SiStripDetCabling;
class DQMStore;
class QReport;
class SiStripHistoPlotter;

class SiStripInformationExtractor {

 public:

  SiStripInformationExtractor();
 ~SiStripInformationExtractor();

  void getSingleModuleHistos(DQMStore * dqm_store, 
       const std::multimap<std::string, std::string>& req_map, xgi::Output * out);
  void getGlobalHistos(DQMStore* dqm_store, 
       const std::multimap<std::string, std::string>& req_map, xgi::Output * out);
  void getHistosFromPath(DQMStore * dqm_store, 
       const std::multimap<std::string, std::string>& req_map, xgi::Output * out);
  void getTrackerMapHistos(DQMStore* dqm_store, 
       const std::multimap<std::string, std::string>& req_map, xgi::Output * out);


  void readModuleAndHistoList(DQMStore* dqm_store,const edm::ESHandle<SiStripDetCabling>& detcabling,xgi::Output * out);
  void plotHistosFromLayout(DQMStore * dqm_store);
  void readSummaryHistoTree(DQMStore* dqm_store, std::string& str_name, 
                xgi::Output * out);
  void readAlarmTree(DQMStore* dqm_store, std::string& str_name, 
                xgi::Output * out);
 
  void readStatusMessage(DQMStore* dqm_store, std::multimap<std::string, std::string>& req_map, xgi::Output * out);
  void readGlobalHistoList(DQMStore* dqm_store, std::string& dname, xgi::Output * out);
  void readLayoutNames(std::multimap<std::string, std::string>& req_map, xgi::Output * out);

  void readQTestSummary(DQMStore* dqm_store, std::string type, const edm::ESHandle<SiStripDetCabling>& detcabling, xgi::Output * out);

  void readNonGeomHistoTree(DQMStore* dqm_store, std::string& fld_name, xgi::Output * out);
  
  void createImages(DQMStore* dqm_store);
  void getIMGCImage(const std::multimap<std::string, std::string>& req_map, xgi::Output * out);
  


 private:

  void readConfiguration();


 
  void getItemList(const std::multimap<std::string, std::string>& req_map,
                   std::string item_name, std::vector<std::string>& items);
  bool hasItem(std::multimap<std::string, std::string>& req_map,
	      std::string item_name);
  std::string getItemValue(const std::multimap<std::string, std::string>& req_map,
	      std::string item_name);
  void printSummaryHistoList(DQMStore* dqm_store, std::ostringstream& str_val);
  void printAlarmList(DQMStore * dqm_store, std::ostringstream& str_val);
  void printNonGeomHistoList(DQMStore * dqm_store, std::ostringstream& str_val);

  void selectImage(std::string& name, int status);
  void selectImage(std::string& name, std::vector<QReport*> & reports);
  void selectColor(std::string& col, int status);
  void selectColor(std::string& col, std::vector<QReport*>& reports);

  
  void setHTMLHeader(xgi::Output * out);
  void setXMLHeader(xgi::Output * out);
  void setPlainHeader(xgi::Output * out);


  SiStripLayoutParser* layoutParser_;

  std::map<std::string, std::vector< std::string > > layoutMap;
  bool  readReference_;
 
  SiStripHistoPlotter* histoPlotter_;
};
#endif
