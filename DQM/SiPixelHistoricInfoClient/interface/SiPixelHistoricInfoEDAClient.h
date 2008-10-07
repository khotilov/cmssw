#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelPerformanceSummary.h"

// #include "DQM/SiPixelHistoricInfoClient/interface/SiPixelHistoricInfoWebInterface.h"


// namespace edm {
//   class ParameterSet;
//   class Event;
//   class EventId;
//   class Timestamp;
// } 

class SiPixelHistoricInfoEDAClient : public edm::EDAnalyzer {
public:
  explicit SiPixelHistoricInfoEDAClient(const edm::ParameterSet&);
 ~SiPixelHistoricInfoEDAClient();

private:
  virtual void beginJob(const edm::EventSetup&);
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void retrievePointersToModuleMEs();
  void fillSummaryObjects() const;
  void writetoDB() const; 
  void savetoFile(std::string) const; 

private: 
  bool printDebug_;
  bool writeHisto_;
  std::string outputDir_; 
  edm::ParameterSet parameterSet_;
  DQMStore* dbe_;

  bool firstEventInRun; 
  int nEventsInRun; 
  std::map< uint32_t, std::vector<MonitorElement*> > ClientPointersToModuleMEs;
  SiPixelPerformanceSummary* performanceSummary;

  // SiPixelHistoricInfoWebInterface* webInterface_;
  // bool defaultWebPageCreated_; 
};

