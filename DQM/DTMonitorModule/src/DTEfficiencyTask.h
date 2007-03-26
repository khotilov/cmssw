#ifndef DTEfficiencyTask_H
#define DTEfficiencyTask_H


/** \class DTEfficiencyTask
 *  DQM Analysis of 4D DT segments, it produces plots about: <br>
 *      - single cell efficiency 
 *  All histos are produced per Layer
 *
 *
 *  $Date: 2007/03/26 17:30:00 $
 *  $Revision: 1.0 $
 *  \author G. Mila - INFN Torino
 */



#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/EDAnalyzer.h>

#include <string>
#include <map>
#include <vector>

class DaqMonitorBEInterface;
class MonitorElement;


class DTEfficiencyTask: public edm::EDAnalyzer{
public:
  /// Constructor
  DTEfficiencyTask(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~DTEfficiencyTask();

  /// BeginJob
  void beginJob(const edm::EventSetup& c);

  /// Endjob
  void endJob();

  // Operations
  void analyze(const edm::Event& event, const edm::EventSetup& setup);

protected:


private:
  DaqMonitorBEInterface* theDbe;

  // Switch for verbosity
  bool debug;
  std::string theRootFileName;
  bool writeHisto;

  // Lable of 4D segments in the event
  std::string theRecHits4DLabel;
  
  edm::ParameterSet parameters;

  // Book a set of histograms for a give chamber
  void bookHistos(DTLayerId lId, int fisrtWire, int lastWire);

  // Fill a set of histograms for a given L 
  void fillHistos(DTLayerId lId, int firstWire, int lastWire);
  void fillHistos(DTLayerId lId, int firstWire, int lastWire, int missingWire);

  std::map<DTLayerId, std::vector<MonitorElement*> > histosPerL;

};
#endif

