#ifndef DTNoiseAnalysisTest_H
#define DTNoiseAnalysisTest_H


/** \class DTNoiseAnalysisTest
 * *
 *  DQM Test Client
 *
 *  $Date: 2008/07/04 10:14:27 $
 *  $Revision: 1.1 $
 *  \author  G. Mila - INFN Torino
 *   
 */


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/ESHandle.h>


#include <iostream>
#include <string>
#include <map>




class DTGeometry;
class DTChamberId;
class DTSuperLayerId;
class DQMStore;
class MonitorElement;

class DTNoiseAnalysisTest: public edm::EDAnalyzer{

public:

  /// Constructor
  DTNoiseAnalysisTest(const edm::ParameterSet& ps);
  
  /// Destructor
  virtual ~DTNoiseAnalysisTest();

protected:

  /// BeginJob
  void beginJob(const edm::EventSetup& c);

  /// Analyze
  void analyze(const edm::Event& e, const edm::EventSetup& c);

  /// book the summary histograms
  void bookHistos();

  /// Get the ME name
  std::string getMEName(const DTChamberId & chID);

  void beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& context) ;

  /// DQM Client Diagnostic
  void endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& c);


private:

  int nevents;
  
  DQMStore* dbe;
  
  // the dt geometry
  edm::ESHandle<DTGeometry> muonGeom;

  // paramaters from cfg
  int noisyCellDef;

  // wheel summary histograms  
  std::map< int, MonitorElement* > noiseHistos;
  std::map< int, MonitorElement* > noisyCellHistos;
  MonitorElement* summaryNoiseHisto;

};

#endif
