#ifndef SiStripRecHitsValid_h
#define SiStripRecHitsValid_h

/* \class SiStripRecHitsValid
 *
 * Analyzer to validate RecHits in the Strip tracker
 *
 * \author Patrizia Azzi, INFN PD
 *
 * \version   1st version May 2006
 *
 ************************************************************/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
//only mine
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//DQM services for histogram
#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//--- for SimHit association
 #include "SimDataFormats/TrackingHit/interface/PSimHit.h"  
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h" 

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h" 
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h" 
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"

#include <string>


using namespace std;
using namespace edm;


class SiStripRecHitsValid : public edm::EDAnalyzer {

 public:
  
  SiStripRecHitsValid(const edm::ParameterSet& conf);
  
  ~SiStripRecHitsValid();

 protected:

  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  void beginJob(const EventSetup& c);
  void endJob();
  
 private: 
  //Back-End Interface
  DaqMonitorBEInterface* dbe_;
  string outputFile_;
  //TIB
  MonitorElement* meNstpRphiTIB[4];
  MonitorElement* meAdcRphiTIB[4];
  MonitorElement* mePosxRphiTIB[4];
  MonitorElement* meErrxRphiTIB[4];
  MonitorElement* meResRphiTIB[4];
  MonitorElement* meNstpSasTIB[4];
  MonitorElement* meAdcSasTIB[4];
  MonitorElement* mePosxSasTIB[4];
  MonitorElement* meErrxSasTIB[4];
  MonitorElement* meResSasTIB[4];
  MonitorElement* mePosxMatchedTIB[2];
  MonitorElement* mePosyMatchedTIB[2];
  MonitorElement* meErrxMatchedTIB[2];
  MonitorElement* meErryMatchedTIB[2];
  MonitorElement* meResxMatchedTIB[2];
  MonitorElement* meResyMatchedTIB[2];
  //TOB
  MonitorElement* meNstpRphiTOB[6];
  MonitorElement* meAdcRphiTOB[6];
  MonitorElement* mePosxRphiTOB[6];
  MonitorElement* meErrxRphiTOB[6];
  MonitorElement* meResRphiTOB[6];
  MonitorElement* meNstpSasTOB[2];
  MonitorElement* meAdcSasTOB[2];
  MonitorElement* mePosxSasTOB[2];
  MonitorElement* meErrxSasTOB[2];
  MonitorElement* meResSasTOB[2];
  MonitorElement* mePosxMatchedTOB[2];
  MonitorElement* mePosyMatchedTOB[2];
  MonitorElement* meErrxMatchedTOB[2];
  MonitorElement* meErryMatchedTOB[2];
  MonitorElement* meResxMatchedTOB[2];
  MonitorElement* meResyMatchedTOB[2];
  //TID
  MonitorElement* meNstpRphiTID[3];
  MonitorElement* meAdcRphiTID[3];
  MonitorElement* mePosxRphiTID[3];
  MonitorElement* meErrxRphiTID[3];
  MonitorElement* meResRphiTID[3];
  MonitorElement* meNstpSasTID[2];
  MonitorElement* meAdcSasTID[2];
  MonitorElement* mePosxSasTID[2];
  MonitorElement* meErrxSasTID[2];
  MonitorElement* meResSasTID[2];
  MonitorElement* mePosxMatchedTID[2];
  MonitorElement* mePosyMatchedTID[2];
  MonitorElement* meErrxMatchedTID[2];
  MonitorElement* meErryMatchedTID[2];
  MonitorElement* meResxMatchedTID[2];
  MonitorElement* meResyMatchedTID[2];
  //TEC
  MonitorElement* meNstpRphiTEC[7];
  MonitorElement* meAdcRphiTEC[7];
  MonitorElement* mePosxRphiTEC[7];
  MonitorElement* meErrxRphiTEC[7];
  MonitorElement* meResRphiTEC[7];
  MonitorElement* meNstpSasTEC[5];
  MonitorElement* meAdcSasTEC[5];
  MonitorElement* mePosxSasTEC[5];
  MonitorElement* meErrxSasTEC[5];
  MonitorElement* meResSasTEC[5];
  MonitorElement* mePosxMatchedTEC[5];
  MonitorElement* mePosyMatchedTEC[5];
  MonitorElement* meErrxMatchedTEC[5];
  MonitorElement* meErryMatchedTEC[5];
  MonitorElement* meResxMatchedTEC[5];
  MonitorElement* meResyMatchedTEC[5];

  std::vector<PSimHit> matched;
  std::pair<LocalPoint,LocalVector> projectHit( const PSimHit& hit, const StripGeomDetUnit* stripDet,
							const BoundPlane& plane);
  edm::ParameterSet conf_;
  const StripTopology* topol;

  static const int MAXHIT = 100;
  float rechitrphix[MAXHIT];
  float rechitrphierrx[MAXHIT];
  float rechitrphiy[MAXHIT];
  float rechitrphiz[MAXHIT];
  float rechitrphiphi[MAXHIT];
  float rechitrphires[MAXHIT];
  int   clusizrphi[MAXHIT];
  float cluchgrphi[MAXHIT];
  float rechitsasx[MAXHIT];
  float rechitsaserrx[MAXHIT];
  float rechitsasy[MAXHIT];
  float rechitsasz[MAXHIT];
  float rechitsasphi[MAXHIT];
  float rechitsasres[MAXHIT];
  int   clusizsas[MAXHIT];
  float cluchgsas[MAXHIT];
  float rechitmatchedx[MAXHIT];
  float rechitmatchedy[MAXHIT];
  float rechitmatchedz[MAXHIT];
  float rechitmatchederrxx[MAXHIT];
  float rechitmatchederrxy[MAXHIT];
  float rechitmatchederryy[MAXHIT];
  float rechitmatchedphi[MAXHIT];
  float rechitmatchedresx[MAXHIT];
  float rechitmatchedresy[MAXHIT];

};

#endif
