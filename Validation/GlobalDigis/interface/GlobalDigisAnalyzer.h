#ifndef GlobalDigisAnalyzer_h
#define GlobalDigisAnalyzer_h

/** \class GlobalDigiAnalyzer
 *  
 *  Class to fill PGlobalDigi object to be inserted into data stream 
 *  containing information about various sub-systems in global coordinates 
 *  with full geometry
 *
 *  $Date: 2007/10/23 23:16:30 $
 *  $Revision: 1.1 $
 *  \author M. Strang SUNY-Buffalo
 */

// framework & common header files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "DataFormats/Common/interface/Provenance.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//DQM services
#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// ecal calorimeter info
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/ESDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "CalibCalorimetry/EcalTrivialCondModules/interface/EcalTrivialConditionRetriever.h"

// hcal calorimeter info
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalDigi/interface/HFDataFrame.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"
//#include "Geometry/Records/interface/IdealGeometryRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

// silicon strip info
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"

// silicon pixel info
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"

// muon DT info
#include "DataFormats/DTDigi/interface/DTDigi.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"

// muon CSC Strip info
#include "DataFormats/CSCDigi/interface/CSCStripDigi.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"

// muon CSC Wire info
#include "DataFormats/CSCDigi/interface/CSCWireDigi.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"

// event info
//#include "SimDataFormats/ValidationFormats/interface/PValidationFormats.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

// helper files
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include <iostream>
#include <stdlib.h>
#include <string>
#include <memory>
#include <vector>
#include <map>

#include "TString.h"

class PGlobalDigi;

class GlobalDigisAnalyzer : public edm::EDAnalyzer
{

 public:
  typedef std::vector<float> FloatVector;
  typedef std::vector<double> DoubleVector;
  typedef std::vector<int> IntVector;
  typedef std::map<uint32_t,float,std::less<uint32_t> > MapType;

  explicit GlobalDigisAnalyzer(const edm::ParameterSet&);
  virtual ~GlobalDigisAnalyzer();
  virtual void beginJob(const edm::EventSetup&);
  virtual void endJob();  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  
 private:

  // production related methods
  void fillECal(const edm::Event&, const edm::EventSetup&);
  //void storeECal(PGlobalDigi&);
  void fillHCal(const edm::Event&, const edm::EventSetup&);
  //void storeHCal(PGlobalDigi&);
  void fillTrk(const edm::Event&, const edm::EventSetup&);
  //void storeTrk(PGlobalDigi&);
  void fillMuon(const edm::Event&, const edm::EventSetup&);
  //void storeMuon(PGlobalDigi&);  

  //void clear();

 private:

  //  parameter information
  std::string fName;
  int verbosity;
  int frequency;
  std::string label;
  bool getAllProvenances;
  bool printProvenanceInfo;

  DaqMonitorBEInterface *dbe;
  std::string outputfile;

  // Electromagnetic info
  // ECal info
 
  MonitorElement *mehEcaln[2];
  MonitorElement *mehEScaln;
  MonitorElement *mehEcalAEE[2];
  MonitorElement *mehEcalSHE[2];
  MonitorElement *mehEcalMaxPos[2];
  MonitorElement *mehEcalMultvAEE[2];
  MonitorElement *mehEcalSHEvAEESHE[2];
  MonitorElement *mehEScalADC[3];

 

  edm::InputTag ECalEBSrc_;
  edm::InputTag ECalEESrc_;
  edm::InputTag ECalESSrc_;

  std::map<int, double, std::less<int> > ECalgainConv_;
  double ECalbarrelADCtoGeV_;
  double ECalendcapADCtoGeV_;

  // HCal info

  MonitorElement *mehHcaln[4];
  MonitorElement *mehHcalAEE[4];
  MonitorElement *mehHcalSHE[4];
  MonitorElement *mehHcalAEESHE[4];
  MonitorElement *mehHcalSHEvAEE[4];

  edm::InputTag HCalSrc_;
  edm::InputTag HCalDigi_;

  // Tracker info
  // SiStrip
  
  MonitorElement *mehSiStripn[19];
  MonitorElement *mehSiStripADC[19];
  MonitorElement *mehSiStripStrip[19];

  edm::InputTag SiStripSrc_;

  // SiPxl

   MonitorElement *mehSiPixeln[7];
  MonitorElement *mehSiPixelADC[7];
  MonitorElement *mehSiPixelRow[7];
  MonitorElement *mehSiPixelCol[7];

  edm::InputTag SiPxlSrc_;

  // Muon info
  // DT

  MonitorElement *mehDtMuonn[4];
  MonitorElement *mehDtMuonLayer[4];
  MonitorElement *mehDtMuonTime[4];
  MonitorElement *mehDtMuonTimevLayer[4];

  edm::InputTag MuDTSrc_;

  // CSC Strip

  MonitorElement *mehCSCStripn;
  MonitorElement *mehCSCStripADC;
  MonitorElement *mehCSCWiren;
  MonitorElement *mehCSCWireTime;

  edm::InputTag MuCSCStripSrc_;
float theCSCStripPedestalSum;
  int theCSCStripPedestalCount;


  edm::InputTag MuCSCWireSrc_;

  // private statistics information
  unsigned int count;

}; // end class declaration

#endif

#ifndef GlobalHitMap
#define GlobalHitMap

// geometry mapping
static const int dTrk             = 1;
static const int sdPxlBrl         = 1;
static const int sdPxlFwd         = 2;
static const int sdSiTIB          = 3;
static const int sdSiTID          = 4;
static const int sdSiTOB          = 5;
static const int sdSiTEC          = 6;

static const int dMuon            = 2;
static const int sdMuonDT         = 1;
static const int sdMuonCSC        = 2;
static const int sdMuonRPC        = 3;
static const int sdMuonRPCRgnBrl  = 0;
static const int sdMuonRPCRgnFwdp = 1;
static const int sdMuonRPCRgnFwdn = -1;

static const int dEcal            = 3;
static const int sdEcalBrl        = 1;
static const int sdEcalFwd        = 2;
static const int sdEcalPS         = 3;
static const int sdEcalTT         = 4;
static const int sdEcalLPnD       = 5;

static const int dHcal            = 4;
static const int sdHcalEmpty      = 0;
static const int sdHcalBrl        = 1;
static const int sdHcalEC         = 2;
static const int sdHcalOut        = 3;
static const int sdHcalFwd        = 4;
static const int sdHcalTT         = 5;
static const int sdHcalCalib      = 6;
static const int sdHcalCompst     = 7;

#endif //PGlobalDigisAnalyzer_h
