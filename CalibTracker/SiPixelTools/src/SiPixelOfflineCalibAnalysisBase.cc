// -*- C++ -*-
//
// Package:    SiPixelOfflineCalibAnalysisBase
// Class:      SiPixelOfflineCalibAnalysisBase
// 
/**\class SiPixelOfflineCalibAnalysisBase SiPixelOfflineCalibAnalysisBase.cc CalibTracker/SiPixelTools/src/SiPixelOfflineCalibAnalysisBase.cc

 Description: Base class for Pixel calibrations

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Evan Klose Friis
//    additions by:  Freya Blekman
//         Created:  Tue Nov  6 17:27:19 CET 2007
// $Id: SiPixelOfflineCalibAnalysisBase.cc,v 1.1 2007/11/19 11:38:25 fblekman Exp $
//
//

#include "CalibTracker/SiPixelTools/interface/SiPixelOfflineCalibAnalysisBase.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h" 

TF1* SiPixelOfflineCalibAnalysisBase::fitFunction_ = NULL;
std::vector<short>  SiPixelOfflineCalibAnalysisBase::vCalValues_(0);
// constructors and destructor
//
SiPixelOfflineCalibAnalysisBase::SiPixelOfflineCalibAnalysisBase(const edm::ParameterSet& iConfig)
{
   siPixelCalibDigiProducer_ = iConfig.getParameter<edm::InputTag>("DetSetVectorSiPixelCalibDigiTag");
   outputFileName_ = iConfig.getParameter<std::string>("outputFileName");
   daqBE_ = &*edm::Service<DaqMonitorBEInterface>();
   folderMaker_ = new SiPixelFolderOrganizer();
}

SiPixelOfflineCalibAnalysisBase::SiPixelOfflineCalibAnalysisBase()
{
   throw cms::Exception("") << "ERROR: Classes derived from SiPixelOfflineCalibAnalysisBase must call SiPixelOfflineCalibAnalysisBase::SiPixelOfflineCalibAnalysisBase(const edm::ParameterSet& iConfig) from their constructor." << std::endl;
}

SiPixelOfflineCalibAnalysisBase::~SiPixelOfflineCalibAnalysisBase()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiPixelOfflineCalibAnalysisBase::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<DetSetVector<SiPixelCalibDigi> > thePlaquettes;
   iEvent.getByLabel(siPixelCalibDigiProducer_, thePlaquettes);

   DetSetVector<SiPixelCalibDigi>::const_iterator digiIter;

   //loop over the plaquettes pulsed in this pattern 
   for(digiIter=thePlaquettes->begin(); digiIter!=thePlaquettes->end(); ++digiIter) 
   {
      uint32_t detId = digiIter->id;
      //check to see if this detID has not been encountered.  If not, run the newDetID (pure virtual) function
      if (detIdsEntered_.find(detId) == detIdsEntered_.end()) 
      {
	 detIdsEntered_.insert(std::make_pair(detId, 0));
	 detIdNames_.insert(std::make_pair(detId, translateDetIdToString(detId)));
	 newDetID(detId);
      }
      DetSet<SiPixelCalibDigi>::const_iterator ipix;
      //loop over pixels pulsed in the current plaquette
      for(ipix=digiIter->data.begin(); ipix!=digiIter->end(); ++ipix)
      {
	 //called derived function to fit this curve
	this->doFits(detId, ipix);
      }
   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
SiPixelOfflineCalibAnalysisBase::beginJob(const edm::EventSetup& iSetup)
{
   //load the calibration information from the database
   iSetup.get<SiPixelCalibConfigurationRcd>().get(calib_);
   iSetup.get<TrackerDigiGeometryRecord>().get( geom_ );

   calibrationMode_ 	= calib_->getCalibrationMode();
   nTriggers_ 		= calib_->getNTriggers();
   vCalValues_		= calib_->getVCalValues();

   edm::LogInfo("SiPixelOfflineCalibAnalysisBase") << "Calibration file loaded. Mode: " << calibrationMode_ << " nTriggers: " << nTriggers_ << " Vcal steps: " << vCalValues_.size() << std::endl;

   //call calibrationSetup virtual function
   this->calibrationSetup(iSetup);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiPixelOfflineCalibAnalysisBase::endJob() {
   edm::LogInfo("SiPixelOfflineCalibAnalysisBase") << "Running end job... output file name is: " << outputFileName_;
   if (!outputFileName_.empty()) 
   {
      edm::LogInfo("SiPixelOfflineCalibAnalysisBase") << "Writing ROOT file to: " << outputFileName_ << std::endl;
      daqBE_->save(outputFileName_);
   }

}

// ------------ helper functions ---------------------------------------------------------

const std::vector<short>* 
SiPixelOfflineCalibAnalysisBase::getVcalValues()
{
   return &vCalValues_;
}

std::string
SiPixelOfflineCalibAnalysisBase::translateDetIdToString(uint32_t detid)
{
   std::map<uint32_t, std::string>::iterator detNameIter = detIdNames_.find(detid);
   if (detNameIter != detIdNames_.end()) {
      return detNameIter->second;
   }
   std::string output = "DetID translation error";
   DetId detId(detid);
   uint32_t detSubId = detId.subdetId();
   if (detSubId > 2 || detSubId < 1)
   {
      edm::LogError("SiPixelOfflineCalibAnalysisBase") << "ERROR: Expected a pixel detector ID (1 - barrel, 2 - forward) but got " << detSubId << std::endl;
      return output;
   }
   if (detSubId == 2)  	//FPIX
   {
      PixelEndcapName nameworker(detid);
      output = nameworker.name();
   } else 		//BPIX
   {
      PixelBarrelName nameworker(detid);
      output = nameworker.name();
   }
   detIdNames_.insert(std::make_pair(detid, output));
   return output;
}

MonitorElement* SiPixelOfflineCalibAnalysisBase::bookDQMHistogram1D(std::string name, std::string title, int nchX, double lowX, double highX)  
{
   return daqBE_->book1D(name, title, nchX, lowX, highX);
}

MonitorElement* SiPixelOfflineCalibAnalysisBase::bookDQMHistogram1D(std::string name, std::string title, int nchX, float *xbinsize)
{
   return daqBE_->book1D(name, title, nchX, xbinsize);
}

MonitorElement* SiPixelOfflineCalibAnalysisBase::bookDQMHistogram2D(std::string name, std::string title, int nchX, double lowX, double highX, int nchY, double lowY, double highY)
{
   return daqBE_->book2D(name, title, nchX, lowX, highX, nchY, lowY, highY);
}

MonitorElement* SiPixelOfflineCalibAnalysisBase::bookDQMHistoPlaquetteSummary2D(std::string name,std::string title,uint32_t detid){
  DetId detId(detid);
  const TrackerGeometry &theTracker(*geom_);
  const PixelGeomDetUnit *theGeomDet = dynamic_cast<const PixelGeomDetUnit*> ( theTracker.idToDet(detId) ); 
  int maxcol = theGeomDet->specificTopology().ncolumns();
  int maxrow = theGeomDet->specificTopology().nrows();
  
  return daqBE_->book2D(name,title,maxcol,0,maxcol,maxrow,0,maxrow);
}

bool
SiPixelOfflineCalibAnalysisBase::setDQMDirectory(std::string dirName)
{
   daqBE_->setCurrentFolder(dirName);
   return daqBE_->dirExists(dirName);
}

bool
SiPixelOfflineCalibAnalysisBase::setDQMDirectory(uint32_t detID)
{
   return folderMaker_->setModuleFolder(detID);
}


// ------------ virtual functions ------------------------------------------------

bool
SiPixelOfflineCalibAnalysisBase::doFits(uint32_t detid, std::vector<SiPixelCalibDigi>::const_iterator ipix)
{
  short row=ipix->row();
  short col=ipix->col();
  std::vector<short> nentries = ipix->getnentries();
  std::vector<uint32_t> sum = ipix->getsum();
  std::vector<uint32_t> sumquares = ipix->getsumsquares();
   //do nothing
   //return false;
   //
   //DEBUG
   std::cout << "Row: " << row << "   Col: " << col << std::endl;
   for (unsigned int i = 0; i < sum.size(); i++)
   {
      std::cout << sum[i] << "  ";
   }
   std::cout << std::endl;
   return false;

}

void 
SiPixelOfflineCalibAnalysisBase::calibrationSetup(const edm::EventSetup&)
{
   //do nothing
}

void 
SiPixelOfflineCalibAnalysisBase::newDetID(uint32_t detid)
{
   //do nothing
   edm::LogInfo("SiPixelOfflineCalibAnalysisBase") << "SiPixelOfflineCalibAnalysisBase - Found new DetID: " << detid << "  Name: " << detIdNames_[detid];
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelOfflineCalibAnalysisBase);
