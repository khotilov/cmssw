// -*- C++ -*-
//
// Package:    SiPixelCalibDigiProducer
// Class:      SiPixelCalibDigiProducer
// 
/**\class SiPixelCalibDigiProducer SiPixelCalibDigiProducer.cc CalibTracker/SiPixelCalibDigiProducer/src/SiPixelCalibDigiProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Wed Oct 31 15:28:52 CET 2007
// $Id: SiPixelCalibDigiProducer.cc,v 1.8 2008/01/23 12:04:04 fblekman Exp $
//
//


// system include files

// user include files
#include "FWCore/Framework/interface/EventSetup.h"






#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h" 

#include "CondFormats/SiPixelObjects/interface/SiPixelFrameConverter.h"
#include "CondFormats/SiPixelObjects/interface/ElectronicIndex.h"
#include "CondFormats/SiPixelObjects/interface/DetectorIndex.h"
#include "CondFormats/SiPixelObjects/interface/LocalPixel.h"

#include "SiPixelCalibDigiProducer.h"
#include <sstream>

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
SiPixelCalibDigiProducer::SiPixelCalibDigiProducer(const edm::ParameterSet& iConfig):
  src_(iConfig.getParameter<edm::InputTag>("src")),
  iEventCounter_(0),
  ignore_non_pattern_(iConfig.getParameter<bool>("ignoreNonPattern")),
  control_pattern_size_(iConfig.getParameter<bool>("checkPatternEachEvent")),
  conf_(iConfig),
  number_of_pixels_per_pattern_(0)

{
   //register your products
  produces< edm::DetSetVector<SiPixelCalibDigi> >();

}


SiPixelCalibDigiProducer::~SiPixelCalibDigiProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

/////////////////////////////////////////////
// function description:
// this function checks where/when the pattern changes so that one can fill the data from the temporary container into the event 
bool 
SiPixelCalibDigiProducer::store()
{
  //  std::cout << "in store() " << std::endl;
  if(iEventCounter_%pattern_repeat_==0){
    //    std::cout << "now at event " << iEventCounter_ <<" where we save the calibration information into the CMSSW digi";
    return 1;
  }
  else if(iEventCounter_==calib_->expectedTotalEvents())
    return 1;
  else
    return 0;
  return 1;
}
////////////////////////////////////////////////////////////////////
// function description:
// fill function, uses maps to keep track of which pixel is where in the local storage container. Called every event to fill and compress the digi data into calibdig format
void
SiPixelCalibDigiProducer::fill(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  std::cout << "in fill() " << std::endl;

  // figure out which calibration point we're on now..
  short icalibpoint = calib_->vcalIndexForEvent(iEventCounter_);
  edm::Handle< edm::DetSetVector<PixelDigi> > pixelDigis;
  iEvent.getByLabel( src_, pixelDigis );

    // loop over the data and store things
  edm::DetSetVector<PixelDigi>::const_iterator digiIter;
  for(digiIter=pixelDigis->begin(); digiIter!=pixelDigis->end(); ++digiIter){// ITERATOR OVER DET IDs
    uint32_t detid = digiIter->id;
    edm::DetSet<PixelDigi>::const_iterator ipix; // ITERATOR OVER DIGI DATA  

    for(ipix = digiIter->data.begin(); ipix!=digiIter->end(); ++ipix){

      // fill in the appropriate location of the temporary data container
      fillPixel(detid,ipix->row(),ipix->column(),icalibpoint,ipix->adc());
    }
  }
}
////////////////////////////////////////////////////
// function description:
// this is the function where we check the cabling map and see if we can assign a fed id to the det ID. 
// returns false if no fed <-> detid association was found
bool SiPixelCalibDigiProducer::checkFED(uint32_t detid){
  //  std::cout << "in checkFED" << std::endl;

  if(detid_to_fedid_[detid])
    return true;
  for(int fedid=0; fedid<=40; ++fedid){
    //    std::cout << " looking at fedid " << fedid << std::endl;
    SiPixelFrameConverter converter(theCablingMap_.product(),fedid);
    if(converter.hasDetUnit(detid)){
      detid_to_fedid_[detid]=fedid;
      edm::LogInfo("SiPixelCalibDigiProducer") << "matched detid " << detid << " to fed " << detid_to_fedid_[detid] << std::endl;
      return true;
    }
  }
  return false;
}

////////////////////////////////////////////////////
// function description:
// this is the function where we look in the maps to find the correct calibration digi container, after which the data is filled.
void SiPixelCalibDigiProducer::fillPixel(uint32_t detid, short row, short col, short ipoint, short adc){
  //  std::cout << " in fillpixel()" << std::endl;
 
  if(!checkFED(detid)){
    edm::LogError("SiPixelCalibDigiProducer") << " was unable to match detid " << detid << " to a FED!" << std::endl;
    return;
  }
  if(!checkPixel(detid,row,col)){
    return;
  }
  // now the check if the pixel exists and fill
  //
  pixelstruct temppixelworker;
  temppixelworker.first=detid;
  temppixelworker.second.first=row;
  temppixelworker.second.second=col;  
  std::map<pixelstruct,SiPixelCalibDigi>::const_iterator ipix = intermediate_data_.find(temppixelworker);
  
  if(ipix == intermediate_data_.end()){
    SiPixelCalibDigi tempdigi(calib_->nVCal());
    tempdigi.setrowcol(row,col);
    intermediate_data_[temppixelworker]=tempdigi;
  }

  intermediate_data_[temppixelworker].fill(ipoint,adc);
  return;

}
//////////////////////////////////////////////////////////////
// function description:
// this function cleans up after everything in the pattern is filled. This involves setting the content of the intermediate_data_ containers to zero and completely emptying the map 
void
SiPixelCalibDigiProducer::clear(){
  //  std::cout << "in clear() " << std::endl;
  // this is where we empty the containers so they can be re-filled
  // the idea: the detPixelMap_ container shrinks/expands as a function 
  // of the number of pixels looked at at one particular time... 
  // unlike the intermediate_data_ container which only expands when 
  // detPixelMap_ becomes bigger than intermedate_data_
  
  // shrink the detPixelMap_
  uint32_t tempsize = intermediate_data_.size();
  if(tempsize>number_of_pixels_per_pattern_){
    edm::LogError("SiPixelCalibDigiProducer") << "Number of pixels in pattern is now: " << tempsize<< ", size is was " << number_of_pixels_per_pattern_ << std::endl;
    number_of_pixels_per_pattern_=tempsize;
  }

  intermediate_data_.erase(intermediate_data_.begin(),intermediate_data_.end());
  intermediate_data_.clear();
}

////////////////////////////////////////////////////
// function description:
// This method gets the pattern from the calib_ (SiPixelCalibConfiguration) object and fills a vector of pairs that is easier to check 
void 
SiPixelCalibDigiProducer::setPattern(){
  //  std::cout << "in setPattern()" << std::endl;
  uint32_t patternnumber = abs(iEventCounter_-1)/pattern_repeat_;
  uint32_t rowpatternnumber = patternnumber/calib_->nColumnPatterns();
  uint32_t colpatternnumber = patternnumber%calib_->nColumnPatterns();
  edm::LogInfo("SiPixelCalibDigiProducer") << " rowpatternnumbers = " << rowpatternnumber << " " << colpatternnumber << " " << patternnumber << std::endl;
  // update currentpattern_
  std::vector<short> calibcols = calib_->getColumnPattern();
  std::vector<short> calibrows = calib_->getRowPattern();
  std::vector<short> temprowvals(0);
  std::vector<short> tempcolvals(0);
  uint32_t nminuscol=0;
  uint32_t nminusrow=0;
  uint32_t npatterns=0;
  for(uint32_t icol=0; icol<calibcols.size(); icol++){
    if(calibcols[icol]==-1){
      nminuscol++;
    }
    else if(nminuscol==colpatternnumber){
      //std::cout << "col " << calibcols[icol] << std::endl;
      short val=calibcols[icol];
      tempcolvals.push_back(val);
    }
    else if (nminuscol> colpatternnumber)
      break;
  }
  for(uint32_t irow=0; irow<calibrows.size(); irow++){
    // std::cout << "row " << irow <<" "<< nminusrow<<" "  << calibrows[irow] << std::endl;
    if(calibrows[irow]==-1)
      nminusrow++;
    else if(nminusrow==rowpatternnumber){
      short val=calibrows[irow];
      temprowvals.push_back(val);
    }
    else if(nminusrow>rowpatternnumber)
      break;
  }
  //now clean up the currentpattern_;
  while(currentpattern_.size()>temprowvals.size()*tempcolvals.size()){
    currentpattern_.erase(currentpattern_.end());
  }
  for(uint32_t irow=0; irow<temprowvals.size(); irow++){
    for(uint32_t icol=0; icol<tempcolvals.size(); icol++){
      std::pair<short,short> pattern(temprowvals[irow],tempcolvals[icol]);
      npatterns++;
      if(npatterns>currentpattern_.size())
	currentpattern_.push_back(pattern);
      else
	currentpattern_[npatterns-1]=pattern;
    }
  }
  
//   std::cout << "summary of created patterns: " ;
//   for(uint32_t i=0; i<currentpattern_.size(); ++i){
//     if(i!=0)
//       std::cout << " - ";
//     std::cout << currentpattern_[i].first << ","<< currentpattern_[i].second ;
   
//   }
//   std::cout << std::endl;

}

////////////////////////////////////////////
// function description:
// produce method. This is the main loop method
void
SiPixelCalibDigiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  std::cout <<"in produce() " << std::endl;
  using namespace edm;
  iEventCounter_++;
  if(iEventCounter_%pattern_repeat_==1)
    setPattern();
  
  fill(iEvent,iSetup); // fill method where the actual looping over the digis is done.

  std::auto_ptr<edm::DetSetVector<SiPixelCalibDigi> > pOut(new edm::DetSetVector<SiPixelCalibDigi>);
 
  // copy things over into pOut if necessary (this is only once per pattern)
  if(store()){
    //    std::cout << "in loop" << std::endl;
    for(std::map<pixelstruct,SiPixelCalibDigi>::const_iterator idet=intermediate_data_.begin(); idet!=intermediate_data_.end();++idet){
      if(!control_pattern_size_)
	if(! checkPixel(idet->first.first,idet->first.second.first,idet->first.second.second))
	  continue;
      uint32_t detid=idet->first.first;
      SiPixelCalibDigi tempdigi=idet->second;
      edm::DetSet<SiPixelCalibDigi> & detSet = pOut->find_or_insert(detid);
      detSet.data.push_back(tempdigi);
    }
    edm::LogInfo("INFO") << "now filling event " << iEventCounter_ << " as pixel pattern changes every " <<  pattern_repeat_ << " events..." << std::endl;
    clear();
  }
  iEvent.put(pOut);
}
//-----------------------------------------------
//  method to check that the pixels are actually valid...
bool SiPixelCalibDigiProducer::checkPixel(uint32_t detid, short row, short col){

  if(!control_pattern_size_ && !store())
    return true;
  
  if( !ignore_non_pattern_ )
    return true;
  
  
  //  std::cout << "Event" << iEventCounter_ << ",now in checkpixel() " << std::endl;
  if(currentpattern_.size()==0)
    setPattern();
  //  uint32_t iroc;
  uint32_t fedid = detid_to_fedid_[detid];

  SiPixelFrameConverter formatter(theCablingMap_.product(),fedid);
  sipixelobjects::DetectorIndex detector ={detid, row, col};
  sipixelobjects::ElectronicIndex cabling;
  
  formatter.toCabling(cabling,detector);
  // cabling should now contain cabling.roc and cabling.dcol  and cabling.pxid

  // however, the coordinates now need to be converted from dcl, pxid to the row,col coordinates used in the calibration info
  sipixelobjects::LocalPixel::DcolPxid loc;
  loc.dcol = cabling.dcol;
  loc.pxid = cabling.pxid;
  sipixelobjects::LocalPixel locpixel(loc);
  currentpair_.first = locpixel.rocRow();
  currentpair_.second = locpixel.rocCol();

  for(uint32_t i=0; i<currentpattern_.size(); ++i){
    //    std::cout << "found pair " << currentpair_.first << "," << currentpair_.second << " calib " << currentpattern_[i].first << ","<< currentpattern_[i].second << " input " << row << "," << col << std::endl;
    if(currentpair_==currentpattern_[i]){
      return true;
    }
  }
  std::ostringstream errorlog;
  errorlog <<  "DETID "<<detid<<", row, col (offline)="<<row<<","<<col<<" row, col (ROC) ="<<currentpair_.first<<","<< currentpair_.second<< " found no match in list of patterns: " ;
  for(uint32_t i=0; i<currentpattern_.size(); ++i){
    if(i!=0 && i!=currentpattern_.size()-1)
      errorlog<<" ";
    errorlog<<"(";
    errorlog<<currentpattern_[i].first;
    errorlog<<",";
    errorlog<<currentpattern_[i].second;
    errorlog<<")";
  }
  edm::LogError("ERROR") << errorlog.str() << std::endl;
  return false;
}


// ------------ method called once each job just before starting event loop  ------------
void 
SiPixelCalibDigiProducer::beginJob(const edm::EventSetup& iSetup)
{
  iSetup.get<SiPixelCalibConfigurationRcd>().get(calib_);
  pattern_repeat_=calib_->getNTriggers()*calib_->nVCal();
  edm::LogInfo("SiPixelCalibDigiProducer") << "got calibration configuration, pattern repeat is " << pattern_repeat_ << ", number of VCal points is " << calib_->nVCal() << std::endl;
  iSetup.get<TrackerDigiGeometryRecord>().get( theGeometry_ );
  iSetup.get<SiPixelFedCablingMapRcd>().get(theCablingMap_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiPixelCalibDigiProducer::endJob() {
}
DEFINE_FWK_MODULE(SiPixelCalibDigiProducer);
