// -*- C++ -*-
//
// Package:    TrackerOfflineValidation
// Class:      TrackerOfflineValidation
// 
/**\class TrackerOfflineValidation TrackerOfflineValidation.cc Alignment/Validator/src/TrackerOfflineValidation.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Erik Butz
//         Created:  Tue Dec 11 14:03:05 CET 2007
// $Id: TrackerOfflineValidation.cc,v 1.26 2009/03/24 16:13:54 jdraeger Exp $
//
//


// system include files
#include <memory>
#include <map>
#include <sstream>
#include <math.h>
#include <utility>
#include <vector>
#include <iostream>

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"
#include "Alignment/OfflineValidation/interface/TkOffTreeVariables.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"
#include "Alignment/CommonAlignment/interface/AlignableComposite.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Alignment/CommonAlignment/interface/Utilities.h"
#include "Alignment/CommonAlignment/interface/AlignableObjectId.h"
#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"

#include "DataFormats/Math/interface/deltaPhi.h"
//
// class declaration
//

class TrackerOfflineValidation : public edm::EDAnalyzer {
public:
  explicit TrackerOfflineValidation(const edm::ParameterSet&);
  ~TrackerOfflineValidation();
  
  enum HistogrammType { XResidual, NormXResidual, 
			XprimeResidual, NormXprimeResidual, 
			YprimeResidual, NormYprimeResidual};
  
private:

  
  struct ModuleHistos{
    ModuleHistos() :  ResHisto(), NormResHisto(), ResXprimeHisto(), NormResXprimeHisto(), 
		      ResYprimeHisto(), NormResYprimeHisto() {} 
    TH1* ResHisto;
    TH1* NormResHisto;
    TH1* ResXprimeHisto;
    TH1* NormResXprimeHisto;
    TH1* ResYprimeHisto;
    TH1* NormResYprimeHisto;
  };


  // container struct to organize collection of histogramms during endJob
  struct SummaryContainer{
    SummaryContainer() : sumXResiduals_(), summaryXResiduals_(), 
			 sumNormXResiduals_(), summaryNormXResiduals_(),
			 sumYResiduals_(), summaryYResiduals_() ,
			 sumNormYResiduals_(), summaryNormYResiduals_() {}
    
    TH1* sumXResiduals_;
    TH1* summaryXResiduals_;
    TH1* sumNormXResiduals_;
    TH1* summaryNormXResiduals_;
    TH1* sumYResiduals_;
    TH1* summaryYResiduals_;
    TH1* sumNormYResiduals_;
    TH1* summaryNormYResiduals_;
  };


  // 
  // ------------- private member function -------------
  // 
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual void checkBookHists(const edm::EventSetup &setup);



  void bookGlobalHists(TFileDirectory &tfd);
  void bookDirHists(TFileDirectory &tfd, const Alignable& ali, const AlignableObjectId &aliobjid);
  void bookHists(TFileDirectory &tfd, const Alignable& ali, align::StructureType type, int i, 
		 const AlignableObjectId &aliobjid);
 
  void collateSummaryHists( TFileDirectory &tfd, const Alignable& ali, int i, 
			    const AlignableObjectId &aliobjid, 
			    std::vector<TrackerOfflineValidation::SummaryContainer > &v_levelProfiles);
  
  void fillTree(TTree &tree, const std::map<int, TrackerOfflineValidation::ModuleHistos> &moduleHist_, 
		TkOffTreeVariables &treeMem, const TrackerGeometry &tkgeom );
  
  TrackerOfflineValidation::SummaryContainer bookSummaryHists(TFileDirectory &tfd, 
							      const Alignable& ali, 
							      align::StructureType type, int i, 
							      const AlignableObjectId &aliobjid); 

  ModuleHistos& getHistStructFromMap(const DetId& detid); 

  bool isBarrel(uint32_t subDetId);
  bool isEndCap(uint32_t subDetId);
  bool isPixel(uint32_t subDetId);
  bool isDetOrDetUnit(align::StructureType type);

  TH1* bookTH1F(bool isTransient, TFileDirectory& tfd, const char* histName, const char* histTitle, 
		int nBinsX, double lowX, double highX);

  void getBinning(uint32_t subDetId, TrackerOfflineValidation::HistogrammType residualtype, 
		  int &nBinsX, double &lowerBoundX, double &upperBoundX);

  void summarizeBinInContainer(int bin, SummaryContainer &targetContainer, 
			       SummaryContainer &sourceContainer);

  void summarizeBinInContainer(int bin, uint32_t subDetId, SummaryContainer &targetContainer, 
			       ModuleHistos &sourceContainer);

  void setSummaryBin(int bin, TH1* targetHist, TH1* sourceHist);
    
  float Fwhm(const TH1* hist) const;
  std::pair<float,float> fitResiduals(TH1 *hist) const; //, float meantmp, float rmstmp);
  float getMedian( const TH1 *hist) const; 
 // From MillePedeAlignmentMonitor: Get Index for Arbitary vector<class> by name
  template <class OBJECT_TYPE>  
  int GetIndex(const std::vector<OBJECT_TYPE*> &vec, const TString &name);

  // ---------- member data ---------------------------


  const edm::ParameterSet parset_;
  edm::ESHandle<TrackerGeometry> tkGeom_;
  const TrackerGeometry *bareTkGeomPtr_; // ugly hack to book hists only once, but check 

  
  // parameters from cfg to steer
  bool lCoorHistOn_;
  bool moduleLevelHistsTransient_;
  bool overlappOn_;
  bool stripYResiduals_;
  bool useFwhm_;
  bool useFit_;
  bool useOverflowForRMS_;
  std::map< std::pair<uint32_t, uint32_t >, TH1*> hOverlappResidual;

  // a vector to keep track which pointers should be deleted at the very end
  std::vector<TH1*> vDeleteObjects_;

  // 
  std::vector<TH1*> vTrackHistos_;
  std::vector<TProfile*> vTrackProfiles_;
  std::vector<TH2*> vTrack2DHistos_;
  
  std::map<int,TrackerOfflineValidation::ModuleHistos> mPxbResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mPxeResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mTibResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mTidResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mTobResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mTecResiduals_;



};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
template <class OBJECT_TYPE>  
int TrackerOfflineValidation::GetIndex(const std::vector<OBJECT_TYPE*> &vec, const TString &name)
{
  int result = 0;
  for (typename std::vector<OBJECT_TYPE*>::const_iterator iter = vec.begin(), iterEnd = vec.end();
       iter != iterEnd; ++iter, ++result) {
    if (*iter && (*iter)->GetName() == name) return result;
  }
  edm::LogError("Alignment") << "@SUB=TrackerOfflineValidation::GetIndex" << " could not find " << name;
  return -1;
}
//
// constructors and destructor
//
TrackerOfflineValidation::TrackerOfflineValidation(const edm::ParameterSet& iConfig)
  : parset_(iConfig), bareTkGeomPtr_(0), lCoorHistOn_(parset_.getParameter<bool>("localCoorHistosOn")),
    moduleLevelHistsTransient_(parset_.getParameter<bool>("moduleLevelHistsTransient")),
    overlappOn_(parset_.getParameter<bool>("overlappOn")), 
    stripYResiduals_(parset_.getParameter<bool>("stripYResiduals")), 
    useFwhm_(parset_.getParameter<bool>("useFwhm")),
    useFit_(parset_.getParameter<bool>("useFit")),
    useOverflowForRMS_(parset_.getParameter<bool>("useOverflowForRMS"))
  
{
   //now do what ever initialization is needed
}


TrackerOfflineValidation::~TrackerOfflineValidation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  for( std::vector<TH1*>::const_iterator it = vDeleteObjects_.begin(), itEnd = vDeleteObjects_.end(); 
       it != itEnd; ++it) delete *it;
    
}


//
// member functions
//


// ------------ method called once each job just before starting event loop  ------------
void
TrackerOfflineValidation::checkBookHists(const edm::EventSetup &es)
{
  es.get<TrackerDigiGeometryRecord>().get( tkGeom_ );
  const TrackerGeometry *newBareTkGeomPtr = &(*tkGeom_);
  if (newBareTkGeomPtr == bareTkGeomPtr_) return; // already booked hists, nothing changed

  if (!bareTkGeomPtr_) { // pointer not yet set: called the first time => book hists
    edm::Service<TFileService> fs;    
    AlignableObjectId aliobjid;
    
    // construct alignable tracker to get access to alignable hierarchy 
    AlignableTracker aliTracker(&(*tkGeom_));
    
    edm::LogInfo("TrackerOfflineValidation") << "There are " << newBareTkGeomPtr->detIds().size()
					     << " detUnits in the Geometry record";
    
    //
    // Book Histogramms for global track quantities
    TFileDirectory trackglobal = fs->mkdir("GlobalTrackVariables");  
    this->bookGlobalHists(trackglobal);
    
    // recursively book histogramms on lowest level
//     this->bookDirHists(static_cast<TFileDirectory&>(*fs), aliTracker, aliobjid);  
    this->bookDirHists(*fs, aliTracker, aliobjid);  
  } else { // histograms booked, but changed TrackerGeometry?
    edm::LogWarning("GeometryChange") << "@SUB=checkBookHists"
				      << "TrackerGeometry changed, but will not re-book hists!";
  }

  bareTkGeomPtr_ = newBareTkGeomPtr;
}


void 
TrackerOfflineValidation::bookGlobalHists(TFileDirectory &tfd )
{

  vTrackHistos_.push_back(tfd.make<TH1F>("h_tracketa",
					 "Track #eta;#eta_{Track};Number of Tracks",
					 90,-3.,3.));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_curvature",
					 "Curvature #kappa;#kappa_{Track};Number of Tracks",
					 100,-.05,.05));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_curvature_pos",
					 "Curvature |#kappa| Positive Tracks;|#kappa_{pos Track}|;Number of Tracks",
					 100,.0,.05));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_curvature_neg",
					 "Curvature |#kappa| Negative Tracks;|#kappa_{neg Track}|;Number of Tracks",
					 100,.0,.05));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_diff_curvature",
					 "Curvature |#kappa| Tracks Difference;|#kappa_{Track}|;# Pos Tracks - # Neg Tracks",
					 100,.0,.05));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_chi2",
					 "#chi^{2};#chi^{2}_{Track};Number of Tracks",
					 500,-0.01,500.));	       
  vTrackHistos_.push_back(tfd.make<TH1F>("h_normchi2",
					 "#chi^{2}/ndof;#chi^{2}/ndof;Number of Tracks",
					 100,-0.01,10.));     
  vTrackHistos_.push_back(tfd.make<TH1F>("h_pt",
					 "p_{T}^{track};p_{T}^{track} [GeV];Number of Tracks",
					 100,0.,2500));           
  vTrackHistos_.push_back(tfd.make<TH1F>("h_ptResolution",
					 "#delta{p_{T}/p_{T}^{track}};#delta_{p_{T}/p_{T}^{track}};Number of Tracks",
					 100,0.,0.5));           

  vTrackProfiles_.push_back(tfd.make<TProfile>("p_d0_vs_phi",
					       "Transverse Impact Parameter vs. #phi;#phi_{Track};#LT d_{0} #GT [cm]",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_dz_vs_phi",
					       "Longitudinal Impact Parameter vs. #phi;#phi_{Track};#LT d_{z} #GT [cm]",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_d0_vs_eta",
					       "Transverse Impact Parameter vs. #eta;#eta_{Track};#LT d_{0} #GT [cm]",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_dz_vs_eta",
					       "Longitudinal Impact Parameter vs. #eta;#eta_{Track};#LT d_{z} #GT [cm]",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_chi2_vs_phi",
					       "#chi^{2} vs. #phi;#phi_{Track};#LT #chi^{2} #GT",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_normchi2_vs_phi",
					       "#chi^{2}/ndof vs. #phi;#phi_{Track};#LT #chi^{2}/ndof #GT",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_chi2_vs_eta",
					       "#chi^{2} vs. #eta;#eta_{Track};#LT #chi^{2} #GT",
					       100,-3.15,3.15));  
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_normchi2_vs_eta",
					       "#chi^{2}/ndof vs. #eta;#eta_{Track};#LT #chi^{2}/ndof #GT",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_kappa_vs_phi",
					       "#kappa vs. #phi;#phi_{Track};#kappa",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_kappa_vs_eta",
					       "#kappa vs. #eta;#eta_{Track};#kappa",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_ptResolution_vs_phi",
					       "#delta_{p_{T}}/p_{T}^{track};#phi^{track};#delta_{p_{T}}/p_{T}^{track}",
					       100, -3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_ptResolution_vs_eta",
					       "#delta_{p_{T}}/p_{T}^{track};#eta^{track};#delta_{p_{T}}/p_{T}^{track}",
					       100, -3.15,3.15));


  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_d0_vs_phi",
					   "Transverse Impact Parameter vs. #phi;#phi_{Track};d_{0} [cm]",
					   100, -3.15, 3.15, 100,-1.,1.) );
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_dz_vs_phi",
					   "Longitudinal Impact Parameter vs. #phi;#phi_{Track};d_{z} [cm]",
					   100, -3.15, 3.15, 100,-100.,100.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_d0_vs_eta",
					   "Transverse Impact Parameter vs. #eta;#eta_{Track};d_{0} [cm]",
					   100, -3.15, 3.15, 100,-1.,1.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_dz_vs_eta",
					   "Longitudinal Impact Parameter vs. #eta;#eta_{Track};d_{z} [cm]",
					   100, -3.15, 3.15, 100,-100.,100.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_chi2_vs_phi",
					   "#chi^{2} vs. #phi;#phi_{Track};#chi^{2}",
					   100, -3.15, 3.15, 500, 0., 500.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_normchi2_vs_phi",
					   "#chi^{2}/ndof vs. #phi;#phi_{Track};#chi^{2}/ndof",
					   100, -3.15, 3.15, 100, 0., 10.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_chi2_vs_eta",
					   "#chi^{2} vs. #eta;#eta_{Track};#chi^{2}",
					   100, -3.15, 3.15, 500, 0., 500.));  
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_normchi2_vs_eta",
					   "#chi^{2}/ndof vs. #eta;#eta_{Track};#chi^{2}/ndof",
					   100,-3.15,3.15, 100, 0., 10.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_kappa_vs_phi",
					   "#kappa vs. #phi;#phi_{Track};#kappa",
					   100,-3.15,3.15, 100, .0,.05));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_kappa_vs_eta",
					   "#kappa vs. #eta;#eta_{Track};#kappa",
					   100,-3.15,3.15, 100, .0,.05));
 

}


void
TrackerOfflineValidation::bookDirHists( TFileDirectory &tfd, const Alignable& ali, const AlignableObjectId &aliobjid)
{
  std::vector<Alignable*> alivec(ali.components());
  for(int i=0, iEnd = ali.components().size();i < iEnd; ++i) {
    std::string structurename  = aliobjid.typeToName((alivec)[i]->alignableObjectId());
    LogDebug("TrackerOfflineValidation") << "StructureName = " << structurename;
    std::stringstream dirname;
    dirname << structurename;
    // add no suffix counter to Strip and Pixel, just aesthetics
    if (structurename != "Strip" && structurename != "Pixel") dirname << "_" << i+1;

    if (structurename.find("Endcap",0) != std::string::npos )   {
      TFileDirectory f = tfd.mkdir((dirname.str()).c_str());
      bookHists(f, *(alivec)[i], ali.alignableObjectId() , i, aliobjid);
      bookDirHists( f, *(alivec)[i], aliobjid);
    } else if( !(this->isDetOrDetUnit( (alivec)[i]->alignableObjectId()) )
	      || alivec[i]->components().size() > 1) {      
      TFileDirectory f = tfd.mkdir((dirname.str()).c_str());
      bookHists(tfd, *(alivec)[i], ali.alignableObjectId() , i, aliobjid);
      bookDirHists( f, *(alivec)[i], aliobjid);
    } else {
      bookHists(tfd, *(alivec)[i], ali.alignableObjectId() , i, aliobjid);
    }
  }
}




void 
TrackerOfflineValidation::bookHists(TFileDirectory &tfd, const Alignable& ali, align::StructureType type, int i, const AlignableObjectId &aliobjid)
{

  TrackerAlignableId aliid;
  const DetId id = ali.id();

  // comparing subdetandlayer to subdetIds gives a warning at compile time
  // -> subdetandlayer could also be pair<uint,uint> but this has to be adapted
  // in AlignableObjId 
  std::pair<int,int> subdetandlayer = aliid.typeAndLayerFromDetId(id);

  align::StructureType subtype = align::invalid;
  
  // are we on or just above det, detunit level respectively?
  if (type == align::AlignableDetUnit )subtype = type;
  else if( this->isDetOrDetUnit(ali.alignableObjectId()) ) subtype = ali.alignableObjectId();
  
  // construct histogramm title and name
  std::stringstream histoname, histotitle, normhistoname, normhistotitle, 
    xprimehistoname, xprimehistotitle, normxprimehistoname, normxprimehistotitle,
    yprimehistoname, yprimehistotitle, normyprimehistoname, normyprimehistotitle;
  
  std::string wheel_or_layer;

  if( this->isEndCap(static_cast<uint32_t>(subdetandlayer.first)) ) wheel_or_layer = "_wheel_";
  else if ( this->isBarrel(static_cast<uint32_t>(subdetandlayer.first)) ) wheel_or_layer = "_layer_";
  else edm::LogWarning("TrackerOfflineValidation") << "@SUB=TrackerOfflineValidation::bookHists" 
						   << "Unknown subdetid: " <<  subdetandlayer.first;     
  
  
  histoname << "h_residuals_subdet_" << subdetandlayer.first 
	    << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  xprimehistoname << "h_xprime_residuals_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  yprimehistoname << "h_yprime_residuals_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  normhistoname << "h_normresiduals_subdet_" << subdetandlayer.first 
		<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  normxprimehistoname << "h_normxprimeresiduals_subdet_" << subdetandlayer.first 
		      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  normyprimehistoname << "h_normyprimeresiduals_subdet_" << subdetandlayer.first 
		      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  histotitle << "Residual for module " << id.rawId() << ";x_{pred} - x_{rec} [cm]";
  normhistotitle << "Normalized Residual for module " << id.rawId() << ";x_{pred} - x_{rec}/#sigma";
  xprimehistotitle << "X' Residual for module " << id.rawId() << ";(x_{pred} - x_{rec})' [cm]";
  normxprimehistotitle << "Normalized X' Residual for module " << id.rawId() << ";(x_{pred} - x_{rec})'/#sigma";
  yprimehistotitle << "Y' Residual for module " << id.rawId() << ";(y_{pred} - y_{rec})' [cm]";
  normyprimehistotitle << "Normalized Y' Residual for module " << id.rawId() << ";(y_{pred} - y_{rec})'/#sigma";
  
  
  if( this->isDetOrDetUnit( subtype ) ) {
    ModuleHistos &histStruct = this->getHistStructFromMap(id);
    int nbins = 0;
    double xmin = 0., xmax = 0.;

    // decide via cfg if hists in local coordinates should be booked 
    if(lCoorHistOn_) {
      this->getBinning(id.subdetId(), XResidual, nbins, xmin, xmax);
      histStruct.ResHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					   histoname.str().c_str(),histotitle.str().c_str(),		     
					   nbins, xmin, xmax);
      this->getBinning(id.subdetId(), NormXResidual, nbins, xmin, xmax);
      histStruct.NormResHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd,
					       normhistoname.str().c_str(),normhistotitle.str().c_str(),
					       nbins, xmin, xmax);
    } 
    this->getBinning(id.subdetId(), XprimeResidual, nbins, xmin, xmax);
    histStruct.ResXprimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					       xprimehistoname.str().c_str(),xprimehistotitle.str().c_str(),
					       nbins, xmin, xmax);
    this->getBinning(id.subdetId(), NormXprimeResidual, nbins, xmin, xmax);
    histStruct.NormResXprimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						   normxprimehistoname.str().c_str(),normxprimehistotitle.str().c_str(),
						   nbins, xmin, xmax);

    if( this->isPixel(subdetandlayer.first) || stripYResiduals_ ) {
      this->getBinning(id.subdetId(), YprimeResidual, nbins, xmin, xmax);
      histStruct.ResYprimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd,
						 yprimehistoname.str().c_str(),yprimehistotitle.str().c_str(),
						 nbins, xmin, xmax);
      this->getBinning(id.subdetId(), NormYprimeResidual, nbins, xmin, xmax);
      histStruct.NormResYprimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						     normyprimehistoname.str().c_str(),normyprimehistotitle.str().c_str(),
						     nbins, xmin, xmax);
    }

  }
  
}


TH1* TrackerOfflineValidation::bookTH1F(bool isTransient, TFileDirectory& tfd, const char* histName, const char* histTitle, 
		int nBinsX, double lowX, double highX)
{
  if(isTransient) {
    vDeleteObjects_.push_back(new TH1F(histName, histTitle, nBinsX, lowX, highX));
    return vDeleteObjects_.back(); // return last element of vector
  }
  else
    return tfd.make<TH1F>(histName, histTitle, nBinsX, lowX, highX);


}


bool TrackerOfflineValidation::isBarrel(uint32_t subDetId)
{
  return (subDetId == StripSubdetector::TIB ||
	  subDetId == StripSubdetector::TOB ||
	  subDetId == PixelSubdetector::PixelBarrel );

}

bool TrackerOfflineValidation::isEndCap(uint32_t subDetId)
{
  return ( subDetId == StripSubdetector::TID ||
	   subDetId == StripSubdetector::TEC ||
	   subDetId == PixelSubdetector::PixelEndcap);
}

bool TrackerOfflineValidation::isPixel(uint32_t subDetId)
{
  return (subDetId == PixelSubdetector::PixelBarrel || subDetId == PixelSubdetector::PixelEndcap);
}


bool TrackerOfflineValidation::isDetOrDetUnit(align::StructureType type)
{
  return ( type == align::AlignableDet || type == align::AlignableDetUnit);
}

void 
TrackerOfflineValidation::getBinning(uint32_t subDetId, 
				     TrackerOfflineValidation::HistogrammType residualType, 
				     int &nBinsX, double &lowerBoundX, double &upperBoundX)
{
  // determine if 
  const bool isPixel = this->isPixel(subDetId);
  
  edm::ParameterSet binningPSet;
  
  switch(residualType) 
    {
    case XResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1XResPixelModules");                
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1XResStripModules");                
      break;
    case NormXResidual : 
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1NormXResPixelModules");             
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1NormXResStripModules");                
      break;
    case XprimeResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1XprimeResPixelModules");                
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1XprimeResStripModules");                
      break;
    case NormXprimeResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1NormXprimeResPixelModules");
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1NormXprimeResStripModules");
      break;
    case YprimeResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1YResPixelModules");                
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1YResStripModules");                
      break; 
    case NormYprimeResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1NormYResPixelModules");             
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1NormYResStripModules");  
      break;
    }
  nBinsX      = binningPSet.getParameter<int32_t>("Nbinx");		       
  lowerBoundX = binningPSet.getParameter<double>("xmin");		       
  upperBoundX = binningPSet.getParameter<double>("xmax");     
  
}

void 
TrackerOfflineValidation::setSummaryBin(int bin, TH1* targetHist, TH1* sourceHist)
{
  if(targetHist && sourceHist) {
    targetHist->SetBinContent(bin, sourceHist->GetMean(1));
    if(useFwhm_) targetHist->SetBinError(bin, Fwhm(sourceHist)/2.);
    else targetHist->SetBinError(bin, sourceHist->GetRMS(1) );
  } else {
    return;
  }

}


void 
TrackerOfflineValidation::summarizeBinInContainer( int bin, SummaryContainer &targetContainer, 
						   SummaryContainer &sourceContainer)
{
  
  
  this->setSummaryBin(bin, targetContainer.summaryXResiduals_, sourceContainer.sumXResiduals_);
  this->setSummaryBin(bin, targetContainer.summaryNormXResiduals_, sourceContainer.sumNormXResiduals_);
  // If no y-residual hists, just returns:
  this->setSummaryBin(bin, targetContainer.summaryYResiduals_, sourceContainer.sumYResiduals_);
  this->setSummaryBin(bin, targetContainer.summaryNormYResiduals_, sourceContainer.sumNormYResiduals_);

}

void 
TrackerOfflineValidation::summarizeBinInContainer( int bin, uint32_t subDetId, 
						   SummaryContainer &targetContainer, 
						   ModuleHistos &sourceContainer)
{

  // takes two summary Containers and sets summaryBins for all histogramms
  this->setSummaryBin(bin, targetContainer.summaryXResiduals_, sourceContainer.ResXprimeHisto);
  this->setSummaryBin(bin, targetContainer.summaryNormXResiduals_, sourceContainer.NormResXprimeHisto);
  if( this->isPixel(subDetId) || stripYResiduals_ ) {
    this->setSummaryBin(bin, targetContainer.summaryYResiduals_, sourceContainer.ResYprimeHisto);
    this->setSummaryBin(bin, targetContainer.summaryNormYResiduals_, sourceContainer.NormResYprimeHisto);
  }
}




TrackerOfflineValidation::ModuleHistos& 
TrackerOfflineValidation::getHistStructFromMap(const DetId& detid)
{

  // get a struct with histogramms from the respective map
  // if no object exist, the reference is automatically created by the map
  // throw exception if non-tracker id is passed
  uint subdetid = detid.subdetId();
  if(subdetid == PixelSubdetector::PixelBarrel ) {
    return mPxbResiduals_[detid.rawId()];
  } else if (subdetid == PixelSubdetector::PixelEndcap) {
    return mPxeResiduals_[detid.rawId()];
  } else if(subdetid  == StripSubdetector::TIB) {
    return mTibResiduals_[detid.rawId()];
  } else if(subdetid  == StripSubdetector::TID) {
    return mTidResiduals_[detid.rawId()];
  } else if(subdetid  == StripSubdetector::TOB) {
    return mTobResiduals_[detid.rawId()];
  } else if(subdetid  == StripSubdetector::TEC) {
    return mTecResiduals_[detid.rawId()];
  } else {
    throw cms::Exception("Geometry Error")
      << "[TrackerOfflineValidation] Error, tried to get reference for non-tracker subdet " << subdetid 
      << " from detector " << detid.det();
    return mPxbResiduals_[0];
  }
  
}


// ------------ method called to for each event  ------------
void
TrackerOfflineValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (useOverflowForRMS_)TH1::StatOverflows(kTRUE);
  this->checkBookHists(iSetup); // check whether hists are are booked and do so if not yet done
  
  //using namespace edm;
  TrackerValidationVariables avalidator_(iSetup,parset_);
  edm::Service<TFileService> fs;
    
  std::vector<TrackerValidationVariables::AVTrackStruct> vTrackstruct;
  avalidator_.fillTrackQuantities(iEvent,vTrackstruct);
  std::vector<TrackerValidationVariables::AVHitStruct> v_hitstruct;
  avalidator_.fillHitQuantities(iEvent,v_hitstruct);
  
  
  for (std::vector<TrackerValidationVariables::AVTrackStruct>::const_iterator it = vTrackstruct.begin(),
	 itEnd = vTrackstruct.end(); it != itEnd; ++it) {
    
    // Fill 1D track histos
    static const int etaindex = this->GetIndex(vTrackHistos_,"h_tracketa");
    vTrackHistos_[etaindex]->Fill(it->eta);
    static const int kappaindex = this->GetIndex(vTrackHistos_,"h_curvature");
    vTrackHistos_[kappaindex]->Fill(it->kappa);
    static const int kappaposindex = this->GetIndex(vTrackHistos_,"h_curvature_pos");
    if(it->charge > 0)
      vTrackHistos_[kappaposindex]->Fill(fabs(it->kappa));
    static const int kappanegindex = this->GetIndex(vTrackHistos_,"h_curvature_neg");
    if(it->charge < 0)
      vTrackHistos_[kappanegindex]->Fill(fabs(it->kappa));
    static const int normchi2index = this->GetIndex(vTrackHistos_,"h_normchi2");
    vTrackHistos_[normchi2index]->Fill(it->normchi2);
    static const int chi2index = this->GetIndex(vTrackHistos_,"h_chi2");
    vTrackHistos_[chi2index]->Fill(it->chi2);
    static const int ptindex = this->GetIndex(vTrackHistos_,"h_pt");
    vTrackHistos_[ptindex]->Fill(it->pt);
    if(it->ptError != 0.) {
      static const int ptResolutionindex = this->GetIndex(vTrackHistos_,"h_ptResolution");
      vTrackHistos_[ptResolutionindex]->Fill(it->ptError/it->pt);
    }
    // Fill track profiles
    static const int d0phiindex = this->GetIndex(vTrackProfiles_,"p_d0_vs_phi");
    vTrackProfiles_[d0phiindex]->Fill(it->phi,it->d0);
    static const int dzphiindex = this->GetIndex(vTrackProfiles_,"p_dz_vs_phi");
    vTrackProfiles_[dzphiindex]->Fill(it->phi,it->dz);
    static const int d0etaindex = this->GetIndex(vTrackProfiles_,"p_d0_vs_eta");
    vTrackProfiles_[d0etaindex]->Fill(it->eta,it->d0);
    static const int dzetaindex = this->GetIndex(vTrackProfiles_,"p_dz_vs_eta");
    vTrackProfiles_[dzetaindex]->Fill(it->eta,it->dz);
    static const int chiphiindex = this->GetIndex(vTrackProfiles_,"p_chi2_vs_phi");
    vTrackProfiles_[chiphiindex]->Fill(it->phi,it->chi2);
    static const int normchiphiindex = this->GetIndex(vTrackProfiles_,"p_normchi2_vs_phi");
    vTrackProfiles_[normchiphiindex]->Fill(it->phi,it->normchi2);
    static const int chietaindex = this->GetIndex(vTrackProfiles_,"p_chi2_vs_eta");
    vTrackProfiles_[chietaindex]->Fill(it->eta,it->chi2);
    static const int normchietaindex = this->GetIndex(vTrackProfiles_,"p_normchi2_vs_eta");
    vTrackProfiles_[normchietaindex]->Fill(it->eta,it->normchi2);
    static const int kappaphiindex = this->GetIndex(vTrackProfiles_,"p_kappa_vs_phi");
    vTrackProfiles_[kappaphiindex]->Fill(it->phi,it->kappa);
    static const int kappaetaindex = this->GetIndex(vTrackProfiles_,"p_kappa_vs_eta");
    vTrackProfiles_[kappaetaindex]->Fill(it->eta,it->kappa);
    static const int ptResphiindex = this->GetIndex(vTrackProfiles_,"p_ptResolution_vs_phi");
    vTrackProfiles_[ptResphiindex]->Fill(it->phi,it->ptError/it->pt);
    static const int ptResetaindex = this->GetIndex(vTrackProfiles_,"p_ptResolution_vs_eta");
    vTrackProfiles_[ptResetaindex]->Fill(it->eta,it->ptError/it->pt);

    // Fill 2D track histos
    static const int d0phiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_d0_vs_phi");
    vTrack2DHistos_[d0phiindex_2d]->Fill(it->phi,it->d0);
    static const int dzphiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_dz_vs_phi");
    vTrack2DHistos_[dzphiindex_2d]->Fill(it->phi,it->dz);
    static const int d0etaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_d0_vs_eta");
    vTrack2DHistos_[d0etaindex_2d]->Fill(it->eta,it->d0);
    static const int dzetaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_dz_vs_eta");
    vTrack2DHistos_[dzetaindex_2d]->Fill(it->eta,it->dz);
    static const int chiphiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_chi2_vs_phi");
    vTrack2DHistos_[chiphiindex_2d]->Fill(it->phi,it->chi2);
    static const int normchiphiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_normchi2_vs_phi");
    vTrack2DHistos_[normchiphiindex_2d]->Fill(it->phi,it->normchi2);
    static const int chietaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_chi2_vs_eta");
    vTrack2DHistos_[chietaindex_2d]->Fill(it->eta,it->chi2);
    static const int normchietaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_normchi2_vs_eta");
    vTrack2DHistos_[normchietaindex_2d]->Fill(it->eta,it->normchi2);
    static const int kappaphiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_kappa_vs_phi");
    vTrack2DHistos_[kappaphiindex_2d]->Fill(it->phi,it->kappa);
    static const int kappaetaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_kappa_vs_eta");
    vTrack2DHistos_[kappaetaindex_2d]->Fill(it->eta,it->kappa);
     
  } // finish loop over track quantities


  // hit quantities: residuals, normalized residuals
  for (std::vector<TrackerValidationVariables::AVHitStruct>::const_iterator it = v_hitstruct.begin(),
  	 itEnd = v_hitstruct.end(); it != itEnd; ++it) {
    DetId detid(it->rawDetId);
    ModuleHistos &histStruct = this->getHistStructFromMap(detid);
    
    // fill histos in local coordinates if set in cf
    if(lCoorHistOn_) {
      histStruct.ResHisto->Fill(it->resX);
      if(it->resErrX != 0) histStruct.NormResHisto->Fill(it->resX/it->resErrX);
    }
    if(it->resXprime != -999.) {
      histStruct.ResXprimeHisto->Fill(it->resXprime);
      if(it->resXprimeErr != 0 && it->resXprimeErr != -999 ) {	
	histStruct.NormResXprimeHisto->Fill(it->resXprime/it->resXprimeErr);
      } 
    }
    if(it->resYprime != -999.) {
      if( this->isPixel(detid.subdetId())  || stripYResiduals_ ) {
	histStruct.ResYprimeHisto->Fill(it->resYprime);
	if(it->resYprimeErr != 0 && it->resYprimeErr != -999. ) {	
	  histStruct.NormResYprimeHisto->Fill(it->resYprime/it->resYprimeErr);
	} 
      }
    }

    
    if(overlappOn_) {
      std::pair<uint32_t,uint32_t> tmp_pair(std::make_pair(it->rawDetId, it->overlapres.first));
      if(it->overlapres.first != 0 ) {
	if( hOverlappResidual[tmp_pair] ) {
	  hOverlappResidual[tmp_pair]->Fill(it->overlapres.second);
	} else if( hOverlappResidual[std::make_pair( it->overlapres.first, it->rawDetId) ]) {
	  hOverlappResidual[std::make_pair( it->overlapres.first, it->rawDetId) ]->Fill(it->overlapres.second);
	} else {
	  TFileDirectory tfd = fs->mkdir("OverlappResiduals");
	  hOverlappResidual[tmp_pair] = tfd.make<TH1F>(Form("hOverlappResidual_%d_%d",tmp_pair.first,tmp_pair.second),
						       "Overlapp Residuals",100,-50,50);
	  hOverlappResidual[tmp_pair]->Fill(it->overlapres.second);
	}
      }
    } // end overlappOn

  }
  if (useOverflowForRMS_) TH1::StatOverflows(kFALSE);  
}



// ------------ method called once each job just after ending the event loop  ------------
void 
TrackerOfflineValidation::endJob()
{
  AlignableTracker aliTracker(&(*tkGeom_));
  edm::Service<TFileService> fs;   
  AlignableObjectId aliobjid;

  TTree *tree = fs->make<TTree>("TkOffVal","TkOffVal");
  TkOffTreeVariables *treeMemPtr = new TkOffTreeVariables;
  // We create branches for all members of 'TkOffTreeVariables' (even if not needed).
  // This works because we have a dictionary for 'TkOffTreeVariables'
  // (see src/classes_def.xml and src/classes.h):
  tree->Branch("TkOffTreeVariables", &treeMemPtr); // address of pointer!

  this->fillTree(*tree, mPxbResiduals_, *treeMemPtr, *tkGeom_);
  this->fillTree(*tree, mPxeResiduals_, *treeMemPtr, *tkGeom_);
  this->fillTree(*tree, mTibResiduals_, *treeMemPtr, *tkGeom_);
  this->fillTree(*tree, mTidResiduals_, *treeMemPtr, *tkGeom_);
  this->fillTree(*tree, mTobResiduals_, *treeMemPtr, *tkGeom_);
  this->fillTree(*tree, mTecResiduals_, *treeMemPtr, *tkGeom_);

  delete treeMemPtr; treeMemPtr = 0;

  static const int kappadiffindex = this->GetIndex(vTrackHistos_,"h_diff_curvature");
  vTrackHistos_[kappadiffindex]->Add(vTrackHistos_[this->GetIndex(vTrackHistos_,"h_curvature_neg")],
				     vTrackHistos_[this->GetIndex(vTrackHistos_,"h_curvature_pos")],-1,1);

  // Collate Information for Subdetectors
  // create summary histogramms recursively
 
  std::vector<TrackerOfflineValidation::SummaryContainer > vTrackerprofiles;
  this->collateSummaryHists((*fs),(aliTracker), 0, aliobjid, vTrackerprofiles);
   
}


void
TrackerOfflineValidation::collateSummaryHists( TFileDirectory &tfd, const Alignable& ali, int i, 
					       const AlignableObjectId &aliobjid, 
					       std::vector< TrackerOfflineValidation::SummaryContainer > &v_levelProfiles)
{
  
  std::vector<Alignable*> alivec(ali.components());
  if( this->isDetOrDetUnit((alivec)[0]->alignableObjectId()) ) return;

  for(int iComp=0, iCompEnd = ali.components().size();iComp < iCompEnd; ++iComp) {
    std::vector< TrackerOfflineValidation::SummaryContainer > v_profiles;        
    std::string structurename  = aliobjid.typeToName((alivec)[iComp]->alignableObjectId());
 
    LogDebug("TrackerOfflineValidation") << "StructureName = " << structurename;
    std::stringstream dirname;
    dirname << structurename;
    
    // add no suffix counter to strip and pixel -> just aesthetics
    if (structurename != "Strip" && structurename != "Pixel") dirname << "_" << iComp+1;
    
    if(  !(this->isDetOrDetUnit( (alivec)[iComp]->alignableObjectId()) )
	 || (alivec)[0]->components().size() > 1 ) {
      TFileDirectory f = tfd.mkdir((dirname.str()).c_str());
      this->collateSummaryHists( f, *(alivec)[iComp], i, aliobjid, v_profiles);
      v_levelProfiles.push_back(this->bookSummaryHists(tfd, *(alivec[iComp]), ali.alignableObjectId(), iComp, aliobjid));
      TH1 *hY = v_levelProfiles[iComp].sumYResiduals_;
      TH1 *hNormY = v_levelProfiles[iComp].sumNormYResiduals_;
      for(uint n = 0; n < v_profiles.size(); ++n) {
	this->summarizeBinInContainer(n+1, v_levelProfiles[iComp], v_profiles[n] );
	v_levelProfiles[iComp].sumXResiduals_->Add(v_profiles[n].sumXResiduals_);
	v_levelProfiles[iComp].sumNormXResiduals_->Add(v_profiles[n].sumNormXResiduals_);
	if (hY)     hY->Add(v_profiles[n].sumYResiduals_);         // only if existing
	if (hNormY) hNormY->Add(v_profiles[n].sumNormYResiduals_); // dito (pxl, stripYResiduals_)
      }
      //add fit values to stat box
      this->fitResiduals(v_levelProfiles[iComp].sumXResiduals_);
      this->fitResiduals(v_levelProfiles[iComp].sumNormXResiduals_);
      if (hY)     this->fitResiduals(hY);     // only if existing (pixel or stripYResiduals_)
      if (hNormY) this->fitResiduals(hNormY); // dito
    } else {
      // nothing to be done for det or detunits
      continue;
    }

  }

}

TrackerOfflineValidation::SummaryContainer 
TrackerOfflineValidation::bookSummaryHists(TFileDirectory &tfd, const Alignable& ali, 
					   align::StructureType type, int i, 
					   const AlignableObjectId &aliobjid)
{
  const uint aliSize = ali.components().size();
  const align::StructureType alitype = ali.alignableObjectId();
  const align::StructureType subtype = ali.components()[0]->alignableObjectId();
  const char *aliTypeName = aliobjid.typeToName(alitype).c_str(); // lifetime of char* OK
  const char *aliSubtypeName = aliobjid.typeToName(subtype).c_str();
  const char *typeName = aliobjid.typeToName(type).c_str();

  const DetId aliDetId = ali.id(); 
  // y residuals only if pixel or specially requested for strip:
  const bool bookResidY = this->isPixel(aliDetId.subdetId()) || stripYResiduals_;

  SummaryContainer sumContainer;
  
  // Book summary hists with one bin per component, 
  // but special case for Det with two DetUnit that we want to summarize one level up 
  // (e.g. in TOBRods with 12 bins for 6 stereo and 6 rphi DetUnit.)
  //    component of ali is not Det or Det with just one components
  const uint subcompSize = ali.components()[0]->components().size();
  if (subtype != align::AlignableDet || subcompSize == 1) { // Det with 1 comp. should not exist anymore...
    const TString title(Form("Summary for substructures in %s %d;%s;",aliTypeName,i,aliSubtypeName));
    sumContainer.summaryXResiduals_ = tfd.make<TH1F>(Form("h_summaryX%s_%d",aliTypeName,i), 
						     title + "#LT #Delta x' #GT",
						     aliSize, 0.5, aliSize+0.5);
    sumContainer.summaryNormXResiduals_ = tfd.make<TH1F>(Form("h_summaryNormX%s_%d",aliTypeName,i), 
							 title + "#LT #Delta x'/#sigma #GT",
							 aliSize,0.5,aliSize+0.5);
    
    if (bookResidY) {
      sumContainer.summaryYResiduals_ = tfd.make<TH1F>(Form("h_summaryY%s_%d",aliTypeName,i), 
						       title + "#LT #Delta y' #GT",
						       aliSize, 0.5, aliSize+0.5);
      sumContainer.summaryNormYResiduals_ = tfd.make<TH1F>(Form("h_summaryNormY%s_%d",aliTypeName,i), 
							   title + "#LT #Delta y'/#sigma #GT",
							   aliSize,0.5,aliSize+0.5);
    }
    
  } else if (subtype == align::AlignableDet && subcompSize > 1) { // fixed: was aliSize before
    if (subcompSize != 2) { // strange... expect only 2 DetUnits in DS layers
      // this 2 is hardcoded factor 2 in binning below and also assummed later on
      edm::LogError("Alignment") << "@SUB=bookSummaryHists"
				 << "Det with " << subcompSize << " components";
    }
    // title contains x-title
    const TString title(Form("Summary for substructures in %s %d;%s;", aliTypeName, i,
			     aliobjid.typeToName(ali.components()[0]->components()[0]->alignableObjectId()).c_str()));
    sumContainer.summaryXResiduals_ 
      = tfd.make<TH1F>(Form("h_summaryX%s_%d", aliTypeName, i), 
		       title + "#LT #Delta x' #GT", (2*aliSize), 0.5, 2*aliSize+0.5);
    sumContainer.summaryNormXResiduals_ 
      = tfd.make<TH1F>(Form("h_summaryNormX%s_%d", aliTypeName, i), 
		       title + "#LT #Delta x'/#sigma #GT", (2*aliSize), 0.5, 2*aliSize+0.5);

    if (bookResidY) {
      sumContainer.summaryYResiduals_ 
	= tfd.make<TH1F>(Form("h_summaryY%s_%d", aliTypeName, i), 
			 title + "#LT #Delta y' #GT", (2*aliSize), 0.5, 2*aliSize+0.5);
      sumContainer.summaryNormYResiduals_ 
	= tfd.make<TH1F>(Form("h_summaryNormY%s_%d", aliTypeName, i), 
			 title + "#LT #Delta y'/#sigma #GT", (2*aliSize), 0.5, 2*aliSize+0.5);
    }

  } else {
    edm::LogError("TrackerOfflineValidation") << "@SUB=TrackerOfflineValidation::bookSummaryHists" 
					      << "No summary histogramm for hierarchy level " 
					      << aliTypeName << " in subdet " << aliDetId.subdetId();
  }

  // Now book hists that just sum up the residual histograms from lower levels.
  // Axis title is copied from lowest level module of structure.
  // Should be safe that y-hists are only touched if non-null pointers...
  int nbins = 0;
  double xmin = 0., xmax = 0.;
  const TString sumTitle(Form("Residual for %s %d in %s;", aliTypeName, i, typeName));
  const ModuleHistos &xTitHists = this->getHistStructFromMap(aliDetId); // for x-axis titles
  this->getBinning(aliDetId.subdetId(), XprimeResidual, nbins, xmin, xmax);
  sumContainer.sumXResiduals_ = tfd.make<TH1F>(Form("h_Xprime_%s_%d", aliTypeName, i),
					       sumTitle + xTitHists.ResXprimeHisto->GetXaxis()->GetTitle(),
					       nbins, xmin, xmax);
  
  this->getBinning(aliDetId.subdetId(), NormXprimeResidual, nbins, xmin, xmax);
  sumContainer.sumNormXResiduals_ = tfd.make<TH1F>(Form("h_NormXprime_%s_%d",aliTypeName,i), 
						   sumTitle + xTitHists.NormResXprimeHisto->GetXaxis()->GetTitle(),
						   nbins, xmin, xmax);
  if (bookResidY) {
    this->getBinning(aliDetId.subdetId(), YprimeResidual, nbins, xmin, xmax);
    sumContainer.sumYResiduals_ = tfd.make<TH1F>(Form("h_Yprime_%s_%d",aliTypeName,i), 
						 sumTitle + xTitHists.ResYprimeHisto->GetXaxis()->GetTitle(),
						 nbins, xmin, xmax);
    
    this->getBinning(aliDetId.subdetId(), NormYprimeResidual, nbins, xmin, xmax);
    sumContainer.sumNormYResiduals_ = tfd.make<TH1F>(Form("h_NormYprime_%s_%d",aliTypeName,i), 
						     sumTitle + xTitHists.NormResYprimeHisto->GetXaxis()->GetTitle(),
						     nbins, xmin, xmax);
  }
  
  // If we are at the lowest level, we already sum up and fill the summary.

  // special case I: For DetUnits and Detwith  only one subcomponent start filling summary histos
  if( (  subtype == align::AlignableDet && subcompSize == 1) || 
      subtype  == align::AlignableDetUnit  
      ) {
    for(uint k = 0; k < aliSize; ++k) {
      DetId detid = ali.components()[k]->id();
      ModuleHistos &histStruct = this->getHistStructFromMap(detid);
      this->summarizeBinInContainer(k+1, detid.subdetId() ,sumContainer, histStruct );
      sumContainer.sumXResiduals_->Add(histStruct.ResXprimeHisto);
      sumContainer.sumNormXResiduals_->Add(histStruct.NormResXprimeHisto);
      if( this->isPixel(detid.subdetId()) || stripYResiduals_ ) {
      	sumContainer.sumYResiduals_->Add(histStruct.ResYprimeHisto);
      	sumContainer.sumNormYResiduals_->Add(histStruct.NormResYprimeHisto);
      }
    }
  } else if( subtype == align::AlignableDet && subcompSize > 1) { // fixed: was aliSize before
    // special case II: Fill summary histos for dets with two detunits 
    for(uint k = 0; k < aliSize; ++k) {
      for(uint j = 0; j < subcompSize; ++j) { // assumes all have same size (as binning does)
	DetId detid = ali.components()[k]->components()[j]->id();
	ModuleHistos &histStruct = this->getHistStructFromMap(detid);	
	this->summarizeBinInContainer(2*k+j+1, detid.subdetId() ,sumContainer, histStruct );
	sumContainer.sumXResiduals_->Add( histStruct.ResXprimeHisto);
	sumContainer.sumNormXResiduals_->Add( histStruct.NormResXprimeHisto);
	if( this->isPixel(detid.subdetId()) || stripYResiduals_ ) {
	  sumContainer.sumYResiduals_->Add( histStruct.ResYprimeHisto);
	  sumContainer.sumNormYResiduals_->Add( histStruct.NormResYprimeHisto);
	}
      }
    }
  }
  
  return sumContainer;
}


float 
TrackerOfflineValidation::Fwhm (const TH1* hist) const
{
  float max = hist->GetMaximum();
  int left = -1, right = -1;
  for(unsigned int i = 1, iEnd = hist->GetNbinsX(); i <= iEnd; ++i) {
    if(hist->GetBinContent(i) < max/2. && hist->GetBinContent(i+1) > max/2. && left == -1) {
      if(max/2. - hist->GetBinContent(i) < hist->GetBinContent(i+1) - max/2.) {
	left = i;
	++i;
      } else {
	left = i+1;
	++i;
      }
    }
    if(left != -1 && right == -1) {
      if(hist->GetBinContent(i) > max/2. && hist->GetBinContent(i+1) < max/2.) {
	if( hist->GetBinContent(i) - max/2. < max/2. - hist->GetBinContent(i+1)) {
	  right = i;
	} else {
	  right = i+1;
	}
	
      }
    }
  }
  return hist->GetXaxis()->GetBinCenter(right) - hist->GetXaxis()->GetBinCenter(left);
}

////////////////////////////////////////////////////////////////////////////////////
void 
TrackerOfflineValidation::fillTree(TTree &tree,
				   const std::map<int, TrackerOfflineValidation::ModuleHistos> &moduleHist_,
				   TkOffTreeVariables &treeMem, const TrackerGeometry &tkgeom)
{
 
  for(std::map<int, TrackerOfflineValidation::ModuleHistos>::const_iterator it = moduleHist_.begin(), 
	itEnd= moduleHist_.end(); it != itEnd;++it ) { 
    treeMem.clear(); // make empty/default
    //variables concerning the tracker components/hierarchy levels
    DetId detId_ = it->first;
    treeMem.moduleId = detId_;
    treeMem.subDetId = detId_.subdetId();
    treeMem.isDoubleSide =0;

    if(treeMem.subDetId == PixelSubdetector::PixelBarrel){
      PXBDetId pxbId(detId_); 
      treeMem.layer = pxbId.layer(); 
      treeMem.rod = pxbId.ladder();
  
    } else if(treeMem.subDetId == PixelSubdetector::PixelEndcap){
      PXFDetId pxfId(detId_); 
      treeMem.layer = pxfId.disk(); 
      treeMem.side = pxfId.side();
      treeMem.blade = pxfId.blade(); 
      treeMem.panel = pxfId.panel();

    } else if(treeMem.subDetId == StripSubdetector::TIB){
      TIBDetId tibId(detId_); 
      treeMem.layer = tibId.layer(); 
      treeMem.side = tibId.string()[0];
      treeMem.rod = tibId.string()[2]; 
      treeMem.outerInner = tibId.string()[1]; 
      treeMem.isStereo = tibId.stereo();
      treeMem.isDoubleSide = tibId.isDoubleSide();
    } else if(treeMem.subDetId == StripSubdetector::TID){
      TIDDetId tidId(detId_); 
      treeMem.layer = tidId.wheel(); 
      treeMem.side = tidId.side();
      treeMem.ring = tidId.ring(); 
      treeMem.outerInner = tidId.module()[0]; 
      treeMem.isStereo = tidId.stereo();
      treeMem.isDoubleSide = tidId.isDoubleSide();
    } else if(treeMem.subDetId == StripSubdetector::TOB){
      TOBDetId tobId(detId_); 
      treeMem.layer = tobId.layer(); 
      treeMem.side = tobId.rod()[0];
      treeMem.rod = tobId.rod()[1]; 
      treeMem.isStereo = tobId.stereo();
      treeMem.isDoubleSide = tobId.isDoubleSide();
    } else if(treeMem.subDetId == StripSubdetector::TEC) {
      TECDetId tecId(detId_); 
      treeMem.layer = tecId.wheel(); 
      treeMem.side  = tecId.side();
      treeMem.ring  = tecId.ring(); 
      treeMem.petal = tecId.petal()[1]; 
      treeMem.outerInner = tecId.petal()[0];
      treeMem.isStereo = tecId.stereo();
      treeMem.isDoubleSide = tecId.isDoubleSide(); 
    }
    
    //variables concerning the tracker geometry
    
    const Surface::PositionType &gPModule = tkgeom.idToDet(detId_)->position();
    treeMem.posPhi = gPModule.phi();
    treeMem.posEta = gPModule.eta();
    treeMem.posR   = gPModule.perp();
    treeMem.posX   = gPModule.x();
    treeMem.posY   = gPModule.y();
    treeMem.posZ   = gPModule.z();
 
    const Surface& surface =  tkgeom.idToDet(detId_)->surface();
    
    //global Orientation of local coordinate system of dets/detUnits   
    LocalPoint  lUDirection(1.,0.,0.), lVDirection(0.,1.,0.), lWDirection(0.,0.,1.);
    GlobalPoint gUDirection = surface.toGlobal(lUDirection),
                gVDirection = surface.toGlobal(lVDirection),
		gWDirection = surface.toGlobal(lWDirection);
    double dR(999.), dPhi(999.), dZ(999.);
    if(treeMem.subDetId==PixelSubdetector::PixelBarrel || treeMem.subDetId==StripSubdetector::TIB || treeMem.subDetId==StripSubdetector::TOB){
      dR = gWDirection.perp() - gPModule.perp();
      dPhi = deltaPhi(gUDirection.phi(),gPModule.phi());
      dZ = gVDirection.z() - gPModule.z();
      if(dZ>=0.)treeMem.rOrZDirection = 1; else treeMem.rOrZDirection = -1;
    }else if(treeMem.subDetId==PixelSubdetector::PixelEndcap){
      dR = gUDirection.perp() - gPModule.perp();
      dPhi = deltaPhi(gVDirection.phi(),gPModule.phi());
      dZ = gWDirection.z() - gPModule.z();
      if(dR>=0.)treeMem.rOrZDirection = 1; else treeMem.rOrZDirection = -1;
    }else if(treeMem.subDetId==StripSubdetector::TID || treeMem.subDetId==StripSubdetector::TEC){
      dR = gVDirection.perp() - gPModule.perp();
      dPhi = deltaPhi(gUDirection.phi(),gPModule.phi());
      dZ = gWDirection.z() - gPModule.z();
      if(dR>=0.)treeMem.rOrZDirection = 1; else treeMem.rOrZDirection = -1;
    }
    if(dR>=0.)treeMem.rDirection = 1; else treeMem.rDirection = -1;
    if(dPhi>=0.)treeMem.phiDirection = 1; else treeMem.phiDirection = -1;
    if(dZ>=0.)treeMem.zDirection = 1; else treeMem.zDirection = -1;
    
    
    //mean and RMS values (extracted from histograms(Xprime on module level)
    treeMem.entries = static_cast<UInt_t>(it->second.ResXprimeHisto->GetEntries());
    treeMem.meanX = it->second.ResXprimeHisto->GetMean();
    treeMem.rmsX  = it->second.ResXprimeHisto->GetRMS();
    //treeMem.sigmaX = Fwhm(it->second.ResXprimeHisto)/2.355;
    if (useFit_) {
      
      //call fit function which returns mean and sigma from the fit
      //for absolute residuals
      std::pair<float,float> fitResult1 = this->fitResiduals(it->second.ResXprimeHisto);
      treeMem.fitMeanX = fitResult1.first;
      treeMem.fitSigmaX = fitResult1.second;
      //for normalized residuals
      std::pair<float,float> fitResult2 = this->fitResiduals(it->second.NormResXprimeHisto);
      treeMem.fitMeanNormX = fitResult2.first;
      treeMem.fitSigmaNormX = fitResult2.second;
    }
    //get median for absolute residuals
    treeMem.medianX   = this->getMedian(it->second.ResXprimeHisto);


    int numberOfBins=it->second.ResXprimeHisto->GetNbinsX();
    treeMem.numberOfUnderflows = it->second.ResXprimeHisto->GetBinContent(0);
    treeMem.numberOfOverflows = it->second.ResXprimeHisto->GetBinContent(numberOfBins+1);
    treeMem.numberOfOutliers =  it->second.ResXprimeHisto->GetBinContent(0)+it->second.ResXprimeHisto->GetBinContent(numberOfBins+1);
    //mean and RMS values (extracted from histograms(normalized Xprime on module level)
    treeMem.meanNormX = it->second.NormResXprimeHisto->GetMean();
    treeMem.rmsNormX = it->second.NormResXprimeHisto->GetRMS();

    double stats[20];
    it->second.NormResXprimeHisto->GetStats(stats);
    // GF  treeMem.chi2PerDofX = stats[3]/(stats[0]-1);
    if (stats[0]) treeMem.chi2PerDofX = stats[3]/stats[0];
    
    treeMem.sigmaNormX = Fwhm(it->second.NormResXprimeHisto)/2.355;
    treeMem.histNameX = it->second.ResXprimeHisto->GetName();
    treeMem.histNameNormX = it->second.NormResXprimeHisto->GetName();
    

    // fill tree variables in local coordinates if set in cfg
    if(lCoorHistOn_) {
      treeMem.meanLocalX = it->second.ResHisto->GetMean();
      treeMem.rmsLocalX = it->second.ResHisto->GetRMS();
      treeMem.meanNormLocalX = it->second.NormResHisto->GetMean();
      treeMem.rmsNormLocalX = it->second.NormResHisto->GetRMS();

      treeMem.histNameLocalX = it->second.ResHisto->GetName();
      treeMem.histNameNormLocalX = it->second.NormResHisto->GetName();
    }

    // mean and RMS values in local y (extracted from histograms(normalized Yprime on module level)
    // might exist in pixel only
    if (it->second.ResYprimeHisto) {//(stripYResiduals_){
      TH1 *h = it->second.ResYprimeHisto;
      treeMem.meanY = h->GetMean();
      treeMem.rmsY  = h->GetRMS();
      
      if (useFit_) { // fit function which returns mean and sigma from the fit
	std::pair<float,float> fitMeanSigma = this->fitResiduals(h);
	treeMem.fitMeanY  = fitMeanSigma.first;
	treeMem.fitSigmaY = fitMeanSigma.second;
      }
      //get median for absolute residuals
      treeMem.medianY   = this->getMedian(h);

      treeMem.histNameY = h->GetName();
    }
    if (it->second.NormResYprimeHisto) {
      TH1 *h = it->second.NormResYprimeHisto;
      treeMem.meanNormY = h->GetMean();
      treeMem.rmsNormY  = h->GetRMS();
      h->GetStats(stats); // stats buffer defined above
      if (stats[0]) treeMem.chi2PerDofY = stats[3]/stats[0];

      if (useFit_) { // fit function which returns mean and sigma from the fit
	std::pair<float,float> fitMeanSigma = this->fitResiduals(h);
	treeMem.fitMeanNormY  = fitMeanSigma.first;
	treeMem.fitSigmaNormY = fitMeanSigma.second;
      }
      treeMem.histNameNormY = h->GetName();
    }
    
    tree.Fill();
  }
}

std::pair<float,float> 
TrackerOfflineValidation::fitResiduals(TH1 *hist) const
{
  std::pair<float,float> fitResult(9999., 9999.);
  if (!hist || hist->GetEntries() < 20) return fitResult;

  float mean  = hist->GetMean();
  float sigma = hist->GetRMS();

  try { // for < CMSSW_2_2_0 since ROOT warnings from fit are converted to exceptions
    // Remove the try/catch for more recent CMSSW!
    // first fit: two RMS around mean
    TF1 func("tmp", "gaus", mean - 2.*sigma, mean + 2.*sigma); 
    if (0 == hist->Fit(&func,"QNR")) { // N: do not blow up file by storing fit!
      mean  = func.GetParameter(1);
      sigma = func.GetParameter(2);
      // second fit: three sigma of first fit around mean of first fit
      func.SetRange(mean - 3.*sigma, mean + 3.*sigma);
      // I: integral gives more correct results if binning is too wide
      // L: Likelihood can treat empty bins correctly (if hist not weighted...)
      if (0 == hist->Fit(&func, "Q0LR")) {
	if (hist->GetFunction(func.GetName())) { // Take care that it is later on drawn:
	  hist->GetFunction(func.GetName())->ResetBit(TF1::kNotDraw);
	}
	fitResult.first = func.GetParameter(1);
	fitResult.second = func.GetParameter(2);
      }
    }
  } catch (cms::Exception const & e) {
    edm::LogWarning("Alignment") << "@SUB=TrackerOfflineValidation::fitResiduals"
				 << "Caught this exception during ROOT fit: "
				 << e.what();
  }
  
  return fitResult;
}
float 
TrackerOfflineValidation::getMedian(const TH1 *histo) const
{

  float median = 999;
  int nbins = histo->GetNbinsX();

 
  //extract median from histogram
  double *x = new double[nbins];
  double *y = new double[nbins];
  for (int j = 0; j < nbins; j++) {
    x[j] = histo->GetBinCenter(j+1);
    y[j] = histo->GetBinContent(j+1);
  }
  median = TMath::Median(nbins, x, y);
  

  delete[] x; x = 0;
  delete [] y; y = 0;  

  return median;

}
//define this as a plug-in
DEFINE_FWK_MODULE(TrackerOfflineValidation);
