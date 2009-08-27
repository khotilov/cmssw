////////////////////////////////////////////////////////////////////////////////
// Package:          TrackingAnalysis/Cosmics
// Class:            HitEff
// Original Author:  Daniele Benedetti-INFN perugia
//   Update for P5:  Keith Ulmer--University of Colorado
//
///////////////////////////////////////////////////////////////////////////////

// system include files
#include <memory>
#include <string>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingAnalysis/Cosmics/interface/HitEff.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "TrackingAnalysis/Cosmics/interface/TrajectoryAtValidHit.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"

#include "AnalysisDataFormats/SiStripClusterInfo/interface/SiStripClusterInfo.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h" 

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


#include "TMath.h"
#include "TH1F.h"

//
// constructors and destructor
//

using namespace std;
HitEff::HitEff(const edm::ParameterSet& conf) : 
  conf_(conf)
{
  layers =conf_.getParameter<int>("Layer");
  DEBUG = conf_.getParameter<bool>("Debug");
}

// Virtual destructor needed.
HitEff::~HitEff() { }

void HitEff::beginJob(const edm::EventSetup& c){

  edm::Service<TFileService> fs;

  traj = fs->make<TTree>("traj","tree of trajectory positions");
  traj->Branch("TrajGlbX",&TrajGlbX,"TrajGlbX/F");
  traj->Branch("TrajGlbY",&TrajGlbY,"TrajGlbY/F");
  traj->Branch("TrajGlbZ",&TrajGlbZ,"TrajGlbZ/F");
  traj->Branch("TrajLocX",&TrajLocX,"TrajLocX/F");
  traj->Branch("TrajLocY",&TrajLocY,"TrajLocY/F");
  traj->Branch("TrajLocErrX",&TrajLocErrX,"TrajLocErrX/F");
  traj->Branch("TrajLocErrY",&TrajLocErrY,"TrajLocErrY/F");
  traj->Branch("TrajLocAngleX",&TrajLocAngleX,"TrajLocAngleX/F");
  traj->Branch("TrajLocAngleY",&TrajLocAngleY,"TrajLocAngleY/F");
  traj->Branch("ClusterLocX",&ClusterLocX,"ClusterLocX/F");
  traj->Branch("ClusterLocY",&ClusterLocY,"ClusterLocY/F");
  traj->Branch("ClusterLocErrX",&ClusterLocErrX,"ClusterLocErrX/F");
  traj->Branch("ClusterLocErrY",&ClusterLocErrY,"ClusterLocErrY/F");
  traj->Branch("ClusterStoN",&ClusterStoN,"ClusterStoN/F");
  traj->Branch("ResX",&ResX,"ResX/F");
  traj->Branch("ResXSig",&ResXSig,"ResXSig/F");
  traj->Branch("ModIsBad",&ModIsBad,"ModIsBad/i");
  traj->Branch("SiStripQualBad",&SiStripQualBad,"SiStripQualBad/i");
  traj->Branch("withinAcceptance",&withinAcceptance, "withinAcceptance/O");
  traj->Branch("Id",&Id,"Id/i");
  traj->Branch("run",&run,"run/i");
  traj->Branch("event",&event,"event/i");
  traj->Branch("layer",&whatlayer,"layer/i");
  traj->Branch("timeInOut",&timeInOut,"timeInOut/F");
  traj->Branch("timeInOutErr",&timeInOutErr,"timeInOutErr/F");
  traj->Branch("timeOutIn",&timeOutIn,"timeOutIn/F");
  traj->Branch("timeOutInErr",&timeOutInErr,"timeOutInErr/F");
  traj->Branch("timeDOF",&timeDOF,"timeDOF/I");

  events = 0;
  EventTrackCKF = 0;
  
}


void HitEff::analyze(const edm::Event& e, const edm::EventSetup& es){

  //  bool DEBUG = false;

  if (DEBUG)  cout << "beginning analyze from HitEff" << endl;

  using namespace edm;
  using namespace reco;
  // Step A: Get Inputs 

  int run_nr = e.id().run();
  int ev_nr = e.id().event();

  //CombinatoriaTrack
  edm::Handle<reco::TrackCollection> trackCollectionCKF;
  edm::InputTag TkTagCKF = conf_.getParameter<edm::InputTag>("combinatorialTracks");
  e.getByLabel(TkTagCKF,trackCollectionCKF);
  
  edm::Handle<std::vector<Trajectory> > TrajectoryCollectionCKF;
  edm::InputTag TkTrajCKF = conf_.getParameter<edm::InputTag>("trajectories");
  e.getByLabel(TkTrajCKF,TrajectoryCollectionCKF);
  
  // Clusters
  // get the SiStripClusters from the event
  edm::Handle< edmNew::DetSetVector<SiStripCluster> > theClusters;
  e.getByLabel("siStripClusters", theClusters);

  //get tracker geometry
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  const TrackerGeometry * tkgeom=&(* tracker);

  //get Cluster Parameter Estimator
  //std::string cpe = conf_.getParameter<std::string>("StripCPE");
  edm::ESHandle<StripClusterParameterEstimator> parameterestimator;
  es.get<TkStripCPERecord>().get("StripCPEfromTrackAngle", parameterestimator); 
  const StripClusterParameterEstimator &stripcpe(*parameterestimator);

  // get the SiStripQuality records
  edm::ESHandle<SiStripQuality> SiStripQuality_;
  try {
    es.get<SiStripQualityRcd>().get("forCluster",SiStripQuality_);
  }
  catch (...) {
    es.get<SiStripQualityRcd>().get(SiStripQuality_);
  }
  
  edm::ESHandle<MagneticField> magFieldHandle;
  es.get<IdealMagneticFieldRecord>().get(magFieldHandle);
  const MagneticField* magField_ = magFieldHandle.product();

  timeInOut = -999.0; timeInOutErr = -999.0; timeOutIn = -999.0; timeOutInErr = -999.0; timeDOF = -999;
  // retrieve the muon time
  Handle<MuonCollection> muH;
  e.getByLabel("muonsWitht0Correction",muH);
  const MuonCollection & muonsT0  =  *muH.product();
  cout << "muonsT0 size = " << muonsT0.size() << endl;
  if(muonsT0.size()<1) return;
  MuonTime mt0 = muonsT0[0].time();
  timeInOut = mt0.timeAtIpInOut; 
  timeInOutErr = mt0.timeAtIpInOutErr;
  timeOutIn = mt0.timeAtIpOutIn; 
  timeOutInErr = mt0.timeAtIpOutInErr;
  timeDOF = mt0.nDof;
    cout << "timeInOut = " << timeInOut << "+-" << timeInOutErr << "    timeOutIn = " << timeOutIn << "+-" << timeOutInErr << " with dof = " << timeDOF << endl;

  events++;
  
  // *************** SiStripCluster Collection
  const edmNew::DetSetVector<SiStripCluster>& input = *theClusters;

  //go through clusters to write out global position of good clusters for the layer understudy for comparison
  // Loop through clusters just to print out locations

  for (edmNew::DetSetVector<SiStripCluster>::const_iterator DSViter = input.begin(); DSViter != input.end(); DSViter++) {
    // DSViter is a vector of SiStripClusters located on a single module
    unsigned int ClusterId = DSViter->id();
    DetId ClusterDetId(ClusterId);
    const StripGeomDetUnit * stripdet=(const StripGeomDetUnit*)tkgeom->idToDetUnit(ClusterDetId);
    
    edmNew::DetSet<SiStripCluster>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiStripCluster>::const_iterator end  =DSViter->end();
    for(edmNew::DetSet<SiStripCluster>::const_iterator iter=begin;iter!=end;++iter) {
      //iter is a single SiStripCluster
      StripClusterParameterEstimator::LocalValues parameters=stripcpe.localParameters(*iter,*stripdet);
      
      const Surface* surface;
      surface = &(tracker->idToDet(ClusterDetId)->surface());
      LocalPoint lp = parameters.first;
      GlobalPoint gp = surface->toGlobal(lp);
      
      uint layer=0;
      unsigned int subid=ClusterDetId.subdetId();
	
      if    (subid==  StripSubdetector::TIB) {
	layer = TIBDetId(ClusterDetId).layer();
      }
      if    (subid==  StripSubdetector::TOB) {
	layer = TOBDetId(ClusterDetId).layer() + 4;
      }
      if    (subid==  StripSubdetector::TID) {
	layer = TIDDetId(ClusterDetId).wheel() + 10;
      }
      if    (subid==  StripSubdetector::TEC) {
	layer = TECDetId(ClusterDetId).wheel() + 13;
      }

      if(DEBUG) cout << "Found hit in cluster collection layer = " << layer << " with id = " << ClusterId << "   local X position = " << lp.x() << " +- " << sqrt(parameters.second.xx()) << "   matched/stereo/rphi = " << ((ClusterId & 0x3)==0) << "/" << ((ClusterId & 0x3)==1) << "/" << ((ClusterId & 0x3)==2) << endl;
    }
  }
  
  // Tracking 
  const   reco::TrackCollection *tracksCKF=trackCollectionCKF.product();
  if (DEBUG)  cout << "number ckf tracks found = " << tracksCKF->size() << endl;
  if (tracksCKF->size() == 1 ){
    if (DEBUG)    cout << "starting checking good single track event" << endl;
    reco::TrackCollection::const_iterator iCKF=trackCollectionCKF.product()->begin();
    EventTrackCKF++;  
    
    const Trajectory traject = *(TrajectoryCollectionCKF.product()->begin());
    std::vector<TrajectoryMeasurement> TMeas=traject.measurements();
    vector<TrajectoryMeasurement>::iterator itm;
    double xloc = 0.;
    double yloc = 0.;
    double xErr = 0.;
    double yErr = 0.;
    double angleX = -999.;
    double angleY = -999.;
    double xglob,yglob,zglob;
    
    for (itm=TMeas.begin();itm!=TMeas.end();itm++){
      ConstReferenceCountingPointer<TransientTrackingRecHit> theInHit;
      theInHit = (*itm).recHit();
      
      if(DEBUG) cout << "theInHit is valid = " << theInHit->isValid() << endl;
      
      uint iidd = theInHit->geographicalId().rawId();

      StripSubdetector strip=StripSubdetector(iidd);
      unsigned int subid=strip.subdetId();
      uint TKlayers = 0;
      if (subid ==  StripSubdetector::TIB) { 
	TIBDetId tibid(iidd);
	TKlayers = tibid.layer();
      }
      if (subid ==  StripSubdetector::TOB) { 
	TOBDetId tobid(iidd);
	TKlayers = tobid.layer() + 4 ; 
      }
      if (subid ==  StripSubdetector::TID) { 
	TIDDetId tidid(iidd);
	TKlayers = tidid.wheel() + 10;
      }
      if (subid ==  StripSubdetector::TEC) { 
	TECDetId tecid(iidd);
	TKlayers = tecid.wheel() + 13 ; 
      }
      if (DEBUG)	cout << "TKlayer from trajectory: " << TKlayers << "  from module = " << iidd <<  "   matched/stereo/rphi = " << ((iidd & 0x3)==0) << "/" << ((iidd & 0x3)==1) << "/" << ((iidd & 0x3)==2) << endl;

      // Make vector of TrajectoryAtValidHits to hold the trajectories
      std::vector<TrajectoryAtValidHit> TMs;
      
      // Make AnalyticalPropagator to use in TAVH constructor
      AnalyticalPropagator propagator(magField_,anyDirection); 

      // for double sided layers check both sensors--if no hit was found on either sensor surface,
      // the trajectory measurements only have one invalid hit entry on the matched surface
      // so get the TrajectoryAtValidHit for both surfaces and include them in the study
      if (isDoubleSided(iidd) &&  ((iidd & 0x3)==0) ) {
	// do hit eff check twice--once for each sensor
	//add a TM for each surface
	TMs.push_back(TrajectoryAtValidHit(*itm,tkgeom, propagator, 1));
	TMs.push_back(TrajectoryAtValidHit(*itm,tkgeom, propagator, 2));
      } else if ( isDoubleSided(iidd) && (!check2DPartner(iidd, TMeas)) ) {
      // if only one hit was found the trajectory measurement is on that sensor surface, and the other surface from
      // the matched layer should be added to the study as well
	TMs.push_back(TrajectoryAtValidHit(*itm,tkgeom, propagator, 1));
	TMs.push_back(TrajectoryAtValidHit(*itm,tkgeom, propagator, 2));
	if (DEBUG) cout << " found a hit with a missing partner" << endl;
      } else {
	//only add one TM for the single surface and the other will be added in the next iteration
	TMs.push_back(TrajectoryAtValidHit(*itm,tkgeom, propagator));
      }
      
      // Modules Constraints
      
      for(std::vector<TrajectoryAtValidHit>::const_iterator TM=TMs.begin();TM!=TMs.end();++TM) {
	
	// --> Get trajectory from combinatedState 
	iidd = TM->monodet_id();
	if (DEBUG) cout << "setting iidd = " << iidd << " before checking efficiency and ";
	
	xloc = TM->localX();
	yloc = TM->localY();
	
	xErr =  TM->localErrorX();
	yErr =  TM->localErrorY();

	angleX = atan( TM->localDxDz() );
	angleY = atan( TM->localDyDz() );

	xglob = TM->globalX();
	yglob = TM->globalY();
	zglob = TM->globalZ();
	withinAcceptance = TM->withinAcceptance();

	float YINTCons = -1.;
	if(TKlayers >= 5 && TKlayers < 11) {
	  YINTCons = 1.0;  // TOB Bonding Region cut
	}
	
	if ( (layers == TKlayers) || (layers == 0) ) {   // Look at the layer not used to reconstruct the track
	  whatlayer = TKlayers;
	  if (DEBUG)	  cout << "Looking at layer under study" << endl;
	  TrajGlbX = 0.0; TrajGlbY = 0.0; TrajGlbZ = 0.0; ModIsBad = 0; Id = 0; SiStripQualBad = 0; 
	  run = 0; event = 0; TrajLocX = 0.0; TrajLocY = 0.0; TrajLocErrX = 0.0; TrajLocErrY = 0.0; 
	  TrajLocAngleX = -999.0; TrajLocAngleY = -999.0;	ResX = 0.0; ResXSig = 0.0;
	  ClusterLocX = 0.0; ClusterLocY = 0.0; ClusterLocErrX = 0.0; ClusterLocErrY = 0.0; ClusterStoN = 0.0;

	  // RPhi RecHit Efficiency 
	  
	  if (input.size() > 0 ) {  
	    if (DEBUG) cout << "Checking clusters with size = " << input.size() << endl;
	    int nClusters = 0;
	    std::vector< std::vector<float> > VCluster_info; //fill with X residual, X residual pull, local X, sig(X), local Y, sig(Y), StoN
	    for (edmNew::DetSetVector<SiStripCluster>::const_iterator DSViter = input.begin(); DSViter != input.end(); DSViter++) {
	      // DSViter is a vector of SiStripClusters located on a single module
	      if (DEBUG)      cout << "the ID from the DSViter = " << DSViter->id() << endl; 
	      unsigned int ClusterId = DSViter->id();
	      if (ClusterId == iidd) {
		cout << "found  (ClusterId == iidd) with ClusterId = " << ClusterId << " and iidd = " << iidd << endl;
		DetId ClusterDetId(ClusterId);
		const StripGeomDetUnit * stripdet=(const StripGeomDetUnit*)tkgeom->idToDetUnit(ClusterDetId);
		
		for(edmNew::DetSet<SiStripCluster>::const_iterator iter=DSViter->begin();iter!=DSViter->end();++iter) {
		  //iter is a single SiStripCluster
		  StripClusterParameterEstimator::LocalValues parameters=stripcpe.localParameters(*iter,*stripdet);
		  float res = (parameters.first.x() - xloc);
      		  float sigma = checkConsistency(parameters , xloc, xErr);
		  // The consistency is probably more accurately measured with the Chi2MeasurementEstimator. To use it
		  // you need a TransientTrackingRecHit instead of the cluster
		  //theEstimator=       new Chi2MeasurementEstimator(30);
		  //const Chi2MeasurementEstimator *theEstimator(100);
		  //theEstimator->estimate(TM->tsos(), TransientTrackingRecHit);

		  //SiStripClusterInfo clusterInfo = SiStripClusterInfo(*iter, es);  
		  // signal to noise from SiStripClusterInfo not working in 225. I'll fix this after the interface
		  // redesign in 300 -ku
		  //float cluster_info[7] = {res, sigma, parameters.first.x(), sqrt(parameters.second.xx()), parameters.first.y(), sqrt(parameters.second.yy()), signal_to_noise};
		  std::vector< float > cluster_info;
		  cluster_info.push_back(res); 
		  cluster_info.push_back(sigma);
		  cluster_info.push_back(parameters.first.x()); 
		  cluster_info.push_back(sqrt(parameters.second.xx()));
		  cluster_info.push_back(parameters.first.y());
		  cluster_info.push_back(sqrt(parameters.second.yy()));
		  cout << "before getting signal over noise" << endl;
		  //cluster_info.push_back( clusterInfo.signalOverNoise() );
		  //cluster_info.push_back( clusterInfo.getSignalOverNoise() );
		  cout << "after getting signal over noise" << endl;
		  VCluster_info.push_back(cluster_info);
		  nClusters++;
		  if (DEBUG) cout << "Have ID match. residual = " << VCluster_info.back()[0] << "  res sigma = " << VCluster_info.back()[1] << endl;
		  cout << "trajectory measurement compatability estimate = " << (*itm).estimate() << endl;
		  cout << "hit position = " << parameters.first.x() << "  hit error = " << sqrt(parameters.second.xx()) << "  trajectory position = " << xloc << "  traj error = " << xErr << endl;
		}
	      }
	    }
	    float FinalResSig = 1000.0;
	    float FinalCluster[7]= {1000.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	    if (nClusters > 0) {
	      if (DEBUG) cout << "found clusters > 0" << endl;
	      if (nClusters > 1) {
		//get the smallest one
		vector< vector<float> >::iterator ires;
		for (ires=VCluster_info.begin(); ires!=VCluster_info.end(); ires++){
		  if ( abs((*ires)[1]) < abs(FinalResSig)) {
		    FinalResSig = (*ires)[1];
		    for (uint i = 0; i<ires->size(); i++) {
		      cout << "filling final cluster. i = " << i << " before fill FinalCluster[i]=" << FinalCluster[i] << " and (*ires)[i] =" << (*ires)[i] << endl;
		      FinalCluster[i] = (*ires)[i];
		      cout << "filling final cluster. i = " << i << " after fill FinalCluster[i]=" << FinalCluster[i] << " and (*ires)[i] =" << (*ires)[i] << endl;
		    }
		  }
		  if (DEBUG) cout << "iresidual = " << (*ires)[0] << "  isigma = " << (*ires)[1] << "  and FinalRes = " << FinalCluster[0] << endl;
		}
	      }
	      else {
		FinalResSig = VCluster_info.at(0)[1];
		for (uint i = 0; i<VCluster_info.at(0).size(); i++) {
		  FinalCluster[i] = VCluster_info.at(0)[i];
		}
	      }
	      nClusters=0;
	      VCluster_info.clear();
	    }
	    
	    if (DEBUG) cout << "Final residual in X = " << FinalCluster[0] << "+-" << (FinalCluster[0]/FinalResSig) << endl;
	    if (DEBUG) cout << "Checking location of trajectory: abs(yloc) = " << abs(yloc) << "  abs(xloc) = " << abs(xloc) << endl;
	    if (DEBUG) cout << "Checking location of cluster hit: yloc = " << FinalCluster[4] << "+-" << FinalCluster[5] << "  xloc = " << FinalCluster[2] << "+-" << FinalCluster[3] << endl;
	    if (DEBUG) cout << "Final cluster signal to noise = " << FinalCluster[6] << endl;
	    
	    // YINTCons is the bonding region absolute cut
	    if ( abs(yloc) < YINTCons ) { 
	      withinAcceptance = false;
	      cout << "setting withinAcceptance false from glued region and value is now = " << withinAcceptance << endl;
	    }

	    // fill ntuple varibles
	    //get global position from module id number iidd
	    TrajGlbX = xglob;
	    TrajGlbY = yglob;
	    TrajGlbZ = zglob;	  
	    Id = iidd;
	    run = run_nr;
	    event = ev_nr;
	    //if ( SiStripQuality_->IsModuleBad(iidd) ) {
	    if ( SiStripQuality_->getBadApvs(iidd)!=0 ) {
	      SiStripQualBad = 1; 
	      if(DEBUG) cout << "strip is bad from SiStripQuality" << endl;
	    } else {
	      SiStripQualBad = 0; 
	      if(DEBUG) cout << "strip is good from SiStripQuality" << endl;
	    }
	    
	    TrajLocX = xloc;
	    TrajLocY = yloc;
	    TrajLocErrX = xErr;
	    TrajLocErrY = yErr;
	    TrajLocAngleX = angleX;
	    TrajLocAngleY = angleY;
	    ResX = FinalCluster[0];
	    ResXSig = FinalResSig;
	    if (FinalResSig != FinalCluster[1]) cout << "Problem with best cluster selection because FinalResSig = " << FinalResSig << " and FinalCluster[1] = " << FinalCluster[1] << endl;
	    ClusterLocX = FinalCluster[2];
	    ClusterLocY = FinalCluster[4];
	    ClusterLocErrX = FinalCluster[3];
	    ClusterLocErrY = FinalCluster[5];
	    ClusterStoN = FinalCluster[6];
	    
	    if (DEBUG)	      cout << "before check good" << endl;
	    
	    if ( FinalResSig < 999.0) {  //could make requirement on track/hit consistency, but for
	      //now take anything with a hit on the module
	      cout << "hit being counted as good" << endl;
	      ModIsBad = 0;
	      traj->Fill();
	    }
	    else {
	      cout << "hit being counted as bad   ######### Invalid RPhi FinalResX " << FinalCluster[0] << " FinalRecHit " << 
		iidd << "   TKlayers  "  <<  TKlayers  << " xloc " <<  xloc << " yloc  " << yloc << " module " << iidd << 
		"   matched/stereo/rphi = " << ((iidd & 0x3)==0) << "/" << ((iidd & 0x3)==1) << "/" << ((iidd & 0x3)==2) << endl;
	      ModIsBad = 1;
	      traj->Fill();
	      
	      cout << " RPhi Error " << sqrt(xErr*xErr + yErr*yErr) << " ErrorX " << xErr << " yErr " <<  yErr <<  endl;
	    } if (DEBUG) cout << "after good location check" << endl;
	  } if (DEBUG) cout << "after list of clusters" << endl;
	} if (DEBUG) cout << "After layers=TKLayers if" << endl;
      } if (DEBUG) cout << "After looping over TrajAtValidHit list" << endl;
    } if (DEBUG) cout << "end TMeasurement loop" << endl;
  }
}

void HitEff::endJob(){
  traj->GetDirectory()->cd();
  traj->Write();  
  
  cout << " Events Analysed             " <<  events          << endl;
  cout << " Number Of Tracked events    " <<  EventTrackCKF   << endl;
}

double HitEff::checkConsistency(StripClusterParameterEstimator::LocalValues parameters, double xx, double xerr) {
  double error = sqrt(parameters.second.xx() + xerr*xerr);
  double separation = abs(parameters.first.x() - xx);
  double consistency = separation/error;
  return consistency;
}

double HitEff::checkConsistency(const SiStripRecHit2D* rechit, double xx, double xerr)
{
  double error = sqrt(rechit->localPositionError().xx() + xerr*xerr);
  double separation = abs(rechit->localPosition().x() - xx);
  double consistency = separation/error;
  return consistency;
}

bool HitEff::isDoubleSided(uint iidd) const {
  StripSubdetector strip=StripSubdetector(iidd);
  unsigned int subid=strip.subdetId();
  uint layer = 0;
  if (subid ==  StripSubdetector::TIB) { 
    TIBDetId tibid(iidd);
    layer = tibid.layer();
    if (layer == 1 || layer == 2) return true;
    else return false;
  }
  else if (subid ==  StripSubdetector::TOB) { 
    TOBDetId tobid(iidd);
    layer = tobid.layer() + 4 ; 
    if (layer == 5 || layer == 6) return true;
    else return false;
  }
  else if (subid ==  StripSubdetector::TID) { 
    TIDDetId tidid(iidd);
    layer = tidid.ring() + 10;
    if (layer == 11 || layer == 12) return true;
    else return false;
  }
  else if (subid ==  StripSubdetector::TEC) { 
    TECDetId tecid(iidd);
    layer = tecid.ring() + 13 ; 
    if (layer == 14 || layer == 15 || layer == 18) return true;
    else return false;
  }
  else
    return false;
}

bool HitEff::check2DPartner(uint iidd, std::vector<TrajectoryMeasurement> traj) {
  uint partner_iidd = 0;
  bool found2DPartner = false;
  // first get the id of the other detector
  if ((iidd & 0x3)==1) partner_iidd = iidd+1;
  if ((iidd & 0x3)==2) partner_iidd = iidd-1;
  // next look in the trajectory measurements for a measurement from that detector
  // loop through trajectory measurements to find the partner_iidd
  for (std::vector<TrajectoryMeasurement>::const_iterator iTM=traj.begin(); iTM!=traj.end(); ++iTM) {
    if (iTM->recHit()->geographicalId().rawId()==partner_iidd) {
      found2DPartner = true;
    }
  }
  return found2DPartner;
}

//  const edmNew::DetSetVector<SiStripCluster>& input = *theClusters;

//define this as a plug-in
DEFINE_FWK_MODULE(HitEff);
