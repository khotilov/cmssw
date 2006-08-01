#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "RecoParticleFlow/PFClusterAlgo/interface/PFClusterAlgo.h"
#include "RecoParticleFlow/PFAlgo/interface/PFBlock.h"
#include "RecoParticleFlow/PFAlgo/interface/PFBlockElement.h"
#include "RecoParticleFlow/PFRootEvent/interface/PFRootEventManager.h"
#include "RecoParticleFlow/PFRootEvent/interface/IO.h"

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TMarker.h>
#include <TH2F.h>
#include <TCutG.h>
#include <TPolyLine.h>
#include <TColor.h>
#include "TGraph.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include "TVector3.h"

#include <iostream>

using namespace std;



PFRootEventManager::PFRootEventManager() {
  graphTrack_.resize(NViews);
}


PFRootEventManager::~PFRootEventManager() {
  // Clear TGraph if needed
  for (unsigned iView = 0; iView < graphTrack_.size(); iView++) {
    if (graphTrack_[iView].size()) {
      for (unsigned iGraph = 0; iGraph < graphTrack_[iView].size(); iGraph++)
	delete graphTrack_[iView][iGraph];
      graphTrack_[iView].clear();
    }
  }
}


PFRootEventManager::PFRootEventManager(const char* file) {  
  options_ = 0;
  readOptions(file);
  iEvent_=0;
  displayHistEtaPhi_=0;
  displayView_.resize(NViews);
  displayHist_.resize(NViews);
  for (unsigned iView = 0; iView < NViews; iView++) {
    displayView_[iView] = 0;
    displayHist_[iView] = 0;
  }
  graphTrack_.resize(NViews);
  maxERecHitEcal_ = -1;
  maxERecHitHcal_ = -1;
}


void PFRootEventManager::reset() { 
  maxERecHitEcal_ = -1;
  maxERecHitHcal_ = -1;  
  clusters_.clear();
  rechits_.clear();
  recTracks_.clear();
}


void PFRootEventManager::readOptions(const char* file, bool refresh) {

  if( !options_ )
    options_ = new IO(file);
  else if( refresh) {
    delete options_;
    options_ = new IO(file);
  }

  options_->GetOpt("root","file", inFileName_);
  
  file_ = TFile::Open(inFileName_.c_str() );
  if(file_->IsZombie() ) {
    return;
  }

  tree_ = (TTree*) file_->Get("Events");
  
  string hitbranchname;
  options_->GetOpt("root","hits_branch", hitbranchname);
  
  hitsBranch_ = tree_->GetBranch(hitbranchname.c_str());
  if(!hitsBranch_) {
    cerr<<"PFRootEventManager::ReadOptions : branch not found : "
	<<hitbranchname<<endl;
  }
  else {
    hitsBranch_->SetAddress(&rechits_);
    // cout << "Set branch " << hitbranchname << endl;
  }

  string recTracksbranchname;
  options_->GetOpt("root","recTracks_branch", recTracksbranchname);

  recTracksBranch_ = tree_->GetBranch(recTracksbranchname.c_str());
  if(!recTracksBranch_) {
    cerr <<"PFRootEventManager::ReadOptions : branch not found : "
	 << recTracksbranchname << endl;
  }
  else {
    recTracksBranch_->SetAddress(&recTracks_);
    // cout << "Set branch " << recTracksbranchname << " " << recTracksBranch_->GetEntries() << endl;
  }    

  vector<int> algos;
  options_->GetOpt("display", "algos", algos);
  algosToDisplay_.clear();
  for(unsigned i=0; i< algos.size(); i++) algosToDisplay_.insert( algos[i] );
  typedef map<int, TCanvas * >::iterator IT;
  for(IT it = displayEtaPhi_.begin(); it!=displayEtaPhi_.end(); it++) {
    if( algosToDisplay_.find(it->first) == algosToDisplay_.end() ) {
      it->second->Close();
      displayEtaPhi_.erase(it);
      if(displayEtaPhi_.empty() ) break;
      it++;
    }    
  }  

  viewSizeEtaPhi_.clear();
  options_->GetOpt("display", "viewsize_etaphi", viewSizeEtaPhi_);
  if(viewSizeEtaPhi_.size() != 2) {
    cerr<<"PFRootEventManager::ReadOptions, bad display/viewsize_etaphi tag...using 700/350"
	<<endl;
    viewSizeEtaPhi_.clear();
    viewSizeEtaPhi_.push_back(700); 
    viewSizeEtaPhi_.push_back(350); 
  }

  viewSize_.clear();
  options_->GetOpt("display", "viewsize_xy", viewSize_);
  if(viewSize_.size() != 2) {
    cerr<<"PFRootEventManager::ReadOptions, bad display/viewsize_xy tag...using 700/350"
	<<endl;
    viewSize_.clear();
    viewSize_.push_back(600); 
    viewSize_.push_back(600); 
  }

  displayColorClusters_ = false;
  options_->GetOpt("display", "color_clusters", displayColorClusters_);

  displayZoomFactor_ = 10;  
  options_->GetOpt("display", "zoom_factor", displayZoomFactor_);

  threshEcalBarrel_ = 0.1;
  options_->GetOpt("clustering", "thresh_Ecal_Barrel", threshEcalBarrel_);
  
  threshSeedEcalBarrel_ = 0.3;
  options_->GetOpt("clustering", "thresh_Seed_Ecal_Barrel", 
		   threshSeedEcalBarrel_);

  threshEcalEndcap_ = 0.2;
  options_->GetOpt("clustering", "thresh_Ecal_Endcap", threshEcalEndcap_);

  threshSeedEcalEndcap_ = 0.8;
  options_->GetOpt("clustering", "thresh_Seed_Ecal_Endcap",
		   threshSeedEcalEndcap_);


  nNeighboursEcal_ = 4;
  options_->GetOpt("clustering", "neighbours_Ecal", nNeighboursEcal_);


  int dcormode = -1;
  options_->GetOpt("clustering", "depthCor_Mode", dcormode);
  
  double dcora = -1;
  options_->GetOpt("clustering", "depthCor_A", dcora);
  double dcorb = -1;
  options_->GetOpt("clustering", "depthCor_B", dcorb);
  double dcorap = -1;
  options_->GetOpt("clustering", "depthCor_A_preshower", dcorap);
  double dcorbp = -1;
  options_->GetOpt("clustering", "depthCor_B_preshower", dcorbp);

  if( dcormode > -0.5 && 
      dcora > -0.5 && 
      dcorb > -0.5 && 
      dcorap > -0.5 && 
      dcorbp > -0.5 )
    reco::PFCluster::setDepthCorParameters( dcormode, 
					    dcora, dcorb, 
					    dcorap, dcorbp);


  

  threshHcalBarrel_ = 0.1;
  options_->GetOpt("clustering", "thresh_Hcal_Barrel", threshHcalBarrel_);
  
  threshSeedHcalBarrel_ = 0.3;
  options_->GetOpt("clustering", "thresh_Seed_Hcal_Barrel", 
		   threshSeedHcalBarrel_);

  threshHcalEndcap_ = 0.2;
  options_->GetOpt("clustering", "thresh_Hcal_Endcap", threshHcalEndcap_);

  threshSeedHcalEndcap_ = 0.8;
  options_->GetOpt("clustering", "thresh_Seed_Hcal_Endcap",
		   threshSeedHcalEndcap_);

  nNeighboursHcal_ = 4;
  options_->GetOpt("clustering", "neighbours_Hcal", nNeighboursHcal_);


  // options for particle flow

  string map_ECAL_eta;
  options_->GetOpt("particle_flow", "resolution_map_ECAL_eta", map_ECAL_eta);
  string map_ECAL_phi;
  options_->GetOpt("particle_flow", "resolution_map_ECAL_phi", map_ECAL_phi);
  string map_ECALec_x;
  options_->GetOpt("particle_flow", "resolution_map_ECALec_x", map_ECALec_x);
  string map_ECALec_y;
  options_->GetOpt("particle_flow", "resolution_map_ECALec_y", map_ECALec_y);
  string map_HCAL_eta;
  options_->GetOpt("particle_flow", "resolution_map_HCAL_eta", map_HCAL_eta);
  string map_HCAL_phi;
  options_->GetOpt("particle_flow", "resolution_map_HCAL_phi", map_HCAL_phi);
  PFBlock::SetResMaps(map_ECAL_eta,
		      map_ECAL_phi, 
		      map_ECALec_x,
		      map_ECALec_y,
		      map_HCAL_eta,
		      map_HCAL_phi);

  double chi2_ECAL_HCAL=0;
  options_->GetOpt("particle_flow", "chi2_ECAL_HCAL", chi2_ECAL_HCAL);
  double chi2_ECAL_PS=0;
  options_->GetOpt("particle_flow", "chi2_ECAL_PS", chi2_ECAL_PS);
  double chi2_HCAL_PS=0;
  options_->GetOpt("particle_flow", "chi2_HCAL_PS", chi2_HCAL_PS);
  double chi2_ECAL_Track=100;
  options_->GetOpt("particle_flow", "chi2_ECAL_Track", chi2_ECAL_Track);
  double chi2_HCAL_Track=100;
  options_->GetOpt("particle_flow", "chi2_HCAL_Track", chi2_HCAL_Track);
  double chi2_PS_Track=0;
  options_->GetOpt("particle_flow", "chi2_PS_Track", chi2_PS_Track);

  PFBlock::SetMaxChi2(chi2_ECAL_HCAL,
		      chi2_ECAL_PS,
		      chi2_HCAL_PS,
		      chi2_ECAL_Track,
		      chi2_HCAL_Track,
		      chi2_PS_Track );
  double nsigma;
  options_->GetOpt("particle_flow", "nsigma_neutral", nsigma);
  PFBlock::SetNsigmaNeutral(nsigma);

  vector<double> ecalib;
  options_->GetOpt("particle_flow", "ecalib", ecalib);
  if(ecalib.size() == 2) {
    PFBlock::SetEcalib(ecalib[0], ecalib[1]);
  }
  else {
    PFBlock::SetEcalib(0, 1);
  }

  reconMethod_ = 3; 
  options_->GetOpt("particle_flow", "recon_method", reconMethod_);

  displayJetColors_ = false;
  options_->GetOpt("display", "jet_colors", displayJetColors_);
  

}


bool PFRootEventManager::processEntry(int entry) {

  reset();

  cout<<"process entry "<< entry << endl;
  
  if(hitsBranch_) hitsBranch_->GetEntry(entry);
  if(recTracksBranch_) recTracksBranch_->GetEntry(entry);

  cout<<"number of rechits : "<<rechits_.size()<<endl;
  cout<<"number of recTracks : "<<recTracks_.size()<<endl;

  clustering(); 
  particleFlow();

  return false;
}


void PFRootEventManager::clustering() {
  
  std::map<unsigned,  reco::PFRecHit* > rechits;
   
  for(unsigned i=0; i<rechits_.size(); i++) {
    rechits_[i].CalculatePositionREP();
    rechits.insert( make_pair(rechits_[i].getDetId(), &rechits_[i] ) );
  }

  for( PFClusterAlgo::IDH ih = rechits.begin(); ih != rechits.end(); ih++) {
    ih->second->findPtrsToNeighbours( rechits );
  }

  PFClusterAlgo clusteralgo; 

  clusteralgo.setThreshEcalBarrel( threshEcalBarrel_ );
  clusteralgo.setThreshSeedEcalBarrel( threshSeedEcalBarrel_ );
  
  clusteralgo.setThreshEcalEndcap( threshEcalEndcap_ );
  clusteralgo.setThreshSeedEcalEndcap( threshSeedEcalEndcap_ );

  clusteralgo.setNNeighboursEcal( nNeighboursEcal_  );
  
  clusteralgo.setThreshHcalBarrel( threshHcalBarrel_ );
  clusteralgo.setThreshSeedHcalBarrel( threshSeedHcalBarrel_ );
  
  clusteralgo.setThreshHcalEndcap( threshHcalEndcap_ );
  clusteralgo.setThreshSeedHcalEndcap( threshSeedHcalEndcap_ );

  clusteralgo.setNNeighboursHcal( nNeighboursHcal_ );

  clusteralgo.init( rechits ); 
  clusteralgo.allClusters();

  clusters_ = clusteralgo.getClusters();
}


void PFRootEventManager::particleFlow() {
  
  // clear stuff from previous call

  for(PFBlock::IT ie = allElements_.begin(); ie!=  allElements_.end(); ie++ ) {
    delete *ie;
  } 
  allElements_.clear();

  allPFBs_.clear();

  // create PFBlockElements from clusters and rectracks

  for(unsigned i=0; i<clusters_.size(); i++) {
   
    if( clusters_[i].getType() != reco::PFCluster::TYPE_PF ) continue;

    int layer = clusters_[i].getLayer();
      
    switch( layer ) {
    case PFLayer::ECAL_BARREL:
    case PFLayer::ECAL_ENDCAP:
      allElements_.insert( new PFBlockElementECAL( & clusters_[i] ) );
      break;
    case PFLayer::PS1:
    case PFLayer::PS2:
      allElements_.insert( new PFBlockElementPS( & clusters_[i] ) );
      break;
    case PFLayer::HCAL_BARREL1:
    case PFLayer::HCAL_BARREL2:
    case PFLayer::HCAL_ENDCAP:
      allElements_.insert( new PFBlockElementHCAL( & clusters_[i] ) );
      break;
    default:
      break;
    }    
  }



  for(unsigned i=0; i<recTracks_.size(); i++) {
    recTracks_[i].calculatePositionREP();
    allElements_.insert( new PFBlockElementTrack( & recTracks_[i] ) );  
  }



  PFBlock::SetAllElements( allElements_ );
  int efbcolor = 2;
  for(PFBlock::IT iele = allElements_.begin(); 
      iele != allElements_.end(); iele++) {

    if( (*iele) -> GetPFBlock() ) continue; // already associated
    
    allPFBs_.push_back( PFBlock() );
    allPFBs_.back().Associate( 0, *iele );
    
    if( displayJetColors_ ) efbcolor = 1;
    
    allPFBs_.back().Finalize(efbcolor, reconMethod_); 
//     cout<<"new eflowblock----------------------"<<endl;
//     cout<<allPFBs_.back()<<endl;
    efbcolor++;
  }

  
  for(unsigned iefb = 0; iefb<allPFBs_.size(); iefb++) {

    switch(reconMethod_) {
    case 1:
      allPFBs_[iefb].ReconstructParticles1();
      break;
    case 2:
      allPFBs_[iefb].ReconstructParticles2();
      break;
    case 3:
      allPFBs_[iefb].ReconstructParticles3();
      break;
    default:
      break;
    }    
    cout<<(allPFBs_[iefb])<<endl;
  }

}


void PFRootEventManager::display(int ientry) {
  
  processEntry(ientry);
  //  displayEtaPhi();
  displayView(RZ);
  displayView(XY);
  displayView(EPE);
  displayView(EPH);
}


void PFRootEventManager::displayView(unsigned viewType) {
  
  // Clear TGraph if needed
  if (graphTrack_[viewType].size()) {
    for (unsigned iGraph = 0; iGraph < graphTrack_[viewType].size(); iGraph++)
      delete graphTrack_[viewType][iGraph];
    graphTrack_[viewType].clear();
  }

  // display or clear canvas
  if(!displayView_[viewType] || !gROOT->GetListOfCanvases()->FindObject(displayView_[viewType]) ) {
    assert(viewSize_.size() == 2);
    switch(viewType) {
    case XY:
      displayView_[viewType] = new TCanvas("displayXY_", "XY view",
					   viewSize_[0], viewSize_[1]);
      break;
    case RZ:
      displayView_[viewType] = new TCanvas("displayRZ_", "RZ view",
					   viewSize_[0], viewSize_[1]);
      break;
    case EPE:
      displayView_[viewType] = new TCanvas("displayEPE_", "eta/phi view, ECAL",
					   viewSize_[0]*2, viewSize_[1]);
      break;
    case EPH:
      displayView_[viewType] = new TCanvas("displayEPH_", "eta/phi view, HCAL",
					   viewSize_[0]*2, viewSize_[1]);
      break;
    }
    displayView_[viewType]->SetGrid(0, 0);
    displayView_[viewType]->SetLeftMargin(0.12);
    displayView_[viewType]->SetBottomMargin(0.12);
    displayView_[viewType]->Draw();  
  } else 
    displayView_[viewType]->Clear();
  displayView_[viewType]->cd();
  
  // Draw support histogram
  double zLow = -500.;
  double zUp  = +500.;
  double rLow = -300.;
  double rUp  = +300.;
  if(!displayHist_[viewType]) {
    switch(viewType) {
    case XY:
      displayHist_[viewType] = new TH2F("hdisplayHist_XY", "", 
					500, rLow, rUp, 500, rLow, rUp);
      displayHist_[viewType]->SetXTitle("X");
      displayHist_[viewType]->SetYTitle("Y");
      break;
    case RZ:
      displayHist_[viewType] = new TH2F("hdisplayHist_RZ", "", 
					500, zLow, zUp, 500, rLow, rUp);
      displayHist_[viewType]->SetXTitle("Z");
      displayHist_[viewType]->SetYTitle("R");
      break;      
    case EPE:
    case EPH:
      if(! displayHist_[EPE] ) {
	displayHist_[EPE] = new TH2F("hdisplayHist_EP", "", 
				     500, -5, 5, 500, -3.5, 3.5);
	displayHist_[EPE]->SetXTitle("#eta");
	displayHist_[EPE]->SetYTitle("#phi");
      }
      displayHist_[EPH] = displayHist_[EPE];
      break;
    default:
      std::cerr << "This kind of view is not implemented" << std::endl;
      break;
    }
    displayHist_[viewType]->SetStats(kFALSE);
  }
  displayHist_[viewType]->Draw();

  switch(viewType) {
  case XY:
    { 
      // Draw ECAL front face
      frontFaceECALXY_.SetX1(0);
      frontFaceECALXY_.SetY1(0);
      // COLIN: should not be hardcoded. 
      // either attribute of this class, or create a 
      // SimpleGeometry or something like this. 
      // NB in PFProduced ECAL radius is set to 130 not 129.
      // same remarks for HCAL
      frontFaceECALXY_.SetR1(129); // <==== BE CAREFUL, ECAL size is hardcoded !!!!
      frontFaceECALXY_.SetR2(129);
      frontFaceECALXY_.SetFillStyle(0);
      frontFaceECALXY_.Draw();
      
      // Draw HCAL front face
      frontFaceHCALXY_.SetX1(0);
      frontFaceHCALXY_.SetY1(0);
      frontFaceHCALXY_.SetR1(183); // <==== BE CAREFUL, HCAL size is hardcoded !!!!
      frontFaceHCALXY_.SetR2(183);
      frontFaceHCALXY_.SetFillStyle(0);
      frontFaceHCALXY_.Draw();
      break;
    }
  case RZ:
    {
      // Draw lines at different etas
      TLine l;
      l.SetLineColor(1);
      l.SetLineStyle(3);
      TLatex etaLeg;
      etaLeg.SetTextSize(0.02);
      float etaMin = -3.;
      float etaMax = +3.;
      float etaBin = 0.2;
      int nEtas = int((etaMax - etaMin)/0.2) + 1;
      for (int iEta = 0; iEta <= nEtas; iEta++) {
	float eta = etaMin + iEta*etaBin;
	float r = 0.9*rUp;
	TVector3 etaImpact;
	etaImpact.SetPtEtaPhi(r, eta, 0.);
	etaLeg.SetTextAlign(21);
	if (eta <= -1.39) {
	  etaImpact.SetXYZ(0.,0.85*zLow*tan(etaImpact.Theta()),0.85*zLow);
	  etaLeg.SetTextAlign(31);
	} else if (eta >= 1.39) {
	  etaImpact.SetXYZ(0.,0.85*zUp*tan(etaImpact.Theta()),0.85*zUp);
	  etaLeg.SetTextAlign(11);
	}
	l.DrawLine(0., 0., etaImpact.Z(), etaImpact.Perp());
	etaLeg.DrawLatex(etaImpact.Z(), etaImpact.Perp(), Form("%2.1f", eta));
      }
      
      frontFaceECALRZ_.SetX1(-303.16);
      frontFaceECALRZ_.SetY1(-129.);
      frontFaceECALRZ_.SetX2(303.16);
      frontFaceECALRZ_.SetY2(129.);
      frontFaceECALRZ_.SetFillStyle(0);
      frontFaceECALRZ_.Draw();
      break;
    }
  default:
    // do nothing for other views
    break;
  }

  // Put highest Pt reconstructed track along phi = pi/2
  double phi0 = 0.;
//   double ptMax = 0.;
//   double pxMax = 0.;
//   double pyMax = 0.;
//   std::vector<reco::PFRecTrack>::iterator itRecTrack;
//   for (itRecTrack = recTracks_.begin(); itRecTrack != recTracks_.end();
//        itRecTrack++) {
//     double pt = itRecTrack->trajectoryPoint(reco::PFTrajectoryPoint::ClosestApproach).momentum().Pt();
//     if (pt > ptMax) {
//       ptMax = pt;
//       phi0 = itRecTrack->trajectoryPoint(reco::PFTrajectoryPoint::ClosestApproach).momentum().Phi();
//       pxMax = itRecTrack->trajectoryPoint(reco::PFTrajectoryPoint::ClosestApproach).momentum().Px();
//       pyMax = itRecTrack->trajectoryPoint(reco::PFTrajectoryPoint::ClosestApproach).momentum().Py();
//     }
//   }


  // display reconstructed objects
  displayView_[viewType]->cd();
  displayRecHits(viewType, phi0);
  displayRecTracks(viewType, phi0);
  displayClusters(viewType, phi0);
}



void PFRootEventManager::displayRecHits(unsigned viewType, double phi0) 
{
  double maxee = getMaxEEcal();
  double maxeh = getMaxEHcal();
  double maxe = maxee>maxeh ? maxee : maxeh;
  
  std::vector<reco::PFRecHit>::iterator itRecHit;
  for (itRecHit = rechits_.begin(); itRecHit != rechits_.end(); 
       itRecHit++) {
    // double maxe = 0;
    double thresh = 0;
    switch(itRecHit->getLayer()) {
    case PFLayer::ECAL_BARREL:
      // maxe = maxee;
      thresh = threshEcalBarrel_;
      break;     
    case PFLayer::ECAL_ENDCAP:
      // maxe = maxee;
      thresh = threshEcalEndcap_;
      break;     
    case PFLayer::HCAL_BARREL1:
    case PFLayer::HCAL_BARREL2:
      // maxe = maxeh;
      thresh = threshHcalBarrel_;
      break;           
    case PFLayer::HCAL_ENDCAP:
      // maxe = maxeh;
      thresh = threshHcalEndcap_;
      break;           
    default:
      cerr<<"manage other layers"<<endl;
      continue;
    }

    if(itRecHit->getEnergy() > thresh )
      displayRecHit(*itRecHit, viewType, maxe, thresh, phi0);
  }
}

void PFRootEventManager::displayRecHit(reco::PFRecHit& rh, unsigned viewType,
				       double maxe, double thresh,
				       double phi0) 
{

  int layer = rh.getLayer();

  // on EPH view, draw only ECAL and preshower
  if(  viewType == EPH && 
       ! (layer == PFLayer::HCAL_BARREL1 || 
	  layer == PFLayer::HCAL_BARREL2 || 
	  layer == PFLayer::HCAL_ENDCAP ) ) return;
  
  // on EPE view, draw only HCAL and preshower
  if(  viewType == EPE && 
       (layer == PFLayer::HCAL_BARREL1 || 
	layer == PFLayer::HCAL_BARREL2 || 
	layer == PFLayer::HCAL_ENDCAP ) ) return;
  
  

  // math::XYZPoint vPhi0(cos(phi0), sin(phi0), 0.);
  if (rh.getEnergy() < thresh ) return;


  double rheta = rh.getPositionREP().Eta();
  double rhphi = rh.getPositionREP().Phi();
  double sign = 1.;
  if (cos(phi0 - rhphi) < 0.) sign = -1.;

  TCutG* cutg = (TCutG*) gROOT->FindObject("CUTG");
  if (cutg) {
    if( !cutg->IsInside(rheta, rhphi) ) return;
  }

  double etaSize[4];
  double phiSize[4];
  double x[5];
  double y[5];
  double z[5];
  double r[5];
  double eta[5];
  double phi[5];
  double xprop[5];
  double yprop[5];
  double etaprop[5];
  double phiprop[5];

  
  const std::vector< math::XYZPoint >& corners = rh.getCornersXYZ();
  
  assert(corners.size() == 4);
  
  double propfact = 0.95; // so that the cells don't overlap ? 
  double ampl = (log(rh.getEnergy() + 1.)/log(maxe + 1.));
  
  for ( unsigned jc=0; jc<4; ++jc ) { 
    phiSize[jc] = rhphi-corners[jc].Phi();
    etaSize[jc] = rheta-corners[jc].Eta();
    if ( phiSize[jc] > 1. ) phiSize[jc] -= 2.*TMath::Pi();  // this is strange...
    if ( phiSize[jc] < -1. ) phiSize[jc]+= 2.*TMath::Pi();
    phiSize[jc] *= propfact;
    etaSize[jc] *= propfact;

    math::XYZPoint cornerposxyz = corners[jc];

//     math::XYZPoint cornerposxyz(corners[jc].X()*vPhi0.Y() - 
// 				corners[jc].Y()*vPhi0.X(), 
// 				corners[jc].X()*vPhi0.X() + 
// 				corners[jc].Y()*vPhi0.Y(), 
// 				corners[jc].Z());
    // cornerposxyz.SetPhi( corners[jc]->Y() - phi0 );

    x[jc] = cornerposxyz.X();
    y[jc] = cornerposxyz.Y();
    z[jc] = cornerposxyz.Z();
    r[jc] = sign*cornerposxyz.Rho();
    eta[jc] = rheta + etaSize[jc];
    phi[jc] = rhphi + phiSize[jc];
    

    // cell area is prop to log(E)
    // not drawn for preshower. 
    // otherwise, drawn for eta/phi view, and for endcaps in xy view
    if( 
       layer != PFLayer::PS1 && 
       layer != PFLayer::PS2 && 
       ( viewType == EPE || viewType == EPH || 
	 ( viewType == XY &&  
	   ( layer == PFLayer::ECAL_ENDCAP && 
	     layer == PFLayer::HCAL_ENDCAP ) ) ) ) {

//        !(layer == PFLayer::ECAL_BARREL || 
// 	 layer == PFLayer::HCAL_BARREL1 || 
// 	 layer == PFLayer::HCAL_BARREL2)
//       ) {
      
//       math::XYZPoint centreXYZrot(rh.getPositionXYZ().X()*vPhi0.Y() - 
// 				  rh.getPositionXYZ().Y()*vPhi0.X(), 
// 				  rh.getPositionXYZ().X()*vPhi0.X() + 
// 				  rh.getPositionXYZ().Y()*vPhi0.Y(), 
// 				  rh.getPositionXYZ().Z());
      // centreXYZrot.SetPhi( fCentre.Y() - phi0 );
      
      math::XYZPoint centreXYZrot = rh.getPositionXYZ();

      math::XYZPoint centertocorner(x[jc] - centreXYZrot.X(), 
				    y[jc] - centreXYZrot.Y(),
				    0 );

      math::XYZPoint centertocornerep(eta[jc] - centreXYZrot.Eta(), 
				      phi[jc] - centreXYZrot.Phi(),
				      0 );
      

      // centertocorner -= centreXYZrot;
      xprop[jc] = centreXYZrot.X() + centertocorner.X()*ampl;
      yprop[jc] = centreXYZrot.Y() + centertocorner.Y()*ampl;

      etaprop[jc] = centreXYZrot.Eta() + centertocornerep.X()*ampl;
      phiprop[jc] = centreXYZrot.Phi() + centertocornerep.Y()*ampl;
    }
  }

  if(layer == PFLayer::ECAL_BARREL  || 
     layer == PFLayer::HCAL_BARREL1 || 
     layer == PFLayer::HCAL_BARREL2 || viewType == RZ) {

    // we are in the barrel. Determining which corners to shift 
    // away from the center to represent the cell energy
    
    // this has to be done also for endcap rz !
    int i1 = -1;
    int i2 = -1;

    if(fabs(phiSize[1]-phiSize[0]) > 0.0001) {
      if (viewType == XY) {
	i1 = 2;
	i2 = 3;
      } else if (viewType == RZ) {
	i1 = 1;
	i2 = 2;
      }
    } else {
      if (viewType == XY) {
	i1 = 1;
	i2 = 2;
      } else if (viewType == RZ) {
	i1 = 2;
	i2 = 3;
      }
    }

    x[i1] *= 1+ampl/2.;
    x[i2] *= 1+ampl/2.;
    y[i1] *= 1+ampl/2.;
    y[i2] *= 1+ampl/2.;
    z[i1] *= 1+ampl/2.;
    z[i2] *= 1+ampl/2.;
    r[i1] *= 1+ampl/2.;
    r[i2] *= 1+ampl/2.;
  }

  
  x[4]=x[0];
  y[4]=y[0]; // closing the polycell
  z[4]=z[0];
  r[4]=r[0]; // closing the polycell
  eta[4]=eta[0];
  phi[4]=phi[0]; // closing the polycell
  

  int color = TColor::GetColor(220, 220, 255);
  if(rh.isSeed() == 1) {
    color = TColor::GetColor(100, 150, 255);
  }

  
  switch( viewType ) {
  case  XY:
    {
      TPolyLine lineSizeXY;
      TPolyLine linePropXY;          
      if(layer == PFLayer::ECAL_BARREL || 
	 layer == PFLayer::HCAL_BARREL1 || 
	 layer == PFLayer::HCAL_BARREL2) {
	lineSizeXY.SetLineColor(color);
	//cout << "x,y " << x[0] << " " << y[0] << endl;
	lineSizeXY.SetFillColor(color);
	lineSizeXY.DrawPolyLine(5,x,y,"f");
      } else {
	//cout << "x,y " << x[0] << " " << y[0] << endl;
	lineSizeXY.SetLineColor(color);
	lineSizeXY.DrawPolyLine(5,x,y);
	
	xprop[4]=xprop[0];
	yprop[4]=yprop[0]; // closing the polycell    
	linePropXY.SetLineColor(color);
	linePropXY.SetFillColor(color);
	linePropXY.DrawPolyLine(5,xprop,yprop,"F");
      }
      break;
    }
  case RZ:
    {
      TPolyLine lineSizeRZ;
      lineSizeRZ.SetLineColor(color);
      lineSizeRZ.SetFillColor(color);
      // cout << "z,r " << z[0] << " " << r[0] << endl;
      lineSizeRZ.DrawPolyLine(5,z,r,"f");
      break;
    }
  case EPE:
  case EPH:
    {
      TPolyLine lineSizeEP;
      TPolyLine linePropEP;          
      
      lineSizeEP.SetLineColor(color);
      lineSizeEP.SetFillColor(color);
      lineSizeEP.DrawPolyLine(5,eta,phi);
      
      etaprop[4]=etaprop[0];
      phiprop[4]=phiprop[0]; // closing the polycell    
      linePropEP.SetLineColor(color);
      linePropEP.SetFillColor(color);
      linePropEP.DrawPolyLine(5,etaprop,phiprop,"F");

      break;
    }
  default:
    break;
    
  }
}

void PFRootEventManager::displayClusters(unsigned viewType, double phi0) {
  
  std::vector<reco::PFCluster>::iterator itCluster;
  for(itCluster = clusters_.begin(); itCluster != clusters_.end(); 
      itCluster++) {
    
    int type = itCluster->getType();
    if(algosToDisplay_.find(type) == algosToDisplay_.end() )
      continue;
    
    TMarker m;

    int color = 4;
    if( displayColorClusters_ ) 
      color = itCluster->getType();

    m.SetMarkerColor(color);
    m.SetMarkerStyle(20);
    
    math::XYZPoint xyzPos = itCluster->getPositionXYZ();

    switch(viewType) {
    case XY:
      m.DrawMarker(xyzPos.X(), xyzPos.Y());
      break;
    case RZ:
      {
	double sign = 1.;
	if (cos(phi0 - itCluster->getPositionXYZ().Phi()) < 0.)
	  sign = -1.;
	m.DrawMarker(xyzPos.Z(), sign*xyzPos.Rho());
	break;
      }
    case EPE:
      if(itCluster->getLayer()<0)
	m.DrawMarker(xyzPos.Eta(), xyzPos.Phi());
      break;
    case EPH:
      if(itCluster->getLayer()>0)
	m.DrawMarker(xyzPos.Eta(), xyzPos.Phi());
      break;
    }      
  }
}

void PFRootEventManager::displayRecTracks(unsigned viewType, double phi0) 
{
  // math::XYZPoint vPhi0(cos(phi0), sin(phi0), 0.);

  std::vector<reco::PFRecTrack>::iterator itRecTrack;
  for (itRecTrack = recTracks_.begin(); itRecTrack != recTracks_.end();
       itRecTrack++) {

    double sign = 1.;

    const reco::PFTrajectoryPoint& tpatecal 
      = itRecTrack->trajectoryPoint(itRecTrack->nTrajectoryMeasurements() +
				    reco::PFTrajectoryPoint::ECALEntrance );
    
    if ( cos(phi0 - tpatecal.momentum().Phi()) < 0.)
      sign = -1.;

    // Check number of measurements with non-zero momentum
    // COLIN: this was a copy ! should pass a reference
    //     std::vector<reco::PFTrajectoryPoint> trajectoryPoints = 
    //       itRecTrack->getTrajectoryPoints();
    
    // COLIN: iterator uninitialized 
    //     std::vector<reco::PFTrajectoryPoint>::iterator itTrajPt;
    //     unsigned nValidPts = 0;
    
    // COLIN: unnecessary loop
    //     for (itTrajPt = trajectoryPoints.begin(); 
    // 	 itTrajPt != trajectoryPoints.end(); itTrajPt++)
    //       if (itTrajPt->getMomentum().P() > 0.) nValidPts++;
    //     if (!nValidPts) continue;


    // reserving space. nb not all trajectory points are valid
    
    const std::vector<reco::PFTrajectoryPoint>& trajectoryPoints = 
      itRecTrack->trajectoryPoints();

    vector<double> xPos;
    xPos.reserve( itRecTrack->nTrajectoryPoints() );
    vector<double> yPos;
    yPos.reserve( itRecTrack->nTrajectoryPoints() );
    
    // COLIN: avoid double* 
//     double* xPos = new double[itTrajPt->getNTrajectoryPoints()];
//     double* yPos = new double[itTrajPt->getNTrajectoryPoints()];
//     unsigned iValid = 0;

    // cout << "Draw a new track " << nValidPts << endl;

    typedef vector<reco::PFTrajectoryPoint>::const_iterator IPT;
    
    for (IPT itTrajPt = trajectoryPoints.begin(); 
	 itTrajPt != trajectoryPoints.end(); itTrajPt++ ) {

      if (itTrajPt->momentum().P() > 0.) {

	// COLIN: this is bugged
	// math::XYZPoint xyzPos(itTrajPt->xyzPosition().X()*vPhi0.Y() - itTrajPt->xyzPosition().Y()*vPhi0.X(), itTrajPt->xyzPosition().X()*vPhi0.X() + itTrajPt->xyzPosition().Y()*vPhi0.Y(), itTrajPt->xyzPosition().Z());
	// xyzPos.SetPhi(xyzPos.Phi()); // <=== Does not work ??? why ???
	
	math::XYZPoint xyzPos = itTrajPt->positionXYZ();
	
	switch(viewType) {
	case XY:
	  xPos.push_back(xyzPos.X());
	  yPos.push_back(xyzPos.Y());
	  // cout << "\t" << itTrajPt->xyzPosition().X() << " " 
	  //     << itTrajPt->xyzPosition().Y() << endl;
	  break;
	case RZ:
	  xPos.push_back(xyzPos.Z());
	  yPos.push_back(sign*xyzPos.Rho());
	  break;
	case EPE:
	case EPH:	 
	  // closest approach is meaningless in eta/phi
	  if( itTrajPt->layer() == reco::PFTrajectoryPoint::ClosestApproach)
	    continue;
	  
	  // 	  cout<<itTrajPt->getPositionXYZ().Eta()<<" "
	  // 	      <<itTrajPt->getPositionXYZ().Phi()<<" "
	  // 	      <<itTrajPt->getPositionXYZ().X()<<" "
	  // 	      <<itTrajPt->getPositionXYZ().Y()<<endl;
	  // 	  cout<<xyzPos.Eta()<<" "<<xyzPos.Phi()<<" "<<xyzPos.X()<<" "<<xyzPos.Y()<<endl;
	  
	  xPos.push_back(xyzPos.Eta());
	  yPos.push_back(xyzPos.Phi());
	  break;
	}
      }
    }  

    graphTrack_[viewType].push_back(new TGraph(xPos.size(), &xPos[0], &yPos[0]));
    int color = 103;
 
    unsigned lastHisto = graphTrack_[viewType].size() - 1;

    graphTrack_[viewType][lastHisto]->SetMarkerColor(color);
    graphTrack_[viewType][lastHisto]->SetMarkerStyle(8);
    graphTrack_[viewType][lastHisto]->SetMarkerSize(0.5);
    graphTrack_[viewType][lastHisto]->SetLineColor(color);
    graphTrack_[viewType][lastHisto]->SetLineStyle(itRecTrack->algoType());
    graphTrack_[viewType][lastHisto]->Draw("pl");

  }
  return;
}


void PFRootEventManager::unZoom() {
  if(displayHistEtaPhi_) {
    displayHistEtaPhi_->GetXaxis()->UnZoom();
    displayHistEtaPhi_->GetYaxis()->UnZoom();
  }
  
  for( unsigned i=0; i<displayHist_.size(); i++) {
    assert( displayHist_[i] );
    displayHist_[i]->GetXaxis()->UnZoom();
    displayHist_[i]->GetYaxis()->UnZoom();
  }

  updateDisplay();
}



void PFRootEventManager::updateDisplay() {

  typedef map<int, TCanvas * >::iterator IT;
  for(IT it = displayEtaPhi_.begin(); it!=displayEtaPhi_.end(); it++) {
    if( gROOT->GetListOfCanvases()->FindObject(it->second) ) {
      it->second->cd(1);
      gPad->Modified();
      it->second->cd(2);
      gPad->Modified();
    }
  }

  for( unsigned i=0; i<displayView_.size(); i++) {
    if( gROOT->GetListOfCanvases()->FindObject(displayView_[i]) )
    displayView_[i]->Modified();
  }
}


void PFRootEventManager::lookForMaxRecHit(bool ecal) {

  // look for the rechit with max e in ecal or hcal
  double maxe = -999;
  reco::PFRecHit* maxrh = 0;

  for(unsigned i=0; i<rechits_.size(); i++) {

    if(ecal &&
       rechits_[i].getLayer() != PFLayer::ECAL_BARREL && 
       rechits_[i].getLayer() != PFLayer::ECAL_ENDCAP )
      continue;

    if(!ecal &&
       rechits_[i].getLayer() != PFLayer::HCAL_BARREL1 && 
       rechits_[i].getLayer() != PFLayer::HCAL_ENDCAP )
      continue;

    double energy = rechits_[i].getEnergy();

    if(energy > maxe ) {
      maxe = energy;
      maxrh = &(rechits_[i]);
    }      
  }
  
  if(!maxrh) return;

  // center view on this rechit


  // get the cell size to set the eta and phi width 
  // of the display window from one of the cells
  
  double phisize = -1;
  double etasize = -1;
  maxrh->size(phisize, etasize);
   
  double etagate = displayZoomFactor_ * etasize;
  double phigate = displayZoomFactor_ * phisize;
  
  double eta =  maxrh->getPositionREP().Eta();
  double phi =  maxrh->getPositionREP().Phi();
  

  if(displayHist_[EPE]) {
    displayHist_[EPE]->GetXaxis()->SetRangeUser(eta-etagate, eta+etagate);
    displayHist_[EPE]->GetYaxis()->SetRangeUser(phi-phigate, phi+phigate);
  }
  if(displayHist_[EPH]) {
    displayHist_[EPH]->GetXaxis()->SetRangeUser(eta-etagate, eta+etagate);
    displayHist_[EPH]->GetYaxis()->SetRangeUser(phi-phigate, phi+phigate);
  }
  
  updateDisplay();
}  





double PFRootEventManager::getMaxEEcal() {
  
  if( maxERecHitEcal_<0 ) {
    double maxeec = getMaxE( PFLayer::ECAL_ENDCAP );
    double maxeb =  getMaxE( PFLayer::ECAL_BARREL );
    maxERecHitEcal_ = maxeec > maxeb ? maxeec:maxeb; 
    // max of both barrel and endcap
  }
  return  maxERecHitEcal_;
}




double PFRootEventManager::getMaxEHcal() {

  if(maxERecHitHcal_ < 0) {
    double maxeec = getMaxE( PFLayer::HCAL_ENDCAP );
    double maxeb =  getMaxE( PFLayer::HCAL_BARREL1 );
    maxERecHitHcal_ =  maxeec>maxeb  ?  maxeec:maxeb;
  }
  return maxERecHitHcal_;
}


double PFRootEventManager::getMaxE(int layer) const {

  double maxe = -9999;

  for( unsigned i=0; i<rechits_.size(); i++) {
    if( rechits_[i].getLayer() != layer ) continue;
    if( rechits_[i].getEnergy() > maxe)
      maxe = rechits_[i].getEnergy();
  }

  return maxe;
}

