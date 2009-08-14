/**_________________________________________________________________
   class:   BeamSpotAnalyzer.cc
   package: RecoVertex/BeamSpotProducer
   


 author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
         Geng-Yuan Jeng, UC Riverside (Geng-Yuan.Jeng@cern.ch)

 version $Id: BeamSpotAnalyzer.cc,v 1.7 2009/03/26 20:04:12 yumiceva Exp $

________________________________________________________________**/


// C++ standard
#include <string>
// CMS
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoVertex/BeamSpotProducer/interface/BeamSpotAnalyzer.h"
#include "RecoVertex/BeamSpotProducer/interface/BSFitter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/BeamSpotObjects/interface/BeamSpotObjects.h"

#include "TMath.h"

BeamSpotAnalyzer::BeamSpotAnalyzer(const edm::ParameterSet& iConfig)
{

  outputfilename_ = iConfig.getUntrackedParameter<std::string>("OutputFileName");
  
//   file_ = TFile::Open(outputfilename_.c_str(),"RECREATE");

//   ftree_ = new TTree("mytree","mytree");
//   ftree_->AutoSave();
  
//   ftree_->Branch("pt",&fpt,"fpt/D");
//   ftree_->Branch("d0",&fd0,"fd0/D");
//   ftree_->Branch("sigmad0",&fsigmad0,"fsigmad0/D");
//   ftree_->Branch("phi0",&fphi0,"fphi0/D");
//   ftree_->Branch("z0",&fz0,"fz0/D");
//   ftree_->Branch("sigmaz0",&fsigmaz0,"fsigmaz0/D");
//   ftree_->Branch("theta",&ftheta,"ftheta/D");
//   ftree_->Branch("eta",&feta,"feta/D");
//   ftree_->Branch("charge",&fcharge,"fcharge/I");
//   ftree_->Branch("chi2",&fchi2,"fchi2/D");
//   ftree_->Branch("ndof",&fndof,"fndof/D");
//   ftree_->Branch("nHit",&fnHit,"fnHit/i");
//   ftree_->Branch("nStripHit",&fnStripHit,"fnStripHit/i");
//   ftree_->Branch("nPixelHit",&fnPixelHit,"fnPixelHit/i");
//   ftree_->Branch("nTIBHit",&fnTIBHit,"fnTIBHit/i");
//   ftree_->Branch("nTOBHit",&fnTOBHit,"fnTOBHit/i");
//   ftree_->Branch("nTIDHit",&fnTIDHit,"fnTIDHit/i");
//   ftree_->Branch("nTECHit",&fnTECHit,"fnTECHit/i");
//   ftree_->Branch("nPXBHit",&fnPXBHit,"fnPXBHit/i");
//   ftree_->Branch("nPXFHit",&fnPXFHit,"fnPXFHit/i");
//   ftree_->Branch("cov",&fcov,"fcov[7][7]/D");
  
   
//   fBSvector.clear();

  //dump to file
  fasciiFileName = outputfilename_.replace(outputfilename_.size()-4,outputfilename_.size(),"txt");
  fasciiFile.open(fasciiFileName.c_str());

  // get parameter
  write2DB_ = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getParameter<bool>("WriteToDB");
  runallfitters_ = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getParameter<bool>("RunAllFitters");
  inputBeamWidth_ = iConfig.getParameter<edm::ParameterSet>("BSAnalyzerParameters").getUntrackedParameter<double>("InputBeamWidth",-1.);
  
  theBeamFitter = new BeamFitter(iConfig);
  theBeamFitter->resetTrkVector();

  ftotal_tracks = 0;
  ftotalevents = 0;
  
}


BeamSpotAnalyzer::~BeamSpotAnalyzer()
{
  delete theBeamFitter;
}


void
BeamSpotAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	
// 	ftree_->SetBranchAddress("theta",&ftheta);
// 	ftree_->SetBranchAddress("pt",&fpt);
// 	ftree_->SetBranchAddress("eta",&feta);
// 	ftree_->SetBranchAddress("charge",&fcharge);
// 	ftree_->SetBranchAddress("chi2",&fchi2);
// 	ftree_->SetBranchAddress("ndof",&fndof);
//      ftree_->SetBranchAddress("d0",&fd0);
// 	ftree_->SetBranchAddress("sigmad0",&fsigmad0);
// 	ftree_->SetBranchAddress("phi0",&fphi0);
// 	ftree_->SetBranchAddress("z0",&fz0);
// 	ftree_->SetBranchAddress("sigmaz0",&fsigmaz0);
// 	ftree_->SetBranchAddress("nHit",&fnHit);
// 	ftree_->SetBranchAddress("nStripHit",&fnStripHit);
// 	ftree_->SetBranchAddress("nPixelHit",&fnPixelHit);
// 	ftree_->SetBranchAddress("nTIBHit",&fnTIBHit);
// 	ftree_->SetBranchAddress("nTOBHit",&fnTOBHit);
// 	ftree_->SetBranchAddress("nTIDHit",&fnTIDHit);
// 	ftree_->SetBranchAddress("nTECHit",&fnTECHit);
// 	ftree_->SetBranchAddress("nPXBHit",&fnPXBHit);
// 	ftree_->SetBranchAddress("nPXFHit",&fnPXFHit);
// 	ftree_->SetBranchAddress("cov",&fcov);

	
// 	  ftree_->Fill();

  theBeamFitter->readEvent(iEvent);
  ftotalevents++;

}



void 
BeamSpotAnalyzer::beginJob(const edm::EventSetup&)
{
}

void 
BeamSpotAnalyzer::endJob() {
  std::cout << "\n-------------------------------------\n" << std::endl;
  std::cout << "\n Total number of events processed: "<< ftotalevents << std::endl;
  std::cout << "\n-------------------------------------\n\n" << std::endl;

  if(theBeamFitter->runFitter()){
    reco::BeamSpot beam_default = theBeamFitter->getBeamSpot();
    
    std::cout << "\n RESULTS OF DEFAULT FIT:" << std::endl;
    std::cout << beam_default << std::endl;
    
    // dump to file
    fasciiFile << "X " << beam_default.x0() << std::endl;
    fasciiFile << "Y " << beam_default.y0() << std::endl;
    fasciiFile << "Z " << beam_default.z0() << std::endl;
    fasciiFile << "sigmaZ " << beam_default.sigmaZ() << std::endl;
    fasciiFile << "dxdz " << beam_default.dxdz() << std::endl;
    fasciiFile << "dydz " << beam_default.dydz() << std::endl;
    if (inputBeamWidth_ > 0 ) {
      fasciiFile << "BeamWidthX " << inputBeamWidth_ << std::endl;
      fasciiFile << "BeamWidthY " << inputBeamWidth_ << std::endl;
    } else {
      fasciiFile << "BeamWidthX " << beam_default.BeamWidthX() << std::endl;
      fasciiFile << "BeamWidthY " << beam_default.BeamWidthY() << std::endl;
    }
	
    for (int i = 0; i<6; ++i) {
      fasciiFile << "Cov("<<i<<",j) ";
      for (int j=0; j<7; ++j) {
	fasciiFile << beam_default.covariance(i,j) << " ";
      }
      fasciiFile << std::endl;
    }
    // beam width error
    if (inputBeamWidth_ > 0 ) {
      fasciiFile << "Cov(6,j) 0 0 0 0 0 0 " << pow(2.e-4,2) << std::endl;
    } else {
      fasciiFile << "Cov(6,j) 0 0 0 0 0 0 " << beam_default.covariance(6,6) << std::endl;
    }
    fasciiFile << "EmittanceX " << beam_default.emittanceX() << std::endl;
    fasciiFile << "EmittanceY " << beam_default.emittanceY() << std::endl;
    fasciiFile << "BetaStar " << beam_default.betaStar() << std::endl;
	
	
    if (write2DB_) {
      std::cout << "\n-------------------------------------\n\n" << std::endl;
      std::cout << " write results to DB..." << std::endl;
      
      BeamSpotObjects *pBSObjects = new BeamSpotObjects();

      //pBSObjects->Put(beam_default);
      pBSObjects->SetPosition(beam_default.position().X(),beam_default.position().Y(),beam_default.position().Z());
      //std::cout << " wrote: x= " << beam_default.position().X() << " y= "<< beam_default.position().Y() << " z= " << beam_default.position().Z() << std::endl;
      pBSObjects->SetSigmaZ(beam_default.sigmaZ());
      pBSObjects->Setdxdz(beam_default.dxdz());
      pBSObjects->Setdydz(beam_default.dydz());
      if (inputBeamWidth_ > 0 ) {
	std::cout << " beam width value forced to be " << inputBeamWidth_ << std::endl;
	pBSObjects->SetBeamWidthX(inputBeamWidth_);
	pBSObjects->SetBeamWidthY(inputBeamWidth_);
      } else {
	// need to fix this
	std::cout << " using default value, 15e-4, for beam width!!!"<<std::endl;
	pBSObjects->SetBeamWidthX(15.0e-4);
	pBSObjects->SetBeamWidthY(15.0e-4);

      }
		
      for (int i = 0; i<7; ++i) {
	for (int j=0; j<7; ++j) {
	  pBSObjects->SetCovariance(i,j,beam_default.covariance(i,j));
	}
      }
      edm::Service<cond::service::PoolDBOutputService> poolDbService;
      if( poolDbService.isAvailable() ) {
	std::cout << "poolDBService available"<<std::endl;
	if ( poolDbService->isNewTagRequest( "BeamSpotObjectsRcd" ) ) {
	  std::cout << "new tag requested" << std::endl;
	  poolDbService->createNewIOV<BeamSpotObjects>( pBSObjects, poolDbService->beginOfTime(),poolDbService->endOfTime(),
							"BeamSpotObjectsRcd"  );
	}
	else {
	  std::cout << "no new tag requested" << std::endl;
	  poolDbService->appendSinceTime<BeamSpotObjects>( pBSObjects, poolDbService->currentTime(),
							   "BeamSpotObjectsRcd" );
	}

      }
    }

    if (runallfitters_) {
      theBeamFitter->runAllFitter();
      
// 	// add new branches
// 	std::cout << " add new branches to output file " << std::endl;
// 	beam_default = myalgo->Fit_d0phi();
// 	file_->cd();
// 	TTree *newtree = new TTree("mytreecorr","mytreecorr");
// 	newtree->Branch("d0phi_chi2",&fd0phi_chi2,"fd0phi_chi2/D");
// 	newtree->Branch("d0phi_d0",&fd0phi_d0,"fd0phi_d0/D");
// 	newtree->SetBranchAddress("d0phi_chi2",&fd0phi_chi2);
// 	newtree->SetBranchAddress("d0phi_d0",&fd0phi_d0);
// 	std::vector<BSTrkParameters>  tmpvector = myalgo->GetData();
	
// 	std::vector<BSTrkParameters>::iterator iparam = tmpvector.begin();
// 	for( iparam = tmpvector.begin() ;
// 		 iparam != tmpvector.end() ; ++iparam) {
// 		fd0phi_chi2 = iparam->d0phi_chi2();
// 		fd0phi_d0   = iparam->d0phi_d0();
// 		newtree->Fill();
// 	}
// 	newtree->Write();
    }

	// let's close everything
// 	file_->cd();
// 	ftree_->Write();
// 	file_->Close();

  }    		
  std::cout << "[BeamSpotAnalyzer] endJob done \n" << std::endl;
}

//define this as a plug-in
DEFINE_ANOTHER_FWK_MODULE(BeamSpotAnalyzer);
