/**_________________________________________________________________
   class:   BeamFitter.cc
   package: RecoVertex/BeamSpotProducer
   


   author: Francisco Yumiceva, Fermilab (yumiceva@fnal.gov)
           Geng-Yuan Jeng, UC Riverside (Geng-Yuan.Jeng@cern.ch)
 
   version $Id: BeamFitter.cc,v 1.49 2010/04/03 10:17:13 jengbou Exp $

________________________________________________________________**/

#include "RecoVertex/BeamSpotProducer/interface/BeamFitter.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/BeamSpotObjects/interface/BeamSpotObjects.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// ----------------------------------------------------------------------
// Useful function:
// ----------------------------------------------------------------------
static char * formatTime(const std::time_t & t)  {
  struct std::tm * ptm;
  ptm = gmtime(&t);
  static char ts[32];
  strftime(ts,sizeof(ts),"%Y.%m.%d %H:%M:%S %Z",ptm);
  return ts;
}

BeamFitter::BeamFitter(const edm::ParameterSet& iConfig)
{

  debug_             = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<bool>("Debug");
  tracksLabel_       = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<edm::InputTag>("TrackCollection");
  writeTxt_          = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<bool>("WriteAscii");
  outputTxt_         = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<std::string>("AsciiFileName");
  appendRunTxt_      = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<bool>("AppendRunToFileName");
  writeDIPTxt_       = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<bool>("WriteDIPAscii");
  outputDIPTxt_      = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<std::string>("DIPFileName");
  saveNtuple_        = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<bool>("SaveNtuple");
  saveBeamFit_       = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<bool>("SaveFitResults");
  isMuon_            = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<bool>("IsMuonCollection");

  trk_MinpT_         = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<double>("MinimumPt");
  trk_MaxEta_        = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<double>("MaximumEta");
  trk_MaxIP_         = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<double>("MaximumImpactParameter");
  trk_MaxZ_          = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<double>("MaximumZ");
  trk_MinNTotLayers_ = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<int>("MinimumTotalLayers");
  trk_MinNPixLayers_ = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<int>("MinimumPixelLayers");
  trk_MaxNormChi2_   = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<double>("MaximumNormChi2");
  trk_Algorithm_     = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<std::vector<std::string> >("TrackAlgorithm");
  trk_Quality_       = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<std::vector<std::string> >("TrackQuality");
  min_Ntrks_         = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<int>("MinimumInputTracks");
  convergence_       = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<double>("FractionOfFittedTrks");
  inputBeamWidth_    = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<double>("InputBeamWidth",-1.);
  
  for (unsigned int j=0;j<trk_Algorithm_.size();j++)
    algorithm_.push_back(reco::TrackBase::algoByName(trk_Algorithm_[j]));
  for (unsigned int j=0;j<trk_Quality_.size();j++)
    quality_.push_back(reco::TrackBase::qualityByName(trk_Quality_[j]));
  
  if (saveNtuple_ || saveBeamFit_) {
    outputfilename_ = iConfig.getParameter<edm::ParameterSet>("BeamFitter").getUntrackedParameter<std::string>("OutputFileName");
    file_ = TFile::Open(outputfilename_.c_str(),"RECREATE");
  }
  if (saveNtuple_) {
    ftree_ = new TTree("mytree","mytree");
    ftree_->AutoSave();
    
    ftree_->Branch("pt",&fpt,"fpt/D");
    ftree_->Branch("d0",&fd0,"fd0/D");
    ftree_->Branch("d0bs",&fd0bs,"fd0bs/D");
    ftree_->Branch("sigmad0",&fsigmad0,"fsigmad0/D");
    ftree_->Branch("phi0",&fphi0,"fphi0/D");
    ftree_->Branch("z0",&fz0,"fz0/D");
    ftree_->Branch("sigmaz0",&fsigmaz0,"fsigmaz0/D");
    ftree_->Branch("theta",&ftheta,"ftheta/D");
    ftree_->Branch("eta",&feta,"feta/D");
    ftree_->Branch("charge",&fcharge,"fcharge/I");
    ftree_->Branch("normchi2",&fnormchi2,"fnormchi2/D");
    ftree_->Branch("nTotLayerMeas",&fnTotLayerMeas,"fnTotLayerMeas/i");
    ftree_->Branch("nStripLayerMeas",&fnStripLayerMeas,"fnStripLayerMeas/i");
    ftree_->Branch("nPixelLayerMeas",&fnPixelLayerMeas,"fnPixelLayerMeas/i");
    ftree_->Branch("nTIBLayerMeas",&fnTIBLayerMeas,"fnTIBLayerMeas/i");
    ftree_->Branch("nTOBLayerMeas",&fnTOBLayerMeas,"fnTOBLayerMeas/i");
    ftree_->Branch("nTIDLayerMeas",&fnTIDLayerMeas,"fnTIDLayerMeas/i");
    ftree_->Branch("nTECLayerMeas",&fnTECLayerMeas,"fnTECLayerMeas/i");
    ftree_->Branch("nPXBLayerMeas",&fnPXBLayerMeas,"fnPXBLayerMeas/i");
    ftree_->Branch("nPXFLayerMeas",&fnPXFLayerMeas,"fnPXFLayerMeas/i");
    ftree_->Branch("cov",&fcov,"fcov[7][7]/D");
    ftree_->Branch("vx",&fvx,"fvx/D");
    ftree_->Branch("vy",&fvy,"fvy/D");
    ftree_->Branch("quality",&fquality,"fquality/O");
    ftree_->Branch("algo",&falgo,"falgo/O");
    ftree_->Branch("run",&frun,"frun/i");
    ftree_->Branch("lumi",&flumi,"flumi/i");
    ftree_->Branch("pvValid",&fpvValid,"fpvValid/O");
    ftree_->Branch("pvx", &fpvx, "fpvx/D");
    ftree_->Branch("pvy", &fpvy, "fpvy/D");
    ftree_->Branch("pvz", &fpvz, "fpvz/D");
  }
  if (saveBeamFit_){
    ftreeFit_ = new TTree("fitResults","fitResults");
    ftreeFit_->AutoSave();
    ftreeFit_->Branch("run",&frunFit,"frunFit/i");
    ftreeFit_->Branch("beginLumi",&fbeginLumiOfFit,"fbeginLumiOfFit/i");
    ftreeFit_->Branch("endLumi",&fendLumiOfFit,"fendLumiOfFit/i");
    ftreeFit_->Branch("beginTime",fbeginTimeOfFit,"fbeginTimeOfFit/C");
    ftreeFit_->Branch("endTime",fendTimeOfFit,"fendTimeOfFit/C");
    ftreeFit_->Branch("x",&fx,"fx/D");
    ftreeFit_->Branch("y",&fy,"fy/D");
    ftreeFit_->Branch("z",&fz,"fz/D");
    ftreeFit_->Branch("sigmaZ",&fsigmaZ,"fsigmaZ/D");
    ftreeFit_->Branch("dxdz",&fdxdz,"fdxdz/D");
    ftreeFit_->Branch("dydz",&fdydz,"fdydz/D");
    ftreeFit_->Branch("xErr",&fxErr,"fxErr/D");
    ftreeFit_->Branch("yErr",&fyErr,"fyErr/D");
    ftreeFit_->Branch("zErr",&fzErr,"fzErr/D");
    ftreeFit_->Branch("sigmaZErr",&fsigmaZErr,"fsigmaZErr/D");
    ftreeFit_->Branch("dxdzErr",&fdxdzErr,"fdxdzErr/D");
    ftreeFit_->Branch("dydzErr",&fdydzErr,"fdydzErr/D");
    ftreeFit_->Branch("widthX",&fwidthX,"fwidthX/D");
    ftreeFit_->Branch("widthY",&fwidthY,"fwidthY/D");
    ftreeFit_->Branch("widthXErr",&fwidthXErr,"fwidthXErr/D");
    ftreeFit_->Branch("widthYErr",&fwidthYErr,"fwidthYErr/D");
  }
  
  fBSvector.clear();
  ftotal_tracks = 0;
  fnTotLayerMeas = fnPixelLayerMeas = fnStripLayerMeas = fnTIBLayerMeas = 0;
  fnTIDLayerMeas = fnTOBLayerMeas = fnTECLayerMeas = fnPXBLayerMeas = fnPXFLayerMeas = 0;
  frun = flumi = -1;
  frunFit = fbeginLumiOfFit = fendLumiOfFit = -1;
  fquality = falgo = true;
  fpvValid = true;
  fpvx = fpvy = fpvz = 0;
  fitted_ = false;
  resetRefTime();
  
  //debug histograms
  h1ntrks = new TH1F("h1ntrks","number of tracks per event",50,0,50);
  h1vz_event = new TH1F("h1vz_event","track Vz", 50, -30, 30 );
  h1cutFlow = new TH1F("h1cutFlow","Cut flow table of track selection", 9, 0, 9);
  h1cutFlow->GetXaxis()->SetBinLabel(1,"No cut");
  h1cutFlow->GetXaxis()->SetBinLabel(2,"Traker hits");
  h1cutFlow->GetXaxis()->SetBinLabel(3,"Pixel hits");
  h1cutFlow->GetXaxis()->SetBinLabel(4,"norm. #chi^{2}");
  h1cutFlow->GetXaxis()->SetBinLabel(5,"algo");
  h1cutFlow->GetXaxis()->SetBinLabel(6,"quality");
  h1cutFlow->GetXaxis()->SetBinLabel(7,"d_{0}");
  h1cutFlow->GetXaxis()->SetBinLabel(8,"z_{0}");
  h1cutFlow->GetXaxis()->SetBinLabel(9,"p_{T}");
  resetCutFlow();

  // Primary vertex fitter
  MyPVFitter = new PVFitter(iConfig);
  MyPVFitter->resetAll();
}

BeamFitter::~BeamFitter() {
  if (saveNtuple_) {
    file_->cd();
    if (fitted_ && h1z) h1z->Write();
    h1ntrks->Write();
    h1vz_event->Write();
    if (h1cutFlow) h1cutFlow->Write();
    ftree_->Write();
  }
  if (saveBeamFit_){
    file_->cd();
    ftreeFit_->Write();
  }
  if (saveNtuple_ || saveBeamFit_)
    file_->Close();
}


void BeamFitter::readEvent(const edm::Event& iEvent)
{

  //open txt files at begin of run
  if (frun == -1) {
    if (writeDIPTxt_)
      fasciiDIP.open(outputDIPTxt_.c_str());
    if (writeTxt_) {
      std::string tmpname = outputTxt_;
      char index[15];
      if (appendRunTxt_) {
          sprintf(index,"%s%i","_Run", int(iEvent.id().run()));
          tmpname.insert(outputTxt_.length()-4,index);
      }
      fasciiFile.open(tmpname.c_str());
      outputTxt_ = tmpname;
    }
  }
  frun = iEvent.id().run();
  const edm::TimeValue_t ftimestamp = iEvent.time().value();
  const std::time_t ftmptime = ftimestamp >> 32;

  if (fbeginLumiOfFit == -1) freftime[0] = freftime[1] = ftmptime;
  if (freftime[0] == 0 || ftmptime < freftime[0]) freftime[0] = ftmptime;
  if (freftime[1] == 0 || ftmptime > freftime[1]) freftime[1] = ftmptime;
  flumi = iEvent.luminosityBlock();
  frunFit = frun;
  if (fbeginLumiOfFit == -1 || fbeginLumiOfFit > flumi) fbeginLumiOfFit = flumi;
  if (fendLumiOfFit == -1 || fendLumiOfFit < flumi) fendLumiOfFit = flumi;

  edm::Handle<reco::TrackCollection> TrackCollection;
  iEvent.getByLabel(tracksLabel_, TrackCollection);

  //------ Primary Vertices
  edm::Handle< edm::View<reco::Vertex> > PVCollection;
  bool hasPVs = false;
  edm::View<reco::Vertex> pv;
  if ( iEvent.getByLabel("offlinePrimaryVertices", PVCollection ) ) {
      pv = *PVCollection;
      hasPVs = true;
  }
  //------

  //------ Beam Spot in current event
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  const reco::BeamSpot *refBS =  0;
  if ( iEvent.getByLabel("offlineBeamSpot",recoBeamSpotHandle) )
      refBS = recoBeamSpotHandle.product();
  //-------
  
  const reco::TrackCollection *tracks = TrackCollection.product();
  
  double eventZ = 0;
  double averageZ = 0;

  for ( reco::TrackCollection::const_iterator track = tracks->begin();
	track != tracks->end();
	++track ) {
    
    if ( ! isMuon_) {
      
      const reco::HitPattern& trkHP = track->hitPattern();
      
      fnPixelLayerMeas = trkHP.pixelLayersWithMeasurement();
      fnStripLayerMeas = trkHP.stripLayersWithMeasurement();
      fnTotLayerMeas = trkHP.trackerLayersWithMeasurement();
      fnPXBLayerMeas = trkHP.pixelBarrelLayersWithMeasurement();
      fnPXFLayerMeas = trkHP.pixelEndcapLayersWithMeasurement();
      fnTIBLayerMeas = trkHP.stripTIBLayersWithMeasurement();
      fnTIDLayerMeas = trkHP.stripTIDLayersWithMeasurement();
      fnTOBLayerMeas = trkHP.stripTOBLayersWithMeasurement();
      fnTECLayerMeas = trkHP.stripTECLayersWithMeasurement();
    } else {
      
      fnTotLayerMeas = track->numberOfValidHits();
      
    }
    
    fpt = track->pt();
    feta = track->eta();
    fphi0 = track->phi();
    fcharge = track->charge();
    fnormchi2 = track->normalizedChi2();
    fd0 = track->d0();
    if (refBS) fd0bs = -1*track->dxy(refBS->position());
    else fd0bs = 0.;
    
    fsigmad0 = track->d0Error();
    fz0 = track->dz();
    fsigmaz0 = track->dzError();
    ftheta = track->theta();
    fvx = track->vx();
    fvy = track->vy();

    for (int i=0; i<5; ++i) {
      for (int j=0; j<5; ++j) {
	fcov[i][j] = track->covariance(i,j);
      }
    }
    
    fquality = true;
    falgo = true;
    
    if (! isMuon_ ) {
      if (quality_.size()!=0) {
          fquality = false;
          for (unsigned int i = 0; i<quality_.size();++i) {
              if(debug_) std::cout << "quality_[" << i << "] = " << track->qualityName(quality_[i]) << std::endl;
              if (track->quality(quality_[i])) {
                  fquality = true;
                  break;
              }
          }
      }
      
      
      // Track algorithm
      
      if (algorithm_.size()!=0) {
	if (std::find(algorithm_.begin(),algorithm_.end(),track->algo())==algorithm_.end())
	  falgo = false;
      }
      
    }
    
    // check if we have a valid PV
    fpvValid = false;

    if ( hasPVs ) {
    
        for ( size_t ipv=0; ipv != pv.size(); ++ipv ) {

            if (! pv[ipv].isFake() ) fpvValid = true;
      
            if ( ipv==0 && !pv[0].isFake() ) { fpvx = pv[0].x(); fpvy = pv[0].y(); fpvz = pv[0].z(); } // fix this later

            
        }
    
    }
    
    
    if (saveNtuple_) ftree_->Fill();
    ftotal_tracks++;

    countPass[0] = ftotal_tracks;
    // Track selection
    if (fnTotLayerMeas >= trk_MinNTotLayers_) { countPass[1] += 1;
      if (fnPixelLayerMeas >= trk_MinNPixLayers_) { countPass[2] += 1;
	if (fnormchi2 < trk_MaxNormChi2_) { countPass[3] += 1;
	  if (falgo) {countPass[4] += 1;
	    if (fquality) { countPass[5] += 1;
	      if (std::abs( fd0 ) < trk_MaxIP_) { countPass[6] += 1;
		if (std::abs( fz0 ) < trk_MaxZ_){ countPass[7] += 1;		  
		  if (fpt > trk_MinpT_) {
		    countPass[8] += 1;
		    if (std::abs( feta ) < trk_MaxEta_
			//&& fpvValid
			) {
		      if (debug_) {
			std::cout << "Selected track quality = " << track->qualityMask();
			std::cout << "; track algorithm = " << track->algoName() << "= TrackAlgorithm: " << track->algo() << std::endl;
		      }
		      BSTrkParameters BSTrk(fz0,fsigmaz0,fd0,fsigmad0,fphi0,fpt,0.,0.);
		      BSTrk.setVx(fvx);
		      BSTrk.setVy(fvy);
		      fBSvector.push_back(BSTrk);
		      averageZ += fz0;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }// track selection

  }// tracks

  averageZ = averageZ/(float)(fBSvector.size());

  for( std::vector<BSTrkParameters>::const_iterator iparam = fBSvector.begin(); iparam != fBSvector.end(); ++iparam) {

    eventZ += fabs( iparam->z0() - averageZ );
    
  }
  
  h1ntrks->Fill( fBSvector.size() );
  h1vz_event->Fill( eventZ/(float)(fBSvector.size() ) ) ;
  for (unsigned int i=0; i < sizeof(countPass)/sizeof(countPass[0]); i++) 
    h1cutFlow->SetBinContent(i+1,countPass[i]);

  MyPVFitter->readEvent(iEvent);

}

bool BeamFitter::runFitter() {
  const char* fbeginTime = formatTime(freftime[0]);
  sprintf(fbeginTimeOfFit,"%s",fbeginTime);
  const char* fendTime = formatTime(freftime[1]);
  sprintf(fendTimeOfFit,"%s",fendTime);
  if (freftime[0] == 0 || freftime[1] == 0 || fbeginLumiOfFit == -1 || fendLumiOfFit == -1) return false;

  bool fit_ok = false;
  // default fit to extract beam spot info
  if(fBSvector.size() > 1 ){
    if(debug_){
      std::cout << "Calculating beam spot..." << std::endl;
      std::cout << "We will use " << fBSvector.size() << " good tracks out of " << ftotal_tracks << std::endl;
    }
    BSFitter *myalgo = new BSFitter( fBSvector );
    myalgo->SetMaximumZ( trk_MaxZ_ );
    myalgo->SetConvergence( convergence_ );
    myalgo->SetMinimumNTrks( min_Ntrks_ );
    if (inputBeamWidth_ > 0 ) myalgo->SetInputBeamWidth( inputBeamWidth_ );

    fbeamspot = myalgo->Fit();

    if ( fbeamspot.type() != 0 && MyPVFitter->runFitter() ) { // Run PVFitter if d0-phi fit returns non zero results

      fbeamspot.setBeamWidthX( MyPVFitter->getWidthX() );
      fbeamspot.setBeamWidthY( MyPVFitter->getWidthY() );
      // to pass errors I have to create another beam spot object
      reco::BeamSpot::CovarianceMatrix matrix;
      for (int j = 0 ; j < 7 ; ++j) {
          for(int k = j ; k < 7 ; ++k) {
              matrix(j,k) = fbeamspot.covariance(j,k);
          }
      }
      matrix(6,6) = MyPVFitter->getWidthXerr() * MyPVFitter->getWidthXerr();

      // get Z and sigmaZ from PV fit
      matrix(2,2) = MyPVFitter->getBeamSpot().covariance(2,2);
      matrix(3,3) = MyPVFitter->getBeamSpot().covariance(3,3);
      
      reco::BeamSpot tmpbs(reco::BeamSpot::Point(fbeamspot.x0(), fbeamspot.y0(),
                                                 MyPVFitter->getBeamSpot().z0() ),
                           MyPVFitter->getBeamSpot().sigmaZ() ,
			   fbeamspot.dxdz(),
			   fbeamspot.dydz(),
			   fbeamspot.BeamWidthX(),
			   matrix,
			   fbeamspot.type() );
      tmpbs.setBeamWidthY( fbeamspot.BeamWidthY() );

      fbeamspot = tmpbs;
      
    }
    else{ // PVFitter returns false
      fbeamspot.setType(reco::BeamSpot::Unknown);
      //       if (fbeamspot.type() == 2)
      // 	fbeamspot.setType(reco::BeamSpot::BadPVFit); // bad PV fit ==> FIXME: Add type to BeamSpot.h?
      //       else
      // 	fbeamspot.setType(reco::BeamSpot::Bad); // bad d0-phi and PV fits ==> FIXME: Add type to BeamSpot.h?
    }

    if(writeTxt_ && fasciiFile.is_open()) dumpTxtFile(outputTxt_,true); // all reaults
    if(writeDIPTxt_ && fasciiDIP.is_open()) dumpTxtFile(outputDIPTxt_,false); // for DQM/DIP

    // retrieve histogram for Vz
    h1z = (TH1F*) myalgo->GetVzHisto();

    delete myalgo;
    if ( fbeamspot.type() != 0 ) {// save all results except for Fake (all 0.)
      fit_ok = true;
      if (saveBeamFit_){
	fx = fbeamspot.x0();
	fy = fbeamspot.y0();
	fz = fbeamspot.z0();
	fsigmaZ = fbeamspot.sigmaZ();
	fdxdz = fbeamspot.dxdz();
	fdydz = fbeamspot.dydz();
	fwidthX = fbeamspot.BeamWidthX();
	fwidthY = fbeamspot.BeamWidthY();
	fxErr = fbeamspot.x0Error();
	fyErr = fbeamspot.y0Error();
	fzErr = fbeamspot.z0Error();
	fsigmaZErr = fbeamspot.sigmaZ0Error();
	fdxdzErr = fbeamspot.dxdzError();
	fdydzErr = fbeamspot.dydzError();
	fwidthXErr = fbeamspot.BeamWidthXError();
	fwidthYErr = fbeamspot.BeamWidthYError();
	ftreeFit_->Fill();
      }
    }
  }
  else{ // tracks <= 1
    reco::BeamSpot tmpbs;
    fbeamspot = tmpbs;
    fbeamspot.setType(reco::BeamSpot::Fake);
    if(debug_) std::cout << "Not enough good tracks selected! No beam fit!" << std::endl;
    if(writeTxt_ && fasciiFile.is_open()) dumpTxtFile(outputTxt_,true);  // all results
    if(writeDIPTxt_ && fasciiDIP.is_open()) dumpTxtFile(outputDIPTxt_,false);// for DQM/DIP
  }
  fitted_ = true;
  return fit_ok;
}

bool BeamFitter::runBeamWidthFitter() {
  bool widthfit_ok = false;
  // default fit to extract beam spot info
  if(fBSvector.size() > 1 ){
    if(debug_){
      std::cout << "Calculating beam spot positions('d0-phi0' method) and width using llh Fit"<< std::endl;
      std::cout << ""<<std::endl;
      std::cout << "We will use " << fBSvector.size() << " good tracks out of " << ftotal_tracks << std::endl;
    }
        BSFitter *myalgo = new BSFitter( fBSvector );
        myalgo->SetMaximumZ( trk_MaxZ_ );
        myalgo->SetConvergence( convergence_ );
        myalgo->SetMinimumNTrks(min_Ntrks_);
        if (inputBeamWidth_ > 0 ) myalgo->SetInputBeamWidth( inputBeamWidth_ );


   myalgo->SetFitVariable(std::string("d*z"));
   myalgo->SetFitType(std::string("likelihood"));
   fbeamWidthFit = myalgo->Fit();

   //Add to .txt file
   if(writeTxt_   ) dumpBWTxtFile(outputTxt_);
   
   delete myalgo;

 if ( fbeamspot.type() != 0 ) // not Fake
      widthfit_ok = true;
  }
  else{
    fbeamspot.setType(reco::BeamSpot::Fake);
    if(debug_) std::cout << "Not enough good tracks selected! No beam fit!" << std::endl;
  }
  return widthfit_ok;
}

void BeamFitter::dumpBWTxtFile(std::string & fileName ){
    std::ofstream outFile;
    outFile.open(fileName.c_str(),std::ios::app);
    outFile<<"-------------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
    outFile<<"Beam width(in cm) from Log-likelihood fit (Here we assume a symmetric beam(SigmaX=SigmaY)!)"<<std::endl;
    outFile<<"   "<<std::endl;
    if (inputBeamWidth_ > 0 ) {
        outFile<< "BeamWidth= " << inputBeamWidth_ << std::endl;
    } else {
        outFile << "BeamWidth =  " <<fbeamWidthFit.BeamWidthX() <<" +/- "<<fbeamWidthFit.BeamWidthXError() << std::endl;
    }
    outFile.close();
}

void BeamFitter::dumpTxtFile(std::string & fileName, bool append){
  std::ofstream outFile;
  if (!append)
    outFile.open(fileName.c_str());
  else
    outFile.open(fileName.c_str(),std::ios::app);

  outFile << "Runnumber " << frun << std::endl;
  outFile << "BeginTimeOfFit " << fbeginTimeOfFit << std::endl;
  outFile << "EndTimeOfFit " << fendTimeOfFit << std::endl;
  outFile << "LumiRange " << fbeginLumiOfFit << " - " << fendLumiOfFit << std::endl;
  outFile << "Type " << fbeamspot.type() << std::endl;
  outFile << "X0 " << fbeamspot.x0() << std::endl;
  outFile << "Y0 " << fbeamspot.y0() << std::endl;
  outFile << "Z0 " << fbeamspot.z0() << std::endl;
  outFile << "sigmaZ0 " << fbeamspot.sigmaZ() << std::endl;
  outFile << "dxdz " << fbeamspot.dxdz() << std::endl;
  outFile << "dydz " << fbeamspot.dydz() << std::endl;
  if (inputBeamWidth_ > 0 ) {
    outFile << "BeamWidthX " << inputBeamWidth_ << std::endl;
    outFile << "BeamWidthY " << inputBeamWidth_ << std::endl;
  } else {
    outFile << "BeamWidthX " << fbeamspot.BeamWidthX() << std::endl;
    outFile << "BeamWidthY " << fbeamspot.BeamWidthY() << std::endl;
  }
	
  for (int i = 0; i<6; ++i) {
    outFile << "Cov("<<i<<",j) ";
    for (int j=0; j<7; ++j) {
      outFile << fbeamspot.covariance(i,j) << " ";
    }
    outFile << std::endl;
  }
  // beam width error
  if (inputBeamWidth_ > 0 ) {
    outFile << "Cov(6,j) 0 0 0 0 0 0 " << "1e-4" << std::endl;
  } else {
    outFile << "Cov(6,j) 0 0 0 0 0 0 " << fbeamspot.covariance(6,6) << std::endl;
  }
  outFile << "EmittanceX " << fbeamspot.emittanceX() << std::endl;
  outFile << "EmittanceY " << fbeamspot.emittanceY() << std::endl;
  outFile << "BetaStar " << fbeamspot.betaStar() << std::endl;
  outFile.close();
}

void BeamFitter::write2DB(){
  BeamSpotObjects *pBSObjects = new BeamSpotObjects();
  
  pBSObjects->SetPosition(fbeamspot.position().X(),fbeamspot.position().Y(),fbeamspot.position().Z());
  //std::cout << " wrote: x= " << fbeamspot.position().X() << " y= "<< fbeamspot.position().Y() << " z= " << fbeamspot.position().Z() << std::endl;
  pBSObjects->SetSigmaZ(fbeamspot.sigmaZ());
  pBSObjects->Setdxdz(fbeamspot.dxdz());
  pBSObjects->Setdydz(fbeamspot.dydz());
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
      pBSObjects->SetCovariance(i,j,fbeamspot.covariance(i,j));
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

void BeamFitter::runAllFitter() {
  if(fBSvector.size()!=0){
    BSFitter *myalgo = new BSFitter( fBSvector );    
    fbeamspot = myalgo->Fit_d0phi();
        
    // iterative
    if(debug_)
      std::cout << " d0-phi Iterative:" << std::endl;
    BSFitter *myitealgo = new BSFitter( fBSvector );
    myitealgo->Setd0Cut_d0phi(4.0);
    reco::BeamSpot beam_ite = myitealgo->Fit_ited0phi();
    if (debug_){
      std::cout << beam_ite << std::endl;
      std::cout << "\n Now run tests of the different fits\n";
    }
    // from here are just tests
    std::string fit_type = "chi2";
    myalgo->SetFitVariable(std::string("z"));
    myalgo->SetFitType(std::string("chi2"));
    reco::BeamSpot beam_fit_z_chi2 = myalgo->Fit();
    if (debug_){
      std::cout << " z Chi2 Fit ONLY:" << std::endl;
      std::cout << beam_fit_z_chi2 << std::endl;
    }
    
    fit_type = "combined";
    myalgo->SetFitVariable("z");
    myalgo->SetFitType("combined");
    reco::BeamSpot beam_fit_z_lh = myalgo->Fit();
    if (debug_){
      std::cout << " z Log-Likelihood Fit ONLY:" << std::endl;
      std::cout << beam_fit_z_lh << std::endl;
    }
    
    myalgo->SetFitVariable("d");
    myalgo->SetFitType("d0phi");
    reco::BeamSpot beam_fit_dphi = myalgo->Fit();
    if (debug_){
      std::cout << " d0-phi0 Fit ONLY:" << std::endl;
      std::cout << beam_fit_dphi << std::endl;
    }
    /*
    myalgo->SetFitVariable(std::string("d*z"));
    myalgo->SetFitType(std::string("likelihood"));
    reco::BeamSpot beam_fit_dz_lh = myalgo->Fit();
    if (debug_){
      std::cout << " Log-Likelihood Fit:" << std::endl;
      std::cout << beam_fit_dz_lh << std::endl;
    }
    
    myalgo->SetFitVariable(std::string("d*z"));
    myalgo->SetFitType(std::string("resolution"));
    reco::BeamSpot beam_fit_dresz_lh = myalgo->Fit();
    if(debug_){
      std::cout << " IP Resolution Fit" << std::endl;
      std::cout << beam_fit_dresz_lh << std::endl;
      
      std::cout << "c0 = " << myalgo->GetResPar0() << " +- " << myalgo->GetResPar0Err() << std::endl;
      std::cout << "c1 = " << myalgo->GetResPar1() << " +- " << myalgo->GetResPar1Err() << std::endl;
    }
    */
  }
  else
    if (debug_) std::cout << "No good track selected! No beam fit!" << std::endl;
}
