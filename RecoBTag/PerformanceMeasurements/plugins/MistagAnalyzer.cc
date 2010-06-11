
#include "RecoBTag/PerformanceMeasurements/interface/MistagAnalyzer.h"

MistagAnalyzer::MistagAnalyzer(const edm::ParameterSet& iConfig): classifier_(iConfig)  
{
//now do what ever initialization is needed

// Flavour identification
  flavourMatchOptionf = iConfig.getParameter<std::string>( "flavourMatchOption" );
  if (flavourMatchOptionf == "fastMC") {
    flavourSourcef = iConfig.getParameter<edm::InputTag>("flavourSource");
  } else if (flavourMatchOptionf == "genParticle") {
    flavourSourcef = iConfig.getParameter<edm::InputTag> ("flavourSource");
  }
  
// Parameters
  minJetPt_  = iConfig.getParameter<double>("MinPt");
  maxJetEta_ = iConfig.getParameter<double>("MaxEta");
  selTagger_ = iConfig.getParameter<int>("selTagger");
  tagCut_    = iConfig.getParameter<double>("tagCut");
  vetoPos_   = iConfig.getParameter<double>("vetoPos");
  ntrackMin_ = iConfig.getParameter<int>("ntrackMin");

  isData_    = iConfig.getParameter<bool>("isData");
  useTrackHistory_       = iConfig.getParameter<bool>("useTrackHistory");
  produceJetProbaTree_ = iConfig.getParameter<bool>("produceJetProbaTree");
  
// Modules
  primaryVertexColl_   = iConfig.getParameter<std::string>("primaryVertexColl");

  jetCorrector_          = iConfig.getParameter<std::string>("jetCorrector");
  CaloJetCollectionTags_ = iConfig.getParameter<std::string>("Jets");
 
  trackCHEModuleName_    = iConfig.getParameter<std::string>("trackCHEModuleName");
  trackCNegHEModuleName_ = iConfig.getParameter<std::string>("trackCNegHEModuleName");
 
  trackCHPModuleName_    = iConfig.getParameter<std::string>("trackCHPModuleName");
  trackCNegHPModuleName_ = iConfig.getParameter<std::string>("trackCNegHPModuleName");
  
  jetPModuleName_        = iConfig.getParameter<std::string>("jetPModuleName");
  jetPPosModuleName_     = iConfig.getParameter<std::string>("jetPPosModuleName");
  jetPNegModuleName_     = iConfig.getParameter<std::string>("jetPNegModuleName");
 
  svtxModuleNameHighEff_       = iConfig.getParameter<std::string>("svtxModuleNameHighEff");
  svtxNegModuleNameHighEff_    = iConfig.getParameter<std::string>("svtxNegModuleNameHighEff");
  
  svtxModuleNameHighPur_       = iConfig.getParameter<std::string>("svtxModuleNameHighPur");
  svtxNegModuleNameHighPur_    = iConfig.getParameter<std::string>("svtxNegModuleNameHighPur");
  
  combinedSvtxModuleName_     = iConfig.getParameter<std::string>("combinedSvtxModuleName");
  combinedSvtxNegModuleName_  = iConfig.getParameter<std::string>("combinedSvtxNegModuleName");
  
  softMuonModuleName_       = iConfig.getParameter<std::string>("softMuonModuleName");
  softMuonNegModuleName_    = iConfig.getParameter<std::string>("softMuonNegModuleName");
  softMuonTagInfoName_      = iConfig.getParameter<std::string>("softMuonTagInfoName");
  
///////////////
// Some Histograms

  hData_All_NJets       = fs->make<TH1F>("hData_All_NJets","nb. of jets",21,-0.5,20.5);
  hData_All_NTracks     = fs->make<TH1F>("hData_All_NTracks","nb. of tracks",20,0.5,20.5);
  hData_All_JetPt       = fs->make<TH1F>("hData_All_JetPt","pt(jet)",21,20.,230.);
  hData_All_JetEta      = fs->make<TH1F>("hData_All_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_NJets           = fs->make<TH1F>("hData_NJets","nb. of jets",21,-0.5,20.5);
  hData_NTracks         = fs->make<TH1F>("hData_NTracks","nb. of tracks",20,0.5,20.5);
  hData_JetPt           = fs->make<TH1F>("hData_JetPt","pt(jet)",21,20.,230.);
  hData_JetEta          = fs->make<TH1F>("hData_JetEta","|#eta(jet)|",25,0.,2.5);
  hData_Tagger          = fs->make<TH1F>("hData_Tagger","Tagger",100,-50.,50.);
  hData_Tagger_TCHE     = fs->make<TH1F>("hData_Tagger_TCHE","Tagger_TCHE",100,-50.,50.);
  hData_Tagger_TCHP     = fs->make<TH1F>("hData_Tagger_TCHP","Tagger_TCHP",100,-50.,50.);
  hData_Tagger_JP       = fs->make<TH1F>("hData_Tagger_JP","Tagger_JP",100,-50.,50.);
  hData_Tagger_SSVHE    = fs->make<TH1F>("hData_Tagger_SSV HE ","Tagger_SSV HE",100,-50.,50.);
  hData_Tagger_SSVHP    = fs->make<TH1F>("hData_Tagger_SSV HP","Tagger_SSV HP",100,-50.,50.);
  hData_Tagger_CSV      = fs->make<TH1F>("hData_Tagger_CSV","Tagger_CSV",100,-50.,50.);
  hData_Tagger_MU       = fs->make<TH1F>("hData_Tagger_MU","Tagger_MU",100,-50.,50.);
  
  hAllFlav_Flavour      = fs->make<TH1F>("hAllFlav_Flavour","Flavour",22,-0.5,21.5);
  hAllFlav_Tagger       = fs->make<TH1F>("hAllFlav_Tagger","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Gam   = fs->make<TH1F>("hAllFlav_Tagger_Gam","Tagger",100,-50.,50.);
  hAllFlav_Tagger_K0s   = fs->make<TH1F>("hAllFlav_Tagger_K0s","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Lam   = fs->make<TH1F>("hAllFlav_Tagger_Lam","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Bwd   = fs->make<TH1F>("hAllFlav_Tagger_Bwd","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Cwd   = fs->make<TH1F>("hAllFlav_Tagger_Cwd","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Tau   = fs->make<TH1F>("hAllFlav_Tagger_Tau","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Int   = fs->make<TH1F>("hAllFlav_Tagger_Int","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Fak   = fs->make<TH1F>("hAllFlav_Tagger_Fak","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Bad   = fs->make<TH1F>("hAllFlav_Tagger_Bad","Tagger",100,-50.,50.);
  hAllFlav_Tagger_Oth   = fs->make<TH1F>("hAllFlav_Tagger_Oth","Tagger",100,-50.,50.);
  
  hLightFlav_Tagger     = fs->make<TH1F>("hLightFlav_Tagger","Tagger",100,-50.,50.);
  hGluonFlav_Tagger     = fs->make<TH1F>("hGluonFlav_Tagger","Tagger",100,-50.,50.);
  hUDSFlav_Tagger       = fs->make<TH1F>("hUDSFlav_Tagger","Tagger",100,-50.,50.);
  hCFlav_Tagger         = fs->make<TH1F>("hCFlav_Tagger","Tagger",100,-50.,50.);
  hBFlav_Tagger         = fs->make<TH1F>("hBFlav_Tagger","Tagger",100,-50.,50.);

//$$
  IPSign_cat0           = fs->make<TH1F>("IPSign_cat0","IP/s",100, -20., 0.);
  IPSign_cat1           = fs->make<TH1F>("IPSign_cat1","IP/s",100, -20., 0.);
  IPSign_cat2           = fs->make<TH1F>("IPSign_cat2","IP/s",100, -20., 0.);
  IPSign_cat3           = fs->make<TH1F>("IPSign_cat3","IP/s",100, -20., 0.);
  IPSign_cat4           = fs->make<TH1F>("IPSign_cat4","IP/s",100, -20., 0.);
  IPSign_cat5           = fs->make<TH1F>("IPSign_cat5","IP/s",100, -20., 0.);
  IPSign_cat6           = fs->make<TH1F>("IPSign_cat6","IP/s",100, -20., 0.);
  IPSign_cat7           = fs->make<TH1F>("IPSign_cat7","IP/s",100, -20., 0.);
  IPSign_cat8           = fs->make<TH1F>("IPSign_cat8","IP/s",100, -20., 0.);
  IPSign_cat9           = fs->make<TH1F>("IPSign_cat9","IP/s",100, -20., 0.);
  TrackProbaNeg         = fs->make<TH1F>("TrackProbaNEG","TrackProbaNEG",100, 0., 1.);
  TrackProbaNeg_Cat0    = fs->make<TH1F>("TrackProbaNEG Cat0","TrackProbaNEG Cat0",100, 0., 1.);
  TrackProbaNeg_Cat1    = fs->make<TH1F>("TrackProbaNEG Cat1","TrackProbaNEG Cat1",100, 0., 1.);
  TrackProbaNeg_Cat2    = fs->make<TH1F>("TrackProbaNEG Cat2","TrackProbaNEG Cat2",100, 0., 1.);
  TrackProbaNeg_Cat3    = fs->make<TH1F>("TrackProbaNEG Cat3","TrackProbaNEG Cat3",100, 0., 1.);
  TrackProbaNeg_Cat4    = fs->make<TH1F>("TrackProbaNEG Cat4","TrackProbaNEG Cat4",100, 0., 1.);
  TrackProbaNeg_Cat5    = fs->make<TH1F>("TrackProbaNEG Cat5","TrackProbaNEG Cat5",100, 0., 1.);
  TrackProbaNeg_Cat6    = fs->make<TH1F>("TrackProbaNEG Cat6","TrackProbaNEG Cat6",100, 0., 1.);
  TrackProbaNeg_Cat7    = fs->make<TH1F>("TrackProbaNEG Cat7","TrackProbaNEG Cat7",100, 0., 1.);
  TrackProbaNeg_Cat8    = fs->make<TH1F>("TrackProbaNEG Cat8","TrackProbaNEG Cat8",100, 0., 1.);
  TrackProbaNeg_Cat9    = fs->make<TH1F>("TrackProbaNEG Cat9","TrackProbaNEG Cat9",100, 0., 1.);
//$$

///////////////
// TTree
  
  smalltree = fs->make<TTree>("ttree", "ttree");
  
  smalltree->Branch("BitTrigger",&BitTrigger,"BitTrigger/I");
  smalltree->Branch("nJet",&nJet,"nJet/I");
  smalltree->Branch("Run",&Run,"Run/I");
  smalltree->Branch("Evt",&Evt,"Evt/I");
  smalltree->Branch("LumiBlock",&LumiBlock,"LumiBlock/I");
  
  if ( produceJetProbaTree_ ) {
     
  //--------------------------------------
  //tracks information 
  //--------------------------------------
  smalltree->Branch("nTrack",	       &nTrack, 	"nTrack/I");
  smalltree->Branch("Track_dxy",       Track_dxy,	"Track_dxy[nTrack]/F");
  smalltree->Branch("Track_dz",        Track_dz,	"Track_dz[nTrack]/F");
  smalltree->Branch("Track_zIP",       Track_zIP,	"Track_zIP[nTrack]/F");
  smalltree->Branch("Track_length",    Track_length,	"Track_length[nTrack]/F");
  smalltree->Branch("Track_dist",      Track_dist,	"Track_dist[nTrack]/F");
  smalltree->Branch("Track_IP2D",      Track_IP2D,	"Track_IP2D[nTrack]/F");
//$$
  smalltree->Branch("Track_IP2Dsig",   Track_IP2Dsig,	"Track_IP2Dsig[nTrack]/F");
//$$
  smalltree->Branch("Track_IP",        Track_IP,	"Track_IP[nTrack]/F");
  smalltree->Branch("Track_IPsig",     Track_IPsig,	"Track_IPsig[nTrack]/F");
  smalltree->Branch("Track_Proba",     Track_Proba,	"Track_[nTrack]Proba/F");
//$$
  smalltree->Branch("Track_p",         Track_p,	        "Track_p[nTrack]/F");
//$$
  smalltree->Branch("Track_pt",        Track_pt,	"Track_pt[nTrack]/F");
  smalltree->Branch("Track_eta",       Track_eta,	"Track_eta[nTrack]/F");
  smalltree->Branch("Track_phi",       Track_phi,	"Track_phi[nTrack]/F");
  smalltree->Branch("Track_chi2",      Track_chi2,	"Track_chi2[nTrack]/F");
  smalltree->Branch("Track_charge",    Track_charge,	"Track_charge[nTrack]/I");
  smalltree->Branch("Track_history",   Track_history,	"Track_history[nTrack]/I");

  smalltree->Branch("Track_nHitStrip", Track_nHitStrip, "Track_nHitStrip[nTrack]/I");
  smalltree->Branch("Track_nHitPixel", Track_nHitPixel, "Track_nHitPixel[nTrack]/I");
  smalltree->Branch("Track_nHitAll",   Track_nHitAll,	"Track_nHitAll[nTrack]/I");
  smalltree->Branch("Track_nHitTIB",   Track_nHitTIB,	"Track_nHitTIB[nTrack]/I");
  smalltree->Branch("Track_nHitTID",   Track_nHitTID,	"Track_nHitTID[nTrack]/I");
  smalltree->Branch("Track_nHitTOB",   Track_nHitTOB,	"Track_nHitTOB[nTrack]/I");
  smalltree->Branch("Track_nHitTEC",   Track_nHitTEC,	"Track_nHitTEC[nTrack]/I");
  smalltree->Branch("Track_nHitPXB",   Track_nHitPXB,	"Track_nHitPXB[nTrack]/I");
  smalltree->Branch("Track_nHitPXF",   Track_nHitPXF,	"Track_nHitPXF[nTrack]/I");
  smalltree->Branch("Track_isHitL1",   Track_isHitL1,	"Track_isHitL1[nTrack]/I");
  
//$$
  smalltree->Branch("Track_PV",        Track_PV,        "Track_PV[nTrack]/I");
  smalltree->Branch("Track_SV",        Track_SV,        "Track_SV[nTrack]/I");
  smalltree->Branch("Track_PVweight",  Track_PVweight,  "Track_PVweight[nTrack]/F");
  smalltree->Branch("Track_SVweight",  Track_SVweight,  "Track_SVweight[nTrack]/F");

  smalltree->Branch("Track_category",  Track_category,  "Track_category[nTrack]/I");
//$$

  //--------------------------------------
  //primary vertex information 
  //--------------------------------------
  smalltree->Branch("nPV"	   ,&nPV	 ,"nPV/I");
  smalltree->Branch("PV_x"	   ,PV_x	 ,"PV_x[nPV]/F");
  smalltree->Branch("PV_y"	   ,PV_y	 ,"PV_y[nPV]/F");
  smalltree->Branch("PV_z"	   ,PV_z	 ,"PV_z[nPV]/F");
  smalltree->Branch("PV_ex"	   ,PV_ex	 ,"PV_ex[nPV]/F");
  smalltree->Branch("PV_ey"	   ,PV_ey	 ,"PV_ey[nPV]/F");
  smalltree->Branch("PV_ez"	   ,PV_ez	 ,"PV_ez[nPV]/F");
  smalltree->Branch("PV_chi2"	   ,PV_chi2	 ,"PV_chi2[nPV]/F");
  smalltree->Branch("PV_ndf"	   ,PV_ndf	 ,"PV_ndf[nPV]/F");
  smalltree->Branch("PV_isgood"    ,PV_isgood	 ,"PV_isgood[nPV]/I");
  smalltree->Branch("PV_isfake"    ,PV_isfake	 ,"PV_isfake[nPV]/I");
  
  //--------------------------------------
  //secondary vertex information 
  //--------------------------------------
  smalltree->Branch("nSV"	   ,&nSV	 ,"nSV/I");
  smalltree->Branch("SV_x"	   ,SV_x	 ,"SV_x[nSV]/F");
  smalltree->Branch("SV_y"	   ,SV_y	 ,"SV_y[nSV]/F");
  smalltree->Branch("SV_z"	   ,SV_z	 ,"SV_z[nSV]/F");
  smalltree->Branch("SV_ex"	   ,SV_ex	 ,"SV_ex[nSV]/F");
  smalltree->Branch("SV_ey"	   ,SV_ey	 ,"SV_ey[nSV]/F");
  smalltree->Branch("SV_ez"	   ,SV_ez	 ,"SV_ez[nSV]/F");
  smalltree->Branch("SV_chi2"	   ,SV_chi2	 ,"SV_chi2[nSV]/F");
  smalltree->Branch("SV_ndf"	   ,SV_ndf	 ,"SV_ndf[nSV]/F");
  smalltree->Branch("SV_flight"    ,SV_flight	 ,"SV_flight[nSV]/F");
  smalltree->Branch("SV_flightErr" ,SV_flightErr ,"SV_flightErr[nSV]/F");
  }
    
  //--------------------------------------
  //jets information 
  //--------------------------------------
  smalltree->Branch("Jet_pt",          Jet_pt	       ,"Jet_pt[nJet]/F");
  smalltree->Branch("Jet_jes",         Jet_jes         ,"Jet_jes[nJet]/F");
  smalltree->Branch("Jet_eta",         Jet_eta         ,"Jet_eta[nJet]/F");
  smalltree->Branch("Jet_phi",         Jet_phi         ,"Jet_phi[nJet]/F");
  smalltree->Branch("Jet_ntracks",     Jet_ntracks     ,"Jet_ntracks[nJet]/I");
  smalltree->Branch("Jet_flavour",     Jet_flavour     ,"Jet_flavour[nJet]/I");
  smalltree->Branch("Jet_Ip1N",        Jet_Ip1N        ,"Jet_Ip1N[nJet]/F");
  smalltree->Branch("Jet_Ip1P",        Jet_Ip1P        ,"Jet_Ip1P[nJet]/F");
  smalltree->Branch("Jet_Ip2N",        Jet_Ip2N        ,"Jet_Ip2N[nJet]/F");
  smalltree->Branch("Jet_Ip2P",        Jet_Ip2P        ,"Jet_Ip2P[nJet]/F");
  smalltree->Branch("Jet_Ip3N",        Jet_Ip3N        ,"Jet_Ip3N[nJet]/F");
  smalltree->Branch("Jet_Ip3P",        Jet_Ip3P        ,"Jet_Ip3P[nJet]/F");
  smalltree->Branch("Jet_ProbaN",      Jet_ProbaN      ,"Jet_ProbaN[nJet]/F");
  smalltree->Branch("Jet_ProbaP",      Jet_ProbaP      ,"Jet_ProbaP[nJet]/F");
  smalltree->Branch("Jet_Proba",       Jet_Proba       ,"Jet_Proba[nJet]/F");
//$$
//   smalltree->Branch("Jet_TkProba",     Jet_TkProba     ,"Jet_TkProba[nJet]/F"); // these are the distribution of JetProba computed by track information
//   smalltree->Branch("Jet_TkProbaN",    Jet_TkProbaN    ,"Jet_TkProbaN[nJet]/F");
//   smalltree->Branch("Jet_TkProbaP",    Jet_TkProbaP    ,"Jet_TkProbaP[nJet]/F");
//$$
  smalltree->Branch("Jet_SvxN",        Jet_SvxN        ,"Jet_SvxN[nJet]/F");
  smalltree->Branch("Jet_Svx",         Jet_Svx         ,"Jet_Svx[nJet]/F");
  smalltree->Branch("Jet_SvxNTracks",  Jet_SvxNTracks  ,"Jet_SvxNTracks[nJet]/I");
  smalltree->Branch("Jet_SvxTracks",   Jet_SvxTracks   ,"Jet_SvxTracks[nJet]/I");
  smalltree->Branch("Jet_SvxNHP",      Jet_SvxNHP      ,"Jet_SvxNHP[nJet]/F");
  smalltree->Branch("Jet_SvxHP",       Jet_SvxHP       ,"Jet_SvxHP[nJet]/F");
  smalltree->Branch("Jet_CombSvxN",    Jet_CombSvxN    ,"Jet_CombSvxN[nJet]/F");
  smalltree->Branch("Jet_CombSvx",     Jet_CombSvx     ,"Jet_CombSvx[nJet]/F");
  smalltree->Branch("Jet_SoftMuN",     Jet_SoftMuN     ,"Jet_SoftMuN[nJet]/F");
  smalltree->Branch("Jet_SoftMu",      Jet_SoftMu      ,"Jet_SoftMu[nJet]/F");
  smalltree->Branch("Jet_hist1",       Jet_hist1       ,"Jet_hist1[nJet]/I");
  smalltree->Branch("Jet_hist2",       Jet_hist2       ,"Jet_hist2[nJet]/I");
  smalltree->Branch("Jet_hist3",       Jet_hist3       ,"Jet_hist3[nJet]/I");
  smalltree->Branch("Jet_histJet",     Jet_histJet     ,"Jet_histJet[nJet]/I");
  smalltree->Branch("Jet_histSvx",     Jet_histSvx     ,"Jet_histSvx[nJet]/I");
  smalltree->Branch("Jet_histMuon",    Jet_histMuon    ,"Jet_histMuon[nJet]/I");
  smalltree->Branch("Jet_mu_nHit",     Jet_mu_nHit     ,"Jet_mu_nHit[nJet]/I");
  smalltree->Branch("Jet_mu_chi2",     Jet_mu_chi2     ,"Jet_mu_chi2[nJet]/F");
  smalltree->Branch("Jet_mu_pt",       Jet_mu_pt       ,"Jet_mu_pt[nJet]/F");
  smalltree->Branch("Jet_mu_ptrel",    Jet_mu_ptrel    ,"Jet_mu_ptrel[nJet]/F");
  smalltree->Branch("Jet_nFirstTrack", Jet_nFirstTrack ,"Jet_nFirstTrack[nJet]/I");
  smalltree->Branch("Jet_nLastTrack",  Jet_nLastTrack  ,"Jet_nLastTrack[nJet]/I"); 
  smalltree->Branch("Jet_nFirstSV",    Jet_nFirstSV    ,"Jet_nFirstSV[nJet]/I");
  smalltree->Branch("Jet_nLastSV",     Jet_nLastSV     ,"Jet_nLastSV[nJet]/I");
    
}

 
MistagAnalyzer::~MistagAnalyzer() 
{
}


static std::vector<std::size_t> sortedIndexes(std::vector<TrackIPTagInfo::TrackIPData > const & values)
{
  std::multimap<float, std::size_t> sortedIdx;
  std::vector<std::size_t> result;
  
  for(size_t i = 0; i < values.size(); ++i)
    sortedIdx.insert( std::pair<float, std::size_t>(values[i].ip3d.significance() , i) );
  
  for(std::multimap<float, std::size_t>::reverse_iterator it = sortedIdx.rbegin(); it != sortedIdx.rend(); ++it)
    result.push_back(it->second);
 
  return result;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void MistagAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  //------------------------------------------------------
  //Event informations
  //------------------------------------------------------
  Run = iEvent.id().run();
  
  if ( !isData_ && Run > 0 ) Run = -Run;
  
  Evt  = iEvent.id().event();
  LumiBlock  = iEvent.luminosityBlock();  
  
  // Tag Jets
  if ( useTrackHistory_ ) classifier_.newEvent(iEvent, iSetup);
  
  const JetCorrector *acorrector = JetCorrector::getJetCorrector(jetCorrector_,iSetup);
  // Calo Jets
  //Handle<reco::CaloJetCollection> jetsColl;
  //iEvent.getByLabel(CaloJetCollectionTags_, jetsColl);
  //const reco::CaloJetCollection recoJets =   *(jetsColl.product());
  
  edm::Handle <edm::View <reco::Jet> > jetsCollHandle;
  iEvent.getByLabel (CaloJetCollectionTags_, jetsCollHandle);
  edm::View<reco::Jet> jetsColl = *jetsCollHandle;


  // initialize flavour identifiers
  edm::Handle<JetFlavourMatchingCollection> jetMC;
 
 
  //------------------------------------------------------
  //TrackCounting taggers
  //------------------------------------------------------
  Handle<std::vector<reco::TrackIPTagInfo> > tagInfo;
  iEvent.getByLabel("impactParameterTagInfos", tagInfo);
  
  edm::Handle<reco::JetTagCollection> jetTags_TCHighEff;
  iEvent.getByLabel(trackCHEModuleName_, jetTags_TCHighEff);
  edm::Handle<reco::JetTagCollection> jetTags_NegTCHighEff;
  iEvent.getByLabel(trackCNegHEModuleName_, jetTags_NegTCHighEff);
  
  edm::Handle<reco::JetTagCollection> jetTags_TCHighPur;
  iEvent.getByLabel(trackCHPModuleName_, jetTags_TCHighPur);
  edm::Handle<reco::JetTagCollection> jetTags_NegTCHighPur;
  iEvent.getByLabel(trackCNegHPModuleName_, jetTags_NegTCHighPur);
  
  //------------------------------------------------------
  //Jet Probability tagger
  //------------------------------------------------------
  edm::Handle<reco::JetTagCollection> jetTags_JP;
  iEvent.getByLabel(jetPModuleName_, jetTags_JP);
  
  edm::Handle<reco::JetTagCollection> jetTags_PosJP;
  iEvent.getByLabel(jetPPosModuleName_, jetTags_PosJP);
  
  edm::Handle<reco::JetTagCollection> jetTags_NegJP;
  iEvent.getByLabel(jetPNegModuleName_, jetTags_NegJP);
  
  //------------------------------------------------------
  //Secondary vertex taggers
  //------------------------------------------------------
  Handle<std::vector<reco::SecondaryVertexTagInfo> > tagInfoSVx;
  iEvent.getByLabel("secondaryVertexTagInfos", tagInfoSVx);
  
  Handle<std::vector<reco::SecondaryVertexTagInfo> > tagInfoNegSVx;
  iEvent.getByLabel("secondaryVertexNegativeTagInfos", tagInfoNegSVx);
    
  edm::Handle<reco::JetTagCollection> jetTags_SvtxHighEff;
  iEvent.getByLabel(svtxModuleNameHighEff_, jetTags_SvtxHighEff);
  edm::Handle<reco::JetTagCollection> jetTags_negSvtxHighEff;
  iEvent.getByLabel(svtxNegModuleNameHighEff_, jetTags_negSvtxHighEff);
  
  edm::Handle<reco::JetTagCollection> jetTags_SvtxHighPur;
  iEvent.getByLabel(svtxModuleNameHighPur_, jetTags_SvtxHighPur);
  edm::Handle<reco::JetTagCollection> jetTags_negSvtxHighPur;
  iEvent.getByLabel(svtxNegModuleNameHighPur_, jetTags_negSvtxHighPur);
  
  edm::Handle<reco::JetTagCollection> jetTags_CombinedSvtx;
  iEvent.getByLabel(combinedSvtxModuleName_, jetTags_CombinedSvtx);
  edm::Handle<reco::JetTagCollection> jetTags_negCombinedSvtx;
  iEvent.getByLabel(combinedSvtxNegModuleName_, jetTags_negCombinedSvtx);
  
  //------------------------------------------------------
  //Soft muon tagger
  //------------------------------------------------------
  
  edm::Handle<reco::JetTagCollection> jetTags_softMu;
  iEvent.getByLabel(softMuonModuleName_, jetTags_softMu);
  edm::Handle<reco::JetTagCollection> jetTags_softMuneg;
  iEvent.getByLabel(softMuonNegModuleName_, jetTags_softMuneg);
  
  edm::Handle<reco::SoftLeptonTagInfoCollection> tagInos_softmuon;
  iEvent.getByLabel(softMuonTagInfoName_, tagInos_softmuon);
  

  //------------------
  //Primary vertex
  //------------------
  Handle<reco::VertexCollection> primaryVertex;
  iEvent.getByLabel(primaryVertexColl_,primaryVertex);
  
  const  reco::Vertex  *pv;
  bool newvertex = false;
  
  bool pvFound = (primaryVertex->size() != 0);
  if ( pvFound ) {
    pv = &(*primaryVertex->begin());
  }
  else { 
    reco::Vertex::Error e;
    e(0,0)=0.0015*0.0015;
    e(1,1)=0.0015*0.0015;
    e(2,2)=15.*15.;
    reco::Vertex::Point p(0,0,0);
    pv=  new reco::Vertex(p,e,1,1,1);
    newvertex = true;
  }
  GlobalPoint Pv_point = GlobalPoint((*pv).x(), (*pv).y(), (*pv).z());
  
  nPV=0;
  for(unsigned int i = 0; i< primaryVertex->size() ; i++){
    PV_x[nPV]    = (*primaryVertex)[i].x();
    PV_y[nPV]    = (*primaryVertex)[i].y();
    PV_z[nPV]    = (*primaryVertex)[i].z();
    PV_ex[nPV]   = (*primaryVertex)[i].xError();
    PV_ey[nPV]   = (*primaryVertex)[i].yError();
    PV_ez[nPV]   = (*primaryVertex)[i].zError();
    PV_chi2[nPV]   = (*primaryVertex)[i].normalizedChi2();
    PV_ndf[nPV]  = (*primaryVertex)[i].ndof();
    PV_isgood[nPV]    = (*primaryVertex)[i].isValid();
    PV_isfake[nPV]    = (*primaryVertex)[i].isFake();
    
    nPV++;
  }


  //------------------------------------------------------
  //Trigger info
  //------------------------------------------------------
  
  TriggerResults tr;
  Handle<TriggerResults> h_trigRes;
  iEvent.getByLabel(InputTag("TriggerResults::HLT"), h_trigRes);
  tr = *h_trigRes;

  BitTrigger = 0;

//$$   if ( TriggerInfo_ ) { 
    vector<string> triggerList;
    Service<service::TriggerNamesService> tns;
    bool foundNames = tns->getTrigPaths(tr,triggerList);
    if ( !foundNames ) cout << "Could not get trigger names!\n";
    if ( tr.size() != triggerList.size() ) cout << "ERROR: length of names and paths not the same: " << triggerList.size() << "," << tr.size() << endl;
  
    for (unsigned int i=0; i< tr.size(); i++) {
    if ( !tr[i].accept() == 1 ) continue;
// Trigger table 2010
      if ( triggerList[i] == "HLT_L1Jet6U"  	    ) BitTrigger +=1 ;
      if ( triggerList[i] == "HLT_L1Jet10U"  	    ) BitTrigger +=2 ;
      if ( triggerList[i] == "HLT_Jet15U"  	    ) BitTrigger +=10 ;	
      if ( triggerList[i] == "HLT_Jet30U"  	    ) BitTrigger +=20 ;	
      if ( triggerList[i] == "HLT_Jet50U"           ) BitTrigger +=40 ;   
      if ( triggerList[i] == "HLT_DiJetAve15U_8E29" ) BitTrigger +=100 ;   
      if ( triggerList[i] == "HLT_DiJetAve30U_8E29" ) BitTrigger +=200 ;   
      if ( triggerList[i] == "HLT_L1MuOpen" ) BitTrigger +=1000 ;   
      if ( triggerList[i] == "HLT_L1Mu"     ) BitTrigger +=2000 ;   
      if ( triggerList[i] == "HLT_L1Mu20"   ) BitTrigger +=4000 ;   
      if ( triggerList[i] == "HLT_L2Mu0"    ) BitTrigger +=10000 ;   
      if ( triggerList[i] == "HLT_L2Mu3"    ) BitTrigger +=20000 ;   
      if ( triggerList[i] == "HLT_L2Mu9"    ) BitTrigger +=40000 ;   
      if ( triggerList[i] == "HLT_Mu3"      ) BitTrigger +=100000 ;   
      if ( triggerList[i] == "HLT_Mu5"      ) BitTrigger +=200000 ;   
      if ( triggerList[i] == "HLT_Mu9"      ) BitTrigger +=400000 ;   
// std::cout << " Run Evt " << Run << " " << Evt << " trigger list " << triggerList[i] << std::endl;
    }
// std::cout << " Run Evt " << Run << " " << Evt << " bit trigger " << BitTrigger << std::endl;
//$$   }

    
  if ( flavourMatchOptionf == "fastMC" && !isData_ ) {
    iEvent.getByLabel(flavourSourcef, jetMC);
    for(JetFlavourMatchingCollection::const_iterator iter =
	  jetMC->begin(); iter != jetMC->end(); iter++)
      flavoursMapf.insert(std::pair <const edm::RefToBase<reco::Jet>, unsigned int>((*iter).first, (unsigned int)((*iter).second).getFlavour()));
  } else if (flavourMatchOptionf == "genParticle" && !isData_) {
    iEvent.getByLabel (flavourSourcef, theJetPartonMapf);
  }
  

  int cap0, cap1, cap2, cap3, cap4, cap5, cap6, cap7, cap8;
  int can0, can1, can2, can3, can4, can5, can6, can7, can8;

  nTrack = 0;
  nSV = 0; 
  int Njets = 0;
  int numjet = 0;
  nJet = 0;

  //*********************************
  // Loop over the jets
  
  //CaloJetCollection::const_iterator jet;
      
  for ( unsigned int ijet = 0; ijet < jetsColl.size(); ijet++) {
    
    double JES = acorrector->correction(jetsColl.at(ijet) );
    double jetpt  = (jetsColl.at(ijet)).pt()  ;
    double jeteta = (jetsColl.at(ijet)).eta() ;
    double ptjet  = jetpt* JES;
    
//$$
    if ( jetpt*JES < 10. || std::fabs( jeteta ) > 2.5 ) continue;
//$$
    if ( jetpt*JES > minJetPt_ && std::fabs( jeteta ) < maxJetEta_ ) Njets++; 

    int ith_tagged = -1;      
    float Flavour  =-1;    
    int flavour  =-1  ; 
    
    if ( !isData_ ) {
      flavour = abs(getMatchedParton(jetsColl.at(ijet)).getFlavour());   
      if ( flavour >= 1 && flavour <= 3 ) flavour = 1;  
      Flavour = flavour;
    }
    
    Jet_flavour[nJet] = Flavour;
    Jet_eta[nJet]     = (jetsColl.at(ijet)).eta();
    Jet_phi[nJet]     = (jetsColl.at(ijet)).phi();
    Jet_pt[nJet]      = (jetsColl.at(ijet)).pt();
    
    Jet_jes[nJet]     = JES;
    
    float etajet = TMath::Abs( (jetsColl.at(ijet)).eta());
    float phijet = (jetsColl.at(ijet)).phi();
    if (phijet < 0.) phijet += 2*TMath::Pi();
     

//*****************************************************************
// Taggers 
//*****************************************************************
   
    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_JP);
    
// Loop on Selected Tracks
//
    int ntagtracks = (*tagInfo)[ith_tagged].probabilities(0).size();
    Jet_ntracks[nJet] = ntagtracks;

    Jet_nFirstTrack[nJet]  = nTrack;
    int k=0;
    
    const edm::RefVector<reco::TrackCollection> & assotracks((*tagInfo)[ith_tagged].selectedTracks());

    if ( produceJetProbaTree_ ) {
      for (unsigned int itt=0; itt < assotracks.size(); itt++) {
	//(*tagInfo)[ith_tagged].probability(itt,0);

//$$
	const reco::Track  &ptrack = (*(assotracks[itt]));
//$$
        
        GlobalPoint maximumClose = (*tagInfo)[ith_tagged].impactParameterData()[k].closestToJetAxis;
        float decayLen = (maximumClose - (Pv_point)).mag(); 
	float distJetAxis = (*tagInfo)[ith_tagged].impactParameterData()[k].distanceToJetAxis.value();
        float deta = (assotracks[itt])->eta() - Jet_eta[nJet];
        float dphi = (assotracks[itt])->phi() - Jet_phi[nJet];
        if ( dphi > TMath::Pi() ) dphi = 2.*TMath::Pi() - dphi;
        float deltaR = TMath::Sqrt(deta*deta + dphi*dphi);
	
      // track selection
//$$
	if ( std::fabs(distJetAxis) < 0.07 && decayLen < 5.0 && deltaR < 0.3 ) {
//$$
	  Track_dxy[nTrack]      = (assotracks[itt])->dxy(pv->position());
	  Track_dz[nTrack]       = (assotracks[itt])->dz(pv->position());
	  Track_zIP[nTrack]      = (assotracks[itt])->dz()-(*pv).z();	
	  Track_length[nTrack]   = decayLen;
	  Track_dist[nTrack]     = distJetAxis;
	  Track_IP2D[nTrack]     = (*tagInfo)[ith_tagged].impactParameterData()[k].ip2d.value();
//$$
	  Track_IP2Dsig[nTrack]  = (*tagInfo)[ith_tagged].impactParameterData()[k].ip2d.significance();
//$$
	  Track_IP[nTrack]       = (*tagInfo)[ith_tagged].impactParameterData()[k].ip3d.value();
	  Track_IPsig[nTrack]    = (*tagInfo)[ith_tagged].impactParameterData()[k].ip3d.significance();
	  Track_Proba[nTrack]    = (*tagInfo)[ith_tagged].probabilities(0)[k];
	  
//$$
	  Track_p[nTrack]        = (assotracks[itt])->p();
//$$
	  Track_pt[nTrack]       = (assotracks[itt])->pt();
	  Track_eta[nTrack]      = (assotracks[itt])->eta();
	  Track_phi[nTrack]      = (assotracks[itt])->phi();
// 	  Track_chi2[nTrack]     = (assotracks[itt])->chi2();    
	  Track_chi2[nTrack]     = (assotracks[itt])->normalizedChi2();   
	  Track_charge[nTrack]   = (assotracks[itt])->charge();
	  
	  Track_nHitAll[nTrack]  = (assotracks[itt])->numberOfValidHits();
	  Track_nHitPixel[nTrack]= (assotracks[itt])->hitPattern().numberOfValidPixelHits();
	  Track_nHitStrip[nTrack]= (assotracks[itt])->hitPattern().numberOfValidStripHits();
	  Track_nHitTIB[nTrack]  = (assotracks[itt])->hitPattern().numberOfValidStripTIBHits();
	  Track_nHitTID[nTrack]  = (assotracks[itt])->hitPattern().numberOfValidStripTIDHits();
	  Track_nHitTOB[nTrack]  = (assotracks[itt])->hitPattern().numberOfValidStripTOBHits();
	  Track_nHitTEC[nTrack]  = (assotracks[itt])->hitPattern().numberOfValidStripTECHits();
	  Track_nHitPXB[nTrack]  = (assotracks[itt])->hitPattern().numberOfValidPixelBarrelHits();
	  Track_nHitPXF[nTrack]  = (assotracks[itt])->hitPattern().numberOfValidPixelEndcapHits();
	  Track_isHitL1[nTrack]  = (assotracks[itt])->hitPattern().hasValidHitInFirstPixelBarrel();
	  
	  Track_PV[nTrack] = 0;
	  Track_SV[nTrack] = 0;
	  Track_PVweight[nTrack] = 0.;
	  Track_SVweight[nTrack] = 0.;
	  
	  Track_history[nTrack] = 0;
          if ( useTrackHistory_ && !isData_ ) {
            TrackCategories::Flags theFlag = classifier_.evaluate( (*tagInfo)[ith_tagged].selectedTracks()[k] ).flags();
	
	    if      ( theFlag[TrackCategories::BWeakDecay] )         Track_history[nTrack] += pow(10, -1 + 1); 
	    else if ( theFlag[TrackCategories::CWeakDecay] )         Track_history[nTrack] += pow(10, -1 + 2); 
	    else if ( theFlag[TrackCategories::TauDecay] )           Track_history[nTrack] += pow(10, -1 + 3); 
	    else if ( theFlag[TrackCategories::ConversionsProcess] ) Track_history[nTrack] += pow(10, -1 + 4); 
	    else if ( theFlag[TrackCategories::KsDecay] )            Track_history[nTrack] += pow(10, -1 + 5); 
	    else if ( theFlag[TrackCategories::LambdaDecay] )        Track_history[nTrack] += pow(10, -1 + 6); 
	    else if ( theFlag[TrackCategories::HadronicProcess] )    Track_history[nTrack] += pow(10, -1 + 7); 
	    else if ( theFlag[TrackCategories::Fake] )               Track_history[nTrack] += pow(10, -1 + 8); 
	    else if ( theFlag[TrackCategories::SharedInnerHits] )    Track_history[nTrack] += pow(10, -1 + 9); 
          }

//$$
// Track categories for Jet Probability calibration
	  if ( Track_IPsig[nTrack] < 0 ) TrackProbaNeg->Fill(-Track_Proba[nTrack]);

	  std::vector<float > theVectOrdered  = getTrackProbabilies(((*tagInfo)[ith_tagged].probabilities(0)), 0);
	  Jet_TkProba[nJet]  = calculProbability( theVectOrdered );
	  std::vector<float > theVectOrdered2 = getTrackProbabilies(((*tagInfo)[ith_tagged].probabilities(0)), 1);
	  Jet_TkProbaP[nJet] = calculProbability( theVectOrdered2 );
	  std::vector<float > theVectOrdered3 = getTrackProbabilies(((*tagInfo)[ith_tagged].probabilities(0)), 2);
	  Jet_TkProbaN[nJet] = calculProbability( theVectOrdered3 );
  
	  Track_category[nTrack] = -1;

	  if ( findCat( &ptrack, *&cat0 ) ) {
	    Track_category[nTrack] = 0;
	    if ( Track_IPsig[nTrack] < 0 ) { 
	      IPSign_cat0->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat0->Fill( -Track_Proba[nTrack] );
	    }
	  }

	  if ( findCat( &ptrack, *&cat1 ) && Track_category[nTrack] != 0 ) {
	    Track_category[nTrack]  = 1;
	    if ( Track_IPsig[nTrack] < 0 ) {
	      IPSign_cat1->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat1->Fill( -Track_Proba[nTrack] );
	    }
	  }
	  
	  if ( findCat( &ptrack, *&cat2 ) && Track_category[nTrack] != 0 && Track_category[nTrack] != 1 ) {
	    Track_category[nTrack] = 2;
	    if ( Track_IPsig[nTrack] < 0 ) {
	      IPSign_cat2->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat2->Fill( -Track_Proba[nTrack] );
	    }
	  }
	  
	  if ( findCat( &ptrack, *&cat3 ) && Track_category[nTrack] != 0 && Track_category[nTrack] != 1 && Track_category[nTrack] != 2 ) {
	    Track_category[nTrack] = 3;
	    if ( Track_IPsig[nTrack] < 0 ) {
	      IPSign_cat3->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat3->Fill( -Track_Proba[nTrack] );
	    }
	  }
	  
	  if ( findCat( &ptrack, *&cat4 ) && Track_category[nTrack] != 0 && Track_category[nTrack] != 1 && Track_category[nTrack] != 2 && Track_category[nTrack] != 3 ) {
	    Track_category[nTrack] = 4;
	    if ( Track_IPsig[nTrack] < 0 ) {
	      IPSign_cat4->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat4->Fill( -Track_Proba[nTrack] );
	    }
	  }
	  
	  if ( findCat( &ptrack, *&cat5 ) && Track_category[nTrack] != 0 && Track_category[nTrack] != 1 && Track_category[nTrack] != 2 && Track_category[nTrack] != 3 && Track_category[nTrack] != 4 ) {
	    Track_category[nTrack] = 5;
	    if ( Track_IPsig[nTrack] < 0 ) {
	      IPSign_cat5->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat5->Fill( -Track_Proba[nTrack] );
	    }
	  }
	  
	  if ( findCat( &ptrack, *&cat6 ) && Track_category[nTrack] != 0 && Track_category[nTrack] != 1 && Track_category[nTrack] != 2 && Track_category[nTrack] != 3 && Track_category[nTrack] != 4 && Track_category[nTrack] != 5 ) {
	    Track_category[nTrack] = 6;
	    if ( Track_IPsig[nTrack] < 0 ) {
	      IPSign_cat6->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat6->Fill( -Track_Proba[nTrack] );
	    }
	  }
	  
	  if ( findCat( &ptrack, *&cat7 ) && Track_category[nTrack] != 0 && Track_category[nTrack] != 1 && Track_category[nTrack] != 2 && Track_category[nTrack] != 3 && Track_category[nTrack] != 4 && Track_category[nTrack] != 5 && Track_category[nTrack] != 6 ) {
	    Track_category[nTrack] = 7;
	    if ( Track_IPsig[nTrack] < 0 ) {
	      IPSign_cat7->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat7->Fill( -Track_Proba[nTrack] );
	    }
	  }
	  
	  if ( findCat( &ptrack, *&cat8 ) && Track_category[nTrack] != 0 && Track_category[nTrack] != 1 && Track_category[nTrack] != 2 && Track_category[nTrack] != 3 && Track_category[nTrack] != 4 && Track_category[nTrack] != 5 && Track_category[nTrack] != 6 && Track_category[nTrack] != 7 ) {
	    Track_category[nTrack] = 8;
	    if ( Track_IPsig[nTrack] < 0 ) {
	      IPSign_cat8->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat8->Fill( -Track_Proba[nTrack] );
	    }
	  }
	  
	  if ( findCat( &ptrack, *&cat9 ) && Track_category[nTrack] != 0 && Track_category[nTrack] != 1 && Track_category[nTrack] != 2 && Track_category[nTrack] != 3 && Track_category[nTrack]!= 4 && Track_category[nTrack] != 5 && Track_category[nTrack] != 6 && Track_category[nTrack] != 7 && Track_category[nTrack] != 8 ) {
	    Track_category[nTrack] = 9;
	    if ( Track_IPsig[nTrack] < 0 ) {
	      IPSign_cat9->Fill( Track_IPsig[nTrack] );
	      TrackProbaNeg_Cat9->Fill( -Track_Proba[nTrack] );
	    }
	  }
//$$

	  nTrack++;
	}
	k++;
      } // end loop on tracks
    }
    Jet_nLastTrack[nJet]   = nTrack;
    
// b-tagger discriminants
    float Proba  = (*jetTags_JP)[ith_tagged].second;
    ith_tagged   = this->TaggedJet(jetsColl.at(ijet),jetTags_PosJP);
    float ProbaP = (*jetTags_PosJP)[ith_tagged].second;
    ith_tagged   = this->TaggedJet(jetsColl.at(ijet),jetTags_NegJP);
    float ProbaN = (*jetTags_NegJP)[ith_tagged].second;
    
    ith_tagged          = this->TaggedJet(jetsColl.at(ijet),jetTags_CombinedSvtx);
    float CombinedSvtx  = (*jetTags_CombinedSvtx)[ith_tagged].second;
    ith_tagged          = this->TaggedJet(jetsColl.at(ijet),jetTags_negCombinedSvtx);
    float CombinedSvtxN = (*jetTags_negCombinedSvtx)[ith_tagged].second;
    
    //reco::TaggingVariableList variables =
    //                   svComputer((*tagInfo)[ith_tagged], (*tagInfoSVx)[ith_tagged]);
    //int vertexCat = variables[reco::btau::vertexCategory];
    //int vertexCat = -1;
    
    ith_tagged    = this->TaggedJet(jetsColl.at(ijet),jetTags_SvtxHighEff);
    float Svtx    = (*jetTags_SvtxHighEff)[ith_tagged].second;
    ith_tagged    = this->TaggedJet(jetsColl.at(ijet),jetTags_negSvtxHighEff);
    float SvtxN   = (*jetTags_negSvtxHighEff)[ith_tagged].second;

    ith_tagged    = this->TaggedJet(jetsColl.at(ijet),jetTags_SvtxHighPur);
    float SvtxHP  = (*jetTags_SvtxHighPur)[ith_tagged].second;
    ith_tagged    = this->TaggedJet(jetsColl.at(ijet),jetTags_negSvtxHighPur);
    float SvtxNHP = (*jetTags_negSvtxHighPur)[ith_tagged].second;
    
    float SoftM  = 0;
    float SoftMN = 0;
    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_softMu);
    if ( (*jetTags_softMu)[ith_tagged].second     > -100000 )   
      SoftM  = (*jetTags_softMu)[ith_tagged].second;
    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_softMuneg);
    if ( (*jetTags_softMuneg)[ith_tagged].second  > -100000 )   
      SoftMN = ((*jetTags_softMuneg)[ith_tagged].second);
    if ( SoftMN > 0 ) SoftMN = -SoftMN;
    
    Jet_ProbaN[nJet]   = ProbaN;	 
    Jet_ProbaP[nJet]   = ProbaP;	  
    Jet_Proba[nJet]    = Proba; 	
    Jet_SvxN[nJet]     = SvtxN;        
    Jet_Svx[nJet]      = Svtx;        
    Jet_SvxNHP[nJet]   = SvtxNHP;        
    Jet_SvxHP[nJet]    = SvtxHP;   	
    Jet_CombSvxN[nJet] = CombinedSvtxN;
    Jet_CombSvx[nJet]  = CombinedSvtx;  
    Jet_SoftMuN[nJet]  = SoftMN;	 
    Jet_SoftMu[nJet]   = SoftM; 	 
    
    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_TCHighEff);
    std::vector<TrackIPTagInfo::TrackIPData>  ipdata = (*tagInfo)[ith_tagged].impactParameterData();
    TrackRefVector tracks( (*tagInfo)[ith_tagged].selectedTracks() );
    std::vector<std::size_t> indexes( sortedIndexes(ipdata) );
    
    TrackCategories::Flags flags1P;
    TrackCategories::Flags flags2P;
    TrackCategories::Flags flags3P;
    TrackCategories::Flags flags1N;
    TrackCategories::Flags flags2N;
    TrackCategories::Flags flags3N;
    int idSize = 0;
    
    if ( useTrackHistory_ && indexes.size() != 0 && !isData_ ) {
      
      idSize = indexes.size();
      flags1P = classifier_.evaluate( tracks[indexes[0]] ).flags();
      if(idSize > 1) flags2P = classifier_.evaluate( tracks[indexes[1]] ).flags();
      if(idSize > 2) flags3P = classifier_.evaluate( tracks[indexes[2]] ).flags();
      flags1N = classifier_.evaluate( tracks[indexes[idSize-1]] ).flags();
      if(idSize > 1) flags2N = classifier_.evaluate( tracks[indexes[idSize-2]] ).flags();
      if(idSize > 2) flags3N = classifier_.evaluate( tracks[indexes[idSize-3]] ).flags();
    }
    
    Jet_Ip1P[nJet]  = -1000;
    Jet_Ip1N[nJet]  =  1000;
    for (std::vector<TrackIPTagInfo::TrackIPData>::const_iterator itipdata = ipdata.begin();
	 itipdata != ipdata.end(); itipdata++) {

      float decayLen = ((*itipdata).closestToJetAxis - (Pv_point)).mag(); 
      float distJetAxis = (*itipdata).distanceToJetAxis.value();
//$$
    if ( abs(distJetAxis) > 0.07 || decayLen > 5.0 ) continue;
//$$	 
      double ip3D = (*itipdata).ip3d.significance();
      if( ip3D > Jet_Ip1P[nJet] ) Jet_Ip1P[nJet] = ip3D;         
      if( ip3D < Jet_Ip1N[nJet] ) Jet_Ip1N[nJet] = ip3D;
    }
    Jet_Ip1N[nJet] = -Jet_Ip1N[nJet];
    
    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_TCHighEff);
    Jet_Ip2P[nJet]   = (*jetTags_TCHighEff)[ith_tagged].second;
    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_TCHighPur);
    Jet_Ip3P[nJet]   = (*jetTags_TCHighPur)[ith_tagged].second;
    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_NegTCHighEff);
    Jet_Ip2N[nJet]   = (*jetTags_NegTCHighEff)[ith_tagged].second;
    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_NegTCHighPur);
    Jet_Ip3N[nJet]   = (*jetTags_NegTCHighPur)[ith_tagged].second;
    
    //*****************************************************************
    //get track histories for 1st, 2nd and 3rd track (TC)
    //*****************************************************************
    Jet_hist1[nJet] = 0;
    Jet_hist2[nJet] = 0;
    Jet_hist3[nJet] = 0;
    
    // Track history
    if ( useTrackHistory_ && indexes.size()!=0 && !isData_ ) {
      if ( flags1P[TrackCategories::BWeakDecay] )         Jet_hist1[nJet] += int(pow(10., -1 + 1)); 
      if ( flags1P[TrackCategories::CWeakDecay] )         Jet_hist1[nJet] += int(pow(10., -1 + 2)); 
      if ( flags1P[TrackCategories::TauDecay] )           Jet_hist1[nJet] += int(pow(10., -1 + 3)); 
      if ( flags1P[TrackCategories::ConversionsProcess] ) Jet_hist1[nJet] += int(pow(10., -1 + 4)); 
      if ( flags1P[TrackCategories::KsDecay] )            Jet_hist1[nJet] += int(pow(10., -1 + 5)); 
      if ( flags1P[TrackCategories::LambdaDecay] )        Jet_hist1[nJet] += int(pow(10., -1 + 6)); 
      if ( flags1P[TrackCategories::HadronicProcess] )    Jet_hist1[nJet] += int(pow(10., -1 + 7)); 
      if ( flags1P[TrackCategories::Fake] )	          Jet_hist1[nJet] += int(pow(10., -1 + 8)); 
      if ( flags1P[TrackCategories::SharedInnerHits] )	  Jet_hist1[nJet] += int(pow(10., -1 + 9)); 
      if ( idSize > 1 ) {
        if ( flags2P[TrackCategories::BWeakDecay] )         Jet_hist2[nJet] += int(pow(10., -1 + 1)); 
        if ( flags2P[TrackCategories::CWeakDecay] )         Jet_hist2[nJet] += int(pow(10., -1 + 2)); 
        if ( flags2P[TrackCategories::TauDecay] )           Jet_hist2[nJet] += int(pow(10., -1 + 3)); 
	if ( flags2P[TrackCategories::ConversionsProcess] ) Jet_hist2[nJet] += int(pow(10., -1 + 4)); 
        if ( flags2P[TrackCategories::KsDecay] )            Jet_hist2[nJet] += int(pow(10., -1 + 5)); 
        if ( flags2P[TrackCategories::LambdaDecay] )        Jet_hist2[nJet] += int(pow(10., -1 + 6)); 
        if ( flags2P[TrackCategories::HadronicProcess] )    Jet_hist2[nJet] += int(pow(10., -1 + 7)); 
        if ( flags2P[TrackCategories::Fake] )	            Jet_hist2[nJet] += int(pow(10., -1 + 8)); 
        if ( flags2P[TrackCategories::SharedInnerHits] )    Jet_hist2[nJet] += int(pow(10., -1 + 9)); 
      }
      if ( idSize > 2 ) {
        if ( flags3P[TrackCategories::BWeakDecay] )         Jet_hist3[nJet] += int(pow(10., -1 + 1)); 
        if ( flags3P[TrackCategories::CWeakDecay] )         Jet_hist3[nJet] += int(pow(10., -1 + 2)); 
        if ( flags3P[TrackCategories::TauDecay] )           Jet_hist3[nJet] += int(pow(10., -1 + 3)); 
        if ( flags3P[TrackCategories::ConversionsProcess] ) Jet_hist3[nJet] += int(pow(10., -1 + 4)); 
        if ( flags3P[TrackCategories::KsDecay] )            Jet_hist3[nJet] += int(pow(10., -1 + 5)); 
        if ( flags3P[TrackCategories::LambdaDecay] )        Jet_hist3[nJet] += int(pow(10., -1 + 6)); 
        if ( flags3P[TrackCategories::HadronicProcess] )    Jet_hist3[nJet] += int(pow(10., -1 + 7)); 
        if ( flags3P[TrackCategories::Fake] )	            Jet_hist3[nJet] += int(pow(10., -1 + 8)); 
        if ( flags3P[TrackCategories::SharedInnerHits] )    Jet_hist3[nJet] += int(pow(10., -1 + 9)); 
      }
      if ( flags1N[TrackCategories::BWeakDecay] )         Jet_hist1[nJet] += 2*int(pow(10., -1 + 1)); 
      if ( flags1N[TrackCategories::CWeakDecay] )         Jet_hist1[nJet] += 2*int(pow(10., -1 + 2)); 
      if ( flags1N[TrackCategories::TauDecay] )           Jet_hist1[nJet] += 2*int(pow(10., -1 + 3)); 
      if ( flags1N[TrackCategories::ConversionsProcess] ) Jet_hist1[nJet] += 2*int(pow(10., -1 + 4)); 
      if ( flags1N[TrackCategories::KsDecay] )            Jet_hist1[nJet] += 2*int(pow(10., -1 + 5)); 
      if ( flags1N[TrackCategories::LambdaDecay] )        Jet_hist1[nJet] += 2*int(pow(10., -1 + 6)); 
      if ( flags1N[TrackCategories::HadronicProcess] )    Jet_hist1[nJet] += 2*int(pow(10., -1 + 7)); 
      if ( flags1N[TrackCategories::Fake] )	          Jet_hist1[nJet] += 2*int(pow(10., -1 + 8)); 
      if ( flags1N[TrackCategories::SharedInnerHits] )	  Jet_hist1[nJet] += 2*int(pow(10., -1 + 9)); 
      if ( idSize > 1 ) {
        if ( flags2N[TrackCategories::BWeakDecay] )         Jet_hist2[nJet] += 2*int(pow(10., -1 + 1)); 
        if ( flags2N[TrackCategories::CWeakDecay] )         Jet_hist2[nJet] += 2*int(pow(10., -1 + 2)); 
        if ( flags2N[TrackCategories::TauDecay] )           Jet_hist2[nJet] += 2*int(pow(10., -1 + 3)); 
        if ( flags2N[TrackCategories::ConversionsProcess] ) Jet_hist2[nJet] += 2*int(pow(10., -1 + 4)); 
        if ( flags2N[TrackCategories::KsDecay] )            Jet_hist2[nJet] += 2*int(pow(10., -1 + 5)); 
        if ( flags2N[TrackCategories::LambdaDecay] )        Jet_hist2[nJet] += 2*int(pow(10., -1 + 6)); 
        if ( flags2N[TrackCategories::HadronicProcess] )    Jet_hist2[nJet] += 2*int(pow(10., -1 + 7)); 
        if ( flags2N[TrackCategories::Fake] )	            Jet_hist2[nJet] += 2*int(pow(10., -1 + 8)); 
        if ( flags2N[TrackCategories::SharedInnerHits] )    Jet_hist2[nJet] += 2*int(pow(10., -1 + 9)); 
      }
      if ( idSize > 2 ) {
        if ( flags3N[TrackCategories::BWeakDecay] )         Jet_hist3[nJet] += 2*int(pow(10., -1 + 1)); 
        if ( flags3N[TrackCategories::CWeakDecay] )         Jet_hist3[nJet] += 2*int(pow(10., -1 + 2)); 
        if ( flags3N[TrackCategories::TauDecay] )           Jet_hist3[nJet] += 2*int(pow(10., -1 + 3)); 
        if ( flags3N[TrackCategories::ConversionsProcess] ) Jet_hist3[nJet] += 2*int(pow(10., -1 + 4)); 
        if ( flags3N[TrackCategories::KsDecay] )            Jet_hist3[nJet] += 2*int(pow(10., -1 + 5)); 
        if ( flags3N[TrackCategories::LambdaDecay] )        Jet_hist3[nJet] += 2*int(pow(10., -1 + 6)); 
        if ( flags3N[TrackCategories::HadronicProcess] )    Jet_hist3[nJet] += 2*int(pow(10., -1 + 7)); 
        if ( flags3N[TrackCategories::Fake] )	            Jet_hist3[nJet] += 2*int(pow(10., -1 + 8)); 
        if ( flags3N[TrackCategories::SharedInnerHits] )    Jet_hist3[nJet] += 2*int(pow(10., -1 + 9)); 
      }
    }
    
    //*****************************************************************
    //get track histories of tracks in jets (for Jet Proba) 
    //*****************************************************************
    Jet_histJet[nJet] = 0;
    
    if(useTrackHistory_ && isData_!=1){
      ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_PosJP);
      TrackRefVector jetProbTracks( (*tagInfo)[ith_tagged].selectedTracks() );  
      
      cap0=0; cap1=0; cap2=0; cap3=0; cap4=0; cap5=0; cap6=0; cap7=0; cap8=0; 
      can0=0; can1=0; can2=0; can3=0; can4=0; can5=0; can6=0; can7=0; can8=0; 

      for (unsigned int i=0; i<jetProbTracks.size(); i++) {
	reco::TrackIPTagInfo::TrackIPData ip = ((*tagInfo)[ith_tagged].impactParameterData())[i];

	if ( ip.ip3d.significance() > 0 ) {
	  TrackCategories::Flags theFlag = classifier_.evaluate( jetProbTracks[i] ).flags();
	  if ( theFlag[TrackCategories::BWeakDecay] )	      cap0 = 1; 
	  if ( theFlag[TrackCategories::CWeakDecay] )	      cap1 = 1; 
	  if ( theFlag[TrackCategories::TauDecay] )	      cap2 = 1;  
	  if ( theFlag[TrackCategories::ConversionsProcess] ) cap3 = 1; 
	  if ( theFlag[TrackCategories::KsDecay] )	      cap4 = 1; 
	  if ( theFlag[TrackCategories::LambdaDecay] )        cap5 = 1; 
	  if ( theFlag[TrackCategories::HadronicProcess] )    cap6 = 1; 
	  if ( theFlag[TrackCategories::Fake] ) 	      cap7 = 1; 
	  if ( theFlag[TrackCategories::SharedInnerHits] )    cap8 = 1; 
	}
	else {
	  TrackCategories::Flags theFlag = classifier_.evaluate( jetProbTracks[i] ).flags();
	  if ( theFlag[TrackCategories::BWeakDecay] )	      can0 = 2; 
	  if ( theFlag[TrackCategories::CWeakDecay] )	      can1 = 2; 
	  if ( theFlag[TrackCategories::TauDecay] )	      can2 = 2;  
	  if ( theFlag[TrackCategories::ConversionsProcess] ) can3 = 2; 
	  if ( theFlag[TrackCategories::KsDecay] )	      can4 = 2; 
	  if ( theFlag[TrackCategories::LambdaDecay] )        can5 = 2; 
	  if ( theFlag[TrackCategories::HadronicProcess] )    can6 = 2; 
	  if ( theFlag[TrackCategories::Fake] ) 	      can7 = 2; 
	  if ( theFlag[TrackCategories::SharedInnerHits] )    can8 = 2; 
	}
      }
      Jet_histJet[nJet] =  cap0+can0 + (cap1+can1)*10 + (cap2+can2)*100 
                        + (cap3+can3)*1000     + (cap4+can4)*10000 
                        + (cap5+can5)*100000   + (cap6+can6)*1000000 
		        + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }
    
    //*****************************************************************
    //get track histories  associated to sec. vertex (for simple SV)
    //*****************************************************************
    Jet_histSvx[nJet] = 0;
    Jet_nFirstSV[nJet]  = nSV;
    
    ith_tagged =    this->TaggedJet(jetsColl.at(ijet),jetTags_SvtxHighEff);
    
    if ( produceJetProbaTree_  && (*tagInfoSVx)[ith_tagged].nVertices() > 0) {
      SV_x[nSV]    = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).x();
      SV_y[nSV]    = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).y();
      SV_z[nSV]    = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).z();
      SV_ex[nSV]   = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).xError();
      SV_ey[nSV]   = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).yError();
      SV_ez[nSV]   = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).zError();
      SV_chi2[nSV] = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).chi2();
      SV_ndf[nSV]  = (*tagInfoSVx)[ith_tagged].secondaryVertex(0).ndof();
      SV_flight[nSV]	= (*tagInfoSVx)[ith_tagged].flightDistance(0).value();
      SV_flightErr[nSV] = (*tagInfoSVx)[ith_tagged].flightDistance(0).error();
//$$      setTracksPV( &(*tagInfoSVx)[ith_tagged].secondaryVertex(0), false );
      nSV++;
    }
    Jet_nLastSV[nJet]  = nSV;    
    
    cap0=0; cap1=0; cap2=0; cap3=0; cap4=0; cap5=0; cap6=0; cap7=0; cap8=0; 
    can0=0; can1=0; can2=0; can3=0; can4=0; can5=0; can6=0; can7=0; can8=0; 

    TrackRefVector  svxPostracks( (*tagInfoSVx)[ith_tagged].vertexTracks(0) );

    if ( useTrackHistory_ && !isData_ ) {
     for (unsigned int i=0; i<svxPostracks.size(); i++) {
	TrackCategories::Flags theFlag = classifier_.evaluate( svxPostracks[i] ).flags();
	if ( theFlag[TrackCategories::BWeakDecay] )         cap0 = 1; 
	if ( theFlag[TrackCategories::CWeakDecay] )         cap1 = 1; 
	if ( theFlag[TrackCategories::TauDecay] )           cap2 = 1;  
	if ( theFlag[TrackCategories::ConversionsProcess] ) cap3 = 1; 
	if ( theFlag[TrackCategories::KsDecay] )            cap4 = 1; 
	if ( theFlag[TrackCategories::LambdaDecay] )        cap5 = 1; 
	if ( theFlag[TrackCategories::HadronicProcess] )    cap6 = 1; 
	if ( theFlag[TrackCategories::Fake] )	            cap7 = 1; 
	if ( theFlag[TrackCategories::SharedInnerHits] )    cap8 = 1; 
      }
    }

    ith_tagged  =  this->TaggedJet(jetsColl.at(ijet),jetTags_negSvtxHighEff);
    TrackRefVector svxNegtracks( (*tagInfoNegSVx)[ith_tagged].vertexTracks(0) );

    if ( useTrackHistory_ && !isData_ ) {
      for (unsigned int i=0; i<svxNegtracks.size(); i++) {
	TrackCategories::Flags theFlag = classifier_.evaluate( svxNegtracks[i] ).flags();
	if ( theFlag[TrackCategories::BWeakDecay] )         can0 = 2; 
	if ( theFlag[TrackCategories::CWeakDecay] )         can1 = 2; 
	if ( theFlag[TrackCategories::TauDecay] )           can2 = 2;  
	if ( theFlag[TrackCategories::ConversionsProcess] ) can3 = 2; 
	if ( theFlag[TrackCategories::KsDecay] )            can4 = 2; 
	if ( theFlag[TrackCategories::LambdaDecay] )        can5 = 2; 
	if ( theFlag[TrackCategories::HadronicProcess] )    can6 = 2; 
	if ( theFlag[TrackCategories::Fake] )	            can7 = 2; 
	if ( theFlag[TrackCategories::SharedInnerHits] )    can8 = 2; 
      }
      Jet_histSvx[nJet] =  cap0+can0 + (cap1+can1)*10 + (cap2+can2)*100 
                        + (cap3+can3)*1000     + (cap4+can4)*10000 
                        + (cap5+can5)*100000   + (cap6+can6)*1000000 
		        + (cap7+can7)*10000000 + (cap8+can8)*100000000;
    }

    Jet_SvxNTracks[nJet] = svxNegtracks.size();	 
    Jet_SvxTracks[nJet]  = svxPostracks.size();	 
   
    //*****************************************************************
    //get track histories of the muon (SoftMuon tagger)
    //*****************************************************************
    Jet_mu_nHit[nJet]  = -1;
    Jet_mu_chi2[nJet]  = -1;
    Jet_mu_pt[nJet]    = -1;
    Jet_mu_ptrel[nJet] = -1;
    Jet_histMuon[nJet] = 0;

    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_softMu);
    if ( (*tagInos_softmuon)[ith_tagged].leptons() !=0 ) {
      Jet_mu_ptrel[nJet] = calculPtRel( (*(*tagInos_softmuon)[ith_tagged].lepton(0)), jetsColl.at(ijet), JES);
      Jet_mu_nHit[nJet]  = (*tagInos_softmuon)[ith_tagged].lepton(0)->hitPattern().numberOfValidHits();
      Jet_mu_chi2[nJet]  = (*tagInos_softmuon)[ith_tagged].lepton(0)->normalizedChi2()  	      ;
      Jet_mu_pt[nJet]    = (*tagInos_softmuon)[ith_tagged].lepton(0)->pt()			      ;
    }
    if ( SoftM > 0 && (*tagInos_softmuon)[ith_tagged].leptons()!=0 ) {
      if ( useTrackHistory_ && !isData_ ) {     
        TrackCategories::Flags theFlagP = classifier_.evaluate( (*tagInos_softmuon)[ith_tagged].lepton(0) ).flags();
        if ( theFlagP[TrackCategories::BWeakDecay] )         Jet_histMuon[nJet] += int(pow(10., -1 + 1)); 
        if ( theFlagP[TrackCategories::CWeakDecay] )         Jet_histMuon[nJet] += int(pow(10., -1 + 2)); 
        if ( theFlagP[TrackCategories::TauDecay] )           Jet_histMuon[nJet] += int(pow(10., -1 + 3));  
        if ( theFlagP[TrackCategories::ConversionsProcess] ) Jet_histMuon[nJet] += int(pow(10., -1 + 4)); 
        if ( theFlagP[TrackCategories::KsDecay] )            Jet_histMuon[nJet] += int(pow(10., -1 + 5)); 
        if ( theFlagP[TrackCategories::LambdaDecay] )        Jet_histMuon[nJet] += int(pow(10., -1 + 6)); 
        if ( theFlagP[TrackCategories::HadronicProcess] )    Jet_histMuon[nJet] += int(pow(10., -1 + 7)); 
        if ( theFlagP[TrackCategories::Fake] )	             Jet_histMuon[nJet] += int(pow(10., -1 + 8)); 
        if ( theFlagP[TrackCategories::SharedInnerHits] )    Jet_histMuon[nJet] += int(pow(10., -1 + 9)); 
      }
    }

    ith_tagged = this->TaggedJet(jetsColl.at(ijet),jetTags_softMuneg);
     //std::cout << "SoftMN " << SoftMN << std::endl;
    if ( SoftMN < 0 && (*tagInos_softmuon)[ith_tagged].leptons()!=0 ) {
      Jet_mu_ptrel[nJet] =  calculPtRel( (*(*tagInos_softmuon)[ith_tagged].lepton(0)), jetsColl.at(ijet), JES);
      Jet_mu_nHit[nJet]  = (*tagInos_softmuon)[ith_tagged].lepton(0)->hitPattern().numberOfValidHits();
      Jet_mu_chi2[nJet]           = (*tagInos_softmuon)[ith_tagged].lepton(0)->normalizedChi2()                ;
      Jet_mu_pt[nJet]             = (*tagInos_softmuon)[ith_tagged].lepton(0)->pt()                            ;
      if ( useTrackHistory_ && !isData_ ) {     
        TrackCategories::Flags theFlagN = classifier_.evaluate( (*tagInos_softmuon)[ith_tagged].lepton(0) ).flags();
        if ( theFlagN[TrackCategories::BWeakDecay] )         Jet_histMuon[nJet] += 2*int(pow(10., -1 + 1)); 
        if ( theFlagN[TrackCategories::CWeakDecay] )         Jet_histMuon[nJet] += 2*int(pow(10., -1 + 2)); 
        if ( theFlagN[TrackCategories::TauDecay] )           Jet_histMuon[nJet] += 2*int(pow(10., -1 + 3));  
        if ( theFlagN[TrackCategories::ConversionsProcess] ) Jet_histMuon[nJet] += 2*int(pow(10., -1 + 4)); 
        if ( theFlagN[TrackCategories::KsDecay] )            Jet_histMuon[nJet] += 2*int(pow(10., -1 + 5)); 
        if ( theFlagN[TrackCategories::LambdaDecay] )        Jet_histMuon[nJet] += 2*int(pow(10., -1 + 6)); 
        if ( theFlagN[TrackCategories::HadronicProcess] )    Jet_histMuon[nJet] += 2*int(pow(10., -1 + 7)); 
        if ( theFlagN[TrackCategories::Fake] )	             Jet_histMuon[nJet] += 2*int(pow(10., -1 + 8)); 
        if ( theFlagN[TrackCategories::SharedInnerHits] )    Jet_histMuon[nJet] += 2*int(pow(10., -1 + 9)); 
      }
    }
        
//$$    
    nJet++;
//$$


    //*********************************
    // Jet analysis
    
    if ( (jetpt * JES) <= minJetPt_ || std::fabs( jeteta ) >= maxJetEta_ ) continue;
    numjet++;
        
    //*****************************************************************
    //define positive and negative tags
    //*****************************************************************
    float varpos = -1000.;
    float varneg = -1000.;
    
    if ( selTagger_ == 0 ) {       // jet proba
      if ( ProbaP > 0 ) varpos = 20.*ProbaP;
      if ( ProbaN > 0 ) varneg = 20.*ProbaN;
      //       if ( ProbaP > tagCut_ ) TagPos = true;
      //       if ( ProbaN > tagCut_ ) TagNeg = true;
    }
    else if ( selTagger_ == 1 ) {     // jet fisrt track
      if ( Jet_Ip1P[nJet] > 0. ) varpos = Jet_Ip1P[nJet];
      if ( Jet_Ip1N[nJet] > 0. ) varneg = Jet_Ip1N[nJet];
    }
    else if ( selTagger_ == 2 ) {     // TC High Eff.
      if ( Jet_Ip2P[nJet] > 0. ) varpos = Jet_Ip2P[nJet];
      if ( Jet_Ip2N[nJet] > 0. ) varneg = Jet_Ip2N[nJet];
    }
    else if ( selTagger_ == 3 ) {     // TC High Pure.
      if ( Jet_Ip3P[nJet] > 0. ) varpos = Jet_Ip3P[nJet];
      if ( Jet_Ip3N[nJet] > 0. ) varneg = Jet_Ip3N[nJet];
    }
    else if ( selTagger_ == 4 ) {    // SSVHE
      if  (Jet_Svx[nJet] > 1)   varpos =  10.*Jet_Svx[nJet] - 10.;
      if  (Jet_SvxN[nJet] < -1) varneg = -10.*Jet_SvxN[nJet] - 10.;
      //       if ( Svtx  > tagCut_ ) TagPos = true;
      //       if (-SvtxN > tagCut_ ) TagNeg = true; 
    }
    else if ( selTagger_ == 5 ) {    // SV combined
      if  (Jet_CombSvx[nJet] > 0)  varpos = 50.*Jet_CombSvx[nJet];
      if  (Jet_CombSvxN[nJet] > 0) varneg = 50.*Jet_CombSvxN[nJet];
      //       if ( CombinedSvtx  > tagCut_ ) TagPos = true;
      //       if (-CombinedSvtxN > tagCut_ ) TagNeg = true;   
    }
    else if ( selTagger_ == 6 ) {    // soft muon ptrel
      if  ( Jet_SoftMu[nJet]  > 0) varpos =  5*Jet_SoftMu[nJet];
      if  ( Jet_SoftMuN[nJet] < 0) varneg = -5*Jet_SoftMuN[nJet];
      //       if  (  SoftM > tagCut_ ) TagPos = true;
      //       if  (  SoftMN> tagCut_ ) TagNeg = true;   
    } 
    else if ( selTagger_ == 7 ) {    // SSVHP
      if  (Jet_SvxHP[nJet] > 1)   varpos =  10.*Jet_SvxHP[nJet] - 10.;
      if  (Jet_SvxNHP[nJet] < -1) varneg = -10.*Jet_SvxNHP[nJet] - 10.;
      //       if ( Svtx  > tagCut_ ) TagPos = true;
      //       if (-SvtxN > tagCut_ ) TagNeg = true; 
    }
    
    //     if ( selTagger_ >= 1 && selTagger_ <= 3) {
    //       if ( varpos > tagCut_ ) TagPos = true;
    //       if ( varneg > tagCut_ ) TagNeg = true;
    //     }
    
    //     // veto on positive tag
    //     if ( Jet_Ip1P[nJet] < vetoPos_ ) Veto = true;
    
    //***************************
    
    int cat1P = 0, cat2P = 0, cat3P = 0, catP = 0;
    int cat1N = 0, cat2N = 0, cat3N = 0, catN = 0;
    
    // Track history
    if ( useTrackHistory_ && indexes.size()!=0 && !isData_ ) {
      flags1P[TrackCategories::ConversionsProcess] ;
           if ( flags1P[TrackCategories::BWeakDecay] )  cat1P = 1; 
      else if ( flags1P[TrackCategories::CWeakDecay] )  cat1P = 2; 
      else if ( flags1P[TrackCategories::TauDecay] )    cat1P = 3; 
      else if ( flags1P[TrackCategories::ConversionsProcess] )  cat1P = 4; 
      else if ( flags1P[TrackCategories::KsDecay] )     cat1P = 5; 
      else if ( flags1P[TrackCategories::LambdaDecay] ) cat1P = 6; 
      else if ( flags1P[TrackCategories::HadronicProcess] ) cat1P = 7; 
      else if ( flags1P[TrackCategories::Fake] )        cat1P = 8; 
      else if ( flags1P[TrackCategories::SharedInnerHits] )         cat1P = 9; 
      if(idSize > 1){
	     if ( flags2P[TrackCategories::BWeakDecay] )  cat2P = 1;
	else if ( flags2P[TrackCategories::CWeakDecay] )  cat2P = 2;
	else if ( flags2P[TrackCategories::TauDecay] )    cat2P = 3;
	else if ( flags2P[TrackCategories::ConversionsProcess] )  cat2P = 4;
	else if ( flags2P[TrackCategories::KsDecay] )     cat2P = 5;
	else if ( flags2P[TrackCategories::LambdaDecay] ) cat2P = 6;
	else if ( flags2P[TrackCategories::HadronicProcess] ) cat2P = 7;
	else if ( flags2P[TrackCategories::Fake] )        cat2P = 8;
	else if ( flags2P[TrackCategories::SharedInnerHits] )         cat2P = 9;
      }
      if(idSize > 2){
	     if ( flags3P[TrackCategories::BWeakDecay] )  cat3P = 1;
	else if ( flags3P[TrackCategories::CWeakDecay] )  cat3P = 2;
	else if ( flags3P[TrackCategories::TauDecay] )    cat3P = 3;
	else if ( flags3P[TrackCategories::ConversionsProcess] )  cat3P = 4;
	else if ( flags3P[TrackCategories::KsDecay] )     cat3P = 5;
	else if ( flags3P[TrackCategories::LambdaDecay] ) cat3P = 6;
	else if ( flags3P[TrackCategories::HadronicProcess] ) cat3P = 7;
	else if ( flags3P[TrackCategories::Fake] )        cat3P = 8;
	else if ( flags3P[TrackCategories::SharedInnerHits] )         cat3P = 9;
      }
           if ( flags1N[TrackCategories::BWeakDecay] )  cat1N = 1;
      else if ( flags1N[TrackCategories::CWeakDecay] )  cat1N = 2;
      else if ( flags1N[TrackCategories::TauDecay] )    cat1N = 3;
      else if ( flags1N[TrackCategories::ConversionsProcess] )  cat1N = 4;
      else if ( flags1N[TrackCategories::KsDecay] )     cat1N = 5;
      else if ( flags1N[TrackCategories::LambdaDecay] ) cat1N = 6;
      else if ( flags1N[TrackCategories::HadronicProcess] ) cat1N = 7;
      else if ( flags1N[TrackCategories::Fake] )        cat1N = 8;
      else if ( flags1N[TrackCategories::SharedInnerHits] )         cat1N = 9;
      if(idSize > 1){
	     if ( flags2N[TrackCategories::BWeakDecay] )  cat2N = 1;
	else if ( flags2N[TrackCategories::CWeakDecay] )  cat2N = 2;
	else if ( flags2N[TrackCategories::TauDecay] )    cat2N = 3;
	else if ( flags2N[TrackCategories::ConversionsProcess] )  cat2N = 4;
	else if ( flags2N[TrackCategories::KsDecay] )     cat2N = 5;
	else if ( flags2N[TrackCategories::LambdaDecay] ) cat2N = 6;
	else if ( flags2N[TrackCategories::HadronicProcess] ) cat2N = 7;
	else if ( flags2N[TrackCategories::Fake] )        cat2N = 8;
	else if ( flags2N[TrackCategories::SharedInnerHits] )         cat2N = 9;
      }
      if(idSize > 2){
	     if ( flags3N[TrackCategories::BWeakDecay] )  cat3N = 1;
	else if ( flags3N[TrackCategories::CWeakDecay] )  cat3N = 2;
	else if ( flags3N[TrackCategories::TauDecay] )    cat3N = 3;
	else if ( flags3N[TrackCategories::ConversionsProcess] )  cat3N = 4;
	else if ( flags3N[TrackCategories::KsDecay] )     cat3N = 5;
	else if ( flags3N[TrackCategories::LambdaDecay] ) cat3N = 6;
	else if ( flags3N[TrackCategories::HadronicProcess] ) cat3N = 7;
	else if ( flags3N[TrackCategories::Fake] )        cat3N = 8;
	else if ( flags3N[TrackCategories::SharedInnerHits] )         cat3N = 9;
      }
      
      if ( selTagger_ == 0 || selTagger_ == 3 ) {
	catP = cat3P;
	catN = cat3N;
      }
      else if ( selTagger_ == 2 ) {
	catP = cat2P;
	catN = cat2N;
      }
      else if ( selTagger_ == 1 ) {
	catP = cat1P;
	catN = cat1N;
      }
    }
    
    //*****************************************************************
    //fill the histograms
    //*****************************************************************

    hData_All_NTracks->Fill( ntagtracks ) ;
    hData_All_JetPt->Fill( ptjet );
    hData_All_JetEta->Fill( etajet );
    
    //*********************************
    // Taggability
    if ( ntagtracks < ntrackMin_ ) continue;

    hData_NTracks->Fill( ntagtracks ) ;
    hData_JetPt->Fill( ptjet );
    hData_JetEta->Fill( etajet );
    
    if ( isData_ != 1 ) {
      hAllFlav_Flavour->Fill( flavour );
    }

    //*********************************
    // Tagging
    
    if ( varneg > 0 ) hData_Tagger->Fill(-varneg );
    if ( varpos > 0 ) hData_Tagger->Fill( varpos );

    if ( Jet_Ip2P[nJet] > 0 )            hData_Tagger_TCHE->Fill( Jet_Ip2P[nJet] );
    if ( Jet_Ip2N[nJet] > 0 )            hData_Tagger_TCHE->Fill(-Jet_Ip2N[nJet] );
    if ( Jet_Ip3P[nJet] > 0 )            hData_Tagger_TCHP->Fill( Jet_Ip3P[nJet] );
    if ( Jet_Ip3N[nJet] > 0 )            hData_Tagger_TCHP->Fill(-Jet_Ip3N[nJet] );
    if ( Jet_ProbaP[nJet] > 0 )          hData_Tagger_JP->Fill( 20.*Jet_ProbaP[nJet] );
    if ( Jet_ProbaN[nJet] > 0 )          hData_Tagger_JP->Fill(-20.*Jet_ProbaN[nJet] );
    if ( Jet_Svx[nJet]  >  1 )           hData_Tagger_SSVHE->Fill( 10.*Jet_Svx[nJet] - 10 );
    if ( Jet_SvxN[nJet] < -1 )           hData_Tagger_SSVHE->Fill( 10*Jet_SvxN[nJet] + 10 );
    if ( Jet_SvxHP[nJet]  >  1 )         hData_Tagger_SSVHP->Fill( 10.*Jet_SvxHP[nJet] - 10 );
    if ( Jet_SvxNHP[nJet] < -1 )         hData_Tagger_SSVHP->Fill( 10*Jet_SvxNHP[nJet] + 10 );
    if ( Jet_CombSvx[nJet] > 0  )        hData_Tagger_CSV->Fill( 50.*Jet_CombSvx[nJet] );
    if ( Jet_CombSvxN[nJet] > 0 )        hData_Tagger_CSV->Fill(-50.*Jet_CombSvxN[nJet] );
    if ( Jet_SoftMu[nJet]  > 0        )  hData_Tagger_MU->Fill( 5*Jet_SoftMu[nJet] );
    if ( Jet_SoftMuN[nJet] < 0        )  hData_Tagger_MU->Fill( 5*Jet_SoftMuN[nJet] );

    if ( isData_ != 1 ) {
      if ( varneg > 0 ) hAllFlav_Tagger->Fill(-varneg );
      if ( varpos > 0 ) hAllFlav_Tagger->Fill( varpos );
      
      if ( varneg > 0 ) {
	if      ( catN == 1 ) hAllFlav_Tagger_Bwd->Fill(-varneg );
	else if ( catN == 2 ) hAllFlav_Tagger_Cwd->Fill(-varneg );
	else if ( catN == 3 ) hAllFlav_Tagger_Tau->Fill(-varneg );
	else if ( catN == 4 ) hAllFlav_Tagger_Gam->Fill(-varneg );
	else if ( catN == 5 ) hAllFlav_Tagger_K0s->Fill(-varneg );
	else if ( catN == 6 ) hAllFlav_Tagger_Lam->Fill(-varneg );
	else if ( catN == 7 ) hAllFlav_Tagger_Int->Fill(-varneg );
	else if ( catN == 8 ) hAllFlav_Tagger_Fak->Fill(-varneg );
	else if ( catN == 9 ) hAllFlav_Tagger_Bad->Fill(-varneg );
	else                  hAllFlav_Tagger_Oth->Fill(-varneg );
      }
      if ( varpos > 0 ) {
	if      ( catP == 1 ) hAllFlav_Tagger_Bwd->Fill( varpos );
	else if ( catP == 2 ) hAllFlav_Tagger_Cwd->Fill( varpos );
	else if ( catP == 3 ) hAllFlav_Tagger_Tau->Fill( varpos );
	else if ( catP == 4 ) hAllFlav_Tagger_Gam->Fill( varpos );
	else if ( catP == 5 ) hAllFlav_Tagger_K0s->Fill( varpos );
	else if ( catP == 6 ) hAllFlav_Tagger_Lam->Fill( varpos );
	else if ( catP == 7 ) hAllFlav_Tagger_Int->Fill( varpos );
	else if ( catP == 8 ) hAllFlav_Tagger_Fak->Fill( varpos );
	else if ( catP == 9 ) hAllFlav_Tagger_Bad->Fill( varpos );
	else                  hAllFlav_Tagger_Oth->Fill( varpos );
      }
      
      if ( flavour == 1 || flavour == 21 ) { // light q+gluon
	if ( varneg > 0 ) hLightFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hLightFlav_Tagger->Fill( varpos );
      } // light q+gluon
      
      if ( flavour == 21 ) { // gluon jets
	if ( varneg > 0 ) hGluonFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hGluonFlav_Tagger->Fill( varpos );
      } // gluon jets
      
      else if ( flavour == 1 ) { // uds jets
	if ( varneg > 0 ) hUDSFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hUDSFlav_Tagger->Fill( varpos );
      } // uds jets
      
      else if ( flavour == 4 ) { // c jets
	if ( varneg > 0 ) hCFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hCFlav_Tagger->Fill( varpos );
      } // c jets
      
      else if ( flavour == 5 ) { // b jets
	if ( varneg > 0 ) hBFlav_Tagger->Fill(-varneg );
	if ( varpos > 0 ) hBFlav_Tagger->Fill( varpos );
      } // b jets
    }
    
  }// end loop on jet
  
//
// Fill TTree
//
  smalltree->Fill();

  hData_All_NJets->Fill( Njets );
  hData_NJets->Fill( numjet );
}


float MistagAnalyzer::calculPtRel(reco::Track theMuon, reco::Jet theJet, double JES )
{
  double pmu = TMath::Sqrt( theMuon.px()*theMuon.px() + theMuon.py()*theMuon.py()  + theMuon.pz()*theMuon.pz() );

//$$  double jetpx = theJet.px() * JES + theMuon.px();
//$$  double jetpy = theJet.py() * JES + theMuon.py();
//$$  double jetpz = theJet.pz() * JES + theMuon.pz();

  double jetpx = theJet.px() * JES ;
  double jetpy = theJet.py() * JES ;
  double jetpz = theJet.pz() * JES ;
  double jetp = TMath::Sqrt(jetpx*jetpx + jetpy*jetpy + jetpz*jetpz);

  double ptrel  = ( jetpx * theMuon.px()  + jetpy * theMuon.py() + jetpz * theMuon.pz() ) / jetp;
  ptrel = TMath::Sqrt( pmu * pmu  - ptrel * ptrel );

  return ptrel;
}

void MistagAnalyzer::setTracksPV( const reco::Vertex *pv, bool isPV ) 
{
  for (reco::Vertex::trackRef_iterator itt = (*pv).tracks_begin(); itt != (*pv).tracks_end(); itt++) {
    for (int i=0; i<nTrack; i++) {
      if ( fabs( (&**itt)->pt() - Track_pt[i] )  < 1.e-5 ) {
	if ( isPV ) {
	  Track_PV[i] = nPV + 1;
	  Track_PVweight[i] = (*pv).trackWeight(*itt);
	}
	else {
	  Track_SV[i] = nSV + 1;
	  Track_SVweight[i] = (*pv).trackWeight(*itt);
//   std::cout << " track " << i << " at SV " << Track_SV[i] << std::endl;
	}
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void
MistagAnalyzer::beginJob() {
//$$
//cat 0
  cat0.etaMax = 2.5;
  cat0.etaMin = 0;
  cat0.nHitsMax= 50;
  cat0.nHitsMin= 8;
  cat0.nPixelHitsMax = 1;
  cat0.nPixelHitsMin = 1;
  cat0.pMax= 5000;
  cat0.pMin= 0;
  cat0.chiMin= 0;
  cat0.chiMax= 5;
  cat0.withFirstPixel = 0;
  
  //cat 1
  cat1.etaMax        = 2.5;
  cat1.etaMin        = 0;
  cat1.nHitsMax      = 50;
  cat1.nHitsMin      = 8;
  cat1.nPixelHitsMax = 8;
  cat1.nPixelHitsMin = 2;
  cat1.pMax          = 5000;
  cat1.pMin          = 0;
  cat1.chiMin        = 2.5;
  cat1.chiMax        = 5;
  cat1.withFirstPixel = 0;
  
  //cat 2
  cat2.etaMax        = 0.8;
  cat2.etaMin        = 0;
  cat2.nHitsMax      = 50;
  cat2.nHitsMin      = 8;
  cat2.nPixelHitsMax = 8;
  cat2.nPixelHitsMin = 3;
  cat2.pMax          = 8;
  cat2.pMin          = 0;
  cat2.chiMin        = 0;
  cat2.chiMax        = 2.5;
  cat2.withFirstPixel = 0;

  //cat 3
  cat3.etaMax        = 1.6;
  cat3.etaMin        = 0.8;
  cat3.nHitsMax      = 50;
  cat3.nHitsMin      = 8;
  cat3.nPixelHitsMax = 8;
  cat3.nPixelHitsMin = 3;
  cat3.pMax          = 8;
  cat3.pMin          = 0;
  cat3.chiMin        = 0;
  cat3.chiMax        = 2.5;
  cat3.withFirstPixel = 0;
  
  //cat 4
  cat4.etaMax        = 2.5;
  cat4.etaMin        = 1.6;
  cat4.nHitsMax      = 50;
  cat4.nHitsMin      = 8;
  cat4.nPixelHitsMax = 8;
  cat4.nPixelHitsMin = 3;
  cat4.pMax          = 8;
  cat4.pMin          = 0;
  cat4.chiMin        = 0;
  cat4.chiMax        = 2.5;
  cat4.withFirstPixel = 0;
  
  //cat 5
  cat5.etaMax        = 2.5;
  cat5.etaMin        = 0;
  cat5.nHitsMax      = 50;
  cat5.nHitsMin      = 8;
  cat5.nPixelHitsMax = 8;
  cat5.nPixelHitsMin = 2;
  cat5.pMax          = 8;
  cat5.pMin          = 0;
  cat5.chiMin        = 0;
  cat5.chiMax        = 2.5;
  cat5.withFirstPixel = 0;
 
  //cat 6
  cat6.etaMax        = 0.8;
  cat6.etaMin        = 0;
  cat6.nHitsMax      = 50;
  cat6.nHitsMin      = 8;
  cat6.nPixelHitsMax = 8;
  cat6.nPixelHitsMin = 3;
  cat6.pMax          = 5000;
  cat6.pMin          = 8;
  cat6.chiMin        = 0;
  cat6.chiMax        = 2.5;
  cat6.withFirstPixel = 0;
  
  //cat 7
  cat7.etaMax        = 1.6;
  cat7.etaMin        = 0.8;
  cat7.nHitsMax      = 50;
  cat7.nHitsMin      = 8;
  cat7.nPixelHitsMax = 8;
  cat7.nPixelHitsMin = 3;
  cat7.pMax          = 5000;
  cat7.pMin          = 8;
  cat7.chiMin        = 0;
  cat7.chiMax        = 2.5;
  cat7.withFirstPixel = 0;
  
  //cat 8
  cat8.etaMax        = 2.5;
  cat8.etaMin        = 1.6;
  cat8.nHitsMax      = 50;
  cat8.nHitsMin      = 8;
  cat8.nPixelHitsMax = 8;
  cat8.nPixelHitsMin = 3;
  cat8.pMax          = 5000;
  cat8.pMin          = 8;
  cat8.chiMin        = 0;
  cat8.chiMax        = 2.5;
  cat8.withFirstPixel = 0;
  
  //cat 9
  cat9.etaMax        = 2.5;
  cat9.etaMin        = 0;
  cat9.nHitsMax      = 50;
  cat9.nHitsMin      = 8;
  cat9.nPixelHitsMax = 2;
  cat9.nPixelHitsMin = 2;
  cat9.pMax          = 5000;
  cat9.pMin          = 8;
  cat9.chiMin        = 0;
  cat9.chiMax        = 2.5;
  cat9.withFirstPixel = 0;
//$$
}


// ------------ method called once each job just after ending the event loop  ------------
void MistagAnalyzer::endJob() {
}


reco::JetFlavour MistagAnalyzer::getMatchedParton(const reco::Jet &jet)
{
  
  reco::JetFlavour jetFlavour;
  if (flavourMatchOptionf == "fastMC") {
    
    //jetFlavour.underlyingParton4Vec(jet.p4());
    
  } else if (flavourMatchOptionf == "hepMC") {
    
    //jetFlavour = jetFlavourIdentifier_.identifyBasedOnPartons(jet);
    
  } else if (flavourMatchOptionf == "genParticle") {
    for ( JetFlavourMatchingCollection::const_iterator j  = theJetPartonMapf->begin();
	  j != theJetPartonMapf->end();
	  j ++ ) {
      RefToBase<Jet> aJet  = (*j).first;       //      const JetFlavour aFlav = (*j).second;
      if ( fabs(aJet->phi() - jet.phi()) < 1.e-5 && fabs(aJet->eta() - jet.eta())< 1.e-5 ){
	// matched
	jetFlavour = reco::JetFlavour (aJet->p4(), math::XYZPoint(0,0,0), (*j).second.getFlavour());
		
	//jetFlavour.flavour((*j).second.getFlavour());
      }
    }
    
    return jetFlavour;
  }
  
  return jetFlavour;
}


//$$
std::vector< float > MistagAnalyzer::getTrackProbabilies(std::vector< float > v, int ipType){
  
  std::vector< float > vectTrackProba;
  
  for (std::vector<float>::const_iterator q = v.begin(); q != v.end(); q++) {
    //positives and negatives tracks
    double p3d = 0;
    if(ipType == 0){
      if ( *q>=0){p3d= (*q)/2.;}else{p3d=1.+(*q)/2.;}
      //if(-log(p3d)> 5) p3d=exp(-5.0);  
      vectTrackProba.push_back(p3d);
    }
    
    //positives tracks only
    if(ipType == 1 && *q >=0 ){
      vectTrackProba.push_back(*q);
    }
    
    //negatives tracks only
    if(ipType == 2 && *q <0){
      vectTrackProba.push_back(-(*q));
    }
  }

  return vectTrackProba;
}


double MistagAnalyzer::calculProbability(std::vector< float > v) {
  
  int ngoodtracks=v.size();
  double SumJet=0.;
  double m_minTrackProb = 0.005;
  for (std::vector<float>::const_iterator q = v.begin(); q != v.end(); q++) {
    SumJet += (*q>m_minTrackProb)?log(*q):log(m_minTrackProb);
  }

  double ProbJet;
  double Loginvlog = 0;
  
  if ( SumJet < 0. ) {
    if ( ngoodtracks >= 2 ) {
      Loginvlog = log(-SumJet);
    }
    double Prob = 1.;
    double lfact = 1.;
    for (int l=1; l!=ngoodtracks; l++) {
      lfact *= l;
      Prob += exp( l*Loginvlog-log(1.*lfact) );
    }
    double LogProb = log(Prob);
    ProbJet = std::min(exp( std::max(LogProb+SumJet,-30.) ), 1.);
  } else {
    ProbJet = 1.;
  }
  return ProbJet;
}


bool MistagAnalyzer::findCat(const reco::Track* track, CategoryFinder& d) {
  //int numcat=-1;
  
  double   p = track->p();
  double eta = track->eta();
  double chi = track->normalizedChi2();
  int   nhit = track->numberOfValidHits();
  int   npix = track->hitPattern().numberOfValidPixelHits();

  bool result = ( p > d.pMin	          &&  p  < d.pMax	      &&
		  fabs(eta) > d.etaMin	  &&  fabs(eta) < d.etaMax    &&
		  nhit >= d.nHitsMin	  &&  nhit <= d.nHitsMax      &&
		  npix >= d.nPixelHitsMin &&  npix <= d.nPixelHitsMax &&
		  chi >= d.chiMin         &&  chi <= d.chiMax );
 
  return result;
}
//$$


int MistagAnalyzer::TaggedJet(reco::Jet calojet, edm::Handle<reco::JetTagCollection > jetTags ) {
  double small = 1.e-5;
  int result = -1;
  
  //std::string moduleLabel = (jetTags).provenance()->moduleLabel();
  //std::string processname = (jetTags).provenance()->processName();
  for (size_t t = 0; t < jetTags->size(); t++) {
    edm::RefToBase<reco::Jet> jet_p = (*jetTags)[t].first;
    if (jet_p.isNull()) {
      continue;
    }
    if (DeltaR<reco::Candidate>()( calojet, *jet_p ) < small) {
      result = (int) t;
    }
    std::string moduleLabel = (jetTags).provenance()->moduleLabel();
    std::string processname = (jetTags).provenance()->processName();
  }
  return result;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MistagAnalyzer);

