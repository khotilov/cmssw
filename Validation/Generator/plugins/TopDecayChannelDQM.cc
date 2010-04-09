#include "FWCore/Utilities/interface/EDMException.h"
#include "Validation/Generator/plugins/TopDecayChannelDQM.h"


/// static const string for status check in  
/// TopCecayChannelDQM::search functions
static const std::string kGenParticles = "genParticles";

// maximal number of daughters 
// to be printed for debugging
static const unsigned int kMAX=5; 


/// constructor
TopDecayChannelDQM::TopDecayChannelDQM(const edm::ParameterSet& cfg):
  log_( cfg.getParameter<unsigned int>( "logEvents"  ) ),
  src_( cfg.getParameter<edm::InputTag>( "src"  ) ), 
  evts_(0)
{
  // register the DQM Service
  store_ = edm::Service<DQMStore>().operator->();
}

/// destructor
TopDecayChannelDQM::~TopDecayChannelDQM()
{
  // free memory
  //delete dqmStore_;
}

/// check the decay chain for different shower types
TopDecayChannelDQM::ShowerType
TopDecayChannelDQM::showerType(const edm::View<reco::GenParticle>& parts) const
{
  for(edm::View<reco::GenParticle>::const_iterator top=parts.begin(); top!=parts.end(); ++top){
    if(top->pdgId()==6 && top->status()==3){
      // check for kHerwig type showers: here the status 3 top quark will 
      // have a single status 2 top quark as daughter, which has again 3 
      // or more status 2 daughters (top3 --> top2 --> X2).
      if( top->numberOfDaughters()==1){
	if( top->begin()->pdgId ()==6 && top->begin()->status()==2 && top->begin()->numberOfDaughters()>1)
	  return kHerwig;
      }
      // check for kPythia type showers: here the status 3 top quark will 
      // have all decay products and a status 2 top quark as daughters; 
      // the status 2 top quark will be w/o further daughters (top3 --> 
      // top2, X2).
      if( top->numberOfDaughters() >1){
	bool containsWBoson=false, containsQuarkDaughter=false;
	for(reco::GenParticle::const_iterator td=top->begin(); td!=top->end(); ++td){
	  if( td->pdgId () < 6 ) containsQuarkDaughter=true;
	  if( td->pdgId ()==24 ) containsWBoson=true;
	}
	if(containsQuarkDaughter && containsWBoson)
	  return kPythia;
      }
    }
  }
  return kNone;
}

/// all that needs to done during the event loop
void
TopDecayChannelDQM::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  // recieve particle listing
  edm::Handle<edm::View<reco::GenParticle> > src; 
  event.getByLabel(src_, src);

  // do the event logging (i.e. print a full listing of all particles) as 
  // long as evts_ is smaller then log_. For log_=0 this kind of logging
  // is switched off.
  if(evts_<log_){ dumpDecayChain(*src); }

  // check every event for corruption; in case throw an exception.
  if(showerType(*src)==kNone){
    throw edm::Exception( edm::errors::LogicError, "particle listing does not correspond to any of the supported listings \n\n");
  }
  bool herwig=(showerType(*src)==kHerwig);

  // major loop in the event; all information about the particle
  // listing that will be monitored is collected here.
  unsigned int iLep=0, iWBoson=0;
  unsigned int iTop=0,iBeauty=0,iElec=0,iMuon=0,iTau=0;
  std::vector<unsigned int> tauDecayType, quarkDecayType;
  for(edm::View<reco::GenParticle>::const_iterator top=src->begin(); top!=src->end(); ++top){
    if( search(top,  6, src_.label()) ){
      ++iTop;
      for(reco::GenParticle::const_iterator td=(!herwig?top->begin():top->begin()->begin()); td!=(!herwig?top->end():top->begin()->end()); ++td){
	if( abs(td->pdgId())<6 ){
	  // keep flavor of the top daughter quark 
	  // for later.
	  quarkDecayType.push_back(abs(td->pdgId()));
	}
	if( search(td,  5, src_.label()) ){
	  ++iBeauty;
	}
	if( search(td, 24, !herwig?src_.label():"test") ){
	  ++iWBoson;
	  for(reco::GenParticle::const_iterator wd=(!herwig?td->begin():td->begin()->begin()->begin()); wd!=(!herwig?td->end():td->begin()->begin()->end()); ++wd){
	    if( abs(wd->pdgId())==11 ){
	      ++iElec;
	    }
	    if( abs(wd->pdgId())==13 ){
	      ++iMuon;
	    }
	    if( abs(wd->pdgId())==15  ){ 
	      // count as iTau if it is leptonic, one-prong
	      // or three-prong and ignore increasing iLep
	      // though else.
	      tauDecayType.push_back(tauDecay(*wd));
	      ++iTau; 
	    }
	  }
	}
      }
    }
  }
  if(evts_<log_){
    edm::LogWarning log("TopDecayChannelDQM");
    log << "----------------------" << "\n"
	<< " iTop    : " << iTop    << "\n"
	<< " iBeauty : " << iBeauty << "\n"
	<< " iElec   : " << iElec   << "\n"
	<< " iMuon   : " << iMuon   << "\n"
	<< " iTau    : " << iTau    << "\n"
	<< "- - - - - - - - - - - " << "\n";
    ++evts_;
  }
  iLep=iElec+iMuon+iTau;

  /* 
   ---
   do the histogram filling
   ---
  */

  bool error=false;
  // *TopDecayChannel*: fill hadronic=0, single leptonic=1, dileptonic=2
  // according to the number of leptons found in the event. If more than
  // two leptons are found in the event the error bin is filled.
  if(iLep<3){
    hists_.find("TopDecayChannel")->second->Fill(iLep);
  }
  else{
    hists_.find("TopDecayChannel")->second->Fill( 3. );
  }
  // *TopDecayWBosons*: fill the number of W bosons found in the decay 
  // chain; no other but the bin for 2 W Bosons should be filled here.
  hists_.find("TopDecayWBosons")->second->Fill(iWBoson);
  // *SemiLeptonType*: fill the lepton flavors for events that decay 
  // single leptonically. If the event is classified to be single lep-
  // tonic but no lepton type found the error bin is filled.
  if(iLep == 1){
    if     ( iElec>0 ) hists_["SemiLeptonType"]->Fill( 0. );
    else if( iMuon>0 ) hists_["SemiLeptonType"]->Fill( 1. );
    else if( iTau >0 ) hists_["SemiLeptonType"]->Fill( 2. );
    else{
      edm::LogWarning("TopDecayWBosons") 
	<< "\n error in checking single leptonic decay type:" 
	<< "\n iLep  == " << iLep  
	<< "\n iElec == " << iElec 
	<< "\n iMuon == " << iMuon 
	<< "\n iTau  == " << iTau  << "\n";
      hists_["SemiLeptonType"]->Fill( 3. );
      error=true;
    }
  }
  // *FullLeptonType*: fill lepton flavors for events that decay di-
  // leptonically. If the event is classified to be di-leptonic but 
  // none of the given combinations is found the error bin is filled.
  if(iLep == 2){
    if     ( iElec==1 && iMuon==1 ) hists_["FullLeptonType"]->Fill( 0. );
    else if( iElec==1 && iTau ==1 ) hists_["FullLeptonType"]->Fill( 1. );
    else if( iMuon==1 && iTau ==1 ) hists_["FullLeptonType"]->Fill( 2. );
    else if( iElec==2             ) hists_["FullLeptonType"]->Fill( 3. );
    else if( iMuon==2             ) hists_["FullLeptonType"]->Fill( 4. );
    else if( iTau ==2             ) hists_["FullLeptonType"]->Fill( 5. );
    else{
      edm::LogWarning("FullLeptonType") 
	<< "\n error in checking di leptonic decay type:" 
	<< "\n iLep  == " << iLep  
	<< "\n iElec == " << iElec 
	<< "\n iMuon == " << iMuon 
	<< "\n iTau  == " << iTau  << "\n";
      hists_["FullLeptonType"]->Fill( 6. );
      error=true;
    }
  }
  // *TauDecayMode*: fill decay mode of tau lepton for any events with 
  // a tau lepton in the decay chain. One-prong, three-prong, leptonic
  // (elec) and leptonic(muon) are checked. If other decay modes than
  // those are found the error bin is filled.
  if(iTau >  0){
    for(unsigned int idx=0; idx<tauDecayType.size(); ++idx){
      if     ( tauDecayType[idx]== 1) hists_["TauDecayMode"]->Fill( 0. );
      else if( tauDecayType[idx]== 3) hists_["TauDecayMode"]->Fill( 1. );
      else if( tauDecayType[idx]==11) hists_["TauDecayMode"]->Fill( 2. );
      else if( tauDecayType[idx]==13) hists_["TauDecayMode"]->Fill( 3. );
      else{
	edm::LogWarning("TauDecayMode") 
	  << "\n un-monitored decay mode found in checking tau decay mode:" 
	  << "\n tauDecayType  == " << tauDecayType[idx]
	  << "\n monitored: 3:3-prong ; 1:1-prong ; 11:electron ; 13:muon \n";
	hists_["TauDecayMode"]->Fill( 4. );
	error=true;
      }
    }
  }
  // *TopDecayQuark*: fill the flavor of the direct top daughter quark. The 
  // first bins are filled if at least one b, s, d quark was found in the 
  // decay chain. As Vtb may be smaller 1 also decays in the other down type
  // quarks are allowed. It may also happen that other quark flavors are 
  // found, when being produced in pairs. If no down type quark is found as 
  // a direct daughter of the top quark the error bin is filled.
  bool containsBeauty=false, containsStrange=false, containsDown=false;
  for(unsigned int idx=0; idx<quarkDecayType.size(); ++idx){
    if     ( quarkDecayType[idx]==5){ containsBeauty =true, hists_["TopDecayQuark"]->Fill( 0. );}
    else if( quarkDecayType[idx]==3){ containsStrange=true, hists_["TopDecayQuark"]->Fill( 1. );}
    else if( quarkDecayType[idx]==1){ containsDown   =true, hists_["TopDecayQuark"]->Fill( 2. );}
    else{
      edm::LogWarning("TopDecayQuark") 
	<< "\n checking top decay associated quark:" << "\n quarkDecayType  == " << quarkDecayType[idx];
      error=true;
    }
  }
  if( !(containsBeauty || containsStrange || containsDown) ){
    hists_["TopDecayQuark"]->Fill( 3. );
  }
  if(error){ dumpDecayChain(*src); };
}

/// count the number of charged particles for tau decays
unsigned int 
TopDecayChannelDQM::countProngs(const reco::Candidate& part) const
{
  // if stable, return 1 or 0
  if(part.status()==1){
    return (part.charge()!=0);
  }
  // if unstable, call recursively on daughters
  int prong =0;
  for(reco::Candidate::const_iterator daughter=part.begin();daughter!=part.end(); ++daughter){
    prong += countProngs(*daughter);
  }
  return prong;
}

/// check tau decay to be leptonic, 1-prong or 3-prong
int
TopDecayChannelDQM::tauDecay(const reco::Candidate& tau) const
{
  bool leptonic = false;
  unsigned int nch = 0;
  unsigned int leptonId = 0;
  // loop on tau decays, check for an elec
  // or muon and count charged particles
  for(reco::Candidate::const_iterator daughter=tau.begin();daughter!=tau.end(); ++daughter){
    // if the tau daughter is again a tau, this means that the particle has 
    // still to be propagated; in that case, return the result of the same 
    // method applied on the on that daughter.
    if(daughter->pdgId()==tau.pdgId()){
      return tauDecay(*daughter);
    }
    // check for leptons
    leptonic |= (abs(daughter->pdgId())==11 || abs(daughter->pdgId())==13);
    if(leptonic && leptonId==0){
      leptonId= abs(daughter->pdgId());
    }
    // count charged particles
    nch += countProngs(*daughter);
  }
  return leptonic ? leptonId : (int)nch;
}

/// search for particle with pdgId (for top)
bool
TopDecayChannelDQM::search(edm::View<reco::GenParticle>::const_iterator& part, int pdgId, const std::string& inputType) const
{
  if(inputType==kGenParticles){
    return (abs(part->pdgId())==pdgId && part->status()==3) ? true : false;
  }
  else{
    return (abs(part->pdgId())==pdgId) ? true : false;
  }
}

/// search for particle with pdgId (overloaded for top daughters)
bool
TopDecayChannelDQM::search(reco::GenParticle::const_iterator& part, int pdgId, const std::string& inputType) const
{
  if(inputType==kGenParticles){
    return (abs(part->pdgId())==pdgId && part->status()==3) ? true : false;
  }
  else{
    return (abs(part->pdgId())==pdgId) ? true : false;
  }
}

void 
TopDecayChannelDQM::dumpDecayChain(const edm::View<reco::GenParticle>& src) const
{
  edm::LogWarning log("DumpDecayChain");
  log << "\n   idx   pdg   stat      px          py         pz             mass          daughter pdg's  "
      << "\n===========================================================================================\n";

  unsigned int idx=0;
  for(edm::View<reco::GenParticle>::const_iterator p=src.begin(); p!=src.end(); ++p, ++idx){
    // loop the top daughters
    log << std::right << std::setw( 5) << idx
	<< std::right << std::setw( 7) << src[idx].pdgId()
	<< std::right << std::setw( 5) << src[idx].status() << "  "
	<< std::right << std::setw(10) << std::setprecision( 6 ) << src[idx].p4().x() << "  "	
	<< std::right << std::setw(10) << std::setprecision( 6 ) << src[idx].p4().y() << "  "	
	<< std::right << std::setw(10) << std::setprecision( 6 ) << src[idx].p4().z() << "  "	
	<< std::right << std::setw(15) << std::setprecision( 6 ) << src[idx].p4().mass() 
	<< "   ";
    // search for potential daughters; if they exits 
    // print the daughter to the screen in the last 
    // column of the table separated by ','
    TString pdgIds;
    unsigned int jdx=0;
    for(reco::GenParticle::const_iterator d=p->begin(); d!=p->end(); ++d, ++jdx){
      if(jdx<kMAX){
	pdgIds+=d->pdgId();
	if(d+1 != p->end()){
	  pdgIds+= ",";
	}
      }
      else{
	pdgIds+="...(";
	pdgIds+= p->numberOfDaughters();
	pdgIds+=")";
	break;
      }
    }
    if(idx>0){
      log << std::setfill( ' ' ) << std::right << std::setw(15) << pdgIds; 
      log << "\n";
    }
    else{
      log << std::setfill( ' ' ) << std::right << std::setw(15) << "-\n";
    }
  }
}

/// all that needs to be done at the beginning of a run
void 
TopDecayChannelDQM::beginJob()
{
  store_->setCurrentFolder("Physics/Top/TopDecayChannelDQM");

  // top decay channel
  hists_["TopDecayChannel"] = store_->book1D("TopDecayChannel" , "TopDecayChannel" , 4, 0., 4.);
  hists_["TopDecayChannel"]->getTH1()->SetOption("TEXT");
  hists_["TopDecayChannel"]->getTH1()->GetXaxis()->SetBinLabel(1 , "Full Hadronic");
  hists_["TopDecayChannel"]->getTH1()->GetXaxis()->SetBinLabel(2 , "Single Lepton");
  hists_["TopDecayChannel"]->getTH1()->GetXaxis()->SetBinLabel(3 , "Di Lepton"    );
  hists_["TopDecayChannel"]->getTH1()->GetXaxis()->SetBinLabel(4 , "Other/Error"  );
  // number of W bosons found in the decay chain
  hists_["TopDecayWBosons"] = store_->book1D("TopDecayWBosons" , "TopDecayWBosons" , 3, 0., 3.);
  hists_["TopDecayWBosons" ]->getTH1()->SetOption("TEXT");
  hists_["TopDecayWBosons" ]->getTH1()->GetXaxis()->SetBinLabel(1 , "0 WBoson(s)" );
  hists_["TopDecayWBosons" ]->getTH1()->GetXaxis()->SetBinLabel(2 , "1 WBoson(s)" );
  hists_["TopDecayWBosons" ]->getTH1()->GetXaxis()->SetBinLabel(3 , "2 WBoson(s)" );
  // lepton flavours in the semi-leptonic decay
  hists_["SemiLeptonType" ] = store_->book1D("SemiLeptonType"  , "SemiLeptonType"  , 4, 0., 4.);
  hists_["SemiLeptonType" ]->getTH1()->SetOption("TEXT");
  hists_["SemiLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(1 , "Electron"     );
  hists_["SemiLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(2 , "Muon"         );
  hists_["SemiLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(3 , "Tau"          );
  hists_["SemiLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(4 , "Other/Error"  );
  // lepton flavours in the full-leptonic decay
  hists_["FullLeptonType" ] = store_->book1D("FullLeptonType"  , "FullLeptonType"  , 7, 0., 7.);
  hists_["FullLeptonType" ]->getTH1()->SetOption("TEXT");
  hists_["FullLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(1 , "Elec-Muon"    );
  hists_["FullLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(2 , "Elec-Tau "    );
  hists_["FullLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(3 , "Muon-Tau "    );
  hists_["FullLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(4 , "Elec-Elec"    );
  hists_["FullLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(5 , "Muon-Muon"    );
  hists_["FullLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(6 , "Tau -Tau "    );
  hists_["FullLeptonType" ]->getTH1()->GetXaxis()->SetBinLabel(7 , "Other/Error"  );
  // prongs in subsequent leptonic/hadronic tau decays 
  hists_["TauDecayMode"   ] = store_->book1D("TauDecayMode"   , "TauDecayMode"    , 5, 0., 5.);
  hists_["TauDecayMode"   ]->getTH1()->SetOption("TEXT");
  hists_["TauDecayMode"   ]->getTH1()->GetXaxis()->SetBinLabel(1 , "One-Prong"    );
  hists_["TauDecayMode"   ]->getTH1()->GetXaxis()->SetBinLabel(2 , "Three-Prong"  );
  hists_["TauDecayMode"   ]->getTH1()->GetXaxis()->SetBinLabel(3 , "Lepton(Elec)" );
  hists_["TauDecayMode"   ]->getTH1()->GetXaxis()->SetBinLabel(4 , "Lepton(Muon)" );
  hists_["TauDecayMode"   ]->getTH1()->GetXaxis()->SetBinLabel(5 , "Other"        );
  // flavors of quarks originating from the top quark
  hists_["TopDecayQuark"  ] = store_->book1D("TopDecayQuark"   , "TopDecayQuark"   , 4, 0., 4.);
  hists_["TopDecayQuark"  ]->getTH1()->SetOption("TEXT");
  hists_["TopDecayQuark"  ]->getTH1()->GetXaxis()->SetBinLabel(1 , "b-Quark"      );
  hists_["TopDecayQuark"  ]->getTH1()->GetXaxis()->SetBinLabel(2 , "s-Quark"      );
  hists_["TopDecayQuark"  ]->getTH1()->GetXaxis()->SetBinLabel(3 , "d-Quark"      );
  hists_["TopDecayQuark"  ]->getTH1()->GetXaxis()->SetBinLabel(4 , "Other/Error"  );
}

/// all that needs to be done at the end of a run 
void 
TopDecayChannelDQM::endJob() 
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopDecayChannelDQM);
