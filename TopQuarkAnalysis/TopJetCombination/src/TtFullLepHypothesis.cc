#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "TopQuarkAnalysis/TopJetCombination/interface/TtFullLepHypothesis.h"

/// default constructor
TtFullLepHypothesis::TtFullLepHypothesis(const edm::ParameterSet& cfg):
  elecs_(cfg.getParameter<edm::InputTag>("electrons")),
  mus_  (cfg.getParameter<edm::InputTag>("muons")),
  jets_ (cfg.getParameter<edm::InputTag>("jets")),
  mets_ (cfg.getParameter<edm::InputTag>("mets")),

  lepton_(0), leptonBar_(0), b_(0), 
  bBar_(0), neutrino_(0), neutrinoBar_(0)
{  
  getMatch_ = false;
  if( cfg.exists("match") ) {
    getMatch_ = true;
    match_ = cfg.getParameter<edm::InputTag>("match");
  }
  if( cfg.exists("jetCorrectionLevel") ) {
    jetCorrectionLevel_ = cfg.getParameter<std::string>("jetCorrectionLevel");
  }
  produces<std::vector<std::pair<reco::CompositeCandidate, std::vector<int> > > >();
  produces<int>("Key");
}

/// default destructor
TtFullLepHypothesis::~TtFullLepHypothesis()
{
  if( lepton_      ) delete lepton_;
  if( leptonBar_   ) delete leptonBar_;
  if( b_           ) delete b_;
  if( bBar_        ) delete bBar_;
  if( neutrino_    ) delete neutrino_;
  if( neutrinoBar_ ) delete neutrinoBar_;
  //if( met_         ) delete met_;
}

/// produce the event hypothesis as CompositeCandidate and Key
void
TtFullLepHypothesis::produce(edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<std::vector<pat::Electron> > elecs;
  evt.getByLabel(elecs_, elecs);
  
  edm::Handle<std::vector<pat::Muon> > mus;
  evt.getByLabel(mus_, mus);  

  edm::Handle<std::vector<pat::Jet> > jets;
  evt.getByLabel(jets_, jets);

  edm::Handle<std::vector<pat::MET> > mets;
  evt.getByLabel(mets_, mets);

  std::vector<std::vector<int> > matchVec;
  if( getMatch_ ) {
    edm::Handle<std::vector<std::vector<int> > > matchHandle;
    evt.getByLabel(match_, matchHandle);;  
    matchVec = *matchHandle;
  }
  else {
    std::vector<int> dummyMatch;
    for(unsigned int i = 0; i < 4; ++i) 
      dummyMatch.push_back( -1 );
    matchVec.push_back( dummyMatch );
  }

  // declare auto_ptr for products
  std::auto_ptr<std::vector<std::pair<reco::CompositeCandidate, std::vector<int> > > >
    pOut( new std::vector<std::pair<reco::CompositeCandidate, std::vector<int> > > );
  std::auto_ptr<int> pKey(new int);
  
  // build and feed out key
  buildKey();  
  *pKey=key();
  evt.put(pKey, "Key");

  // go through given vector of jet combinations
  unsigned int idMatch = 0;  
  typedef std::vector<std::vector<int> >::iterator MatchVecIterator;
  for(MatchVecIterator match = matchVec.begin(); match != matchVec.end(); ++match) {        
    // reset pointers
    resetCandidates();
    // build hypothesis
    buildHypo(evt, elecs, mus, jets, mets, *match, idMatch++);        
    pOut->push_back( std::make_pair(hypo(), *match) );    
  }
  // feed out hyps and matches
  evt.put(pOut);
}

/// reset candidate pointers before hypo build process
void
TtFullLepHypothesis::resetCandidates()
{
  lepton_     = 0;
  leptonBar_  = 0;
  b_          = 0;
  bBar_       = 0;
  neutrino_   = 0;
  neutrinoBar_= 0;
  //met_        = 0;
}

/// return event hypothesis
reco::CompositeCandidate
TtFullLepHypothesis::hypo()
{ 
  // check for sanity of the hypothesis  
  if( !lepton_ || !leptonBar_ || !b_ || !bBar_ ){
    return reco::CompositeCandidate();
  } 
       
  if( key()==TtFullLeptonicEvent::kGenMatch && (!recNu || !recNuBar) ){
    edm::LogWarning("TtFullHypothesis") << "no neutrinos for gen match" << std::endl; 
    return reco::CompositeCandidate();
  }  
  if( key()==TtFullLeptonicEvent::kKinSolution && (!neutrino_ || !neutrinoBar_) ){
    edm::LogWarning("TtFullHypothesis") << "no neutrinos for kin solution" << std::endl;  
    return reco::CompositeCandidate();    
  }
          
  // setup transient references
  reco::CompositeCandidate hyp, Top, WPlus, TopBar, WMinus;

  AddFourMomenta addFourMomenta; 
     
  // build up the top branch
  WPlus.addDaughter(*leptonBar_, TtFullLepDaughter::LepBar);
  if(key()==TtFullLeptonicEvent::kKinSolution)
    WPlus.addDaughter(*neutrino_, TtFullLepDaughter::Nu);  
  else if(key()==TtFullLeptonicEvent::kGenMatch)
    WPlus.addDaughter(*recNu, TtFullLepDaughter::Nu);     
  addFourMomenta.set(WPlus);
  Top.addDaughter(WPlus, TtFullLepDaughter::WPlus);
  Top.addDaughter(*b_,TtFullLepDaughter::B);
  addFourMomenta.set(Top);   
   
  // build up the anti top branch
  WMinus.addDaughter(*lepton_, TtFullLepDaughter::Lep); 
  if(key()==TtFullLeptonicEvent::kKinSolution)
    WMinus.addDaughter(*neutrinoBar_, TtFullLepDaughter::NuBar);  
  else if(key()==TtFullLeptonicEvent::kGenMatch)  
    WMinus.addDaughter(*recNuBar, TtFullLepDaughter::NuBar);   
  addFourMomenta.set(WMinus);
  TopBar.addDaughter(WMinus, TtFullLepDaughter::WMinus);
  TopBar.addDaughter(*bBar_, TtFullLepDaughter::BBar);
  addFourMomenta.set(TopBar);   

  // build ttbar hypothesis
  hyp.addDaughter(Top, TtFullLepDaughter::Top);
  hyp.addDaughter(TopBar, TtFullLepDaughter::TopBar);
  addFourMomenta.set( hyp );
     
  // the four momentum of the met is not added to the hypothesis
  // because it is allready included through the neutrinos     
  //hyp.addDaughter(*met_, TtFullLepDaughter::Met);
  return hyp;
}

/// helper function to contruct the proper correction level string for corresponding quarkType, for unknown quarkTypes or if JetCorrectionLevel_ was not configured an emty string is returned 
std::string
TtFullLepHypothesis::jetCorrectionLevel(const std::string& quarkType)
{
  std::string level;
  // jetCorrectionLevel was not configured
  if(jetCorrectionLevel_.empty())
    return level;

  // quarkType is unknown
  if( !(quarkType=="lightQuark" || quarkType=="bQuark") )
    return level;

  // combine correction level; start with a ':' even if 
  // there is no flavor tag to be added, as it is needed
  // by setCandidate to disentangle the correction tag 
  // from a potential flavor tag, which can be empty
  level=jetCorrectionLevel_+":";
  if( level=="had:" || level=="ue:" || level=="part:" ){
    if(quarkType=="lightQuark"){level+="uds";}
    if(quarkType=="bQuark"){level+="b";  }
  }
  return level;
}

/// use one object in a jet collection to set a ShallowClonePtrCandidate with proper jet corrections
void 
TtFullLepHypothesis::setCandidate(const edm::Handle<std::vector<pat::Jet> >& handle, const int& idx, reco::ShallowClonePtrCandidate*& clone, const std::string& correctionLevel)
{
  // disentangle the correction from the potential flavor tag 
  // by the separating ':'; the flavor tag can be empty though
  std::string step   = correctionLevel.substr(0,correctionLevel.find(":"));
  std::string flavor = correctionLevel.substr(1+correctionLevel.find(":"));
  edm::Ptr<pat::Jet> ptr = edm::Ptr<pat::Jet>(handle, idx);
  clone = new reco::ShallowClonePtrCandidate( ptr, ptr->charge(), ptr->correctedJet(step, flavor).p4(), ptr->vertex() );
}
