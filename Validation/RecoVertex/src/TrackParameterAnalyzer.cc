#include "Validation/RecoVertex/interface/TrackParameterAnalyzer.h"
#include <string>
#include <vector>
#include "SimGeneral/HepPDT/interface/HepPDTable.h"
#include "SimGeneral/HepPDT/interface/HepParticleData.h"


//
//
// constants, enums and typedefs
//
const double fBfield=4.06;

//
// static data member definitions
//

//
// constructors and destructor
//
TrackParameterAnalyzer::TrackParameterAnalyzer(const edm::ParameterSet& iConfig) :
  simG4_( iConfig.getParameter<edm::InputTag>( "simG4" ) ) {
   //now do whatever initialization is needed

  recoTrackProducer_   = iConfig.getUntrackedParameter<std::string>("recoTrackProducer");
  // open output file to store histograms}
  outputFile_   = iConfig.getUntrackedParameter<std::string>("outputFile");
  rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE"); 
  verbose_= iConfig.getUntrackedParameter<bool>("verbose", false);
}


TrackParameterAnalyzer::~TrackParameterAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete rootFile_;
}



//
// member functions
//
void TrackParameterAnalyzer::beginJob(edm::EventSetup const&){
  rootFile_->cd();
  h1_pull0_ = new TH1F("pull0","pull kappa",100,-25.,25.);
  h1_pull1_ = new TH1F("pull1","pull theta",100,-25.,25.);
  h1_pull2_ = new TH1F("pull2","pull phi  ",100,-25.,25.);
  h1_pull3_ = new TH1F("pull3","pull dca  ",100,-25.,25.);
  h1_pull4_ = new TH1F("pull4","pull zdca ",100,-25.,25.);
  h1_Beff_  = new TH1F("Beff", "Beff",2000,-10.,10.);
  h2_dvsphi_ = new TH2F("dvsphi","dvsphi",360,-3.14159,3.14159,100,-0.1,0.1);
}


void TrackParameterAnalyzer::endJob() {
  rootFile_->cd();
  h1_pull0_->Write();
  h1_pull1_->Write();
  h1_pull2_->Write();
  h1_pull3_->Write();
  h1_pull4_->Write();
  h1_Beff_->Write();
  h2_dvsphi_->Write();
}

// helper function
bool TrackParameterAnalyzer::match(const PerigeeTrajectoryParameters::ParameterVector  &a, 
				   const PerigeeTrajectoryParameters::ParameterVector  &b){
  double dtheta=a(1)-b(1);
  double dphi  =a(2)-b(2);
  if (dphi>M_PI){ dphi-=M_2_PI; }else if(dphi<-M_PI){dphi+=M_2_PI;}
  return ( (fabs(dtheta)<0.1) && (fabs(dphi)<0.1) );
}

// ------------ method called to produce the data  ------------
void
TrackParameterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  
   Handle<edm::SimVertexContainer> simVtcs;
   iEvent.getByLabel( simG4_, simVtcs);
   if(verbose_){
     std::cout << "SimVertex " << simVtcs->size() << std::endl;
     for(edm::SimVertexContainer::const_iterator v=simVtcs->begin();
	 v!=simVtcs->end(); ++v){
       std::cout << "simvtx "
		 << std::setw(10) << std::setprecision(3)
		 << v->position().x() << " "
		 << v->position().y() << " "
		 << v->position().z() << " "
		 << v->parentIndex() << " "
		 << v->noParent() << " "
		 << std::endl;
     }
   }
   
   // get the simulated tracks, extract perigee parameters
   Handle<SimTrackContainer> simTrks;
   iEvent.getByLabel( simG4_, simTrks);
   std::vector<PerigeeTrajectoryParameters::ParameterVector > tsim;
   for(edm::SimTrackContainer::const_iterator t=simTrks->begin();
       t!=simTrks->end(); ++t){
     if (t->noVertex()){
       std::cout << "simtrk  has no vertex" << std::endl;
       return;
     }else{
       // get the vertex position
       HepLorentzVector v=(*simVtcs)[t->vertIndex()].position();
       int pdgCode=t->type();
       double Q=HepPDT::theTable().getParticleData(pdgCode)->charge();
       HepLorentzVector p=t->momentum();
       if(verbose_){
	 std::cout << ".simtrk "
		   << std::setw(10) << std::setprecision(3)
		   << t->genpartIndex() << " "
		   << t->vertIndex() << " "
		   << t->type() << " "
		   << Q << " " 
		   << p.perp() << " " 
		   << " vx="  << v.x() << " vy=" << v.y() << " vz=" << v.z()   << " " 
		   << std::endl;
       }
       double kappa=-Q*0.002998*fBfield/p.perp();
       double D0=v.x()*sin(p.phi())-v.y()*cos(p.phi())-0.5*kappa*(v.x()*v.x()+v.y()*v.y());
       double q=sqrt(1.-2.*kappa*D0);
       double s0=(v.x()*cos(p.phi())+v.y()*sin(p.phi()))/q;
       double s1;
       if (fabs(kappa*s0)>0.001){
	 s1=asin(kappa*s0)/kappa;
       }else{
	 double ks02=(kappa*s0)*(kappa*s0);
	 s1=s0*(1.+ks02/6.+3./40.*ks02*ks02+5./112.*pow(ks02,3));
       }
       PerigeeTrajectoryParameters::ParameterVector par;
       par[reco::TrackBase::i_transverseCurvature] = kappa;
       par[reco::TrackBase::i_theta] = p.theta();
       par[reco::TrackBase::i_phi0] = p.phi()-asin(kappa*s0);
       par[reco::TrackBase::i_d0] = 2.*D0/(1.+q);
       par[reco::TrackBase::i_dz] = v.z()-s1/tan(p.theta());
       tsim.push_back(par);
     }
   }

   // simtrack parameters are in now tsim
   // loop over tracks and try to match them to simulated tracks


   Handle<reco::TrackCollection> recTracks;
   iEvent.getByLabel(recoTrackProducer_, recTracks);
   for(reco::TrackCollection::const_iterator t=recTracks->begin();
       t!=recTracks->end(); ++t){
     PerigeeTrajectoryParameters::ParameterVector  p = t->parameters();
     PerigeeTrajectoryError::CovarianceMatrix c = t->covariance();
     for(std::vector<PerigeeTrajectoryParameters::ParameterVector>::const_iterator s=tsim.begin();
	 s!=tsim.end(); ++s){
       if (match(*s,p)){
	 //std::cout << "match found" << std::endl;
	 h1_pull0_->Fill((p(0)-(*s)(0))/sqrt(c(0,0)));
	 h1_pull1_->Fill((p(1)-(*s)(1))/sqrt(c(1,1)));
	 h1_pull2_->Fill((p(2)-(*s)(2))/sqrt(c(2,2)));
	 h1_pull3_->Fill((p(3)-(*s)(3))/sqrt(c(3,3)));
	 h1_pull4_->Fill((p(4)-(*s)(4))/sqrt(c(4,4)));
	 h1_Beff_->Fill(p(0)/(*s)(0)*fBfield);
	 h2_dvsphi_->Fill(p(2),p(3));
       }
     }
   }




}

//define this as a plug-in
//DEFINE_FWK_MODULE(TrackParameterAnalyzer)
