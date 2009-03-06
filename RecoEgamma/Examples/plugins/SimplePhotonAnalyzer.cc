/**\class PhotonSimpleAnalyzer
 **
 ** $Date: 2008/12/08 22:33:47 $
 ** $Revision: 1.18 $
 ** \author Nancy Marinelli, U. of Notre Dame, US
*/

#include "RecoEgamma/Examples/plugins/SimplePhotonAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"
//
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
//
#include "CLHEP/Units/PhysicalConstants.h"
#include "TFile.h"

//========================================================================
SimplePhotonAnalyzer::SimplePhotonAnalyzer( const edm::ParameterSet& ps )
//========================================================================
{
  photonCollectionProducer_ = ps.getParameter<std::string>("phoProducer");
  photonCollection_ = ps.getParameter<std::string>("photonCollection");

  barrelEcalHits_   = ps.getParameter<edm::InputTag>("barrelEcalHits");
  endcapEcalHits_   = ps.getParameter<edm::InputTag>("endcapEcalHits");



  mcProducer_ = ps.getParameter<std::string>("mcProducer");
  //mcCollection_ = ps.getParameter<std::string>("mcCollection");
  vertexProducer_ = ps.getParameter<std::string>("primaryVertexProducer");
  sample_ = ps.getParameter<int>("sample");


}


//========================================================================
SimplePhotonAnalyzer::~SimplePhotonAnalyzer()
//========================================================================
{


}

//========================================================================
void
SimplePhotonAnalyzer::beginJob() {
//========================================================================

  edm::Service<TFileService> fs;

  float hiE=0;
  float loE=0;
  float hiEt=0;
  float loEt=0;
  float dPhi=0;
  float loRes=0;
  float hiRes=0;
  if ( sample_ ==1 ) {
    loE=0.;
    hiE=30.;
    loEt=0.;
    hiEt=30.;
    dPhi=0.2;
    loRes=0.;
    hiRes=1.2;
  } else if ( sample_ ==2 ) {
    loE=0.;
    hiE=200.;
    loEt=0.;
    hiEt=50.;
    dPhi=0.05;
    loRes=0.7;
    hiRes=1.2;
  }  else if (sample_==4) {
    loE=0.;
    hiE=6000.;
    loEt=0.;
    hiEt=1200.;
    dPhi=0.05;
    loRes=0.7;
    hiRes=1.2;
  }


  effEta_ = fs->make<TProfile> ("effEta"," Photon reconstruction efficiency",50,-2.5,2.5);
  effPhi_ = fs->make<TProfile> ("effPhi"," Photon reconstruction efficiency",80, -3.14, 3.14);

  h1_scEta_ = fs->make<TH1F>("scEta"," SC Eta ",40,-3., 3.);
  h1_scPhi_ = fs->make<TH1F>("scPhi"," SC Phi ",40,-3.14, 3.14);
  h1_deltaEtaSC_ = fs->make<TH1F>("deltaEtaSC"," SC Eta minus Generated photon Eta  ",100,-0.02, 0.02);
  h1_deltaPhiSC_ = fs->make<TH1F>("deltaPhiSC"," SC Phi minus Generated photon Phi ",100,-dPhi, dPhi);
  h1_deltaEta_ = fs->make<TH1F>("deltaEta"," Reco photon Eta minus Generated photon Eta  ",100,-0.2, 0.2);
  h1_deltaPhi_ = fs->make<TH1F>("deltaPhi","Reco photon Phi minus Generated photon Phi ",100,-dPhi, dPhi);
  h1_pho_Eta_ = fs->make<TH1F>("phoEta","Photon  Eta ",40,-3., 3.);
  h1_pho_Phi_ = fs->make<TH1F>("phoPhi","Photon  Phi ",40,-3.14, 3.14);


  h1_scEt_ = fs->make<TH1F>("scEt"," SC Et ",100,loEt,hiEt);
  h1_scE_ = fs->make<TH1F>("scE"," SC Energy ",100,loE,hiE);
  h1_pho_E_ = fs->make<TH1F>("phoE","Photon Energy ",100,loE,hiE);

  h1_e5x5_unconvBarrel_ = fs->make<TH1F>("e5x5_unconvBarrelOverEtrue"," Photon rec/true energy if R9>0.93 Barrel ",100,loRes, hiRes);
  h1_e5x5_unconvEndcap_ = fs->make<TH1F>("e5x5_unconvEndcapOverEtrue"," Photon rec/true energy if R9>0.93 Endcap ",100,loRes, hiRes);
  h1_ePho_convBarrel_ = fs->make<TH1F>("ePho_convBarrelOverEtrue"," Photon rec/true energy if R9<=0.93 Barrel ",100,loRes, hiRes);
  h1_ePho_convEndcap_ = fs->make<TH1F>("ePho_convEndcapOverEtrue"," Photon rec/true energy if R9<=0.93 Endcap ",100,loRes, hiRes);


 //
  h1_recEoverTrueEBarrel_ = fs->make<TH1F>("recEoverTrueEBarrel"," Reco photon Energy over Generated photon Energy: Barrel ",100,loRes, hiRes);
  h1_recEoverTrueEEndcap_ = fs->make<TH1F>("recEoverTrueEEndcap"," Reco photon Energy over Generated photon Energy: Endcap ",100,loRes, hiRes);
  h1_recESCoverTrueEBarrel_ = fs->make<TH1F>("recESCoverTrueEBarrel"," Reco SC Energy over Generated photon Energy: Barrel ",100,loRes, hiRes);
  h1_recESCoverTrueEEndcap_ = fs->make<TH1F>("recESCoverTrueEEndcap"," Reco SC Energy over Generated photon Energy: Endcap ",100,loRes, hiRes);

  //

  h1_pho_R9Barrel_ = fs->make<TH1F>("phoR9Barrel","Photon  3x3 energy / SuperCluster energy : Barrel ",100,0.,1.2);
  h1_pho_R9Endcap_ = fs->make<TH1F>("phoR9Endcap","Photon  3x3 energy / SuperCluster energy : Endcap ",100,0.,1.2);



}


//========================================================================
void
SimplePhotonAnalyzer::analyze( const edm::Event& evt, const edm::EventSetup& es ) {
//========================================================================

  using namespace edm; // needed for all fwk related classes
  edm::LogInfo("PhotonAnalyzer") << "Analyzing event number: " << evt.id() << "\n";


 // get the  calo topology  from the event setup:
  edm::ESHandle<CaloTopology> pTopology;
  es.get<CaloTopologyRecord>().get(theCaloTopo_);
  const CaloTopology *topology = theCaloTopo_.product();


  // Get the  corrected  photon collection (set in the configuration) which also contains infos about conversions

  Handle<reco::PhotonCollection> photonHandle;
  evt.getByLabel(photonCollectionProducer_, photonCollection_ , photonHandle);
  const reco::PhotonCollection photonCollection = *(photonHandle.product());


  // Get the primary event vertex
  Handle<reco::VertexCollection> vertexHandle;
  evt.getByLabel(vertexProducer_, vertexHandle);
  reco::VertexCollection vertexCollection = *(vertexHandle.product());
  math::XYZPoint vtx(0.,0.,0.);
  if (vertexCollection.size()>0) {
    std::cout << " I am using the Primary vertex position " << std::endl;
    vtx = vertexCollection.begin()->position();
  }

  /// Get the MC truth
  Handle< HepMCProduct > hepProd ;
  evt.getByLabel( mcProducer_.c_str(),  hepProd ) ;
  const HepMC::GenEvent * myGenEvent = hepProd->GetEvent();

  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) {
    if ( !( (*p)->pdg_id() == 22 && (*p)->status()==1 )  )  continue;

    // single primary photons or photons from Higgs or RS Graviton
    HepMC::GenParticle* mother = 0;
    if ( (*p)->production_vertex() )  {
      if ( (*p)->production_vertex()->particles_begin(HepMC::parents) !=
           (*p)->production_vertex()->particles_end(HepMC::parents))
	mother = *((*p)->production_vertex()->particles_begin(HepMC::parents));
    }
    if ( ((mother == 0) || ((mother != 0) && (mother->pdg_id() == 25))
	  || ((mother != 0) && (mother->pdg_id() == 22)))) {

      float minDelta=10000.;
      std::vector<reco::Photon> localPhotons;
      int index=0;
      int iMatch=-1;

      float phiPho=(*p)->momentum().phi();
      float etaPho=(*p)->momentum().eta();
      etaPho = etaTransformation(etaPho, (*p)->production_vertex()->position().z()/10. );



      // loop over corrected  Photon candidates
      for( reco::PhotonCollection::const_iterator  iPho = photonCollection.begin(); iPho != photonCollection.end(); iPho++) {

	/////  Set event vertex
	reco::Photon localPho = reco::Photon(*iPho);
	//	localPho.setVertex(vtx);
	localPhotons.push_back(localPho);

	/// Match reconstructed photon candidates with the nearest generated photonPho;
	float phiClu=localPho.phi();
	float etaClu=localPho.eta();
	float deltaPhi = phiClu-phiPho;
	float deltaEta = etaClu-etaPho;

	if ( deltaPhi > pi )  deltaPhi -= twopi;
	if ( deltaPhi < -pi) deltaPhi += twopi;
	deltaPhi=pow(deltaPhi,2);
	deltaEta=pow(deltaEta,2);
	float delta = sqrt( deltaPhi+deltaEta);
	if ( delta<0.1 && delta < minDelta ) {
	  minDelta=delta;
	  iMatch=index;

	}
	index++;
      } // End loop over uncorrected photons

      double wt=0.;
      if ( iMatch>-1 ) wt=1.;

      effEta_ ->Fill ( etaPho, wt);
      effPhi_ ->Fill ( phiPho, wt);

      /// Plot kinematic disctributions for matched photons
      if (iMatch>-1) {


	bool  phoIsInBarrel=false;
	bool  phoIsInEndcap=false;
	if ( fabs(localPhotons[iMatch].superCluster()->position().eta() ) <  1.479 ) {
	  phoIsInBarrel=true;
	} else {
	  phoIsInEndcap=true;
	}
	edm::Handle<EcalRecHitCollection>   ecalRecHitHandle;


        if (  phoIsInBarrel ) {
	  // Get handle to rec hits ecal barrel
	  evt.getByLabel(barrelEcalHits_, ecalRecHitHandle);
	  if (!ecalRecHitHandle.isValid()) {
	    edm::LogError("PhotonProducer") << "Error! Can't get the product "<<barrelEcalHits_.label();
	    return;
	  }

	}  else if ( phoIsInEndcap ) {

	  // Get handle to rec hits ecal encap
	  evt.getByLabel(endcapEcalHits_, ecalRecHitHandle);
	  if (!ecalRecHitHandle.isValid()) {
	    edm::LogError("PhotonProducer") << "Error! Can't get the product "<<endcapEcalHits_.label();
	    return;
	  }


	}
	const EcalRecHitCollection ecalRecHitCollection = *(ecalRecHitHandle.product());
	float e5x5= EcalClusterTools::e5x5( *( localPhotons[iMatch].superCluster()->seed()), &ecalRecHitCollection, &(*topology));
	float e3x3=   EcalClusterTools::e3x3(  *(    localPhotons[iMatch].superCluster()->seed()  ), &ecalRecHitCollection, &(*topology));
	float r9 =e3x3/(  localPhotons[iMatch].superCluster()->rawEnergy()+  localPhotons[iMatch].superCluster()->preshowerEnergy());




	h1_scE_->Fill( localPhotons[iMatch].superCluster()->energy() );
	h1_scEt_->Fill( localPhotons[iMatch].superCluster()->energy()/cosh(localPhotons[iMatch].superCluster()->position().eta()) );
	h1_scEta_->Fill( localPhotons[iMatch].superCluster()->position().eta() );
	h1_scPhi_->Fill( localPhotons[iMatch].superCluster()->position().phi() );

	float trueEta=  (*p)->momentum().eta() ;
	trueEta = etaTransformation(trueEta, (*p)->production_vertex()->position().z()/10. );
	h1_deltaEtaSC_ -> Fill(  localPhotons[iMatch].superCluster()->eta()- trueEta  );
	h1_deltaPhiSC_ -> Fill(  localPhotons[iMatch].phi()- (*p)->momentum().phi()  );

	h1_pho_E_->Fill( localPhotons[iMatch].energy() );
	h1_pho_Eta_->Fill( localPhotons[iMatch].eta() );
	h1_pho_Phi_->Fill( localPhotons[iMatch].phi() );

	h1_deltaEta_ -> Fill(  localPhotons[iMatch].eta()- (*p)->momentum().eta()  );
	h1_deltaPhi_ -> Fill(  localPhotons[iMatch].phi()- (*p)->momentum().phi()  );




	if ( phoIsInBarrel ) {
	  h1_recEoverTrueEBarrel_ -> Fill (localPhotons[iMatch].energy() / (*p)->momentum().e() );
	  h1_recESCoverTrueEBarrel_ -> Fill (localPhotons[iMatch].superCluster()->energy() / (*p)->momentum().e() );
	  h1_pho_R9Barrel_->Fill( r9 );



	  if ( r9 > 0.93 )
	    h1_e5x5_unconvBarrel_ -> Fill (  e5x5/ (*p)->momentum().e() );
	  else
	    h1_ePho_convBarrel_ -> Fill (  localPhotons[iMatch].energy()/ (*p)->momentum().e() );

	} else {
	  h1_recEoverTrueEEndcap_ -> Fill (localPhotons[iMatch].energy() / (*p)->momentum().e() );
	  h1_recESCoverTrueEEndcap_ -> Fill (localPhotons[iMatch].energy() / (*p)->momentum().e() );
	  h1_pho_R9Endcap_->Fill( r9 );


	  if ( r9 > 0.93 )
	    h1_e5x5_unconvEndcap_ -> Fill ( e5x5/ (*p)->momentum().e() );
          else
	    h1_ePho_convEndcap_ -> Fill (  localPhotons[iMatch].energy()/ (*p)->momentum().e() );

	}


      }




    } // End loop over MC particles

  }


}


float SimplePhotonAnalyzer::etaTransformation(  float EtaParticle , float Zvertex)  {

  //---Definitions
  const float PI    = 3.1415927;
  //UNUSED const float TWOPI = 2.0*PI;

  //---Definitions for ECAL
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479;

  //---ETA correction

  float Theta = 0.0  ;
  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0) Theta = Theta+PI ;
  float ETA = - log(tan(0.5*Theta));

  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - Zvertex ;
      float RR = Zlen/sinh(EtaParticle);
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+PI ;
      ETA = - log(tan(0.5*Theta));
    }
  //---Return the result
  return ETA;
  //---end
}




//========================================================================
void
SimplePhotonAnalyzer::endJob() {
//========================================================================



}
