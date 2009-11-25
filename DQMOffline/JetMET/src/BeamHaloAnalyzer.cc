#include "DQMOffline/JetMET/interface/BeamHaloAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//author : Ronny Remington, University of Florida
//date : 11/11/09

using namespace edm;
using namespace std;
using namespace reco;

int Phi_To_iPhi(float phi) 
{
  phi = phi < 0 ? phi + 2.*TMath::Pi() : phi ;
  float phi_degrees = phi * (360.) / ( 2. * TMath::Pi() ) ;
  int iPhi = (int) ( ( phi_degrees/5. ) + 1.);
   
  return iPhi < 73 ? iPhi : 73 ;
}


BeamHaloAnalyzer::BeamHaloAnalyzer( const edm::ParameterSet& iConfig)
{
  OutputFileName = iConfig.getParameter<std::string>("OutputFile"); 
 
  //Get Input Tags

  //Digi Level 
  IT_L1MuGMTReadout = iConfig.getParameter<edm::InputTag>("L1MuGMTReadoutLabel");
  
  //RecHit Level
  IT_CSCRecHit   = iConfig.getParameter<edm::InputTag>("CSCRecHitLabel");
  IT_EBRecHit    = iConfig.getParameter<edm::InputTag>("EBRecHitLabel");
  IT_EERecHit    = iConfig.getParameter<edm::InputTag>("EERecHitLabel");
  IT_ESRecHit    = iConfig.getParameter<edm::InputTag>("ESRecHitLabel");
  IT_HBHERecHit  = iConfig.getParameter<edm::InputTag>("HBHERecHitLabel");
  IT_HFRecHit    = iConfig.getParameter<edm::InputTag>("HFRecHitLabel");
  IT_HORecHit    = iConfig.getParameter<edm::InputTag>("HORecHitLabel");

  //Higher Level Reco 
  IT_CSCSegment = iConfig.getParameter<edm::InputTag>("CSCSegmentLabel");  
  IT_CosmicStandAloneMuon = iConfig.getParameter<edm::InputTag>("CosmicStandAloneMuonLabel"); 
  IT_BeamHaloMuon = iConfig.getParameter<edm::InputTag>("BeamHaloMuonLabel");
  IT_CollisionMuon = iConfig.getParameter<edm::InputTag>("CollisionMuonLabel");
  IT_CollisionStandAloneMuon  = iConfig.getParameter<edm::InputTag>("CollisionStandAloneMuonLabel"); 
  IT_met = iConfig.getParameter<edm::InputTag>("metLabel");
  IT_CaloTower = iConfig.getParameter<edm::InputTag>("CaloTowerLabel");
  IT_SuperCluster = iConfig.getParameter<edm::InputTag>("SuperClusterLabel");
  IT_Photon = iConfig.getParameter<edm::InputTag>("PhotonLabel") ;
  
  //Halo Data
  IT_CSCHaloData = iConfig.getParameter<edm::InputTag> ("CSCHaloDataLabel");
  IT_EcalHaloData = iConfig.getParameter<edm::InputTag>("EcalHaloDataLabel");
  IT_HcalHaloData = iConfig.getParameter<edm::InputTag>("HcalHaloDataLabel");
  IT_GlobalHaloData = iConfig.getParameter<edm::InputTag>("GlobalHaloDataLabel");
  IT_BeamHaloSummary = iConfig.getParameter<edm::InputTag>("BeamHaloSummaryLabel");
  FolderName = iConfig.getParameter<std::string>("folderName");
}


void BeamHaloAnalyzer::beginJob(const edm::EventSetup& iSetup){}

void BeamHaloAnalyzer::beginRun(const edm::Run&, const edm::EventSetup& iSetup){
  
  dqm = edm::Service<DQMStore>().operator->();
  if( dqm ) {
  

    // EcalHaloData
    dqm->setCurrentFolder(FolderName+"/EcalHaloData");
    ME["EcalHaloData_PhiWedgeMultiplicity"] = dqm->book1D("EcalHaloData_PhiWedgeMultiplicity","",20, -0.5, 19.5);
    ME["EcalHaloData_PhiWedgeEnergy"]       = dqm->book1D("EcalHaloData_PhiWedgeEnergy","", 50,-0.5,199.5);
    ME["EcalHaloData_PhiWedgeConstituents"] = dqm->book1D("EcalHaloData_PhiWedgeConstituents","",20,-0.5, 19.5);
    ME["EcalHaloData_PhiWedgeMinTime"]      = dqm->book1D("EcalHaloData_PhiWedgeMinTime","", 100, -225.0, 225.0);
    ME["EcalHaloData_PhiWedgeMaxTime"]      = dqm->book1D("EcalHaloData_PhiWedgeMaxTime","", 100, -225.0, 225.0);
    ME["EcalHaloData_PhiWedgeiPhi"]         = dqm->book1D("EcalHaloData_PhiWedgeiPhi","", 72, 0.5, 72.5) ;
    ME["EcalHaloData_PhiWedgePlusZDirectionConfidence"] = dqm->book1D("EcalHaloData_PlusZDirectionConfidence","",  50, 0., 1.0);
    ME["EcalHaloData_PhiWedgeMinVsMaxTime"] = dqm->book2D("EcalHaloData_PhiWedgeMinVsMaxTime","", 50,-100.0, 100.0, 50, -100.0, 100.0);
    ME["EcalHaloData_SuperClusterShowerShapes"]  = dqm->book2D("EcalHaloData_SuperClusterShowerShapes","", 25,0.0, TMath::Pi(), 25,0.0, 2.0);
    
    // HcalHaloData
    dqm->setCurrentFolder(FolderName+"/HcalHaloData");    
    ME["HcalHaloData_PhiWedgeMultiplicity"] = dqm->book1D("HcalHaloData_PhiWedgeMultiplicity","", 20, -0.5, 19.5);
    ME["HcalHaloData_PhiWedgeEnergy"]       = dqm->book1D("HcalHaloData_PhiWedgeEnergy", "", 50,-0.5,199.5);
    ME["HcalHaloData_PhiWedgeConstituents"] = dqm->book1D("HcalHaloData_PhiWedgeConstituents","", 20,-0.5, 19.5);
    ME["HcalHaloData_PhiWedgeiPhi"]         = dqm->book1D("HcalHaloData_PhiWedgeiPhi","", 72, 0.5,72.5);
    ME["HcalHaloData_PhiWedgeMinTime"]      = dqm->book1D("HcalHaloData_PhiWedgeMinTime", "", 50, -100.0, 100.0);
    ME["HcalHaloData_PhiWedgeMaxTime"]      = dqm->book1D("HcalHaloData_PhiWedgeMaxTime", "", 50, -100.0, 100.0);
    ME["HcalHaloData_PhiWedgePlusZDirectionConfidence"] = dqm->book1D("HcalHaloData_PlusZDirectionConfidence","",  50, 0., 1.0);
    ME["HcalHaloData_PhiWedgeMinVsMaxTime"] = dqm->book2D("HcalHaloData_PhiWedgeMinVsMaxTime","" , 50,-100.0, 100.0, 50, -100.0, 100.0);
    
    // CSCHaloData
    dqm->setCurrentFolder(FolderName+"/CSCHaloData");
    ME["CSCHaloData_TrackMultiplicity"]  = dqm->book1D("CSCHaloData_TrackMultiplicity", "", 15, -0.5, 14.5);
    ME["CSCHaloData_TrackMultiplicityMEPlus"]  = dqm->book1D("CSCHaloData_TrackMultiplicityMEPlus", "", 15, -0.5, 14.5);
    ME["CSCHaloData_TrackMultiplicityMEMinus"]  = dqm->book1D("CSCHaloData_TrackMultiplicityMEMinus", "", 15, -0.5, 14.5);
    ME["CSCHaloData_InnerMostTrackHitXY"]  = dqm->book2D("CSCHaloData_InnerMostTrackHitXY","", 100,-700,700,100, -700,700);
    ME["CSCHaloData_InnerMostTrackHitR"]  = dqm->book1D("CSCHaloData_InnerMostTrackHitR", "", 400, -0.5, 799.5);
    ME["CSCHaloData_InnerMostTrackHitiPhi"]  = dqm->book1D("CSCHaloData_InnerMostTrackHitiPhi","", 72, 0.5, 72.5);
    ME["CSCHaloData_HaloTriggersMEPlus"]  = dqm->book1D("CSCHaloData_HaloTriggersMEPlus", "", 10, -0.5, 9.5);
    ME["CSCHaloData_HaloTriggersMEMinus"]  = dqm->book1D("CSCHaloData_HaloTriggersMEMinus", "" , 10, -0.5, 9.5);
    ME["CSCHaloData_HaloTriggers"]  = dqm->book1D("CSCHaloData_HaloTriggers", "", 10, -0.5, 9.5);
    ME["CSCHaloData_CaloPointingHaloTrackMultiplicity"] = dqm->book1D("CSCHaloData_CaloPointingHaloTrackMultiplicity","", 10, -0.5, 9.5 );

    // GlobalHaloData
    dqm->setCurrentFolder(FolderName+"/GlobalHaloData");
    ME["GlobalHaloData_MExCorrection"]  = dqm->book1D("GlobalHaloData_MExCorrection", "" , 200, -200., 200.);
    ME["GlobalHaloData_MEyCorrection"]  = dqm->book1D("GlobalHaloData_MEyCorrection", "" , 200, -200., 200.);
    ME["GlobalHaloData_SumEtCorrection"] = dqm->book1D("GlobalHaloData_SumEtCorrection", "" , 200, -0.5, 399.5);
    ME["GlobalHaloData_HaloCorrectedMET"] = dqm->book1D("GlobalHaloData_HaloCorrectedMET", "" , 500, -0.5, 1999.5);
    ME["GlobalHaloData_RawMETMinusHaloCorrectedMET"] = dqm->book1D("GlobalHaloData_RawMETMinusHaloCorrectedMET","" , 250, -500., 500.);
    ME["GlobalHaloData_RawMETOverSumEt"]  = dqm->book1D("GlobalHaloData_RawMETOverSumEt","" , 100, 0.0, 1.0);
    ME["GlobalHaloData_MatchedHcalPhiWedgeMultiplicity"] = dqm->book1D("GlobalHaloData_MatchedHcalPhiWedgeMultiplicity","", 15, -0.5, 14.5);    
    ME["GlobalHaloData_MatchedHcalPhiWedgeEnergy"]       = dqm->book1D("GlobalHaloData_MatchedHcalPhiWedgeEnergy", "", 50,-0.5,199.5);
    ME["GlobalHaloData_MatchedHcalPhiWedgeConstituents"] = dqm->book1D("GlobalHaloData_MatchedHcalPhiWedgeConstituents","", 20,-0.5, 19.5);
    ME["GlobalHaloData_MatchedHcalPhiWedgeiPhi"]         = dqm->book1D("GlobalHaloData_MatchedHcalPhiWedgeiPhi","", 72, 0.5,72.5);
    ME["GlobalHaloData_MatchedHcalPhiWedgeMinTime"]      = dqm->book1D("GlobalHaloData_MatchedHcalPhiWedgeMinTime", "", 50, -100.0, 100.0);
    ME["GlobalHaloData_MatchedHcalPhiWedgeMaxTime"]      = dqm->book1D("GlobalHaloData_MatchedHcalPhiWedgeMaxTime", "", 50, -100.0, 100.0);
    ME["GlobalHaloData_MatchedEcalPhiWedgeMultiplicity"] = dqm->book1D("GlobalHaloData_MatchedEcalPhiWedgeMultiplicity","", 15, -0.5, 14.5);
    ME["GlobalHaloData_MatchedEcalPhiWedgeEnergy"]       = dqm->book1D("GlobalHaloData_MatchedEcalPhiWedgeEnergy", "", 50,-0.5,199.5);
    ME["GlobalHaloData_MatchedEcalPhiWedgeConstituents"] = dqm->book1D("GlobalHaloData_MatchedEcalPhiWedgeConstituents","", 20,-0.5, 19.5);
    ME["GlobalHaloData_MatchedEcalPhiWedgeiPhi"]         = dqm->book1D("GlobalHaloData_MatchedEcalPhiWedgeiPhi","", 72, 0.5,72.5);
    ME["GlobalHaloData_MatchedEcalPhiWedgeMinTime"]      = dqm->book1D("GlobalHaloData_MatchedEcalPhiWedgeMinTime", "", 50, -100.0, 100.0);
    ME["GlobalHaloData_MatchedEcalPhiWedgeMaxTime"]      = dqm->book1D("GlobalHaloData_MatchedEcalPhiWedgeMaxTime", "", 50, -100.0, 100.0);

    // BeamHaloSummary 
    dqm->setCurrentFolder(FolderName+"/BeamHaloSummary");
    ME["BeamHaloSummary_Id"] = dqm->book1D("BeamHaloSumamry_Id", "", 11, 0.5,11.5);
    ME["BeamHaloSummary_Id"] ->setBinLabel(1,"CSC Loose");
    ME["BeamHaloSummary_Id"] ->setBinLabel(2,"CSC Tight");
    ME["BeamHaloSummary_Id"] ->setBinLabel(3,"Ecal Loose");
    ME["BeamHaloSummary_Id"] ->setBinLabel(4,"Ecal Tight");
    ME["BeamHaloSummary_Id"] ->setBinLabel(5,"Hcal Loose");
    ME["BeamHaloSummary_Id"] ->setBinLabel(6,"Hcal Tight");
    ME["BeamHaloSummary_Id"] ->setBinLabel(7,"Global Loose");
    ME["BeamHaloSummary_Id"] ->setBinLabel(8,"Global Tight");
    ME["BeamHaloSummary_Id"] ->setBinLabel(9,"Event Loose");
    ME["BeamHaloSummary_Id"] ->setBinLabel(10,"Event Tight");
    ME["BeamHaloSummary_Id"] ->setBinLabel(11,"Nothing");

    ME["BeamHaloSummary_BXN"] = dqm->book2D("BeamHaloSummary_BXN", "",11, 0.5, 11.5, 4000, -0.5,3999.5);
    ME["BeamHaloSummary_BXN"] ->setBinLabel(1,"CSC Loose");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(2,"CSC Tight");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(3,"Ecal Loose");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(4,"Ecal Tight");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(5,"Hcal Loose");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(6,"Hcal Tight");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(7,"Global Loose");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(8,"Global Tight");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(9,"Event Loose");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(10,"Event Tight");
    ME["BeamHaloSummary_BXN"] ->setBinLabel(11,"Nothing");

    // Extra
    dqm->setCurrentFolder(FolderName+"/ExtraHaloData");
    ME["Extra_CSCActivityWithMET"]= dqm->book2D("Extra_CSCActivityWithMET", "", 4, 0.5, 4.5, 4, 0.5, 4.5);
    ME["Extra_CSCActivityWithMET"]->setBinLabel(1,"Track",1);
    ME["Extra_CSCActivityWithMET"]->setBinLabel(1,"Track",2);
    ME["Extra_CSCActivityWithMET"]->setBinLabel(2, "Segments",1);
    ME["Extra_CSCActivityWithMET"]->setBinLabel(2, "Segments",2);
    ME["Extra_CSCActivityWithMET"]->setBinLabel(3, "RecHits", 1);
    ME["Extra_CSCActivityWithMET"]->setBinLabel(3, "RecHits", 2);
    ME["Extra_CSCActivityWithMET"]->setBinLabel(4, "Nothing", 1);
    ME["Extra_CSCActivityWithMET"]->setBinLabel(4, "Nothing", 2);
    ME["Extra_HcalToF"]  = dqm->book2D("HcalToF","" , 83,-41.5,41.5 , 1000, -125., 125.); 
    ME["Extra_HcalToF_HaloId"]  = dqm->book2D("HcalToF","", 83,-41.5,41.5 , 1000, -125., 125.); 
    ME["Extra_EcalToF"]  = dqm->book2D("EcalToF","",  171,-85.5,85.5 , 2000, -225., 225.); 
    ME["Extra_EcalToF_HaloId"]  = dqm->book2D("EcalToF","",  171,-85.5,85.5 , 2000, -225., 225.); 


  }
}

void BeamHaloAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  EventID TheEvent = iEvent.id();
  int BXN = iEvent.bunchCrossing() ;
  //  int TheEventNumber = TheEvent.event();
  
  //Get CSC Geometry
  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  //Get CaloGeometry
  edm::ESHandle<CaloGeometry> TheCaloGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCaloGeometry);

  // Get BeamHalo Muon Collection from Cosmic Muon Reconstruction 
  /*
  edm::Handle<reco::MuonCollection> TheBeamHaloMuons;
  iEvent.getByLabel(IT_BeamHaloMuon,TheBeamHaloMuons);
  if( TheBeamHaloMuons.isValid() )
    {
      for( reco::MuonCollection::const_iterator iMuon = TheBeamHaloMuons->begin() ; iMuon != TheBeamHaloMuons->end(); iMuon++)
        {
        }
    }
  */

  //Get Collision Stand-Alone Muons
  /*
  edm::Handle<reco::TrackCollection> TheSAMuons;
  iEvent.getByLabel( IT_CollisionStandAloneMuon, TheSAMuons);
  if( TheSAMuons.isValid() )
    {
      for( reco::TrackCollection::const_iterator sa = TheSAMuons->begin() ; sa != TheSAMuons->end() ; sa++ )
	{
	}
    }
  */

  //Get Stand-alone Muons from Cosmic Muon Reconstruction
  edm::Handle< reco::TrackCollection > TheCosmics;
  iEvent.getByLabel(IT_CosmicStandAloneMuon, TheCosmics);
  bool CSCTrackPlus = false; bool CSCTrackMinus = false;
  if( TheCosmics.isValid() )
    {
      for( reco::TrackCollection::const_iterator cosmic = TheCosmics ->begin() ; cosmic != TheCosmics->end() ; cosmic++ )
	{
	  if( !CSCTrackPlus || !CSCTrackMinus)
            {
              if( cosmic->eta() > 0 || cosmic->outerPosition().z() > 0  || cosmic->innerPosition().z() > 0 ) CSCTrackPlus = true ;
	      else if( cosmic->eta() < 0 || cosmic->outerPosition().z() < 0 || cosmic->innerPosition().z() < 0) CSCTrackMinus = true;
            }
	}
    }

  //Get CSC Segments
  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByLabel(IT_CSCSegment, TheCSCSegments);

  // Group segments according to endcaps
  std::vector< CSCSegment> vCSCSegments_Plus;
  std::vector< CSCSegment> vCSCSegments_Minus;

  bool CSCSegmentPlus = false; 
  bool CSCSegmentMinus=false;
  if( TheCSCSegments.isValid() ) 
    {
      for(CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end(); iSegment++) 
	{
	  const std::vector<CSCRecHit2D> vCSCRecHits = iSegment->specificRecHits();
	  CSCDetId iDetId  = (CSCDetId)(*iSegment).cscDetId();
	  
	  if ( iDetId.endcap() == 1 ) vCSCSegments_Plus.push_back( *iSegment );
	  else vCSCSegments_Minus.push_back( *iSegment );
	}      
    }
  
  // Are there segments on the plus/minus side?  
  if( vCSCSegments_Plus.size() ) CSCSegmentPlus = true;
  if( vCSCSegments_Minus.size() ) CSCSegmentMinus = true;

  //Get CSC RecHits
  Handle<CSCRecHit2DCollection> TheCSCRecHits;
  iEvent.getByLabel(IT_CSCRecHit, TheCSCRecHits);
  bool CSCRecHitPlus = false; 
  bool CSCRecHitMinus = false;
  if( TheCSCRecHits.isValid() )
    {
      for(CSCRecHit2DCollection::const_iterator iCSCRecHit = TheCSCRecHits->begin();   iCSCRecHit != TheCSCRecHits->end(); iCSCRecHit++ )
	{
	  DetId TheDetUnitId(iCSCRecHit->geographicalId());
	  const GeomDetUnit *TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);
	  LocalPoint TheLocalPosition = iCSCRecHit->localPosition();
	  const BoundPlane& TheSurface = TheUnit->surface();
	  GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);

	  //Are there hits on the plus/minus side?
	  if ( TheGlobalPosition.z() > 0 ) CSCRecHitPlus = true;
	  else CSCRecHitMinus = true;
	}
    }
  //Get L1MuGMT 
  //edm::Handle < L1MuGMTReadoutCollection > TheL1GMTReadout ;
  //iEvent.getByLabel (IT_L1MuGMTReadout, TheL1GMTReadout);
  
  //Get  EB RecHits
  edm::Handle<EBRecHitCollection> TheEBRecHits;
  iEvent.getByLabel(IT_EBRecHit, TheEBRecHits);
  int EBHits=0;
  if( TheEBRecHits.isValid() )
    {
      for( EBRecHitCollection::const_iterator iEBRecHit = TheEBRecHits->begin() ; iEBRecHit != TheEBRecHits->end(); iEBRecHit++)
	{
	  if( iEBRecHit->energy() < 0.5 ) continue;
	  DetId id = DetId( iEBRecHit->id() ) ;
	  EBDetId EcalId ( id.rawId() );
	  int ieta = EcalId.ieta() ;
	  ME["Extra_EcalToF"] ->Fill(ieta, iEBRecHit->time() );
	  EBHits++;
	}
    }
  
  //Get EE RecHits
  /*
  edm::Handle<EERecHitCollection> TheEERecHits;
  iEvent.getByLabel(IT_EERecHit, TheEERecHits);
  if( TheEERecHits.isValid() )
    {
      for( EERecHitCollection::const_iterator iEERecHit = TheEERecHits->begin() ; iEERecHit != TheEERecHits->end(); iEERecHit++)
	{
	  if( iEERecHit->energy() < 0.2 ) continue;
	}
    }
  */

  //Get ES RecHits
  /*
  edm::Handle<ESRecHitCollection> TheESRecHits;
  iEvent.getByLabel(IT_ESRecHit, TheESRecHits);
  if( TheESRecHits.isValid() )
    {
      for( ESRecHitCollection::const_iterator iESRecHit = TheESRecHits->begin() ; iESRecHit != TheESRecHits->end(); iESRecHit++)
	{
	  if( iESRecHit->energy() < 0.2 ) continue;
	}
    }
  */

  //Get HB/HE RecHits
  edm::Handle<HBHERecHitCollection> TheHBHERecHits;
  iEvent.getByLabel(IT_HBHERecHit, TheHBHERecHits);
  if( TheHBHERecHits.isValid() )
    {
      for( HBHERecHitCollection::const_iterator iHBHERecHit = TheHBHERecHits->begin(); iHBHERecHit != TheHBHERecHits->end(); iHBHERecHit++)  
	{
	  if( iHBHERecHit->energy() < 1.) continue;
	  HcalDetId id = HcalDetId( iHBHERecHit->id() );
	  ME["Extra_HcalToF"]->Fill( id.ieta(), iHBHERecHit->time() ) ;
	}
    }

  //Get HF RecHits
  /*
  edm::Handle<HFRecHitCollection> TheHFRecHits;
  iEvent.getByLabel(IT_HFRecHit, TheHFRecHits);
  if( TheHFRecHits.isValid() )
    {
      for( HFRecHitCollection::const_iterator iHFRecHit = TheHFRecHits->begin(); iHFRecHit != TheHFRecHits->end(); iHFRecHit++ )
	{
	}
    }
  */

  //Get HO RecHits
  /*
  edm::Handle<HORecHitCollection> TheHORecHits;
  iEvent.getByLabel(IT_HORecHit, TheHORecHits);
  if( TheHORecHits.isValid() )
    {
      for(HORecHitCollection::const_iterator iHORecHit = TheHORecHits->begin(); iHORecHit != TheHORecHits->end(); iHORecHit++ )
	{
	}
    }
  */

  //Get ECAL Barrel SuperClusters                                                                                                                         
  /*
  edm::Handle<reco::SuperClusterCollection> TheSuperClusters;
  iEvent.getByLabel(IT_SuperCluster, TheSuperClusters);
  if( TheSuperClusters.isValid() )
    {
      for( SuperClusterCollection::const_iterator iSCluster = TheSuperClusters->begin() ; iSCluster != TheSuperClusters->end() ; iSCluster++ )
        {
        }
    }
  */
  
  // Get Photons
  /*
  edm::Handle<reco::PhotonCollection> ThePhotons;
  iEvent.getByLabel(IT_Photon, ThePhotons);
  if(ThePhotons.isValid())
    {      
      for(reco::PhotonCollection::const_iterator iPhoton = ThePhotons->begin() ; iPhoton != ThePhotons->end() ; iPhoton++ )
	{
	}
    }
  */
 
  //Get CaloTowers
  /*
  edm::Handle<edm::View<Candidate> > TheCaloTowers;
  iEvent.getByLabel(IT_CaloTower,TheCaloTowers);
  for( edm::View<Candidate>::const_iterator iCandidate = TheCaloTowers->begin() ; iCandidate != TheCaloTowers->end() ; iCandidate++ )
    {
      const Candidate* c = &(*iCandidate);
      if ( c )
	{
	//  const CaloTower* iTower = dynamic_cast<const CaloTower*> (c);
	}
    }
  */
  
  //Get MET
  edm::Handle< reco::CaloMETCollection > TheCaloMET;
  iEvent.getByLabel(IT_met, TheCaloMET);

  //Get CSCHaloData
  edm::Handle<reco::CSCHaloData> TheCSCDataHandle;
  iEvent.getByLabel(IT_CSCHaloData,TheCSCDataHandle);
  int TheHaloOrigin = 0;
  if (TheCSCDataHandle.isValid())
    {
      const CSCHaloData CSCData = (*TheCSCDataHandle.product());
      ME["CSCHaloData_CaloPointingHaloTrackMultiplicity"]->Fill(CSCData.GetCSCTrackImpactPositions().size());                                                      

      if( CSCData.NumberOfHaloTriggers(1) && !CSCData.NumberOfHaloTriggers(-1) ) TheHaloOrigin = 1;
      else if ( CSCData.NumberOfHaloTriggers(-1) && !CSCData.NumberOfHaloTriggers(1) ) TheHaloOrigin = -1 ;

      for( std::vector<GlobalPoint>::const_iterator i=CSCData.GetCSCTrackImpactPositions().begin();  i != CSCData.GetCSCTrackImpactPositions().end() ; i++ )   
	{                                                                                                                                                      
	  //ME["CSCHaloData_InnerMostTrackHitXY"]->Fill( i->x(), i->y() );
	  ME["CSCHaloData_InnerMostTrackHitR"]  ->Fill( TMath::Sqrt( i->x()*i->x() + i->y()*i->y() ));
	  ME["CSCHaloData_InnerMostTrackHitiPhi"]  ->Fill( Phi_To_iPhi( i->phi())); 
	}                          
      ME["CSCHaloData_HaloTriggersMEPlus"]   -> Fill ( CSCData.NumberOfHaloTriggers(1) );
      ME["CSCHaloData_HaloTriggersMEMinus"]  -> Fill ( CSCData.NumberOfHaloTriggers(-1));
      ME["CSCHaloData_HaloTriggers"]  -> Fill ( CSCData.NumberOfHaloTriggers());
      ME["CSCHaloData_TrackMultiplicityMEPlus"] ->Fill ( CSCData.NumberOfHaloTracks(1) );
      ME["CSCHaloData_TrackMultiplicityMEMinus"] ->Fill ( CSCData.NumberOfHaloTracks(-1) );
      ME["CSCHaloData_TrackMultiplicity"]->Fill( CSCData.GetTracks().size() );
    }

  //Get EcalHaloData 
  edm::Handle<reco::EcalHaloData> TheEcalHaloData;
  iEvent.getByLabel(IT_EcalHaloData, TheEcalHaloData );
  if( TheEcalHaloData.isValid() ) 
    {
      const EcalHaloData EcalData = (*TheEcalHaloData.product()); 
      std::vector<PhiWedge> EcalWedges = EcalData.GetPhiWedges();                                                                                              
      for(std::vector<PhiWedge>::const_iterator iWedge = EcalWedges.begin() ; iWedge != EcalWedges.end(); iWedge ++ )                                  
	{                                                                                                                                                     
	  ME["EcalHaloData_PhiWedgeEnergy"]->Fill( iWedge->Energy() );
	  ME["EcalHaloData_PhiWedgeConstituents"]->Fill( iWedge->NumberOfConstituents() ) ;
	  ME["EcalHaloData_PhiWedgeMinTime"]     ->Fill( iWedge->MinTime() );
	  ME["EcalHaloData_PhiWedgeMaxTime"]     ->Fill( iWedge->MaxTime() );
	  cout << "Min,Max Time       " <<  iWedge->MinTime() << "         " << iWedge->MaxTime() << endl;
	  
	  ME["EcalHaloData_PhiWedgePlusZDirectionConfidence"]->Fill( iWedge->PlusZDirectionConfidence() );
	  ME["EcalHaloData_PhiWedgeMinVsMaxTime"]->Fill(iWedge->MinTime() , iWedge->MaxTime() ) ;
	  ME["EcalHaloData_PhiWedgeiPhi"]->Fill(iWedge->iPhi() ) ;
	}      

      ME["EcalHaloData_PhiWedgeMultiplicity"]->Fill( EcalWedges.size() );

      edm::ValueMap<float> vm_Angle = EcalData.GetShowerShapesAngle();
      edm::ValueMap<float> vm_Roundness = EcalData.GetShowerShapesRoundness();
      //Access selected SuperClusters
      for(unsigned int n = 0 ; n < EcalData.GetSuperClusters().size() ; n++ )
	{
	  edm::Ref<SuperClusterCollection> cluster(EcalData.GetSuperClusters(), n );
	  float angle = vm_Angle[cluster];
	  float roundness = vm_Roundness[cluster];
	  ME["EcalHaloData_SuperClusterShowerShapes"]->Fill(angle, roundness);
	}
    }

  //Get HcalHaloData
  edm::Handle<reco::HcalHaloData> TheHcalHaloData;
  iEvent.getByLabel(IT_HcalHaloData ,TheHcalHaloData );
  if( TheHcalHaloData.isValid( ) )
    {
      const HcalHaloData HcalData = (*TheHcalHaloData.product());                                                                
      std::vector<PhiWedge> HcalWedges = HcalData.GetPhiWedges();                                                                                   
      ME["HcalHaloData_PhiWedgeMultiplicity"] ->Fill( HcalWedges.size() );
      for(std::vector<PhiWedge>::const_iterator iWedge = HcalWedges.begin() ; iWedge != HcalWedges.end(); iWedge ++ )                               
	{                                                                                                                                        
	  ME["HcalHaloData_PhiWedgeEnergy"]       ->Fill( iWedge->Energy() );
	  ME["HcalHaloData_PhiWedgeConstituents"] ->Fill( iWedge->NumberOfConstituents() );
	  ME["HcalHaloData_PhiWedgeiPhi"]         ->Fill( iWedge->iPhi() );
	  ME["HcalHaloData_PhiWedgeMinTime"]      ->Fill( iWedge->MinTime() );
	  ME["HcalHaloData_PhiWedgeMaxTime"]      ->Fill( iWedge->MaxTime() );
	  ME["HcalHaloData_PhiWedgePlusZDirectionConfidence"] ->Fill( iWedge->PlusZDirectionConfidence() );
	  ME["HcalHaloData_PhiWedgeMinVsMaxTime"]  ->Fill( iWedge->MinTime() , iWedge->MaxTime() );
	}
    }
  
  //Get GlobalHaloData
  edm::Handle<reco::GlobalHaloData> TheGlobalHaloData;
  iEvent.getByLabel(IT_GlobalHaloData, TheGlobalHaloData );
  if( TheGlobalHaloData.isValid() ) 
    {
      const GlobalHaloData GlobalData =(*TheGlobalHaloData.product());                                                           
      if( TheCaloMET.isValid() ) 
	{
	  // Get Raw Uncorrected CaloMET
	  const CaloMETCollection *calometcol = TheCaloMET.product();
	  const CaloMET *RawMET = &(calometcol->front());

	  // Get BeamHalo Corrected CaloMET 
	  const CaloMET CorrectedMET = GlobalData.GetCorrectedCaloMET(*RawMET);
	  ME["GlobalHaloData_MExCorrection"]  ->Fill( GlobalData.DeltaMEx() );
	  ME["GlobalHaloData_MEyCorrection"]  ->Fill( GlobalData.DeltaMEy() );
	  //	  ME["GlobalHaloData_SumEtCorrection"] ->Fill( GlobalData.DeltaSumEt() );
	  ME["GlobalHaloData_HaloCorrectedMET"]->Fill(CorrectedMET.pt() );
	  ME["GlobalHaloData_RawMETMinusHaloCorrectedMET"] ->Fill( RawMET->pt() - CorrectedMET.pt() );
	  if( RawMET->sumEt() )
	    ME["GlobalHaloData_RawMETOverSumEt"] ->Fill( RawMET->pt() / RawMET->sumEt() ); 
	  
	}                

      // Get Matched Hcal Phi Wedges
      std::vector<PhiWedge> HcalWedges = GlobalData.GetMatchedHcalPhiWedges();
      ME["GlobalHaloData_MatchedHcalPhiWedgeMultiplicity"] ->Fill(HcalWedges.size());
      // Loop over Matched Hcal Phi Wedges
      for( std::vector<PhiWedge>::const_iterator iWedge = HcalWedges.begin() ; iWedge != HcalWedges.end() ; iWedge ++ )
	{
	  ME["GlobalHaloData_MatchedHcalPhiWedgeEnergy"]       ->Fill( iWedge->Energy() );
	  ME["GlobalHaloData_MatchedHcalPhiWedgeConstituents"] ->Fill( iWedge->NumberOfConstituents());
	  ME["GlobalHaloData_MatchedHcalPhiWedgeiPhi"]         ->Fill( iWedge->iPhi() );
	  ME["GlobalHaloData_MatchedHcalPhiWedgeMinTime"]      ->Fill( iWedge->MinTime() );
	  ME["GlobalHaloData_MatchedHcalPhiWedgeMaxTime"]      ->Fill( iWedge->MaxTime() );
	  if( TheHBHERecHits.isValid() )
	    {
	      for( HBHERecHitCollection::const_iterator iHBHERecHit = TheHBHERecHits->begin(); iHBHERecHit != TheHBHERecHits->end(); iHBHERecHit++)  
		{
		  HcalDetId id = HcalDetId( iHBHERecHit->id() ) ;
		  int iphi = id.iphi() ;
		  if( iphi != iWedge->iPhi() ) continue;
		  if( iHBHERecHit->energy() < 1.0) continue;  // Otherwise there are thousands of hits per event (even with negative energies)
		  
		  float time = iHBHERecHit->time();
		  int ieta = id.ieta();
		  ME["Extra_HcalToF_HaloId"] ->Fill( ieta, time );
		}
	    }
	}

      // Get Matched Hcal Phi Wedges
      std::vector<PhiWedge> EcalWedges = GlobalData.GetMatchedEcalPhiWedges();
      ME["GlobalHaloData_MatchedEcalPhiWedgeMultiplicity"] ->Fill(EcalWedges.size());
      for( std::vector<PhiWedge>::const_iterator iWedge = EcalWedges.begin() ; iWedge != EcalWedges.end() ; iWedge ++ )
	{
	  ME["GlobalHaloData_MatchedEcalPhiWedgeEnergy"]       ->Fill(iWedge->Energy());
	  ME["GlobalHaloData_MatchedEcalPhiWedgeConstituents"] ->Fill(iWedge->NumberOfConstituents());
	  ME["GlobalHaloData_MatchedEcalPhiWedgeiPhi"]         ->Fill(iWedge->iPhi());
	  ME["GlobalHaloData_MatchedEcalPhiWedgeMinTime"]      ->Fill(iWedge->MinTime());
	  ME["GlobalHaloData_MatchedEcalPhiWedgeMaxTime"]      ->Fill(iWedge->MaxTime());

	  if( TheEBRecHits.isValid() ) 
	    {
	      for( EBRecHitCollection::const_iterator iEBRecHit = TheEBRecHits->begin() ; iEBRecHit != TheEBRecHits->end(); iEBRecHit++ )
		{
		  if( iEBRecHit->energy() < 0.5 ) continue;
		  DetId id = DetId( iEBRecHit->id() ) ;
		  EBDetId EcalId ( id.rawId() );
		  int iPhi = EcalId.iphi() ;
		  iPhi = (iPhi-1)/5 + 1;
		  if( iPhi != iWedge->iPhi() ) continue;
		  ME["Extra_EcalToF_HaloId"] ->Fill(EcalId.ieta(), iEBRecHit->time() );
		}
	    }
	}
    }


  // Get BeamHaloSummary 
  edm::Handle<BeamHaloSummary> TheBeamHaloSummary ;
  iEvent.getByLabel(IT_BeamHaloSummary, TheBeamHaloSummary) ;
  if( TheBeamHaloSummary.isValid() ) 
    {
      const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );
      if( TheSummary.CSCLooseHaloId() ) 
	{
	  ME["BeamHaloSummary_Id"] ->Fill(1);
	  ME["BeamHaloSummary_BXN"] -> Fill( 1, BXN );
	}
      if( TheSummary.CSCTightHaloId() ) 
	{
	  ME["BeamHaloSummary_Id"] ->Fill(2);
	  ME["BeamHaloSummary_BXN"] -> Fill( 2, BXN );
	}
      if( TheSummary.EcalLooseHaloId() )
	{
	  ME["BeamHaloSummary_Id"] ->Fill(3);
	  ME["BeamHaloSummary_BXN"] -> Fill( 3, BXN );
	}
      if( TheSummary.EcalTightHaloId() ) 
	{
	  ME["BeamHaloSummary_Id"] ->Fill(4);
	  ME["BeamHaloSummary_BXN"] -> Fill( 4, BXN );
	}
      if( TheSummary.HcalLooseHaloId() ) 
	{
	  ME["BeamHaloSummary_Id"] ->Fill(5);
	  ME["BeamHaloSummary_BXN"] -> Fill( 5, BXN );
	}
      if( TheSummary.HcalTightHaloId() ) 
	{
	  ME["BeamHaloSummary_Id"] ->Fill(6);
	  ME["BeamHaloSummary_BXN"] -> Fill( 6, BXN );
	}
      if( TheSummary.GlobalLooseHaloId()) 
	{
	  ME["BeamHaloSummary_Id"] ->Fill(7);
	  ME["BeamHaloSummary_BXN"] -> Fill( 7, BXN );
	}
      if( TheSummary.GlobalTightHaloId() )
	{
	  ME["BeamHaloSummary_Id"] ->Fill(8);	
	  ME["BeamHaloSummary_BXN"] -> Fill( 8, BXN );
	}
      if( TheSummary.LooseId() ) 
	{
	  ME["BeamHaloSummary_Id"] ->Fill(9);
	  ME["BeamHaloSummary_BXN"] -> Fill( 9, BXN );
	}
      if( TheSummary.TightId() )
	{
	  ME["BeamHaloSummary_Id"] ->Fill(10);
	  ME["BeamHaloSummary_BXN"] -> Fill( 10, BXN );
	}
      if( !TheSummary.EcalLooseHaloId()  && !TheSummary.HcalLooseHaloId() && !TheSummary.CSCLooseHaloId() && !TheSummary.GlobalLooseHaloId() )
	{
	  ME["BeamHaloSummary_Id"] ->Fill(11);
	  ME["BeamHaloSummary_BXN"] -> Fill( 11, BXN );
	}
    }

  if( TheCaloMET.isValid() )
    {
      const CaloMETCollection *calometcol = TheCaloMET.product();
      const CaloMET *calomet = &(calometcol->front());
      
      //Fill CSC Activity Plot 
      if( calomet->pt() > 15.0 ) 
	{
	  if( TheHaloOrigin > 0 )
	    {
	      if( CSCTrackPlus && CSCTrackMinus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(1,1);
	      else if( CSCTrackPlus && CSCSegmentMinus) 
		ME["Extra_CSCActivityWithMET"]->Fill(1,2);
	      else if( CSCTrackPlus && CSCRecHitMinus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(1,3);
	      else if( CSCTrackPlus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(1,4);
	      else if( CSCSegmentPlus && CSCTrackMinus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(2,1);
	      else if( CSCSegmentPlus && CSCSegmentMinus )
		ME["Extra_CSCActivityWithMET"]-> Fill(2,2);
	      else if( CSCSegmentPlus && CSCRecHitMinus   )
		ME["Extra_CSCActivityWithMET"]-> Fill(2,3);
	      else if( CSCSegmentPlus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(2,4 );
	      else if( CSCRecHitPlus && CSCTrackMinus  ) 
		ME["Extra_CSCActivityWithMET"]->Fill(3,1);
	      else if( CSCRecHitPlus && CSCSegmentMinus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(3,2);
	      else if( CSCRecHitPlus && CSCRecHitMinus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(3,3);
	      else if( CSCRecHitPlus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(3,4);
	      else 
		ME["Extra_CSCActivityWithMET"]->Fill(4,4);
	    }
	  else if( TheHaloOrigin < 0 )
	    {
	      if( CSCTrackMinus && CSCTrackPlus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(1,1);
	      else if( CSCTrackMinus && CSCSegmentPlus)
		ME["Extra_CSCActivityWithMET"]->Fill(1,2);
	      else if( CSCTrackMinus && CSCRecHitPlus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(1,3);
	      else if( CSCTrackMinus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(1,4);
	      else if( CSCSegmentMinus && CSCTrackPlus) 
		ME["Extra_CSCActivityWithMET"]->Fill(2,1);
	      else if( CSCSegmentMinus && CSCSegmentPlus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(2,2 );
	      else if( CSCSegmentMinus && CSCRecHitPlus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(2,3);
	      else if( CSCSegmentMinus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(2,4);
	      else if( CSCRecHitMinus && CSCTrackPlus )
		ME["Extra_CSCActivityWithMET"]->Fill(3,1 );
	      else if( CSCRecHitMinus && CSCSegmentPlus )
		ME["Extra_CSCActivityWithMET"]->Fill(3,2 );
	      else if( CSCRecHitMinus && CSCRecHitPlus ) 
		ME["Extra_CSCActivityWithMET"]->Fill(3,3);
	      else if( CSCRecHitMinus )
		ME["Extra_CSCActivityWithMET"]->Fill(3,4);
	      else ME["Extra_CSCActivityWithMET"]->Fill(4,4);
	    }
	}
    }
  
}

void BeamHaloAnalyzer::endJob()
{

}

BeamHaloAnalyzer::~BeamHaloAnalyzer(){
}

//DEFINE_FWK_MODULE(CMSEventAnalyzer);




