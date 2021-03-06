// Class:      L25and3TauEfficiencyAnalyzer
// Original Author:  Eduardo Luiggi, modified by Sho Maruyama
//         Created:  Fri Apr  4 16:37:44 CDT 2008
// $Id: L25and3TauEfficiencyAnalyzer.cc,v 1.18 2012/02/02 07:00:09 mkortela Exp $
#include "ElectroWeakAnalysis/TauTriggerEfficiency/interface/L25and3TauEfficiencyAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

#include "RecoTauTag/TauTagTools/interface/TauTagTools.h"

using namespace edm;
using namespace reco;
using namespace std;
L25and3TauEfficiencyAnalyzer::L25and3TauEfficiencyAnalyzer(){}

L25and3TauEfficiencyAnalyzer::L25and3TauEfficiencyAnalyzer(const edm::ParameterSet& iConfig){}

L25and3TauEfficiencyAnalyzer::~L25and3TauEfficiencyAnalyzer(){}

void L25and3TauEfficiencyAnalyzer::Setup(const edm::ParameterSet& iConfig,TTree* l25tree){
  l25JetSource = iConfig.getParameter<InputTag>("l25JetSource");
  l25PtCutSource = iConfig.getParameter<InputTag>("l25PtCutSource");
  l3IsoSource  = iConfig.getParameter<InputTag>("l3IsoSource");
  l25MatchingCone = iConfig.getParameter<double>("l25MatchingCone");
  HLTPFTau = iConfig.getParameter<bool>("HLTPFTau");

  edm::ParameterSet qualityCuts = iConfig.getParameter<edm::ParameterSet>("l3IsoQualityCuts");
  filterMinTrackPt = qualityCuts.getParameter<double>("minTrackPt");
  filterMinPixelHits = qualityCuts.getParameter<unsigned>("minTrackPixelHits");
  filterMinTrackerHits = qualityCuts.getParameter<unsigned>("minTrackHits");
  filterMaxIP = qualityCuts.getParameter<double>("maxTransverseImpactParameter");
  filterMaxChi2 = qualityCuts.getParameter<double>("maxTrackChi2");
  filterMaxDeltaZ = qualityCuts.getParameter<double>("maxDeltaZ");
  filterMinGammaEt = qualityCuts.getParameter<double>("minGammaEt");

  l25tree->Branch("l25Eta", &l25Eta);
  l25tree->Branch("l25Phi", &l25Phi);
  l25tree->Branch("l25Et", &l25Et);
  l25tree->Branch("l25Pt", &l25Pt);
  l25tree->Branch("l25PtLdgLoose",&l25PtLdgLoose);
  l25tree->Branch("l25EtaLdgLoose",&l25EtaLdgLoose);
  l25tree->Branch("l25PhiLdgLoose",&l25PhiLdgLoose);
  l25tree->Branch("l25PtLdgJetDRLoose",&l25PtLdgJetDRLoose);
  l25tree->Branch("l25InvPt", &l25InvPt);
  l25tree->Branch("matchedToL2Jet", &matchedToL2Jet);
  l25tree->Branch("foundTracksInJet", &foundTracksInJet);
  l25tree->Branch("tauTagHasTracks", &tauTagHasTracks);
  l25tree->Branch("matchedToHLTPFTau", &matchedToHLTPFTau);
  l25tree->Branch("l25PFTauLeadTrkIsValid", &l25PFTauLeadTrkIsValid);
  l25tree->Branch("leadDepth1", &leadDepth1);
  l25tree->Branch("leadDepth2", &leadDepth2);
  l25tree->Branch("minDR", &minDR);
  l25tree->Branch("bareEt", &bareEt);
  l25tree->Branch("l25DefDisc_Trk5_IsoPtMin1_Ntrk0",&l25DefDisc_Trk5_IsoPtMin1_Ntrk0);
  l25tree->Branch("l25DefDisc_Trk5_IsoPtMin1_Ntrk1",&l25DefDisc_Trk5_IsoPtMin1_Ntrk1);
  l25tree->Branch("l25DefDisc_Trk5_IsoPtMin1_5_Ntrk0",&l25DefDisc_Trk5_IsoPtMin1_5_Ntrk0);
  l25tree->Branch("l25DefDisc_Trk5_IsoPtMin1_5_Ntrk1",&l25DefDisc_Trk5_IsoPtMin1_5_Ntrk1);
  l25tree->Branch("l25DefDisc_Trk5_IsoPtMin2_Ntrk0",&l25DefDisc_Trk5_IsoPtMin2_Ntrk0);
  l25tree->Branch("l25IsoPtSum",&l25IsoPtSum);
  l25tree->Branch("l25IsoTrkNHits",&l25IsoTrkNHits);
  l25tree->Branch("l25IsoTrkChi2",&l25IsoTrkChi2);
  l25tree->Branch("l25IsoTrkPt",&l25IsoTrkPt);
  l25tree->Branch("l25TrkIsoPtSum",&l25TrkIsoPtSum);
  l25tree->Branch("l25EcalIsoEtSum",&l25EcalIsoEtSum);
  l25tree->Branch("l25TrkIsoPtMax",&l25TrkIsoPtMax);
  l25tree->Branch("l25EcalIsoEtMax",&l25EcalIsoEtMax);
  l25tree->Branch("l25TrkIsoPtMaxAll",&l25TrkIsoPtMaxAll);
  l25tree->Branch("l25EcalIsoEtMaxAll",&l25EcalIsoEtMaxAll);
  l25tree->Branch("l25NTrksIso",&l25NTrksIso);
  l25tree->Branch("l25NGammaIso",&l25NGammaIso);
  l25tree->Branch("l25NTrksIsoAll",&l25NTrksIsoAll);
  l25tree->Branch("l25NGammaIsoAll",&l25NGammaIsoAll);
  l25tree->Branch("l25Prong", &l25Prong);
  l25tree->Branch("l25PFTauEt", &l25PFTauEt);
  l25tree->Branch("l25PFTauPt", &l25PFTauPt);
  l25tree->Branch("primaryVertexIsValid", &primaryVertexIsValid);
  NMatchedToL2 = 0;
  NJetsWithTracks = 0;
}

void L25and3TauEfficiencyAnalyzer::fill(const edm::Event& iEvent, const reco::Particle& tau) {
  fill(iEvent,tau.p4());
}

void L25and3TauEfficiencyAnalyzer::fill(const edm::Event& iEvent, const reco::Candidate& tau) {
  fill(iEvent,tau.p4());
}

void L25and3TauEfficiencyAnalyzer::fill(const edm::Event& iEvent, const LorentzVector& tau) {
  // PF specific quantities are moved to TTEff analyzer.
  //this was originally the calotau method, but since it only uses direction from the calojet
  //make this method the LorentzVector one and add a CaloTau stub that calls this method - gfball
  
  //modified to collect extra discriminator parameters
  l25Et = 0;
  l25Phi = -10.;
  l25Eta = -10.;
  l25Pt = 0;
  l25InvPt = 0;
  l25PtLdgJetDRLoose = 0;
  l25EtaLdgLoose = 0;
  l25PhiLdgLoose = 0;
  l25PtLdgLoose = 0;
  
  l25DefDisc_Trk5_IsoPtMin1_Ntrk0 = 0;
  l25DefDisc_Trk5_IsoPtMin1_Ntrk1 = 0;
  l25DefDisc_Trk5_IsoPtMin1_5_Ntrk0 = 0;
  l25DefDisc_Trk5_IsoPtMin1_5_Ntrk1 = 0;
  l25DefDisc_Trk5_IsoPtMin2_Ntrk0 = 0;
 
  l25IsoPtSum = 0;
  l25EcalIsoEtSum = 0;
  l25TrkIsoPtMax = 0;
  l25EcalIsoEtMax = 0;
  l25TrkIsoPtMaxAll = 0;
  l25EcalIsoEtMaxAll = 0;
  l25NTrksIso = 0;
  l25NGammaIso = 0;
  l25NTrksIsoAll = 0;
  l25NGammaIsoAll = 0;
  l25Prong = 0;
  l25PFTauEt = 0;
  l25PFTauPt = 0;
  matchedToL2Jet = false;
  foundTracksInJet = false;
  matchedToHLTPFTau = false;
  tauTagHasTracks = false;
  l25PFTauLeadTrkIsValid = false;
  primaryVertexIsValid = false;

  if(HLTPFTau == true){
    Handle<PFTauTagInfoCollection> tags;
    iEvent.getByLabel(l25JetSource, tags);
    
    Handle<PFTauCollection> ptJets;
    iEvent.getByLabel(l25PtCutSource, ptJets); // Leading Pt Cut > X GeV/c applied, check HLT Config file
    
    Handle<TrackCollection> theHLTTracks;
    iEvent.getByLabel("hltPFJetCtfWithMaterialTracks", theHLTTracks);
    
    Handle<TrackCollection> theHQTracks;
    iEvent.getByLabel("hltPFlowTrackSelectionHighPurity", theHQTracks);
    
    Handle<JetTracksAssociationCollection> theJetTracks;
    iEvent.getByLabel("hltPFTauJetTracksAssociator", theJetTracks);
    
    Handle<VertexCollection> thePrimaryVertices;
    iEvent.getByLabel("hltPixelVertices", thePrimaryVertices);
    const reco::Vertex& theHltPrimaryVertex = (*thePrimaryVertices->begin());
    math::XYZPoint theVertexPosition = math::XYZPoint(0.,0.,0.);
    if(thePrimaryVertices->size() > 0) {
      theVertexPosition = math::XYZPoint(theHltPrimaryVertex.position().x(),
                                         theHltPrimaryVertex.position().y(),
                                         theHltPrimaryVertex.position().z());  
      primaryVertexIsValid = true;
    }
    //Handle<PFTauCollection> isoJets;
    //iEvent.getByLabel(l3IsoSource, isoJets);
    
    float minPFTauTagDeltaR = l25MatchingCone;
    PFTauTagInfo theMatchedPFTauTag;
    unsigned int nTauTags = 0;
    if(tags.isValid()){
      for(unsigned int j = 0; j < tags->size(); j++){ // bare L2 Taus
	nTauTags++;
	if(deltaR(tau, *(tags->at(j).pfjetRef())) < minPFTauTagDeltaR){ // dr < l25MatchingCone
	  minPFTauTagDeltaR = deltaR(tau, *(tags->at(j).pfjetRef()));
	  theMatchedPFTauTag = tags->at(j);
	  matchedToL2Jet = true;
	  NMatchedToL2++;
	}
      }
    }  
    else {
      edm::LogWarning("TTEffAnalyzer") <<"No L25Jetsource in event!"<<std::endl;
    }
    if(matchedToL2Jet){ //Found the closest pfTauTagHLT
      if(!theMatchedPFTauTag.pfjetRef().isNonnull()) cout <<"Invalid PFTauTag Ref!"<<std::endl;
      l25Eta  = theMatchedPFTauTag.pfjetRef()->eta();	   
      l25Phi  = theMatchedPFTauTag.pfjetRef()->phi();	   
      l25Et   = theMatchedPFTauTag.pfjetRef()->et();	   
      const PFCandidateRefVector chargedHads = theMatchedPFTauTag.PFChargedHadrCands(); 								 
      for(unsigned int k = 0; k < chargedHads.size(); k++){											 
        if(chargedHads.at(k)->pt() > 0.0 && deltaR(*(theMatchedPFTauTag.pfjetRef()),*(chargedHads.at(k)) ) < 0.8)  foundTracksInJet = true;
      } // charged hadron loop
    }
    float minPFTauDeltaR = l25MatchingCone;
    PFTau theMatchedHLTPFtau;
    if(ptJets.isValid()){ // Leading Pt Cut > X GeV/c applied, check HLT Config file
      for(unsigned int j = 0; j < ptJets->size(); j++){
 	if(deltaR(tau, ptJets->at(j)) < minPFTauDeltaR){ // dr < l25MatchingCone
           matchedToHLTPFTau = true;
	   minPFTauDeltaR = deltaR(tau, ptJets->at(j));
	   theMatchedHLTPFtau = ptJets->at(j);
	}
      }
    }
    else {
      edm::LogWarning("TTEffAnalyzer") <<"No L25PtJetsource in event!"<<std::endl;
    }
    if(matchedToHLTPFTau){ //Got matching pftauHLT
      const PFCandidateRef theLeadPFChargedCand = theMatchedHLTPFtau.leadPFChargedHadrCand();
      if(theLeadPFChargedCand.isNonnull()){
        NJetsWithTracks++;
        l25PFTauLeadTrkIsValid = 1;
        l25Pt = theLeadPFChargedCand->pt();
        l25InvPt = 1./l25Pt;
        l25Prong = theMatchedHLTPFtau.signalPFChargedHadrCands().size();
        l25PFTauEt = theMatchedHLTPFtau.et();
        l25PFTauPt = theMatchedHLTPFtau.pt();
        l25IsoPtSum = theMatchedHLTPFtau.isolationPFChargedHadrCandsPtSum();
        l25EcalIsoEtSum = theMatchedHLTPFtau.isolationPFGammaCandsEtSum();
	l25NTrksIsoAll = theMatchedHLTPFtau.isolationPFChargedHadrCands().size();
        l25NGammaIsoAll = theMatchedHLTPFtau.isolationPFGammaCands().size();

        for(size_t i=0; i<theMatchedHLTPFtau.isolationPFChargedHadrCands().size(); ++i) {
          l25TrkIsoPtMaxAll = std::max(l25TrkIsoPtMaxAll, static_cast<float>(theMatchedHLTPFtau.isolationPFChargedHadrCands()[i]->pt()));
        }
        for(size_t i=0; i<theMatchedHLTPFtau.isolationPFGammaCands().size(); ++i) {
          l25EcalIsoEtMaxAll = std::max(l25EcalIsoEtMaxAll, static_cast<float>(theMatchedHLTPFtau.isolationPFGammaCands()[i]->et()));
        }

        if(primaryVertexIsValid) {
          // the number of candidates after quality filtering
          reco::PFCandidateRefVector cands;
          cands = TauTagTools::filteredPFChargedHadrCands(theMatchedHLTPFtau.isolationPFChargedHadrCands(),
                                                          filterMinTrackPt,
                                                          filterMinPixelHits,
                                                          filterMinTrackerHits,
                                                          filterMaxIP,
                                                          filterMaxChi2,
                                                          filterMaxDeltaZ,
                                                          theHltPrimaryVertex,
                                                          theHltPrimaryVertex.position().z());
          l25NTrksIso = cands.size();
        
          cands = TauTagTools::filteredPFGammaCands(theMatchedHLTPFtau.isolationPFGammaCands(),
                                                    filterMinGammaEt);
          l25NGammaIso = cands.size();

          // Maximum pt/et of candidates after quality filtering
          cands = TauTagTools::filteredPFChargedHadrCands(theMatchedHLTPFtau.isolationPFChargedHadrCands(),
                                                          0,
                                                          filterMinPixelHits,
                                                          filterMinTrackerHits,
                                                          filterMaxIP,
                                                          filterMaxChi2,
                                                          filterMaxDeltaZ,
                                                          theHltPrimaryVertex,
                                                          theHltPrimaryVertex.position().z());
          for(size_t i=0; i<cands.size(); ++i) {
            l25TrkIsoPtMax = std::max(l25TrkIsoPtMax, static_cast<float>(cands[i]->pt()));
          }

          cands = TauTagTools::filteredPFGammaCands(theMatchedHLTPFtau.isolationPFGammaCands(), 0);
          for(size_t i=0; i<cands.size(); ++i) {
            l25EcalIsoEtMax = std::max(l25EcalIsoEtMax, static_cast<float>(cands[i]->pt()));
          }
        }
      }
      
      else{
/*
      	cout.precision(3);									    
        cout << "Number of TauTag Jets\t" << nTauTags << endl;
        cout << "TauTagInfoTracks \n";
        const PFCandidateRefVector chargedHads = theMatchedPFTauTag.PFChargedHadrCands(); 								 
        for(unsigned int trkIt = 0; trkIt < chargedHads.size(); trkIt++){											 
     	  cout << trkIt << "\t Trk Pt " << chargedHads.at(trkIt)->pt()
     	       << "\t DeltaR " << deltaR(*(theMatchedPFTauTag.pfjetRef()),*(chargedHads.at(trkIt)))
    	       << "\t Chi2 " << chargedHads.at(trkIt)->trackRef()->chi2()
     	       << "\t Chi2/NdF " << chargedHads.at(trkIt)->trackRef()->normalizedChi2()
     	       << "\t Hits " << chargedHads.at(trkIt)->trackRef()->numberOfValidHits()
     	       << "\t PxlHits " << chargedHads.at(trkIt)->trackRef()->hitPattern().numberOfValidPixelHits()
     	       << "\t Vtx " << theVertexPosition
     	       << "\t dxy " << chargedHads.at(trkIt)->trackRef()->dxy(theVertexPosition)
     	       << "\t dz " << chargedHads.at(trkIt)->trackRef()->dz(theVertexPosition)
     	       << "\n";
        }
	
      	//if no tracks with pt > 1 and DR < 0.2 Look at all tracks in jet
      	cout << "hltPFJetCtfWithMaterialTracks \n";						    
      	for(unsigned int trkIt = 0; trkIt < theHLTTracks->size(); trkIt++){			    
     	  cout << trkIt << "\t Trk Pt " << theHLTTracks->at(trkIt).pt()
     	       << "\t DeltaR " << deltaR(*(theMatchedPFTauTag.pfjetRef()),theHLTTracks->at(trkIt))
     	       << "\t Chi2 " << theHLTTracks->at(trkIt).chi2()
     	       << "\t Chi2/NdF " << theHLTTracks->at(trkIt).normalizedChi2()
     	       << "\t Hits " << theHLTTracks->at(trkIt).numberOfValidHits()
     	       << "\t PxlHits " << theHLTTracks->at(trkIt).hitPattern().numberOfValidPixelHits()
     	       << "\t Vtx " << theVertexPosition
     	       << "\t dxy " << theHLTTracks->at(trkIt).dxy(theVertexPosition)
     	       << "\t dz " << theHLTTracks->at(trkIt).dz(theVertexPosition)
     	       << "\n";
      	}
      	cout << "hltPFlowTrackSelectionHighPurity \n";
      	for(unsigned int trkIt = 0; trkIt < theHQTracks->size(); trkIt++){			    
     	  cout << trkIt << "\t Trk Pt " << theHQTracks->at(trkIt).pt()
     	       << "\t DeltaR " << deltaR(*(theMatchedPFTauTag.pfjetRef()),theHQTracks->at(trkIt))
     	       << "\t Chi2 " << theHQTracks->at(trkIt).chi2()
     	       << "\t Chi2/NdF " << theHQTracks->at(trkIt).normalizedChi2()
     	       << "\t Hits " << theHQTracks->at(trkIt).numberOfValidHits()
     	       << "\t PxlHits " << theHQTracks->at(trkIt).hitPattern().numberOfValidPixelHits()
     	       << "\t Vtx " << theVertexPosition
     	       << "\t dxy " << theHQTracks->at(trkIt).dxy(theVertexPosition)
     	       << "\t dz " << theHQTracks->at(trkIt).dz(theVertexPosition)
     	       << "\n";
      	}
*/
      }
 
/*
    
      if(isoJets.isValid()){
        for(unsigned int j = 0; j < isoJets->size(); j++){
          if(deltaR(tau, isoJets->at(j)) < l25MatchingCone){ // dr < l25MatchingCone
            if(l25Depth < 4){
              l25Depth = 4; // iso match
              break;
            }
          }
        }
      }

      else {
        //edm::LogWarning("TTEffAnalyzer") <<"No L3IsoJetsource in event!"<<std::endl;
      }
    
*/   
    }
  }
  else {

    Handle<IsolatedTauTagInfoCollection> tags;
    iEvent.getByLabel(l25JetSource, tags);
/*    
    Handle<CaloJetCollection> ptJets;
    Handle<CaloJetCollection> isoJets;

    iEvent.getByLabel(l25PtCutSource, ptJets); // Leading Pt Cut > X GeV/c applied, check HLT Config file
    iEvent.getByLabel(l3IsoSource, isoJets);
*/    
    //float minDeltaR = l25MatchingCone;
    //IsolatedTauTagInfo theMatchedTauTag;
    if(tags.isValid()){
      for(unsigned int j = 0; j < tags->size(); j++){ // bare L2 Taus
       // find the closest matched tauTagInfo object to offline tau
	if(deltaR(tau, *(tags->at(j).jet())) < l25MatchingCone){ // dr < l25MatchingCone
	  matchedToL2Jet = true; // L2 match
	  if(tags->at(j).allTracks().size() > 0) foundTracksInJet = true;					         
          l25Eta  = tags->at(j).jet()->eta();							         
          l25Phi  = tags->at(j).jet()->phi();							         
      	  l25Et   = tags->at(j).jet()->et();
          const TrackRef leadTrk = tags->at(j).leadingSignalTrack(0.2,1.0);// track finding cone = 0.2   
           if(leadTrk.isNonnull()){  		 					     
             l25PFTauLeadTrkIsValid = true; 	 							     
             l25Pt = leadTrk->pt();  		 								  
             l25InvPt = 1./leadTrk->pt();
             
	     l25DefDisc_Trk5_IsoPtMin1_Ntrk0 = tags->at(j).discriminator(0.2,0.15,0.5,5.,1.,0,0.2);
             l25DefDisc_Trk5_IsoPtMin1_Ntrk1 = tags->at(j).discriminator(0.2,0.15,0.5,5.,1.,1,0.2);	   
             l25DefDisc_Trk5_IsoPtMin1_5_Ntrk0 = tags->at(j).discriminator(0.2,0.15,0.5,5.,1.5,0,0.2);	    
             l25DefDisc_Trk5_IsoPtMin1_5_Ntrk1 = tags->at(j).discriminator(0.2,0.15,0.5,5.,1.5,1,0.2);	    
             l25DefDisc_Trk5_IsoPtMin2_Ntrk0 = tags->at(j).discriminator(0.2,0.15,0.5,5.,2.0,0,0.2);	 
             l25IsoPtSum = 0.;  											  
             const TrackRefVector theSignalTracks = tags->at(j).tracksInCone(leadTrk->momentum(), 0.15, 1.0);	     
             const TrackRefVector theIsoTracks = tags->at(j).tracksInCone(leadTrk->momentum(), 0.5, 1.0);  		
             l25NTrksIso = theIsoTracks.size() - theSignalTracks.size();						  
             														
             for(TrackRefVector::const_iterator isoIt = theIsoTracks.begin(); isoIt != theIsoTracks.end(); ++isoIt){	  
               if(deltaR(leadTrk->momentum(), (*isoIt)->momentum()) > 0.15 && deltaR(leadTrk->momentum(), (*isoIt)->momentum()) < 0.5){
             	 l25IsoPtSum += (*isoIt)->pt(); 									  
             	 l25IsoTrkChi2 = (*isoIt)->chi2();									  
             	 l25IsoTrkPt = (*isoIt)->pt();  									  
             	 l25IsoTrkNHits = (*isoIt)->numberOfValidHits();							  
               }													  
             }  													  
	   }	 								  
        } 
      }  
    }// non empty collection
    else{
      edm::LogWarning("TTEffAnalyzer") <<"No L25Jetsource in event!"<<std::endl;
    }
/*    
    if(l25Depth > 0 && theMatchedTauTag.jet().isNonnull()){ // found a matched hlt object; now get the info
      // Use the matched tauTagInfo object to extract the needed info  
      l25Eta  = theMatchedTauTag.jet()->eta();       		       
      l25Phi  = theMatchedTauTag.jet()->phi();       		       
      l25Et   = theMatchedTauTag.jet()->et();	     		       
      const TrackRef leadTrk = theMatchedTauTag.leadingSignalTrack(0.2,1.0);// track finding cone = 0.2
      if(leadTrk.isNonnull()){ 									
        if(l25Depth < 2) l25Depth = 2; 										
        l25Pt = leadTrk->pt();  								      		     
        l25InvPt = 1./leadTrk->pt();							      			     
        													     
        // get the tracks in iso region around leading track 							     
        l25IsoPtSum = 0.;										   	     
        const TrackRefVector theSignalTracks = theMatchedTauTag.tracksInCone(leadTrk->momentum(), 0.15, 1.0); 	   	
        const TrackRefVector theIsoTracks = theMatchedTauTag.tracksInCone(leadTrk->momentum(), 0.5, 1.0);     	   	   
        l25NTrksIso = theIsoTracks.size() - theSignalTracks.size();					   	     
        													   
        for(TrackRefVector::const_iterator isoIt = theIsoTracks.begin(); isoIt != theIsoTracks.end(); ++isoIt){      
          if(deltaR(leadTrk->momentum(), (*isoIt)->momentum()) > 0.15 && deltaR(leadTrk->momentum(), (*isoIt)->momentum()) < 0.5){
            l25IsoPtSum += (*isoIt)->pt();								   	     
            l25IsoTrkChi2 = (*isoIt)->chi2();								   	     
            l25IsoTrkPt = (*isoIt)->pt();								   	     
            l25IsoTrkNHits = (*isoIt)->numberOfValidHits();						   	     
          }													     
        }													     
        //evaluate a series of different discriminator parameters				      		     
        // MatchConeSize, sigCone, isoCone, ltPt, pt min, nTracksIsoRing, ltDz  		      		     
        if(theMatchedTauTag.discriminator(0.2,0.15,0.5,5.,1.,0,0.2))l25DefDisc_Trk5_IsoPtMin1_Ntrk0=1;     	     
        if(theMatchedTauTag.discriminator(0.2,0.15,0.5,5.,1.,1,0.2))l25DefDisc_Trk5_IsoPtMin1_Ntrk1=1;     	     
        if(theMatchedTauTag.discriminator(0.2,0.15,0.5,5.,1.5,0,0.2))l25DefDisc_Trk5_IsoPtMin1_5_Ntrk0=1;  	     
        if(theMatchedTauTag.discriminator(0.2,0.15,0.5,5.,1.5,1,0.2))l25DefDisc_Trk5_IsoPtMin1_5_Ntrk1=1;  	     
        if(theMatchedTauTag.discriminator(0.2,0.15,0.5,5.,2.0,0,0.2))l25DefDisc_Trk5_IsoPtMin2_Ntrk0=1;    	     
        													     
        //Find lead track within a very large cone...								     
        const TrackRef leadTrkLoose = theMatchedTauTag.leadingSignalTrack(100.,1.0);// Loose Track finding  	     
        if(leadTrkLoose.isNonnull()){								        	     
          l25PtLdgLoose = leadTrkLoose->pt();							        	     
          l25EtaLdgLoose = leadTrkLoose->eta(); 						        	     
          l25PhiLdgLoose = leadTrkLoose->phi(); 						        	     
          double dphi = fabs(l25PhiLdgLoose-l25Phi);						        	     
          if(dphi>2*acos(-1.))dphi=2*acos(-1.)-dphi;						        	     
          double deta = fabs(l25EtaLdgLoose-l25Eta);						        	     
          l25PtLdgJetDRLoose = sqrt(dphi*dphi+deta*deta);					        	     
	}
      }// good lead cand 												
*/    


/*      
      if(ptJets.isValid()){ // Leading Pt Cut > X GeV/c applied, check HLT Config file
        for(unsigned int j = 0; j < ptJets->size(); j++){
          if(deltaR(tau, ptJets->at(j) ) < l25MatchingCone){ // dr < l25MatchingCone
            if(l25Depth < 3){
              l25Depth = 3; // lead pt cut match
              break;
            }
          }// pf and l25 tau match dr < l25MatchingCone
        }// for jet loop
      }// non empty collection
      else {
        //edm::LogWarning("TTEffAnalyzer") <<"No L25PtJetsource in event!"<<std::endl;
      }
      if(isoJets.isValid()){
        for(unsigned int j = 0; j < isoJets->size(); j++){
          if(deltaR(tau, isoJets->at(j)) < l25MatchingCone){ // dr < l25MatchingCone
            if(l25Depth < 4){
              l25Depth = 4; // iso match
              break;
            }
          }
        }
      }
      else {
        //edm::LogWarning("TTEffAnalyzer") <<"No L3IsoJetsource in event!"<<std::endl;
      } 
    }
*/    
  }
}// tau ends

void L25and3TauEfficiencyAnalyzer::beginJob(const edm::EventSetup&){} 
void L25and3TauEfficiencyAnalyzer::analyze(const edm::Event&, const edm::EventSetup&){}
void L25and3TauEfficiencyAnalyzer::endJob(){
  cout << "\n\tTracksInJets/MatchedJets " << (double)NJetsWithTracks/(double)NMatchedToL2 << endl;
} 
