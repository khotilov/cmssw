 //
// Original Authors: Fabian Stoeckli: fabian.stoeckli@cern.ch
//                   Nicholas Wardle: nckw@cern.ch
//                   Rishi Patel rpatel@cern.ch
//

#include "RecoParticleFlow/PFProducer/interface/PFPhotonAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementSuperCluster.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TFile.h>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace reco;


PFPhotonAlgo::PFPhotonAlgo(std::string mvaweightfile,  
			   double mvaConvCut, 
			   bool useReg,
			   std::string mvaWeightFilePFClusCorr, 
			   std::string mvaWeightFilePFPhoCorr, 
			   std::string X0_Map,
			   const reco::Vertex& primary,
			   const boost::shared_ptr<PFEnergyCalibration>& thePFEnergyCalibration,
                           double sumPtTrackIsoForPhoton,
                           double sumPtTrackIsoSlopeForPhoton
			   ) : 
  isvalid_(false), 
  verbosityLevel_(Silent), 
  MVACUT(mvaConvCut),
  useReg_(useReg),
  thePFEnergyCalibration_(thePFEnergyCalibration),
  sumPtTrackIsoForPhoton_(sumPtTrackIsoForPhoton),
  sumPtTrackIsoSlopeForPhoton_(sumPtTrackIsoSlopeForPhoton)
{  
    primaryVertex_=primary;  
    //Book MVA  
    tmvaReader_ = new TMVA::Reader("!Color:Silent");  
    tmvaReader_->AddVariable("del_phi",&del_phi);  
    tmvaReader_->AddVariable("nlayers", &nlayers);  
    tmvaReader_->AddVariable("chi2",&chi2);  
    tmvaReader_->AddVariable("EoverPt",&EoverPt);  
    tmvaReader_->AddVariable("HoverPt",&HoverPt);  
    tmvaReader_->AddVariable("track_pt", &track_pt);  
    tmvaReader_->AddVariable("STIP",&STIP);  
    tmvaReader_->AddVariable("nlost", &nlost);  
    tmvaReader_->BookMVA("BDT",mvaweightfile.c_str());  
    if(useReg_)
      {
	tmvaLCRegReader_=new TMVA::Reader("!Color:Silent");
	tmvaLCRegReader_->AddVariable("genZ", &VtxZ_);
	tmvaLCRegReader_->AddVariable("PFPhoCluseta[0]/abs(PFPhoCluseta[0])", &EB);
	tmvaLCRegReader_->AddVariable("abs(PFPhoCluseta[0])", &ClusEta_);
	tmvaLCRegReader_->AddVariable("PFPhoClusphi[0]", &ClusPhi_);
	tmvaLCRegReader_->AddVariable("log(PFPhoClusE[0])", &logPFClusE_);
	tmvaLCRegReader_->AddVariable("PFPhoClusE3x3[0]/PFPhoClusE[0]", &ClusR9_);
	tmvaLCRegReader_->AddVariable("PFPhoClusE5x5[0]/PFPhoClusE[0]", &Clus5x5ratio_);   
	tmvaLCRegReader_->AddVariable("PFCrysPhiCrack[0]", &PFCrysPhiCrack_);
	tmvaLCRegReader_->AddVariable("PFCrysEtaCrack[0]", &PFCrysEtaCrack_);
	tmvaLCRegReader_->AddVariable("PFCrysEta[0]", &CrysEta_);
	tmvaLCRegReader_->AddVariable("PFCrysPhi[0]", &CrysPhi_);
	tmvaLCRegReader_->BookMVA("BDTG_LocCorr",mvaWeightFilePFClusCorr.c_str()); 
	
	tmvaGCRegReader_=new TMVA::Reader("!Color:Silent");   
	tmvaGCRegReader_->AddVariable("PFPhoEta", &PFPhoEta_);
	tmvaGCRegReader_->AddVariable("PFPhoEtCorr", &PFPhoEt_);
	tmvaGCRegReader_->AddVariable("PFPhoR9", &PFPhoR9_);
	tmvaGCRegReader_->AddVariable("PFPhoPhi", &PFPhoPhi_);
	tmvaGCRegReader_->AddVariable("SCEtaWidth", &SCEtaWidth_);
	tmvaGCRegReader_->AddVariable("SCPhiWidth", &SCPhiWidth_);
	tmvaGCRegReader_->AddVariable("X0_inner", &x0inner_);
	tmvaGCRegReader_->AddVariable("X0_middle", &x0middle_);
	tmvaGCRegReader_->AddVariable("X0_outer", &x0outer_);
	tmvaGCRegReader_->AddVariable("nPFClus", &nPFClus_);
	tmvaGCRegReader_->AddVariable("Rconv", &RConv_);
	tmvaGCRegReader_->AddVariable("PFPhoClusLowE", &LowClusE_);
	tmvaGCRegReader_->AddVariable("dEtaLow", &dEta_);
	tmvaGCRegReader_->AddVariable("dPhiLow", &dPhi_);
	tmvaGCRegReader_->AddVariable("excluded", &excluded_);
	tmvaGCRegReader_->AddVariable("Mustache_Et_out/(Mustache_Et_in+Mustache_Et_out)", &Mustache_EtRatio_);
	
	tmvaGCRegReader_->BookMVA("BDTG_GCorr",mvaWeightFilePFPhoCorr.c_str());
    //Material Map
	TFile *XO_File = new TFile(X0_Map.c_str(),"READ");
	X0_sum=(TH2D*)XO_File->Get("TrackerSum");
	X0_inner = (TH2D*)XO_File->Get("Inner");
	X0_middle = (TH2D*)XO_File->Get("Middle");
	X0_outer = (TH2D*)XO_File->Get("Outer");
      }
}

void PFPhotonAlgo::RunPFPhoton(const reco::PFBlockRef&  blockRef,
			       std::vector<bool>& active,
			       std::auto_ptr<PFCandidateCollection> &pfCandidates,
			       std::vector<reco::PFCandidatePhotonExtra>& pfPhotonExtraCandidates,
			       std::vector<reco::PFCandidate> 
			       &tempElectronCandidates
){
  
  //std::cout<<" calling RunPFPhoton "<<std::endl;
  
  /*      For now we construct the PhotonCandidate simply from 
	  a) adding the CORRECTED energies of each participating ECAL cluster
	  b) build the energy-weighted direction for the Photon
  */


  // define how much is printed out for debugging.
  // ... will be setable via CFG file parameter
  verbosityLevel_ = Chatty;          // Chatty mode.
  

  // loop over all elements in the Block
  const edm::OwnVector< reco::PFBlockElement >&          elements         = blockRef->elements();
  edm::OwnVector< reco::PFBlockElement >::const_iterator ele              = elements.begin();
  std::vector<bool>::const_iterator                      actIter          = active.begin();
  PFBlock::LinkData                                      linkData         = blockRef->linkData();
  bool                                                   isActive         = true;


  if(elements.size() != active.size()) {
    // throw excpetion...
    //std::cout<<" WARNING: Size of collection and active-vectro don't agree!"<<std::endl;
    return;
  }
  
  // local vecotr to keep track of the indices of the 'elements' for the Photon candidate
  // once we decide to keep the candidate, the 'active' entriesd for them must be set to false
  std::vector<unsigned int> elemsToLock;
  elemsToLock.resize(0);
  
  for( ; ele != elements.end(); ++ele, ++actIter ) {

    // if it's not a SuperCluster, go to the next element
    if( !( ele->type() == reco::PFBlockElement::SC ) ) continue;
    
    // Photon kienmatics, will be updated for each identified participating element
    float photonEnergy_        =   0.;
    float photonX_             =   0.;
    float photonY_             =   0.;
    float photonZ_             =   0.;
    float RawEcalEne           =   0.;

    // Total pre-shower energy
    float ps1TotEne      = 0.;
    float ps2TotEne      = 0.;
    
    bool hasConvTrack=false;  
    bool hasSingleleg=false;  
    std::vector<unsigned int> AddClusters(0);  
    std::vector<unsigned int> IsoTracks(0);  
    std::multimap<unsigned int, unsigned int>ClusterAddPS1;  
    std::multimap<unsigned int, unsigned int>ClusterAddPS2;
    std::vector<reco::TrackRef>singleLegRef;
    std::vector<float>MVA_values(0);
    reco::ConversionRefVector ConversionsRef_;
    isActive = *(actIter);
    //cout << " Found a SuperCluster.  Energy " ;
    const reco::PFBlockElementSuperCluster *sc = dynamic_cast<const reco::PFBlockElementSuperCluster*>(&(*ele));
    //std::cout << sc->superClusterRef()->energy () << " Track/Ecal/Hcal Iso " << sc->trackIso()<< " " << sc->ecalIso() ;
    //std::cout << " " << sc->hcalIso() <<std::endl;
    if (!(sc->fromPhoton()))continue;

    // check the status of the SC Element... 
    // ..... I understand it should *always* be active, since PFElectronAlgo does not touch this (yet?) RISHI: YES
    if( !isActive ) {
      //std::cout<<" SuperCluster is NOT active.... "<<std::endl;
      continue;
    }
    elemsToLock.push_back(ele-elements.begin()); //add SC to elements to lock
    // loop over its constituent ECAL cluster
    std::multimap<double, unsigned int> ecalAssoPFClusters;
    blockRef->associatedElements( ele-elements.begin(), 
				  linkData,
				  ecalAssoPFClusters,
				  reco::PFBlockElement::ECAL,
				  reco::PFBlock::LINKTEST_ALL );
    
    // loop over the ECAL clusters linked to the iEle 
    if( ! ecalAssoPFClusters.size() ) {
      // This SC element has NO ECAL elements asigned... *SHOULD NOT HAPPEN*
      //std::cout<<" Found SC element with no ECAL assigned "<<std::endl;
      continue;
    }
    
    // This is basically CASE 2
    // .... we loop over all ECAL cluster linked to each other by this SC
    for(std::multimap<double, unsigned int>::iterator itecal = ecalAssoPFClusters.begin(); 
	itecal != ecalAssoPFClusters.end(); ++itecal) { 
      
      // to get the reference to the PF clusters, this is needed.
      reco::PFClusterRef clusterRef = elements[itecal->second].clusterRef();	
      
      // from the clusterRef get the energy, direction, etc
      //      float ClustRawEnergy = clusterRef->energy();
      //      float ClustEta = clusterRef->position().eta();
      //      float ClustPhi = clusterRef->position().phi();

      // initialize the vectors for the PS energies
      vector<double> ps1Ene(0);
      vector<double> ps2Ene(0);
      double ps1=0;  
      double ps2=0;  
      hasSingleleg=false;  
      hasConvTrack=false;

      /*
      cout << " My cluster index " << itecal->second 
	   << " energy " <<  ClustRawEnergy
	   << " eta " << ClustEta
	   << " phi " << ClustPhi << endl;
      */
      // check if this ECAL element is still active (could have been eaten by PFElectronAlgo)
      // ......for now we give the PFElectron Algo *ALWAYS* Shot-Gun on the ECAL elements to the PFElectronAlgo
      
      if( !( active[itecal->second] ) ) {
	//std::cout<< "  .... this ECAL element is NOT active anymore. Is skipped. "<<std::endl;
	continue;
      }
      
      // ------------------------------------------------------------------------------------------
      // TODO: do some tests on the ECAL cluster itself, deciding to use it or not for the Photons
      // ..... ??? Do we need this?
      if ( false ) {
	// Check if there are a large number tracks that do not pass pre-ID around this ECAL cluster
	bool useIt = true;
	int mva_reject=0;  
	bool isClosest=false;  
	std::multimap<double, unsigned int> Trackscheck;  
	blockRef->associatedElements( itecal->second,  
				      linkData,  
				      Trackscheck,  
				      reco::PFBlockElement::TRACK,  
				      reco::PFBlock::LINKTEST_ALL);  
	for(std::multimap<double, unsigned int>::iterator track = Trackscheck.begin();  
	    track != Trackscheck.end(); ++track) {  
	   
	  // first check if is it's still active  
	  if( ! (active[track->second]) ) continue;  
	  hasSingleleg=EvaluateSingleLegMVA(blockRef,  primaryVertex_, track->second);  
	  //check if it is the closest linked track  
	  std::multimap<double, unsigned int> closecheck;  
	  blockRef->associatedElements(track->second,  
				       linkData,  
				       closecheck,  
				       reco::PFBlockElement::ECAL,  
				       reco::PFBlock::LINKTEST_ALL);  
	  if(closecheck.begin()->second ==itecal->second)isClosest=true;  
	  if(!hasSingleleg)mva_reject++;  
	}  
	 
	if(mva_reject>0 &&  isClosest)useIt=false;  
	//if(mva_reject==1 && isClosest)useIt=false;
	if( !useIt ) continue;    // Go to next ECAL cluster within SC
      }
      // ------------------------------------------------------------------------------------------
      
      // We decided to keep the ECAL cluster for this Photon Candidate ...
      elemsToLock.push_back(itecal->second);
      
      // look for PS in this Block linked to this ECAL cluster      
      std::multimap<double, unsigned int> PS1Elems;
      std::multimap<double, unsigned int> PS2Elems;
      //PS Layer 1 linked to ECAL cluster
      blockRef->associatedElements( itecal->second,
				    linkData,
				    PS1Elems,
				    reco::PFBlockElement::PS1,
				    reco::PFBlock::LINKTEST_ALL );
      //PS Layer 2 linked to the ECAL cluster
      blockRef->associatedElements( itecal->second,
				    linkData,
				    PS2Elems,
				    reco::PFBlockElement::PS2,
				    reco::PFBlock::LINKTEST_ALL );
      
      // loop over all PS1 and compute energy
      for(std::multimap<double, unsigned int>::iterator iteps = PS1Elems.begin();
	  iteps != PS1Elems.end(); ++iteps) {

	// first chekc if it's still active
	if( !(active[iteps->second]) ) continue;
	
	//Check if this PS1 is not closer to another ECAL cluster in this Block          
	std::multimap<double, unsigned int> ECALPS1check;  
	blockRef->associatedElements( iteps->second,  
				      linkData,  
				      ECALPS1check,  
				      reco::PFBlockElement::ECAL,  
				      reco::PFBlock::LINKTEST_ALL );  
	if(itecal->second==ECALPS1check.begin()->second)//then it is closest linked  
	  {
	    reco::PFClusterRef ps1ClusterRef = elements[iteps->second].clusterRef();
	    ps1Ene.push_back( ps1ClusterRef->energy() );
	    ps1=ps1+ps1ClusterRef->energy(); //add to total PS1
	    // incativate this PS1 Element
	    elemsToLock.push_back(iteps->second);
	  }
      }
      for(std::multimap<double, unsigned int>::iterator iteps = PS2Elems.begin();
	  iteps != PS2Elems.end(); ++iteps) {

	// first chekc if it's still active
	if( !(active[iteps->second]) ) continue;
	
	// Check if this PS2 is not closer to another ECAL cluster in this Block:
	std::multimap<double, unsigned int> ECALPS2check;  
	blockRef->associatedElements( iteps->second,  
				      linkData,  
				      ECALPS2check,  
				      reco::PFBlockElement::ECAL,  
				      reco::PFBlock::LINKTEST_ALL );  
	if(itecal->second==ECALPS2check.begin()->second)//is closest linked  
	  {
	    reco::PFClusterRef ps2ClusterRef = elements[iteps->second].clusterRef();
	    ps2Ene.push_back( ps2ClusterRef->energy() );
	    ps2=ps2ClusterRef->energy()+ps2; //add to total PS2
	    // incativate this PS2 Element
	    elemsToLock.push_back(iteps->second);
	  }
      }
            
      // loop over the HCAL Clusters linked to the ECAL cluster (CASE 6)
      std::multimap<double, unsigned int> hcalElems;
      blockRef->associatedElements( itecal->second,linkData,
				    hcalElems,
				    reco::PFBlockElement::HCAL,
				    reco::PFBlock::LINKTEST_ALL );

      for(std::multimap<double, unsigned int>::iterator ithcal = hcalElems.begin();
	  ithcal != hcalElems.end(); ++ithcal) {

	if ( ! (active[ithcal->second] ) ) continue; // HCAL Cluster already used....
	
	// TODO: Decide if this HCAL cluster is to be used
	// .... based on some Physics
	// .... To we need to check if it's closer to any other ECAL/TRACK?

	bool useHcal = false;
	if ( !useHcal ) continue;
	//not locked
	//elemsToLock.push_back(ithcal->second);
      }

      // This is entry point for CASE 3.
      // .... we loop over all Tracks linked to this ECAL and check if it's labeled as conversion
      // This is the part for looping over all 'Conversion' Tracks
      std::multimap<double, unsigned int> convTracks;
      blockRef->associatedElements( itecal->second,
				    linkData,
				    convTracks,
				    reco::PFBlockElement::TRACK,
				    reco::PFBlock::LINKTEST_ALL);
      for(std::multimap<double, unsigned int>::iterator track = convTracks.begin();
	  track != convTracks.end(); ++track) {

	// first check if is it's still active
	if( ! (active[track->second]) ) continue;
	
	// check if it's a CONV track
	const reco::PFBlockElementTrack * trackRef = dynamic_cast<const reco::PFBlockElementTrack*>((&elements[track->second])); 	
	
	//Check if track is a Single leg from a Conversion  
	mvaValue=-999;  
	hasSingleleg=EvaluateSingleLegMVA(blockRef,  primaryVertex_, track->second);  

	// Daniele; example for mvaValues, do the same for single leg trackRef and convRef
	//          
	// 	if(hasSingleleg)
	// 	  mvaValues.push_back(mvaValue);

	//If it is not then it will be used to check Track Isolation at the end  
	if(!hasSingleleg)  
	  {  
	    bool included=false;  
	    //check if this track is already included in the vector so it is linked to an ECAL cluster that is already examined  
	    for(unsigned int i=0; i<IsoTracks.size(); i++)  
	      {if(IsoTracks[i]==track->second)included=true;}  
	    if(!included)IsoTracks.push_back(track->second);  
	  }  
	//For now only Pre-ID tracks that are not already identified as Conversions  
	if(hasSingleleg &&!(trackRef->trackType(reco::PFBlockElement::T_FROM_GAMMACONV)))  
	  {  
	    elemsToLock.push_back(track->second);
	    
	    reco::TrackRef t_ref=elements[track->second].trackRef();
	    bool matched=false;
	    for(unsigned int ic=0; ic<singleLegRef.size(); ic++)
	      if(singleLegRef[ic]==t_ref)matched=true;
	    
	    if(!matched){
	      singleLegRef.push_back(t_ref);
	      MVA_values.push_back(mvaValue);
	    }
	    //find all the clusters linked to this track  
	    std::multimap<double, unsigned int> moreClusters;  
	    blockRef->associatedElements( track->second,  
					  linkData,  
					  moreClusters,  
					  reco::PFBlockElement::ECAL,  
					  reco::PFBlock::LINKTEST_ALL);  
	     
	    float p_in=sqrt(elements[track->second].trackRef()->innerMomentum().x() * elements[track->second].trackRef()->innerMomentum().x() +  
			    elements[track->second].trackRef()->innerMomentum().y()*elements[track->second].trackRef()->innerMomentum().y()+  
			    elements[track->second].trackRef()->innerMomentum().z()*elements[track->second].trackRef()->innerMomentum().z());  
	    float linked_E=0;  
	    for(std::multimap<double, unsigned int>::iterator clust = moreClusters.begin();  
		clust != moreClusters.end(); ++clust)  
	      {  
		if(!active[clust->second])continue;  
		//running sum of linked energy  
		linked_E=linked_E+elements[clust->second].clusterRef()->energy();  
		//prevent too much energy from being added  
		if(linked_E/p_in>1.5)break;  
		bool included=false;  
		//check if these ecal clusters are already included with the supercluster  
		for(std::multimap<double, unsigned int>::iterator cluscheck = ecalAssoPFClusters.begin();  
		    cluscheck != ecalAssoPFClusters.end(); ++cluscheck)  
		  {  
		    if(cluscheck->second==clust->second)included=true;  
		  }  
		if(!included)AddClusters.push_back(clust->second);//Add to a container of clusters to be Added to the Photon candidate  
	      }  
	  }

	// Possibly need to be more smart about them (CASE 5)
	// .... for now we simply skip non id'ed tracks
	if( ! (trackRef->trackType(reco::PFBlockElement::T_FROM_GAMMACONV) ) ) continue;  
	hasConvTrack=true;  
	elemsToLock.push_back(track->second);
	//again look at the clusters linked to this track  
	//if(elements[track->second].convRef().isNonnull())
	//{	    
	//  ConversionsRef_.push_back(elements[track->second].convRef());
	//}
	std::multimap<double, unsigned int> moreClusters;  
	blockRef->associatedElements( track->second,  
				      linkData,  
				      moreClusters,  
				      reco::PFBlockElement::ECAL,  
				      reco::PFBlock::LINKTEST_ALL);
	
	float p_in=sqrt(elements[track->second].trackRef()->innerMomentum().x() * elements[track->second].trackRef()->innerMomentum().x() +  
			elements[track->second].trackRef()->innerMomentum().y()*elements[track->second].trackRef()->innerMomentum().y()+  
			elements[track->second].trackRef()->innerMomentum().z()*elements[track->second].trackRef()->innerMomentum().z());  
	float linked_E=0;  
	for(std::multimap<double, unsigned int>::iterator clust = moreClusters.begin();  
	    clust != moreClusters.end(); ++clust)  
	  {  
	    if(!active[clust->second])continue;  
	    linked_E=linked_E+elements[clust->second].clusterRef()->energy();  
	    if(linked_E/p_in>1.5)break;  
	    bool included=false;  
	    for(std::multimap<double, unsigned int>::iterator cluscheck = ecalAssoPFClusters.begin();  
		cluscheck != ecalAssoPFClusters.end(); ++cluscheck)  
	      {  
		if(cluscheck->second==clust->second)included=true;  
	      }  
	    if(!included)AddClusters.push_back(clust->second);//again only add if it is not already included with the supercluster  
	  }
	
	// we need to check for other TRACKS linked to this conversion track, that point possibly no an ECAL cluster not included in the SC
	// .... This is basically CASE 4.
	
	std::multimap<double, unsigned int> moreTracks;
	blockRef->associatedElements( track->second,
				      linkData,
				      moreTracks,
				      reco::PFBlockElement::TRACK,
				      reco::PFBlock::LINKTEST_ALL);
	
	for(std::multimap<double, unsigned int>::iterator track2 = moreTracks.begin();
	    track2 != moreTracks.end(); ++track2) {
	  
	  // first check if is it's still active
	  if( ! (active[track2->second]) ) continue;
	  //skip over the 1st leg already found above  
	  if(track->second==track2->second)continue;	  
	  // check if it's a CONV track
	  const reco::PFBlockElementTrack * track2Ref = dynamic_cast<const reco::PFBlockElementTrack*>((&elements[track2->second])); 	
	  if( ! (track2Ref->trackType(reco::PFBlockElement::T_FROM_GAMMACONV) ) ) continue;  // Possibly need to be more smart about them (CASE 5)
	  elemsToLock.push_back(track2->second);
	  // so it's another active conversion track, that is in the Block and linked to the conversion track we already found
	  // find the ECAL cluster linked to it...
	  std::multimap<double, unsigned int> convEcal;
	  blockRef->associatedElements( track2->second,
					linkData,
					convEcal,
					reco::PFBlockElement::ECAL,
					reco::PFBlock::LINKTEST_ALL);
	  float p_in=sqrt(elements[track->second].trackRef()->innerMomentum().x()*elements[track->second].trackRef()->innerMomentum().x()+
			  elements[track->second].trackRef()->innerMomentum().y()*elements[track->second].trackRef()->innerMomentum().y()+  
			  elements[track->second].trackRef()->innerMomentum().z()*elements[track->second].trackRef()->innerMomentum().z());  
	  
	  
	  float linked_E=0;
	  for(std::multimap<double, unsigned int>::iterator itConvEcal = convEcal.begin();
	      itConvEcal != convEcal.end(); ++itConvEcal) {
	    
	    if( ! (active[itConvEcal->second]) ) continue;
	    bool included=false;  
	    for(std::multimap<double, unsigned int>::iterator cluscheck = ecalAssoPFClusters.begin();  
		cluscheck != ecalAssoPFClusters.end(); ++cluscheck)  
	      {  
		if(cluscheck->second==itConvEcal->second)included=true;  
	      }
	    linked_E=linked_E+elements[itConvEcal->second].clusterRef()->energy();
	    if(linked_E/p_in>1.5)break;
	    if(!included){AddClusters.push_back(itConvEcal->second);
	    }
	    
	    // it's still active, so we have to add it.
	    // CAUTION: we don't care here if it's part of the SC or not, we include it anyways
	    
	    // loop over the HCAL Clusters linked to the ECAL cluster (CASE 6)
	    std::multimap<double, unsigned int> hcalElems_conv;
	    blockRef->associatedElements( itecal->second,linkData,
					  hcalElems_conv,
					  reco::PFBlockElement::HCAL,
					  reco::PFBlock::LINKTEST_ALL );
	    
	    for(std::multimap<double, unsigned int>::iterator ithcal2 = hcalElems_conv.begin();
		ithcal2 != hcalElems_conv.end(); ++ithcal2) {
	      
	      if ( ! (active[ithcal2->second] ) ) continue; // HCAL Cluster already used....
	      
	      // TODO: Decide if this HCAL cluster is to be used
	      // .... based on some Physics
	      // .... To we need to check if it's closer to any other ECAL/TRACK?
	      
	      bool useHcal = true;
	      if ( !useHcal ) continue;
	      
	      //elemsToLock.push_back(ithcal2->second);

	    } // end of loop over HCAL clusters linked to the ECAL cluster from second CONVERSION leg
	    
	  } // end of loop over ECALs linked to second T_FROM_GAMMACONV
	  
	} // end of loop over SECOND conversion leg

	// TODO: Do we need to check separatly if there are HCAL cluster linked to the track?
	
      } // end of loop over tracks
      
            
      // Calibrate the Added ECAL energy
      float addedCalibEne=0;
      float addedRawEne=0;
      std::vector<double>AddedPS1(0);
      std::vector<double>AddedPS2(0);  
      double addedps1=0;  
      double addedps2=0;  
      for(unsigned int i=0; i<AddClusters.size(); i++)  
	{  
	  std::multimap<double, unsigned int> PS1Elems_conv;  
	  std::multimap<double, unsigned int> PS2Elems_conv;  
	  blockRef->associatedElements(AddClusters[i],  
				       linkData,  
				       PS1Elems_conv,  
				       reco::PFBlockElement::PS1,  
				       reco::PFBlock::LINKTEST_ALL );  
	  blockRef->associatedElements( AddClusters[i],  
					linkData,  
					PS2Elems_conv,  
					reco::PFBlockElement::PS2,  
					reco::PFBlock::LINKTEST_ALL );  
	   
	  for(std::multimap<double, unsigned int>::iterator iteps = PS1Elems_conv.begin();  
	      iteps != PS1Elems_conv.end(); ++iteps)  
	    {  
	      if(!active[iteps->second])continue;  
	      std::multimap<double, unsigned int> PS1Elems_check;  
	      blockRef->associatedElements(iteps->second,  
					   linkData,  
					   PS1Elems_check,  
					   reco::PFBlockElement::ECAL,  
					   reco::PFBlock::LINKTEST_ALL );  
	      if(PS1Elems_check.begin()->second==AddClusters[i])  
		{  
		   
		  reco::PFClusterRef ps1ClusterRef = elements[iteps->second].clusterRef();  
		  AddedPS1.push_back(ps1ClusterRef->energy());  
		  addedps1=addedps1+ps1ClusterRef->energy();  
		  elemsToLock.push_back(iteps->second);  
		}  
	    }  
	   
	  for(std::multimap<double, unsigned int>::iterator iteps = PS2Elems_conv.begin();  
	      iteps != PS2Elems_conv.end(); ++iteps) {  
	    if(!active[iteps->second])continue;  
	    std::multimap<double, unsigned int> PS2Elems_check;  
	    blockRef->associatedElements(iteps->second,  
					 linkData,  
					 PS2Elems_check,  
					 reco::PFBlockElement::ECAL,  
					 reco::PFBlock::LINKTEST_ALL );  
	     
	    if(PS2Elems_check.begin()->second==AddClusters[i])  
	      {  
		reco::PFClusterRef ps2ClusterRef = elements[iteps->second].clusterRef();  
		AddedPS2.push_back(ps2ClusterRef->energy());  
		addedps2=addedps2+ps2ClusterRef->energy();  
		elemsToLock.push_back(iteps->second);  
	      }  
	  }  
	  reco::PFClusterRef AddclusterRef = elements[AddClusters[i]].clusterRef();  
	  addedRawEne = AddclusterRef->energy()+addedRawEne;  
	  addedCalibEne = thePFEnergyCalibration_->energyEm(*AddclusterRef,AddedPS1,AddedPS2,false)+addedCalibEne;  
	  AddedPS2.clear(); 
	  AddedPS1.clear();  
	  elemsToLock.push_back(AddClusters[i]);  
	}  
      AddClusters.clear();
      float EE=thePFEnergyCalibration_->energyEm(*clusterRef,ps1Ene,ps2Ene,false)+addedCalibEne; 
      if(useReg_){
	if(clusterRef->layer()==PFLayer::ECAL_BARREL){
	  float LocCorr=EvaluateLCorrMVA(clusterRef);
	  EE=LocCorr*clusterRef->energy()+addedCalibEne;
	  //cout<<"LocCorr "<<LocCorr<<endl;
	  //cout<<"Clust E "<<clusterRef->energy()<<"Locally Corrected "<<clusterRef->energy()*LocCorr<<endl;
	}
	else EE = thePFEnergyCalibration_->energyEm(*clusterRef,ps1Ene,ps2Ene,false)+addedCalibEne; 
      }
      //cout<<"Original Energy "<<EE<<"Added Energy "<<addedCalibEne<<endl;
      
      photonEnergy_ +=  EE;
      RawEcalEne    +=  clusterRef->energy()+addedRawEne;
      photonX_      +=  EE * clusterRef->position().X();
      photonY_      +=  EE * clusterRef->position().Y();
      photonZ_      +=  EE * clusterRef->position().Z();	        
      ps1TotEne     +=  ps1+addedps1;
      ps2TotEne     +=  ps2+addedps2;
    } // end of loop over all ECAL cluster within this SC
    AddFromElectron_.clear();
    float Elec_energy=0;
    float Elec_rawEcal=0;
    float Elec_totPs1=0;
    float Elec_totPs2=0;
    float ElectronX=0;
    float ElectronY=0;
    float ElectronZ=0;  
    std::vector<double>AddedPS1(0);
    std::vector<double>AddedPS2(0);
    
    EarlyConversion(    
   	    tempElectronCandidates,
    	    sc
    	    );   
    
    if(AddFromElectron_.size()>0)
      {	
	//collect elements from early Conversions that are reconstructed as Electrons
	
	for(std::vector<unsigned int>::const_iterator it = 
	      AddFromElectron_.begin();
	    it != AddFromElectron_.end(); ++it)
	  {
	    
	    if(elements[*it].type()== reco::PFBlockElement::ECAL)
	      {
	        //cout<<"Cluster ind "<<*it<<endl;
		AddedPS1.clear();
		AddedPS2.clear();
		unsigned int index=*it;
		reco::PFClusterRef clusterRef = 
		elements[index].clusterRef();
		//match to PS1 and PS2 to this cluster for calibration
		Elec_rawEcal=Elec_rawEcal+
		  elements[index].clusterRef()->energy();
		std::multimap<double, unsigned int> PS1Elems;  
		std::multimap<double, unsigned int> PS2Elems;  
		
		blockRef->associatedElements(index,  
					     linkData,  
					     PS1Elems,  					     reco::PFBlockElement::PS1,  
					     reco::PFBlock::LINKTEST_ALL );  
		blockRef->associatedElements( index,  
					      linkData,  
					      PS2Elems,  
					      reco::PFBlockElement::PS2,  
					      reco::PFBlock::LINKTEST_ALL );
		
		
		for(std::multimap<double, unsigned int>::iterator iteps = 
		      PS1Elems.begin();  
		    iteps != PS1Elems.end(); ++iteps) 
		  {  
		    std::multimap<double, unsigned int> Clustcheck;  		    	    blockRef->associatedElements( iteps->second,  								   linkData,  
															  Clustcheck,  
															  reco::PFBlockElement::ECAL,  
															  reco::PFBlock::LINKTEST_ALL );
		    if(Clustcheck.begin()->second==index)
		      {
			AddedPS1.push_back(elements[iteps->second].clusterRef()->energy());
			Elec_totPs1=Elec_totPs1+elements[iteps->second].clusterRef()->energy();
		      }
		  }
		
		for(std::multimap<double, unsigned int>::iterator iteps = 
		      PS2Elems.begin();  
		    iteps != PS2Elems.end(); ++iteps) 
		  {  
		    std::multimap<double, unsigned int> Clustcheck;  		    	    blockRef->associatedElements( iteps->second,  								   linkData,  
															  Clustcheck,  
															  reco::PFBlockElement::ECAL,  
															  reco::PFBlock::LINKTEST_ALL );
		    if(Clustcheck.begin()->second==index)
		      {
			AddedPS2.push_back(elements[iteps->second].clusterRef()->energy());
			Elec_totPs2=Elec_totPs2+elements[iteps->second].clusterRef()->energy();
		      }
		  }
		
	      //energy calibration 
		float EE=thePFEnergyCalibration_->
		  energyEm(*clusterRef,AddedPS1,AddedPS2,false);
		if(useReg_){
		  EE=0;
		  if(clusterRef->layer()==PFLayer::ECAL_BARREL){
		    float LocCorr=EvaluateLCorrMVA(clusterRef);
		    EE=LocCorr*clusterRef->energy();
		    //cout<<"LocCorr "<<LocCorr<<endl;
		    //cout<<"Clust E "<<clusterRef->energy()<<"Locally Corrected "<<clusterRef->energy()*LocCorr<<endl;		  
		  }
		  else EE=thePFEnergyCalibration_->energyEm(*clusterRef,AddedPS1,AddedPS2,false);
		}
		Elec_energy    += EE;
		ElectronX      +=  EE * clusterRef->position().X();
		ElectronY      +=  EE * clusterRef->position().Y();
		ElectronZ      +=  EE * clusterRef->position().Z();
		
	      }
	  }
	
      }
    
    //std::cout<<"Added Energy to Photon "<<Elec_energy<<" to "<<photonEnergy_<<std::endl;   
    photonEnergy_ +=  Elec_energy;
      RawEcalEne    +=  Elec_rawEcal;
      photonX_      +=  ElectronX;
      photonY_      +=  ElectronY;
      photonZ_      +=  ElectronZ;	        
      ps1TotEne     +=  Elec_totPs1;
      ps2TotEne     +=  Elec_totPs2;
    
    // we've looped over all ECAL clusters, ready to generate PhotonCandidate
    if( ! (photonEnergy_ > 0.) ) continue;    // This SC is not a Photon Candidate
    float sum_track_pt=0;
    //Now check if there are tracks failing isolation outside of the Jurassic isolation region  
    for(unsigned int i=0; i<IsoTracks.size(); i++)sum_track_pt=sum_track_pt+elements[IsoTracks[i]].trackRef()->pt();  
    


    math::XYZVector photonPosition(photonX_,
				   photonY_,
				   photonZ_);

    math::XYZVector photonDirection=photonPosition.Unit();
    
    math::XYZTLorentzVector photonMomentum(photonEnergy_* photonDirection.X(),
					   photonEnergy_* photonDirection.Y(),
					   photonEnergy_* photonDirection.Z(),
					   photonEnergy_           );

    if(sum_track_pt>(sumPtTrackIsoForPhoton_ + sumPtTrackIsoSlopeForPhoton_ * photonMomentum.pt()))
      {
	//cout<<"Hit Continue "<<endl;
	match_ind.clear(); //release the matched Electron candidates
	continue;
      }

	//THIS SC is not a Photon it fails track Isolation
    //if(sum_track_pt>(2+ 0.001* photonMomentum.pt()))
    //continue;//THIS SC is not a Photon it fails track Isolation

    /*
    std::cout<<" Created Photon with energy = "<<photonEnergy_<<std::endl;
    std::cout<<"                         pT = "<<photonMomentum.pt()<<std::endl;
    std::cout<<"                     RawEne = "<<RawEcalEne<<std::endl;
    std::cout<<"                          E = "<<photonMomentum.e()<<std::endl;
    std::cout<<"                        eta = "<<photonMomentum.eta()<<std::endl;
    std::cout<<"             TrackIsolation = "<< sum_track_pt <<std::endl;
    */

    reco::PFCandidate photonCand(0,photonMomentum, reco::PFCandidate::gamma);
    
    photonCand.setPs1Energy(ps1TotEne);
    photonCand.setPs2Energy(ps2TotEne);
    photonCand.setEcalEnergy(RawEcalEne,photonEnergy_);
    photonCand.setHcalEnergy(0.,0.);
    photonCand.set_mva_nothing_gamma(1.);  
    photonCand.setSuperClusterRef(sc->superClusterRef());
    math::XYZPoint v(primaryVertex_.x(), primaryVertex_.y(), primaryVertex_.z());
    photonCand.setVertex( v  );
    if(hasConvTrack || hasSingleleg)photonCand.setFlag( reco::PFCandidate::GAMMA_TO_GAMMACONV, true);
    int matches=match_ind.size();
    int count=0;
    for ( std::vector<reco::PFCandidate>::const_iterator ec=tempElectronCandidates.begin();   ec != tempElectronCandidates.end(); ++ec ){
      for(int i=0; i<matches; i++)
	{
	  if(count==match_ind[i])photonCand.addDaughter(*ec);
	  count++;
	}
    }
    //photonCand.setPositionAtECALEntrance(math::XYZPointF(photonMom_.position()));
    // set isvalid_ to TRUE since we've found at least one photon candidate
    isvalid_ = true;
    // push back the candidate into the collection ...
    //Add Elements from Electron
	for(std::vector<unsigned int>::const_iterator it = 
	      AddFromElectron_.begin();
	    it != AddFromElectron_.end(); ++it)photonCand.addElementInBlock(blockRef,*it);


    // ... and lock all elemts used
    for(std::vector<unsigned int>::const_iterator it = elemsToLock.begin();
	it != elemsToLock.end(); ++it)
      {
	if(active[*it])
	  {
	    photonCand.addElementInBlock(blockRef,*it);
	    if( elements[*it].type() == reco::PFBlockElement::TRACK  )
	      {
		if(elements[*it].convRef().isNonnull())
		  {
		    //make sure it is not stored already as the partner track
		    bool matched=false;
		    for(unsigned int ic = 0; ic < ConversionsRef_.size(); ic++)
		      {
			if(ConversionsRef_[ic]==elements[*it].convRef())matched=true;
		      }
		    if(!matched)ConversionsRef_.push_back(elements[*it].convRef());
		  }
	      }
	  }
	active[*it] = false;	
      }
    //Do Global Corrections here:
    if(useReg_){
    float GCorr=EvaluateGCorrMVA(photonCand);
    cout<<"GCorr "<<GCorr<<endl;
    math::XYZTLorentzVector photonCorrMomentum(GCorr*photonEnergy_* photonDirection.X(),
					       GCorr*photonEnergy_* photonDirection.Y(),
					       GCorr*photonEnergy_* photonDirection.Z(),
					       GCorr * photonEnergy_           );
    photonCand.setP4(photonCorrMomentum);
    }
    // here add the extra information
    PFCandidatePhotonExtra myExtra(sc->superClusterRef());

    //    Daniele example for mvaValues
    //    do the same for single leg trackRef and convRef
    for(unsigned int ic = 0; ic < MVA_values.size(); ic++)
      {
	myExtra.addSingleLegConvMva(MVA_values[ic]);
	myExtra.addSingleLegConvTrackRef(singleLegRef[ic]);
	//cout<<"Single Leg Tracks "<<singleLegRef[ic]->pt()<<" MVA "<<MVA_values[ic]<<endl;
      }
    for(unsigned int ic = 0; ic < ConversionsRef_.size(); ic++)
      {
	myExtra.addConversionRef(ConversionsRef_[ic]);
	//cout<<"Conversion Pairs "<<ConversionsRef_[ic]->pairMomentum()<<endl;
      }
    pfPhotonExtraCandidates.push_back(myExtra);
    pfCandidates->push_back(photonCand);
    // ... and reset the vector
    elemsToLock.resize(0);
    hasConvTrack=false;
    hasSingleleg=false;
  } // end of loops over all elements in block
  
  return;

}

std::vector<int>PFPhotonAlgo::getPFMustacheClus(int nClust, std::vector<float>& ClustEt, std::vector<float>& ClustEta, std::vector<float>& ClustPhi){
  float etmax = 0;
  int imax = -1;
  float phot_eta_maxcl = 0.0;
  float phot_phi_maxcl = 0.0;
  std::vector<int> included(0);
  Float_t fExcluded = 0.0;
 
  float deta, dphi;
  float upper_cut, lower_cut;
  float b_upper, b_lower;
  float a_upper, a_lower;
  float curv_low, curv_up;
  float midpoint;
 
 //loop over clusters
 
  for(int k=0; k<nClust; k++){
    
    //search for highest Et cluster, set phi and eta
    
    if(etmax < ClustEt[k]){
      imax = k;
      etmax = ClustEt[k];
      phot_eta_maxcl = ClustEta[k];
      phot_phi_maxcl= ClustPhi[k];
    }//end search for highest Et cluster
    //cout<<"Here"<<endl;
  }//end loop over clusters
  
  for(int k=0; k<nClust; k++){
    
    deta = 0.0;
    dphi = 0.0; 
    upper_cut = 0.0;
    lower_cut = 0.0;
    b_upper = 0.0;
    b_lower = 0.0;
    a_upper = 0.0;
    a_lower = 0.0;
    curv_low = 0.0;
    curv_up = 0.0;
    midpoint = 0.0;  
    float w00 = -0.00571429;
    float w01 = -0.002;
    float w10 = 0.0135714;
    float w11 = 0.001;
    float p00 = -0.107537;
    float p01 = 0.590969;
    float p02 = -0.076494;
    float p10 = -0.0268843;
    float p11 = 0.147742;
    float p12 = -0.0191235;
    
    
    deta = sin(phot_eta_maxcl)*(ClustEta[k]-phot_eta_maxcl);	
    dphi = ClustPhi[k]-phot_phi_maxcl;
   
    //2 parabolas (upper and lower) 
    //of the form: y = a*x*x + b      
    
    //b comes from a fit to the width
    //and has a slight dependence on Et on the upper edge
    
    b_lower = w00*sin(phot_eta_maxcl)*phot_eta_maxcl + w01 / sqrt(log10(ClustEt[k])+1.1);
    b_upper = w10*sin(phot_eta_maxcl)*phot_eta_maxcl + w11  / sqrt(log10(ClustEt[k])+1.1);
    //cout<<"upper_b "<< b_upper<<" lower b "<<b_lower <<endl;
    //here make an adjustment to the width for the offset from 0.
    midpoint = b_upper - (b_upper-b_lower)/2.;
    b_lower = b_lower - midpoint;
    b_upper = b_upper - midpoint;
    
    //the curvature comes from a parabolic 
    //fit for many slices in eta given a 
    //slice -0.1 < log10(Et) < 0.1
    curv_up = p00*pow(phot_eta_maxcl*sin(phot_eta_maxcl),2)+p01*phot_eta_maxcl*sin(phot_eta_maxcl)+p02;
    curv_low = p10*pow(phot_eta_maxcl*sin(phot_eta_maxcl),2)+p11*phot_eta_maxcl*sin(phot_eta_maxcl)+p12;
    
    //solving for the curviness given the width of this particular point
    a_lower = (1/(4*curv_low))-fabs(b_lower);
    a_upper = (1/(4*curv_up))-fabs(b_upper);
    //cout<<"upper_a "<< a_upper<<" lower a "<<a_lower <<endl;
    upper_cut =(1./(4.*a_upper))*pow(dphi,2)+b_upper;
    lower_cut =(1./(4.*a_lower))*pow(dphi,2)+b_lower;
    //cout<<"upper_cut "<< upper_cut<<" lower cut "<<lower_cut<<endl;
    //remove the Not in the function to make it included clusters
    if ((deta < upper_cut && deta > lower_cut)){
      included.push_back(k);
    }
    
    
  }
  
 return included;
 
}

float PFPhotonAlgo::EvaluateGCorrMVA(reco::PFCandidate photon){
  float BDTG=1;
  PFPhoEta_=photon.eta();
  PFPhoPhi_=photon.phi();
  PFPhoEt_=photon.pt();
  //recalculate R9 from sum PFClusterEnergy and E3x3 from Highest PFCluster Energy
  SCPhiWidth_=photon.superClusterRef()->phiWidth();
  SCEtaWidth_=photon.superClusterRef()->etaWidth();  
  //get from track with min R;
  RConv_=130;
  float ClustSumEt=0;
  std::vector<float>Clust_E(0);
  std::vector<float>Clust_Et(0);
  std::vector<float>Clust_Eta(0);
  std::vector<float>Clust_Phi(0);
  //Multimap to sort clusters by energy
  std::multimap<float, int>Clust;
  PFCandidate::ElementsInBlocks eleInBlocks = photon.elementsInBlocks();
  for(unsigned i=0; i<eleInBlocks.size(); i++)
    {
      PFBlockRef blockRef = eleInBlocks[i].first;
      unsigned indexInBlock = eleInBlocks[i].second;
      const edm::OwnVector< reco::PFBlockElement >&  elements=eleInBlocks[i].first->elements();
      const reco::PFBlockElement& element = elements[indexInBlock];
      if(element.type()==reco::PFBlockElement::TRACK){
	float R=sqrt(element.trackRef()->innerPosition().X()*element.trackRef()->innerPosition().X()+element.trackRef()->innerPosition().Y()*element.trackRef()->innerPosition().Y());
	if(RConv_>R)RConv_=R;
      }
      
      if(element.type()==reco::PFBlockElement::ECAL){
	reco::PFClusterRef ClusterRef = element.clusterRef();
	Clust_E.push_back(ClusterRef->energy());
	Clust_Et.push_back(ClusterRef->pt());	
	ClustSumEt=ClustSumEt+ClusterRef->pt();
	Clust_Eta.push_back(ClusterRef->eta());
	Clust_Phi.push_back(ClusterRef->phi());
      }
      
      if(element.type()==reco::PFBlockElement::GSF)
	{
	  //elements[indexInBlock].GsftrackRef();
	  // RConv_=sqrt(element.GsftrackRef()->innerPosition().X()*element.GsftrackRef()->innerPosition().X() + element.GsftrackRef()->innerPosition().Y()*element.GsftrackRef()->innerPosition().Y());
	}
      
    }
  
  nPFClus_=Clust_Et.size();
  std::vector<int>included(0);
  included=getPFMustacheClus(nPFClus_, Clust_Et, Clust_Eta, Clust_Phi);
  excluded_=nPFClus_-included.size();
  //order the clusters by energy
  float Mustache_Et=0;
  float ClusSum=0;
  for(unsigned int i=0; i<included.size(); ++i)
    {
      
      Clust.insert(make_pair(Clust_E[i], i));
      Mustache_Et=Mustache_Et+Clust_Et[i];
      ClusSum=ClusSum+Clust_E[i];
    }
  std::multimap<float, int>::reverse_iterator it;
  it=Clust.rbegin();
  int max_c=(*it).second;
  it=Clust.rend();
  int min_c=(*it).second;
  if(nPFClus_>1)LowClusE_=Clust_E[min_c];
  else LowClusE_=0;
  if(nPFClus_>1){
    dEta_=fabs(Clust_Eta[max_c]-Clust_Eta[min_c]);
    dPhi_=acos(cos(Clust_Phi[max_c]-Clust_Phi[min_c]));
  }
  else{
    dEta_=0;
    dPhi_=0;
  }
  Mustache_EtRatio_=(Mustache_Et-ClustSumEt)/ClustSumEt;
  
  float dRmin=999;
  float SCphi=photon.superClusterRef()->position().phi();
  float SCeta=photon.superClusterRef()->position().eta();
  for(unsigned i=0; i<eleInBlocks.size(); i++)
    {
      PFBlockRef blockRef = eleInBlocks[i].first;
      unsigned indexInBlock = eleInBlocks[i].second;
      const edm::OwnVector< reco::PFBlockElement >&  elements=eleInBlocks[i].first->elements();
      const reco::PFBlockElement& element = elements[indexInBlock];
      if(element.type()==reco::PFBlockElement::ECAL){
	reco::PFClusterRef ClusterRef = element.clusterRef();
	float eta=ClusterRef->position().eta();
	float phi=ClusterRef->position().phi();
	float dR=deltaR(SCeta, SCphi, eta, phi);
	if(dR<dRmin){
	  dRmin=dR;
	  fill5x5Map(ClusterRef);
	  PFPhoR9_=e3x3_/ClusSum;
	}
      }
    } 
  //fill Material Map:
  int ix = X0_sum->GetXaxis()->FindBin(PFPhoEta_);
  int iy = X0_sum->GetYaxis()->FindBin(PFPhoPhi_);
  x0inner_= X0_inner->GetBinContent(ix,iy);
  x0middle_=X0_middle->GetBinContent(ix,iy);
  x0outer_=X0_outer->GetBinContent(ix,iy);
  BDTG=tmvaGCRegReader_->EvaluateRegression("BDTG_GCorr")[0];
  //  cout<<"BDTG Parameters X0"<<x0inner_<<", "<<x0middle_<<", "<<x0outer_<<endl;
  // cout<<"Et, Eta, Phi "<<PFPhoEt_<<", "<<PFPhoEta_<<", "<<PFPhoPhi_<<endl;
  // cout<<"PFPhoR9 "<<PFPhoR9_<<endl;
  // cout<<"R "<<RConv_<<endl;
  
  return BDTG;

}

float PFPhotonAlgo::EvaluateLCorrMVA(reco::PFClusterRef clusterRef ){
  float BDTG=1;
  
  GetCrysCoordinates(clusterRef);
  fill5x5Map(clusterRef);
  VtxZ_=primaryVertex_.z();
  ClusPhi_=clusterRef->position().phi(); 
  ClusEta_=fabs(clusterRef->position().eta());
  EB=fabs(clusterRef->position().eta())/clusterRef->position().eta();
  logPFClusE_=log(clusterRef->energy());
  /*
  cout<<"BDTG Parameters "<<" Crys Eta, Phi "<<CrysEta_<<", "<<CrysPhi_<<endl;
  cout<<"BDTG Parameters "<<" Crys Eta, Phi Index "<<CrysIEta_<<", "<<CrysIPhi_<<endl;
  cout<<"BDTG Parameters "<<" Clus Eta, Phi "<<ClusEta_<<", "<<ClusPhi_<<endl;
  cout<<"BDTG Parameters "<<" EB "<<EB<<endl;
  cout<<"BDTG Parameters "<<" Z "<<VtxZ_<<endl;
  cout<<"BDTG Parameters "<<" log E "<<logPFClusE_<<endl;
  cout<<"BDTG Parameters "<<" R9 & 5x5 "<<ClusR9_<<", "<<Clus5x5ratio_<<endl;
  */
  BDTG=tmvaLCRegReader_->EvaluateRegression("BDTG_LocCorr")[0];
  return BDTG;
  
}


void PFPhotonAlgo::GetCrysCoordinates(reco::PFClusterRef clusterRef){
  float PFSeedEta=99;
  float PFSeedPhi=99;
  float PFSeedTheta=99;
  double PFSeedE=0;
  unsigned int SeedDetId=-1;
  float seedPhi=0;
  float seedEta=0;
  DetId idseed;
  const std::vector< reco::PFRecHitFraction >& PFRecHits=
    clusterRef->recHitFractions();
  for ( std::vector< reco::PFRecHitFraction >::const_iterator it = PFRecHits.begin();
	it != PFRecHits.end(); ++it){
    const PFRecHitRef& RefPFRecHit = it->recHitRef();
    unsigned index=it-PFRecHits.begin();
    float frac=clusterRef->hitsAndFractions()[index].second;
    float E= RefPFRecHit->energy()* frac;
    if(E>PFSeedE){
      SeedDetId=RefPFRecHit.index();
      PFSeedE=E;  
      PFSeedEta=RefPFRecHit->positionREP().eta(); 
      PFSeedPhi=RefPFRecHit->positionREP().phi();
      PFSeedTheta=RefPFRecHit->positionREP().theta();
      RefPFRecHit->positionREP().theta();
      idseed = RefPFRecHit->detId();
    }
  }
  EBDetId EBidSeed=EBDetId(idseed.rawId());
  CrysIEta_=EBidSeed.ieta();
  CrysIPhi_=EBidSeed.iphi();
  
  //Crystal Coordinates:
  double Pi=3.14159265358979323846;
  float Phi=clusterRef->position().phi(); 
  float Eta=clusterRef->position().eta();
  double Theta = -(clusterRef->position().theta())+0.5* Pi;
  double PhiCentr = TVector2::Phi_mpi_pi(PFSeedPhi);
  double PhiWidth = (Pi/180.);
  double PhiCry = (TVector2::Phi_mpi_pi(Phi-PhiCentr))/PhiWidth;
  double ThetaCentr = -PFSeedTheta+0.5*Pi;
  double ThetaWidth = (Pi/180.)*cos(ThetaCentr);
  
  //cout<<"Clust Theta "<<Theta<<" Crys Theta "<<ThetaCentr<<endl;
  //cout<<" Width "<< ThetaWidth<<endl;
  double EtaCry = (Theta-ThetaCentr)/ThetaWidth; 
  CrysEta_=EtaCry;
  CrysPhi_=PhiCry;
 
  //check Module and crack:
  int iphi=CrysIPhi_;
  int phimod=iphi%20;
  if(phimod>1)PFCrysPhiCrack_=2;
  else PFCrysPhiCrack_=phimod; //should be 0, 1
  
  if(abs(CrysIEta_)==1 || abs(CrysIEta_)==2 )
    PFCrysEtaCrack_=abs(CrysIEta_);
  if(abs(CrysIEta_)>2 && abs(CrysIEta_)<24)
    PFCrysEtaCrack_=3;
  if(abs(CrysIEta_)==24)
    PFCrysEtaCrack_=4;
  if(abs(CrysIEta_)==25)
	PFCrysEtaCrack_=5;
      if(abs(CrysIEta_)==26)
	PFCrysEtaCrack_=6;
      if(abs(CrysIEta_)==27)
		PFCrysEtaCrack_=7;
      if(abs(CrysIEta_)>27 &&  abs(CrysIEta_)<44)
	PFCrysEtaCrack_=8;
      if(abs(CrysIEta_)==44)
		PFCrysEtaCrack_=9;
      if(abs(CrysIEta_)==45)
	PFCrysEtaCrack_=10;
      if(abs(CrysIEta_)==46)
		PFCrysEtaCrack_=11;
      if(abs(CrysIEta_)==47)
	PFCrysEtaCrack_=12;
      if(abs(CrysIEta_)>47 &&  abs(CrysIEta_)<64)
	PFCrysEtaCrack_=13;
      if(abs(CrysIEta_)==64)
	PFCrysEtaCrack_=14;
      if(abs(CrysIEta_)==65)
	PFCrysEtaCrack_=15;
      if(abs(CrysIEta_)==66)
	PFCrysEtaCrack_=16;
      if(abs(CrysIEta_)==67)
	PFCrysEtaCrack_=17;
      if(abs(CrysIEta_)>67 &&  abs(CrysIEta_)<84)
	PFCrysEtaCrack_=18;
      if(abs(CrysIEta_)==84)
	PFCrysEtaCrack_=19;
      if(abs(CrysIEta_)==85)
	PFCrysEtaCrack_=20;
      
}

void PFPhotonAlgo::fill5x5Map(reco::PFClusterRef clusterRef){

  float PFSeedEta=99;
  float PFSeedPhi=99;
  double PFSeedE=0;
  unsigned int SeedDetId=-1;
  float seedPhi=0;
  float seedEta=0;
  DetId idseed;
  const std::vector< reco::PFRecHitFraction >& PFRecHits=
    clusterRef->recHitFractions();
  for ( std::vector< reco::PFRecHitFraction >::const_iterator it = PFRecHits.begin();
	it != PFRecHits.end(); ++it){
    const PFRecHitRef& RefPFRecHit = it->recHitRef();
    unsigned index=it-PFRecHits.begin();
    //DetId id=eclusterRef->hitsAndFractions()[index].first;
    float frac=clusterRef->hitsAndFractions()[index].second;
    float E= RefPFRecHit->energy()* frac;
    if(E>PFSeedE){
      SeedDetId=RefPFRecHit.index();
      PFSeedE=E;  
      PFSeedEta=RefPFRecHit->positionREP().eta(); 
      PFSeedPhi=RefPFRecHit->positionREP().phi(); 
      idseed = RefPFRecHit->detId();
    }
  }
  
  
  //initialize 5x5 map
  for(int i=0; i<5; ++i)
    for(int j=0; j<5; ++j)e5x5Map[i][j]=0;
  float E3x3=0;
  float E5x5=0;
  int count=0;
  for ( std::vector< reco::PFRecHitFraction >::const_iterator it = PFRecHits.begin();
	it != PFRecHits.end(); ++it){
    unsigned index=it-PFRecHits.begin();
    const PFRecHitRef& RefPFRecHit = it->recHitRef();
    float frac=clusterRef->hitsAndFractions()[index].second;
    DetId id = RefPFRecHit->detId();
    if(idseed.subdetId()==EcalBarrel){
      int deta=EBDetId::distanceEta(id,idseed);
      int dphi=EBDetId::distancePhi(id,idseed);	
      EBDetId EBidSeed=EBDetId(idseed.rawId());
      if(abs(dphi)<=1 && abs(deta)<=1)
	{
	  E3x3=E3x3+(RefPFRecHit->energy()*frac);
	}
      if(abs(dphi)<=2 && abs(deta)<=2)
	{
	  //center the array on [2][2] to be the Seed
	  //deta and deta are unsigned so you want to recompute 
	  //to get top bottom left right 
	  EBDetId EBid=EBDetId(id.rawId());
	  //note this is actually opposite to make it consistent to the lazy tools left right which inverts them
	  int i=EBidSeed.ieta()-EBid.ieta();
	  
	  int j=EBid.iphi()-EBidSeed.iphi();
	  int iEta=i+2;
	  int iPhi=j+2;
	  
	  e5x5Map[iEta][iPhi]=RefPFRecHit->energy()*frac;
	  E5x5=E5x5+(RefPFRecHit->energy()*frac);
	}
    }
    if(idseed.subdetId()==EcalEndcap){
      //dx and dy are unsigned so you want to recompute 
      //to get top bottom left right 
      int dx=EEDetId::distanceX(id,idseed);
      int dy=EEDetId::distanceY(id,idseed);
      EEDetId EEidSeed=EEDetId(idseed.rawId());
      //dx and dy are unsigned so you want to recompute 
      //to get top bottom left right 
      if(abs(dx)<=1 && abs(dy)<=1)
	{
	  E3x3=E3x3+(RefPFRecHit->energy()*frac);
	}
      if(abs(dx)<=2 && abs(dy)<=2)
	{
	  EEDetId EEid=EEDetId(id.rawId());
	  int i=EEid.ix()-EEidSeed.ix();
	  int j=EEid.iy()-EEidSeed.iy();
	  //center the array on [2][2] to be the Seed
	  int ix=i+2;
	  int iy=j+2;
	  e5x5Map[ix][iy]=RefPFRecHit->energy()*frac;
	  E5x5=E5x5+(RefPFRecHit->energy()*frac);
	}
    }
    
  }
  e3x3_=E3x3;
  ClusR9_=E3x3/clusterRef->energy();
  Clus5x5ratio_=E5x5/clusterRef->energy();
  //cout<<"E5x5 "<<E5x5<<" Clus Energy "<< clusterRef->energy();
}




bool PFPhotonAlgo::EvaluateSingleLegMVA(const reco::PFBlockRef& blockref, const reco::Vertex& primaryvtx, unsigned int track_index)  
{  
  bool convtkfound=false;  
  const reco::PFBlock& block = *blockref;  
  const edm::OwnVector< reco::PFBlockElement >& elements = block.elements();  
  //use this to store linkdata in the associatedElements function below  
  PFBlock::LinkData linkData =  block.linkData();  
  //calculate MVA Variables  
  chi2=elements[track_index].trackRef()->chi2()/elements[track_index].trackRef()->ndof();  
  nlost=elements[track_index].trackRef()->trackerExpectedHitsInner().numberOfLostHits();  
  nlayers=elements[track_index].trackRef()->hitPattern().trackerLayersWithMeasurement();  
  track_pt=elements[track_index].trackRef()->pt();  
  STIP=elements[track_index].trackRefPF()->STIP();  
   
  float linked_e=0;  
  float linked_h=0;  
  std::multimap<double, unsigned int> ecalAssoTrack;  
  block.associatedElements( track_index,linkData,  
			    ecalAssoTrack,  
			    reco::PFBlockElement::ECAL,  
			    reco::PFBlock::LINKTEST_ALL );  
  std::multimap<double, unsigned int> hcalAssoTrack;  
  block.associatedElements( track_index,linkData,  
			    hcalAssoTrack,  
			    reco::PFBlockElement::HCAL,  
			    reco::PFBlock::LINKTEST_ALL );  
  if(ecalAssoTrack.size() > 0) {  
    for(std::multimap<double, unsigned int>::iterator itecal = ecalAssoTrack.begin();  
	itecal != ecalAssoTrack.end(); ++itecal) {  
      linked_e=linked_e+elements[itecal->second].clusterRef()->energy();  
    }  
  }  
  if(hcalAssoTrack.size() > 0) {  
    for(std::multimap<double, unsigned int>::iterator ithcal = hcalAssoTrack.begin();  
	ithcal != hcalAssoTrack.end(); ++ithcal) {  
      linked_h=linked_h+elements[ithcal->second].clusterRef()->energy();  
    }  
  }  
  EoverPt=linked_e/elements[track_index].trackRef()->pt();  
  HoverPt=linked_h/elements[track_index].trackRef()->pt();  
  GlobalVector rvtx(elements[track_index].trackRef()->innerPosition().X()-primaryvtx.x(),  
		    elements[track_index].trackRef()->innerPosition().Y()-primaryvtx.y(),  
		    elements[track_index].trackRef()->innerPosition().Z()-primaryvtx.z());  
  double vtx_phi=rvtx.phi();  
  //delta Phi between conversion vertex and track  
  del_phi=fabs(deltaPhi(vtx_phi, elements[track_index].trackRef()->innerMomentum().Phi()));  
  mvaValue = tmvaReader_->EvaluateMVA("BDT");  
  if(mvaValue > MVACUT)convtkfound=true;  
  return convtkfound;  
}

//Recover Early Conversions reconstructed as PFelectrons
void PFPhotonAlgo::EarlyConversion(    
				   //std::auto_ptr< reco::PFCandidateCollection > 
				   //&pfElectronCandidates_,
				   std::vector<reco::PFCandidate>& 
				   tempElectronCandidates,
				   const reco::PFBlockElementSuperCluster* sc
				   ){
  //step 1 check temp electrons for clusters that match Photon Supercluster:
  // permElectronCandidates->clear();
  int count=0;
  for ( std::vector<reco::PFCandidate>::const_iterator ec=tempElectronCandidates.begin();   ec != tempElectronCandidates.end(); ++ec ) 
    {
      //      bool matched=false;
      int mh=ec->gsfTrackRef()->trackerExpectedHitsInner().numberOfLostHits();
      if(mh==0)continue;//Case where missing hits greater than zero
      
      reco::GsfTrackRef gsf=ec->gsfTrackRef();
      //some hoopla to get Electron SC ref
      
      if(gsf->extra().isAvailable() && gsf->extra()->seedRef().isAvailable()) 
	{
	  reco::ElectronSeedRef seedRef=  gsf->extra()->seedRef().castTo<reco::ElectronSeedRef>();
	  if(seedRef.isAvailable() && seedRef->isEcalDriven()) 
	    {
	      reco::SuperClusterRef ElecscRef = seedRef->caloCluster().castTo<reco::SuperClusterRef>();
	      
	      if(ElecscRef.isNonnull()){
		//finally see if it matches:
		reco::SuperClusterRef PhotscRef=sc->superClusterRef();
		
		if(PhotscRef==ElecscRef)
		  {
		    match_ind.push_back(count);
		    //  matched=true; 
		    //cout<<"Matched Electron with Index "<<count<<" This is the electron "<<*ec<<endl;
		    //find that they have the same SC footprint start to collect Clusters and tracks and these will be passed to PFPhoton
		    reco::PFCandidate::ElementsInBlocks eleInBlocks = ec->elementsInBlocks();
		    for(unsigned i=0; i<eleInBlocks.size(); i++) 
		      {
			reco::PFBlockRef blockRef = eleInBlocks[i].first;
			unsigned indexInBlock = eleInBlocks[i].second;	 
			//const edm::OwnVector< reco::PFBlockElement >&  elements=eleInBlocks[i].first->elements();
			//const reco::PFBlockElement& element = elements[indexInBlock];  		
			
			AddFromElectron_.push_back(indexInBlock);	       	
		      }		    
		  }		
	      }
	    }	  
	}           
      count++;
    }
}
