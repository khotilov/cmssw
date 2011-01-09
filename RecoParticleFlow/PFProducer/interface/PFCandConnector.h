#ifndef PFProducer_PFCandConnector_H_
#define PFProducer_PFCandConnector_H_

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
// \author : M. Gouzevitch
// \date : May 2010

/// Based on a class from : V. Roberfroid, February 2008

class PFCandConnector {
    
    public :
       
       PFCandConnector( ) { 
	 pfC_ = std::auto_ptr<reco::PFCandidateCollection>(new reco::PFCandidateCollection); 
	 debug_ = false;
	 bCorrect_ = false;
	 bCalibPrimary_ =  false;

	 fConst_.push_back(1), fConst_.push_back(0);
	 fNorm_.push_back(0), fNorm_.push_back(0);
	 fExp_.push_back(0);

	 dptRel_PrimaryTrack_ = 0.;
         dptRel_MergedTrack_ = 0.;
	 ptErrorSecondary_ = 0.;
       }
       
       void setParameters(const edm::ParameterSet& iCfgCandConnector){
	 
	 bool bCorrect, bCalibPrimary;
	 double dptRel_PrimaryTrack, dptRel_MergedTrack, ptErrorSecondary;
	 std::vector<double> nuclCalibFactors;

	 /// Flag to apply the correction procedure for nuclear interactions
	 bCorrect = iCfgCandConnector.getParameter<bool>("bCorrect");
	 /// Flag to calibrate the reconstructed nuclear interactions with primary or merged tracks
	 bCalibPrimary =  iCfgCandConnector.getParameter<bool>("bCalibPrimary");

	 if(iCfgCandConnector.exists("dptRel_PrimaryTrack")) dptRel_PrimaryTrack = iCfgCandConnector.getParameter<double>("dptRel_PrimaryTrack");
	 else { std::cout << "dptRel_PrimaryTrack doesn't exist. Setting a default safe value 0" << std::endl; dptRel_PrimaryTrack = 0;}

	 if(iCfgCandConnector.exists("dptRel_MergedTrack"))  dptRel_MergedTrack = iCfgCandConnector.getParameter<double>("dptRel_MergedTrack");
	 else { std::cout << "dptRel_MergedTrack doesn't exist. Setting a default safe value 0" << std::endl; dptRel_MergedTrack = 0;}

	 if(iCfgCandConnector.exists("ptErrorSecondary"))    ptErrorSecondary = iCfgCandConnector.getParameter<double>("ptErrorSecondary");
	 else { std::cout << "ptErrorSecondary doesn't exist. Setting a default safe value 0" << std::endl; ptErrorSecondary = 0;}

	 if(iCfgCandConnector.exists("nuclCalibFactors"))    nuclCalibFactors = iCfgCandConnector.getParameter<std::vector<double> >("nuclCalibFactors");  
	 else { std::cout << "nuclear calib factors doesn't exist the factor would not be applyed" << std::endl; }

	 setParameters(bCorrect, bCalibPrimary, dptRel_PrimaryTrack, dptRel_MergedTrack, ptErrorSecondary, nuclCalibFactors);

       }


       void setParameters(bool bCorrect, bool bCalibPrimary, double dptRel_PrimaryTrack, double dptRel_MergedTrack, double ptErrorSecondary, std::vector<double> nuclCalibFactors){

	 bCorrect_ = bCorrect;
	 bCalibPrimary_ =  bCalibPrimary;
	 dptRel_PrimaryTrack_ = dptRel_PrimaryTrack;
	 dptRel_MergedTrack_ = 	dptRel_MergedTrack;
	 ptErrorSecondary_ = ptErrorSecondary;

	 if (nuclCalibFactors.size() == 5) {
	   fConst_[0] =  nuclCalibFactors[0]; 
	   fConst_[1] =  nuclCalibFactors[1]; 
	   
	   fNorm_[0] = nuclCalibFactors[2]; 
	   fNorm_[1] = nuclCalibFactors[3]; 
	   
	   fExp_[0] = nuclCalibFactors[4];
	 } else {
	   std::cout << "Wrong calibration factors for nuclear interactions. The calibration procedure would not be applyed." << std::endl;
	   bCalibPrimary_ =  false;
	 }

	 std::string sCorrect = bCorrect_ ? "On" : "Off";
	 std::cout << " ====================== The PFCandConnector is switched " << sCorrect.c_str() << " ==================== " << std::endl;
	 std::string sCalibPrimary = bCalibPrimary_ ? "used for calibration" : "not used for calibration";
	 if (bCorrect_) std::cout << "Primary Tracks are " << sCalibPrimary.c_str() << std::endl;
	 if (bCorrect_ && bCalibPrimary_) std::cout << "Under the condition that the precision on the Primary track is better than " << dptRel_PrimaryTrack_ << " % "<< std::endl;
	 if (bCorrect_ && bCalibPrimary_) std::cout << "      and on merged tracks better than " <<  dptRel_MergedTrack_ << " % "<< std::endl;
	 if (bCorrect_ && bCalibPrimary_) std::cout << "      and secondary tracks in some cases more precise than " << ptErrorSecondary_ << " GeV"<< std::endl;
	 if (bCorrect_ && bCalibPrimary_) std::cout << "factor = (" << fConst_[0] << " + " << fConst_[1] << "*cFrac) - (" 
				       << fNorm_[0] << " - " << fNorm_[1] << "cFrac)*exp( "
				       << -1*fExp_[0] << "*pT )"<< std::endl;
	 std::cout << " =========================================================== " << std::endl;	 

       }

       void setDebug( bool debug ) {debug_ = debug;}

       

       std::auto_ptr<reco::PFCandidateCollection> connect(std::auto_ptr<reco::PFCandidateCollection>& pfCand);

       
 
    private :

       /// Analyse nuclear interactions where a primary or merged track is present
       void analyseNuclearWPrim(std::auto_ptr<reco::PFCandidateCollection>&, unsigned int);

       /// Analyse nuclear interactions where a secondary track is present
       void analyseNuclearWSec(std::auto_ptr<reco::PFCandidateCollection>&, unsigned int);

       bool isPrimaryNucl( const reco::PFCandidate& pf ) const;

       bool isSecondaryNucl( const reco::PFCandidate& pf ) const;

       /// Return a calibration factor for a reconstructed nuclear interaction
       double rescaleFactor( const double pt, const double cFrac ) const;

       /// Collection of primary PFCandidates to be transmitted to the Event
       std::auto_ptr<reco::PFCandidateCollection> pfC_;
       /// A mask to define the candidates which shall not be transmitted
       std::vector<bool> bMask_;

       /// Parameters
       bool debug_;
       bool bCorrect_;

       /// Calibration parameters for the reconstructed nuclear interactions
       bool bCalibPrimary_;
       std::vector< double > fConst_;
       std::vector< double > fNorm_;
       std::vector< double > fExp_;

       // Maximal accepatble uncertainty on primary tracks to usem them as MC truth for calibration
       double dptRel_PrimaryTrack_;
       double dptRel_MergedTrack_;
       double ptErrorSecondary_;

       /// Useful constants
       static const double pion_mass2;
       static const reco::PFCandidate::Flags fT_TO_DISP_;
       static const reco::PFCandidate::Flags fT_FROM_DISP_;
};

#endif
