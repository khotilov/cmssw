
// Producer for validation histograms for CaloJet objects
// F. Ratnikov, Sept. 7, 2006
// Modified by J F Novak July 10, 2008
// $Id: JPTJetTester.cc,v 1.14 2011/06/30 15:16:00 kovitang Exp $

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

//#include "DataFormats/JetReco/interface/CaloJet.h"
//#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
// #include "DataFormats/METReco/interface/GenMET.h"
// #include "DataFormats/METReco/interface/GenMETCollection.h"
// #include "DataFormats/METReco/interface/MET.h"
// #include "DataFormats/METReco/interface/METCollection.h"

#include "RecoJets/JetProducers/interface/JetMatchingTools.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "JPTJetTester.h"

#include <cmath>

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

using namespace edm;
using namespace reco;
using namespace std;

namespace {
  bool is_B (const reco::Jet& fJet) {return fabs (fJet.eta()) < 1.3;}
  bool is_E (const reco::Jet& fJet) {return fabs (fJet.eta()) >= 1.3 && fabs (fJet.eta()) < 3.;}
  bool is_F (const reco::Jet& fJet) {return fabs (fJet.eta()) >= 3.;}
}

JPTJetTester::JPTJetTester(const edm::ParameterSet& iConfig)
  : mInputCollection (iConfig.getParameter<edm::InputTag>( "src" )),
    mInputGenCollection (iConfig.getParameter<edm::InputTag>( "srcGen" )),
    mOutputFile (iConfig.getUntrackedParameter<std::string>("outputFile", "")),
    mMatchGenPtThreshold (iConfig.getParameter<double>("genPtThreshold")),
    mGenEnergyFractionThreshold (iConfig.getParameter<double>("genEnergyFractionThreshold")),
    mReverseEnergyFractionThreshold (iConfig.getParameter<double>("reverseEnergyFractionThreshold")),
    mRThreshold (iConfig.getParameter<double>("RThreshold")),
     JetCorrectionService  (iConfig.getParameter<std::string>  ("JetCorrectionService"  )),
    mTurnOnEverything (iConfig.getUntrackedParameter<std::string>("TurnOnEverything",""))
{
    numberofevents
    = mEta = mEtaFineBin = mPhi = mPhiFineBin = mE = mE_80 
    = mP = mP_80  = mPt = mPt_80 
    = mMass = mMass_80 
      //    = mConstituents = mConstituents_80
    = mEtaFirst = mPhiFirst  = mPtFirst = mPtFirst_80 = mPtFirst_3000
    = mMjj = mMjj_3000 = mDelEta = mDelPhi = mDelPt 
    /*  
    = mMaxEInEmTowers = mMaxEInHadTowers 
    = mHadEnergyInHO = mHadEnergyInHB = mHadEnergyInHF = mHadEnergyInHE 
    = mHadEnergyInHO_80 = mHadEnergyInHB_80 = mHadEnergyInHE_80 
    = mHadEnergyInHO_3000 = mHadEnergyInHB_3000 = mHadEnergyInHE_3000 
    = mEmEnergyInEB = mEmEnergyInEE = mEmEnergyInHF 
    = mEmEnergyInEB_80 = mEmEnergyInEE_80
    = mEmEnergyInEB_3000 = mEmEnergyInEE_3000
    = mEnergyFractionHadronic = mEnergyFractionEm 
    = mHFLong = mHFTotal = mHFLong_80 = mHFLong_3000 = mHFShort = mHFShort_80 = mHFShort_3000
    = mN90
    */
    //
    //= mCaloMEx = mCaloMEx_3000 = mCaloMEy = mCaloMEy_3000 = mCaloMETSig = mCaloMETSig_3000
    //= mCaloMET = mCaloMET_3000 =  mCaloMETPhi = mCaloSumET  = mCaloSumET_3000   
    = mHadTiming = mEmTiming 
    = mNJetsEtaC = mNJetsEtaF = mNJets1 = mNJets2
      //= mAllGenJetsPt = mMatchedGenJetsPt = mAllGenJetsEta = mMatchedGenJetsEta 
      //= mGenJetMatchEnergyFraction = mReverseMatchEnergyFraction = mRMatch
      = mDeltaEta = mDeltaPhi //= mEScale = mlinEScale = mDeltaE
    = mHadEnergyProfile = mEmEnergyProfile = mJetEnergyProfile 
    /*
    = mHadJetEnergyProfile = mEMJetEnergyProfile
    */
    = mEScale_pt10 = mEScaleFineBin
      //= mpTScaleB_s = mpTScaleE_s = mpTScaleF_s 
    = mpTScaleB_d = mpTScaleE_d = mpTScaleF_d
      = mpTScalePhiB_d = mpTScalePhiE_d = mpTScalePhiF_d
      //= mpTScale_60_120_s = mpTScale_200_300_s = mpTScale_600_900_s = mpTScale_2700_3500_s
    = mpTScale_60_120_d = mpTScale_200_300_d = mpTScale_600_900_d = mpTScale_2700_3500_d
      
    = mpTScale1DB_60_120    = mpTScale1DE_60_120    = mpTScale1DF_60_120 
    = mpTScale1DB_200_300   = mpTScale1DE_200_300   = mpTScale1DF_200_300 
    = mpTScale1DB_600_900   = mpTScale1DE_600_900   = mpTScale1DF_600_900 
    = mpTScale1DB_2700_3500 = mpTScale1DE_2700_3500 = mpTScale1DF_2700_3500
	/*
    = mpTScale1D_60_120 = mpTScale1D_200_300 = mpTScale1D_600_900 = mpTScale1D_2700_3500
    = mHBEne = mHBTime = mHEEne = mHETime = mHFEne = mHFTime = mHOEne = mHOTime
    = mEBEne = mEBTime = mEEEne = mEETime 
      */
    = mPthat_80 = mPthat_3000

      //Corr Jet
    = mCorrJetPt = mCorrJetPt =mCorrJetPt_80  =mCorrJetEta =mCorrJetPhi =mpTRatio =mpTResponse
      = mpTRatioB_d = mpTRatioE_d = mpTRatioF_d
      = mpTRatio_60_120_d = mpTRatio_200_300_d = mpTRatio_600_900_d = mpTRatio_2700_3500_d
      = mpTResponseB_d = mpTResponseE_d = mpTResponseF_d
      = mpTResponse_60_120_d = mpTResponse_200_300_d = mpTResponse_600_900_d = mpTResponse_2700_3500_d
      = mpTResponse_30_d
    = 0;
  
  DQMStore* dbe = &*edm::Service<DQMStore>();
  if (dbe) {
    dbe->setCurrentFolder("JetMET/RecoJetsV/JPTJetTask_" + mInputCollection.label());
    //
    numberofevents    = dbe->book1D("numberofevents","numberofevents", 3, 0 , 2);
    //
    mEta              = dbe->book1D("Eta", "Eta", 120, -6, 6); 
    mEtaFineBin       = dbe->book1D("EtaFineBin_Pt10", "EtaFineBin_Pt10", 600, -6, 6); 
    /*
    mEtaFineBin1p     = dbe->book1D("EtaFineBin1p_Pt10", "EtaFineBin1p_Pt10", 100, 0, 1.3); 
    mEtaFineBin2p     = dbe->book1D("EtaFineBin2p_Pt10", "EtaFineBin2p_Pt10", 100, 1.3, 3); 
    mEtaFineBin3p     = dbe->book1D("EtaFineBin3p_Pt10", "EtaFineBin3p_Pt10", 100, 3, 5); 
    mEtaFineBin1m     = dbe->book1D("EtaFineBin1m_Pt10", "EtaFineBin1m_Pt10", 100, -1.3, 0); 
    mEtaFineBin2m     = dbe->book1D("EtaFineBin2m_Pt10", "EtaFineBin2m_Pt10", 100, -3, -1.3); 
    mEtaFineBin3m     = dbe->book1D("EtaFineBin3m_Pt10", "EtaFineBin3m_Pt10", 100, -5, -3); 
    */
    //
    mPhi              = dbe->book1D("Phi", "Phi", 70, -3.5, 3.5); 
    mPhiFineBin       = dbe->book1D("PhiFineBin_Pt10", "PhiFineBin_Pt10", 350, -3.5, 3.5); 
    //
    mE                = dbe->book1D("E", "E", 100, 0, 500); 
    mE_80             = dbe->book1D("E_80", "E_80", 100, 0, 5000); 
    //
    mP                = dbe->book1D("P", "P", 100, 0, 500); 
    mP_80             = dbe->book1D("P_80", "P_80", 100, 0, 5000); 
    //
    mPt               = dbe->book1D("Pt", "Pt", 100, 0, 150); 
    mPt_80            = dbe->book1D("Pt_80", "Pt_80", 100, 0, 4000);
    //
    mMass             = dbe->book1D("Mass", "Mass", 100, 0, 200); 
    mMass_80          = dbe->book1D("Mass_80", "Mass_80", 100, 0, 500); 
    //
    //    mConstituents     = dbe->book1D("Constituents", "# of Constituents", 100, 0, 100); 
    //    mConstituents_80  = dbe->book1D("Constituents_80", "# of Constituents_80", 40, 0, 40); 
    //
    mEtaFirst         = dbe->book1D("EtaFirst", "EtaFirst", 120, -6, 6); 
    mPhiFirst         = dbe->book1D("PhiFirst", "PhiFirst", 70, -3.5, 3.5);      
    mPtFirst          = dbe->book1D("PtFirst", "PtFirst", 100, 0, 50); 
    mPtFirst_80       = dbe->book1D("PtFirst_80", "PtFirst_80", 100, 0, 140);
    mPtFirst_3000     = dbe->book1D("PtFirst_3000", "PtFirst_3000", 100, 0, 4000);
    //
    mMjj              = dbe->book1D("Mjj", "Mjj", 100, 0, 2000); 
    mMjj_3000         = dbe->book1D("Mjj_3000", "Mjj_3000", 100, 0, 10000); 
    mDelEta           = dbe->book1D("DelEta", "DelEta", 100, -.5, .5); 
    mDelPhi           = dbe->book1D("DelPhi", "DelPhi", 100, -.5, .5); 
    mDelPt            = dbe->book1D("DelPt", "DelPt", 100, -1, 1); 
    //

    /*
    mMaxEInEmTowers   = dbe->book1D("MaxEInEmTowers", "MaxEInEmTowers", 100, 0, 100); 
    mMaxEInHadTowers  = dbe->book1D("MaxEInHadTowers", "MaxEInHadTowers", 100, 0, 100); 
    mHadEnergyInHO    = dbe->book1D("HadEnergyInHO", "HadEnergyInHO", 100, 0, 10); 
    mHadEnergyInHB    = dbe->book1D("HadEnergyInHB", "HadEnergyInHB", 100, 0, 50); 
    mHadEnergyInHF    = dbe->book1D("HadEnergyInHF", "HadEnergyInHF", 100, 0, 50); 
    mHadEnergyInHE    = dbe->book1D("HadEnergyInHE", "HadEnergyInHE", 100, 0, 100); 
    //
    mHadEnergyInHO_80    = dbe->book1D("HadEnergyInHO_80", "HadEnergyInHO_80", 100, 0, 50); 
    mHadEnergyInHB_80    = dbe->book1D("HadEnergyInHB_80", "HadEnergyInHB_80", 100, 0, 200); 
    mHadEnergyInHE_80    = dbe->book1D("HadEnergyInHE_80", "HadEnergyInHE_80", 100, 0, 1000); 
    mHadEnergyInHO_3000  = dbe->book1D("HadEnergyInHO_3000", "HadEnergyInHO_3000", 100, 0, 500); 
    mHadEnergyInHB_3000  = dbe->book1D("HadEnergyInHB_3000", "HadEnergyInHB_3000", 100, 0, 3000); 
    mHadEnergyInHE_3000  = dbe->book1D("HadEnergyInHE_3000", "HadEnergyInHE_3000", 100, 0, 2000); 
    //
    mEmEnergyInEB     = dbe->book1D("EmEnergyInEB", "EmEnergyInEB", 100, 0, 50); 
    mEmEnergyInEE     = dbe->book1D("EmEnergyInEE", "EmEnergyInEE", 100, 0, 50); 
    mEmEnergyInHF     = dbe->book1D("EmEnergyInHF", "EmEnergyInHF", 120, -20, 100); 
    mEmEnergyInEB_80  = dbe->book1D("EmEnergyInEB_80", "EmEnergyInEB_80", 100, 0, 200); 
    mEmEnergyInEE_80  = dbe->book1D("EmEnergyInEE_80", "EmEnergyInEE_80", 100, 0, 1000); 
    mEmEnergyInEB_3000= dbe->book1D("EmEnergyInEB_3000", "EmEnergyInEB_3000", 100, 0, 3000); 
    mEmEnergyInEE_3000= dbe->book1D("EmEnergyInEE_3000", "EmEnergyInEE_3000", 100, 0, 2000); 
    mEnergyFractionHadronic = dbe->book1D("EnergyFractionHadronic", "EnergyFractionHadronic", 120, -0.1, 1.1); 
    mEnergyFractionEm = dbe->book1D("EnergyFractionEm", "EnergyFractionEm", 120, -0.1, 1.1); 
    //
    mHFTotal          = dbe->book1D("HFTotal", "HFTotal", 100, 0, 500);
    mHFTotal_80       = dbe->book1D("HFTotal_80", "HFTotal_80", 100, 0, 3000);
    mHFTotal_3000     = dbe->book1D("HFTotal_3000", "HFTotal_3000", 100, 0, 6000);
    mHFLong           = dbe->book1D("HFLong", "HFLong", 100, 0, 500);
    mHFLong_80        = dbe->book1D("HFLong_80", "HFLong_80", 100, 0, 200);
    mHFLong_3000      = dbe->book1D("HFLong_3000", "HFLong_3000", 100, 0, 1500);
    mHFShort          = dbe->book1D("HFShort", "HFShort", 100, 0, 500);
    mHFShort_80       = dbe->book1D("HFShort_80", "HFShort_80", 100, 0, 200);
    mHFShort_3000     = dbe->book1D("HFShort_3000", "HFShort_3000", 100, 0, 1500);
    //
    mN90              = dbe->book1D("N90", "N90", 50, 0, 50); 
    */
    //
    mGenEta           = dbe->book1D("GenEta", "GenEta", 120, -6, 6);
    mGenPhi           = dbe->book1D("GenPhi", "GenPhi", 70, -3.5, 3.5);
    mGenPt            = dbe->book1D("GenPt", "GenPt", 100, 0, 150);
    mGenPt_80         = dbe->book1D("GenPt_80", "GenPt_80", 100, 0, 1500);
    //
    mGenEtaFirst      = dbe->book1D("GenEtaFirst", "GenEtaFirst", 100, -5, 5);
    mGenPhiFirst      = dbe->book1D("GenPhiFirst", "GenPhiFirst", 70, -3.5, 3.5);
    //
    /*
    mCaloMEx          = dbe->book1D("CaloMEx","CaloMEx",200,-150,150);
    mCaloMEx_3000     = dbe->book1D("CaloMEx_3000","CaloMEx_3000",100,-500,500);
    mCaloMEy          = dbe->book1D("CaloMEy","CaloMEy",200,-150,150);
    mCaloMEy_3000     = dbe->book1D("CaloMEy_3000","CaloMEy_3000",100,-500,500);
    mCaloMETSig       = dbe->book1D("CaloMETSig","CaloMETSig",100,0,15);
    mCaloMETSig_3000  = dbe->book1D("CaloMETSig_3000","CaloMETSig_3000",100,0,50);
    mCaloMET          = dbe->book1D("CaloMET","CaloMET",100,0,150);
    mCaloMET_3000     = dbe->book1D("CaloMET_3000","CaloMET_3000",100,0,1000);
    mCaloMETPhi       = dbe->book1D("CaloMETPhi","CaloMETPhi",70, -3.5, 3.5);
    mCaloSumET        = dbe->book1D("CaloSumET","CaloSumET",100,0,500);
    mCaloSumET_3000   = dbe->book1D("CaloSumET_3000","CaloSumET_3000",100,3000,8000);
    */
    //
    mHadTiming        = dbe->book1D("HadTiming", "HadTiming", 75, -50, 100);
    mEmTiming         = dbe->book1D("EMTiming", "EMTiming", 75, -50, 100);
    //
    mNJetsEtaC        = dbe->book1D("NJetsEtaC_Pt10", "NJetsEtaC_Pt10", 15, 0, 15);
    mNJetsEtaF        = dbe->book1D("NJetsEtaF_Pt10", "NJetsEtaF_Pt10", 15, 0, 15);
    //
    mNJets1           = dbe->bookProfile("NJets1", "NJets1", 100, 0, 200,  100, 0, 50, "s");
    mNJets2           = dbe->bookProfile("NJets2", "NJets2", 100, 0, 4000, 100, 0, 50, "s");
    //
    /*
    mHBEne     = dbe->book1D( "HBEne",  "HBEne", 1000, -20, 100 );
    mHBTime    = dbe->book1D( "HBTime", "HBTime", 200, -200, 200 );
    mHEEne     = dbe->book1D( "HEEne",  "HEEne", 1000, -20, 100 );
    mHETime    = dbe->book1D( "HETime", "HETime", 200, -200, 200 );
    mHOEne     = dbe->book1D( "HOEne",  "HOEne", 1000, -20, 100 );
    mHOTime    = dbe->book1D( "HOTime", "HOTime", 200, -200, 200 );
    mHFEne     = dbe->book1D( "HFEne",  "HFEne", 1000, -20, 100 );
    mHFTime    = dbe->book1D( "HFTime", "HFTime", 200, -200, 200 );
    mEBEne     = dbe->book1D( "EBEne",  "EBEne", 1000, -20, 100 );
    mEBTime    = dbe->book1D( "EBTime", "EBTime", 200, -200, 200 );
    mEEEne     = dbe->book1D( "EEEne",  "EEEne", 1000, -20, 100 );
    mEETime    = dbe->book1D( "EETime", "EETime", 200, -200, 200 );
    */
    //
    mPthat_80            = dbe->book1D("Pthat_80", "Pthat_80", 100, 0.0, 1000.0); 
    mPthat_3000          = dbe->book1D("Pthat_3000", "Pthat_3000", 100, 1000.0, 4000.0); 

    //Corr
    mCorrJetPt  = dbe->book1D("CorrPt", "CorrPt", 100, 0, 150);
    mCorrJetPt_80 = dbe->book1D("CorrPt_80", "CorrPt_80", 100, 0, 4000);
    mCorrJetEta = dbe->book1D("CorrEta", "CorrEta", 120, -6, 6);
    mCorrJetPhi = dbe->book1D("CorrPhi", "CorrPhi", 70, -3.5, 3.5);

    //
    double log10PtMin = 0.5; //=3.1622766
    double log10PtMax = 3.75; //=5623.41325
    int log10PtBins = 26; 
    double etaMin = -6.;
    double etaMax = 6.;
    int etaBins = 50;

    //double linPtMin = 5;
    //double linPtMax = 155;
    //int linPtBins = 15;

    int log10PtFineBins = 50;
    /*
    mAllGenJetsPt = dbe->book1D("GenJetLOGpT", "GenJet LOG(pT_gen)", 
				log10PtBins, log10PtMin, log10PtMax);
    mMatchedGenJetsPt = dbe->book1D("MatchedGenJetLOGpT", "MatchedGenJet LOG(pT_gen)", 
				    log10PtBins, log10PtMin, log10PtMax);
    mAllGenJetsEta = dbe->book2D("GenJetEta", "GenJet Eta vs LOG(pT_gen)", 
				 log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax);
    mMatchedGenJetsEta = dbe->book2D("MatchedGenJetEta", "MatchedGenJet Eta vs LOG(pT_gen)", 
				     log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax);
    */
    //
    if (mTurnOnEverything.compare("yes")==0) {
      /*
      mHadEnergyProfile = dbe->bookProfile2D("HadEnergyProfile", "HadEnergyProfile", 82, -41, 41, 73, 0, 73, 100, 0, 10000, "s");
      mEmEnergyProfile  = dbe->bookProfile2D("EmEnergyProfile", "EmEnergyProfile", 82, -41, 41, 73, 0, 73, 100, 0, 10000, "s");
      */
    }
    /*
    mJetEnergyProfile = dbe->bookProfile2D("JetEnergyProfile", "JetEnergyProfile", 50, -5, 5, 36, -3.1415987, 3.1415987, 100, 0, 10000, "s");
    mHadJetEnergyProfile = dbe->bookProfile2D("HadJetEnergyProfile", "HadJetEnergyProfile", 50, -5, 5, 36, -3.1415987, 3.1415987, 100, 0, 10000, "s");
    mEMJetEnergyProfile = dbe->bookProfile2D("EMJetEnergyProfile", "EMJetEnergyProfile", 50, -5, 5, 36, -3.1415987, 3.1415987, 100, 0, 10000, "s");
    */

    //
    if (mTurnOnEverything.compare("yes")==0) {
      /*
      mGenJetMatchEnergyFraction  = dbe->book3D("GenJetMatchEnergyFraction", "GenJetMatchEnergyFraction vs LOG(pT_gen) vs eta", 
						log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax, 101, 0, 1.01);
      mReverseMatchEnergyFraction  = dbe->book3D("ReverseMatchEnergyFraction", "ReverseMatchEnergyFraction vs LOG(pT_gen) vs eta", 
					       log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax, 101, 0, 1.01);
      mRMatch  = dbe->book3D("RMatch", "delta(R)(Gen-Calo) vs LOG(pT_gen) vs eta", 
			     log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax, 60, 0, 3);
      */
      mDeltaEta = dbe->book3D("DeltaEta", "DeltaEta vs LOG(pT_gen) vs eta", 
			      log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax, 100, -1, 1);
      mDeltaPhi = dbe->book3D("DeltaPhi", "DeltaPhi vs LOG(pT_gen) vs eta", 
			      log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax, 100, -1, 1);
      /*
      mEScale = dbe->book3D("EScale", "EnergyScale vs LOG(pT_gen) vs eta", 
			    log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax, 100, 0, 2);
      mlinEScale = dbe->book3D("linEScale", "EnergyScale vs LOG(pT_gen) vs eta", 
			       linPtBins, linPtMin, linPtMax, etaBins, etaMin, etaMax, 100, 0, 2);
      mDeltaE = dbe->book3D("DeltaE", "DeltaE vs LOG(pT_gen) vs eta", 
			    log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax, 2000, -200, 200);
      */
      //
      mEScale_pt10 = dbe->book3D("EScale_pt10", "EnergyScale vs LOG(pT_gen) vs eta", 
				 log10PtBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax, 100, 0, 2);
      mEScaleFineBin = dbe->book3D("EScaleFineBins", "EnergyScale vs LOG(pT_gen) vs eta", 
				   log10PtFineBins, log10PtMin, log10PtMax, etaBins, etaMin, etaMax, 100, 0, 2);
    }
    /*
    mpTScaleB_s = dbe->bookProfile("pTScaleB_s", "pTScale_s_0<|eta|<1.3",
				    log10PtBins, log10PtMin, log10PtMax, 0, 2, "s");
    mpTScaleE_s = dbe->bookProfile("pTScaleE_s", "pTScale_s_1.3<|eta|<3.0",
				    log10PtBins, log10PtMin, log10PtMax, 0, 2, "s");
    mpTScaleF_s = dbe->bookProfile("pTScaleF_s", "pTScale_s_3.0<|eta|<5.0",
				    log10PtBins, log10PtMin, log10PtMax, 0, 2, "s");
    */
    mpTScaleB_d = dbe->bookProfile("pTScaleB_d", "pTScale_d_0<|eta|<1.3",
				   log10PtBins, log10PtMin, log10PtMax, 0, 2, " ");
    mpTScaleE_d = dbe->bookProfile("pTScaleE_d", "pTScale_d_1.3<|eta|<3.0",
				   log10PtBins, log10PtMin, log10PtMax, 0, 2, " ");
    mpTScaleF_d = dbe->bookProfile("pTScaleF_d", "pTScale_d_3.0<|eta|<6.0",
				   log10PtBins, log10PtMin, log10PtMax, 0, 2, " ");
    mpTScalePhiB_d = dbe->bookProfile("pTScalePhiB_d", "pTScalePhi_d_0<|eta|<1.3",
				   log10PtBins, log10PtMin, log10PtMax, 0, 2, " ");
    mpTScalePhiE_d = dbe->bookProfile("pTScalePhiE_d", "pTScalePhi_d_1.3<|eta|<3.0",
				   log10PtBins, log10PtMin, log10PtMax, 0, 2, " ");
    mpTScalePhiF_d = dbe->bookProfile("pTScalePhiF_d", "pTScalePhi_d_3.0<|eta|<6.0",
				   log10PtBins, log10PtMin, log10PtMax, 0, 2, " ");
    /*
    mpTScale_60_120_s    = dbe->bookProfile("pTScale_60_120_s", "pTScale_s_60<pT<120",
					  etaBins, etaMin, etaMax, 0., 2., "s");
    mpTScale_200_300_s   = dbe->bookProfile("pTScale_200_300_s", "pTScale_s_200<pT<300",
					  etaBins, etaMin, etaMax, 0., 2., "s");
    mpTScale_600_900_s   = dbe->bookProfile("pTScale_600_900_s", "pTScale_s_600<pT<900",
					  etaBins, etaMin, etaMax, 0., 2., "s");
    mpTScale_2700_3500_s = dbe->bookProfile("pTScale_2700_3500_s", "pTScale_s_2700<pt<3500",
    etaBins, etaMin, etaMax, 0., 2., "s");
*/
    mpTScale_60_120_d    = dbe->bookProfile("pTScale_60_120_d", "pTScale_d_60<pT<120",
					  etaBins, etaMin, etaMax, 0., 2., " ");
    mpTScale_200_300_d   = dbe->bookProfile("pTScale_200_300_d", "pTScale_d_200<pT<300",
					  etaBins, etaMin, etaMax, 0., 2., " ");
    mpTScale_600_900_d   = dbe->bookProfile("pTScale_600_900_d", "pTScale_d_600<pT<900",
					  etaBins, etaMin, etaMax, 0., 2., " ");
    mpTScale_2700_3500_d = dbe->bookProfile("pTScale_2700_3500_d", "pTScale_d_2700<pt<3500",
                                          etaBins, etaMin, etaMax, 0., 2., " ");
    
    mpTScale1DB_60_120 = dbe->book1D("pTScale1DB_60_120", "pTScale_distribution_for_0<|eta|<1.3_60_120",
				   100, 0, 2);
    mpTScale1DE_60_120 = dbe->book1D("pTScale1DE_60_120", "pTScale_distribution_for_1.3<|eta|<3.0_60_120",
				   50, 0, 2);
    mpTScale1DF_60_120 = dbe->book1D("pTScale1DF_60_120", "pTScale_distribution_for_3.0<|eta|<6.0_60_120",
				   50, 0, 2);

    mpTScale1DB_200_300 = dbe->book1D("pTScale1DB_200_300", "pTScale_distribution_for_0<|eta|<1.3_200_300",
				   100, 0, 2);
    mpTScale1DE_200_300 = dbe->book1D("pTScale1DE_200_300", "pTScale_distribution_for_1.3<|eta|<3.0_200_300",
				   50, 0, 2);
    mpTScale1DF_200_300 = dbe->book1D("pTScale1DF_200_300", "pTScale_distribution_for_3.0<|eta|<6.0_200_300",
				   50, 0, 2);

    mpTScale1DB_600_900 = dbe->book1D("pTScale1DB_600_900", "pTScale_distribution_for_0<|eta|<1.3_600_900",
				   100, 0, 2);
    mpTScale1DE_600_900 = dbe->book1D("pTScale1DE_600_900", "pTScale_distribution_for_1.3<|eta|<3.0_600_900",
				   50, 0, 2);
    mpTScale1DF_600_900 = dbe->book1D("pTScale1DF_600_900", "pTScale_distribution_for_3.0<|eta|<6.0_600_900",
				   50, 0, 2);

    mpTScale1DB_2700_3500 = dbe->book1D("pTScale1DB_2700_3500", "pTScale_distribution_for_0<|eta|<1.3_2700_3500",
				   100, 0, 2);
    mpTScale1DE_2700_3500 = dbe->book1D("pTScale1DE_2700_3500", "pTScale_distribution_for_1.3<|eta|<3.0_2700_3500",
				   50, 0, 2);
    mpTScale1DF_2700_3500 = dbe->book1D("pTScale1DF_2700_3500", "pTScale_distribution_for_3.0<|eta|<6.0_2700_3500",
				   50, 0, 2);
/*
    mpTScale1D_60_120    = dbe->book1D("pTScale1D_60_120", "pTScale_distribution_for_60<pT<120",
					    100, 0, 2);
    mpTScale1D_200_300    = dbe->book1D("pTScale1D_200_300", "pTScale_distribution_for_200<pT<300",
					    100, 0, 2);
    mpTScale1D_600_900    = dbe->book1D("pTScale1D_600_900", "pTScale_distribution_for_600<pT<900",
					    100, 0, 2);
    mpTScale1D_2700_3500 = dbe->book1D("pTScale1D_2700_3500", "pTScale_distribution_for_2700<pt<3500",
					    100, 0, 2);
    */
 ///////////Corr profile//////////////
    mpTRatio = dbe->bookProfile("pTRatio", "pTRatio",
                                log10PtBins, log10PtMin, log10PtMax, 100, 0.,5., " ");
    mpTRatioB_d = dbe->bookProfile("pTRatioB_d", "pTRatio_d_0<|eta|<1.3",
                                   log10PtBins, log10PtMin, log10PtMax, 0, 5, " ");
    mpTRatioE_d = dbe->bookProfile("pTRatioE_d", "pTRatio_d_1.3<|eta|<3.0",
                                   log10PtBins, log10PtMin, log10PtMax, 0, 5, " ");
    mpTRatioF_d = dbe->bookProfile("pTRatioF_d", "pTRatio_d_3.0<|eta|<6.0",
                                   log10PtBins, log10PtMin, log10PtMax, 0, 5, " ");
    mpTRatio_60_120_d    = dbe->bookProfile("pTRatio_60_120_d", "pTRatio_d_60<pT<120",
                                          etaBins, etaMin, etaMax, 0., 5., " ");
    mpTRatio_200_300_d   = dbe->bookProfile("pTRatio_200_300_d", "pTRatio_d_200<pT<300",
                                          etaBins, etaMin, etaMax, 0., 5., " ");
    mpTRatio_600_900_d   = dbe->bookProfile("pTRatio_600_900_d", "pTRatio_d_600<pT<900",
                                          etaBins, etaMin, etaMax, 0., 5., " ");
    mpTRatio_2700_3500_d = dbe->bookProfile("pTRatio_2700_3500_d", "pTRatio_d_2700<pt<3500",
                                          etaBins, etaMin, etaMax, 0., 5., " ");
    mpTResponse = dbe->bookProfile("pTResponse", "pTResponse",
				log10PtBins, log10PtMin, log10PtMax, 100, 0.8,1.2, " ");
    mpTResponseB_d = dbe->bookProfile("pTResponseB_d", "pTResponse_d_0<|eta|<1.3",
				      log10PtBins, log10PtMin, log10PtMax, 0.8, 1.2, " ");
    mpTResponseE_d = dbe->bookProfile("pTResponseE_d", "pTResponse_d_1.3<|eta|<3.0",
                                   log10PtBins, log10PtMin, log10PtMax, 0.8, 1.2, " ");
    mpTResponseF_d = dbe->bookProfile("pTResponseF_d", "pTResponse_d_3.0<|eta|<6.0",
                                   log10PtBins, log10PtMin, log10PtMax, 0.8, 1.2, " ");
    mpTResponse_60_120_d    = dbe->bookProfile("pTResponse_60_120_d", "pTResponse_d_60<pT<120",
                                          etaBins, etaMin, etaMax, 0.8, 1.2, " ");
    mpTResponse_200_300_d   = dbe->bookProfile("pTResponse_200_300_d", "pTResponse_d_200<pT<300",
                                          etaBins, etaMin, etaMax, 0.8, 1.2, " ");
    mpTResponse_600_900_d   = dbe->bookProfile("pTResponse_600_900_d", "pTResponse_d_600<pT<900",
                                          etaBins, etaMin, etaMax, 0.8, 1.2, " ");
    mpTResponse_2700_3500_d = dbe->bookProfile("pTResponse_2700_3500_d", "pTResponse_d_2700<pt<3500",
					       etaBins, etaMin, etaMax, 0.8, 1.2, " ");
    mpTResponse_30_d = dbe->bookProfile("pTResponse_30_d", "pTResponse_d_pt>30",
					       etaBins, etaMin, etaMax, 0.8, 1.2, " ");
  } // if (dbe)

  if (mOutputFile.empty ()) {
    LogInfo("OutputInfo") << " JPTJet histograms will NOT be saved";
  } 
  else {
    LogInfo("OutputInfo") << " JPTJethistograms will be saved to file:" << mOutputFile;
  }
}
   
JPTJetTester::~JPTJetTester()
{
}

void JPTJetTester::beginJob(){
}

void JPTJetTester::endJob() {
 if (!mOutputFile.empty() && &*edm::Service<DQMStore>()) edm::Service<DQMStore>()->save (mOutputFile);
}


void JPTJetTester::analyze(const edm::Event& mEvent, const edm::EventSetup& mSetup)
{
  double countsfornumberofevents = 1;
  numberofevents->Fill(countsfornumberofevents);
  // *********************************
  // *** Get pThat
  // *********************************
if (!mEvent.isRealData()){
  edm::Handle<HepMCProduct> evt;
  mEvent.getByLabel("generator", evt);
  if (evt.isValid()) {
  HepMC::GenEvent * myGenEvent = new HepMC::GenEvent(*(evt->GetEvent()));
  
  double pthat = myGenEvent->event_scale();

  mPthat_80->Fill(pthat);
  mPthat_3000->Fill(pthat);

  delete myGenEvent; 
  }
}
  // ***********************************
  // *** Get CaloMET
  // ***********************************

  const CaloMET *calomet;
  edm::Handle<CaloMETCollection> calo;
  mEvent.getByLabel("met", calo);
  if (!calo.isValid()) {
    edm::LogInfo("OutputInfo") << " failed to retrieve data required by MET Task";
    edm::LogInfo("OutputInfo") << " MET Task cannot continue...!";
  } else {
    const CaloMETCollection *calometcol = calo.product();
    calomet = &(calometcol->front());
    /*
    double caloSumET = calomet->sumEt();
    double caloMETSig = calomet->mEtSig();
    double caloMET = calomet->pt();
    double caloMEx = calomet->px();
    double caloMEy = calomet->py();
    double caloMETPhi = calomet->phi();

    mCaloMEx->Fill(caloMEx);
    mCaloMEx_3000->Fill(caloMEx);
    mCaloMEy->Fill(caloMEy);
    mCaloMEy_3000->Fill(caloMEy);
    mCaloMET->Fill(caloMET);
    mCaloMET_3000->Fill(caloMET);
    mCaloMETPhi->Fill(caloMETPhi);
    mCaloSumET->Fill(caloSumET);
    mCaloSumET_3000->Fill(caloSumET);
    mCaloMETSig->Fill(caloMETSig);
    mCaloMETSig_3000->Fill(caloMETSig);
    */
  }

  // ***********************************
  // *** Get the CaloTower collection
  // ***********************************
  Handle<CaloTowerCollection> caloTowers;
  mEvent.getByLabel( "towerMaker", caloTowers );
  if (caloTowers.isValid()) {
  for( CaloTowerCollection::const_iterator cal = caloTowers->begin(); cal != caloTowers->end(); ++ cal ){

    //To compensate for the index
    if (mTurnOnEverything.compare("yes")==0) {
      /*
      if (cal->ieta() >> 0 ){mHadEnergyProfile->Fill (cal->ieta()-1, cal->iphi(), cal->hadEnergy());
      mEmEnergyProfile->Fill (cal->ieta()-1, cal->iphi(), cal->emEnergy());}
      mHadEnergyProfile->Fill (cal->ieta(), cal->iphi(), cal->hadEnergy());
      mEmEnergyProfile->Fill (cal->ieta(), cal->iphi(), cal->emEnergy());
      */  
  }

    mHadTiming->Fill (cal->hcalTime());
    mEmTiming->Fill (cal->ecalTime());    
  }
  }
  
  // ***********************************
  // *** Get the RecHits collection
  // ***********************************
  try {
    std::vector<edm::Handle<HBHERecHitCollection> > colls;
    mEvent.getManyByType(colls);
    std::vector<edm::Handle<HBHERecHitCollection> >::iterator i;
    for (i=colls.begin(); i!=colls.end(); i++) {
      for (HBHERecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) {
        //      std::cout << *j << std::endl;
	/*
        if (j->id().subdet() == HcalBarrel) {
          mHBEne->Fill(j->energy()); 
          mHBTime->Fill(j->time()); 
        }
        if (j->id().subdet() == HcalEndcap) {
          mHEEne->Fill(j->energy()); 
          mHETime->Fill(j->time()); 
        }
	*/
      }
    }
  } catch (...) {
    edm::LogInfo("OutputInfo") << " No HB/HE RecHits.";
  }
  
  try {
    std::vector<edm::Handle<HFRecHitCollection> > colls;
    mEvent.getManyByType(colls);
    std::vector<edm::Handle<HFRecHitCollection> >::iterator i;
    for (i=colls.begin(); i!=colls.end(); i++) {
      for (HFRecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) {
        //      std::cout << *j << std::endl;
	/*
        if (j->id().subdet() == HcalForward) {
          mHFEne->Fill(j->energy()); 
          mHFTime->Fill(j->time()); 
        }
	*/
      }
    }
  } catch (...) {
    edm::LogInfo("OutputInfo") << " No HF RecHits.";
  }

  try {
    std::vector<edm::Handle<HORecHitCollection> > colls;
    mEvent.getManyByType(colls);
    std::vector<edm::Handle<HORecHitCollection> >::iterator i;
    for (i=colls.begin(); i!=colls.end(); i++) {
      for (HORecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) {
	/*
        if (j->id().subdet() == HcalOuter) {
          mHOEne->Fill(j->energy()); 
          mHOTime->Fill(j->time()); 
        }
	*/
      }
    }
  } catch (...) {
    edm::LogInfo("OutputInfo") << " No HO RecHits.";
  }
  try {
    std::vector<edm::Handle<EBRecHitCollection> > colls;
    mEvent.getManyByType(colls);
    std::vector<edm::Handle<EBRecHitCollection> >::iterator i;
    for (i=colls.begin(); i!=colls.end(); i++) {
      for (EBRecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) {
        //      if (j->id() == EcalBarrel) {
	//mEBEne->Fill(j->energy()); 
	//mEBTime->Fill(j->time()); 
	//    }
        //      std::cout << *j << std::endl;
        //      std::cout << j->id() << std::endl;
      }
    }
  } catch (...) {
    edm::LogInfo("OutputInfo") << " No EB RecHits.";
  }

  try {
    std::vector<edm::Handle<EERecHitCollection> > colls;
    mEvent.getManyByType(colls);
    std::vector<edm::Handle<EERecHitCollection> >::iterator i;
    for (i=colls.begin(); i!=colls.end(); i++) {
      for (EERecHitCollection::const_iterator j=(*i)->begin(); j!=(*i)->end(); j++) {
        //      if (j->id().subdet() == EcalEndcap) {
	//mEEEne->Fill(j->energy()); 
	//mEETime->Fill(j->time()); 
	//    }
	//      std::cout << *j << std::endl;
      }
    }
  } catch (...) {
    edm::LogInfo("OutputInfo") << " No EE RecHits.";
  }

  //***********************************
  //*** Get the Jet collection
  //***********************************
  math::XYZTLorentzVector p4tmp[2];
  Handle<JPTJetCollection> jptJets;
  mEvent.getByLabel(mInputCollection, jptJets);
  if (!jptJets.isValid()) return;
  JPTJetCollection::const_iterator jet = jptJets->begin ();
  int jetIndex = 0;
  int nJet = 0;
  int nJetF = 0;
  int nJetC = 0;
  for (; jet != jptJets->end (); jet++, jetIndex++) {
    if (jet->pt() > 10.) {
      if (fabs(jet->eta()) > 1.3) 
	nJetF++;
      else 
	nJetC++;	  
    }
    if (jet->pt() > 10.) {
      if (mEta) mEta->Fill (jet->eta());
      if (mEtaFineBin) mEtaFineBin->Fill (jet->eta());
      //if (mEtaFineBin1p) mEtaFineBin1p->Fill (jet->eta());
      //if (mEtaFineBin2p) mEtaFineBin2p->Fill (jet->eta());
      //if (mEtaFineBin3p) mEtaFineBin3p->Fill (jet->eta());
      //if (mEtaFineBin1m) mEtaFineBin1m->Fill (jet->eta());
      //if (mEtaFineBin2m) mEtaFineBin2m->Fill (jet->eta());
      //if (mEtaFineBin3m) mEtaFineBin3m->Fill (jet->eta());
      if (mPhiFineBin) mPhiFineBin->Fill (jet->phi());
    }
    if (mPhi) mPhi->Fill (jet->phi());
    if (mE) mE->Fill (jet->energy());
    if (mE_80) mE_80->Fill (jet->energy());
    if (mP) mP->Fill (jet->p());
    if (mP_80) mP_80->Fill (jet->p());
    if (mPt) mPt->Fill (jet->pt());
    if (mPt_80) mPt_80->Fill (jet->pt());
    if (mMass) mMass->Fill (jet->mass());
    if (mMass_80) mMass_80->Fill (jet->mass());
    //    if (mConstituents) mConstituents->Fill (jet->nConstituents());
    //    if (mConstituents_80) mConstituents_80->Fill (jet->nConstituents());
    if (jet == jptJets->begin ()) { // first jet
      if (mEtaFirst) mEtaFirst->Fill (jet->eta());
      if (mPhiFirst) mPhiFirst->Fill (jet->phi());
      if (mPtFirst) mPtFirst->Fill (jet->pt());
      if (mPtFirst_80) mPtFirst_80->Fill (jet->pt());
      if (mPtFirst_3000) mPtFirst_3000->Fill (jet->pt());
    }
    if (jetIndex == 0) {
      nJet++;
      p4tmp[0] = jet->p4();     
    }
    if (jetIndex == 1) {
      nJet++;
      p4tmp[1] = jet->p4();     
    }

    /*
    if (mMaxEInEmTowers) mMaxEInEmTowers->Fill ( (jet->getCaloJetRef()->maxEInEmTowers());
    if (mMaxEInHadTowers) mMaxEInHadTowers->Fill (jet->getCaloJetRef()->maxEInHadTowers());
    if (mHadEnergyInHO) mHadEnergyInHO->Fill (jet->getCaloJetRef()->hadEnergyInHO());
    if (mHadEnergyInHO_80)   mHadEnergyInHO_80->Fill (jet->getCaloJetRef()->hadEnergyInHO());
    if (mHadEnergyInHO_3000) mHadEnergyInHO_3000->Fill (jet->getCaloJetRef()->hadEnergyInHO());
    if (mHadEnergyInHB) mHadEnergyInHB->Fill (jet->getCaloJetRef()->hadEnergyInHB());
    if (mHadEnergyInHB_80)   mHadEnergyInHB_80->Fill (jet->getCaloJetRef()->hadEnergyInHB());
    if (mHadEnergyInHB_3000) mHadEnergyInHB_3000->Fill (jet->getCaloJetRef()->hadEnergyInHB());
    if (mHadEnergyInHF) mHadEnergyInHF->Fill (jet->getCaloJetRef()->hadEnergyInHF());
    if (mHadEnergyInHE) mHadEnergyInHE->Fill (jet->getCaloJetRef()->hadEnergyInHE());
    if (mHadEnergyInHE_80)   mHadEnergyInHE_80->Fill (jet->getCaloJetRef()->hadEnergyInHE());
    if (mHadEnergyInHE_3000) mHadEnergyInHE_3000->Fill (jet->getCaloJetRef()->hadEnergyInHE());
    if (mEmEnergyInEB) mEmEnergyInEB->Fill (jet->getCaloJetRef()->emEnergyInEB());
    if (mEmEnergyInEB_80)   mEmEnergyInEB_80->Fill (jet->getCaloJetRef()->emEnergyInEB());
    if (mEmEnergyInEB_3000) mEmEnergyInEB_3000->Fill (jet->getCaloJetRef()->emEnergyInEB());
    if (mEmEnergyInEE) mEmEnergyInEE->Fill (jet->getCaloJetRef()->emEnergyInEE());
    if (mEmEnergyInEE_80)   mEmEnergyInEE_80->Fill (jet->getCaloJetRef()->emEnergyInEE());
    if (mEmEnergyInEE_3000) mEmEnergyInEE_3000->Fill (jet->getCaloJetRef()->emEnergyInEE());
    if (mEmEnergyInHF) mEmEnergyInHF->Fill (jet->getCaloJetRef()->emEnergyInHF());
    if (mEnergyFractionHadronic) mEnergyFractionHadronic->Fill (jet->getCaloJetRef()->energyFractionHadronic());
    if (mEnergyFractionEm) mEnergyFractionEm->Fill (jet->getCaloJetRef()->emEnergyFraction());

    if (mHFTotal)      mHFTotal->Fill (jet->getCaloJetRef()->hadEnergyInHF()+jet->getCaloJetRef()->emEnergyInHF());
    if (mHFTotal_80)   mHFTotal_80->Fill (jet->getCaloJetRef()->hadEnergyInHF()+jet->getCaloJetRef()->emEnergyInHF());
    if (mHFTotal_3000) mHFTotal_3000->Fill (jet->getCaloJetRef()->hadEnergyInHF()+jet->getCaloJetRef()->emEnergyInHF());
    if (mHFLong)       mHFLong->Fill (jet->getCaloJetRef()->hadEnergyInHF()*0.5+jet->getCaloJetRef()->emEnergyInHF());
    if (mHFLong_80)    mHFLong_80->Fill (jet->getCaloJetRef()->hadEnergyInHF()*0.5+jet->getCaloJetRef()->emEnergyInHF());
    if (mHFLong_3000)  mHFLong_3000->Fill (jet->getCaloJetRef()->hadEnergyInHF()*0.5+jet->getCaloJetRef()->emEnergyInHF());
    if (mHFShort)      mHFShort->Fill (jet->getCaloJetRef()->hadEnergyInHF()*0.5);
    if (mHFShort_80)   mHFShort_80->Fill (jet->getCaloJetRef()->hadEnergyInHF()*0.5);
    if (mHFShort_3000) mHFShort_3000->Fill (jet->getCaloJetRef()->hadEnergyInHF()*0.5);
    */

    /*
    if (mN90) mN90->Fill (jet->getCaloJetRef()->n90());
    */
    //mJetEnergyProfile->Fill (jet->eta(), jet->phi(), jet->energy());
    /*
    mHadJetEnergyProfile->Fill (jet->eta(), jet->phi(), 
				jet->getCaloJetRef()->hadEnergyInHO()+
				jet->getCaloJetRef()->hadEnergyInHB()+
				jet->getCaloJetRef()->hadEnergyInHF()+
				jet->getCaloJetRef()->hadEnergyInHE());
    mEMJetEnergyProfile->Fill (jet->eta(), jet->phi(), 
			       jet->getCaloJetRef()->emEnergyInEB()+
			       jet->getCaloJetRef()->emEnergyInEE()+
			       jet->getCaloJetRef()->emEnergyInHF());
    */
  }

  if (mNJetsEtaC) mNJetsEtaC->Fill( nJetC );
  if (mNJetsEtaF) mNJetsEtaF->Fill( nJetF );

  if (nJet == 2) {
    if (mMjj) mMjj->Fill( (p4tmp[0]+p4tmp[1]).mass() );
    if (mMjj_3000) mMjj_3000->Fill( (p4tmp[0]+p4tmp[1]).mass() );
  }

 // Correction jets
  const JetCorrector* corrector = JetCorrector::getJetCorrector (JetCorrectionService,mSetup);

  for (JPTJetCollection::const_iterator jet = jptJets->begin(); jet !=jptJets ->end(); jet++) 
  {
 
      JPTJet  correctedJet = *jet;
      //double scale = corrector->correction(jet->p4()); 
      //double scale = corrector->correction(*jet);
      double scale = corrector->correction(*jet,mEvent,mSetup);
      correctedJet.scaleEnergy(scale); 
      mCorrJetPt->Fill(correctedJet.pt());
      mCorrJetPt_80->Fill(correctedJet.pt());  
      if (correctedJet.pt()>10) mCorrJetEta->Fill(correctedJet.eta());
      mCorrJetPhi->Fill(correctedJet.phi());
      mpTRatio->Fill(log10(jet->pt()),correctedJet.pt()/jet->pt());

       if (fabs(jet->eta())<1.3) {
           mpTRatioB_d->Fill(log10(jet->pt()), correctedJet.pt()/jet->pt());
      }

     if (fabs(jet->eta())>1.3 && fabs(jet->eta())<3.0) {
        mpTRatioE_d->Fill (log10(jet->pt()), correctedJet.pt()/jet->pt());
     }
     if (fabs(jet->eta())>3.0 && fabs(jet->eta())<6.0) {
        mpTRatioF_d->Fill (log10(jet->pt()), correctedJet.pt()/jet->pt());
    }
     if (jet->pt()>60.0 && jet->pt()<120.0) {
    mpTRatio_60_120_d->Fill (jet->eta(),correctedJet.pt()/jet->pt());
  }
   if (jet->pt()>200.0 && jet->pt()<300.0) {
    mpTRatio_200_300_d->Fill (jet->eta(),correctedJet.pt()/jet->pt());
  }
   if (jet->pt()>600.0 && jet->pt()<900.0) {
    mpTRatio_600_900_d->Fill (jet->eta(),correctedJet.pt()/jet->pt());
  }
   if (jet->pt()>2700.0 && jet->pt()<3500.0) {
    mpTRatio_2700_3500_d->Fill (jet->eta(),correctedJet.pt()/jet->pt());
  }

  }

  // Count Jets above Pt cut
  for (int istep = 0; istep < 100; ++istep) {
    int     njet = 0;
    float ptStep = (istep * (200./100.));

    for ( JPTJetCollection::const_iterator jpt = jptJets->begin(); jpt != jptJets->end(); ++ jpt ) {
      if ( jpt->pt() > ptStep ) njet++;
    }
    mNJets1->Fill( ptStep, njet );
  }

  for (int istep = 0; istep < 100; ++istep) {
    int     njet = 0;
    float ptStep = (istep * (4000./100.));
    for ( JPTJetCollection::const_iterator jpt = jptJets->begin(); jpt != jptJets->end(); ++ jpt ) {
      if ( jpt->pt() > ptStep ) njet++;
    }
    mNJets2->Fill( ptStep, njet );
  }

if (!mEvent.isRealData()){
  // Gen jet analysis
  Handle<GenJetCollection> genJets;
  mEvent.getByLabel(mInputGenCollection, genJets);
  if (!genJets.isValid()) return;
  GenJetCollection::const_iterator gjet = genJets->begin ();
  int gjetIndex = 0;
  for (; gjet != genJets->end (); gjet++, gjetIndex++) {
    if (mGenEta) mGenEta->Fill (gjet->eta());
    if (mGenPhi) mGenPhi->Fill (gjet->phi());
    if (mGenPt) mGenPt->Fill (gjet->pt());
    if (mGenPt_80) mGenPt_80->Fill (gjet->pt());
    if (gjet == genJets->begin ()) { // first jet
      if (mGenEtaFirst) mGenEtaFirst->Fill (gjet->eta());
      if (mGenPhiFirst) mGenPhiFirst->Fill (gjet->phi());
    }
  }


  // now match JPTJets to GenJets
  JetMatchingTools jetMatching (mEvent);
  if (!(mInputGenCollection.label().empty())) {
    //    Handle<GenJetCollection> genJets;
    //    mEvent.getByLabel(mInputGenCollection, genJets);

    std::vector <std::vector <const reco::GenParticle*> > genJetConstituents (genJets->size());
    std::vector <std::vector <const reco::GenParticle*> > jptJetConstituents (jptJets->size());
    //     if (mRThreshold > 0) { 
    //     }
    //     else {
    //       for (unsigned iGenJet = 0; iGenJet < genJets->size(); ++iGenJet) {
    // 	genJetConstituents [iGenJet] = jetMatching.getGenParticles ((*genJets) [iGenJet]);
    //       }
    
    //       for (unsigned iJPTJet = 0; iJPTJet < jptJets->size(); ++iJPTJet) {
    // 	jptJetConstituents [iJPTJet] = jetMatching.getGenParticles ((*(*jptJets) [iJPTJet].getCaloJetRef()), false);
    //       }
    //     }
    
    for (unsigned iGenJet = 0; iGenJet < genJets->size(); ++iGenJet) {               //****************************************************************
    //for (unsigned iGenJet = 0; iGenJet < 1; ++iGenJet) {                           // only FIRST Jet !!!!
      const GenJet& genJet = (*genJets) [iGenJet];
      double genJetPt = genJet.pt();

      //std::cout << iGenJet <<". Genjet: pT = " << genJetPt << "GeV" << std::endl;  //  *****************************************************

      if (fabs(genJet.eta()) > 6.) continue; // out of detector 
      if (genJetPt < mMatchGenPtThreshold) continue; // no low momentum 
      //double logPtGen = log10 (genJetPt);
      //mAllGenJetsPt->Fill (logPtGen);
      //mAllGenJetsEta->Fill (logPtGen, genJet.eta());
      if (jptJets->size() <= 0) continue; // no JPTJets - nothing to match
      if (mRThreshold > 0) {
	unsigned iJPTJetBest = 0;
	double deltaRBest = 999.;
	for (unsigned iJPTJet = 0; iJPTJet < jptJets->size(); ++iJPTJet) {
	  double dR = deltaR (genJet.eta(), genJet.phi(), (*jptJets) [iJPTJet].eta(), (*jptJets) [iJPTJet].phi());
	  if (deltaRBest < mRThreshold && dR < mRThreshold && genJet.pt() > 5.) {
	    /*
	    std::cout << "Yet another matched jet for GenJet pt=" << genJet.pt()
		      << " previous JPTJet pt/dr: " << (*jptJets) [iJPTJetBest].pt() << '/' << deltaRBest
		      << " new JPTJet pt/dr: " << (*jptJets) [iJPTJet].pt() << '/' << dR
		      << std::endl;
	    */
	  }
	  if (dR < deltaRBest) {
	    iJPTJetBest = iJPTJet;
	    deltaRBest = dR;
	  }
	}
	if (mTurnOnEverything.compare("yes")==0) {
	  //mRMatch->Fill (logPtGen, genJet.eta(), deltaRBest);
	}
	if (deltaRBest < mRThreshold) { // Matched
	  fillMatchHists (genJet, (*jptJets) [iJPTJetBest]);
	}

		///////////pT Response///////////////
	double CorrdeltaRBest = 999.;
	double CorrJetPtBest = 0;
	for (JPTJetCollection::const_iterator jet = jptJets->begin(); jet !=jptJets ->end(); jet++) {
	  JPTJet  correctedJet = *jet;
	  //double scale = corrector->correction(jet->p4());
	  //double scale = corrector->correction(*jet);
          double scale = corrector->correction(*jet,mEvent,mSetup);
	  correctedJet.scaleEnergy(scale);
	  double CorrJetPt = correctedJet.pt();
	  double CorrdR = deltaR (genJet.eta(), genJet.phi(), correctedJet.eta(), correctedJet.phi());
	  if (CorrdR < CorrdeltaRBest) {
	    CorrdeltaRBest = CorrdR;
	    CorrJetPtBest = CorrJetPt;
	  }
	}
	if (deltaRBest < mRThreshold) { // Matched
	  mpTResponse->Fill(log10(genJet.pt()),CorrJetPtBest/genJet.pt());
	  
	  if (fabs(genJet.eta())<1.3) {
	    mpTResponseB_d->Fill(log10(genJet.pt()), CorrJetPtBest/genJet.pt());
	  }	
	  
	  if (fabs(genJet.eta())>1.3 && fabs(genJet.eta())<3.0) {
	    mpTResponseE_d->Fill (log10(genJet.pt()), CorrJetPtBest/genJet.pt());   
	  }
	  
	  if (fabs(genJet.eta())>3.0 && fabs(genJet.eta())<6.0) {
        mpTResponseF_d->Fill (log10(genJet.pt()), CorrJetPtBest/genJet.pt());
	  }
	  if (genJet.pt()>60.0 && genJet.pt()<120.0) {
	    mpTResponse_60_120_d->Fill (genJet.eta(),CorrJetPtBest/genJet.pt());
	  }
	  if (genJet.pt()>200.0 && genJet.pt()<300.0) {
	    mpTResponse_200_300_d->Fill (genJet.eta(),CorrJetPtBest/genJet.pt());
	  }
	  if (genJet.pt()>600.0 && genJet.pt()<900.0) {
	    mpTResponse_600_900_d->Fill (genJet.eta(),CorrJetPtBest/genJet.pt());
	  }
	  if (genJet.pt()>2700.0 && genJet.pt()<3500.0) {
	    mpTResponse_2700_3500_d->Fill (genJet.eta(),CorrJetPtBest/genJet.pt());
	  }
	  if (genJet.pt()>30.0) {
	    mpTResponse_30_d->Fill (genJet.eta(),CorrJetPtBest/genJet.pt());
	  }
	}
	///////////////////////////////////
      }
      //       else {
      // 	unsigned iJPTJetBest = 0;
      // 	double energyFractionBest = 0.;
      // 	for (unsigned iJPTJet = 0; iJPTJet < jptJets->size(); ++iJPTJet) {
      // 	  double energyFraction = jetMatching.overlapEnergyFraction (genJetConstituents [iGenJet], 
      // 								     jptJetConstituents [iJPTJet]);
      // 	  if (energyFraction > energyFractionBest) {
      // 	    iJPTJetBest = iJPTJet;
      // 	    energyFractionBest = energyFraction;
      // 	  }
      // 	}
      // 	if (mTurnOnEverything.compare("yes")==0) {
      // 	  mGenJetMatchEnergyFraction->Fill (logPtGen, genJet.eta(), energyFractionBest);
      // 	}
      // 	if (energyFractionBest > mGenEnergyFractionThreshold) { // good enough
      // 	  double reverseEnergyFraction = jetMatching.overlapEnergyFraction (jptJetConstituents [iJPTJetBest], 
      // 									    genJetConstituents [iGenJet]);
      // 	  if (mTurnOnEverything.compare("yes")==0) {
      // 	    mReverseMatchEnergyFraction->Fill (logPtGen, genJet.eta(), reverseEnergyFraction);
      // 	  }
      // 	  if (reverseEnergyFraction > mReverseEnergyFractionThreshold) { // Matched
      // 	    fillMatchHists (genJet, (*jptJets) [iJPTJetBest]);
      // 	  }
      // 	}
      //       } // else - mRThreshold
    }
  }
}

}////Gen Close

void JPTJetTester::fillMatchHists (const reco::GenJet& fGenJet, const reco::JPTJet& fJPTJet) {
  double logPtGen = log10 (fGenJet.pt());
  double PtGen = fGenJet.pt();
  double PtJpt = fJPTJet.pt();
  //mMatchedGenJetsPt->Fill (logPtGen);
  //mMatchedGenJetsEta->Fill (logPtGen, fGenJet.eta());

  double PtThreshold = 10.;

  if (mTurnOnEverything.compare("yes")==0) {
    mDeltaEta->Fill (logPtGen, fGenJet.eta(), fJPTJet.eta()-fGenJet.eta());
    mDeltaPhi->Fill (logPtGen, fGenJet.eta(), fJPTJet.phi()-fGenJet.phi());
    //mEScale->Fill (logPtGen, fGenJet.eta(), fJPTJet.energy()/fGenJet.energy());
    //mlinEScale->Fill (fGenJet.pt(), fGenJet.eta(), fJPTJet.energy()/fGenJet.energy());
    //mDeltaE->Fill (logPtGen, fGenJet.eta(), fJPTJet.energy()-fGenJet.energy());

    mEScaleFineBin->Fill (logPtGen, fGenJet.eta(), fJPTJet.energy()/fGenJet.energy());
  
    if (fGenJet.pt()>PtThreshold) {
      mEScale_pt10->Fill (logPtGen, fGenJet.eta(), fJPTJet.energy()/fGenJet.energy());

    }

  }
  if (fJPTJet.pt() > PtThreshold) {
    mDelEta->Fill (fGenJet.eta()-fJPTJet.eta());
    mDelPhi->Fill (fGenJet.phi()-fJPTJet.phi());
    mDelPt->Fill  ((fGenJet.pt()-fJPTJet.pt())/fGenJet.pt());
  }

  if (fabs(fGenJet.eta())<1.3) {

    //mpTScaleB_s->Fill (log10(PtGen), PtJpt/PtGen);
    mpTScaleB_d->Fill (log10(PtGen), PtJpt/PtGen);
    mpTScalePhiB_d->Fill (fGenJet.phi(), PtJpt/PtGen);
    
    if (PtGen>60.0 && PtGen<120.0) {
      mpTScale1DB_60_120->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    if (PtGen>200.0 && PtGen<300.0) {
      mpTScale1DB_200_300->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    if (PtGen>600.0 && PtGen<900.0) {
      mpTScale1DB_600_900->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    if (PtGen>2700.0 && PtGen<3500.0) {
      mpTScale1DB_2700_3500->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    
  }

  if (fabs(fGenJet.eta())>1.3 && fabs(fGenJet.eta())<3.0) {
    //mpTScaleE_s->Fill (log10(PtGen), PtJpt/PtGen);
    mpTScaleE_d->Fill (log10(PtGen), PtJpt/PtGen);
    mpTScalePhiE_d->Fill (fGenJet.phi(), PtJpt/PtGen);
    
    if (PtGen>60.0 && PtGen<120.0) {
      mpTScale1DE_60_120->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    if (PtGen>200.0 && PtGen<300.0) {
      mpTScale1DE_200_300->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    if (PtGen>600.0 && PtGen<900.0) {
      mpTScale1DE_600_900->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    if (PtGen>2700.0 && PtGen<3500.0) {
      mpTScale1DE_2700_3500->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    
  }

  if (fabs(fGenJet.eta())>3.0 && fabs(fGenJet.eta())<6.0) {

    //mpTScaleF_s->Fill (log10(PtGen), PtJpt/PtGen);
    mpTScaleF_d->Fill (log10(PtGen), PtJpt/PtGen);
    mpTScalePhiF_d->Fill (fGenJet.phi(), PtJpt/PtGen);
    
    if (PtGen>60.0 && PtGen<120.0) {
      mpTScale1DF_60_120->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    if (PtGen>200.0 && PtGen<300.0) {
      mpTScale1DF_200_300->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    if (PtGen>600.0 && PtGen<900.0) {
      mpTScale1DF_600_900->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    if (PtGen>2700.0 && PtGen<3500.0) {
      mpTScale1DF_2700_3500->Fill (fJPTJet.pt()/fGenJet.pt());
    }
    
  }

  if (fGenJet.pt()>60.0 && fGenJet.pt()<120.0) {
    //mpTScale_60_120_s->Fill (fGenJet.eta(),fJPTJet.pt()/fGenJet.pt());
    mpTScale_60_120_d->Fill (fGenJet.eta(),fJPTJet.pt()/fGenJet.pt());
    //mpTScale1D_60_120->Fill (fJPTJet.pt()/fGenJet.pt());
  }

  if (fGenJet.pt()>200.0 && fGenJet.pt()<300.0) {
    //mpTScale_200_300_s->Fill (fGenJet.eta(),fJPTJet.pt()/fGenJet.pt());
    mpTScale_200_300_d->Fill (fGenJet.eta(),fJPTJet.pt()/fGenJet.pt());
    //mpTScale1D_200_300->Fill (fJPTJet.pt()/fGenJet.pt());
  }

  if (fGenJet.pt()>600.0 && fGenJet.pt()<900.0) {
    //mpTScale_600_900_s->Fill (fGenJet.eta(),fJPTJet.pt()/fGenJet.pt());
    mpTScale_600_900_d->Fill (fGenJet.eta(),fJPTJet.pt()/fGenJet.pt());
    //mpTScale1D_600_900->Fill (fJPTJet.pt()/fGenJet.pt());
  }

  if (fGenJet.pt()>2700.0 && fGenJet.pt()<3500.0) {
    //mpTScale_2700_3500_s->Fill (fGenJet.eta(),fJPTJet.pt()/fGenJet.pt());
    mpTScale_2700_3500_d->Fill (fGenJet.eta(),fJPTJet.pt()/fGenJet.pt());
    //mpTScale1D_2700_3500->Fill (fJPTJet.pt()/fGenJet.pt());
  }



}

double JPTJetTester::getSumPt(const reco::TrackRefVector& tracks){

  double sumpt = 0.;
  
  for (reco::TrackRefVector::const_iterator itrack = tracks.begin(); itrack != tracks.end(); ++itrack){
    const reco::Track& track = **itrack;
    sumpt += track.pt(); 
  }
  
  return sumpt;

}
