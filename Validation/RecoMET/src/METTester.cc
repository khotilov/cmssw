#include "Validation/RecoMET/interface/METTester.h"
// author: Mike Schmitt, University of Florida
// first version 8/24/2006
// modification: Bobby Scurlock
// date:  03.11.2006
// note:  added RMS(METx) vs SumET capability
// modification: Rick Cavanaugh
// date:  05.11.2006
// note:  cleaned up constructor and beginJob, removed int conv. warning
//        added configuration params
// modification: Mike Schmitt
// date:  02.28.2007
// note:  code rewrite. Now uses STL map for monitoring element container. 
// modification: Bobby Scurlock
// date:  04.03.2007
// note:  Eliminated automated resolution fitting. This is now done in a ROOT script.

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
//#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"

#include <vector>
#include <utility>
#include <ostream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <memory>
#include "DQMServices/Core/interface/DQMStore.h"

METTester::METTester(const edm::ParameterSet& iConfig)
{

  outputFile_              = iConfig.getUntrackedParameter<std::string>("OutputFile");
  inputMETLabel_           = iConfig.getUntrackedParameter<edm::InputTag>("InputMETLabel");
  METType_                 = iConfig.getUntrackedParameter<std::string>("METType");
  
  if (outputFile_.size() > 0)
    edm::LogInfo("OutputInfo") << " MET Task histograms will be saved to '" << outputFile_.c_str() << "'";
  else edm::LogInfo("OutputInfo") << " MET Task histograms will NOT be saved";
  
}

void METTester::beginJob(const edm::EventSetup& iSetup)
{
  
  // get ahold of back-end interface
  dbe_ = edm::Service<DQMStore>().operator->();
  
  if (dbe_) {
    TString dirName = "RecoMETV/METTask/";
    TString label(inputMETLabel_.label());
    dirName += label;
    dbe_->setCurrentFolder((string)dirName);
    
    if (METType_ == "CaloMET")
      { 
	// CaloMET Histograms
	me["hNevents"]                = dbe_->book1D("METTask_Nevents","METTask_Nevents",1,0,1);
	me["hCaloMEx"]                = dbe_->book1D("METTask_CaloMEx","METTask_CaloMEx",2001,-500,501);
	me["hCaloMEy"]                = dbe_->book1D("METTask_CaloMEy","METTask_CaloMEy",2001,-500,501);
	me["hCaloEz"]                 = dbe_->book1D("METTask_CaloEz","METTask_CaloEz",2001,-500,501);
	me["hCaloMETSig"]             = dbe_->book1D("METTask_CaloMETSig","METTask_CaloMETSig",51,0,51);
	me["hCaloMET"]                = dbe_->book1D("METTask_CaloMET","METTask_CaloMET",2001,0,2001);
	me["hCaloMETPhi"]             = dbe_->book1D("METTask_CaloMETPhi","METTask_CaloMETPhi",80,-4,4);
	me["hCaloSumET"]              = dbe_->book1D("METTask_CaloSumET","METTask_CaloSumET",4001,0,4001);
	me["hCaloMaxEtInEmTowers"]    = dbe_->book1D("METTask_CaloMaxEtInEmTowers","METTask_CaloMaxEtInEmTowers",4001,0,4001);
	me["hCaloMaxEtInHadTowers"]   = dbe_->book1D("METTask_CaloMaxEtInHadTowers","METTask_CaloMaxEtInHadTowers",4001,0,4001);
	me["hCaloEtFractionHadronic"] = dbe_->book1D("METTask_CaloEtFractionHadronic","METTask_CaloEtFractionHadronic",100,0,1);
	me["hCaloEmEtFraction"]       = dbe_->book1D("METTask_CaloEmEtFraction","METTask_CaloEmEtFraction",100,0,1);
	me["hCaloHadEtInHB"]          = dbe_->book1D("METTask_CaloHadEtInHB","METTask_CaloHadEtInHB",4001,0,4001);
	me["hCaloHadEtInHO"]          = dbe_->book1D("METTask_CaloHadEtInHO","METTask_CaloHadEtInHO",4001,0,4001);
	me["hCaloHadEtInHE"]          = dbe_->book1D("METTask_CaloHadEtInHE","METTask_CaloHadEtInHE",4001,0,4001);
	me["hCaloHadEtInHF"]          = dbe_->book1D("METTask_CaloHadEtInHF","METTask_CaloHadEtInHF",4001,0,4001);
	me["hCaloHadEtInEB"]          = dbe_->book1D("METTask_CaloHadEtInEB","METTask_CaloHadEtInEB",4001,0,4001);
	me["hCaloHadEtInEE"]          = dbe_->book1D("METTask_CaloHadEtInEE","METTask_CaloHadEtInEE",4001,0,4001);
	me["hCaloEmEtInHF"]           = dbe_->book1D("METTask_CaloEmEtInHF","METTask_CaloEmEtInHF",4001,0,4001);
	me["hCaloEmEtInEE"]           = dbe_->book1D("METTask_CaloEmEtInEE","METTask_CaloEmEtInEE",4001,0,4001);
	me["hCaloEmEtInEB"]           = dbe_->book1D("METTask_CaloEmEtInEB","METTask_CaloEmEtInEB",4001,0,4001);
      }

    else if (METType_ == "GenMET")
      {
	// GenMET Histograms
	me["hNevents"]                = dbe_->book1D("METTask_Nevents","METTask_Nevents",1,0,1);
	me["hGenMEx"]                 = dbe_->book1D("METTask_GenMEx","METTask_GenMEx",2001,-500,501);
	me["hGenMEy"]                 = dbe_->book1D("METTask_GenMEy","METTask_GenMEy",2001,-500,501);
	me["hGenEz"]                  = dbe_->book1D("METTask_GenEz","METTask_GenEz",2001,-500,501);
	me["hGenMETSig"]              = dbe_->book1D("METTask_GenMETSig","METTask_GenMETSig",51,0,51);
	me["hGenMET"]                 = dbe_->book1D("METTask_GenMET","METTask_GenMET",2001,0,2001);
	me["hGenMETPhi"]              = dbe_->book1D("METTask_GenMETPhi","METTask_GenMETPhi",80,-4,4);
	me["hGenSumET"]               = dbe_->book1D("METTask_GenSumET","METTask_GenSumET",4001,0,4001);
	me["hGenEmEnergy"]            = dbe_->book1D("METTask_GenEmEnergy","METTask_GenEmEnergy",4001,0,4001);
	me["hGenHadEnergy"]           = dbe_->book1D("METTask_GenHadEnergy","METTask_GenHadEnergy",4001,0,4001);
	me["hGenInvisibleEnergy"]     = dbe_->book1D("METTask_GenInvisibleEnergy","METTask_GenInvisibleEnergy",4001,0,4001);
	me["hGenAuxiliaryEnergy"]     = dbe_->book1D("METTask_GenAuxiliaryEnergy","METTask_GenAuxiliaryEnergy",4001,0,4001);
      }
    else if (METType_ == "MET")
      {
	// MET Histograms
	me["hNevents"]                = dbe_->book1D("METTask_Nevents","METTask_Nevents",1,0,1);
	me["hMEx"]                = dbe_->book1D("METTask_MEx","METTask_MEx",2001,-500,501);
	me["hMEy"]                = dbe_->book1D("METTask_MEy","METTask_MEy",2001,-500,501);
	me["hEz"]                 = dbe_->book1D("METTask_Ez","METTask_Ez",2001,-500,501);
	me["hMETSig"]             = dbe_->book1D("METTask_METSig","METTask_METSig",51,0,51);
	me["hMET"]                = dbe_->book1D("METTask_MET","METTask_MET",2001,0,2001);
	me["hMETPhi"]             = dbe_->book1D("METTask_METPhi","METTask_METPhi",80,-4,4);
	me["hSumET"]              = dbe_->book1D("METTask_SumET","METTask_SumET",4001,0,4001);
      }
    
    else
      {
	edm::LogInfo("OutputInfo") << " METType not correctly specified!'" << outputFile_.c_str();
      }
  }
}

void METTester::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  if (METType_ == "CaloMET")
    { 
      const CaloMET *calomet;
      // Get CaloMET
      edm::Handle<CaloMETCollection> calo;
      iEvent.getByLabel(inputMETLabel_, calo);
      if (!calo.isValid()) {
	edm::LogInfo("OutputInfo") << " failed to retrieve data required by MET Task";
	edm::LogInfo("OutputInfo") << " MET Task cannot continue...!";
	return;
      } else {
	const CaloMETCollection *calometcol = calo.product();
	calomet = &(calometcol->front());
      }
      // ==========================================================
      // Reconstructed MET Information
      double caloSumET = calomet->sumEt();
      double caloMETSig = calomet->mEtSig();
      double caloEz = calomet->e_longitudinal();
      double caloMET = calomet->pt();
      double caloMEx = calomet->px();
      double caloMEy = calomet->py();
      double caloMETPhi = calomet->phi();
      double caloMaxEtInEMTowers = calomet->maxEtInEmTowers();
      double caloMaxEtInHadTowers = calomet->maxEtInHadTowers();
      double caloEtFractionHadronic = calomet->etFractionHadronic();
      double caloEmEtFraction = calomet->emEtFraction();
      double caloHadEtInHB = calomet->hadEtInHB();
      double caloHadEtInHO = calomet->hadEtInHO();
      double caloHadEtInHE = calomet->hadEtInHE();
      double caloHadEtInHF = calomet->hadEtInHF();
      double caloEmEtInEB = calomet->emEtInEB();
      double caloEmEtInEE = calomet->emEtInEE();
      double caloEmEtInHF = calomet->emEtInHF();

      me["hCaloMEx"]->Fill(caloMEx);
      me["hCaloMEy"]->Fill(caloMEy);
      me["hCaloMET"]->Fill(caloMET);
      me["hCaloMETPhi"]->Fill(caloMETPhi);
      me["hCaloSumET"]->Fill(caloSumET);
      me["hCaloMETSig"]->Fill(caloMETSig);
      me["hCaloEz"]->Fill(caloEz);
      me["hCaloMaxEtInEmTowers"]->Fill(caloMaxEtInEMTowers);
      me["hCaloMaxEtInHadTowers"]->Fill(caloMaxEtInHadTowers);
      me["hCaloEtFractionHadronic"]->Fill(caloEtFractionHadronic);
      me["hCaloEmEtFraction"]->Fill(caloEmEtFraction);
      me["hCaloHadEtInHB"]->Fill(caloHadEtInHB);
      me["hCaloHadEtInHO"]->Fill(caloHadEtInHO);
      me["hCaloHadEtInHE"]->Fill(caloHadEtInHE);
      me["hCaloHadEtInHF"]->Fill(caloHadEtInHF);
      me["hCaloEmEtInEB"]->Fill(caloEmEtInEB);
      me["hCaloEmEtInEE"]->Fill(caloEmEtInEE);
      me["hCaloEmEtInHF"]->Fill(caloEmEtInHF);
    }
  
  else if (METType_ == "GenMET")
    {
      const GenMET *genmet;
      // Get Generated MET
      edm::Handle<GenMETCollection> gen;
      iEvent.getByLabel(inputMETLabel_, gen);
      if (!gen.isValid()) {
	edm::LogInfo("OutputInfo") << " failed to retrieve data required by MET Task";
	edm::LogInfo("OutputInfo") << " MET Task cannot continue...!";
	return;
      } else {
	const GenMETCollection *genmetcol = gen.product();
	genmet = &(genmetcol->front());
      }    
      
      // ==========================================================
      // Genenerated MET Information  
      double genSumET = genmet->sumEt();
      double genMETSig = genmet->mEtSig();
      double genEz = genmet->e_longitudinal();
      double genMET = genmet->pt();
      double genMEx = genmet->px();
      double genMEy = genmet->py();
      double genMETPhi = genmet->phi();
      double genEmEnergy = genmet->emEnergy();
      double genHadEnergy = genmet->hadEnergy();
      double genInvisibleEnergy= genmet->invisibleEnergy();
      double genAuxiliaryEnergy= genmet->auxiliaryEnergy();
      
      me["hNevents"]->Fill(0);
      me["hGenMEx"]->Fill(genMEx);
      me["hGenMEy"]->Fill(genMEy);
      me["hGenMET"]->Fill(genMET);
      me["hGenMETPhi"]->Fill(genMETPhi);
      me["hGenSumET"]->Fill(genSumET);
      me["hGenMETSig"]->Fill(genMETSig);
      me["hGenEz"]->Fill(genEz);
      me["hGenEmEnergy"]->Fill(genEmEnergy);
      me["hGenHadEnergy"]->Fill(genHadEnergy);
      me["hGenInvisibleEnergy"]->Fill(genInvisibleEnergy);
      me["hGenAuxiliaryEnergy"]->Fill(genAuxiliaryEnergy);
    }

  else if (METType_ == "MET")
    {
      const MET *met;
      // Get Generated MET
      edm::Handle<METCollection> hmetcol;
      iEvent.getByLabel(inputMETLabel_, hmetcol);
      if (!hmetcol.isValid()) {
	edm::LogInfo("OutputInfo") << " failed to retrieve data required by MET Task";
	edm::LogInfo("OutputInfo") << " MET Task cannot continue...!";
	return;
      } else {
	const METCollection *metcol = hmetcol.product();
	met = &(metcol->front());
      }   

      // Reconstructed MET Information
      double SumET = met->sumEt();
      double METSig = met->mEtSig();
      double Ez = met->e_longitudinal();
      double MET = met->pt();
      double MEx = met->px();
      double MEy = met->py();
      double METPhi = met->phi();

      me["hMEx"]->Fill(MEx);
      me["hMEy"]->Fill(MEy);
      me["hMET"]->Fill(MET);
      me["hMETPhi"]->Fill(METPhi);
      me["hSumET"]->Fill(SumET);
      me["hMETSig"]->Fill(METSig);
      me["hEz"]->Fill(Ez);
    }


}

void METTester::endJob() 
{

  // Store the DAQ Histograms
  if (outputFile_.size() > 0 && dbe_)
    dbe_->save(outputFile_);

}
