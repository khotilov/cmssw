#ifndef TauAnalysis_Core_GenericEventDump_h
#define TauAnalysis_Core_GenericEventDump_h

/** \class GenericEventDump
 *
 * Base-class for print-out of event level information
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.10 $
 *
 * $Id: GenericEventDump.h,v 1.10 2009/12/17 11:03:33 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include "TauAnalysis/Core/interface/EventDumpBase.h"

#include <TMath.h>

#include <vector>
#include <string>

class GenericEventDump : public EventDumpBase
{
 public:
  // constructor 
  explicit GenericEventDump(const edm::ParameterSet&);
  
  // destructor
  virtual ~GenericEventDump();

  // derrived-class method for print-out of event level information
  virtual void analyze(const edm::Event&, const edm::EventSetup&, 
		       const EventDumpBase::filterResults_type&, const EventDumpBase::filterResults_type&, double);

 protected:
//--- function to count types of particles faking reconstructed electrons,
//    muons and tau-jets
  void countFakeParticles(const edm::Event&);

//--- print functions to be used by derrived classes
  virtual void printEventHeaderInfo(const edm::Event&, double) const;
  virtual void printEventTriggerInfo(const edm::Event&) const;

  virtual void printElectronInfo(const edm::Event&) const;
  virtual void printMuonInfo(const edm::Event&) const;
  virtual void printTauInfo(const edm::Event&) const;

  virtual void printDiTauCandidateInfo(const edm::Event&) const = 0;
  template<typename T1, typename T2>
  void printDiTauCandidateInfoImp(const edm::Event& evt) const
  { 
    if ( !outputStream_ ) {
      edm::LogError ("printDiTauCandidateInfoImp") << " outputStream = NULL --> skipping !!";
      return;
    }
    
    if ( diTauCandidateSource_.label() != "" ) {
      typedef std::vector<CompositePtrCandidateT1T2MEt<T1,T2> > CompositePtrCandidateCollection;
      edm::Handle<CompositePtrCandidateCollection> diTauCandidates;
      evt.getByLabel(diTauCandidateSource_, diTauCandidates);
      
      unsigned iDiTauCandidate = 0;
      for ( typename CompositePtrCandidateCollection::const_iterator diTauCandidate = diTauCandidates->begin(); 
	    diTauCandidate != diTauCandidates->end(); ++diTauCandidate ) {
	*outputStream_ << "DiTauCandidate(" << iDiTauCandidate << "):" << std::endl;
	*outputStream_ << " Pt = " << diTauCandidate->pt() << std::endl;
	*outputStream_ << " theta = " << diTauCandidate->theta()*180./TMath::Pi() 
		       << " (eta = " << diTauCandidate->eta() << ")" << std::endl;
	*outputStream_ << " phi = " << diTauCandidate->phi()*180./TMath::Pi() << std::endl;
	*outputStream_ << " Leg1" << std::endl;
	*outputStream_ << "  Pt = " << diTauCandidate->leg1()->pt() << std::endl;
	*outputStream_ << "  theta = " << diTauCandidate->leg1()->theta()*180./TMath::Pi() 
		       << " (eta = " << diTauCandidate->leg1()->eta() << ")" << std::endl;
	*outputStream_ << "  phi = " << diTauCandidate->leg1()->phi()*180./TMath::Pi() << std::endl;
	*outputStream_ << "  pdgId = " << diTauCandidate->leg1()->pdgId() << std::endl;
	*outputStream_ << " Leg2" << std::endl;
	*outputStream_ << "  Pt = " << diTauCandidate->leg2()->pt() << std::endl;
	*outputStream_ << "  theta = " << diTauCandidate->leg2()->theta()*180./TMath::Pi() 
		       << " (eta = " << diTauCandidate->leg2()->eta() << ")" << std::endl;
	*outputStream_ << "  phi = " << diTauCandidate->leg2()->phi()*180./TMath::Pi() << std::endl;
	*outputStream_ << "  pdgId = " << diTauCandidate->leg2()->pdgId() << std::endl;
	*outputStream_ << " dPhi(Leg1,Leg2) = " << diTauCandidate->dPhi12()*180./TMath::Pi() << std::endl;
	*outputStream_ << " M(visible) = " << diTauCandidate->p4Vis().mass() << std::endl;
	*outputStream_ << " Mt(Leg1+MET) = " << diTauCandidate->mt1MET() << std::endl;
	*outputStream_ << " Mt(Leg2+MET) = " << diTauCandidate->mt2MET() << std::endl;
	*outputStream_ << " M(CDF method) = " << diTauCandidate->p4CDFmethod().mass() << std::endl;
	*outputStream_ << " M(collinear Approx.) = " << diTauCandidate->p4CollinearApprox().mass() << std::endl;
	*outputStream_ << "  x1 = " << diTauCandidate->x1CollinearApprox() 
		       << " (gen. = " << diTauCandidate->x1gen() << ")" << std::endl;
	*outputStream_ << "  x2 = " << diTauCandidate->x2CollinearApprox() 
		       << " (gen. = " << diTauCandidate->x2gen() << ")" << std::endl;
	std::string collinearApproxStatus = ( diTauCandidate->collinearApproxIsValid() ) ? "valid" : "invalid";
	*outputStream_ << " (collinear Approx. " << collinearApproxStatus << ")" << std::endl;
	++iDiTauCandidate;
      }

      *outputStream_ << std::endl;
    }
  }

  virtual void printMissingEtInfo(const edm::Event&) const;

  virtual void printJetInfo(const edm::Event&) const;

//--- configuration parameters
  edm::InputTag l1GtReadoutRecordSource_;
  edm::InputTag l1GtObjectMapRecordSource_;
  edm::InputTag hltResultsSource_;

  typedef std::vector<std::string> vstring;
  vstring l1BitsToPrint_;
  vstring hltPathsToPrint_;

  edm::InputTag genParticleSource_;
  edm::InputTag genJetSource_;
  edm::InputTag genTauJetSource_;
  edm::InputTag genEventInfoSource_;

  edm::InputTag patElectronSource_;
  edm::InputTag patMuonSource_;
  edm::InputTag patTauSource_;
  edm::InputTag patJetSource_;

  edm::InputTag diTauCandidateSource_;

  edm::InputTag patCaloMEtSource_;
  edm::InputTag patPFMEtSource_;
  edm::InputTag genMEtSource_;

  std::vector<int> skipPdgIdsGenParticleMatch_;

  edm::InputTag recoTrackSource_;
  edm::InputTag recoVertexSource_;

  edm::InputTag pfChargedHadronSource_;
  edm::InputTag pfGammaSource_;
  edm::InputTag pfNeutralHadronSource_;
  edm::InputTag pfCandidateSource_;

//--- count different types of particles faking reconstructed electrons,
//    muons and tau-jets
  unsigned numRecoElectronsMatchingGenMuons_;
  unsigned numRecoElectronsMatchingGenElectrons_;
  unsigned numRecoElectronsMatchingGenTauJets_;
  unsigned numRecoElectronsMatchingGenBottomQuarks_;
  unsigned numRecoElectronsMatchingGenCharmQuarks_;
  unsigned numRecoElectronsMatchingGenGluons_;
  unsigned numRecoElectronsMatchingGenLightQuarks_;
  unsigned numRecoElectronsUndeterminedGenMatch_;

  unsigned numRecoMuonsMatchingGenMuons_;
  unsigned numRecoMuonsMatchingGenElectrons_;
  unsigned numRecoMuonsMatchingGenTauJets_;
  unsigned numRecoMuonsMatchingGenBottomQuarks_;
  unsigned numRecoMuonsMatchingGenCharmQuarks_;
  unsigned numRecoMuonsMatchingGenGluons_;
  unsigned numRecoMuonsMatchingGenLightQuarks_;
  unsigned numRecoMuonsUndeterminedGenMatch_;

  unsigned numRecoTauJetsMatchingGenMuons_;
  unsigned numRecoTauJetsMatchingGenElectrons_;
  unsigned numRecoTauJetsMatchingGenTauJets_;
  unsigned numRecoTauJetsMatchingGenBottomQuarks_;
  unsigned numRecoTauJetsMatchingGenCharmQuarks_;
  unsigned numRecoTauJetsMatchingGenGluons_;
  unsigned numRecoTauJetsMatchingGenLightQuarks_;
  unsigned numRecoTauJetsUndeterminedGenMatch_;
};

#endif       

