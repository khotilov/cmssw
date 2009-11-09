#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "MuonAnalysis/MomentumScaleCalibration/interface/MomentumScaleCorrector.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/ResolutionFunction.h"

//
// class declaration
//
class DistortedMuonProducer : public edm::EDProducer {
   public:
      explicit DistortedMuonProducer(const edm::ParameterSet&);
      ~DistortedMuonProducer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      edm::InputTag muonTag_;
      edm::InputTag genMatchMapTag_;
      std::vector<double> etaBinEdges_;

      bool useDBForMomentumCorrections_;
      std::string dbScaleLabel_;
      std::string dbDataResolutionLabel_;
      std::string dbMCResolutionLabel_;

      std::vector<double> momentumScaleShift_;
      std::vector<double> uncertaintyOnOneOverPt_; // in [1/GeV]
      std::vector<double> relativeUncertaintyOnPt_;

      std::vector<double> efficiencyRatioOverMC_;

      std::auto_ptr<MomentumScaleCorrector> momCorrector_;
      std::auto_ptr<ResolutionFunction> momResolutionData_;
      std::auto_ptr<ResolutionFunction> momResolutionMC_;
};

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Random/RandGauss.h>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CondFormats/RecoMuonObjects/interface/MuScleFitDBobject.h"
#include "CondFormats/DataRecord/interface/MuScleFitDBobjectRcd.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/BaseFunction.h"
#include "MuonAnalysis/MomentumScaleCalibration/interface/Functions.h"

/////////////////////////////////////////////////////////////////////////////////////
DistortedMuonProducer::DistortedMuonProducer(const edm::ParameterSet& pset) {

  // What is being produced
      produces<std::vector<reco::Muon> >();

  // Input products
      muonTag_ = pset.getUntrackedParameter<edm::InputTag> ("MuonTag", edm::InputTag("muons"));
      genMatchMapTag_ = pset.getUntrackedParameter<edm::InputTag> ("GenMatchMapTag", edm::InputTag("genMatchMap"));
      useDBForMomentumCorrections_ = pset.getUntrackedParameter<bool> ("UseDBForMomentumCorrections", false);
      dbScaleLabel_ = pset.getUntrackedParameter<std::string> ("DBScaleLabel", "scale");
      dbDataResolutionLabel_ = pset.getUntrackedParameter<std::string> ("DBDataResolutionLabel", "datareso");
      dbMCResolutionLabel_ = pset.getUntrackedParameter<std::string> ("DBMCResolutionLabel", "mcreso");

  // Eta edges
      std::vector<double> defEtaEdges;
      defEtaEdges.push_back(-999999.);
      defEtaEdges.push_back(999999.);
      etaBinEdges_ = pset.getUntrackedParameter<std::vector<double> > ("EtaBinEdges",defEtaEdges);
      unsigned int ninputs_expected = etaBinEdges_.size()-1;

  // Distortions in muon momentum
      std::vector<double> defDistortion;
      defDistortion.push_back(0.);

      momentumScaleShift_ = pset.getUntrackedParameter<std::vector<double> > ("MomentumScaleShift",defDistortion);
      if (momentumScaleShift_.size()==1 && ninputs_expected>1) {
            for (unsigned int i=1; i<ninputs_expected; i++){ momentumScaleShift_.push_back(momentumScaleShift_[0]);}
      }

      uncertaintyOnOneOverPt_ = pset.getUntrackedParameter<std::vector<double> > ("UncertaintyOnOneOverPt",defDistortion); // in [1/GeV]
      if (uncertaintyOnOneOverPt_.size()==1 && ninputs_expected>1) {
            for (unsigned int i=1; i<ninputs_expected; i++){ uncertaintyOnOneOverPt_.push_back(uncertaintyOnOneOverPt_[0]);}
      }

      relativeUncertaintyOnPt_ = pset.getUntrackedParameter<std::vector<double> > ("RelativeUncertaintyOnPt",defDistortion);
      if (relativeUncertaintyOnPt_.size()==1 && ninputs_expected>1) {
            for (unsigned int i=1; i<ninputs_expected; i++){ relativeUncertaintyOnPt_.push_back(relativeUncertaintyOnPt_[0]);}
      }

  // Data/MC efficiency ratios
      std::vector<double> defEfficiencyRatio;
      defEfficiencyRatio.push_back(1.);
      efficiencyRatioOverMC_ = pset.getUntrackedParameter<std::vector<double> > ("EfficiencyRatioOverMC",defEfficiencyRatio);
      if (efficiencyRatioOverMC_.size()==1 && ninputs_expected>1) {
            for (unsigned int i=1; i<ninputs_expected; i++){ efficiencyRatioOverMC_.push_back(efficiencyRatioOverMC_[0]);}
      }

  // Send a warning if there are inconsistencies in vector sizes !!
      bool effWrong = efficiencyRatioOverMC_.size()!=ninputs_expected;
      bool momWrong =    momentumScaleShift_.size()!=ninputs_expected 
                      || uncertaintyOnOneOverPt_.size()!=ninputs_expected 
                      || relativeUncertaintyOnPt_.size()!=ninputs_expected;
      momWrong = momWrong && (!useDBForMomentumCorrections_);
      if ( effWrong and momWrong) {
           edm::LogError("") << "WARNING: DistortedMuonProducer : Size of some parameters do not match the EtaBinEdges vector!!";
      }

} 

/////////////////////////////////////////////////////////////////////////////////////
DistortedMuonProducer::~DistortedMuonProducer(){
}

/////////////////////////////////////////////////////////////////////////////////////
void DistortedMuonProducer::beginJob(const edm::EventSetup& iSetup) {
      if (useDBForMomentumCorrections_) {
            edm::ESHandle<MuScleFitDBobject> dbObject1;
            iSetup.get<MuScleFitDBobjectRcd>().get(dbScaleLabel_,dbObject1);
            momCorrector_.reset(new MomentumScaleCorrector(dbObject1.product()));

            LogTrace("") << ">>> Using database for momentum scale corrections !!";

            edm::ESHandle<MuScleFitDBobject> dbObject2;
            iSetup.get<MuScleFitDBobjectRcd>().get(dbDataResolutionLabel_, dbObject2);
            momResolutionData_.reset(new ResolutionFunction(dbObject2.product()));

            edm::ESHandle<MuScleFitDBobject> dbObject3;
            iSetup.get<MuScleFitDBobjectRcd>().get(dbMCResolutionLabel_, dbObject3);
            momResolutionMC_.reset(new ResolutionFunction(dbObject3.product()));

            LogTrace("") << ">>> Using database for momentum resolution corrections !!";

      }
      
}

/////////////////////////////////////////////////////////////////////////////////////
void DistortedMuonProducer::endJob(){
}

/////////////////////////////////////////////////////////////////////////////////////
void DistortedMuonProducer::produce(edm::Event& ev, const edm::EventSetup& iSetup) {

      if (ev.isRealData()) return;

      // Muon collection
      edm::Handle<edm::View<reco::Muon> > muonCollection;
      if (!ev.getByLabel(muonTag_, muonCollection)) {
            edm::LogError("") << ">>> Muon collection does not exist !!!";
            return;
      }

      edm::Handle<reco::GenParticleMatch> genMatchMap;
      if (!ev.getByLabel(genMatchMapTag_, genMatchMap)) {
            edm::LogError("") << ">>> Muon-GenParticle match map does not exist !!!";
            return;
      }
  
      unsigned int muonCollectionSize = muonCollection->size();

      std::auto_ptr<reco::MuonCollection> newmuons (new reco::MuonCollection);

      for (unsigned int i=0; i<muonCollectionSize; i++) {
            edm::RefToBase<reco::Muon> mu = muonCollection->refAt(i);

            double ptgen = mu->pt();
            double etagen = mu->eta();
            double phigen = mu->phi();
            int chrgen = mu->charge();
            reco::GenParticleRef gen = (*genMatchMap)[mu];
            if( !gen.isNull()) {
                  ptgen = gen->pt();
                  etagen = gen->eta();
                  phigen = gen->phi();
                  chrgen = gen->charge();
                  LogTrace("") << ">>> Muon-GenParticle match found; ptmu= " << mu->pt() << ", ptgen= " << ptgen;
            } else {
                  LogTrace("") << ">>> MUON-GENPARTICLE MATCH NOT FOUND!!!";
            }

            // Initialize parameters
            double effRatio = 0.;
            double shift = 0.;
            double sigma1 = 0.;
            double sigma2 = 0.;

            // Find out which eta bin should be used
            unsigned int nbins = etaBinEdges_.size()-1;
            unsigned int etaBin = nbins;
            if (etagen>etaBinEdges_[0] && etagen<etaBinEdges_[nbins]) {
                  for (unsigned int j=1; j<=nbins; ++j) {
                        if (etagen>etaBinEdges_[j]) continue;
                        etaBin = j-1;
                        break;
                  }
            }
            if (etaBin<nbins) {
                  LogTrace("") << ">>> etaBin: " << etaBin << ", for etagen =" << etagen;
            } else {
                  // Muon is rejected if outside the considered eta range
                  LogTrace("") << ">>> Muon outside eta range: reject it; etagen = " << etagen;
                  continue;
            }

            // Set shift
            if (useDBForMomentumCorrections_) {
                  shift = ( (*momCorrector_)(*mu) - mu->pt() ) / ptgen;
            } else {
                  shift = momentumScaleShift_[etaBin];
            }
            LogTrace("") << "\tmomentumScaleShift= " << shift*100 << " [%]"; 

            // Set resolutions
            if (useDBForMomentumCorrections_) {
                  sigma1 = 0.;
                  sigma2 = pow(momResolutionData_->sigmaPt(*mu),2) -
                              pow(momResolutionMC_->sigmaPt(*mu),2);
                  if (sigma2>0.) sigma2 = sqrt(sigma2)/ptgen; else sigma2 = 0.;
            } else {
                  sigma1 = uncertaintyOnOneOverPt_[etaBin];
                  sigma2 = relativeUncertaintyOnPt_[etaBin];
            }

            LogTrace("") << "\tuncertaintyOnOneOverPt= " << sigma1 << " [1/GeV]"; 
            LogTrace("") << "\trelativeUncertaintyOnPt= " << sigma2*100 << " [%]"; 

            // Set efficiency ratio
            effRatio = efficiencyRatioOverMC_[etaBin];
            LogTrace("") << "\tefficiencyRatioOverMC= " << effRatio;

            // Reject muons according to efficiency ratio
            double rndf = CLHEP::RandFlat::shoot();
            if (rndf>effRatio) continue;

            // Gaussian Random numbers for smearing
            double rndg1 = CLHEP::RandGauss::shoot();
            double rndg2 = CLHEP::RandGauss::shoot();
            
            // New muon
            double ptmu = mu->pt();
            ptmu += ptgen * ( shift + sigma1*rndg1*ptgen + sigma2*rndg2);
            reco::Muon* newmu = mu->clone();
            newmu->setP4 (
                  reco::Particle::PolarLorentzVector (
                        ptmu, mu->eta(), mu->phi(), mu->mass()
                  )
            );
            newmuons->push_back(*newmu);

      }

      ev.put(newmuons);
}

DEFINE_FWK_MODULE(DistortedMuonProducer);
