#ifndef L1TauEfficiencyAnalyzer_H
#define L1TauEfficiencyAnalyzer_H

/** \class L1TauEfficiencyAnalyzer
 *
 *
 *  This class is an EDFilter for selecting events
 *  passing chosen HLT bits
 *
 *  \author Sami Lehti  -  HIP Helsinki
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"


#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include <TTree.h>
#include <TFile.h>

#include <string>

using namespace edm;
using namespace l1extra;
using namespace reco;
using namespace std;

class L1TauEfficiencyAnalyzer {

    public:
	typedef math::XYZTLorentzVector LorentzVector;

        L1TauEfficiencyAnalyzer();
        ~L1TauEfficiencyAnalyzer();

        void Setup(const edm::ParameterSet&,TTree *l1tree);
        void fill(const edm::Event&, const reco::PFTau&);
        void fill(const edm::Event&, const reco::CaloTau&);
        void fill(const edm::Event&, const LorentzVector&);

   private:
        // Input parameters
	InputTag PFTauCollectionSource;  
	InputTag L1extraTauJetSource;  
        InputTag L1extraCentralJetSource;
        InputTag L1bitInfoSource;
	double   jetMatchingCone;

        // Output tree and related variables
        TTree *l1tree;

        float jetPt, jetEt, jetEta, jetPhi;
        char hasL1Jet, hasL1TauJet, hasL1CenJet;
        char hasTauVeto, hasEmTauVeto, hasHadTauVeto, hasIsolationVeto, hasSumEtBelowThres, hasMaxEt, hasSoft, hasHard;
};
#endif
