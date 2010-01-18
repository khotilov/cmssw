#include "DQM/Physics/src/EwkElecDQM.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" // I guess this is the right one??
// also need Fwd.h file ???
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "DataFormats/GeometryVector/interface/Phi.h"

#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/Common/interface/View.h"
  
using namespace edm;
using namespace std;
using namespace reco;

EwkElecDQM::EwkElecDQM( const ParameterSet & cfg ) :
      // Input collections
      trigTag_(cfg.getUntrackedParameter<edm::InputTag> ("TrigTag", edm::InputTag("TriggerResults::HLT"))),
      //      muonTag_(cfg.getUntrackedParameter<edm::InputTag> ("MuonTag", edm::InputTag("muons"))),
      elecTag_(cfg.getUntrackedParameter<edm::InputTag> ("ElecTag", edm::InputTag("gsfElectrons"))), 
      //      metTag_(cfg.getUntrackedParameter<edm::InputTag> ("METTag", edm::InputTag("met"))),
      //      metIncludesMuons_(cfg.getUntrackedParameter<bool> ("METIncludesMuons", false)),
      //      jetTag_(cfg.getUntrackedParameter<edm::InputTag> ("JetTag", edm::InputTag("sisCone5CaloJets"))),

      // Main cuts 
      //      muonTrig_(cfg.getUntrackedParameter<std::string> ("MuonTrig", "HLT_Mu9")),
      elecTrig_(cfg.getUntrackedParameter<std::string> ("ElecTrig", "HLT_Ele10_SW_L1R")), // NOT HLT_Ele10_LW_L1R ???
      //      ptCut_(cfg.getUntrackedParameter<double>("PtCut", 25.)),
      ptCut_(cfg.getUntrackedParameter<double>("PtCut", 10.)),
      //      etaCut_(cfg.getUntrackedParameter<double>("EtaCut", 2.1)),
      etaCut_(cfg.getUntrackedParameter<double>("EtaCut", 2.4)),
      sieieCutBarrel_(cfg.getUntrackedParameter<double>("SieieBarrel", 0.01)),
      sieieCutEndcap_(cfg.getUntrackedParameter<double>("SieieEndcap", 0.028)),
      detainCutBarrel_(cfg.getUntrackedParameter<double>("DetainBarrel", 0.0071)),
      detainCutEndcap_(cfg.getUntrackedParameter<double>("DetainEndcap", 0.0066)),
      //      isRelativeIso_(cfg.getUntrackedParameter<bool>("IsRelativeIso", true)),
      //      isCombinedIso_(cfg.getUntrackedParameter<bool>("IsCombinedIso", false)),
      //      isoCut03_(cfg.getUntrackedParameter<double>("IsoCut03", 0.1)),
      ecalIsoCutBarrel_(cfg.getUntrackedParameter<double>("EcalIsoCutBarrel", 5.7)),
      ecalIsoCutEndcap_(cfg.getUntrackedParameter<double>("EcalIsoCutEndcap", 5.0)),
      hcalIsoCutBarrel_(cfg.getUntrackedParameter<double>("HcalIsoCutBarrel", 8.1)),
      hcalIsoCutEndcap_(cfg.getUntrackedParameter<double>("HcalIsoCutEndcap", 3.4)),
      trkIsoCutBarrel_(cfg.getUntrackedParameter<double>("TrkIsoCutBarrel", 7.2)),
      trkIsoCutEndcap_(cfg.getUntrackedParameter<double>("TrkIsoCutEndcap", 5.1))//,
      //      mtMin_(cfg.getUntrackedParameter<double>("MtMin", 50.)),
      //      mtMax_(cfg.getUntrackedParameter<double>("MtMax", 200.)),
      //      metMin_(cfg.getUntrackedParameter<double>("MetMin", -999999.)),
      //      metMax_(cfg.getUntrackedParameter<double>("MetMax", 999999.)),
      //      acopCut_(cfg.getUntrackedParameter<double>("AcopCut", 2.)),

      // Muon quality cuts
      //      dxyCut_(cfg.getUntrackedParameter<double>("DxyCut", 0.2)),
      //      normalizedChi2Cut_(cfg.getUntrackedParameter<double>("NormalizedChi2Cut", 10.)),
      //      trackerHitsCut_(cfg.getUntrackedParameter<int>("TrackerHitsCut", 11)),
      //      isAlsoTrackerMuon_(cfg.getUntrackedParameter<bool>("IsAlsoTrackerMuon", true)),

      // Z rejection
      //      ptThrForZ1_(cfg.getUntrackedParameter<double>("PtThrForZ1", 20.)),
      //      ptThrForZ2_(cfg.getUntrackedParameter<double>("PtThrForZ2", 10.)),

      // Top rejection
      //      eJetMin_(cfg.getUntrackedParameter<double>("EJetMin", 999999.)),
      //      nJetMax_(cfg.getUntrackedParameter<int>("NJetMax", 999999))
{
}

void EwkElecDQM::beginRun(const Run& r, const EventSetup&) {
      nall = 0;
      nsel = 0;

      nrec = 0; 
      neid = 0;
      niso = 0; 
//       nhlt = 0; 
//       nmet = 0;
}


void EwkElecDQM::beginJob(const EventSetup &) {
      theDbe = Service<DQMStore>().operator->();
      theDbe->setCurrentFolder("Physics/EwkElecDQM");

      init_histograms();
}

void EwkElecDQM::init_histograms() {

      char chtitle[256] = "";
      for (int i=0; i<2; ++i) {
//             snprintf(chtitle, 255, "Muon transverse momentum (global muon) [GeV]");
//             pt_before_ = theDbe->book1D("PT_BEFORECUTS",chtitle,100,0.,100.);
//             pt_after_ = theDbe->book1D("PT_LASTCUT",chtitle,100,0.,100.);

            snprintf(chtitle, 255, "Electron transverse momentum [GeV]");
            pt_before_ = theDbe->book1D("PT_BEFORECUTS",chtitle,100,0.,100.);
            pt_after_ = theDbe->book1D("PT_LASTCUT",chtitle,100,0.,100.);

            snprintf(chtitle, 255, "Electron pseudo-rapidity");
            eta_before_ = theDbe->book1D("ETA_BEFORECUTS",chtitle,50,-2.5,2.5);
            eta_after_ = theDbe->book1D("ETA_LASTCUT",chtitle,50,-2.5,2.5);

            snprintf(chtitle, 255, "Electron #sigma_{i#etai#eta} (barrel)");
            sieiebarrel_before_ = theDbe->book1D("SIEIEBARREL_BEFORECUTS",chtitle,70,0.,0.07);
            sieiebarrel_after_ = theDbe->book1D("SIEIEBARREL_LASTCUT",chtitle,70,0.,0.07);

            snprintf(chtitle, 255, "Electron #sigma_{i#etai#eta} (endcap)");
            sieieendcap_before_ = theDbe->book1D("SIEIEENDCAP_BEFORECUTS",chtitle,70,0.,0.07);
            sieieendcap_after_ = theDbe->book1D("SIEIEENDCAP_LASTCUT",chtitle,70,0.,0.07);

            snprintf(chtitle, 255, "Electron #Delta#eta_{in} (barrel)");
            detainbarrel_before_ = theDbe->book1D("DETAINBARREL_BEFORECUTS",chtitle,40,-0.02,0.02);
            detainbarrel_after_ = theDbe->book1D("DETAINBARREL_LASTCUT",chtitle,40,-0.02,0.02);

            snprintf(chtitle, 255, "Electron #Delta#eta_{in} (endcap)");
            detainendcap_before_ = theDbe->book1D("DETAINENDCAP_BEFORECUTS",chtitle,40,-0.02,0.02);
            detainendcap_after_ = theDbe->book1D("DETAINENDCAP_LASTCUT",chtitle,40,-0.02,0.02);

//             snprintf(chtitle, 255, "Muon transverse distance to beam spot [cm]");
//             dxy_before_ = theDbe->book1D("DXY_BEFORECUTS",chtitle,100,-0.5,0.5);
//             dxy_after_ = theDbe->book1D("DXY_LASTCUT",chtitle,100,-0.5,0.5);

//             snprintf(chtitle, 255, "Normalized Chi2, inner track fit");
//             chi2_before_ = theDbe->book1D("CHI2_BEFORECUTS",chtitle,100,0.,100.);
//             chi2_after_ = theDbe->book1D("CHI2_LASTCUT",chtitle,100,0.,100.);

//             snprintf(chtitle, 255, "Number of hits, inner track");
//             nhits_before_ = theDbe->book1D("NHITS_BEFORECUTS",chtitle,40,-0.5,39.5);
//             nhits_after_ = theDbe->book1D("NHITS_LASTCUT",chtitle,40,-0.5,39.5);

//             snprintf(chtitle, 255, "Tracker-muon flag (for global muons)");
//             tkmu_before_ = theDbe->book1D("TKMU_BEFORECUTS",chtitle,2,-0.5,1.5);
//             tkmu_after_ = theDbe->book1D("TKMU_LASTCUT",chtitle,2,-0.5,1.5);

//             if (isRelativeIso_) {
//                   if (isCombinedIso_) {
//                         snprintf(chtitle, 255, "Relative (combined) isolation variable");
//                   } else {
//                         snprintf(chtitle, 255, "Relative (tracker) isolation variable");
//                   }
//                   iso_before_ = theDbe->book1D("ISO_BEFORECUTS",chtitle,100, 0., 1.);
//                   iso_after_ = theDbe->book1D("ISO_LASTCUT",chtitle,100, 0., 1.);
//             } else {
//                   if (isCombinedIso_) {
//                         snprintf(chtitle, 255, "Absolute (combined) isolation variable [GeV]");
//                   } else {
//                         snprintf(chtitle, 255, "Absolute (tracker) isolation variable [GeV]");
//                   }
//                   iso_before_ = theDbe->book1D("ISO_BEFORECUTS",chtitle,100, 0., 20.);
//                   iso_after_ = theDbe->book1D("ISO_LASTCUT",chtitle,100, 0., 20.);
//             }

            snprintf(chtitle, 255, "Absolute electron ECAL isolation variable (barrel) [GeV]");
            ecalisobarrel_before_ = theDbe->book1D("ECALISOBARREL_BEFORECUTS",chtitle,50,0.,50.);
            ecalisobarrel_after_ = theDbe->book1D("ECALISOBARREL_LASTCUT",chtitle,50,0.,50.);

            snprintf(chtitle, 255, "Absolute electron ECAL isolation variable (endcap) [GeV]");
            ecalisoendcap_before_ = theDbe->book1D("ECALISOENDCAP_BEFORECUTS",chtitle,50,0.,50.);
            ecalisoendcap_after_ = theDbe->book1D("ECALISOENDCAP_LASTCUT",chtitle,50,0.,50.);

            snprintf(chtitle, 255, "Absolute electron HCAL isolation variable (barrel) [GeV]");
            hcalisobarrel_before_ = theDbe->book1D("HCALISOBARREL_BEFORECUTS",chtitle,50,0.,50.);
            hcalisobarrel_after_ = theDbe->book1D("HCALISOBARREL_LASTCUT",chtitle,50,0.,50.);

            snprintf(chtitle, 255, "Absolute electron HCAL isolation variable (endcap) [GeV]");
            hcalisoendcap_before_ = theDbe->book1D("HCALISOENDCAP_BEFORECUTS",chtitle,50,0.,50.);
            hcalisoendcap_after_ = theDbe->book1D("HCALISOENDCAP_LASTCUT",chtitle,50,0.,50.);

            snprintf(chtitle, 255, "Absolute electron ECAL isolation variable (barrel) [GeV]");
            trkisobarrel_before_ = theDbe->book1D("TRKISOBARREL_BEFORECUTS",chtitle,50,0.,50.);
            trkisobarrel_after_ = theDbe->book1D("TRKISOBARREL_LASTCUT",chtitle,50,0.,50.);

            snprintf(chtitle, 255, "Absolute electron track isolation variable (endcap) [GeV]");
            trkisoendcap_before_ = theDbe->book1D("TRKISOENDCAP_BEFORECUTS",chtitle,50,0.,50.);
            trkisoendcap_after_ = theDbe->book1D("TRKISOENDCAP_LASTCUT",chtitle,50,0.,50.);

//             snprintf(chtitle, 255, "Trigger response (bit %s)", muonTrig_.data());
//             trig_before_ = theDbe->book1D("TRIG_BEFORECUTS",chtitle,2,-0.5,1.5);
//             trig_after_ = theDbe->book1D("TRIG_LASTCUT",chtitle,2,-0.5,1.5);

            snprintf(chtitle, 255, "Trigger response (bit %s)", elecTrig_.data());
            trig_before_ = theDbe->book1D("TRIG_BEFORECUTS",chtitle,2,-0.5,1.5);
            trig_after_ = theDbe->book1D("TRIG_LASTCUT",chtitle,2,-0.5,1.5);

            snprintf(chtitle, 255, "Di-electron invariant mass [GeV]");
            invmass_before_ = theDbe->book1D("INVMASS_BEFORECUTS",chtitle,100,0.,200.);
            invmass_after_ = theDbe->book1D("INVMASS_AFTERCUTS",chtitle,100,0.,200.);

	    snprintf(chtitle, 255, "Number of electrons in event");
	    nelectrons_before_ = theDbe->book1D("NELECTRONS_BEFORECUTS",chtitle,10,-0.5,9.5);
	    nelectrons_after_ = theDbe->book1D("NELECTRONS_AFTERCUTS",chtitle,10,-0.5,9.5);

//             snprintf(chtitle, 255, "Transverse mass (%s) [GeV]", metTag_.label().data());
//             mt_before_ = theDbe->book1D("MT_BEFORECUTS",chtitle,150,0.,300.);
//             mt_after_ = theDbe->book1D("MT_LASTCUT",chtitle,150,0.,300.);

//             snprintf(chtitle, 255, "Missing transverse energy (%s) [GeV]", metTag_.label().data());
//             met_before_ = theDbe->book1D("MET_BEFORECUTS",chtitle,100,0.,200.);
//             met_after_ = theDbe->book1D("MET_LASTCUT",chtitle,100,0.,200.);

//             snprintf(chtitle, 255, "MU-MET (%s) acoplanarity", metTag_.label().data());
//             acop_before_ = theDbe->book1D("ACOP_BEFORECUTS",chtitle,50,0.,M_PI);
//             acop_after_ = theDbe->book1D("ACOP_LASTCUT",chtitle,50,0.,M_PI);

//             snprintf(chtitle, 255, "Z rejection: number of muons above %.2f GeV", ptThrForZ1_);
//             nz1_before_ = theDbe->book1D("NZ1_BEFORECUTS",chtitle,10,-0.5,9.5);
//             nz1_after_ = theDbe->book1D("NZ1_LASTCUT",chtitle,10,-0.5,9.5);

//             snprintf(chtitle, 255, "Z rejection: number of muons above %.2f GeV", ptThrForZ2_);
//             nz2_before_ = theDbe->book1D("NZ2_BEFORECUTS",chtitle,10,-0.5,9.5);
//             nz2_after_ = theDbe->book1D("NZ2_LASTCUT",chtitle,10,-0.5,9.5);

//             snprintf(chtitle, 255, "Number of jets (%s) above %.2f GeV", jetTag_.label().data(), eJetMin_);
//             njets_before_ = theDbe->book1D("NJETS_BEFORECUTS",chtitle,10,-0.5,9.5);
//             njets_after_ = theDbe->book1D("NJETS_LASTCUT",chtitle,10,-0.5,9.5);

      }
}


void EwkElecDQM::endJob() {
}

void EwkElecDQM::endRun(const Run& r, const EventSetup&) {

  // overall
  double all = nall;
  double esel = nsel/all;
  LogVerbatim("") << "\n>>>>>> SELECTION SUMMARY BEGIN >>>>>>>>>>>>>>>";
  LogVerbatim("") << "Total number of events analyzed: " << nall << " [events]";
  LogVerbatim("") << "Total number of events selected: " << nsel << " [events]";
  LogVerbatim("") << "Overall efficiency:             " << "(" << setprecision(4) << esel*100. <<" +/- "<< setprecision(2) << sqrt(esel*(1-esel)/all)*100. << ")%";
  
  double erec = nrec/all;
  double eeid = neid/all;
  double eiso = niso/all;
//   double ehlt = nhlt/all;
//   double emet = nmet/all;
  
  // general reconstruction step??
  double num = nrec;
  double eff = erec;
  double err = sqrt(eff*(1-eff)/all);
  LogVerbatim("") << "Passing Pt/Eta/Quality cuts:    " << num << " [events], (" << setprecision(4) << eff*100. <<" +/- "<< setprecision(2) << err*100. << ")%";

  // electron ID step  
  num = neid;
  eff = eeid;
  err = sqrt(eff*(1-eff)/all);
  double effstep = 0.;
  double errstep = 0.;
  if (nrec>0) effstep = eeid/erec;
  if (nrec>0) errstep = sqrt(effstep*(1-effstep)/nrec);
  LogVerbatim("") << "Passing eID cuts:         " << num << " [events], (" << setprecision(4) << eff*100. <<" +/- "<< setprecision(2) << err*100. << ")%, to previous step: (" <<  setprecision(4) << effstep*100. << " +/- "<< setprecision(2) << errstep*100. <<")%";
  
  // isolation step  
  num = niso;
  eff = eiso;
  err = sqrt(eff*(1-eff)/all);
  effstep = 0.;
  errstep = 0.;
  if (neid>0) effstep = eiso/eeid;
  if (neid>0) errstep = sqrt(effstep*(1-effstep)/neid);
  LogVerbatim("") << "Passing isolation cuts:         " << num << " [events], (" << setprecision(4) << eff*100. <<" +/- "<< setprecision(2) << err*100. << ")%, to previous step: (" <<  setprecision(4) << effstep*100. << " +/- "<< setprecision(2) << errstep*100. <<")%";
  
//   // trigger step
//   num = nhlt;
//   eff = ehlt;
//   err = sqrt(eff*(1-eff)/all);
//   effstep = 0.;
//   errstep = 0.;
//   if (niso>0) effstep = ehlt/eiso;
//   if (niso>0) errstep = sqrt(effstep*(1-effstep)/niso);
//   LogVerbatim("") << "Passing HLT criteria:           " << num << " [events], (" << setprecision(4) << eff*100. <<" +/- "<< setprecision(2) << err*100. << ")%, to previous step: (" <<  setprecision(4) << effstep*100. << " +/- "<< setprecision(2) << errstep*100. <<")%";
  
  // trigger step
  num = nsel;
  eff = esel;
  err = sqrt(eff*(1-eff)/all);
  effstep = 0.;
  errstep = 0.;
  if (niso>0) effstep = esel/eiso;
  if (niso>0) errstep = sqrt(effstep*(1-effstep)/niso);
  LogVerbatim("") << "Passing HLT criteria:           " << num << " [events], (" << setprecision(4) << eff*100. <<" +/- "<< setprecision(2) << err*100. << ")%, to previous step: (" <<  setprecision(4) << effstep*100. << " +/- "<< setprecision(2) << errstep*100. <<")%";
  
//   // met/acoplanarity cuts 
//   num = nmet;
//   eff = emet;
//   err = sqrt(eff*(1-eff)/all);
//   effstep = 0.;
//   errstep = 0.;
//   if (nhlt>0) effstep = emet/ehlt;
//   if (nhlt>0) errstep = sqrt(effstep*(1-effstep)/nhlt);
//   LogVerbatim("") << "Passing MET/acoplanarity cuts:  " << num << " [events], (" << setprecision(4) << eff*100. <<" +/- "<< setprecision(2) << err*100. << ")%, to previous step: (" <<  setprecision(4) << effstep*100. << " +/- "<< setprecision(2) << errstep*100. <<")%";
  
//   // Z/top selection cuts ALSO LAST STEP so "sel" for "selection"
//   num = nsel;
//   eff = esel;
//   err = sqrt(eff*(1-eff)/all);
//   effstep = 0.;
//   errstep = 0.;
//   if (nmet>0) effstep = esel/emet;
//   if (nmet>0) errstep = sqrt(effstep*(1-effstep)/nmet);
//   LogVerbatim("") << "Passing Z/top rejection cuts:   " << num << " [events], (" << setprecision(4) << eff*100. <<" +/- "<< setprecision(2) << err*100. << ")%, to previous step: (" <<  setprecision(4) << effstep*100. << " +/- "<< setprecision(2) << errstep*100. <<")%";
  
  LogVerbatim("") << ">>>>>> SELECTION SUMMARY END   >>>>>>>>>>>>>>>\n";
}

void EwkElecDQM::analyze (const Event & ev, const EventSetup &) {
      
      // Reset global event selection flags
      bool rec_sel = false;
      bool eid_sel = false;
      bool iso_sel = false;
//       bool hlt_sel = false;
//       bool met_sel = false;
      bool all_sel = false;

//       // Muon collection
//       Handle<View<Muon> > muonCollection;
//       if (!ev.getByLabel(muonTag_, muonCollection)) {
//             LogWarning("") << ">>> Muon collection does not exist !!!";
//             return;
//       }
//       unsigned int muonCollectionSize = muonCollection->size();

      // Electron collection
      Handle<View<GsfElectron> > electronCollection;
      if (!ev.getByLabel(elecTag_, electronCollection)) {
            LogWarning("") << ">>> Electron collection does not exist !!!";
            return;
      }
      unsigned int electronCollectionSize = electronCollection->size();

      // Beam spot
      Handle<reco::BeamSpot> beamSpotHandle;
      if (!ev.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle)) {
            LogWarning("") << ">>> No beam spot found !!!";
            return;
      }
  
//       // MET
//       double met_px = 0.;
//       double met_py = 0.;
//       Handle<View<MET> > metCollection;
//       if (!ev.getByLabel(metTag_, metCollection)) {
//             LogWarning("") << ">>> MET collection does not exist !!!";
//             return;
//       }
//       const MET& met = metCollection->at(0);
//       met_px = met.px();
//       met_py = met.py();
//       if (!metIncludesMuons_) {
//             for (unsigned int i=0; i<muonCollectionSize; i++) {
//                   const Muon& mu = muonCollection->at(i);
//                   if (!mu.isGlobalMuon()) continue;
//                   met_px -= mu.px();
//                   met_py -= mu.py();
//             }
//       }
//       double met_et = sqrt(met_px*met_px+met_py*met_py);
//       LogTrace("") << ">>> MET, MET_px, MET_py: " << met_et << ", " << met_px << ", " << met_py << " [GeV]";
//       met_before_->Fill(met_et);

      // Trigger
      Handle<TriggerResults> triggerResults;
      TriggerNames trigNames;
      if (!ev.getByLabel(trigTag_, triggerResults)) {
            LogWarning("") << ">>> TRIGGER collection does not exist !!!";
            return;
      }
      trigNames.init(*triggerResults);
      bool trigger_fired = false;
      /*
      for (unsigned int i=0; i<triggerResults->size(); i++) {
            if (triggerResults->accept(i)) {
                  LogTrace("") << "Accept by: " << i << ", Trigger: " << trigNames.triggerName(i);
            }
      }
      */

      // the following gives error on CRAFT08 data where itrig1=19 (vector index out of range)
      /*
      int itrig1 = trigNames.triggerIndex(muonTrig_);
      if (triggerResults->accept(itrig1)) trigger_fired = true;
      */
      //suggested replacement: lm250909
      for (unsigned int i=0; i<triggerResults->size(); i++) {
        std::string trigName = trigNames.triggerName(i);
	if ( trigName == elecTrig_ && triggerResults->accept(i)) 
	  {
	    trigger_fired = true;
	  }
      }


      LogTrace("") << ">>> Trigger bit: " << trigger_fired << " (" << elecTrig_ << ")";
      trig_before_->Fill(trigger_fired);

//       // Loop to reject/control Z->mumu is done separately
//       unsigned int nmuonsForZ1 = 0;
//       unsigned int nmuonsForZ2 = 0;
//       for (unsigned int i=0; i<muonCollectionSize; i++) {
//             const Muon& mu = muonCollection->at(i);
//             if (!mu.isGlobalMuon()) continue;
//             double pt = mu.pt();
//             if (pt>ptThrForZ1_) nmuonsForZ1++;
//             if (pt>ptThrForZ2_) nmuonsForZ2++;
//       }
//       LogTrace("") << "> Z rejection: muons above " << ptThrForZ1_ << " [GeV]: " << nmuonsForZ1;
//       LogTrace("") << "> Z rejection: muons above " << ptThrForZ2_ << " [GeV]: " << nmuonsForZ2;
//       nz1_before_->Fill(nmuonsForZ1);
//       nz2_before_->Fill(nmuonsForZ2);
      
//       // Jet collection
//       Handle<View<Jet> > jetCollection;
//       if (!ev.getByLabel(jetTag_, jetCollection)) {
//             LogError("") << ">>> JET collection does not exist !!!";
//             return;
//       }
//       unsigned int jetCollectionSize = jetCollection->size();
//       int njets = 0;
//       for (unsigned int i=0; i<jetCollectionSize; i++) {
//             const Jet& jet = jetCollection->at(i);
//             if (jet.et()>eJetMin_) njets++;
//       }
//       LogTrace("") << ">>> Total number of jets: " << jetCollectionSize;
//       LogTrace("") << ">>> Number of jets above " << eJetMin_ << " [GeV]: " << njets;
//       njets_before_->Fill(njets);

      // Start counting
      nall++;

      // Histograms per event should be done only once, so keep track of them
      bool hlt_hist_done = false;
      //bool minv_hist_done = false;
//       bool met_hist_done = false;
//       bool nz1_hist_done = false;
//       bool nz2_hist_done = false;
//       bool njets_hist_done = false;

      // Central selection criteria
      //const int NFLAGS = 13; // number of individual selection criteria
      const int NFLAGS = 8; // number of individual selection criteria
      // 0: pt cut           | rec
      // 1: eta cut          | rec
      // 2: sieie            | eid
      // 3: detain           | eid
      // 4: ecal iso         | iso
      // 5: hcal iso         | iso
      // 6: trk iso          | iso
      // 7: trigger fired    | hlt/all
      bool electron_sel[NFLAGS];

      // for invariant mass calculation
      // keep track of highest-pt electrons for initial (RECO) electrons 
      // and "good" electrons (passing all cuts)
      // first array dimension is for first or second good electron
      // second array dimension is for relevant quantities of good electron
      //    [0]: 1 for electron found or 0 for not found (if 0, no other quantities filled)
      //    [1]: mSqr
      //    [2]: E
      //    [3]: px
      //    [4]: py
      //    [5]: pz
      // inv mass = sqrt(m_1^2 + m_2^2 + 2*(E_1*E_2 - (px1*px2 + py1*py2 + pz1+pz2) ) )
      double electron[2][6];
      double goodElectron[2][6];
      nGoodElectrons = 0;
      for (unsigned int i = 0; i < 2; i++)
	{
	  for (unsigned int j = 0; j < 6; j++)
	    {
	      electron[i][j] = 0.;
	      goodElectron[i][j] = 0.;
	    }
	}

      for (unsigned int i=0; i<electronCollectionSize; i++) 
	{
	  for (int j=0; j<NFLAGS; ++j) 
	    {
	      electron_sel[j] = false;
            }
	  
	  const GsfElectron& elec = electronCollection->at(i);
	  //if (!mu.isGlobalMuon()) continue;
	  //if (mu.globalTrack().isNull()) continue;
	  //if (mu.innerTrack().isNull()) continue;
	  
	  LogTrace("") << "> elec: processing electron number " << i << "...";
	  //reco::TrackRef gm = mu.globalTrack();
	  //reco::TrackRef tk = mu.innerTrack();
	  // should have stuff for electron track?
	  
	  if (i < 2)
	    {
	      electron[i][0] = 1.;
	      electron[i][1] = elec.massSqr();
	      electron[i][2] = elec.energy();
	      electron[i][3] = elec.px();
	      electron[i][4] = elec.py();
	      electron[i][5] = elec.pz();
	    }

	  // Pt,eta cuts
	  double pt = elec.pt();
	  double eta = elec.eta();
	  LogTrace("") << "\t... pt, eta: " << pt << " [GeV], " << eta;;
	  if (pt>ptCut_) electron_sel[0] = true; 
	  if (fabs(eta)<etaCut_) electron_sel[1] = true; 
	  
	  bool isBarrel = false;
	  bool isEndcap = false;
	  if (eta < 1.4442 && eta > -1.4442)
	    {
	      isBarrel = true;
	    }
	  else if ((eta > 1.56 && eta < 2.4) || (eta < -1.56 && eta > -2.4))
	    {
	      isEndcap = true;
	    }
	  
	  //             // d0, chi2, nhits quality cuts
	  //             double dxy = tk->dxy(beamSpotHandle->position());
	  //             double normalizedChi2 = gm->normalizedChi2();
	  //             double trackerHits = tk->numberOfValidHits();
	  //             LogTrace("") << "\t... dxy, normalizedChi2, trackerHits, isTrackerMuon?: " << dxy << " [cm], " << normalizedChi2 << ", " << trackerHits << ", " << mu.isTrackerMuon();
	  //             if (fabs(dxy)<dxyCut_) muon_sel[2] = true; 
	  //             if (normalizedChi2<normalizedChi2Cut_) muon_sel[3] = true; 
	  //             if (trackerHits>=trackerHitsCut_) muon_sel[4] = true; 
	  //             if (mu.isTrackerMuon()) muon_sel[5] = true; 
	  
	  pt_before_->Fill(pt);
	  eta_before_->Fill(eta);
	  //             dxy_before_->Fill(dxy);
	  //             chi2_before_->Fill(normalizedChi2);
	  //             nhits_before_->Fill(trackerHits);
	  //             tkmu_before_->Fill(mu.isTrackerMuon());
	  
	  
	  // Electron ID cuts
	  double sieie = (double) elec.sigmaIetaIeta();
	  double detain = (double) elec.deltaEtaSuperClusterTrackAtVtx(); // think this is detain
	  if (sieie < sieieCutBarrel_ && isBarrel) electron_sel[2] = true; 
	  if (sieie < sieieCutEndcap_ && isEndcap) electron_sel[2] = true; 
	  if (detain < detainCutBarrel_ && isBarrel) electron_sel[3] = true; 
	  if (detain < detainCutEndcap_ && isEndcap) electron_sel[3] = true; 
	  if (isBarrel)
	    {
	      LogTrace("") << "\t... sieie value " << sieie << " (barrel), pass? " << electron_sel[2]; 
	      LogTrace("") << "\t... detain value " << detain << " (barrel), pass? " << electron_sel[3]; 
	    }
	  else if (isEndcap)
	    {
	      LogTrace("") << "\t... sieie value " << sieie << " (endcap), pass? " << electron_sel[2]; 
	      LogTrace("") << "\t... detain value " << detain << " (endcap), pass? " << electron_sel[2]; 
	    }
	  
	  if (isBarrel) 
	    {
	      sieiebarrel_before_->Fill(sieie);
	      detainbarrel_before_->Fill(detain);
	    }
	  else if (isEndcap) 
	    {
	      sieieendcap_before_->Fill(sieie);
	      detainendcap_before_->Fill(detain);
	    }
	  
	  
	  
	  // Isolation cuts
	  //double isovar = mu.isolationR03().sumPt;
	  double ecalisovar = elec.dr03EcalRecHitSumEt();  // picked one set! 
	  double hcalisovar = elec.dr03HcalTowerSumEt();   // try others if 
	  double trkisovar = elec.dr04TkSumPt();           // doesn't work 
	  //if (isCombinedIso_) {
	  //isovar += mu.isolationR03().emEt;
	  //isovar += mu.isolationR03().hadEt;
	  //}
	  //if (isRelativeIso_) isovar /= pt;
	  if (ecalisovar<ecalIsoCutBarrel_ && isBarrel) electron_sel[4] = true; 
	  if (ecalisovar<ecalIsoCutEndcap_ && isEndcap) electron_sel[4] = true; 
	  if (hcalisovar<hcalIsoCutBarrel_ && isBarrel) electron_sel[5] = true; 
	  if (hcalisovar<hcalIsoCutEndcap_ && isEndcap) electron_sel[5] = true; 
	  if (trkisovar<trkIsoCutBarrel_ && isBarrel) electron_sel[6] = true;   
	  if (trkisovar<trkIsoCutEndcap_ && isEndcap) electron_sel[6] = true;   
	  if (isBarrel)
	    {
	      LogTrace("") << "\t... ecal isolation value " << ecalisovar << " (barrel), pass? " << electron_sel[4]; 
	      LogTrace("") << "\t... hcal isolation value " << hcalisovar << " (barrel), pass? " << electron_sel[5];
	      LogTrace("") << "\t... trk isolation value " << trkisovar << " (barrel), pass? " << electron_sel[6];
	    }
	  else if (isEndcap)
	    {
	      LogTrace("") << "\t... ecal isolation value " << ecalisovar << " (endcap), pass? " << electron_sel[4]; 
	      LogTrace("") << "\t... hcal isolation value " << hcalisovar << " (endcap), pass? " << electron_sel[5];
	      LogTrace("") << "\t... trk isolation value " << trkisovar << " (endcap), pass? " << electron_sel[6];
	    }
	  
	  //iso_before_->Fill(isovar);
	  if (isBarrel) 
	    {
	      ecalisobarrel_before_->Fill(ecalisovar);
	      hcalisobarrel_before_->Fill(hcalisovar);
	      trkisobarrel_before_->Fill(trkisovar);
	    }
	  else if (isEndcap) 
	    {
	      ecalisoendcap_before_->Fill(ecalisovar);
	      hcalisoendcap_before_->Fill(hcalisovar);
	      trkisoendcap_before_->Fill(trkisovar);
	    }
	  
	  
	  // HLT 
	  if (trigger_fired) electron_sel[7] = true; 
	  
	  
	  //             // MET/MT cuts
	  //             double w_et = met_et+mu.pt();
	  //             double w_px = met_px+mu.px();
	  //             double w_py = met_py+mu.py();
	  
	  //             double massT = w_et*w_et - w_px*w_px - w_py*w_py;
	  //             massT = (massT>0) ? sqrt(massT) : 0;
	  
	  //             LogTrace("") << "\t... W mass, W_et, W_px, W_py: " << massT << ", " << w_et << ", " << w_px << ", " << w_py << " [GeV]";
	  //             if (massT>mtMin_ && massT<mtMax_) muon_sel[8] = true; 
	  //             mt_before_->Fill(massT);
	  //             if (met_et>metMin_ && met_et<metMax_) muon_sel[9] = true; 
	  
	  //             // Acoplanarity cuts
	  //             Geom::Phi<double> deltaphi(mu.phi()-atan2(met_py,met_px));
	  //             double acop = deltaphi.value();
	  //             if (acop<0) acop = - acop;
	  //             acop = M_PI - acop;
	  //             LogTrace("") << "\t... acoplanarity: " << acop;
	  //             if (acop<acopCut_) muon_sel[10] = true; 
	  //             acop_before_->Fill(acop);
	  
	  //             // Remaining flags (from global event information)
	  //             if (nmuonsForZ1<1 || nmuonsForZ2<2) muon_sel[11] = true; 
	  //             if (njets<=nJetMax_) muon_sel[12] = true; 
	  
	  // Collect necessary flags "per electron"
	  int flags_passed = 0;
	  bool rec_sel_this = true;
	  bool eid_sel_this = true;
	  bool iso_sel_this = true;
	  bool all_sel_this = true;
	  for (int j=0; j<NFLAGS; ++j) 
	    {
	      if (electron_sel[j]) flags_passed += 1;
	      if (j<2 && !electron_sel[j]) rec_sel_this = false;
	      if (j<4 && !electron_sel[j]) eid_sel_this = false;
	      if (j<7 && !electron_sel[j]) iso_sel_this = false;
	      if (!electron_sel[j]) all_sel_this = false;
	    }
	  
	  if (all_sel_this)
	    {
	      if (nGoodElectrons < 2)
		{
		  goodElectron[nGoodElectrons][0] = 1.;
		  goodElectron[nGoodElectrons][1] = elec.massSqr();
		  goodElectron[nGoodElectrons][2] = elec.energy();
		  goodElectron[nGoodElectrons][3] = elec.px();
		  goodElectron[nGoodElectrons][4] = elec.py();
		  goodElectron[nGoodElectrons][5] = elec.pz();
		}
	      nGoodElectrons++;
	    }

	  //             // "rec" => pt,eta and quality cuts are satisfied
	  //             if (rec_sel_this) rec_sel = true;
	  //             // "iso" => "rec" AND "muon is isolated"
	  //             if (iso_sel_this) iso_sel = true;
	  //             // "hlt" => "iso" AND "event is triggered"
	  //             if (hlt_sel_this) hlt_sel = true;
	  //             // "met" => "hlt" AND "MET/MT and acoplanarity cuts"
	  //             if (met_sel_this) met_sel = true;
	  //             // "all" => "met" AND "Z/top rejection cuts"
	  //             if (all_sel_this) all_sel = true;
	  
	  // "rec" => pt,eta cuts are satisfied
	  if (rec_sel_this) rec_sel = true;
	  // "eid" => "rec" AND "electron passes ID"
	  if (eid_sel_this) iso_sel = true;
	  // "iso" => "eid" AND "electron is isolated"
	  if (iso_sel_this) iso_sel = true;
	  // "all" => "iso" AND "event is triggered"
	  if (all_sel_this) all_sel = true;

	  // Do N-1 histograms now (and only once for global event quantities)
	  if (flags_passed >= (NFLAGS-1)) 
	    {
	      if (!electron_sel[0] || flags_passed==NFLAGS) 
		{
		  pt_after_->Fill(pt);
		}
	      if (!electron_sel[1] || flags_passed==NFLAGS) 
		{
		  eta_after_->Fill(eta);
		}
	      if (!electron_sel[2] || flags_passed==NFLAGS) 
		{
		  if (isBarrel)
		    {
		      sieiebarrel_after_->Fill(sieie);
		    }
		  else if (isEndcap)
		    {
		      sieieendcap_after_->Fill(sieie);
		    }
		}
	      if (!electron_sel[3] || flags_passed==NFLAGS) 
		{
		  if (isBarrel)
		    {
		      detainbarrel_after_->Fill(detain);
		    }
		  else if (isEndcap)
		    {
		      detainendcap_after_->Fill(detain);
		    }
		}
	      if (!electron_sel[4] || flags_passed==NFLAGS) 
		{
		  if (isBarrel)
		    {
		      ecalisobarrel_after_->Fill(ecalisovar);
		    }
		  else if (isEndcap)
		    {
		      ecalisoendcap_after_->Fill(ecalisovar);
		    }
		}
	      if (!electron_sel[5] || flags_passed==NFLAGS) 
		{
		  if (isBarrel)
		    {
		      hcalisobarrel_after_->Fill(hcalisovar);
		    }
		  else if (isEndcap)
		    {
		      hcalisoendcap_after_->Fill(hcalisovar);
		    }
		}
	      if (!electron_sel[6] || flags_passed==NFLAGS) 
		{
		  if (isBarrel)
		    {
		      trkisobarrel_after_->Fill(trkisovar);
		    }
		  else if (isEndcap)
		    {
		      trkisoendcap_after_->Fill(trkisovar);
		    }
		}
	      // 		if (!electron_sel[3] || flags_passed==NFLAGS) 
	      // 		  {
	      // 		    detain_after_->Fill(detain);
	      // 		  }
	      // 		if (!electron_sel[4] || flags_passed==NFLAGS) 
	      // 		  {
	      // 		    ecaliso_after_->Fill(trackerHits);
	      // 		  }
	      // 		if (!electron_sel[5] || flags_passed==NFLAGS) 
	      // 		  {
	      // 		    tkelectr_after_->Fill(electr.isTrackerElectron());
	      // 		  }
	      // 		if (!electron_sel[6] || flags_passed==NFLAGS) 
	      // 		  {
	      // 		    iso_after_->Fill(isovar);
	      // 		  }
	      if (!electron_sel[7] || flags_passed==NFLAGS) 
		{
		  if (!hlt_hist_done) 
		    {
		      trig_after_->Fill(trigger_fired);
		    }
		}
	      hlt_hist_done = true;
	      //                   if (!muon_sel[8] || flags_passed==NFLAGS) 
	      //                         mt_after_->Fill(massT);
	      //                   if (!muon_sel[9] || flags_passed==NFLAGS) 
	      //                         if (!met_hist_done) met_after_->Fill(met_et);
	      //                         met_hist_done = true;
	      //                   if (!muon_sel[10] || flags_passed==NFLAGS) 
	      //                         acop_after_->Fill(acop);
	      //                   if (!muon_sel[11] || flags_passed==NFLAGS) 
	      //                         if (!nz1_hist_done) nz1_after_->Fill(nmuonsForZ1);
	      //                         nz1_hist_done = true;
	      //                   if (!muon_sel[11] || flags_passed==NFLAGS) 
	      //                         if (!nz2_hist_done) nz2_after_->Fill(nmuonsForZ2);
	      //                         nz2_hist_done = true;
	      //                   if (!muon_sel[12] || flags_passed==NFLAGS) 
	      //                         if (!njets_hist_done) njets_after_->Fill(njets);
	      //                         njets_hist_done = true;

            } // end N-1 histos block
	  
	} // end loop through electrons

      // inv mass = sqrt(m_1^2 + m_2^2 + 2*(E_1*E_2 - (px1*px2 + py1*py2 + pz1+pz2) ) )
      double invMass;

      nelectrons_before_->Fill(electronCollectionSize);
      if (electronCollectionSize > 1)
	{
	  invMass = sqrt(electron[0][1] + electron[1][1] + 2*(electron[0][2]*electron[1][2] - (electron[0][3]*electron[1][3] + electron[0][4]*electron[1][4] + electron[0][5]*electron[1][5]) ) );
	  invmass_before_->Fill(invMass);
	}

      nelectrons_after_->Fill(nGoodElectrons);
      if (nGoodElectrons > 1)
	{
	  invMass = sqrt(goodElectron[0][1] + goodElectron[1][1] + 2*(goodElectron[0][2]*goodElectron[1][2] - (goodElectron[0][3]*goodElectron[1][3] + goodElectron[0][4]*goodElectron[1][4] + goodElectron[0][5]*goodElectron[1][5]) ) );
	  invmass_after_->Fill(invMass);
	}

      // Collect final flags
      if (rec_sel) nrec++;
      if (eid_sel) neid++;
      if (iso_sel) niso++;
      //      if (hlt_sel) nhlt++;
      //      if (met_sel) nmet++;

      if (all_sel) {
            nsel++;
            LogTrace("") << ">>>> Event ACCEPTED";
      } else {
            LogTrace("") << ">>>> Event REJECTED";
      }

      return;

}
