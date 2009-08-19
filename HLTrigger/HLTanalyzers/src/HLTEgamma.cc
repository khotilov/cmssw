#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <cstdlib>
#include <cstring>

#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "HLTrigger/HLTanalyzers/interface/HLTEgamma.h"



#include "RecoEgamma/EgammaTools/interface/ECALPositionCalculator.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "HLTMessages.h"

static const size_t kMaxEl     = 10000;
static const size_t kMaxPhot   = 10000;
static const size_t kMaxhPhot  =   500;
static const size_t kMaxhEle   =   500;
static const size_t kMaxhEleLW =   500;

HLTEgamma::HLTEgamma() {
}

/*  Setup the analysis to put the branch-variables size_to the tree. */
void HLTEgamma::setup(const edm::ParameterSet& pSet, TTree* HltTree)
{
  elpt              = new float[kMaxEl];
  elphi             = new float[kMaxEl];
  eleta             = new float[kMaxEl];
  elet              = new float[kMaxEl];
  ele               = new float[kMaxEl];
  eleId             = new int[kMaxEl];// RL  + 2*RT + 4*L +  4*T 
  
  photonpt          = new float[kMaxPhot];
  photonphi         = new float[kMaxPhot];
  photoneta         = new float[kMaxPhot];
  photonet          = new float[kMaxPhot];
  photone           = new float[kMaxPhot];

  hphotet           = new float[kMaxhPhot];
  hphoteta          = new float[kMaxhPhot];
  hphotphi          = new float[kMaxhPhot];
  hphoteiso         = new float[kMaxhPhot];
  hphothiso         = new float[kMaxhPhot];
  hphottiso         = new float[kMaxhPhot];
  hphotl1iso        = new int[kMaxhPhot];
  hphotClusShap     = new float[kMaxhPhot];

  heleet            = new float[kMaxhEle];
  heleeta           = new float[kMaxhEle];
  helephi           = new float[kMaxhEle];
  heleE             = new float[kMaxhEle];
  helep             = new float[kMaxhEle];
  helehiso          = new float[kMaxhEle];
  heletiso          = new float[kMaxhEle];
  helel1iso         = new int[kMaxhEle];
  helePixelSeeds    = new int[kMaxhEle];
  heleNewSC         = new int[kMaxhEle];
  heleClusShap      = new float[kMaxhEle];
  heleDeta          = new float[kMaxhEle];
  heleDphi          = new float[kMaxhEle];

  heleetLW          = new float[kMaxhEleLW];
  heleetaLW         = new float[kMaxhEleLW];
  helephiLW         = new float[kMaxhEleLW];
  heleELW           = new float[kMaxhEleLW];
  helepLW           = new float[kMaxhEleLW];
  helehisoLW        = new float[kMaxhEleLW];
  heletisoLW        = new float[kMaxhEleLW];
  helel1isoLW       = new int[kMaxhEleLW];
  helePixelSeedsLW  = new int[kMaxhEleLW];
  heleNewSCLW       = new int[kMaxhEleLW];
  heleClusShapLW    = new float[kMaxhEleLW];
  heleDetaLW        = new float[kMaxhEleLW];
  heleDphiLW        = new float[kMaxhEleLW];

  nele      = 0;
  nphoton   = 0;
  nhltgam   = 0;
  nhltele   = 0;
  nhlteleLW = 0;

  // Egamma-specific branches of the tree
  HltTree->Branch("NrecoElec",          & nele,             "NrecoElec/I");
  HltTree->Branch("recoElecPt",         elpt,               "recoElecPt[NrecoElec]/F");
  HltTree->Branch("recoElecPhi",        elphi,              "recoElecPhi[NrecoElec]/F");
  HltTree->Branch("recoElecEta",        eleta,              "recoElecEta[NrecoElec]/F");
  HltTree->Branch("recoElecEt",         elet,               "recoElecEt[NrecoElec]/F");
  HltTree->Branch("recoElecE",          ele,                "recoElecE[NrecoElec]/F");
  HltTree->Branch("recoElecEleID",      eleId,              "recoElecEleID[NrecoElec]/I");

  HltTree->Branch("NrecoPhot",          &nphoton,           "NrecoPhot/I");
  HltTree->Branch("recoPhotPt",         photonpt,           "recoPhotPt[NrecoPhot]/F");
  HltTree->Branch("recoPhotPhi",        photonphi,          "recoPhotPhi[NrecoPhot]/F");
  HltTree->Branch("recoPhotEta",        photoneta,          "recoPhotEta[NrecoPhot]/F");
  HltTree->Branch("recoPhotEt",         photonet,           "recoPhotEt[NrecoPhot]/F");
  HltTree->Branch("recoPhotE",          photone,            "recoPhotE[NrecoPhot]/F");

  HltTree->Branch("NohPhot",            & nhltgam,          "NohPhot/I");
  HltTree->Branch("ohPhotEt",           hphotet,            "ohPhotEt[NohPhot]/F");
  HltTree->Branch("ohPhotEta",          hphoteta,           "ohPhotEta[NohPhot]/F");
  HltTree->Branch("ohPhotPhi",          hphotphi,           "ohPhotPhi[NohPhot]/F");
  HltTree->Branch("ohPhotEiso",         hphoteiso,          "ohPhotEiso[NohPhot]/F");
  HltTree->Branch("ohPhotHiso",         hphothiso,          "ohPhotHiso[NohPhot]/F");
  HltTree->Branch("ohPhotTiso",         hphottiso,          "ohPhotTiso[NohPhot]/F");
  HltTree->Branch("ohPhotL1iso",        hphotl1iso,         "ohPhotL1iso[NohPhot]/I");
  HltTree->Branch("ohPhotClusShap",     hphotClusShap,      "ohPhotClusShap[NohPhot]/F");

  HltTree->Branch("NohEle",             & nhltele,          "NohEle/I");
  HltTree->Branch("ohEleEt",            heleet,             "ohEleEt[NohEle]/F");
  HltTree->Branch("ohEleEta",           heleeta,            "ohEleEta[NohEle]/F");
  HltTree->Branch("ohElePhi",           helephi,            "ohElePhi[NohEle]/F");
  HltTree->Branch("ohEleE",             heleE,              "ohEleE[NohEle]/F");
  HltTree->Branch("ohEleP",             helep,              "ohEleP[NohEle]/F");
  HltTree->Branch("ohEleHiso",          helehiso,           "ohEleHiso[NohEle]/F");
  HltTree->Branch("ohEleTiso",          heletiso,           "ohEleTiso[NohEle]/F");
  HltTree->Branch("ohEleL1iso",         helel1iso,          "ohEleLiso[NohEle]/I");
  HltTree->Branch("ohElePixelSeeds",    helePixelSeeds,     "ohElePixelSeeds[NohEle]/I");
  HltTree->Branch("ohEleNewSC",         heleNewSC,          "ohEleNewSC[NohEle]/I");
  HltTree->Branch("ohEleClusShap",      heleClusShap,       "ohEleClusShap[NohEle]/F");
  HltTree->Branch("ohEleDeta",          heleDeta,           "ohEleDeta[NohEle]/F");
  HltTree->Branch("ohEleDphi",          heleDphi,           "ohEleDphi[NohEle]/F");

  HltTree->Branch("NohEleLW",           & nhlteleLW,        "NohEleLW/I");
  HltTree->Branch("ohEleEtLW",          heleetLW,           "ohEleEtLW[NohEleLW]/F");
  HltTree->Branch("ohEleEtaLW",         heleetaLW,          "ohEleEtaLW[NohEleLW]/F");
  HltTree->Branch("ohElePhiLW",         helephiLW,          "ohElePhiLW[NohEleLW]/F");
  HltTree->Branch("ohEleELW",           heleELW,            "ohEleELW[NohEleLW]/F");
  HltTree->Branch("ohElePLW",           helepLW,            "ohElePLW[NohEleLW]/F");
  HltTree->Branch("ohEleHisoLW",        helehisoLW,         "ohEleHisoLW[NohEleLW]/F");
  HltTree->Branch("ohEleTisoLW",        heletisoLW,         "ohEleTisoLW[NohEleLW]/F");
  HltTree->Branch("ohEleL1isoLW",       helel1isoLW,        "ohEleLisoLW[NohEleLW]/I");
  HltTree->Branch("ohElePixelSeedsLW",  helePixelSeedsLW,   "ohElePixelSeedsLW[NohEleLW]/I");
  HltTree->Branch("ohEleNewSCLW",       heleNewSCLW,        "ohEleNewSCLW[NohEleLW]/I");
  HltTree->Branch("ohEleClusShapLW",    heleClusShapLW,     "ohEleClusShapLW[NohEleLW]/F");
  HltTree->Branch("ohEleDetaLW",        heleDetaLW,         "ohEleDetaLW[NohEleLW]/F");
  HltTree->Branch("ohEleDphiLW",        heleDphiLW,         "ohEleDphiLW[NohEleLW]/F");
}

void HLTEgamma::clear(void)
{
  std::memset(elpt,             '\0', kMaxEl     * sizeof(float));
  std::memset(elphi,            '\0', kMaxEl     * sizeof(float));
  std::memset(eleta,            '\0', kMaxEl     * sizeof(float));
  std::memset(elet,             '\0', kMaxEl     * sizeof(float));
  std::memset(ele,              '\0', kMaxEl     * sizeof(float));
  std::memset(ele,              '\0', kMaxEl     * sizeof(int));

  std::memset(photonpt,         '\0', kMaxPhot   * sizeof(float));
  std::memset(photonphi,        '\0', kMaxPhot   * sizeof(float));
  std::memset(photoneta,        '\0', kMaxPhot   * sizeof(float));
  std::memset(photonet,         '\0', kMaxPhot   * sizeof(float));
  std::memset(photone,          '\0', kMaxPhot   * sizeof(float));

  std::memset(hphotet,          '\0', kMaxhPhot  * sizeof(float));
  std::memset(hphoteta,         '\0', kMaxhPhot  * sizeof(float));
  std::memset(hphotphi,         '\0', kMaxhPhot  * sizeof(float));
  std::memset(hphoteiso,        '\0', kMaxhPhot  * sizeof(float));
  std::memset(hphothiso,        '\0', kMaxhPhot  * sizeof(float));
  std::memset(hphottiso,        '\0', kMaxhPhot  * sizeof(float));
  std::memset(hphotl1iso,       '\0', kMaxhPhot  * sizeof(int));
  std::memset(hphotClusShap,    '\0', kMaxhPhot  * sizeof(float));

  std::memset(heleet,           '\0', kMaxhEle   * sizeof(float));
  std::memset(heleeta,          '\0', kMaxhEle   * sizeof(float));
  std::memset(helephi,          '\0', kMaxhEle   * sizeof(float));
  std::memset(heleE,            '\0', kMaxhEle   * sizeof(float));
  std::memset(helep,            '\0', kMaxhEle   * sizeof(float));
  std::memset(helehiso,         '\0', kMaxhEle   * sizeof(float));
  std::memset(heletiso,         '\0', kMaxhEle   * sizeof(float));
  std::memset(helel1iso,        '\0', kMaxhEle   * sizeof(int));
  std::memset(helePixelSeeds,   '\0', kMaxhEle   * sizeof(int));
  std::memset(heleNewSC,        '\0', kMaxhEle   * sizeof(int));
  std::memset(heleClusShap,     '\0', kMaxhEle  * sizeof(float));
  std::memset(heleDeta,         '\0', kMaxhEle  * sizeof(float));
  std::memset(heleDphi,         '\0', kMaxhEle  * sizeof(float));

  std::memset(heleetLW,         '\0', kMaxhEleLW * sizeof(float));
  std::memset(heleetaLW,        '\0', kMaxhEleLW * sizeof(float));
  std::memset(helephiLW,        '\0', kMaxhEleLW * sizeof(float));
  std::memset(heleELW,          '\0', kMaxhEleLW * sizeof(float));
  std::memset(helepLW,          '\0', kMaxhEleLW * sizeof(float));
  std::memset(helehisoLW,       '\0', kMaxhEleLW * sizeof(float));
  std::memset(heletisoLW,       '\0', kMaxhEleLW * sizeof(float));
  std::memset(helel1isoLW,      '\0', kMaxhEleLW * sizeof(int));
  std::memset(helePixelSeedsLW, '\0', kMaxhEleLW * sizeof(int));
  std::memset(heleNewSCLW,      '\0', kMaxhEleLW * sizeof(int));
  std::memset(heleClusShapLW,   '\0', kMaxhEleLW  * sizeof(float));
  std::memset(heleDetaLW,       '\0', kMaxhEleLW  * sizeof(float));
  std::memset(heleDphiLW,       '\0', kMaxhEleLW  * sizeof(float));

  nele      = 0;
  nphoton   = 0;
  nhltgam   = 0;
  nhltele   = 0;
  nhlteleLW = 0;
}

/* **Analyze the event** */
void HLTEgamma::analyze(const edm::Handle<reco::GsfElectronCollection>         & electrons,
                        const edm::Handle<reco::PhotonCollection>              & photons,
                        const edm::Handle<reco::ElectronCollection>            & electronIsoHandle,
                        const edm::Handle<reco::ElectronCollection>            & electronIsoHandleLW,
                        const edm::Handle<reco::ElectronCollection>            & electronNonIsoHandle,
                        const edm::Handle<reco::ElectronCollection>            & electronNonIsoHandleLW,
                        const edm::Handle<reco::ElectronIsolationMap>          & NonIsoTrackEleIsolMap,
                        const edm::Handle<reco::ElectronIsolationMap>          & NonIsoTrackEleIsolMapLW,
                        const edm::Handle<reco::ElectronIsolationMap>          & TrackEleIsolMap,
                        const edm::Handle<reco::ElectronIsolationMap>          & TrackEleIsolMapLW,
                        const edm::Handle<reco::ElectronSeedCollection>   & L1IsoPixelSeedsMap,
                        const edm::Handle<reco::ElectronSeedCollection>   & L1IsoPixelSeedsMapLW,
                        const edm::Handle<reco::ElectronSeedCollection>   & L1NonIsoPixelSeedsMap,
                        const edm::Handle<reco::ElectronSeedCollection>   & L1NonIsoPixelSeedsMapLW,
                        const edm::Handle<reco::RecoEcalCandidateCollection>   & recoIsolecalcands,
                        const edm::Handle<reco::RecoEcalCandidateCollection>   & recoNonIsolecalcands,
                        const edm::Handle<reco::RecoEcalCandidateIsolationMap> & EcalIsolMap,
                        const edm::Handle<reco::RecoEcalCandidateIsolationMap> & EcalNonIsolMap,
                        const edm::Handle<reco::RecoEcalCandidateIsolationMap> & HcalEleIsolMap,
                        const edm::Handle<reco::RecoEcalCandidateIsolationMap> & HcalEleNonIsolMap,
                        const edm::Handle<reco::RecoEcalCandidateIsolationMap> & HcalIsolMap,
                        const edm::Handle<reco::RecoEcalCandidateIsolationMap> & HcalNonIsolMap,
                        const edm::Handle<reco::RecoEcalCandidateIsolationMap> & TrackIsolMap,
                        const edm::Handle<reco::RecoEcalCandidateIsolationMap> & TrackNonIsolMap,
			EcalClusterLazyTools& lazyTools,
			const edm::ESHandle<MagneticField>& theMagField,
			reco::BeamSpot::Point & BSPosition, 
			std::vector<edm::Handle<edm::ValueMap<float> > > & eIDValueMap,
                        TTree* HltTree)
{
  // reset the tree variables
  clear();

  if (electrons.isValid()) {
    GsfElectronCollection myelectrons( electrons->begin(), electrons->end() );
    nele = myelectrons.size();
    std::sort(myelectrons.begin(), myelectrons.end(), EtGreater());
    int iel = 0;
    for (GsfElectronCollection::const_iterator i = myelectrons.begin(); i != myelectrons.end(); i++) {
      elpt[iel]  = i->pt();
      elphi[iel] = i->phi();
      eleta[iel] = i->eta();
      elet[iel]  = i->et();
      ele[iel]   = i->energy();
      iel++;
    }
  } else {
    nele = 0;
  }

// Use of electron-identification (= RobustLoose or RobustTight or Loose or Tight)
// variable, which does not have the Et-sorting
//     nele = electrons->size();
//     for (int iel = 0; iel < nele; iel++) {
//       edm::Ref<reco::GsfElectronCollection> electronRef(electrons,iel);
//       elpt[iel]  = electronRef->pt();
//       elphi[iel] = electronRef->phi();
//       eleta[iel] = electronRef->eta();
//       elet[iel]  = electronRef->et();
//       ele[iel]   = electronRef->energy();
//       int eleidtmp=0;
//       bool foundEleId = false;
//       int pos[4]={1,2,4,8};
//       for(int id=0; id<4; id++){
// 	if(eIDValueMap[id].isValid()){
// 	  foundEleId = true;
// 	  const edm::ValueMap<float> & eIDmap  = * eIDValueMap[id] ;
// 	  eleidtmp += pos[id]*int(eIDmap[electronRef]); 
// 	}
//       }
//       if (!foundEleId){eleidtmp = -1;}
//       eleId[iel]=eleidtmp;
//     }
//   } else {
//     nele = 0;
//   }

  if (photons.isValid()) {
    PhotonCollection myphotons(* photons);
    nphoton = myphotons.size();
    std::sort(myphotons.begin(), myphotons.end(), EtGreater());
    int ipho = 0;
    for (PhotonCollection::const_iterator i = myphotons.begin(); i!= myphotons.end(); i++) {
      photonpt[ipho] = i->pt();
      photonphi[ipho] = i->phi();
      photoneta[ipho] = i->eta();
      photonet[ipho] = i->et();
      photone[ipho] = i->energy();
      ipho++;
    }
  } else {
    nphoton = 0;
  }

  /////////////////////////////// Open-HLT Egammas ///////////////////////////////

  /////    Open-HLT photons /////////////////
  std::vector<OpenHLTPhoton> theHLTPhotons;
  MakeL1IsolatedPhotons(
      theHLTPhotons,
      recoIsolecalcands,
      EcalIsolMap,
      HcalIsolMap,
      TrackIsolMap,
      lazyTools);
  MakeL1NonIsolatedPhotons(
      theHLTPhotons,
      recoNonIsolecalcands,
      EcalNonIsolMap,
      HcalNonIsolMap,
      TrackNonIsolMap,
      lazyTools);

  std::sort(theHLTPhotons.begin(), theHLTPhotons.end(), EtGreater());
  nhltgam = theHLTPhotons.size();

  for (int u = 0; u < nhltgam; u++) {
    hphotet[u]    = theHLTPhotons[u].Et;
    hphoteta[u]   = theHLTPhotons[u].eta;
    hphotphi[u]   = theHLTPhotons[u].phi;
    hphoteiso[u]  = theHLTPhotons[u].ecalIsol;
    hphothiso[u]  = theHLTPhotons[u].hcalIsol;
    hphottiso[u]  = theHLTPhotons[u].trackIsol;
    hphotl1iso[u] = theHLTPhotons[u].L1Isolated;
    hphotClusShap[u] = theHLTPhotons[u].clusterShape;
  }

  /////    Open-HLT electrons /////////////////
  std::vector<OpenHLTElectron> theHLTElectrons;
  MakeL1IsolatedElectrons(
      theHLTElectrons,
      electronIsoHandle,
      recoIsolecalcands,
      HcalEleIsolMap,
      L1IsoPixelSeedsMap,
      TrackEleIsolMap,
      lazyTools,
      theMagField,
      BSPosition);
  MakeL1NonIsolatedElectrons(
      theHLTElectrons,
      electronNonIsoHandle,
      recoNonIsolecalcands,
      HcalEleNonIsolMap,
      L1NonIsoPixelSeedsMap,
      NonIsoTrackEleIsolMap,
      lazyTools,
      theMagField,
      BSPosition);

  std::sort(theHLTElectrons.begin(), theHLTElectrons.end(), EtGreater());
  nhltele = theHLTElectrons.size();

  for (int u = 0; u < nhltele; u++) {
    heleet[u]         = theHLTElectrons[u].Et;
    heleeta[u]        = theHLTElectrons[u].eta;
    helephi[u]        = theHLTElectrons[u].phi;
    heleE[u]          = theHLTElectrons[u].E;
    helep[u]          = theHLTElectrons[u].p;
    helehiso[u]       = theHLTElectrons[u].hcalIsol;
    helePixelSeeds[u] = theHLTElectrons[u].pixelSeeds;
    heletiso[u]       = theHLTElectrons[u].trackIsol;
    helel1iso[u]      = theHLTElectrons[u].L1Isolated;
    heleNewSC[u]      = theHLTElectrons[u].newSC;
    heleClusShap[u]   = theHLTElectrons[u].clusterShape;
    heleDeta[u]       = theHLTElectrons[u].Deta;
    heleDphi[u]       = theHLTElectrons[u].Dphi;
  }

  /////    Open-HLT electrons  LW /////////////////
  std::vector<OpenHLTElectron> theHLTElectronsLargeWindows;
  MakeL1IsolatedElectrons(
      theHLTElectronsLargeWindows,
      electronIsoHandleLW,
      recoIsolecalcands,
      HcalEleIsolMap,
      L1IsoPixelSeedsMapLW,
      TrackEleIsolMapLW,
      lazyTools,
      theMagField,
      BSPosition);
  MakeL1NonIsolatedElectrons(
      theHLTElectronsLargeWindows,
      electronNonIsoHandleLW,
      recoNonIsolecalcands,
      HcalEleNonIsolMap,
      L1NonIsoPixelSeedsMapLW,
      NonIsoTrackEleIsolMapLW,
      lazyTools,
      theMagField,
      BSPosition);

  std::sort(theHLTElectronsLargeWindows.begin(), theHLTElectronsLargeWindows.end(), EtGreater());
  nhlteleLW = theHLTElectronsLargeWindows.size();

  for (int u = 0; u < nhlteleLW; u++) {
    heleetLW[u]         = theHLTElectronsLargeWindows[u].Et;
    heleetaLW[u]        = theHLTElectronsLargeWindows[u].eta;
    helephiLW[u]        = theHLTElectronsLargeWindows[u].phi;
    heleELW[u]          = theHLTElectronsLargeWindows[u].E;
    helepLW[u]          = theHLTElectronsLargeWindows[u].p;
    helehisoLW[u]       = theHLTElectronsLargeWindows[u].hcalIsol;
    helePixelSeedsLW[u] = theHLTElectronsLargeWindows[u].pixelSeeds;
    heletisoLW[u]       = theHLTElectronsLargeWindows[u].trackIsol;
    helel1isoLW[u]      = theHLTElectronsLargeWindows[u].L1Isolated;
    heleNewSCLW[u]      = theHLTElectronsLargeWindows[u].newSC;
    heleClusShapLW[u]   = theHLTElectronsLargeWindows[u].clusterShape;
    heleDetaLW[u]       = theHLTElectronsLargeWindows[u].Deta;
    heleDphiLW[u]       = theHLTElectronsLargeWindows[u].Dphi;

  }
}

void HLTEgamma::MakeL1IsolatedPhotons(
    std::vector<OpenHLTPhoton> & theHLTPhotons,
    const edm::Handle<reco::RecoEcalCandidateCollection>   & recoIsolecalcands,
    const edm::Handle<reco::RecoEcalCandidateIsolationMap> & EcalIsolMap,
    const edm::Handle<reco::RecoEcalCandidateIsolationMap> & HcalIsolMap,
    const edm::Handle<reco::RecoEcalCandidateIsolationMap> & TrackIsolMap,
    EcalClusterLazyTools& lazyTools )
{
  // Iterator to the isolation-map
  reco::RecoEcalCandidateIsolationMap::const_iterator mapi;

  if (recoIsolecalcands.isValid()) {
    // loop over SuperCluster and fill the HLTPhotons
  

    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoIsolecalcands->begin();
         recoecalcand!= recoIsolecalcands->end(); recoecalcand++) {

      OpenHLTPhoton pho;
      pho.ecalIsol   = -999;
      pho.hcalIsol   = -999;
      pho.trackIsol  = -999;
      pho.clusterShape = -999;
      pho.L1Isolated = true;
      pho.Et         = recoecalcand->et();
      pho.eta        = recoecalcand->eta();
      pho.phi        = recoecalcand->phi();

      //Get the cluster shape
      std::vector<float> vCov = lazyTools.covariances( *(recoecalcand->superCluster()->seed()) );
      //      std::vector<float> vCov = lazyTools.localCovariances( *(recoecalcand->superCluster()->seed()) );
      double sigmaee = sqrt(vCov[0]);
      float EtaSC = fabs(recoecalcand->eta());
      if(EtaSC > 1.479 ) {//Endcap
        sigmaee = sigmaee - 0.02*(EtaSC - 2.3);
      }
      pho.clusterShape = sigmaee;

      // Method to get the reference to the candidate
      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoIsolecalcands, distance(recoIsolecalcands->begin(), recoecalcand));

      // First/Second member of the Map: Ref-to-Candidate(mapi)/Isolation(->val)
      // fill the ecal Isolation
      if (EcalIsolMap.isValid()) {
        mapi = (*EcalIsolMap).find(ref);
        if (mapi !=(*EcalIsolMap).end()) { pho.ecalIsol = mapi->val;}
      }
      // fill the hcal Isolation
      if (HcalIsolMap.isValid()) {
        mapi = (*HcalIsolMap).find(ref);
        if (mapi !=(*HcalIsolMap).end()) { pho.hcalIsol = mapi->val;}
      }
      // fill the track Isolation
      if (TrackIsolMap.isValid()) {
        mapi = (*TrackIsolMap).find(ref);
        if (mapi !=(*TrackIsolMap).end()) { pho.trackIsol = mapi->val;}
      }

      // store the photon into the vector
      theHLTPhotons.push_back(pho);
    }
  }
}

void HLTEgamma::MakeL1NonIsolatedPhotons(
    std::vector<OpenHLTPhoton> & theHLTPhotons,
    const edm::Handle<reco::RecoEcalCandidateCollection>   & recoNonIsolecalcands,
    const edm::Handle<reco::RecoEcalCandidateIsolationMap> & EcalNonIsolMap,
    const edm::Handle<reco::RecoEcalCandidateIsolationMap> & HcalNonIsolMap,
    const edm::Handle<reco::RecoEcalCandidateIsolationMap> & TrackNonIsolMap,
    EcalClusterLazyTools& lazyTools )
{
  reco::RecoEcalCandidateIsolationMap::const_iterator mapi;

  if (recoNonIsolecalcands.isValid()) {
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoNonIsolecalcands->begin();
         recoecalcand!= recoNonIsolecalcands->end(); recoecalcand++) {
      // loop over SuperCluster and fill the HLTPhotons
      OpenHLTPhoton pho;
      pho.ecalIsol   = -999;
      pho.hcalIsol   = -999;
      pho.trackIsol  = -999;
      pho.clusterShape = -999;
      pho.L1Isolated = false;
      pho.Et         = recoecalcand->et();
      pho.eta        = recoecalcand->eta();
      pho.phi        = recoecalcand->phi();

      //Get the cluster shape
      std::vector<float> vCov = lazyTools.covariances( *(recoecalcand->superCluster()->seed()) );
      //      std::vector<float> vCov = lazyTools.localCovariances( *(recoecalcand->superCluster()->seed()) );
      double sigmaee = sqrt(vCov[0]);
      float EtaSC = fabs(recoecalcand->eta());
      if(EtaSC > 1.479 ) {//Endcap
        sigmaee = sigmaee - 0.02*(EtaSC - 2.3);
      }
      pho.clusterShape = sigmaee;

      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoNonIsolecalcands, distance(recoNonIsolecalcands->begin(), recoecalcand));

      // fill the ecal Isolation
      if (EcalNonIsolMap.isValid()) {
        mapi = (*EcalNonIsolMap).find(ref);
        if (mapi !=(*EcalNonIsolMap).end()) { pho.ecalIsol = mapi->val;}
      }
      // fill the hcal Isolation
      if (HcalNonIsolMap.isValid()) {
        mapi = (*HcalNonIsolMap).find(ref);
        if (mapi !=(*HcalNonIsolMap).end()) { pho.hcalIsol = mapi->val;}
      }
      // fill the track Isolation
      if (TrackNonIsolMap.isValid()) {
        mapi = (*TrackNonIsolMap).find(ref);
        if (mapi !=(*TrackNonIsolMap).end()) { pho.trackIsol = mapi->val;}
      }

      // store the photon into the vector
      theHLTPhotons.push_back(pho);
    }
  }
}

void HLTEgamma::MakeL1IsolatedElectrons(
    std::vector<OpenHLTElectron> & theHLTElectrons,
    const edm::Handle<reco::ElectronCollection>            & electronIsoHandle,
    const edm::Handle<reco::RecoEcalCandidateCollection>   & recoIsolecalcands,
    const edm::Handle<reco::RecoEcalCandidateIsolationMap> & HcalEleIsolMap,
    const edm::Handle<reco::ElectronSeedCollection>   & L1IsoPixelSeedsMap,
    const edm::Handle<reco::ElectronIsolationMap>          & TrackEleIsolMap,
    EcalClusterLazyTools& lazyTools,
    const edm::ESHandle<MagneticField>& theMagField,
    reco::BeamSpot::Point & BSPosition )
{
  // if there are electrons, then the isolation maps and the SC should be in the event; if not it is an error

  if (recoIsolecalcands.isValid()) {
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoIsolecalcands->begin();
         recoecalcand!= recoIsolecalcands->end(); recoecalcand++) {
      // get the ref to the SC:
      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoIsolecalcands, distance(recoIsolecalcands->begin(), recoecalcand));
      reco::SuperClusterRef recrSC = ref->superCluster();
      //reco::SuperClusterRef recrSC = recoecalcand->superCluster();

      OpenHLTElectron ele;
      ele.hcalIsol   = -999;
      ele.trackIsol  = -999;
      ele.L1Isolated = true;
      ele.p          = -999;
      ele.pixelSeeds = -999;
      ele.newSC      = true;
      ele.clusterShape = -999;
      ele.Dphi = 700; 
      ele.Deta = 700;
      ele.Et         = recoecalcand->et();
      ele.eta        = recoecalcand->eta();
      ele.phi        = recoecalcand->phi();
      ele.E          = recrSC->energy();
      //Get the cluster shape
      std::vector<float> vCov = lazyTools.covariances( *(recrSC->seed()) );
      //      std::vector<float> vCov = lazyTools.localCovariances( *(recrSC->seed()) );
      double sigmaee = sqrt(vCov[0]);
      float EtaSC = fabs(recoecalcand->eta());
      if(EtaSC > 1.479 ) {//Endcap
        sigmaee = sigmaee - 0.02*(EtaSC - 2.3);
      }
      ele.clusterShape = sigmaee;

      // fill the hcal Isolation
      if (HcalEleIsolMap.isValid()) {
        //reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( reco::RecoEcalCandidateRef(recoIsolecalcands, distance(recoIsolecalcands->begin(), recoecalcand)) );
        reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( ref );
        if (mapi !=(*HcalEleIsolMap).end()) { ele.hcalIsol = mapi->val; }
      }
      // look if the SC has associated pixelSeeds
      int nmatch = 0;

      if (L1IsoPixelSeedsMap.isValid()) {
        for (reco::ElectronSeedCollection::const_iterator it = L1IsoPixelSeedsMap->begin();
             it != L1IsoPixelSeedsMap->end(); it++) {
	  edm::RefToBase<reco::CaloCluster> caloCluster = it->caloCluster() ;
 	  reco::SuperClusterRef scRef = caloCluster.castTo<reco::SuperClusterRef>() ;
          if (&(*recrSC) ==  &(*scRef)) { nmatch++; }
        }
      }

      ele.pixelSeeds = nmatch;

      // look if the SC was promoted to an electron:
      if (electronIsoHandle.isValid()) {
        bool FirstElectron = true;
        reco::ElectronRef electronref;
        for (reco::ElectronCollection::const_iterator iElectron = electronIsoHandle->begin();
             iElectron != electronIsoHandle->end(); iElectron++) {
          // 1) find the SC from the electron
          electronref = reco::ElectronRef(electronIsoHandle, iElectron - electronIsoHandle->begin());
          const reco::SuperClusterRef theClus = electronref->superCluster(); // SC from the electron;
          if (&(*recrSC) ==  &(*theClus)) {     // ref is the RecoEcalCandidateRef corresponding to the electron
            if (FirstElectron) {                // the first electron is stored in ele, keeping the ele.newSC = true
              FirstElectron = false;
              ele.p = electronref->track()->momentum().R();
	      float deta=-100, dphi=-100;
              CalculateDetaDphi(theMagField,BSPosition , electronref , deta, dphi);
              ele.Dphi=dphi; ele.Deta=deta;
              // fill the track Isolation
              if (TrackEleIsolMap.isValid()) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find(electronref);
                if (mapTr != (*TrackEleIsolMap).end()) { ele.trackIsol = mapTr->val; }
              }
            }
            else {
              // FirstElectron is false, i.e. the SC of this electron is common to another electron.
              // A new  OpenHLTElectron is inserted in the theHLTElectrons vector setting newSC = false
              OpenHLTElectron ele2;
              ele2.hcalIsol  = ele.hcalIsol;
              ele2.trackIsol = -999;
	      ele2.Dphi = 700; 
	      ele2.Deta = 700;
              ele2.Et  = ele.Et;
              ele2.eta = ele.eta;
              ele2.phi = ele.phi;
              ele2.E   = ele.E;
              ele2.L1Isolated = ele.L1Isolated;
              ele2.pixelSeeds = ele.pixelSeeds;
	      ele2.clusterShape = ele.clusterShape;
              ele2.newSC = false;
              ele2.p = electronref->track()->momentum().R();
	      float deta=-100, dphi=-100;
              CalculateDetaDphi(theMagField,BSPosition , electronref , deta, dphi);
              ele2.Dphi=dphi; ele2.Deta=deta;
              // fill the track Isolation
              if (TrackEleIsolMap.isValid()) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele2.trackIsol = mapTr->val;}
              }
              theHLTElectrons.push_back(ele2);
            }
          }
        } // end of loop over electrons
      } // end of if (electronIsoHandle) {

      //store the electron into the vector
      theHLTElectrons.push_back(ele);
    } // end of loop over ecalCandidates
  } // end of if (recoIsolecalcands) {
}


void HLTEgamma::MakeL1NonIsolatedElectrons(
    std::vector<OpenHLTElectron> & theHLTElectrons,
    const edm::Handle<reco::ElectronCollection>            & electronNonIsoHandle,
    const edm::Handle<reco::RecoEcalCandidateCollection>   & recoNonIsolecalcands,
    const edm::Handle<reco::RecoEcalCandidateIsolationMap> & HcalEleIsolMap,
    const edm::Handle<reco::ElectronSeedCollection>   & L1NonIsoPixelSeedsMap,
    const edm::Handle<reco::ElectronIsolationMap>          & TrackEleIsolMap,
    EcalClusterLazyTools& lazyTools,
    const edm::ESHandle<MagneticField>& theMagField,
    reco::BeamSpot::Point & BSPosition  )
{
  // if there are electrons, then the isolation maps and the SC should be in the event; if not it is an error

  if (recoNonIsolecalcands.isValid()) {
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoNonIsolecalcands->begin();
         recoecalcand!= recoNonIsolecalcands->end(); recoecalcand++) {
      //get the ref to the SC:
      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoNonIsolecalcands, distance(recoNonIsolecalcands->begin(), recoecalcand));
      reco::SuperClusterRef recrSC = ref->superCluster();
      //reco::SuperClusterRef recrSC = recoecalcand->superCluster();

      OpenHLTElectron ele;
      ele.hcalIsol   = -999;
      ele.trackIsol  = -999;
      ele.L1Isolated = false;
      ele.p          = -999;
      ele.pixelSeeds = -999;
      ele.newSC      = true;
      ele.clusterShape = -999;
      ele.Dphi = 700; 
      ele.Deta = 700;
      ele.Et         = recoecalcand->et();
      ele.eta        = recoecalcand->eta();
      ele.phi        = recoecalcand->phi();
      ele.E          = recrSC->energy();
      //Get the cluster shape
      std::vector<float> vCov = lazyTools.covariances( *(recrSC->seed()) );
      //      std::vector<float> vCov = lazyTools.localCovariances( *(recrSC->seed()) );
      double sigmaee = sqrt(vCov[0]);
      float EtaSC = fabs(recoecalcand->eta());
      if(EtaSC > 1.479 ) {//Endcap
        sigmaee = sigmaee - 0.02*(EtaSC - 2.3);
      }
      ele.clusterShape = sigmaee;

      // fill the hcal Isolation
      if (HcalEleIsolMap.isValid()) {
        // reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( reco::RecoEcalCandidateRef(recoNonIsolecalcands, distance(recoNonIsolecalcands->begin(), recoecalcand)) );
        reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( ref );
        if (mapi !=(*HcalEleIsolMap).end()) {ele.hcalIsol = mapi->val;}
      }
      // look if the SC has associated pixelSeeds
      int nmatch = 0;

      if (L1NonIsoPixelSeedsMap.isValid()) {
        for (reco::ElectronSeedCollection::const_iterator it = L1NonIsoPixelSeedsMap->begin();
             it != L1NonIsoPixelSeedsMap->end(); it++) {
	  edm::RefToBase<reco::CaloCluster> caloCluster = it->caloCluster() ;
 	  reco::SuperClusterRef scRef = caloCluster.castTo<reco::SuperClusterRef>() ;
          if (&(*recrSC) == &(*scRef)) { nmatch++;}
        }
      }

      ele.pixelSeeds = nmatch;

      // look if the SC was promoted to an electron:
      if (electronNonIsoHandle.isValid()) {
        bool FirstElectron = true;
        reco::ElectronRef electronref;
        for (reco::ElectronCollection::const_iterator iElectron = electronNonIsoHandle->begin(); 
             iElectron != electronNonIsoHandle->end();iElectron++) {
          // 1) find the SC from the electron
          electronref = reco::ElectronRef(electronNonIsoHandle, iElectron - electronNonIsoHandle->begin());
          const reco::SuperClusterRef theClus = electronref->superCluster(); //SC from the electron;
          if (&(*recrSC) ==  &(*theClus)) { // ref is the RecoEcalCandidateRef corresponding to the electron
            if (FirstElectron) { //the first electron is stored in ele, keeping the ele.newSC = true
              FirstElectron = false;
              ele.p = electronref->track()->momentum().R();
	      float deta=-100, dphi=-100;
              CalculateDetaDphi(theMagField,BSPosition , electronref , deta, dphi);
              ele.Dphi=dphi; ele.Deta=deta;

              // fill the track Isolation
              if (TrackEleIsolMap.isValid()) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele.trackIsol = mapTr->val;}
              }
            } else {
              // FirstElectron is false, i.e. the SC of this electron is common to another electron.
              // A new OpenHLTElectron is inserted in the theHLTElectrons vector setting newSC = false
              OpenHLTElectron ele2;
              ele2.hcalIsol   = ele.hcalIsol;
              ele2.trackIsol  =-999;
	      ele2.Dphi = 700; 
	      ele2.Deta = 700;
              ele2.Et         = ele.Et;
              ele2.eta        = ele.eta;
              ele2.phi        = ele.phi;
              ele2.E          = ele.E;
              ele2.L1Isolated = ele.L1Isolated;
              ele2.pixelSeeds = ele.pixelSeeds;
	      ele2.clusterShape = ele.clusterShape;
              ele2.newSC      = false;
              ele2.p          = electronref->track()->momentum().R();
	      float deta=-100, dphi=-100;
              CalculateDetaDphi(theMagField,BSPosition , electronref , deta, dphi);
              ele2.Dphi=dphi; ele2.Deta=deta;

              // fill the track Isolation
              if (TrackEleIsolMap.isValid()) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele2.trackIsol = mapTr->val;}
              }
              theHLTElectrons.push_back(ele2);
            }
          }
        } // end of loop over electrons
      } // end of if (electronNonIsoHandle) {

      // store the electron into the vector
      theHLTElectrons.push_back(ele);
    } // end of loop over ecalCandidates
  } // end of if (recoNonIsolecalcands) {
}

void HLTEgamma::CalculateDetaDphi(const edm::ESHandle<MagneticField>& theMagField, 
                                  reco::BeamSpot::Point & BSPosition, 
                                  const reco::ElectronRef eleref, 
                                  float& deltaeta, 
                                  float& deltaphi )
{
  
  const reco::SuperClusterRef theClus = eleref->superCluster();
  const math::XYZVector trackMom =  eleref->track()->momentum();

    math::XYZPoint SCcorrPosition(theClus->x()-BSPosition.x(), theClus->y()-BSPosition.y() , theClus->z()-eleref->track()->vz() );
    deltaeta = SCcorrPosition.eta()-eleref->track()->eta();

    ECALPositionCalculator posCalc;
    const math::XYZPoint vertex(BSPosition.x(),BSPosition.y(),eleref->track()->vz());

    float phi1= posCalc.ecalPhi(theMagField.product(),trackMom,vertex,1);
    float phi2= posCalc.ecalPhi(theMagField.product(),trackMom,vertex,-1);

    float deltaphi1=fabs( phi1 - theClus->position().phi() );
    if(deltaphi1>6.283185308) deltaphi1 -= 6.283185308;
    if(deltaphi1>3.141592654) deltaphi1 = 6.283185308-deltaphi1;

    float deltaphi2=fabs( phi2 - theClus->position().phi() );
    if(deltaphi2>6.283185308) deltaphi2 -= 6.283185308;
    if(deltaphi2>3.141592654) deltaphi2 = 6.283185308-deltaphi2;

    deltaphi = deltaphi1;
    if(deltaphi2<deltaphi1){ deltaphi = deltaphi2;}

}
