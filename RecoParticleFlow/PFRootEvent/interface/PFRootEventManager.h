#ifndef Demo_PFRootEvent_PFRootEventManager_h
#define Demo_PFRootEvent_PFRootEventManager_h

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPoint.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "RecoParticleFlow/PFAlgo/interface/PFBlock.h"

#include <TObject.h>
#include "TEllipse.h"
#include "TBox.h"

#include <string>
#include <map>
#include <set>
#include <vector>

class TTree;
class TBranch;
class TFile;
class TCanvas;
class TH2F;
class TGraph;
class IO;

class PFBlockElement;


/// \brief ROOT interface to particle flow package
/*!
  This base class allows to perform clustering and particle flow from 
  ROOT CINT (or any program). It is designed to support analysis and 
  developpement. Users should feel free to create their own PFRootEventManager,
  inheriting from this base class. Just reimplement the ProcessEntry function

  An example:

  \code
  gSystem->Load("libFWCoreFWLite.so");
  gSystem->Load("libRecoParticleFlowPFRootEvent.so");
  AutoLibraryLoader::enable();
  gSystem->Load("libCintex.so");
  ROOT::Cintex::Cintex::Enable();

  PFRootEventManager em("pfRootEvent.opt");
  int i=0;
  em.display(i++);
  \endcode
  
  pfRootEvent.opt is an option file (see IO class):
  \verbatim
  root file test.root

  root hits_branch  recoPFRecHits_pfcluster__Demo.obj
  root recTracks_branch  recoPFRecTracks_pf_PFRecTrackCollection_Demo.obj

  display algos 1 

  display  viewsize_etaphi 600 400
  display  viewsize_xy     400 400

  display  color_clusters		1

  clustering thresh_Ecal_Barrel           0.2
  clustering thresh_Seed_Ecal_Barrel      0.3
  clustering thresh_Ecal_Endcap           0.2
  clustering thresh_Seed_Ecal_Endcap      0.9
  clustering neighbours_Ecal		4

  clustering depthCor_Mode          1
  clustering depthCor_A 		  0.89
  clustering depthCor_B 		  7.3
  clustering depthCor_A_preshower   0.89
  clustering depthCor_B_preshower   4.0

  clustering thresh_Hcal_Barrel           1.0
  clustering thresh_Seed_Hcal_Barrel      1.4
  clustering thresh_Hcal_Endcap           1.0
  clustering thresh_Seed_Hcal_Endcap      1.4
  clustering neighbours_Hcal		4
  \endverbatim
  
  \author Colin Bernet, Renaud Bruneliere
  \date July 2006
*/
class PFRootEventManager {

 public:

  /// viewport definition
  enum View_t { XY = 0, RZ = 1, EPE = 2, EPH = 3, NViews = 4 };

  /// default constructor
  PFRootEventManager();

  /// \param is an option file, see IO
  PFRootEventManager(const char* file);

  /// destructor
  virtual ~PFRootEventManager();
  
  /// reset before next event
  void reset();

  /// parse option file
  void readOptions(const char* file, bool refresh=true);

  /// process one entry 
  virtual bool processEntry(int entry);

  /// performs clustering 
  void clustering();

  /// performs particle flow
  void particleFlow();

  /// display one entry 
  void display(int ientry);

  /// display x/y or r/z
  void displayView(unsigned viewType);

  /// display reconstructed calorimeter hits in x/y or r/z view
  void displayRecHits(unsigned viewType, double phi0 = 0.);

  /// display a reconstructed calorimeter hit in x/y or r/z view
  void displayRecHit(reco::PFRecHit& rh, unsigned viewType,
		     double maxe, double thresh, double phi0 = 0.);

  /// display clusters in x/y or r/z view
  void displayClusters(unsigned viewType, double phi0 = 0.);

  /// display reconstructed tracks in x/y or r/z view
  void displayRecTracks(unsigned viewType, double phi0 = 0.);

  /// unzooms all support histograms
  void unZoom();

  /// updates all displays
  void updateDisplay();

  /// look for ecal
  void lookForMaxRecHit(bool ecal);

  /// look for hcal 

  /// finds max rechit energy in a given layer 
  double getMaxE(int layer) const;

  /// max rechit energy in ecal 
  double getMaxEEcal();

  /// max rechit energy in hcal 
  double getMaxEHcal();

  /// current event
  int        iEvent_;
  
  /// options file parser 
  IO*        options_;      
  
  /// input tree  
  TTree*     tree_;          
  
  /// cluster branch  
  TBranch*   hitsBranch_;          
  
  /// reconstructed tracks branch  
  TBranch*   recTracksBranch_;          
  
  /// rechits
  std::vector<reco::PFRecHit> rechits_;

  /// clusters
  std::vector<reco::PFCluster> clusters_;

  /// reconstructed tracks
  std::vector<reco::PFRecTrack> recTracks_;
  
  /// all pfblock elements
  std::set< PFBlockElement* > allElements_; 

  /// all pfblocks  
  std::vector< PFBlock >     allPFBs_;

  /// input file
  TFile*     file_; 

  /// input file name
  std::string     inFileName_;   

  /// output file
  TFile*     outFile_;   

  /// output filename
  std::string     outFileName_;   


  /// canvases for eta/phi display, one per algo
  /// each is split in 2 : HCAL, ECAL
  std::map<int, TCanvas* > displayEtaPhi_;        

  /// algos to display
  std::set<int>            algosToDisplay_;  

  /// display pad xy size for eta/phi view
  std::vector<int>         viewSizeEtaPhi_;        

  /// display cluster color ? (then color = cluster type )
  bool displayColorClusters_;

  /// size of view in number of cells when centering on a rechit
  double displayZoomFactor_;

  /// support histogram for eta/phi display
  TH2F*                    displayHistEtaPhi_;

  /// vector of canvas for x/y or r/z display
  std::vector<TCanvas*> displayView_;

  /// display pad xy size for (x,y) or (r,z) display
  std::vector<int>      viewSize_;     
 
  /// support histogram for x/y or r/z display. COLIN: why a vector ?
  std::vector<TH2F*>    displayHist_;

  /// ECAL in XY view. COLIN: should be attribute ?
  TEllipse frontFaceECALXY_;

  /// ECAL in RZ view. COLIN: should be attribute ?
  TBox     frontFaceECALRZ_;

  /// HCAL in XY view. COLIN: should be attribute ?
  TEllipse frontFaceHCALXY_;

  /// vector of TGraph used to represent the track in XY or RZ view. COLIN: should be attribute ?  
  std::vector< std::vector<TGraph*> > graphTrack_;

  /// max rechit energy in ecal
  double                   maxERecHitEcal_;

  /// max rechit energy in hcal
  double                   maxERecHitHcal_;


  //----------------- clustering parameters ---------------------

  /// ecal barrel threshold
  double threshEcalBarrel_;

  /// ecal barrel seed threshold
  double threshSeedEcalBarrel_;

  /// ecal endcap threshold
  double threshEcalEndcap_;

  /// ecal endcap seed threshold
  double threshSeedEcalEndcap_;

  /// ecal number of neighbours
  int    nNeighboursEcal_;


  /// hcal barrel threshold
  double threshHcalBarrel_;

  /// hcal barrel seed threshold
  double threshSeedHcalBarrel_;

  /// hcal endcap threshold
  double threshHcalEndcap_;

  /// hcal endcap seed threshold
  double threshSeedHcalEndcap_;

  /// hcal number of neighbours
  int    nNeighboursHcal_;


  // particle flow ------------------------------------------
  bool   displayJetColors_;
  int    reconMethod_;
};
#endif
