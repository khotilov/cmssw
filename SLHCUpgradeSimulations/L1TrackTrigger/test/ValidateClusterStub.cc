//////////////////////////
//  Analyzer by Nicola  //
//    Oct 2011 @ PD     //
//////////////////////////

/////////////////////////
//       HEADERS       //
/////////////////////////

////////////////
// CLASS HEADER
// No more necessary in the current "no *.h file" implementation

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"
//
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/Common/interface/DetSetVector.h"
//
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
//
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing.h"
////////////////////////
// FAST SIMULATION STUFF
#include "FastSimulation/Particle/interface/RawParticle.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
//
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"
#include "DataFormats/SiPixelDetId/interface/StackedTrackerDetId.h"

////////////////
// PHYSICS TOOLS
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "RecoTauTag/TauTagTools/interface/GeneratorTau.h"
//
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementVector.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"

///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH1.h>

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class TTree;
class TFile;
class TH1D;
class TH2D;
class TGraph;
class RectangularPixelTopology;
class TransientInitialStateEstimator;
class MagneticField;
class TrackerGeometry;
class TrajectoryStateOnSurface;
class PTrajectoryStateOnDet;
//
class ValidateClusterStub : public edm::EDAnalyzer
{
  /// Public methods
  public:
    /// Constructor/destructor
    explicit ValidateClusterStub(const edm::ParameterSet& iConfig);
    virtual ~ValidateClusterStub();
    // Typical methods used on Loops over events
    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

    /// Some Type definitions

  /// Protected methods only internally used
  protected:
                     
  /// Private methods and variables
  private:

    /// SimTrack and SimVertex
    TH2D* hSimVtx_XY;
    TH2D* hSimVtx_RZ;

    TH1D* hSimTrk_Pt;
    TH1D* hSimTrk_Eta_Pt10;
    TH1D* hSimTrk_Phi_Pt10;

    TH1D* hSimTrk_Pt_Aux;
    TH1D* hSimTrk_Eta_Aux;
    TH1D* hSimTrk_Phi_Aux;

    /// Global positions of L1TkClusters
    TH2D* hCluster_Barrel_XY;
    TH2D* hCluster_Barrel_XY_Zoom;
    TH2D* hCluster_Endcap_Fw_XY;
    TH2D* hCluster_Endcap_Bw_XY;
    TH2D* hCluster_RZ;
    TH2D* hCluster_Endcap_Fw_RZ_Zoom;
    TH2D* hCluster_Endcap_Bw_RZ_Zoom;

    TH1D* hCluster_Mem;
    TH1D* hCluster_Barrel_Stack;
    TH1D* hCluster_Endcap_Stack;
    TH2D* hCluster_Gen; /// For both Stack Members
    TH2D* hCluster_PID;
    TH2D* hCluster_W;
    
    /// Global positions of L1TkStubs
    TH2D* hStub_Barrel_XY;
    TH2D* hStub_Barrel_XY_Zoom;
    TH2D* hStub_Endcap_Fw_XY;
    TH2D* hStub_Endcap_Bw_XY;
    TH2D* hStub_RZ;
    TH2D* hStub_Endcap_Fw_RZ_Zoom;
    TH2D* hStub_Endcap_Bw_RZ_Zoom;

    TH1D* hStub_Barrel_Stack;
    TH1D* hStub_Endcap_Stack;
    TH1D* hStub_Gen;
    TH1D* hStub_PID;
    TH2D* hStub_Barrel_W;
    TH2D* hStub_Barrel_O;
    TH2D* hStub_Endcap_W;
    TH2D* hStub_Endcap_O;


    /// Denominator for Stub Prod Eff
    std::map< unsigned int, TH1D* > mapCluLayer_hSimTrk_Pt;
    std::map< unsigned int, TH1D* > mapCluLayer_hSimTrk_Eta_Pt10;
    std::map< unsigned int, TH1D* > mapCluLayer_hSimTrk_Phi_Pt10;
    std::map< unsigned int, TH1D* > mapCluDisk_hSimTrk_Pt;
    std::map< unsigned int, TH1D* > mapCluDisk_hSimTrk_Eta_Pt10;
    std::map< unsigned int, TH1D* > mapCluDisk_hSimTrk_Phi_Pt10;
    /// Numerator for Stub Prod Eff
    std::map< unsigned int, TH1D* > mapStubLayer_hSimTrk_Pt;
    std::map< unsigned int, TH1D* > mapStubLayer_hSimTrk_Eta_Pt10;
    std::map< unsigned int, TH1D* > mapStubLayer_hSimTrk_Phi_Pt10;
    std::map< unsigned int, TH1D* > mapStubDisk_hSimTrk_Pt;
    std::map< unsigned int, TH1D* > mapStubDisk_hSimTrk_Eta_Pt10;
    std::map< unsigned int, TH1D* > mapStubDisk_hSimTrk_Phi_Pt10;

    /// Comparison of Stubs to SimTracks
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_InvPt_SimTrk_InvPt;
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_Pt_SimTrk_Pt;
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_Eta_SimTrk_Eta;
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_Phi_SimTrk_Phi;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_InvPt_SimTrk_InvPt;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_Pt_SimTrk_Pt;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_Eta_SimTrk_Eta;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_Phi_SimTrk_Phi;

    /// Residuals
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_InvPtRes_SimTrk_Eta;
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_PtRes_SimTrk_Eta;
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_EtaRes_SimTrk_Eta;
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_PhiRes_SimTrk_Eta;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_InvPtRes_SimTrk_Eta;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_PtRes_SimTrk_Eta;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_EtaRes_SimTrk_Eta;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_PhiRes_SimTrk_Eta;

    /// Stub Width vs Pt
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_W_SimTrk_Pt;
    std::map< unsigned int, TH2D* > mapStubLayer_hStub_W_SimTrk_InvPt;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_W_SimTrk_Pt;
    std::map< unsigned int, TH2D* > mapStubDisk_hStub_W_SimTrk_InvPt;

    /// Containers of parameters passed by python
    /// configuration file
    edm::ParameterSet config;

    bool testedGeometry;
    bool DebugMode;
};

//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
ValidateClusterStub::ValidateClusterStub(edm::ParameterSet const& iConfig) : 
  config(iConfig)
{
  /// Insert here what you need to initialize
  DebugMode = iConfig.getParameter< bool >("DebugMode");
}

/////////////
// DESTRUCTOR
ValidateClusterStub::~ValidateClusterStub()
{
  /// Insert here what you need to delete
  /// when you close the class instance
}  

//////////
// END JOB
void ValidateClusterStub::endJob()//edm::Run& run, const edm::EventSetup& iSetup
{
  /// Things to be done at the exit of the event Loop
  std::cerr << " ValidateClusterStub::endJob" << std::endl;
  /// End of things to be done at the exit from the event Loop
}

////////////
// BEGIN JOB
void ValidateClusterStub::beginJob()
{
  /// Initialize all slave variables
  /// mainly histogram ranges and resolution
  testedGeometry = false;

  std::ostringstream histoName;
  std::ostringstream histoTitle;

  /// Things to be done before entering the event Loop
  std::cerr << " ValidateClusterStub::beginJob" << std::endl;

  /// Book histograms etc
  edm::Service<TFileService> fs;

  /// Prepare for LogXY Plots
  int NumBins = 200;
  double MinPt = 0.0;
  double MaxPt = 100.0;

  double* BinVec = new double[NumBins+1];
  for ( int iBin = 0; iBin < NumBins + 1; iBin++ )
  {
    double temp = pow( 10, (- NumBins + iBin)/(MaxPt - MinPt)  );
    BinVec[ iBin ] = temp;
  }

  /// SimTrack and SimVertex
  hSimVtx_XY      = fs->make<TH2D>( "hSimVtx_XY",       "SimVtx y vs. x",           140, -0.7, 0.7, 140, -0.7, 0.7 );
  hSimVtx_XY->Sumw2();

  hSimVtx_RZ      = fs->make<TH2D>( "hSimVtx_RZ",       "SimVtx #rho vs. z",           200, -50, 50, 140, 0, 0.7 );
  hSimVtx_RZ->Sumw2();

  hSimTrk_Pt       = fs->make<TH1D>( "hSimTrk_Pt",         "SimTrk p_{T}",                    100, 0, 50 );
  hSimTrk_Eta_Pt10 = fs->make<TH1D>( "hSimTrk_Eta_Pt10",   "SimTrk #eta (p_{T} > 10 GeV/c)",  180, -M_PI, M_PI );
  hSimTrk_Phi_Pt10 = fs->make<TH1D>( "hSimTrk_Phi_Pt10",   "SimTrk #phi (p_{T} > 10 GeV/c)",  180, -M_PI, M_PI );
  hSimTrk_Pt->Sumw2();
  hSimTrk_Eta_Pt10->Sumw2();
  hSimTrk_Phi_Pt10->Sumw2();

  hSimTrk_Pt_Aux  = fs->make<TH1D>( "hSimTrk_Pt_Aux",  "SimTrk p_{T}, Aux", 100, 0, 50 );
  hSimTrk_Eta_Aux = fs->make<TH1D>( "hSimTrk_Eta_Aux", "SimTrk #eta, Aux",  180, -M_PI, M_PI );
  hSimTrk_Phi_Aux = fs->make<TH1D>( "hSimTrk_Phi_Aux", "SimTrk #phi, Aux",  180, -M_PI, M_PI );
  hSimTrk_Pt_Aux->Sumw2();
  hSimTrk_Eta_Aux->Sumw2();
  hSimTrk_Phi_Aux->Sumw2();

  /// Global position of L1TkCluster
  hCluster_Barrel_XY          = fs->make<TH2D>( "hCluster_Barrel_XY",         "L1TkCluster Barrel y vs. x",              960, -120, 120, 960, -120, 120 );
  hCluster_Barrel_XY_Zoom     = fs->make<TH2D>( "hCluster_Barrel_XY_Zoom",    "L1TkCluster Barrel y vs. x",              960, 30, 60, 960, -15, 15 );
  hCluster_Endcap_Fw_XY       = fs->make<TH2D>( "hCluster_Endcap_Fw_XY",      "L1TkCluster Forward Endcap y vs. x",      960, -120, 120, 960, -120, 120 );
  hCluster_Endcap_Bw_XY       = fs->make<TH2D>( "hCluster_Endcap_Bw_XY",      "L1TkCluster Backward Endcap y vs. x",     960, -120, 120, 960, -120, 120 );
  hCluster_RZ                 = fs->make<TH2D>( "hCluster_RZ",                "L1TkCluster #rho vs. z",                  900, -300, 300, 480, 0, 120 );
  hCluster_Endcap_Fw_RZ_Zoom  = fs->make<TH2D>( "hCluster_Endcap_Fw_RZ_Zoom", "L1TkCluster Forward Endcap #rho vs. z",   960, 140, 170, 960, 30, 60 );
  hCluster_Endcap_Bw_RZ_Zoom  = fs->make<TH2D>( "hCluster_Endcap_Bw_RZ_Zoom", "L1TkCluster Backward Endcap #rho vs. z",  960, -170, -140, 960, 70, 100 );
  hCluster_Barrel_XY->Sumw2();
  hCluster_Barrel_XY_Zoom->Sumw2();
  hCluster_Endcap_Fw_XY->Sumw2();
  hCluster_Endcap_Bw_XY->Sumw2();
  hCluster_RZ->Sumw2();
  hCluster_Endcap_Fw_RZ_Zoom->Sumw2();
  hCluster_Endcap_Bw_RZ_Zoom->Sumw2();

  hCluster_Mem   = fs->make<TH1D>("hCluster_Mem", "L1TkCluster Stack Member",      2, -0.5, 1.5 );
  hCluster_Barrel_Stack = fs->make<TH1D>("hCluster_Barrel_Stack", "L1TkCluster Stack",           12, -0.5, 11.5 );
  hCluster_Endcap_Stack = fs->make<TH1D>("hCluster_Endcap_Stack", "L1TkCluster Stack",           12, -0.5, 11.5 );
  hCluster_Gen   = fs->make<TH2D>("hCluster_Gen", "L1TkCluster Genuine (Member)",  2, -0.5, 1.5, 2, -0.5, 1.5 );
  hCluster_PID   = fs->make<TH2D>("hCluster_PID", "L1TkCluster PID (Member)",      501, -250.5, 250.5, 2, -0.5, 1.5 );
  hCluster_W     = fs->make<TH2D>("hCluster_W", "L1TkCluster Width (Member)",      10, -0.5, 9.5, 2, -0.5, 1.5 );
  hCluster_Mem->Sumw2();
  hCluster_Barrel_Stack->Sumw2();
  hCluster_Endcap_Stack->Sumw2();
  hCluster_Gen->Sumw2();
  hCluster_PID->Sumw2();
  hCluster_W->Sumw2();

  /// Global position of L1TkStub
  hStub_Barrel_XY          = fs->make<TH2D>( "hStub_Barrel_XY",         "L1TkStub Barrel y vs. x",              960, -120, 120, 960, -120, 120 );
  hStub_Barrel_XY_Zoom     = fs->make<TH2D>( "hStub_Barrel_XY_Zoom",    "L1TkStub Barrel y vs. x",              960, 30, 60, 960, -15, 15 );
  hStub_Endcap_Fw_XY       = fs->make<TH2D>( "hStub_Endcap_Fw_XY",      "L1TkStub Forward Endcap y vs. x",      960, -120, 120, 960, -120, 120 );
  hStub_Endcap_Bw_XY       = fs->make<TH2D>( "hStub_Endcap_Bw_XY",      "L1TkStub Backward Endcap y vs. x",     960, -120, 120, 960, -120, 120 );
  hStub_RZ                 = fs->make<TH2D>( "hStub_RZ",                "L1TkStub #rho vs. z",                  900, -300, 300, 480, 0, 120 );
  hStub_Endcap_Fw_RZ_Zoom  = fs->make<TH2D>( "hStub_Endcap_Fw_RZ_Zoom", "L1TkStub Forward Endcap #rho vs. z",   960, 140, 170, 960, 30, 60 );
  hStub_Endcap_Bw_RZ_Zoom  = fs->make<TH2D>( "hStub_Endcap_Bw_RZ_Zoom", "L1TkStub Backward Endcap #rho vs. z",  960, -170, -140, 960, 70, 100 );
  hStub_Barrel_XY->Sumw2();
  hStub_Barrel_XY_Zoom->Sumw2();
  hStub_Endcap_Fw_XY->Sumw2();
  hStub_Endcap_Bw_XY->Sumw2();
  hStub_RZ->Sumw2();
  hStub_Endcap_Fw_RZ_Zoom->Sumw2();
  hStub_Endcap_Bw_RZ_Zoom->Sumw2();

  hStub_Barrel_Stack     = fs->make<TH1D>("hStub_Barrel_Stack",    "L1TkStub Stack",                          12, -0.5, 11.5 );
  hStub_Endcap_Stack     = fs->make<TH1D>("hStub_Endcap_Stack",    "L1TkStub Stack",                          12, -0.5, 11.5 );
  hStub_Gen       = fs->make<TH1D>("hStub_Gen",      "L1TkStub Genuine",                        2, -0.5, 1.5 );
  hStub_PID       = fs->make<TH1D>("hStub_PID",      "L1TkStub PID",                            501, -250.5, 250.5 );
  hStub_Barrel_W  = fs->make<TH2D>("hStub_Barrel_W", "L1TkStub Post-Corr Displacement (Layer)", 21, -10.5, 10.5, 12, -0.5, 11.5 );
  hStub_Barrel_O  = fs->make<TH2D>("hStub_Barrel_O", "L1TkStub Offset (Layer)",                 21, -10.5, 10.5, 12, -0.5, 11.5 );
  hStub_Endcap_W  = fs->make<TH2D>("hStub_Endcap_W", "L1TkStub Post-Corr Displacement (Layer)", 21, -10.5, 10.5, 12, -0.5, 11.5 );
  hStub_Endcap_O  = fs->make<TH2D>("hStub_Endcap_O", "L1TkStub Offset (Layer)",                 21, -10.5, 10.5, 12, -0.5, 11.5 );
  hStub_Barrel_Stack->Sumw2();
  hStub_Endcap_Stack->Sumw2();
  hStub_Gen->Sumw2();
  hStub_PID->Sumw2();
  hStub_Barrel_W->Sumw2();
  hStub_Barrel_O->Sumw2();
  hStub_Endcap_W->Sumw2();
  hStub_Endcap_O->Sumw2();

  /// Stub Production Efficiency and comparison to SimTrack
  for ( unsigned int stackIdx = 0; stackIdx < 12; stackIdx++ )
  {
    /// BARREL

    /// Denominators
    histoName.str("");  histoName << "hSimTrk_Pt_Clu_L" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk p_{T}, Cluster, Barrel Stack " << stackIdx;
    mapCluLayer_hSimTrk_Pt[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                         100, 0, 50 );
    histoName.str("");  histoName << "hSimTrk_Eta_Pt10_Clu_L" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk #eta (p_{T} > 10 GeV/c), Cluster, Barrel Stack " << stackIdx;
    mapCluLayer_hSimTrk_Eta_Pt10[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                               180, -M_PI, M_PI );
    histoName.str("");  histoName << "hSimTrk_Phi_Pt10_Clu_L" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk #phi (p_{T} > 10 GeV/c), Cluster, Barrel Stack " << stackIdx;
    mapCluLayer_hSimTrk_Phi_Pt10[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                               180, -M_PI, M_PI );
    mapCluLayer_hSimTrk_Pt[ stackIdx ]->Sumw2();
    mapCluLayer_hSimTrk_Eta_Pt10[ stackIdx ]->Sumw2();
    mapCluLayer_hSimTrk_Phi_Pt10[ stackIdx ]->Sumw2();

    /// Numerators GeV/c
    histoName.str("");  histoName << "hSimTrk_Pt_Stub_L" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk p_{T}, Stub, Barrel Stack " << stackIdx;
    mapStubLayer_hSimTrk_Pt[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                          100, 0, 50 );
    histoName.str("");  histoName << "hSimTrk_Eta_Pt10_Stub_L" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk #eta (p_{T} > 10 GeV/c), Stub, Barrel Stack " << stackIdx;
    mapStubLayer_hSimTrk_Eta_Pt10[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                180, -M_PI, M_PI );
    histoName.str("");  histoName << "hSimTrk_Phi_Pt10_Stub_L" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk #phi (p_{T} > 10 GeV/c), Stub, Barrel Stack " << stackIdx;
    mapStubLayer_hSimTrk_Phi_Pt10[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                180, -M_PI, M_PI );
    mapStubLayer_hSimTrk_Pt[ stackIdx ]->Sumw2();
    mapStubLayer_hSimTrk_Eta_Pt10[ stackIdx ]->Sumw2();
    mapStubLayer_hSimTrk_Phi_Pt10[ stackIdx ]->Sumw2();

    /// Comparison to SimTrack
    histoName.str("");  histoName << "hStub_InvPt_SimTrk_InvPt_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub p_{T}^{-1} vs. SimTrk p_{T}^{-1}, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_InvPt_SimTrk_InvPt[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                        200, 0.0, 0.8,
                                                                        200, 0.0, 0.8 );
    mapStubLayer_hStub_InvPt_SimTrk_InvPt[ stackIdx ]->GetXaxis()->Set( NumBins, BinVec );
    mapStubLayer_hStub_InvPt_SimTrk_InvPt[ stackIdx ]->GetYaxis()->Set( NumBins, BinVec );
    mapStubLayer_hStub_InvPt_SimTrk_InvPt[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_Pt_SimTrk_Pt_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub p_{T} vs. SimTrk p_{T}, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_Pt_SimTrk_Pt[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                  100, 0, 50,
                                                                  100, 0, 50 );
    mapStubLayer_hStub_Pt_SimTrk_Pt[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_Eta_SimTrk_Eta_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub #eta vs. SimTrk #eta, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_Eta_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                    180, -M_PI, M_PI,
                                                                    180, -M_PI, M_PI );
    mapStubLayer_hStub_Eta_SimTrk_Eta[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_Phi_SimTrk_Phi_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub #phi vs. SimTrk #phi, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_Phi_SimTrk_Phi[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                    180, -M_PI, M_PI,
                                                                    180, -M_PI, M_PI );
    mapStubLayer_hStub_Phi_SimTrk_Phi[ stackIdx ]->Sumw2();

    /// Residuals
    histoName.str("");  histoName << "hStub_InvPtRes_SimTrk_Eta_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub p_{T}^{-1} - SimTrk p_{T}^{-1} vs. SimTrk #eta, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_InvPtRes_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                         180, -M_PI, M_PI,
                                                                         100, -2.0, 2.0 );
    mapStubLayer_hStub_InvPtRes_SimTrk_Eta[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_PtRes_SimTrk_Eta_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub p_{T} - SimTrk p_{T} vs. SimTrk #eta, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_PtRes_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                      180, -M_PI, M_PI,
                                                                      100, -40, 40 );
    mapStubLayer_hStub_PtRes_SimTrk_Eta[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_EtaRes_SimTrk_Eta_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub #eta - SimTrk #eta vs. SimTrk #eta, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_EtaRes_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                       180, -M_PI, M_PI,
                                                                       100, -2, 2 );
    mapStubLayer_hStub_EtaRes_SimTrk_Eta[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_PhiRes_SimTrk_Eta_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub #phi - SimTrk #phi vs. SimTrk #eta, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_PhiRes_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                       180, -M_PI, M_PI,
                                                                       100, -0.5, 0.5 );
    mapStubLayer_hStub_PhiRes_SimTrk_Eta[ stackIdx ]->Sumw2();

    /// Stub Width vs. Pt
    histoName.str("");  histoName << "hStub_W_SimTrk_Pt_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub Width vs. SimTrk p_{T}, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_W_SimTrk_Pt[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                 200, 0, 50,
                                                                 41, -10.25, 10.25 );
    mapStubLayer_hStub_W_SimTrk_Pt[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_W_SimTrk_InvPt_L" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub Width vs. SimTrk p_{T}^{-1}, Barrel Stack " << stackIdx;
    mapStubLayer_hStub_W_SimTrk_InvPt[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                    200, 0, 0.8,
                                                                    41, -10.25, 10.25 );
    mapStubLayer_hStub_W_SimTrk_InvPt[ stackIdx ]->GetXaxis()->Set( NumBins, BinVec );
    mapStubLayer_hStub_W_SimTrk_InvPt[ stackIdx ]->Sumw2();

    /// ENDCAP

    /// Denominators
    histoName.str("");  histoName << "hSimTrk_Pt_Clu_D" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk p_{T}, Cluster, Endcap Stack " << stackIdx;
    mapCluDisk_hSimTrk_Pt[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                        100, 0, 50 );
    histoName.str("");  histoName << "hSimTrk_Eta_Pt10_Clu_D" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk #eta (p_{T} > 10 GeV/c), Cluster, Endcap Stack " << stackIdx;
    mapCluDisk_hSimTrk_Eta_Pt10[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                              180, -M_PI, M_PI );
    histoName.str("");  histoName << "hSimTrk_Phi_Pt10_Clu_D" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk #phi (p_{T} > 10 GeV/c), Cluster, Endcap Stack " << stackIdx;
    mapCluDisk_hSimTrk_Phi_Pt10[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                              180, -M_PI, M_PI );
    mapCluDisk_hSimTrk_Pt[ stackIdx ]->Sumw2();
    mapCluDisk_hSimTrk_Eta_Pt10[ stackIdx ]->Sumw2();
    mapCluDisk_hSimTrk_Phi_Pt10[ stackIdx ]->Sumw2();

    /// Numerators GeV/c
    histoName.str("");  histoName << "hSimTrk_Pt_Stub_D" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk p_{T}, Stub, Endcap Stack " << stackIdx;
    mapStubDisk_hSimTrk_Pt[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                         100, 0, 50 );
    histoName.str("");  histoName << "hSimTrk_Eta_Pt10_Stub_D" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk #eta (p_{T} > 10 GeV/c), Stub, Endcap Stack " << stackIdx;
    mapStubDisk_hSimTrk_Eta_Pt10[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                               180, -M_PI, M_PI );
    histoName.str("");  histoName << "hSimTrk_Phi_Pt10_Stub_D" << stackIdx;
    histoTitle.str(""); histoTitle << "SimTrk #phi (p_{T} > 10 GeV/c), Stub, Endcap Stack " << stackIdx;
    mapStubDisk_hSimTrk_Phi_Pt10[ stackIdx ] = fs->make<TH1D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                               180, -M_PI, M_PI );
    mapStubDisk_hSimTrk_Pt[ stackIdx ]->Sumw2();
    mapStubDisk_hSimTrk_Eta_Pt10[ stackIdx ]->Sumw2();
    mapStubDisk_hSimTrk_Phi_Pt10[ stackIdx ]->Sumw2();

    /// Comparison to SimTrack
    histoName.str("");  histoName << "hStub_InvPt_SimTrk_InvPt_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub p_{T}^{-1} vs. SimTrk p_{T}^{-1}, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_InvPt_SimTrk_InvPt[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                       200, 0.0, 0.8,
                                                                       200, 0.0, 0.8 );
    mapStubDisk_hStub_InvPt_SimTrk_InvPt[ stackIdx ]->GetXaxis()->Set( NumBins, BinVec );
    mapStubDisk_hStub_InvPt_SimTrk_InvPt[ stackIdx ]->GetYaxis()->Set( NumBins, BinVec );
    mapStubDisk_hStub_InvPt_SimTrk_InvPt[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_Pt_SimTrk_Pt_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub p_{T} vs. SimTrk p_{T}, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_Pt_SimTrk_Pt[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                 100, 0, 50,
                                                                 100, 0, 50 );
    mapStubDisk_hStub_Pt_SimTrk_Pt[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_Eta_SimTrk_Eta_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub #eta vs. SimTrk #eta, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_Eta_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                   180, -M_PI, M_PI,
                                                                   180, -M_PI, M_PI );
    mapStubDisk_hStub_Eta_SimTrk_Eta[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_Phi_SimTrk_Phi_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub #phi vs. SimTrk #phi, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_Phi_SimTrk_Phi[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                   180, -M_PI, M_PI,
                                                                   180, -M_PI, M_PI );
    mapStubDisk_hStub_Phi_SimTrk_Phi[ stackIdx ]->Sumw2();

    /// Residuals
    histoName.str("");  histoName << "hStub_InvPtRes_SimTrk_Eta_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub p_{T}^{-1} - SimTrk p_{T}^{-1} vs. SimTrk #eta, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_InvPtRes_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                        180, -M_PI, M_PI,
                                                                        100, -2.0, 2.0 );
    mapStubDisk_hStub_InvPtRes_SimTrk_Eta[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_PtRes_SimTrk_Eta_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub p_{T} - SimTrk p_{T} vs. SimTrk #eta, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_PtRes_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                     180, -M_PI, M_PI,
                                                                     100, -40, 40 );
    mapStubDisk_hStub_PtRes_SimTrk_Eta[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_EtaRes_SimTrk_Eta_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub #eta - SimTrk #eta vs. SimTrk #eta, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_EtaRes_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                      180, -M_PI, M_PI,
                                                                      100, -2, 2 );
    mapStubDisk_hStub_EtaRes_SimTrk_Eta[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_PhiRes_SimTrk_Eta_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub #phi - SimTrk #phi vs. SimTrk #eta, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_PhiRes_SimTrk_Eta[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                      180, -M_PI, M_PI,
                                                                      100, -0.5, 0.5 );
    mapStubDisk_hStub_PhiRes_SimTrk_Eta[ stackIdx ]->Sumw2();

    /// Stub Width vs. Pt
    histoName.str("");  histoName << "hStub_W_SimTrk_Pt_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub Width vs. SimTrk p_{T}, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_W_SimTrk_Pt[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                200, 0, 50,
                                                                41, -10.25, 10.25 );
    mapStubDisk_hStub_W_SimTrk_Pt[ stackIdx ]->Sumw2();

    histoName.str("");  histoName << "hStub_W_SimTrk_InvPt_D" << stackIdx;
    histoTitle.str(""); histoTitle << "Stub Width vs. SimTrk p_{T}^{-1}, Endcap Stack " << stackIdx;
    mapStubDisk_hStub_W_SimTrk_InvPt[ stackIdx ] = fs->make<TH2D>( histoName.str().c_str(),  histoTitle.str().c_str(),
                                                                   200, 0, 0.8,
                                                                   41, -10.25, 10.25 );
    mapStubDisk_hStub_W_SimTrk_InvPt[ stackIdx ]->GetXaxis()->Set( NumBins, BinVec );
    mapStubDisk_hStub_W_SimTrk_InvPt[ stackIdx ]->Sumw2();
  }

  /// End of things to be done before entering the event Loop
}

//////////
// ANALYZE
void ValidateClusterStub::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /// Geometry handles etc
  edm::ESHandle< TrackerGeometry >                             GeometryHandle;
  edm::ESHandle< StackedTrackerGeometry >         StackedGeometryHandle;
  const StackedTrackerGeometry*                   theStackedGeometry;
  StackedTrackerGeometry::StackContainerIterator  StackedTrackerIterator;

  /// Geometry setup
  /// Set pointers to Geometry
  iSetup.get< TrackerDigiGeometryRecord >().get(GeometryHandle);
  /// Set pointers to Stacked Modules
  iSetup.get< StackedTrackerGeometryRecord >().get(StackedGeometryHandle);
  theStackedGeometry = StackedGeometryHandle.product(); /// Note this is different 
                                                        /// from the "global" geometry

  /// Magnetic Field
  edm::ESHandle< MagneticField > magneticFieldHandle;
  iSetup.get< IdealMagneticFieldRecord >().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();

  /// Sim Tracks and Vtx
  edm::Handle< edm::SimTrackContainer >  SimTrackHandle;
  edm::Handle< edm::SimVertexContainer > SimVtxHandle;
  //iEvent.getByLabel( "famosSimHits", SimTrackHandle );
  //iEvent.getByLabel( "famosSimHits", SimVtxHandle );
  iEvent.getByLabel( "g4SimHits", SimTrackHandle );
  iEvent.getByLabel( "g4SimHits", SimVtxHandle );

  /// Track Trigger
  edm::Handle< L1TkCluster_PixelDigi_Collection > PixelDigiL1TkClusterHandle;
  edm::Handle< L1TkStub_PixelDigi_Collection >    PixelDigiL1TkStubHandle;
  edm::Handle< L1TkStub_PixelDigi_Collection >    PixelDigiL1TkFailedStubHandle;
  iEvent.getByLabel( "L1TkClustersFromPixelDigis",             PixelDigiL1TkClusterHandle );
  iEvent.getByLabel( "L1TkStubsFromPixelDigis", "StubsPass",   PixelDigiL1TkStubHandle );
  iEvent.getByLabel( "L1TkStubsFromPixelDigis", "StubsFail",   PixelDigiL1TkFailedStubHandle );

  ///////////////////////////////////
  /// COLLECT CLUSTER INFORMATION ///
  ///////////////////////////////////

  /// Maps to store SimTrack information
  std::map< unsigned int, std::vector< edm::Ptr< SimTrack > > > simTrackPerLayer;
  std::map< unsigned int, std::vector< edm::Ptr< SimTrack > > > simTrackPerDisk;

  /// Go on only if there are L1TkCluster from PixelDigis
  if ( PixelDigiL1TkClusterHandle->size() > 0 )
  {
    /// Loop over L1TkClusters
    L1TkCluster_PixelDigi_Collection::const_iterator iterL1TkCluster;
    for ( iterL1TkCluster = PixelDigiL1TkClusterHandle->begin();
          iterL1TkCluster != PixelDigiL1TkClusterHandle->end();
          ++iterL1TkCluster )
    {
      StackedTrackerDetId detIdClu( iterL1TkCluster->getDetId() );
      unsigned int memberClu = iterL1TkCluster->getStackMember();
      bool genuineClu     = iterL1TkCluster->isGenuine();
      int partClu         = iterL1TkCluster->findType();
      unsigned int widClu = iterL1TkCluster->findWidth();
      GlobalPoint posClu  = theStackedGeometry->findAverageGlobalPosition( &(*iterL1TkCluster) );
      
      hCluster_RZ->Fill( posClu.z(), posClu.perp() );

      if ( detIdClu.isBarrel() )
      {
        hCluster_Barrel_Stack->Fill( detIdClu.iLayer() );
        hCluster_Barrel_XY->Fill( posClu.x(), posClu.y() );
        hCluster_Barrel_XY_Zoom->Fill( posClu.x(), posClu.y() );
      }
      else if ( detIdClu.isEndcap() )
      {
        hCluster_Endcap_Stack->Fill( detIdClu.iDisk() );
        if ( posClu.z() > 0 )
        {
          hCluster_Endcap_Fw_XY->Fill( posClu.x(), posClu.y() );
          hCluster_Endcap_Fw_RZ_Zoom->Fill( posClu.z(), posClu.perp() );
        }
        else
        {
          hCluster_Endcap_Bw_XY->Fill( posClu.x(), posClu.y() );
          hCluster_Endcap_Bw_RZ_Zoom->Fill( posClu.z(), posClu.perp() );
        }
      }

      hCluster_Mem->Fill( memberClu );
      hCluster_Gen->Fill( genuineClu, memberClu );
      hCluster_PID->Fill( partClu, memberClu );
      hCluster_W->Fill( widClu, memberClu );

      /// Store Track information in maps, skip if the Cluster is not good
      if ( !genuineClu ) continue;

      edm::Ptr< SimTrack > simTrackPtr = iterL1TkCluster->getSimTrackPtrs().at(0);

      if ( detIdClu.isBarrel() )
      {
        if ( simTrackPerLayer.find( detIdClu.iLayer() ) == simTrackPerLayer.end() )
        {
          std::vector< edm::Ptr< SimTrack > > tempVec;
          simTrackPerLayer.insert( make_pair( detIdClu.iLayer(), tempVec ) );
        }
        simTrackPerLayer[detIdClu.iLayer()].push_back( simTrackPtr );
      }
      else if ( detIdClu.isEndcap() )
      {
        if ( simTrackPerDisk.find( detIdClu.iDisk() ) == simTrackPerDisk.end() )
        {
          std::vector< edm::Ptr< SimTrack > > tempVec;
          simTrackPerDisk.insert( make_pair( detIdClu.iDisk(), tempVec ) );
        }
        simTrackPerDisk[detIdClu.iDisk()].push_back( simTrackPtr );
      }
    } /// End of Loop over L1TkClusters
  } /// End of if ( PixelDigiL1TkClusterHandle->size() > 0 )

  /// Clean the maps for SimTracks and fill histograms
  std::map< unsigned int, std::vector< edm::Ptr< SimTrack > > >::iterator iterSimTrackPerLayer;
  std::map< unsigned int, std::vector< edm::Ptr< SimTrack > > >::iterator iterSimTrackPerDisk;

  for ( iterSimTrackPerLayer = simTrackPerLayer.begin();
        iterSimTrackPerLayer != simTrackPerLayer.end();
        ++iterSimTrackPerLayer )
  {
    /// Remove duplicates, if any
    std::vector< edm::Ptr< SimTrack > > tempVec = iterSimTrackPerLayer->second;
    std::sort( tempVec.begin(), tempVec.end() );
    tempVec.erase( std::unique( tempVec.begin(), tempVec.end() ), tempVec.end() );

    /// Loop over the SimTracks in this piece of the map
    for ( unsigned int i = 0; i < tempVec.size(); i++ )
    {
      if ( tempVec.at(i).isNull() ) continue;
      SimTrack thisSimTrack = *(tempVec.at(i));
      mapCluLayer_hSimTrk_Pt[ iterSimTrackPerLayer->first ]->Fill( thisSimTrack.momentum().pt() );
      if ( thisSimTrack.momentum().pt() > 10.0 )
      {
        mapCluLayer_hSimTrk_Eta_Pt10[ iterSimTrackPerLayer->first ]->Fill( thisSimTrack.momentum().eta() );
        mapCluLayer_hSimTrk_Phi_Pt10[ iterSimTrackPerLayer->first ]->Fill( thisSimTrack.momentum().phi() > M_PI ?
                                                                           thisSimTrack.momentum().phi() - 2*M_PI :
                                                                           thisSimTrack.momentum().phi() );    
      }
    }
  }

  for ( iterSimTrackPerDisk = simTrackPerDisk.begin();
        iterSimTrackPerDisk != simTrackPerDisk.end();
        ++iterSimTrackPerDisk )
  {
    /// Remove duplicates, if any
    std::vector< edm::Ptr< SimTrack > > tempVec = iterSimTrackPerDisk->second;
    std::sort( tempVec.begin(), tempVec.end() );
    tempVec.erase( std::unique( tempVec.begin(), tempVec.end() ), tempVec.end() );

    /// Loop over the SimTracks in this piece of the map
    for ( unsigned int i = 0; i < tempVec.size(); i++ )
    {
      if ( tempVec.at(i).isNull() ) continue;
      SimTrack thisSimTrack = *(tempVec.at(i));
      mapCluDisk_hSimTrk_Pt[ iterSimTrackPerDisk->first ]->Fill( thisSimTrack.momentum().pt() );
      if ( thisSimTrack.momentum().pt() > 10.0 )
      {
        mapCluDisk_hSimTrk_Eta_Pt10[ iterSimTrackPerDisk->first ]->Fill( thisSimTrack.momentum().eta() );
        mapCluDisk_hSimTrk_Phi_Pt10[ iterSimTrackPerDisk->first ]->Fill( thisSimTrack.momentum().phi() > M_PI ?
                                                                         thisSimTrack.momentum().phi() - 2*M_PI :
                                                                         thisSimTrack.momentum().phi() );    
      }
    }
  }

  ////////////////////////////////
  /// COLLECT STUB INFORMATION ///
  ////////////////////////////////

  /// Maps to store SimTrack information
  std::map< unsigned int, std::vector< edm::Ptr< SimTrack > > > simTrackPerStubLayer;
  std::map< unsigned int, std::vector< edm::Ptr< SimTrack > > > simTrackPerStubDisk;

  /// Go on only if there are L1TkStub from PixelDigis
  if ( PixelDigiL1TkStubHandle->size() > 0 )
  {
    /// Loop over L1TkStubs
    L1TkStub_PixelDigi_Collection::const_iterator iterL1TkStub;
    for ( iterL1TkStub = PixelDigiL1TkStubHandle->begin();
          iterL1TkStub != PixelDigiL1TkStubHandle->end();
          ++iterL1TkStub )
    {
      StackedTrackerDetId detIdStub( iterL1TkStub->getDetId() );

      bool genuineStub    = iterL1TkStub->isGenuine();
      int partStub        = iterL1TkStub->findType();
      double displStub    = iterL1TkStub->getTriggerDisplacement();
      double offsetStub   = iterL1TkStub->getTriggerOffset();
      GlobalPoint posStub = theStackedGeometry->findGlobalPosition( &(*iterL1TkStub) );
      
      hStub_RZ->Fill( posStub.z(), posStub.perp() );

      if ( detIdStub.isBarrel() )
      {
        hStub_Barrel_XY->Fill( posStub.x(), posStub.y() );
        hStub_Barrel_XY_Zoom->Fill( posStub.x(), posStub.y() );
        hStub_Barrel_Stack->Fill( detIdStub.iLayer() );
      }
      else if ( detIdStub.isEndcap() )
      {
        hStub_Endcap_Stack->Fill( detIdStub.iDisk() );
        if (posStub.z() > 0) 
        {
          hStub_Endcap_Fw_XY->Fill( posStub.x(), posStub.y() );
          hStub_Endcap_Fw_RZ_Zoom->Fill( posStub.z(), posStub.perp() );
        }
        else
        {
          hStub_Endcap_Bw_XY->Fill( posStub.x(), posStub.y() );
          hStub_Endcap_Bw_RZ_Zoom->Fill( posStub.z(), posStub.perp() );
        }
      }

      hStub_Gen->Fill( genuineStub );
      hStub_PID->Fill( partStub );

      /// Store Track information in maps, skip if the Cluster is not good
      if ( !genuineStub ) continue;

      edm::Ptr< SimTrack > simTrackPtr = iterL1TkStub->getSimTrackPtr();

      if ( detIdStub.isBarrel() )
      {
        if ( simTrackPerStubLayer.find( detIdStub.iLayer() ) == simTrackPerStubLayer.end() ) {
          std::vector< edm::Ptr< SimTrack > > tempVec;
          simTrackPerStubLayer.insert( make_pair( detIdStub.iLayer(), tempVec ) );
        }
        simTrackPerStubLayer[detIdStub.iLayer()].push_back( simTrackPtr );

        hStub_Barrel_W->Fill( displStub - offsetStub, detIdStub.iLayer() );
        hStub_Barrel_O->Fill( offsetStub, detIdStub.iLayer() );
      }
      else if ( detIdStub.isEndcap() )
      {
        if ( simTrackPerStubDisk.find( detIdStub.iDisk() ) == simTrackPerStubDisk.end() ) {
          std::vector< edm::Ptr< SimTrack > > tempVec;
          simTrackPerStubDisk.insert( make_pair( detIdStub.iDisk(), tempVec ) );
        }
        simTrackPerStubDisk[detIdStub.iDisk()].push_back( simTrackPtr );

        hStub_Endcap_W->Fill( displStub - offsetStub, detIdStub.iDisk() );
        hStub_Endcap_O->Fill( offsetStub, detIdStub.iDisk() );
      }
      
      /// Compare to SimTrack

      if ( simTrackPtr.isNull() ) continue; /// This prevents to fill the vector if the SimTrack is not found
      SimTrack thisSimTrack = *simTrackPtr;

      double simPt = thisSimTrack.momentum().pt();
      double simEta = thisSimTrack.momentum().eta();
      double simPhi = thisSimTrack.momentum().phi();
      double recPt = theStackedGeometry->findRoughPt( mMagneticFieldStrength, &(*iterL1TkStub) );
      double recEta = theStackedGeometry->findGlobalDirection( &(*iterL1TkStub) ).eta();
      double recPhi = theStackedGeometry->findGlobalDirection( &(*iterL1TkStub) ).phi();

      if ( simPhi > M_PI )
        simPhi -= 2*M_PI;
      if ( recPhi > M_PI )
        recPhi -= 2*M_PI;

      if ( detIdStub.isBarrel() )
      {
        mapStubLayer_hStub_InvPt_SimTrk_InvPt[ detIdStub.iLayer() ]->Fill( 1./simPt, 1./recPt );
        mapStubLayer_hStub_Pt_SimTrk_Pt[ detIdStub.iLayer() ]->Fill( simPt, recPt );
        mapStubLayer_hStub_Eta_SimTrk_Eta[ detIdStub.iLayer() ]->Fill( simEta, recEta );
        mapStubLayer_hStub_Phi_SimTrk_Phi[ detIdStub.iLayer() ]->Fill( simPhi, recPhi );

        mapStubLayer_hStub_InvPtRes_SimTrk_Eta[ detIdStub.iLayer() ]->Fill( simEta, 1./recPt - 1./simPt );
        mapStubLayer_hStub_PtRes_SimTrk_Eta[ detIdStub.iLayer() ]->Fill( simEta, recPt - simPt );
        mapStubLayer_hStub_EtaRes_SimTrk_Eta[ detIdStub.iLayer() ]->Fill( simEta, recEta - simEta );
        mapStubLayer_hStub_PhiRes_SimTrk_Eta[ detIdStub.iLayer() ]->Fill( simEta, recPhi - simPhi );

        mapStubLayer_hStub_W_SimTrk_Pt[ detIdStub.iLayer() ]->Fill( simPt, displStub - offsetStub );
        mapStubLayer_hStub_W_SimTrk_InvPt[ detIdStub.iLayer() ]->Fill( 1./simPt, displStub - offsetStub );
      }
      else if ( detIdStub.isEndcap() )
      {
        mapStubDisk_hStub_InvPt_SimTrk_InvPt[ detIdStub.iDisk() ]->Fill( 1./simPt, 1./recPt );
        mapStubDisk_hStub_Pt_SimTrk_Pt[ detIdStub.iDisk() ]->Fill( simPt, recPt );
        mapStubDisk_hStub_Eta_SimTrk_Eta[ detIdStub.iDisk() ]->Fill( simEta, recEta );
        mapStubDisk_hStub_Phi_SimTrk_Phi[ detIdStub.iDisk() ]->Fill( simPhi, recPhi );

        mapStubDisk_hStub_InvPtRes_SimTrk_Eta[ detIdStub.iDisk() ]->Fill( simEta, 1./recPt - 1./simPt );
        mapStubDisk_hStub_PtRes_SimTrk_Eta[ detIdStub.iDisk() ]->Fill( simEta, recPt - simPt );
        mapStubDisk_hStub_EtaRes_SimTrk_Eta[ detIdStub.iDisk() ]->Fill( simEta, recEta - simEta );
        mapStubDisk_hStub_PhiRes_SimTrk_Eta[ detIdStub.iDisk() ]->Fill( simEta, recPhi - simPhi );

        mapStubDisk_hStub_W_SimTrk_Pt[ detIdStub.iDisk() ]->Fill( simPt, displStub - offsetStub );
        mapStubDisk_hStub_W_SimTrk_InvPt[ detIdStub.iDisk() ]->Fill( 1./simPt, displStub - offsetStub );
      }
    } /// End of loop over L1TkStubs
  } /// End of if ( PixelDigiL1TkStubHandle->size() > 0 )

  /// Clean the maps for SimTracks and fill histograms
  std::map< unsigned int, std::vector< edm::Ptr< SimTrack > > >::iterator iterSimTrackPerStubLayer;
  std::map< unsigned int, std::vector< edm::Ptr< SimTrack > > >::iterator iterSimTrackPerStubDisk;

  for ( iterSimTrackPerStubLayer = simTrackPerStubLayer.begin();
        iterSimTrackPerStubLayer != simTrackPerStubLayer.end();
        ++iterSimTrackPerStubLayer ) 
  {
    /// Remove duplicates, if any
    std::vector< edm::Ptr< SimTrack > > tempVec = iterSimTrackPerStubLayer->second;
    std::sort( tempVec.begin(), tempVec.end() );
    tempVec.erase( std::unique( tempVec.begin(), tempVec.end() ), tempVec.end() );

    /// Loop over the SimTracks in this piece of the map
    for ( unsigned int i = 0; i < tempVec.size(); i++ )
    {
      if ( tempVec.at(i).isNull() ) continue;
      SimTrack thisSimTrack = *(tempVec.at(i));
      mapStubLayer_hSimTrk_Pt[ iterSimTrackPerStubLayer->first ]->Fill( thisSimTrack.momentum().pt() );
      if ( thisSimTrack.momentum().pt() > 10.0 )
      {
        mapStubLayer_hSimTrk_Eta_Pt10[ iterSimTrackPerStubLayer->first ]->Fill( thisSimTrack.momentum().eta() );
        mapStubLayer_hSimTrk_Phi_Pt10[ iterSimTrackPerStubLayer->first ]->Fill( thisSimTrack.momentum().phi() > M_PI ?
                                                                                thisSimTrack.momentum().phi() - 2*M_PI :
                                                                                thisSimTrack.momentum().phi() );    
      }
    }
  }

  for ( iterSimTrackPerStubDisk = simTrackPerStubDisk.begin();
        iterSimTrackPerStubDisk != simTrackPerStubDisk.end();
        ++iterSimTrackPerStubDisk ) 
  {
    /// Remove duplicates, if any
    std::vector< edm::Ptr< SimTrack > > tempVec = iterSimTrackPerStubDisk->second;
    std::sort( tempVec.begin(), tempVec.end() );
    tempVec.erase( std::unique( tempVec.begin(), tempVec.end() ), tempVec.end() );

    /// Loop over the SimTracks in this piece of the map
    for ( unsigned int i = 0; i < tempVec.size(); i++ )
    {
      if ( tempVec.at(i).isNull() ) continue;
      SimTrack thisSimTrack = *(tempVec.at(i));
      mapStubDisk_hSimTrk_Pt[ iterSimTrackPerStubDisk->first ]->Fill( thisSimTrack.momentum().pt() );
      if ( thisSimTrack.momentum().pt() > 10.0 )
      {
        mapStubDisk_hSimTrk_Eta_Pt10[ iterSimTrackPerStubDisk->first ]->Fill( thisSimTrack.momentum().eta() );
        mapStubDisk_hSimTrk_Phi_Pt10[ iterSimTrackPerStubDisk->first ]->Fill( thisSimTrack.momentum().phi() > M_PI ?
                                                                              thisSimTrack.momentum().phi() - 2*M_PI :
                                                                              thisSimTrack.momentum().phi() );    
      }
    }
  }

  /// //////////////////////////
  /// SPECTRUM OF SIM TRACKS ///
  /// WITHIN PRIMARY VERTEX  ///
  /// CONSTRAINTS            ///
  /// //////////////////////////

  /// Go on only if there are SimTracks
  if ( SimTrackHandle->size() != 0 )
  {
    /// Loop over SimTracks
    edm::SimTrackContainer::const_iterator iterSimTracks;
    for ( iterSimTracks = SimTrackHandle->begin();
          iterSimTracks != SimTrackHandle->end();
          ++iterSimTracks )
    {
      /// Get the corresponding vertex
      int vertexIndex = iterSimTracks->vertIndex();
      const SimVertex& theSimVertex = (*SimVtxHandle)[vertexIndex];

      /// Assume perfectly round beamspot
      /// Correct and get the correct SimTrack Vertex position wrt beam center
      math::XYZTLorentzVectorD trkVtxPos = theSimVertex.position();

      /// First of all, check beamspot and correction
      hSimVtx_XY->Fill( trkVtxPos.x(), trkVtxPos.y() );
      hSimVtx_RZ->Fill( trkVtxPos.z(), trkVtxPos.rho() );

      /// Here we have only tracks form primary vertices
      /// Check Pt spectrum and pseudorapidity for over-threshold tracks
      hSimTrk_Pt->Fill( iterSimTracks->momentum().pt() );
      if ( iterSimTracks->momentum().pt() > 10.0 ) {
        hSimTrk_Eta_Pt10->Fill( iterSimTracks->momentum().eta() );
        hSimTrk_Phi_Pt10->Fill( iterSimTracks->momentum().phi() > M_PI ?
                                iterSimTracks->momentum().phi() - 2*M_PI :
                                iterSimTracks->momentum().phi() );
      }
    } /// End of Loop over SimTracks
  } /// End of if ( SimTrackHandle->size() != 0 )
} /// End of analyze()

///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(ValidateClusterStub);

