#ifndef EcalABAnalyzer_h_
#define EcalABAnalyzer_h_

// $Id: EcalABAnalyzer.h

#include <memory>
#include <vector>
#include <map>

#include <FWCore/Framework/interface/EDAnalyzer.h>

using namespace std;

class TShapeAnalysis;
class TAPDPulse;
class TMom;

// Define geometrical constants 
// NOT the same for "EB" and "EE"
//
//     "EB"       "EE"
//
//      0          0
//   1     2    1     2
//   3     4
//   5     6
//   7     8 
//
// 

// "EB" geometry
#define NCRYSEB    1700  // Number of crystals per EB supermodule

// "EE" geometry 
#define NCRYSEE    830   // Number of crystals per EE supermodule

class EcalABAnalyzer: public edm::EDAnalyzer{  

 public:
  
  explicit EcalABAnalyzer(const edm::ParameterSet& iConfig);  
  ~EcalABAnalyzer();
  
  
  virtual void analyze( const edm::Event & e, const  edm::EventSetup& c);
  virtual void beginJob(edm::EventSetup const&);
  virtual void endJob();
  
  
  enum VarCol   { iBlue, iRed, nColor }; 
  
 private:
  
  int iEvent;
  
  
  // Framework parameters
  
  unsigned int  _nsamples;
  unsigned int  _presample;
  unsigned int  _firstsample;
  unsigned int  _lastsample;
  unsigned int  _timingcutlow;
  unsigned int  _timingcuthigh;
  unsigned int  _timingquallow;
  unsigned int  _timingqualhigh;
  double        _ratiomincutlow;
  double        _ratiomincuthigh;
  double        _ratiomaxcutlow;
  double        _presamplecut;
  unsigned int  _niter ;
  double        _alpha ;
  double        _beta;
  unsigned int  _nevtmax;
  double        _noise;
  double        _chi2cut;
  string        _ecalPart;
  int           _fedid;
  double        _qualpercent;
  int           _debug;
  
  TAPDPulse *APDPulse;
  TMom *Delta01;
  TMom *Delta12;
  
  string  resdir_;
  string  digiCollection_;
  string  digiProducer_;
  string  eventHeaderCollection_;
  string  eventHeaderProducer_;
  
  // Output file names
  
  string  alphafile;
  string  alphainitfile;

  TShapeAnalysis *shapana;
  unsigned int nevtAB[NCRYSEB];
  
  //  Define geometrical constants
  //  Default values correspond to "EB" geometry (1700 crystals)
  
  unsigned int nCrys;
  bool doesABTreeExist;
    
  bool          _fitab;
  // Identify run type
  
  int runType;
  int runNum;
  int fedID;
  int dccID;
  int side;
  int lightside;
  int iZ;
  
  
  // Temporary root files and trees
  
  std::vector<int> colors;
  std::map<int, int> channelMapEE;
  std::vector<int> dccMEM;
  std::vector<int> modules;
  
  
  // Declaration of leaves types for temporary trees
  
  int             phi, eta;
  int             event ;
  int             color ;
  double          adc[10];
  int             adcG[10];
  int         channelIteratorEE;
  
  
  int iEta[NCRYSEB],iPhi[NCRYSEB];
  int iTowerID[NCRYSEB],iChannelID[NCRYSEB], idccID[NCRYSEB], iside[NCRYSEB];
  
  // Quality Checks variables and flags
  
  int nEvtBadGain[NCRYSEB];
  int nEvtBadTiming[NCRYSEB];
  int nEvtTot[NCRYSEB];

  bool          wasGainOK[NCRYSEB];
  bool          wasTimingOK[NCRYSEB];
  
  bool          isGainOK;
  bool          isTimingOK;
  
};

#endif


