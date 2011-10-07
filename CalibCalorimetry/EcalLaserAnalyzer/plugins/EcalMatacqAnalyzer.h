// $Id: EcalMatacqAnalyzer.h

#include <memory>
#include <FWCore/Framework/interface/EDAnalyzer.h>

class TFile;
class TTree;
class TMTQ;
class TF1;
class TH1;

#define N_samples 2560
#define N_channels 1
#define NSIDES     2     // Number of sides
#define NCOL       2     // Number of sides
#define FFT2_SIZE   2048  // Number of bins used for FFT
#define FFT_SIZE   1048  // Number of bins used for FFT
#define FFT_START  850   // Keep signal starting at 850 ns

class EcalMatacqAnalyzer: public edm::EDAnalyzer{  

 public:
  
  explicit EcalMatacqAnalyzer(const edm::ParameterSet& iConfig);  
  ~EcalMatacqAnalyzer();
  
  
  virtual void analyze( const edm::Event & e, const  edm::EventSetup& c);
  virtual void beginJob();
  virtual void endJob();
  
  enum VarCol  { iBlue, iRed, nColor }; 
  enum VarSide { iSide0, iSide1, nSide }; 
  
 private:
    
  int iEvent;

//
// Framework parameters
//

  double       _presample;
  unsigned int _nsamplesaftmax;
  unsigned int _noiseCut;
  unsigned int _parabnbefmax;
  unsigned int _parabnaftmax;
  unsigned int _thres;
  unsigned int _lowlev;
  unsigned int _highlev;
  unsigned int _nevlasers;
  unsigned int _timebefmax;
  unsigned int _timeaftmax;
  double       _cutwindow;
  unsigned int _nsamplesshape;
  unsigned int _presampleshape;
  unsigned int _slide;
  int          _fedid;
  int          _debug;

  std::string resdir_;
  std::string digiCollection_;
  std::string digiProducer_;
  std::string eventHeaderCollection_;
  std::string eventHeaderProducer_;

  std::string outfile;
  std::string sampfile;


  // Identify run type

  unsigned int nSides;
  int lightside;
  int runType;
  int runNum;
  int  dccID;
  int  fedID;

  // Count Laser Events
  int laserEvents[2];
  int matacqEvents[2];
  bool isThereMatacq;

  //Declaration of leaves types
  
  int   event ;
  int   color ;
  double  matacq[N_samples]; 
  int  maxsamp; 
  int  nsamples; 
  double tt;
  // vector<int> vernier;

  TFile *sampFile;
  TTree *tree;

  
  TMTQ *MTQ[nColor][nSide];
  TTree *meanTree[nColor];

  std::vector<int> colors;
  std::vector<int> sides;

  TFile *outFile;
  int status;
  double peak, sigma, fit, ampl, trise, fwhm, fw20, fw80, ped, pedsig, ttrig, sliding;
  TTree* mtqShape;


};


