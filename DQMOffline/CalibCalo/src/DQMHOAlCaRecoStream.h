
#include "FWCore/Framework/interface/Frameworkfwd.h"

class DQMHOAlCaRecoStream : public edm::EDAnalyzer {
   public:
      explicit DQMHOAlCaRecoStream(const edm::ParameterSet&);
      ~DQMHOAlCaRecoStream();

   private:

      DQMStore* dbe_; 

      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void endRun(const edm::Run& r, const edm::EventSetup& c);


  MonitorElement* hMuonMultipl;
  MonitorElement* hMuonMom;
  MonitorElement* hMuonEta;
  MonitorElement* hMuonPhi;

  MonitorElement* hDirCosine;
  MonitorElement* hHOTime;
  
  MonitorElement* hSigRing[5];
  //  MonitorElement* hSigRingm1;
  //  MonitorElement* hSigRing00;
  //  MonitorElement* hSigRingp1;
  //  MonitorElement* hSigRingp2;

  MonitorElement* hPedRing[5];
  //  MonitorElement* hPedRingm1;
  //  MonitorElement* hPedRing00;
  //  MonitorElement* hPedRingp1;
  //  MonitorElement* hPedRingp2;

  MonitorElement* hSignal3x3[9];

  int Nevents;
  int Nmuons;

  std::string theRootFileName;
  std::string folderName_;
  double m_sigmaValue;

  double m_lowRadPosInMuch;
  double m_highRadPosInMuch;

  int    m_nbins;
  double m_lowEdge;
  double m_highEdge;

  bool saveToFile_;

      // ----------member data ---------------------------

};
