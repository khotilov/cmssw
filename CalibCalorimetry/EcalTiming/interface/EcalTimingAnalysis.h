#ifndef EcalTimingAnalysis_H
#define EcalTimingAnalysis_H
/**\class EcalTimingAnalysis

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>

*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"

#include <string>
#include "TProfile.h"

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH3F.h"
#include <vector>

//#include<fstream>
//#include<map>
//#include<stl_pair>

//

class EcalTimingAnalysis : public edm::EDAnalyzer {
   public:
      explicit EcalTimingAnalysis( const edm::ParameterSet& );
      ~EcalTimingAnalysis();


      virtual void analyze(  edm::Event const&,   edm::EventSetup const& );
      virtual void beginJob(edm::EventSetup const&);
      virtual void endJob();
 private:
      
      std::string rootfile_;
      std::string hitCollection_;
      std::string hitCollectionEE_;
      std::string hitProducer_;
      std::string hitProducerEE_;
      std::string digiProducer_;
      float ampl_thr_;
      int min_num_ev_;
      int sm_;
      std::string txtFileName_;
      std::string txtFileForChGroups_;
      //std::string pndiodeProducer_;
      std::vector<double> sMAves_;
      std::vector<double> sMCorr_;
      
      TProfile* amplProfileConv_[54][4];
      TProfile* absoluteTimingConv_[54][4];

      TProfile* amplProfileAll_[54][4];
      TProfile* absoluteTimingAll_[54][4];
      
      TProfile* Chi2ProfileConv_[54][4];
      
      TProfile* relativeTimingBlueConv_[54];

      TGraphErrors* ttTiming_[54];
      TGraphErrors* ttTimingAll_;
      TGraphErrors* ttTimingRel_[54];
      TGraphErrors* ttTimingAllRel_;
      TGraphErrors* ttTimingAllSMChng_;
      
      TGraph* lasershiftVsTime_[54];
      TH2F* lasershiftVsTimehist_[54];
      TH1F* lasershiftLM_[54];
      TH1F* lasershift_;
      
      TH3F* ttTimingEtaPhi_;
      TH3F* ttTimingEtaPhiRel_;
      TH3F* chTimingEtaPhi_;
      TH3F* chTimingEtaPhiRel_;
	  
      TH3F* chTimingEtaPhiEEP_;
      TH3F* chTimingEtaPhiEEM_;
      TH3F* chTimingEtaPhiRelEEP_;
      TH3F* chTimingEtaPhiRelEEM_;
      
      TH1F* timeCry1[54]; 
      TH1F* timeCry2[54]; 
      TH1F* timeRelCry1[54]; 
      TH1F* timeRelCry2[54]; 
      
      TH1F* aveRelXtalTime_;
      TH1F* aveRelXtalTimebyDCC_[54];
      TH2F* aveRelXtalTimeVsAbsTime_;
      
      double timerunstart_;
      double timerunlength_;
	  
      TH1F* lasersPerEvt;

      const EcalElectronicsMapping* ecalElectronicsMap_;
 
      int ievt_;
      int numGoodEvtsPerSM_[54];
	  
      static const int numXtals = 15480;
  
      //Allows for running the job on a file
      bool fromfile_;
      std::string fromfilename_;   

};



#endif
