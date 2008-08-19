#ifndef _WBMHTMLGENERATOR_H_
#define _WBMHTMLGENERATOR_H_

// ROOT Schema Headers
#include "RecoLuminosity/ROOTSchema/interface/ROOTFileReader.h"
#include "RecoLuminosity/TCPReceiver/interface/TimeStamp.h"
#include "RecoLuminosity/ROOTSchema/interface/FileToolKit.h"

#include <vector>

namespace HCAL_HLX{

   class HTMLGenerator: public ROOTFileReader, private TimeStamp, public FileToolKit{
      public:
	 HTMLGenerator();
    
	 ~HTMLGenerator();
    
	 void CreateWebPage();

	 void SetOutputDir(std::string outDir){ outputDir_ = outDir; }
	 void SetHistoBins(unsigned int NBins, double XMin, double XMax);
    
      private:
    
	 std::string GetRunDir();
	 std::string GetLSDir();
	 std::string GetHLXDir(const unsigned short int &wedgeNum);
	 std::string GetHLXPicDir(const unsigned short int &wedgeNum);

	 void GenerateIndexPage();
    
	 void GenerateRunPage();
	 void GenerateSectionPage();
    
	 void GenerateAveragePlots();
	 void GenerateAveragePage();
    
	 void GenerateComparePlots();
	 void GenerateComparePage();
    
	 void GenerateHLXPage(const unsigned short int &wedgeNum);
	 void GenerateHLXPlots(const unsigned short int &wedgeNum);

//       void GenerateEtSumPage();

	 void GenerateHistoGroupPage(const std::string &histoType);
	 void GenerateLumiPage();

	 std::string fileName_;
	 std::string outputFile_;
	 std::string outputDir_;
	 std::string plotExt_;
	 std::string baseURL;
    
	 mode_t writeMode;

	 std::vector< std::string > HistoNames;
    
	 unsigned int  previousRunNumber;
    
	 unsigned int NBins_;
	 double XMin_;
	 double XMax_;

	 int iEta_[36];
	 int iPhi_[36];

	 HCAL_HLX::LUMI_SECTION lumiSection_;

	 std::string HLXToHFMap_[36];


   };
  
}

#endif
