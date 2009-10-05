#include "CalibTracker/SiStripQuality/interface/SiStripBadAPVAlgorithmFromClusterOccupancy.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"


SiStripBadAPVAlgorithmFromClusterOccupancy::SiStripBadAPVAlgorithmFromClusterOccupancy(const edm::ParameterSet& iConfig):
  lowoccupancy_(0),
  highoccupancy_(100),
  absolutelow_(0),
  numberiterations_(2),
  Nevents_(0),
  occupancy_(0),
  OutFileName_("Occupancy.root")
  {
    minNevents_=Nevents_*occupancy_;
  }

SiStripBadAPVAlgorithmFromClusterOccupancy::~SiStripBadAPVAlgorithmFromClusterOccupancy(){
  LogTrace("SiStripBadAPVAlgorithmFromClusterOccupancy")<<"[SiStripBadAPVAlgorithmFromClusterOccupancy::~SiStripBadAPVAlgorithmFromClusterOccupancy] "<<std::endl;
}

void SiStripBadAPVAlgorithmFromClusterOccupancy::extractBadAPVs(SiStripQuality* siStripQuality,HistoMap& DM){

  LogTrace("SiStripBadAPVAlgorithmFromClusterOccupancy")<<"[SiStripBadAPVAlgorithmFromClusterOccupancy::extractBadAPVs] "<<std::endl;

  if (WriteOutputFile_==true){
  f = new TFile(OutFileName_.c_str(),"RECREATE");
  f->cd();

  apvtree = new TTree("moduleOccupancy","tree");

  apvtree->Branch("DetRawId",                &detrawid,                "DetRawId/I");
  apvtree->Branch("SubDetId",                &subdetid,                "SubDetId/I");
  apvtree->Branch("Layer_Ring",              &layer_ring,              "Layer_Ring/I");
  apvtree->Branch("Disc",                    &disc,                    "Disc/I");
  apvtree->Branch("IsBack",                  &isback,                  "IsBack/I");
  apvtree->Branch("IsExternalString",        &isexternalstring,        "IsExternalString/I");
  apvtree->Branch("IsZMinusSide",            &iszminusside,            "IsZMinusSide/I");
  apvtree->Branch("RodStringPetal",          &rodstringpetal,          "RodStringPetal/I");
  apvtree->Branch("IsStereo",                &isstereo,                "IsStereo/I");
  apvtree->Branch("ModuleNumber",            &module_number,           "ModuleNumber/I");
  apvtree->Branch("NumberOfStrips",          &number_strips,           "NumberOfStrips/I");
  apvtree->Branch("APVGlobalPositionX",      &global_position_x,       "APVGlobalPositionX/F");
  apvtree->Branch("APVGlobalPositionY",      &global_position_y,       "APVGlobalPositionY/F");
  apvtree->Branch("APVGlobalPositionZ",      &global_position_z,       "APVGlobalPositionZ/F");
  apvtree->Branch("APVNumber",               &apv_number,              "APVNumber/I");
  apvtree->Branch("APVAbsoluteOccupancy",    &apvAbsoluteOccupancy,    "apvAbsoluteOccupancy/D");
  apvtree->Branch("APVMedianOccupancy",      &apvMedianOccupancy,      "apvMedianOccupancy/D");
  }

  HistoMap::iterator it=DM.begin();
  HistoMap::iterator itEnd=DM.end();
  std::vector<unsigned int> badStripList;
  uint32_t detid;
  for (;it!=itEnd;++it){

    Apv APV;

    for (int apv=0; apv<6; apv++)
      {
	APV.apvMedian[apv]        = 0;
	apvabsoluteOccupancy[apv] = 0;

	for (int strip=0; strip<128; strip++)
	  {
	    stripOccupancy[apv][strip] = 0;
	    stripWeight[apv][strip]    = 0;
	  }
      }

    pHisto phisto;
    phisto._th1f = it->second.get();
    phisto._NEntries = (int)phisto._th1f->GetEntries();
    phisto._NBins = phisto._th1f->GetNbinsX();

    number_strips  = (int)phisto._NBins;
    number_apvs    = number_strips/128;
    APV.numberApvs = number_apvs;

    for (int apv=0; apv<number_apvs; apv++)
      {
	for (int strip=0; strip<128; strip++)
	  {
	    stripOccupancy[apv][strip] = phisto._th1f->GetBinContent((apv*128)+strip+1); // Remember: Bin=0 is underflow bin!
	    stripWeight[apv][strip]    = 1;
	    apvabsoluteOccupancy[apv] += phisto._th1f->GetBinContent((apv*128)+strip+1); // Remember: Bin=0 is underflow bin!
	  }
      }

    for (int apv=0; apv<number_apvs; apv++)
      {
	APV.apvMedian[apv] = TMath::Median(128,stripOccupancy[apv],stripWeight[apv]);
      }

    detid=it->first;
    DetId detectorId=DetId(detid);

    if (edm::isDebugEnabled())
      LogTrace("SiStripBadAPV") << "Analyzing detid " << detid<< std::endl;

    detrawid     = detid;
    APV.detrawId = detrawid;
    subdetid     = detectorId.subdetId();
    if (SiStripDetId(detrawid).stereo() !=0 ) isstereo = 1; // It's a stereo module
    else                                      isstereo = 0; // It's an rphi module
    switch (detectorId.subdetId())
      {
      case StripSubdetector::TIB :
	layer_ring = TIBDetId(detrawid).layer();
	disc       = -1;
	isback     = -1;
	if (TIBDetId(detrawid).isExternalString()) isexternalstring = 1;
	else                                       isexternalstring = 0;
	if (TIBDetId(detrawid).isZMinusSide()) iszminusside = 1;
	else                                   iszminusside = 0;
	rodstringpetal     = TIBDetId(detrawid).stringNumber();
	module_number      = TIBDetId(detrawid).moduleNumber();
	APV.modulePosition = module_number;

	if      (layer_ring == 1) medianValues_TIB_Layer1.push_back(APV);
	else if (layer_ring == 2) medianValues_TIB_Layer2.push_back(APV);
	else if (layer_ring == 3) medianValues_TIB_Layer3.push_back(APV);
	else if (layer_ring == 4) medianValues_TIB_Layer4.push_back(APV);
	break;

      case StripSubdetector::TID :
	layer_ring = TIDDetId(detrawid).ring();
	disc       = TIDDetId(detrawid).wheel();
	if (TIDDetId(detrawid).isBackRing()) isback = 1;
	else                                 isback = 0;
	if (TIDDetId(detrawid).isZMinusSide()) iszminusside = 1;
	else                                   iszminusside = 0;
	isexternalstring   = -1;
	rodstringpetal     = -1;
	module_number      = TIDDetId(detrawid).moduleNumber();
	APV.modulePosition = layer_ring;

	if (iszminusside==0)
	  {
	    if      (disc==1) medianValues_TIDPlus_Disc1.push_back(APV);
	    else if (disc==2) medianValues_TIDPlus_Disc2.push_back(APV);
	    else if (disc==3) medianValues_TIDPlus_Disc3.push_back(APV);
	  }
	else if (iszminusside==1)
	  {
	    if      (disc==1) medianValues_TIDMinus_Disc1.push_back(APV);
	    else if (disc==2) medianValues_TIDMinus_Disc2.push_back(APV);
	    else if (disc==3) medianValues_TIDMinus_Disc3.push_back(APV);
	  }
	break;

      case StripSubdetector::TOB :
	layer_ring = TOBDetId(detrawid).layer();
	disc       = -1;
	isback     = -1;
	if (TOBDetId(detrawid).isZMinusSide()) iszminusside = 1;
	else                                   iszminusside = 0;
	isexternalstring   = -1;
	rodstringpetal     = TOBDetId(detrawid).rodNumber();
	module_number      = TOBDetId(detrawid).moduleNumber();
	APV.modulePosition = module_number;

	if      (layer_ring == 1) medianValues_TOB_Layer1.push_back(APV);
	else if (layer_ring == 2) medianValues_TOB_Layer2.push_back(APV);
	else if (layer_ring == 3) medianValues_TOB_Layer3.push_back(APV);
	else if (layer_ring == 4) medianValues_TOB_Layer4.push_back(APV);
	else if (layer_ring == 5) medianValues_TOB_Layer5.push_back(APV);
	else if (layer_ring == 6) medianValues_TOB_Layer6.push_back(APV);
	break;

      case StripSubdetector::TEC :
	layer_ring = TECDetId(detrawid).ring();
	disc       = TECDetId(detrawid).wheel();
	if (TECDetId(detrawid).isBackPetal()) isback = 1;
	else                                  isback = 0;
	if (TECDetId(detrawid).isZMinusSide()) iszminusside = 1;
	else                                   iszminusside = 0;
	isexternalstring   = -1;
	rodstringpetal     = TECDetId(detrawid).petalNumber();
	module_number      = TECDetId(detrawid).moduleNumber();
	APV.modulePosition = layer_ring;

	if (iszminusside==0)
	  {
	    if      (disc==1) medianValues_TECPlus_Disc1.push_back(APV);
	    else if (disc==2) medianValues_TECPlus_Disc2.push_back(APV);
	    else if (disc==3) medianValues_TECPlus_Disc3.push_back(APV);
	    else if (disc==4) medianValues_TECPlus_Disc4.push_back(APV);
	    else if (disc==5) medianValues_TECPlus_Disc5.push_back(APV);
	    else if (disc==6) medianValues_TECPlus_Disc6.push_back(APV);
	    else if (disc==7) medianValues_TECPlus_Disc7.push_back(APV);
	    else if (disc==8) medianValues_TECPlus_Disc8.push_back(APV);
	    else if (disc==9) medianValues_TECPlus_Disc9.push_back(APV);
	  }
	else if (iszminusside==1)
	  {
	    if      (disc==1) medianValues_TECMinus_Disc1.push_back(APV);
	    else if (disc==2) medianValues_TECMinus_Disc2.push_back(APV);
	    else if (disc==3) medianValues_TECMinus_Disc3.push_back(APV);
	    else if (disc==4) medianValues_TECMinus_Disc4.push_back(APV);
	    else if (disc==5) medianValues_TECMinus_Disc5.push_back(APV);
	    else if (disc==6) medianValues_TECMinus_Disc6.push_back(APV);
	    else if (disc==7) medianValues_TECMinus_Disc7.push_back(APV);
	    else if (disc==8) medianValues_TECMinus_Disc8.push_back(APV);
	    else if (disc==9) medianValues_TECMinus_Disc9.push_back(APV);
	  }
	break;

      default :
	std::cout << "### Detector does not belong to TIB, TID, TOB or TEC !? ###" << std::endl;
	std::cout << "### DetRawId: " << detrawid << " ###" << std::endl;
      }

    const StripGeomDetUnit*  theStripDet = dynamic_cast<const StripGeomDetUnit*>( (TkGeom->idToDet(detectorId)) );
    const StripTopology* theStripTopol   = dynamic_cast<const StripTopology*>( &(theStripDet->specificTopology()) );

    for (int apv=0; apv<number_apvs; apv++)
      {
	apv_number           = apv+1;
	apvMedianOccupancy   = APV.apvMedian[apv];
	apvAbsoluteOccupancy = apvabsoluteOccupancy[apv];

	LocalPoint  pos_strip_local  = theStripTopol->localPosition((apv*128));
        GlobalPoint pos_strip_global = (TkGeom->idToDet(detectorId))->surface().toGlobal(pos_strip_local);

        global_position_x = pos_strip_global.x();
        global_position_y = pos_strip_global.y();
        global_position_z = pos_strip_global.z();

	if (WriteOutputFile_==true) apvtree->Fill();
      }

  } // end loop on modules

  // Calculate Mean and RMS for each Layer
  CalculateMeanAndRMS(medianValues_TIB_Layer1,MeanAndRms_TIB_Layer1,numberiterations_);
  CalculateMeanAndRMS(medianValues_TIB_Layer2,MeanAndRms_TIB_Layer2,numberiterations_);
  CalculateMeanAndRMS(medianValues_TIB_Layer3,MeanAndRms_TIB_Layer3,numberiterations_);
  CalculateMeanAndRMS(medianValues_TIB_Layer4,MeanAndRms_TIB_Layer4,numberiterations_);

  CalculateMeanAndRMS(medianValues_TOB_Layer1,MeanAndRms_TOB_Layer1,numberiterations_);
  CalculateMeanAndRMS(medianValues_TOB_Layer2,MeanAndRms_TOB_Layer2,numberiterations_);
  CalculateMeanAndRMS(medianValues_TOB_Layer3,MeanAndRms_TOB_Layer3,numberiterations_);
  CalculateMeanAndRMS(medianValues_TOB_Layer4,MeanAndRms_TOB_Layer4,numberiterations_);
  CalculateMeanAndRMS(medianValues_TOB_Layer5,MeanAndRms_TOB_Layer5,numberiterations_);
  CalculateMeanAndRMS(medianValues_TOB_Layer6,MeanAndRms_TOB_Layer6,numberiterations_);

  CalculateMeanAndRMS(medianValues_TIDPlus_Disc1,MeanAndRms_TIDPlus_Disc1,numberiterations_);
  CalculateMeanAndRMS(medianValues_TIDPlus_Disc2,MeanAndRms_TIDPlus_Disc2,numberiterations_);
  CalculateMeanAndRMS(medianValues_TIDPlus_Disc3,MeanAndRms_TIDPlus_Disc3,numberiterations_);
  CalculateMeanAndRMS(medianValues_TIDMinus_Disc1,MeanAndRms_TIDMinus_Disc1,numberiterations_);
  CalculateMeanAndRMS(medianValues_TIDMinus_Disc2,MeanAndRms_TIDMinus_Disc2,numberiterations_);
  CalculateMeanAndRMS(medianValues_TIDMinus_Disc3,MeanAndRms_TIDMinus_Disc3,numberiterations_);

  CalculateMeanAndRMS(medianValues_TECPlus_Disc1,MeanAndRms_TECPlus_Disc1,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECPlus_Disc2,MeanAndRms_TECPlus_Disc2,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECPlus_Disc3,MeanAndRms_TECPlus_Disc3,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECPlus_Disc4,MeanAndRms_TECPlus_Disc4,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECPlus_Disc5,MeanAndRms_TECPlus_Disc5,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECPlus_Disc6,MeanAndRms_TECPlus_Disc6,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECPlus_Disc7,MeanAndRms_TECPlus_Disc7,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECPlus_Disc8,MeanAndRms_TECPlus_Disc8,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECPlus_Disc9,MeanAndRms_TECPlus_Disc9,numberiterations_);

  CalculateMeanAndRMS(medianValues_TECMinus_Disc1,MeanAndRms_TECMinus_Disc1,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECMinus_Disc2,MeanAndRms_TECMinus_Disc2,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECMinus_Disc3,MeanAndRms_TECMinus_Disc3,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECMinus_Disc4,MeanAndRms_TECMinus_Disc4,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECMinus_Disc5,MeanAndRms_TECMinus_Disc5,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECMinus_Disc6,MeanAndRms_TECMinus_Disc6,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECMinus_Disc7,MeanAndRms_TECMinus_Disc7,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECMinus_Disc8,MeanAndRms_TECMinus_Disc8,numberiterations_);
  CalculateMeanAndRMS(medianValues_TECMinus_Disc9,MeanAndRms_TECMinus_Disc9,numberiterations_);

  pQuality=siStripQuality;
  badStripList.clear();

  // Analyze the APV Occupancy for hot APVs
  AnalyzeOccupancy(siStripQuality,medianValues_TIB_Layer1,MeanAndRms_TIB_Layer1,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TIB_Layer2,MeanAndRms_TIB_Layer2,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TIB_Layer3,MeanAndRms_TIB_Layer3,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TIB_Layer4,MeanAndRms_TIB_Layer4,badStripList);

  AnalyzeOccupancy(siStripQuality,medianValues_TOB_Layer1,MeanAndRms_TOB_Layer1,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TOB_Layer2,MeanAndRms_TOB_Layer2,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TOB_Layer3,MeanAndRms_TOB_Layer3,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TOB_Layer4,MeanAndRms_TOB_Layer4,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TOB_Layer5,MeanAndRms_TOB_Layer5,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TOB_Layer6,MeanAndRms_TOB_Layer6,badStripList);

  AnalyzeOccupancy(siStripQuality,medianValues_TIDPlus_Disc1,MeanAndRms_TIDPlus_Disc1,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TIDPlus_Disc2,MeanAndRms_TIDPlus_Disc2,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TIDPlus_Disc3,MeanAndRms_TIDPlus_Disc3,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TIDMinus_Disc1,MeanAndRms_TIDMinus_Disc1,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TIDMinus_Disc2,MeanAndRms_TIDMinus_Disc2,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TIDMinus_Disc3,MeanAndRms_TIDMinus_Disc3,badStripList);

  AnalyzeOccupancy(siStripQuality,medianValues_TECPlus_Disc1,MeanAndRms_TECPlus_Disc1,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECPlus_Disc2,MeanAndRms_TECPlus_Disc2,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECPlus_Disc3,MeanAndRms_TECPlus_Disc3,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECPlus_Disc4,MeanAndRms_TECPlus_Disc4,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECPlus_Disc5,MeanAndRms_TECPlus_Disc5,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECPlus_Disc6,MeanAndRms_TECPlus_Disc6,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECPlus_Disc7,MeanAndRms_TECPlus_Disc7,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECPlus_Disc8,MeanAndRms_TECPlus_Disc8,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECPlus_Disc9,MeanAndRms_TECPlus_Disc9,badStripList);

  AnalyzeOccupancy(siStripQuality,medianValues_TECMinus_Disc1,MeanAndRms_TECMinus_Disc1,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECMinus_Disc2,MeanAndRms_TECMinus_Disc2,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECMinus_Disc3,MeanAndRms_TECMinus_Disc3,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECMinus_Disc4,MeanAndRms_TECMinus_Disc4,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECMinus_Disc5,MeanAndRms_TECMinus_Disc5,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECMinus_Disc6,MeanAndRms_TECMinus_Disc6,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECMinus_Disc7,MeanAndRms_TECMinus_Disc7,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECMinus_Disc8,MeanAndRms_TECMinus_Disc8,badStripList);
  AnalyzeOccupancy(siStripQuality,medianValues_TECMinus_Disc9,MeanAndRms_TECMinus_Disc9,badStripList);

  siStripQuality->fillBadComponents();

  if (WriteOutputFile_==true){
  f->cd();
  apvtree->Write();
  f->Close();
  }

  LogTrace("SiStripBadAPV") << ss.str() << std::endl;
}


void SiStripBadAPVAlgorithmFromClusterOccupancy::CalculateMeanAndRMS(std::vector<Apv> a, std::pair<double,double>* MeanRMS, int number_iterations)
{
  Double_t tot[7], tot2[7];
  Double_t n[7];

  Double_t Mean[7] = {0};
  Double_t Rms[7]  = {1000,1000,1000,1000,1000,1000,1000};

  int Moduleposition;

  for (int i=0; i<number_iterations; i++)
    {
      for (int j=0; j<7; j++)
	{
	  n[j]    = 0;
	  tot[j]  = 0;
	  tot2[j] = 0;
	}

      for (uint32_t it=0; it<a.size(); it++)
	{
	  Moduleposition = (a[it].modulePosition)-1;

	  for (int apv=0; apv<a[it].numberApvs; apv++)
	    {
	      if (i>0)
		{
		  if (a[it].apvMedian[apv]<(Mean[Moduleposition]-3*Rms[Moduleposition]) || (a[it].apvMedian[apv]>(Mean[Moduleposition]+5*Rms[Moduleposition])))
		    {
		      continue;
		    }
		}
	      tot[Moduleposition]  += a[it].apvMedian[apv];
	      tot2[Moduleposition] += (a[it].apvMedian[apv])*(a[it].apvMedian[apv]);
	      n[Moduleposition]++;
	    }
	}

      for (int j=0; j<7; j++)
	{
	  if (n[j]!=0)
	    {
	      Mean[j] = tot[j]/n[j];
	      Rms[j]  = TMath::Sqrt(TMath::Abs(tot2[j]/n[j] -Mean[j]*Mean[j]));
	    }
	}
    }

  for (int j=0; j<7; j++)
    {
      MeanRMS[j] = std::make_pair(Mean[j],Rms[j]);
    }

}

void SiStripBadAPVAlgorithmFromClusterOccupancy::AnalyzeOccupancy(SiStripQuality* quality, std::vector<Apv>& medianValues, std::pair<double,double>* MeanAndRms, std::vector<unsigned int>& BadStripList)
{
  int Moduleposition;

  for (uint32_t it=0; it<medianValues.size(); it++)
    {
      Moduleposition = (medianValues[it].modulePosition)-1;

      for (int apv=0; apv<medianValues[it].numberApvs; apv++)
	{
	  if (medianValues[it].apvMedian[apv] > minNevents_)
	    {
	      if ((medianValues[it].apvMedian[apv]>(MeanAndRms[Moduleposition].first+highoccupancy_*MeanAndRms[Moduleposition].second)) && (medianValues[it].apvMedian[apv]>absolutelow_))
		BadStripList.push_back(pQuality->encode((apv*128),128,0));
	    }
	  else if (medianValues[it].apvMedian[apv]<(MeanAndRms[Moduleposition].first-lowoccupancy_*MeanAndRms[Moduleposition].second))
	    BadStripList.push_back(pQuality->encode((apv*128),128,0));
	}
      if (BadStripList.begin()!=BadStripList.end())
	{
	  quality->compact(medianValues[it].detrawId,BadStripList);
	  SiStripQuality::Range range(BadStripList.begin(),BadStripList.end());
	  quality->put(medianValues[it].detrawId,range);
	}
      BadStripList.clear();
    }
}

void SiStripBadAPVAlgorithmFromClusterOccupancy::setMinNumOfEvents()
{
  minNevents_=occupancy_*Nevents_;
}
