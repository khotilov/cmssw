// Original Author:  Loic Quertenmont

#ifndef HSCP_ANALYSIS_SAMPLE
#define HSCP_ANALYSIS_SAMPLE

class stSample{
   public:
   std::string CMSSW;
   int         Type;
   std::string Name;
   std::string FileName;
   std::string Legend;
   std::string Pileup;
   double      Mass;
   double      XSec;
   bool        MakePlot;
   float       WNC0;//weight for signal event with 0 Charged HSCP
   float       WNC1;//weight for signal event with 1 Charged HSCP
   float       WNC2;//weight for signal event with 2 Charged HSCP

   //These weights are used to reweight the gluino samples for different scenario of gluinobll fraction (f).  
   //The weights are independent on the mass and are computed as a reweighting factor w.r.t a f=10% gluino ball sample.
   //Weights have been computed in 2010 to be : 
   //for f=00% --> [0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015]
   //for f=10% --> [0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015]
   //for f=50% --> [0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015]

   stSample(){};
   stSample(std::string CMSSW_, int Type_, std::string Name_, std::string FileName_, std::string Legend_, std::string Pileup_, double Mass_, double XSec_, bool MakePlot_, float WNC0_, float WNC1_, float WNC2_){CMSSW=CMSSW_; Type=Type_; Name=Name_; FileName=FileName_; Legend=Legend_; Pileup=Pileup_; Mass=Mass_; XSec=XSec_; MakePlot=MakePlot_; WNC0=WNC0_; WNC1=WNC1_; WNC2=WNC2_;}

   int readFromFile(FILE* pFile){ 
      char line[4096];
      char* toReturn = fgets(line, 4096, pFile);
      if(!toReturn)return EOF;
      char* pch=strtok(line,","); int Arg=0; string tmp; int temp;
      while (pch!=NULL){
         switch(Arg){
            case  0: tmp = pch;  CMSSW    = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  1: sscanf(pch, "%d", &Type); break;
            case  2: tmp = pch;  Name     = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  3: tmp = pch;  FileName = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  4: tmp = pch;  Legend   = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  5: tmp = pch;  Pileup   = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  6: sscanf(pch, "%lf", &Mass); break;
            case  7: sscanf(pch, "%lf", &XSec); break;
            case  8: sscanf(pch, "%d", &temp); MakePlot=(temp>0); break;
            case  9: sscanf(pch, "%f", &WNC0); break;
            case 10: sscanf(pch, "%f", &WNC1); break;
            case 11: sscanf(pch, "%f", &WNC2); break;
            default:return EOF;break;
         }
         pch=strtok(NULL,",");Arg++;
      }
      return 0;
   }

   void print(FILE* pFile=stdout){
      fprintf(pFile, "%-9s, %1d, %-30s, %-40s, %-60s, %-8s, % 5g, %+12.10E, %1i, %5.3f, %5.3f, %5.3f\n", (string("\"")+CMSSW+"\"").c_str(), Type, (string("\"")+Name+"\"").c_str(), (string("\"")+FileName+"\"").c_str(), (string("\"")+Legend+"\"").c_str(), (string("\"")+Pileup+"\"").c_str(), Mass, XSec, MakePlot, WNC0, WNC1, WNC2);
   }

   std::string ModelName(){
      char MassStr[255];sprintf(MassStr, "%.0f",Mass);
      string toReturn=Name; toReturn.erase(toReturn.find(MassStr), string(MassStr).length());
      return toReturn;
   }

   std::string ModelLegend(){
      char MassStr[255];sprintf(MassStr, "%.0f",Mass);
      string toReturn=Legend; toReturn.erase(toReturn.find(MassStr), string(MassStr).length());
      if(toReturn.find(" GeV/#font[12]{c}^{2}")!=string::npos)toReturn.erase(toReturn.find(" GeV/#font[12]{c}^{2}"), string(" GeV/#font[12]{c}^{2}").length());
      return toReturn;
   }

   double GetFGluinoWeight(int NChargedHSCP){ 
      if(Type!=2)return 1;
      double Weight;
      switch(NChargedHSCP){
         case 0: Weight=WNC0; break;
         case 1: Weight=WNC1; break;
         case 2: Weight=WNC2; break;
         default: return 1.0;
      }
      if(Weight<0)return 1.0;
      return Weight;
   }

};

void GetSampleDefinition(std::vector<stSample>& samples, std::string sampleTxtFile="Analysis_Samples.txt"){
      FILE* pFile = fopen(sampleTxtFile.c_str(),"r");
         if(!pFile){printf("Can't open %s\n","Analysis_Samples.txt"); return;}
         stSample newSample;      
         while(newSample.readFromFile(pFile)!=EOF){samples.push_back(newSample);}
      fclose(pFile);
}

void GetInputFiles(stSample sample, std::string BaseDirectory_, std::vector<std::string>& inputFiles, int period=0){
   std::vector<string> fileNames;
   char* tmp = (char*)sample.FileName.c_str();
   char* pch=strtok(tmp,";");
   while (pch!=NULL){
      fileNames.push_back(pch);
      pch=strtok(NULL,",");
   }

   for(unsigned int f=0;f<fileNames.size();f++){  
      if(sample.Type>=2){ //MC Signal
        if (period==0) inputFiles.push_back(BaseDirectory_ + fileNames[f] + "_BX0.root");
        if (period==1) inputFiles.push_back(BaseDirectory_ + fileNames[f] + "_BX1.root");
      }else{ //Data or MC Background
        inputFiles.push_back(BaseDirectory_ + fileNames[f] + ".root");
      }
   }
}

int JobIdToIndex(string JobId, const std::vector<stSample>& samples){
   for(unsigned int s=0;s<samples.size();s++){
      if(samples[s].Name==JobId)return s;
   }return -1;
}

void keepOnlySamplesOfTypeX(std::vector<stSample>& samples, int TypeX){
   for(unsigned int s=0;s<samples.size();s++){if(samples[s].Type!=TypeX){samples.erase(samples.begin()+s);s--;} }
}

void keepOnlyTheXtoYSamples(std::vector<stSample>& samples, unsigned int X, unsigned int Y){
   std::vector<stSample> tmp; for(unsigned int s=X;s<=Y && s<samples.size();s++){tmp.push_back(samples[s]);} samples.clear();
   for(unsigned int s=0;s<tmp.size();s++)samples.push_back(tmp[s]);
}

void keepOnlySamplesOfNameX(std::vector<stSample>& samples, string NameX){
   for(unsigned int s=0;s<samples.size();s++){if(samples[s].Name!=NameX){samples.erase(samples.begin()+s);s--;} }
}  





//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Genertic code related to samples processing in FWLITE --> functions below will be loaded only if FWLITE compiler variable is defined

#ifdef FWLITE

unsigned long GetInitialNumberOfMCEvent(const vector<string>& fileNames);
double GetSampleWeight  (const double& IntegratedLuminosityInPb=-1, const double& IntegratedLuminosityInPbBeforeTriggerChange=-1, const double& CrossSection=0, const double& MCEvents=0, int period=0);
double GetSampleWeightMC(const double& IntegratedLuminosityInPb, const std::vector<string> fileNames, const double& XSection, const double& SampleSize, double MaxEvent);
double GetPUWeight      (const fwlite::ChainEvent& ev, const std::string& pileup, double &PUSystFactor);

// loop on all the lumi blocks for an EDM file in order to count the number of events that are in a sample
// this is useful to determine how to normalize the events (compute weight)
unsigned long GetInitialNumberOfMCEvent(const vector<string>& fileNames)
{
   unsigned long Total = 0;
   for(unsigned int f=0;f<fileNames.size();f++){
      TFile *file;
      size_t place=fileNames[f].find("dcache");
      if(place!=string::npos){
 	 string name=fileNames[f];
         name.replace(place, 7, "dcap://cmsgridftp.fnal.gov:24125");
         file = new TDCacheFile (name.c_str());
      }else{
         file = TFile::Open(fileNames[f].c_str());
      }

      fwlite::LuminosityBlock ls( file );
      for(ls.toBegin(); !ls.atEnd(); ++ls){
         fwlite::Handle<edm::MergeableCounter> nEventsTotalCounter;
         nEventsTotalCounter.getByLabel(ls,"nEventsBefSkim");
         if(!nEventsTotalCounter.isValid()){printf("Invalid nEventsTotalCounterH\n");continue;}
         Total+= nEventsTotalCounter->value;
      }
   }
   return Total;
}

// compute event weight for the signal samples based on number of events in the sample, cross section and intergated luminosity
double GetSampleWeight(const double& IntegratedLuminosityInPb, const double& IntegratedLuminosityInPbBeforeTriggerChange, const double& CrossSection, const double& MCEvents, int period){
  double Weight = 1.0;
  if(IntegratedLuminosityInPb>=IntegratedLuminosityInPbBeforeTriggerChange && IntegratedLuminosityInPb>0){
    double NMCEvents = MCEvents;
    //if(MaxEntry>0)NMCEvents=std::min(MCEvents,(double)MaxEntry);
    if      (period==0)Weight = (CrossSection * IntegratedLuminosityInPbBeforeTriggerChange) / NMCEvents;
    else if (period==1)Weight = (CrossSection * (IntegratedLuminosityInPb-IntegratedLuminosityInPbBeforeTriggerChange)) / NMCEvents;
  }
  return Weight;
}

// compute event weight for the MC background samples based on number of events in the sample, cross section and intergated luminosity
double GetSampleWeightMC(const double& IntegratedLuminosityInPb, const std::vector<string> fileNames, const double& XSection, const double& SampleSize, double MaxEvent){
  double Weight = 1.0;
   unsigned long InitNumberOfEvents = GetInitialNumberOfMCEvent(fileNames); 
   double SampleEquivalentLumi = InitNumberOfEvents / XSection;
   if(MaxEvent<0)MaxEvent=SampleSize;
   printf("GetSampleWeight MC: IntLumi = %6.2E  SampleLumi = %6.2E --> EventWeight = %6.2E --> ",IntegratedLuminosityInPb,SampleEquivalentLumi, IntegratedLuminosityInPb/SampleEquivalentLumi);
//   printf("Sample NEvent = %6.2E   SampleEventUsed = %6.2E --> Weight Rescale = %6.2E\n",SampleSize, MaxEvent, SampleSize/MaxEvent);
   Weight = (IntegratedLuminosityInPb/SampleEquivalentLumi) * (SampleSize/MaxEvent);
   printf("FinalWeight = %6.2f\n",Weight);
   return Weight;
}

// compute a weight to make sure that the pileup distribution in Signal/Background MC is compatible to the pileup distribution in data
double GetPUWeight(const fwlite::ChainEvent& ev, const std::string& pileup, double &PUSystFactor, edm::LumiReWeighting& LumiWeightsMC, reweight::PoissonMeanShifter& PShift){
   fwlite::Handle<std::vector<PileupSummaryInfo> > PupInfo;
   PupInfo.getByLabel(ev, "addPileupInfo");
   if(!PupInfo.isValid()){printf("PileupSummaryInfo Collection NotFound\n");return 1.0;}
   double PUWeight_thisevent=1;
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   int npv = -1;
   if(pileup=="S4"){
      float sum_nvtx = 0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         npv = PVI->getPU_NumInteractions();
         sum_nvtx += float(npv);
      }
      float ave_nvtx = sum_nvtx/3.;
      PUWeight_thisevent = LumiWeightsMC.weight( ave_nvtx );
      PUSystFactor = PShift.ShiftWeight( ave_nvtx );
   }else if(pileup=="S3"){
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         int BX = PVI->getBunchCrossing();
         if(BX == 0) {
            npv = PVI->getPU_NumInteractions();
            continue;
         }
      }
      PUWeight_thisevent = LumiWeightsMC.weight( npv );
      PUSystFactor = PShift.ShiftWeight( npv );
   }
   return PUWeight_thisevent;
}

#endif //end FWLITE block



//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define Theoretical XSection and error band for 7TeVStops and 7TeVGluino

//NLO Stop@7TeV
double THXSEC7TeV_Stop_Mass [] = {100.0,120.0,140.0,160.0,180.0,200.0,220.0,240.0,260.0,280.0,300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0
,390.0,165.0,76.8,38.8,21.1,12.1,7.21,4.46,2.84,1.86,1.24,0.844,0.585,0.412,0.294,0.212,0.155,0.114,0.0847,0.0634,0.0479,0.0365,0.0279,0.0215,0.0166,0.0129,0.0101,0.00796,0.0063,0.00499,0.00397,0.00317,0.00253,0.00203,0.00164,0.00132,0.00107,0.000864,0.000702,0.000572,0.000467,0.000381,0.000312,0.000255,0.00021,0.000173 };
double THXSEC7TeV_Stop_Cen  [] = {390.0,165.0,76.8,38.8,21.1,12.1,7.21,4.46,2.84,1.86,1.24,0.844,0.585,0.412,0.294,0.212,0.155,0.114,0.0847,0.0634,0.0479,0.0365,0.0279,0.0215,0.0166,0.0129,0.0101,0.00796,0.0063,0.00499,0.00397,0.00317,0.00253,0.00203,0.00164,0.00132,0.00107,0.000864,0.000702,0.000572,0.000467,0.000381,0.000312,0.000255,0.00021,0.000173 };
double THXSEC7TeV_Stop_Low  [] = {336.687212239,142.161704016,65.8513147385,33.1253989627,17.918668408,10.2276057109,6.06740619265,3.73576507352,2.36659235722,1.54221509094,1.02378918186,0.693358458913,0.478780586501,0.335646071783,0.238501707841,0.171183802392,0.124608032654,0.0912510120893,0.0674833984023,0.0502792202566,0.037816978755,0.0286883112529,0.0218301431028,0.0167475917045,0.012877462018,0.00995875547206,0.00775938400598,0.00608185167256,0.00479071471604,0.00377411064367,0.00298741802917,0.00237365916816,0.00188189555314,0.00150193564851,0.00120631066392,0.000965374807573,0.000778619696638,0.000624245761673,0.000504320503585,0.00040804137369,0.000330700238198,0.000267790633317,0.000217758112762,0.000176730791724,0.000144344874883,0.000117997311653 };
double THXSEC7TeV_Stop_High [] = {455.337270601,192.484749389,89.6237634573,45.38,24.7104320442,14.2154509779,8.49882260582,5.27246023469,3.3687269963,2.21436345753,1.48113888098,1.01267602748,0.704377399955,0.497929290742,0.356909806888,0.258306503513,0.189623679693,0.14001964183,0.104443008835,0.0785510793723,0.059587810898,0.0455979172605,0.035030746117,0.0271154457368,0.021046423859,0.0164444810585,0.0129238101605,0.0102301624614,0.00812681856605,0.00646997801663,0.00517242916474,0.004153874823,0.00333274574668,0.00268995484779,0.00218320499759,0.00176632187855,0.00143974932708,0.00116964830458,0.000955623561151,0.000783143384745,0.000643472508176,0.000528336426727,0.000435441662818,0.000358165927707,0.000296784169801,0.000246047049368 };

//NLO Stop@8TeV
double THXSEC8TeV_Stop_Mass [] = { 300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0,1020.0,1040.0,1060.0,1080.0,1100.0,1120.0,1140.0,1160.0,1180.0,1200.0,1220.0,1240.0,1260.0,1280.0,1300.0,1320.0,1340.0,1360.0,1380.0,1400.0,1420.0,1440.0,1460.0,1480.0,1500.0,1520.0,1540.0,1560.0,1580.0,1600.0,1620.0,1640.0,1660.0,1680.0,1700.0,1720.0,1740.0,1760.0,1780.0,1800.0,1820.0,1840.0,1860.0,1880.0,1900.0,1920.0,1940.0,1960.0,1980.0 };
double THXSEC8TeV_Stop_Cen  [] = { 103.0,70.6,49.3,35.1,25.3,18.5,13.6,10.2,7.66,5.82,4.46,3.44,2.66,2.08,1.63,1.29,1.02,0.816,0.653,0.525,0.424,0.344,0.279,0.228,0.186,0.153,0.126,0.104,0.0856,0.0709,0.0588,0.0489,0.0408,0.0341,0.0285,0.0239,0.02,0.0168,0.0141,0.0119,0.01,0.00845,0.00714,0.00605,0.00513,0.00435,0.00369,0.00314,0.00267,0.00227,0.00193,0.00164,0.0014,0.00119,0.00102,0.000866,0.000739,0.000631,0.000539,0.00046,0.000393,0.000336,0.000287,0.000245,0.000209,0.000179,0.000153,0.000131,0.000112,9.52e-05,8.13e-05,6.94e-05,5.93e-05,5.06e-05,4.31e-05,3.68e-05,3.14e-05,2.67e-05,2.28e-05,1.94e-05,1.65e-05,1.4e-05,1.19e-05,1.02e-05,8.63e-06};
double THXSEC8TeV_Stop_Low  [] = { 85.4522933503,58.3024212534,40.5404072161,28.7311959371,20.6110330718,14.9981486143,10.9817116552,8.2045980396,6.13599856717,4.64639599616,3.54679219469,2.7264626978,2.09915349508,1.63562861604,1.27699308705,1.00685323305,0.792356354919,0.632441971506,0.503862668979,0.403411581822,0.324940255218,0.262909295502,0.212319587364,0.173005424976,0.140570743272,0.115316873907,0.0947326706086,0.0778816775066,0.0638144532538,0.0526220125758,0.0434789865444,0.0359951306958,0.029944567911,0.0249059836578,0.020754521583,0.0173207451689,0.0144267991917,0.0120799394065,0.0100863120603,0.00847183872318,0.00708717142879,0.00595886341748,0.00501035858415,0.00422415500896,0.00356494485684,0.00300727877119,0.00253860291188,0.00214977087027,0.00181901909611,0.00153605196417,0.00129944484705,0.00109821151599,0.000931476154739,0.000787499214937,0.000671494422707,0.000566354149908,0.000479894886148,0.000407420453834,0.000345671282251,0.000292733637163,0.000248262306276,0.000210814299455,0.000178718470665,0.000151342195875,0.000127869581852,0.000108804207517,9.23062817067e-05,7.82498615939e-05,6.63400478633e-05,5.58615231572e-05,4.72597016726e-05,3.99537211035e-05,3.38756321908e-05,2.86208917333e-05,2.41297344895e-05,2.0403820875e-05,1.72341338273e-05,1.44941525431e-05,1.22590511119e-05,1.03016340378e-05,8.66483254755e-06,7.26134469144e-06,6.10042739754e-06,5.16703392541e-06,4.32379602699e-06 };
double THXSEC8TeV_Stop_High [] = { 123.042915069,84.6204861909,59.3621530586,42.4270342231,30.720245281,22.5429742994,16.6400543991,12.521,9.44012836998,7.19888313862,5.53783839898,4.28718343406,3.32951446035,2.61303391779,2.0546626767,1.63189165427,1.29573309725,1.03959612686,0.835270344823,0.673688637441,0.546222752456,0.444771492296,0.362012997158,0.297010039438,0.243115317971,0.200703894146,0.165936031728,0.137421289844,0.113579295087,0.094359434471,0.0785768709139,0.0656067780802,0.054914187159,0.0460700308168,0.0386624800886,0.0325504144943,0.0273588316481,0.0230737980716,0.019440405649,0.0164822382849,0.0139146663585,0.0118012857472,0.0100182943312,0.00852079529061,0.00725778032384,0.00618194507243,0.00526737821502,0.00450197140212,0.00384501422077,0.00328618755105,0.00280677982274,0.00239852281744,0.00205691112089,0.00175754634874,0.00151421701532,0.00129260450299,0.00110862278184,0.000951522205324,0.000817402023602,0.000701118190558,0.000602322798422,0.000517825498222,0.000444861065588,0.000382177349458,0.000327902228753,0.00028253307733,0.000243007105712,0.00020936089293,0.000180051621273,0.000154016576997,0.000132434560814,0.000113806344946,9.78639305045e-05,8.40542824159e-05,7.21176753055e-05,6.20064018569e-05,5.32790981711e-05,4.56654057754e-05,3.92622063114e-05,3.3654042914e-05,2.88430383571e-05,2.46542444561e-05,2.11238854891e-05,1.82236973231e-05,1.55451905173e-05 };

//NLO Gluino@7TeV
double THXSEC7TeV_Gluino_Mass [] = {300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0,1020.0,1040.0,1060.0,1080.0,1100.0,1120.0,1140.0,1160.0,1180.0,1200.0,1220.0,1240.0,1260.0,1280.0,1300.0,1320.0,1340.0,1360.0,1380.0,1400.0,1420.0,1440.0,1460.0,1480.0,1500.0,1520.0,1540.0,1560.0,1580.0,1600.0,1620.0,1640.0,1660.0,1680.0,1700.0,1720.0,1740.0,1760.0,1780.0,1800.0,1820.0,1840.0,1860.0,1880.0,1900.0,1920.0,1940.0,1960.0,1980.0 };
double THXSEC7TeV_Gluino_Cen  [] = {65.8,44.8,31.0,21.7,15.5,11.2,8.15,6.01,4.47,3.36,2.54,1.93,1.48,1.14,0.888,0.693,0.543,0.428,0.339,0.269,0.214,0.172,0.137,0.111,0.0895,0.0725,0.0588,0.0479,0.0391,0.032,0.0262,0.0215,0.0176,0.0145,0.0119,0.00987,0.00815,0.00674,0.00559,0.00464,0.00386,0.0032,0.00266,0.00222,0.00185,0.00154,0.00128,0.00107,0.000893,0.000746,0.000623,0.00052,0.000435,0.000363,0.000304,0.000254,0.000212,0.000177,0.000148,0.000124,0.000103,8.64e-05,7.22e-05,6.02e-05,5.02e-05,4.18e-05,3.49e-05,2.9e-05,2.42e-05,2.02e-05,1.67e-05,1.39e-05,1.16e-05,9.59e-06,7.95e-06,6.59e-06,5.45e-06,4.51e-06,3.73e-06,3.08e-06,2.54e-06,2.09e-06,1.72e-06,1.42e-06,1.16e-06 };
double THXSEC7TeV_Gluino_Low  [] = {53.9915787502,36.3679605186,25.04773106,17.5484531801,12.5029969909,8.9137230832,6.48814044136,4.71554753537,3.51079991399,2.63456118649,1.96607018967,1.50525499617,1.14834208312,0.886682415379,0.67892705475,0.528365261709,0.415476340911,0.322928423644,0.255343774802,0.202059219533,0.161164211822,0.127363153734,0.102023478594,0.0818182928751,0.065393518095,0.0529567854047,0.0430807023909,0.0345084584098,0.0280914623477,0.0229311955562,0.0184911518907,0.0152076608485,0.0124075251946,0.0102192110979,0.00832945628418,0.00686771000809,0.00568217384168,0.00462615538234,0.00382900127407,0.00313534410297,0.00259805136043,0.00215484050605,0.00180048617066,0.001467341696,0.00122061808,0.00100694617882,0.000838421855019,0.000684769467067,0.000571644050552,0.000476594263719,0.000392527368983,0.000326777849928,0.000268212922309,0.000224036696952,0.000184330609415,0.0001530784447,0.000126064184114,0.000105014721299,8.62123946791e-05,7.21482688124e-05,5.88983434237e-05,4.93030990383e-05,4.04422743416e-05,3.31101477016e-05,2.75456714875e-05,2.25762430633e-05,1.87694070612e-05,1.53040540155e-05,1.24350948334e-05,1.04092345492e-05,8.49269279945e-06,6.97519707277e-06,5.71827736139e-06,4.62718381962e-06,3.82358391349e-06,3.10935813704e-06,2.50561363539e-06,2.0258883204e-06,1.67348322755e-06,1.34494877209e-06,1.08865400522e-06,8.98317665716e-07,7.1490062492e-07,5.71185862116e-07,4.67630704264e-07 };
double THXSEC7TeV_Gluino_High [] = {79.5668288977,54.2211352869,37.5922538006,26.5343616416,19.075,13.786,10.1121682798,7.45102358218,5.57973835763,4.22381192237,3.19868564961,2.44746309321,1.88138665258,1.45913982451,1.13478647663,0.892932349164,0.704893627141,0.556473798188,0.440065907012,0.352253724644,0.282471833163,0.226566443289,0.182846908818,0.147392902942,0.119981370949,0.0972614292036,0.0794933605454,0.0647760279165,0.0531432044804,0.0435131063759,0.0359897509829,0.0297064480657,0.0243363202323,0.0201945890666,0.0166888696477,0.013839213608,0.0115485526692,0.00955917250779,0.00797993717028,0.00667993235001,0.00555403311772,0.00464008175518,0.00388720169471,0.00326472742026,0.00273449708667,0.00230578041856,0.00191402087824,0.00161285546003,0.00135345358155,0.00113904856138,0.000958125352577,0.000805385272244,0.000678241538475,0.000570293065787,0.000481037075751,0.000404734213278,0.000340132555915,0.000285662391546,0.000240024886815,0.000203182671651,0.000169737973538,0.000143208791006,0.000120394572153,0.000101748910812,8.53401846951e-05,7.16178115223e-05,6.0214226498e-05,5.05737327547e-05,4.23835072134e-05,3.58341122051e-05,2.98807306886e-05,2.50681230531e-05,2.11198589408e-05,1.76419321064e-05,1.47213338727e-05,1.23406750758e-05,1.02636214127e-05,8.59048670095e-06,7.18944995696e-06,5.96042280205e-06,4.97224153032e-06,4.14424242035e-06,3.4282369471e-06,2.85290331058e-06,2.35628618585e-06 };

//NLO Gluino@8TeV
double THXSEC8TeV_Gluino_Mass [] = {300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0,1020.0,1040.0,1060.0,1080.0,1100.0,1120.0,1140.0,1160.0,1180.0,1200.0,1220.0,1240.0,1260.0,1280.0,1300.0,1320.0,1340.0,1360.0,1380.0,1400.0,1420.0,1440.0,1460.0,1480.0,1500.0,1520.0,1540.0,1560.0,1580.0,1600.0,1620.0,1640.0,1660.0,1680.0,1700.0,1720.0,1740.0,1760.0,1780.0,1800.0,1820.0,1840.0,1860.0,1880.0,1900.0,1920.0,1940.0,1960.0,1980.0};
double THXSEC8TeV_Gluino_Cen  [] = {103.0,70.6,49.3,35.1,25.3,18.5,13.6,10.2,7.66,5.82,4.46,3.44,2.66,2.08,1.63,1.29,1.02,0.816,0.653,0.525,0.424,0.344,0.279,0.228,0.186,0.153,0.126,0.104,0.0856,0.0709,0.0588,0.0489,0.0408,0.0341,0.0285,0.0239,0.02,0.0168,0.0141,0.0119,0.01,0.00845,0.00714,0.00605,0.00513,0.00435,0.00369,0.00314,0.00267,0.00227,0.00193,0.00164,0.0014,0.00119,0.00102,0.000866,0.000739,0.000631,0.000539,0.00046,0.000393,0.000336,0.000287,0.000245,0.000209,0.000179,0.000153,0.000131,0.000112,9.52e-05,8.13e-05,6.94e-05,5.93e-05,5.06e-05,4.31e-05,3.68e-05,3.14e-05,2.67e-05,2.28e-05,1.94e-05,1.65e-05,1.4e-05,1.19e-05,1.02e-05,8.63e-06};
double THXSEC8TeV_Gluino_Low  [] = {85.4522933503,58.3024212534,40.5404072161,28.7311959371,20.6110330718,14.9981486143,10.9817116552,8.2045980396,6.13599856717,4.64639599616,3.54679219469,2.7264626978,2.09915349508,1.63562861604,1.27699308705,1.00685323305,0.792356354919,0.632441971506,0.503862668979,0.403411581822,0.324940255218,0.262909295502,0.212319587364,0.173005424976,0.140570743272,0.115316873907,0.0947326706086,0.0778816775066,0.0638144532538,0.0526220125758,0.0434789865444,0.0359951306958,0.029944567911,0.0249059836578,0.020754521583,0.0173207451689,0.0144267991917,0.0120799394065,0.0100863120603,0.00847183872318,0.00708717142879,0.00595886341748,0.00501035858415,0.00422415500896,0.00356494485684,0.00300727877119,0.00253860291188,0.00214977087027,0.00181901909611,0.00153605196417,0.00129944484705,0.00109821151599,0.000931476154739,0.000787499214937,0.000671494422707,0.000566354149908,0.000479894886148,0.000407420453834,0.000345671282251,0.000292733637163,0.000248262306276,0.000210814299455,0.000178718470665,0.000151342195875,0.000127869581852,0.000108804207517,9.23062817067e-05,7.82498615939e-05,6.63400478633e-05,5.58615231572e-05,4.72597016726e-05,3.99537211035e-05,3.38756321908e-05,2.86208917333e-05,2.41297344895e-05,2.0403820875e-05,1.72341338273e-05,1.44941525431e-05,1.22590511119e-05,1.03016340378e-05,8.66483254755e-06,7.26134469144e-06,6.10042739754e-06,5.16703392541e-06,4.32379602699e-06};
double THXSEC8TeV_Gluino_High [] = {123.042915069,84.6204861909,59.3621530586,42.4270342231,30.720245281,22.5429742994,16.6400543991,12.521,9.44012836998,7.19888313862,5.53783839898,4.28718343406,3.32951446035,2.61303391779,2.0546626767,1.63189165427,1.29573309725,1.03959612686,0.835270344823,0.673688637441,0.546222752456,0.444771492296,0.362012997158,0.297010039438,0.243115317971,0.200703894146,0.165936031728,0.137421289844,0.113579295087,0.094359434471,0.0785768709139,0.0656067780802,0.054914187159,0.0460700308168,0.0386624800886,0.0325504144943,0.0273588316481,0.0230737980716,0.019440405649,0.0164822382849,0.0139146663585,0.0118012857472,0.0100182943312,0.00852079529061,0.00725778032384,0.00618194507243,0.00526737821502,0.00450197140212,0.00384501422077,0.00328618755105,0.00280677982274,0.00239852281744,0.00205691112089,0.00175754634874,0.00151421701532,0.00129260450299,0.00110862278184,0.000951522205324,0.000817402023602,0.000701118190558,0.000602322798422,0.000517825498222,0.000444861065588,0.000382177349458,0.000327902228753,0.00028253307733,0.000243007105712,0.00020936089293,0.000180051621273,0.000154016576997,0.000132434560814,0.000113806344946,9.78639305045e-05,8.40542824159e-05,7.21176753055e-05,6.20064018569e-05,5.32790981711e-05,4.56654057754e-05,3.92622063114e-05,3.3654042914e-05,2.88430383571e-05,2.46542444561e-05,2.11238854891e-05,1.82236973231e-05,1.55451905173e-05};



//LO GMSB Stau@7TeV
double THXSEC7TeV_GMStau_Mass [] = {100,126,156,200,247,308,370,432,494};
double THXSEC7TeV_GMStau_Cen  [] = {1.3398,0.274591,0.0645953,0.0118093,0.00342512,0.00098447,0.000353388,0.000141817,6.17749e-05};
double THXSEC7TeV_GMStau_Low  [] = {1.18163,0.242982,0.0581651,0.0109992,0.00324853,0.00093519,0.000335826,0.000134024,5.83501e-05};
double THXSEC7TeV_GMStau_High [] = {1.48684,0.304386,0.0709262,0.012632,0.00358232,0.00102099,0.000366819,0.000147665,6.45963e-05};

//LO PairProduced Stau@7TeV
double THXSEC7TeV_PPStau_Mass [] = {100,126,156,200,247,308};
double THXSEC7TeV_PPStau_Cen  [] = {0.038200,0.0161,0.007040,0.002470,0.001010,0.000353};
double THXSEC7TeV_PPStau_Low  [] = {0.037076,0.0155927,0.0067891,0.00237277,0.00096927,0.000335308};
double THXSEC7TeV_PPStau_High [] = {0.0391443,0.016527,0.00723151,0.00253477,0.00103844,0.000363699};

#endif

