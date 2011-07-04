
#ifndef HSCP_ANALYSIS_SAMPLE
#define HSCP_ANALYSIS_SAMPLE


#define SID_GL200     0
#define SID_GL300     1
#define SID_GL400     2
#define SID_GL500     3
#define SID_GL600     4
#define SID_GL900     5
#define SID_GL200N    6
#define SID_GL300N    7
#define SID_GL400N    8
#define SID_GL500N    9
#define SID_GL600N   10
#define SID_GL900N   11
#define SID_ST130    12
#define SID_ST200    13
#define SID_ST300    14
#define SID_ST500    15
#define SID_ST800    16
#define SID_ST130N   17
#define SID_ST200N   18
#define SID_ST300N   19
#define SID_ST500N   20
#define SID_ST800N   21
#define SID_GS100    22
#define SID_GS126    23
#define SID_GS156    24
#define SID_GS200    25
#define SID_GS247    26
#define SID_GS308    27
#define SID_PS100    28
#define SID_PS126    29
#define SID_PS156    30
#define SID_PS200    31
#define SID_PS247    32
#define SID_PS308    33
#define SID_DS121    34
#define SID_DS182    35
#define SID_DS242    36
#define SID_DS302    37


int                  RunningPeriods = 2;
double               IntegratedLuminosity = 802.112979531; //715.041013013; //705.273820; //342.603275; //204.160928; //191.04;
double               IntegratedLuminosityBeforeTriggerChange = 324.787493137; // Total luminosity taken before RPC L1 trigger change (went into effect on run 165970)
float                Event_Weight = 1;
int                  MaxEntry = -1;


class stSignal{
   public:
   std::string Type;
   std::string Name;
   std::string Legend;
   double Mass;
   double XSec;
   bool   MakePlot;

   stSignal(); 
   stSignal(std::string Type_, std::string Name_, std::string Legend_, double Mass_, bool MakePlot_, double XSec_){Type=Type_; Name=Name_; Legend=Legend_; Mass=Mass_; MakePlot=MakePlot_; XSec=XSec_;}
};


void GetSignalDefinition(std::vector<stSignal>& signals){

// signals.push_back(stSignal("Gluino", "Gluino200"    , "#tilde{g} 200"                 , 200,  1, 606.000000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino300"    , "#tilde{g} 300"                 , 300,  1,  57.200000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino400"    , "#tilde{g} 400"                 , 400,  1,   8.980000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino500"    , "#tilde{g} 500"                 , 500,  1,   1.870000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino600"    , "#tilde{g} 600"                 , 600,  1,   0.465000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino700"    , "#tilde{g} 700"                 , 700,  1,   0.130000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino800"    , "#tilde{g} 800"                 , 800,  1,   0.039600) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino900"    , "#tilde{g} 900"                 , 900,  1,   0.012800) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino1000"   , "#tilde{g} 1000"                ,1000,  1,   0.004300) ); //NLO

   signals.push_back(stSignal("Gluino", "Gluino300N"   , "#tilde{g} 300 CS"              , 300,  1,  57.200000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino400N"   , "#tilde{g} 400 CS"              , 400,  1,   8.980000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino500N"   , "#tilde{g} 500 CS"              , 500,  1,   1.870000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino600N"   , "#tilde{g} 600 CS"              , 600,  1,   0.465000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino700N"   , "#tilde{g} 700 CS"              , 700,  1,   0.130000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino800N"   , "#tilde{g} 800 CS"              , 800,  1,   0.039600) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino900N"   , "#tilde{g} 900 CS"              , 900,  1,   0.012800) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino1000N"  , "#tilde{g} 1000 CS"             ,1000,  1,   0.004300) ); //NLO

   signals.push_back(stSignal("Gluino", "Gluino600Z"   , "#tilde{g} 600 Z2"              , 600,  1,   0.465000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino700Z"   , "#tilde{g} 700 Z2"              , 700,  1,   0.130000) ); //NLO
   signals.push_back(stSignal("Gluino", "Gluino800Z"   , "#tilde{g} 800 Z2"              , 800,  1,   0.039600) ); //NLO

   signals.push_back(stSignal("Stop"  , "Stop130"      , "#tilde{t}_{1} 130"             , 130,  1, 120.000000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop200"      , "#tilde{t}_{1} 200"             , 200,  1,  13.000000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop300"      , "#tilde{t}_{1} 300"             , 300,  1,   1.310000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop400"      , "#tilde{t}_{1} 400"             , 400,  1,   0.218000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop500"      , "#tilde{t}_{1} 500"             , 500,  0,   0.047800) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop600"      , "#tilde{t}_{1} 600"             , 600,  1,   0.012500) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop700"      , "#tilde{t}_{1} 700"             , 700,  1,   0.003560) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop800"      , "#tilde{t}_{1} 800"             , 800,  1,   0.001140) ); //NLO

   signals.push_back(stSignal("Stop"  , "Stop130N"     , "#tilde{t}_{1} 130 CS"          , 130,  1, 120.000000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop200N"     , "#tilde{t}_{1} 200 CS"          , 200,  1,  13.000000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop300N"     , "#tilde{t}_{1} 300 CS"          , 300,  1,   1.310000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop400N"     , "#tilde{t}_{1} 400 CS"          , 400,  1,   0.218000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop500N"     , "#tilde{t}_{1} 500 CS"          , 500,  0,   0.047800) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop600N"     , "#tilde{t}_{1} 600 CS"          , 600,  1,   0.012500) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop700N"     , "#tilde{t}_{1} 700 CS"          , 700,  1,   0.003560) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop800N"     , "#tilde{t}_{1} 800 CS"          , 800,  1,   0.001140) ); //NLO

   signals.push_back(stSignal("Stop"  , "Stop300Z"     , "#tilde{t}_{1} 300 Z2"          , 300,  1,   1.310000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop400Z"     , "#tilde{t}_{1} 400 Z2"          , 400,  1,   0.218000) ); //NLO
   signals.push_back(stSignal("Stop"  , "Stop500Z"     , "#tilde{t}_{1} 500 Z2"          , 500,  0,   0.047800) ); //NLO

   signals.push_back(stSignal("Stau"  , "GMStau100"    , "GMSB #tilde{#tau}_{1} 100"     , 100,  1,   1.326000) ); //LO
   signals.push_back(stSignal("Stau"  , "GMStau126"    , "GMSB #tilde{#tau}_{1} 126"     , 126,  1,   0.330000) ); //LO
   signals.push_back(stSignal("Stau"  , "GMStau156"    , "GMSB #tilde{#tau}_{1} 156"     , 156,  0,   0.105000) ); //LO
   signals.push_back(stSignal("Stau"  , "GMStau200"    , "GMSB #tilde{#tau}_{1} 200"     , 200,  1,   0.025000) ); //LO
   signals.push_back(stSignal("Stau"  , "GMStau247"    , "GMSB #tilde{#tau}_{1} 247"     , 247,  0,   0.008000) ); //LO
   signals.push_back(stSignal("Stau"  , "GMStau308"    , "GMSB #tilde{#tau}_{1} 308"     , 308,  1,   0.002000) ); //LO
/*   signals.push_back(stSignal("Stau"  , "PPStau100"    , "Pair #tilde{#tau}_{1} 100"     , 100,  1,   0.032000) ); //LO
   signals.push_back(stSignal("Stau"  , "PPStau126"    , "Pair #tilde{#tau}_{1} 126"     , 126,  0,   0.014000) ); //LO
   signals.push_back(stSignal("Stau"  , "PPStau156"    , "Pair #tilde{#tau}_{1} 156"     , 156,  0,   0.006000) ); //LO
   signals.push_back(stSignal("Stau"  , "PPStau200"    , "Pair #tilde{#tau}_{1} 200"     , 200,  1,   0.002200) ); //LO
   signals.push_back(stSignal("Stau"  , "PPStau247"    , "Pair #tilde{#tau}_{1} 247"     , 247,  0,   0.000900) ); //LO
   signals.push_back(stSignal("Stau"  , "PPStau308"    , "Pair #tilde{#tau}_{1} 308"     , 308,  1,   0.000300) ); //LO
   signals.push_back(stSignal("Stau"  , "DCStau121"    , "DICHAMP #tilde{#tau}_{1} 121"  , 121,  1,   0.450000) ); //LO
   signals.push_back(stSignal("Stau"  , "DCStau182"    , "DICHAMP #tilde{#tau}_{1} 182"  , 182,  0,   0.083000) ); //LO
   signals.push_back(stSignal("Stau"  , "DCStau242"    , "DICHAMP #tilde{#tau}_{1} 242"  , 242,  0,   0.022800) ); //LO
   signals.push_back(stSignal("Stau"  , "DCStau302"    , "DICHAMP #tilde{#tau}_{1} 302"  , 302,  1,   0.007700) ); //LO
*/
}

struct stMC{
   std::string Name;
   double ILumi;
   double MaxPtHat;
   double MaxEvent;

   stMC();
   stMC(std::string Name_, double ILumi_, double MaxPtHat_, int MaxEvent_){Name = Name_; ILumi = ILumi_; MaxPtHat = MaxPtHat_; MaxEvent = MaxEvent_;}
};

void GetMCDefinition(std::vector<stMC>& MC){
//   MC.push_back(stMC("MC_MB"   , 0.000754    , 30, -1 ));
//   MC.push_back(stMC("MC_QCD30", 0.07708     , 80, -1 ));
//   MC.push_back(stMC("MC_QCD80", 3.1700      , -1, -1 ));
}


void GetInputFiles(std::vector<std::string>& inputFiles, std::string SampleName, int period=0){
  //std::string BaseDirectory = "/storage/data/cms/users/quertenmont/HSCP/CMSSW_4_2_3/11_06_28/";
   std::string BaseDirectory = "dcap://cmsdca.fnal.gov:24125/pnfs/cms/WAX/11/store/user/farrell3/EDMFiles/";

   if(SampleName=="Data"){
        inputFiles.push_back(BaseDirectory + "Data_RunA_160000_163250.root");
        inputFiles.push_back(BaseDirectory + "Data_RunA_163251_163500.root");
        inputFiles.push_back(BaseDirectory + "Data_RunA_165001_166033.root");                                                                                
        inputFiles.push_back(BaseDirectory + "Data_RunA_166034_166500.root"); 
        inputFiles.push_back(BaseDirectory + "Data_RunA_166501_166893.root");
        inputFiles.push_back(BaseDirectory + "Data_RunA_166894_167151.root");

//   }else if(SampleName=="MC_MB"){
//      inputFiles.push_back(BaseDirectory + "MC_MB.root");
//   }else if(SampleName=="MC_PPMUX"){
//      inputFiles.push_back(BaseDirectory + "MC_PPMUX.root");
//   }else if(SampleName=="MC_QCD30"){
//      inputFiles.push_back(BaseDirectory + "MC_QCD30.root");
//   }else if(SampleName=="MC_QCD80"){
//      if(rand()%2==0){
//         inputFiles.push_back(BaseDirectory + "MC_QCD80.root");
//      }else{
//         inputFiles.push_back(BaseDirectory + "MC_QCD80_B.root");
//      }
   }else{
     if (period==0) inputFiles.push_back(BaseDirectory + SampleName + ".root");
     if (period==1) inputFiles.push_back(BaseDirectory + SampleName + "_RPCL1Extended.root");
   }
}

#endif
