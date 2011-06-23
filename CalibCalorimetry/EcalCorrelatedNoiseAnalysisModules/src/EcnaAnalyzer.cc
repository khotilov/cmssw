// -*- C++ -*-
//
// Package:    EcalCorrelatedNoiseAnalysisModules
// Class:      EcnaAnalyzer
// // class EcnaAnalyzer
// EcnaAnalyzer.cc
// CalibCalorimetry/EcalCorrelatedNoiseAnalysisModules/src/EcnaAnalyzer.cc

// Description: <one line class summary>

// Implementation:
//     <Notes on implementation>

//
// Original Author:  Bernard Fabbro
//         Created:  Fri Jun  2 10:27:01 CEST 2006
// $Id: EcnaAnalyzer.cc,v 1.2 2007/08/10 14:52:31 ferriff Exp $
//
//          Update: 01/04/2010  

// CMSSW include files

#include <signal.h>

#include "CalibCalorimetry/EcalCorrelatedNoiseAnalysisModules/interface/EcnaAnalyzer.h"
#include "CalibCalorimetry/EcalCorrelatedNoiseAnalysisAlgos/interface/TEcnaParEcal.h"
#include <DataFormats/Provenance/interface/Timestamp.h>

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EcnaAnalyzer::EcnaAnalyzer(const edm::ParameterSet& pSet) : 
  verbosity_(pSet.getUntrackedParameter("verbosity", 1U)),
  nChannels_(0), iEvent_(0)
{
  //now do what ever initialization is needed

  using namespace edm;
  using namespace std;

  std::cout << "*EcnaAnalyzer-constructor> Parameter initialization." << endl;

  fgMaxCar = (Int_t)512;
  fTTBELL = '\007';
  fOutcomeError = kFALSE;

  fMyEBEcal      = new TEcnaParEcal("EB");
  fMyEBNumbering = new TEcnaNumbering("EB", fMyEBEcal);

  fMyEEEcal      = new TEcnaParEcal("EE");
  fMyEENumbering = new TEcnaNumbering("EE", fMyEEEcal);

  //========================================================================================== 
  //.................................. Get parameter values from python file
  eventHeaderProducer_   = pSet.getParameter<std::string>("eventHeaderProducer");
  digiProducer_          = pSet.getParameter<std::string>("digiProducer");

  eventHeaderCollection_ = pSet.getParameter<std::string>("eventHeaderCollection");
  EBdigiCollection_      = pSet.getParameter<std::string>("EBdigiCollection");
  EEdigiCollection_      = pSet.getParameter<std::string>("EEdigiCollection");

  sAnalysisName_  = pSet.getParameter<std::string>("sAnalysisName");
  sNbOfSamples_   = pSet.getParameter<std::string>("sNbOfSamples");
  sFirstReqEvent_ = pSet.getParameter<std::string>("sFirstReqEvent");
  sLastReqEvent_  = pSet.getParameter<std::string>("sLastReqEvent");
  sReqNbOfEvts_   = pSet.getParameter<std::string>("sReqNbOfEvts");
  sStexName_      = pSet.getParameter<std::string>("sStexName");
  sStexNumber_    = pSet.getParameter<std::string>("sStexNumber");

  fAnalysisName  = sAnalysisName_.Data();
  fNbOfSamples   = atoi(sNbOfSamples_.Data());
  fFirstReqEvent = atoi(sFirstReqEvent_.Data());
  fLastReqEvent  = atoi(sLastReqEvent_.Data());
  fReqNbOfEvts   = atoi(sReqNbOfEvts_.Data());
  fStexName      = sStexName_.Data();
  fStexNumber    = atoi(sStexNumber_.Data());

  //------------------------------- ERRORS in requested evts numbers
  if( fFirstReqEvent < 1 )
    {fOutcomeError = AnalysisOutcome("ERR_FNEG");}

  if( (fLastReqEvent >= fFirstReqEvent) && (fReqNbOfEvts > fLastReqEvent - fFirstReqEvent + 1 ) )
    {fOutcomeError = AnalysisOutcome("ERR_LREQ");}

  if( fOutcomeError == kTRUE )return;
  //===========================================================================================

  fRunTypeCounter = 0;
  fMaxRunTypeCounter = 25;
  fRunTypeCounter = new Int_t[fMaxRunTypeCounter];
  for(Int_t i=0; i<fMaxRunTypeCounter; i++){fRunTypeCounter[i] = 0;}

  fMgpaGainCounter = 0;
  fMaxMgpaGainCounter = 4;  // Because chozen gain = 0,1,2,3
  fMgpaGainCounter = new Int_t[fMaxMgpaGainCounter];
  for(Int_t i=0; i<fMaxMgpaGainCounter; i++){fMgpaGainCounter[i] = 0;}

  fFedIdCounter = 0;
  fMaxFedIdCounter = 54;
  fFedIdCounter = new Int_t[fMaxFedIdCounter];
  for(Int_t i=0; i<fMaxFedIdCounter; i++){fFedIdCounter[i] = 0;}

  fEvtNumber = 0;
  fEvtNumberMemo = -1;
  fRecNumber = 0;
  fBadBuildCnaRun = 0;

  fDeeDS5Memo1 = 0;
  fDeeDS5Memo2 = 0;

  fCurrentEventNumber = 0;
  fNbOfSelectedEvents = 0;

  fMemoCutOK = 0;
  fTreatedFedOrder = 0;
  fNbOfTreatedStexs = 0;

  //-------------- Fed
  if( fStexName == "SM"  ){fMaxFedUnitCounter = fMyEBEcal->MaxSMInEB();}  // EB: FED Unit = SM
  if( fStexName == "Dee" ){fMaxFedUnitCounter = fMyEEEcal->MaxDSInEE();}  // EE: FED Unit = Data Sector

  fFedDigiOK = new Int_t[fMaxFedUnitCounter];
  for(Int_t i=0; i<fMaxFedUnitCounter; i++){fFedDigiOK[i] = 0;}

  fFedNbOfTreatedEvents = new Int_t[fMaxFedUnitCounter];
  for(Int_t i=0; i<fMaxFedUnitCounter; i++){fFedNbOfTreatedEvents[i] = 0;}

  fFedStatus = new Int_t[fMaxFedUnitCounter];
  for(Int_t i=0; i<fMaxFedUnitCounter; i++){fFedStatus[i] = 0;}

  fFedStatusOrder = new Int_t[fMaxFedUnitCounter];
  for(Int_t i=0; i<fMaxFedUnitCounter; i++){fFedStatusOrder[i] = 0;}

  fDeeNumberString = new TString[fMaxFedUnitCounter];
  for(Int_t i=0; i<fMaxFedUnitCounter; i++){fDeeNumberString[i] = "SM";}

  if( fStexName == "Dee" )
    { 
      fDeeNumberString[ 0] = "Sector1 Dee4";
      fDeeNumberString[ 1] = "Sector2 Dee4";
      fDeeNumberString[ 2] = "Sector3 Dee4";
      fDeeNumberString[ 3] = "Sector4 Dee4";
      fDeeNumberString[ 4] = "Sector5 Dee4-Dee3";
      fDeeNumberString[ 5] = "Sector6 Dee3";
      fDeeNumberString[ 6] = "Sector7 Dee3";
      fDeeNumberString[ 7] = "Sector8 Dee3";
      fDeeNumberString[ 8] = "Sector9 Dee3";
      fDeeNumberString[ 9] = "Sector1 Dee1";
      fDeeNumberString[10] = "Sector2 Dee1";
      fDeeNumberString[11] = "Sector3 Dee1";
      fDeeNumberString[12] = "Sector4 Dee1";
      fDeeNumberString[13] = "Sector5 Dee1-Dee2";
      fDeeNumberString[14] = "Sector6 Dee2";
      fDeeNumberString[15] = "Sector7 Dee2";
      fDeeNumberString[16] = "Sector8 Dee2";
      fDeeNumberString[17] = "Sector9 Dee2";
    }
  //............................... arrays fSMFromFedDcc and fESFromFedTcc
  //
  //  FED-TCC:   1   2   3   4   5   6   7   8   9
  //      Dee:   3   3   3   4   4   4   4  4-3  3
  //       DS:   7   8   9   1   2   3   4   5   6
  //       ES:   7   8   9   1   2   3   4   5   6  (ES = DS)
  //
  //  FED-TCC:  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
  //       SM:  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36
  //       SM:  -1  -2  -3  -4  -5  -6  -7  -8  -9 -10 -11 -12 -13 -14 -15 -16 -17 -18
  //
  //  FED-TCC:  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45
  //       SM:   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
  //
  //  FED-TCC:  46  47  48  49  50  51  52  53  54
  //      Dee:   2   2   2   1   1   1   1  1-2  2
  //       DS:   7   8   9   1   2   3   4   5   6
  //       ES:  16  17  18  10  11  12  13  14  15  (ES = DS + 9)
  
  Int_t MaxSMAndDS = fMyEBEcal->MaxSMInEB()+fMyEEEcal->MaxDSInEE();

  fSMFromFedTcc = new Int_t[MaxSMAndDS];
  fESFromFedTcc = new Int_t[MaxSMAndDS];
  for(Int_t nFedTcc= 1; nFedTcc<=MaxSMAndDS; nFedTcc++){fESFromFedTcc[nFedTcc-1] = -1;}

  for(Int_t nFedTcc= 1;nFedTcc<= 3;nFedTcc++){fESFromFedTcc[nFedTcc-1]=nFedTcc+ 6;} // Dee3, ES 7,8,9
  for(Int_t nFedTcc= 4;nFedTcc<= 9;nFedTcc++){fESFromFedTcc[nFedTcc-1]=nFedTcc- 3;} // Dee4, ES 1,2,3,4,5; Dee3, DS 5,6

  for(Int_t nFedTcc=10;nFedTcc<=27;nFedTcc++){fSMFromFedTcc[nFedTcc-1]=nFedTcc+ 9;} // EB-  SM 19 to 36
  for(Int_t nFedTcc=28;nFedTcc<=45;nFedTcc++){fSMFromFedTcc[nFedTcc-1]=nFedTcc-27;} // EB+  SM  1 to 18

  for(Int_t nFedTcc=46;nFedTcc<=48;nFedTcc++){fESFromFedTcc[nFedTcc-1]=nFedTcc-30;} // Dee2, ES 16,17,18
  for(Int_t nFedTcc=49;nFedTcc<=54;nFedTcc++){fESFromFedTcc[nFedTcc-1]=nFedTcc-39;} // Dee1, ES 10,11,12,13,14; Dee2, ES 14,15

  //............................... Nb of treated events for "AdcPeg12" and "AdcSPeg12" analysis
  //-------------- Stex
  if( fStexName == "SM"  ){fMaxTreatedStexCounter = fMyEBEcal->MaxSMInEB();}  // EB: Stex = SM
  if( fStexName == "Dee" ){fMaxTreatedStexCounter = fMyEEEcal->MaxDeeInEE();} // EE: Stex = Dee

  fStexNbOfTreatedEvents = new Int_t[fMaxTreatedStexCounter];
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++){fStexNbOfTreatedEvents[i] = 0;}

  fTimeFirst = new time_t[fMaxTreatedStexCounter];
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++){fTimeFirst[i] = 0;}
  fTimeLast = new time_t[fMaxTreatedStexCounter];
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++){fTimeLast[i] = 0;}

  fMemoDateFirstEvent = new Int_t[fMaxTreatedStexCounter];;
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++){fMemoDateFirstEvent[i] = 0;}

  Int_t MaxCar = fgMaxCar; 
  fDateFirst = new TString[fMaxTreatedStexCounter];
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++)
    {
      fDateFirst[i].Resize(MaxCar);
      fDateFirst[i] = "*1st event date not found*";
    }

  MaxCar = fgMaxCar; 
  fDateLast = new TString[fMaxTreatedStexCounter];
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++)
    {
      fDateLast[i].Resize(MaxCar);
      fDateLast[i] = "*last event date not found*";
    }

  fStexStatus = new Int_t[fMaxTreatedStexCounter];
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++){fStexStatus[i] = 0;}

  fStexDigiOK = new Int_t[fMaxTreatedStexCounter];
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++){fStexDigiOK[i] = 0;}

  fNbOfTreatedFedsInDee = new Int_t[fMaxTreatedStexCounter];
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++){fNbOfTreatedFedsInDee[i] = 0;}

  fNbOfTreatedFedsInStex = new Int_t[fMaxTreatedStexCounter];
  for(Int_t i=0; i<fMaxTreatedStexCounter; i++){fNbOfTreatedFedsInStex[i] = 0;}

  //----------------------------------- Analysis name codes ------------------------------------------
  //   
  //                  AnalysisName  RunType         Gain    DBLS (Dynamic BaseLine Substraction)
  //
  //                  AdcPed1       fPEDESTAL_STD      3    No
  //                  AdcPed6       fPEDESTAL_STD      2    No
  //                  AdcPed12      fPEDESTAL_STD      1    No
  //
  //                  AdcPeg12      fPEDESTAL_GAP      1    No
  //
  //                  AdcLaser      fLASER_STD         0    No
  //                  AdcPes12      fPEDSIM            0    No
  //
  //
  //                  AdcSPed1      fPEDESTAL_STD      3    Yes
  //                  AdcSPed6      fPEDESTAL_STD      2    Yes
  //                  AdcSPed12     fPEDESTAL_STD      1    Yes
  //
  //                  AdcSPeg12     fPEDESTAL_GAP      1    Yes
  //
  //                  AdcSLaser     fLASER_STD         0    Yes
  //                  AdcSPes12     fPEDSIM            0    Yes
  //
  //--------------------------------------------------------------------------------------------------

  //................ Run type list
  fLASER_STD    =  4;
  fPEDESTAL_STD =  9;
  fPEDESTAL_GAP = 18;
  fPEDSIM       = 24;
 
  //................ Chozen run type from analysis name
  fChozenRunTypeNumber = fPEDESTAL_STD;   // default
  if( fAnalysisName == "AdcPed1"  || fAnalysisName == "AdcPed6"  || fAnalysisName == "AdcPed12" ||
      fAnalysisName == "AdcSPed1" || fAnalysisName == "AdcSPed6" || fAnalysisName == "AdcSPed12" )
    {fChozenRunTypeNumber = fPEDESTAL_STD;}
  if( fAnalysisName == "AdcPeg12" || fAnalysisName == "AdcSPeg12" )
    {fChozenRunTypeNumber = fPEDESTAL_GAP;}
  if( fAnalysisName == "AdcLaser"  || fAnalysisName == "AdcSLaser"  ){fChozenRunTypeNumber = fLASER_STD;}
  if( fAnalysisName == "AdcPes12 " || fAnalysisName == "AdcSPes12 " ){fChozenRunTypeNumber = fPEDSIM;}

  //................ Gains from analysis name
  fChozenGainNumber = 0;   // default => event always accepted if fChozenGainNumber = 0 ( see USER's Analysis cut in ::analyze(...) )
  if( fAnalysisName == "AdcPed1"   || fAnalysisName == "AdcSPed1"  ){fChozenGainNumber = 3;}
  if( fAnalysisName == "AdcPed6"   || fAnalysisName == "AdcSPed6"  ){fChozenGainNumber = 2;}
  if( fAnalysisName == "AdcPed12"  || fAnalysisName == "AdcSPed12" ){fChozenGainNumber = 1;}
  if( fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcSPeg12" ){fChozenGainNumber = 0;}
  if( fAnalysisName == "AdcLaser"  || fAnalysisName == "AdcSLaser" ){fChozenGainNumber = 0;}
  if( fAnalysisName == "AdcPes12 " || fAnalysisName == "AdcSPes12 "){fChozenGainNumber = 0;}

  //............... Flag for Dynamic BaseLine Substraction from analysis name
  fDynBaseLineSub = "no";   // default
  if( fAnalysisName == "AdcPed1"   || fAnalysisName == "AdcPed6"  || fAnalysisName == "AdcPed12" ||
      fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcLaser" || fAnalysisName == "AdcPes12 " )
    {fDynBaseLineSub = "no";}
  if( fAnalysisName == "AdcSPed1"  || fAnalysisName == "AdcSPed6"  || fAnalysisName == "AdcSPed12" ||
      fAnalysisName == "AdcSPeg12" || fAnalysisName == "AdcSLaser" || fAnalysisName == "AdcSPes12 " )
    {fDynBaseLineSub = "yes";}

  //....................... Index range for ECNA init and for loop on BuildEventDistribution calls
  if( fStexNumber == 0 )
    {
      if( fStexName == "SM"  )
	{
	  fSMIndexBegin   = 0;             fSMIndexStop   = fMyEBEcal->MaxSMInEB();
	  fStexIndexBegin = fSMIndexBegin; fStexIndexStop = fSMIndexStop;
	  fDeeIndexBegin  = 0;             fDeeIndexStop  = 0;
	}
      if( fStexName == "Dee"  )
	{
	  fSMIndexBegin   = 0;              fSMIndexStop   = 0;
	  fDeeIndexBegin  = 0;              fDeeIndexStop  = fMyEEEcal->MaxDeeInEE();
	  fStexIndexBegin = fDeeIndexBegin; fStexIndexStop = fDeeIndexStop;
	}
    }
  else
    {
      if( fStexName == "SM"  )
	{
	  fSMIndexBegin   = fStexNumber-1; fSMIndexStop   = fStexNumber;
	  fStexIndexBegin = fSMIndexBegin; fStexIndexStop = fSMIndexStop;
	  fDeeIndexBegin  = 0;             fDeeIndexStop  = 0;
	}
      if( fStexName == "Dee" )
	{
	  fSMIndexBegin   = 0;              fSMIndexStop   = 0;
	  fDeeIndexBegin  = fStexNumber-1;  fDeeIndexStop  = fStexNumber;
	  fStexIndexBegin = fDeeIndexBegin; fStexIndexStop = fDeeIndexStop;
	}
    }

  //......... DATA DEPENDENT PARAMETERS
  fRunNumber  = 0;

  fMyCnaEBSM  = 0;
  fMyCnaEEDee = 0;

  fRunTypeNumber  = -1;
  fMgpaGainNumber = -1;

  fFedId  = -1;
  fFedTcc = -1;

  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> fAnalysisName        = " << fAnalysisName << endl;
  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> fNbOfSamples         = " << fNbOfSamples << endl;
  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> fFirstReqEvent       = " << fFirstReqEvent << endl;
  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> fLastReqEvent        = " << fLastReqEvent << endl;
  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> fReqNbOfEvts         = " << fReqNbOfEvts << endl;
  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> fStexName            = " << fStexName << endl;
  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> fStexNumber          = " << fStexNumber << endl;
  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> fChozenRunTypeNumber = " << fChozenRunTypeNumber << endl;
  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> fChozenGainNumber    = " << fChozenGainNumber  << endl << endl;

  std::cout << "*EcnaAnalyzer::EcnaAnalyzer-constructor> Init done. " << endl;
}
// end of constructor

EcnaAnalyzer::~EcnaAnalyzer()
{
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

  using namespace std;
  //..................................... format numerical values
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << setprecision(3) << setw(6);
  cout.setf(ios::dec, ios::basefield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.setf(ios::left, ios::adjustfield);
  cout.setf(ios::right, ios::adjustfield);

  std::cout << "EcnaAnalyzer::~EcnaAnalyzer()> destructor is going to be executed." << endl;

  if( fOutcomeError == kTRUE )return;

  //-------------------------------------------------------------------------------

  //....................................................... EB (SM)
  if( fMyCnaEBSM == 0 && fStexName == "SM" )
    {
      std::cout << endl << "!EcnaAnalyzer-destructor> **** ERROR **** fMyCnaEBSM = " << fMyCnaEBSM
		<< ". !===> ECNA HAS NOT BEEN INITIALIZED. Last event run type = " << runtype(fRunTypeNumber)
		<< ", last event fFedId(+601) = " << fFedId+601 << endl 
		<< ", last event Mgpa gain = " << gainvalue(fMgpaGainNumber) << endl << endl;
    }
  else
    {
      for(Int_t iSM = fSMIndexBegin; iSM < fSMIndexStop; iSM++)
	{
	  if( fMyCnaEBSM[iSM] != 0 )
	    {
	      //........................................ register dates 1 and 2
	      fMyCnaEBSM[iSM]->StartStopDate(fDateFirst[iSM], fDateLast[iSM]);
	      fMyCnaEBSM[iSM]->StartStopTime(fTimeFirst[iSM], fTimeLast[iSM]);
	      
	      //........................................ get the sample values in array
	      fMyCnaEBSM[iSM]->GetReadyToCompute();
	      fMyCnaEBSM[iSM]->SampleValues();
	      
	      //........................................ write the sample values in .root file
	      if( fMyCnaEBSM[iSM]->WriteRootFile() == kFALSE )
		{
		  std::cout << "!EcnaAnalyzer-destructor> PROBLEM with write ROOT file for SM" << iSM+1
			    << fTTBELL << endl;
		}
	    }
	  else
	    {
	      std::cout << "*EcnaAnalyzer-destructor> Calculations and writing on file already done for SM "
			<< iSM+1 << endl;
	    }
	}
      delete fMyCnaEBSM;
    }
  //....................................................... EE (Dee)

  if( fMyCnaEEDee == 0 && fStexName == "Dee" )
    {
      std::cout << endl << "!EcnaAnalyzer-destructor> **** ERROR **** fMyCnaEEDee = " << fMyCnaEEDee
		<< ". !===> ECNA HAS NOT BEEN INITIALIZED. Last event run type = " << runtype(fRunTypeNumber)
		<< ", last event fFedId(+601) = " << fFedId+601 << endl 
		<< ", last event Mgpa gain = " << gainvalue(fMgpaGainNumber) << endl << endl;
    }
  else
    {
      for(Int_t iDee = fDeeIndexBegin; iDee < fDeeIndexStop; iDee++)
	{
	  if( fMyCnaEEDee[iDee] != 0 )
	    {
	      //........................................ register dates 1 and 2
	      fMyCnaEEDee[iDee]->StartStopDate(fDateFirst[iDee], fDateLast[iDee]);
	      fMyCnaEEDee[iDee]->StartStopTime(fTimeFirst[iDee], fTimeLast[iDee]);
	      
	      //........................................ get the sample values in array
	      fMyCnaEEDee[iDee]->GetReadyToCompute();
	      fMyCnaEEDee[iDee]->SampleValues();
	      
	      //........................................ write the sample values in .root file
	      if(fMyCnaEEDee[iDee]->WriteRootFile() == kFALSE )
		{
		  std::cout << "!EcnaAnalyzer-destructor> PROBLEM with write ROOT file for Dee" << iDee+1
			    << fTTBELL << endl;
		}
	    }
	  else
	    {
	      std::cout << "*EcnaAnalyzer-destructor> Calculations and writing on file already done for Dee "
			<< iDee+1 << endl;
	    }
	}
      delete fMyCnaEEDee;
    }
  std::cout <<endl;

  //-----------------------------------------------------------------------------------
  std::cout << "*EcnaAnalyzer-destructor> Numbers of events with ERROR(S) returned by BuildEventDistributions(): "
       << fBadBuildCnaRun << endl;

  std::cout << endl<< "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " <<  endl;

  std::cout << "*EcnaAnalyzer-destructor> Run types seen in event headers before selection:" << endl;

  for(Int_t i=0; i<fMaxRunTypeCounter; i++)
    {
      std::cout << " => " << setw(10) << fRunTypeCounter[i]
		<< " event header(s) with run type " << runtype(i) << endl; 
    }

  std::cout << endl<< "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " <<  endl;

  std::cout << "*EcnaAnalyzer-destructor> Mgpa gains seen in event headers before selection:" << endl;

  for(Int_t i=0; i<fMaxMgpaGainCounter; i++)
    {
      std::cout << " => " << setw(10) << fMgpaGainCounter[i]
		<< " event header(s) with gain " << gainvalue(i) << endl; 
    }

  std::cout << endl<< "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " <<  endl;

  std::cout << "*EcnaAnalyzer-destructor> Numbers of selected events for each FED:" << endl;

  for(Int_t i=0; i<fMaxFedIdCounter; i++)
    {
	  std::cout << " => FedId " << i+601 << ": "
		    << setw(10) << fFedIdCounter[i] << " events" << endl;
    }

  std::cout << endl<< "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " <<  endl;

  if( fStexNumber == 0 )
    {
      std::cout << "*EcnaAnalyzer-destructor> fDateFirst = " << fDateFirst[0] << endl
		<< "                          fDateLast  = " << fDateLast[fMaxTreatedStexCounter-1] << endl << endl;
    }
  if( fStexNumber > 0 )
    {
      std::cout << "*EcnaAnalyzer-destructor> fDateFirst = " << fDateFirst[fStexNumber-1] << endl
		<< "                          fDateLast  = " << fDateLast[fStexNumber-1] << endl << endl;
    }

  std::cout << endl<< "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " <<  endl;

  //std::cout << "*EcnaAnalyzer-destructor> fTimeFirst = " << fTimeFirst[0] << endl
  //     << "                          fTimeLast  = " << fTimeLast[fMaxTreatedStexCounter-1] << endl << endl;
  //
  //std::cout << endl<< "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " <<  endl;

  Int_t n0 =0; CheckMsg(n0);

  delete fMyEENumbering;
  delete fMyEEEcal;

  std::cout << "*EcnaAnalyzer-destructor> End of execution." << endl;
}
// end of destructor


//
// member functions
//

// ------------ method called to produce the data  ------------
void EcnaAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  //..................................... format numerical values
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << setprecision(3) << setw(6);
  cout.setf(ios::dec, ios::basefield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.setf(ios::left, ios::adjustfield);
  cout.setf(ios::right, ios::adjustfield);

  using namespace edm;

  fRecNumber++; 

  Int_t iFreq = (fLastReqEvent - fFirstReqEvent + 1)/5;
  if( iFreq <= 0 ){iFreq = 10000;}

  Int_t MaxSMAndDS = fMyEBEcal->MaxSMInEB()+fMyEEEcal->MaxDSInEE();

  //********************************************* EVENT TREATMENT ********************************
  Handle<EcalRawDataCollection> pEventHeader;
  const EcalRawDataCollection* myEventHeader = 0;
  try{
    iEvent.getByLabel(eventHeaderProducer_, eventHeaderCollection_, pEventHeader);
    myEventHeader = pEventHeader.product();
  }catch (std::exception& ex ){
    std::cerr << "Error! can't get the product " << eventHeaderCollection_.c_str() << std::endl;
  }
  //........... Decode myEventHeader infos
  for(EcalRawDataCollection::const_iterator headerItr=myEventHeader->begin();
      headerItr != myEventHeader->end();++headerItr)
    {
      //===> fRunNumber, fRunTypeNumber, fMgpaGainNumber, fFedId, fEvtNumber
      //     will be used in AnalysisOutcome(...) below
      fRunNumber = (Int_t)headerItr->getRunNumber();
      if( fRunNumber <= 0 ){fRunNumber = (Int_t)iEvent.id().run();}
      fRunTypeNumber  = (Int_t)headerItr->getRunType();
      fMgpaGainNumber = (Int_t)headerItr->getMgpaGain();
      fFedId = (Int_t)headerItr->fedId() - 601;   // 1st Fed = 601, FedId = Fed number - 1 
      fEvtNumber = (Int_t)headerItr->getLV1(); 
      if( fEvtNumber <= 0 ){fEvtNumber = (Int_t)iEvent.id().event();}

      if( fEvtNumber != fEvtNumberMemo )
	{
	  fEvtNumberMemo = fEvtNumber;

	  //============================================
	  //  cmsRun INTERRUPTION if analysis complete
	  //  or if fCurrentEventNumber >= LastReqEvent
	  //============================================
	  if( AnalysisOutcome("EVT") == kTRUE ){return;}

	  // no interruption => event has to be analyzed

	  fCurrentEventNumber++;

	  if( fRecNumber ==   1 || fRecNumber ==   50 || fRecNumber == 100 ||
	      fRecNumber == 500 || fRecNumber == 1000 || fRecNumber%iFreq == 0 ){Int_t n1 =1; CheckMsg(n1);}
  
	  if( fCurrentEventNumber < fFirstReqEvent )return; // skip events before fFirstReqEvent
	}

      //.................. Increment Run type and MgpaGain counters
      if( fRunTypeNumber  >= 0 && fRunTypeNumber  < fMaxRunTypeCounter  ){fRunTypeCounter[fRunTypeNumber]++;}
      if( fMgpaGainNumber >= 0 && fMgpaGainNumber < fMaxMgpaGainCounter ){fMgpaGainCounter[fMgpaGainNumber]++;}

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User's analysis cut %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if( !(
	    fRunNumber > 0 &&
	    fRunTypeNumber == fChozenRunTypeNumber &&
	    ( fMgpaGainNumber == fChozenGainNumber || fChozenGainNumber == 0 )
	    )
	  ) return;
      
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if( fMemoCutOK == 0 ){fMemoCutOK = 1;}

      //---- Accelerating selection with "FED-TCC" number [ from headerItr->getDccInTCCCommand() ]
      //     Arrays fSMFromFedTcc[] and fESFromFedTcc[] are initialised in Init()

      if( fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcSPeg12" )
	{
	  fFedTcc = (Int_t)headerItr->getDccInTCCCommand();

	  if( fFedTcc >= 1 && fFedTcc <= MaxSMAndDS )
	    {
	      if( fStexName == "SM" )
		{
		  if( fFedTcc < 10 || fFedTcc > 45 )return;
		  
		  if( fSMFromFedTcc[fFedTcc-1] >= 1 && fSMFromFedTcc[fFedTcc-1] <= fMyEBEcal->MaxSMInEB() &&
		      fStexNbOfTreatedEvents[fSMFromFedTcc[fFedTcc-1]-1] >= fReqNbOfEvts )return;
		}
	      
	      if( fStexName == "Dee" )
		{
		  if( fFedTcc >= 10 && fFedTcc <= 45 )return;
		  
		  if( fESFromFedTcc[fFedTcc-1] >= 1 && fESFromFedTcc[fFedTcc-1] <= fMyEEEcal->MaxDSInEE() &&
		      fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1] >= fReqNbOfEvts )return;
		}
	    } // end of if( fFedTcc >= 1 && fFedTcc <= MaxSMAndDS )
	} // end of if( fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcSPeg12" )
      
      //.................. Increment FedId counters     
      if( fFedId >= 0 && fFedId < fMaxFedIdCounter ){fFedIdCounter[fFedId]++;}
      
    } // end of for(EcalRawDataCollection::const_iterator headerItr=myEventHeader->begin();
      // headerItr != myEventHeader->end();++headerItr)
  
  if( fMemoCutOK == 0 )return;   // return if no event passed the user's analysis cut

  //========================== SELECTED EVENTS ================================

  fNbOfSelectedEvents++;
  if( fNbOfSelectedEvents == 1 ){Int_t n2 = 2; CheckMsg(n2);}

  //============================ Ecna init for the pointers array =================================
  //.................................................................. EB (SM)
  if( fMyCnaEBSM == 0 && fStexName == "SM" )
    {
      fMyCnaEBSM = new TEcnaRun*[fMyEBEcal->MaxSMInEB()];

      for(Int_t i0SM = 0; i0SM < fMyEBEcal->MaxSMInEB(); i0SM++)
	{fMyCnaEBSM[i0SM] = 0;}

      //for(Int_t iSM = 0; iSM < fMyEBEcal->MaxSMInEB(); iSM++)
      //{
      // if( (fStexNumber > 0 && iSM == fStexNumber-1) || (fStexNumber == 0) )
      //  {
      //   fMyCnaEBSM[iSM] = new TEcnaRun("EB", fNbOfSamples);
      //   fMyCnaEBSM[iSM]->GetReadyToReadData(fAnalysisName,  fRunNumber,
      //					  fFirstReqEvent, fLastReqEvent, fReqNbOfEvts,
      //					  iSM+1,          fRunTypeNumber);
      //   std::cout << "*EcnaAnalyzer::analyze(...)> ********* INIT ECNA EB ********* " << endl
      //	   << "                                   fAnalysisName = " << fAnalysisName << endl
      //	   << "                                      fRunNumber = " << fRunNumber << endl
      //	   << "                                  fFirstReqEvent = " << fFirstReqEvent << endl
      //	   << "                                   fLastReqEvent = " << fLastReqEvent << endl
      //	   << "                                    fReqNbOfEvts = " << fReqNbOfEvts << endl
      //	   << "                                              SM = " << iSM+1 << endl
      //           << "                                        run type = " << runtype(fRunTypeNumber) << endl;
      // }
      //  if( fStexNumber > 0 && iSM != fStexNumber-1 ){fMyCnaEBSM[iSM] = new TEcnaRun();}
      //}

    }
  //.................................................................. EE (Dee)
  if( fMyCnaEEDee == 0 && fStexName == "Dee" )
    {
      fMyCnaEEDee = new TEcnaRun*[fMyEEEcal->MaxDeeInEE()];
      for(Int_t iDee = 0; iDee < fMyEEEcal->MaxDeeInEE(); iDee++)
	{fMyCnaEEDee[iDee] = 0;}

      
      //for(Int_t iDee = 0; iDee < fMyEEEcal->MaxDeeInEE(); iDee++)
      //{ 
      // if( (fStexNumber > 0 && iDee == fStexNumber-1) || (fStexNumber == 0) )
      //   {
      //      fMyCnaEEDee[iDee] = new TEcnaRun("EE", fNbOfSamples);
      //      fMyCnaEEDee[iDee]->GetReadyToReadData(fAnalysisName,  fRunNumber,
      // 					    fFirstReqEvent, fLastReqEvent, fReqNbOfEvts,
      // 					    iDee+1,         fRunTypeNumber);
      //      
      //     std::cout << "*EcnaAnalyzer::analyze(...)> ********* INIT ECNA EE ********* " << endl
      // 	   << "                                   fAnalysisName = " << fAnalysisName << endl
      // 	   << "                                      fRunNumber = " << fRunNumber << endl
      // 	   << "                                  fFirstReqEvent = " << fFirstReqEvent << endl
      // 	   << "                                   fLastReqEvent = " << fLastReqEvent << endl
      // 	   << "                                    fReqNbOfEvts = " << fReqNbOfEvts << endl
      // 	   << "                                             Dee = " << iDee+1 << endl
      // 	   << "                                        run type = " << runtype(fRunTypeNumber) << endl;
      //     }
      //   if( fStexNumber > 0 && iDee != fStexNumber-1 ){fMyCnaEEDee[iDee] = new TEcnaRun();}
      // 	}
    }
  
  //============================ EVENT TREATMENT ==============================

  Int_t MaxNbOfStex = 0;
  if( fStexName == "SM"  ){MaxNbOfStex = fMyEBEcal->MaxSMInEB();}
  if( fStexName == "Dee" ){MaxNbOfStex = fMyEEEcal->MaxDeeInEE();}

  if( (fStexNumber > 0 && fNbOfTreatedStexs == 0) || (fStexNumber == 0 && fNbOfTreatedStexs < MaxNbOfStex) )
    {
      //================================================================= Record type EB (SM)
      if( fStexName == "SM" && fSMIndexBegin < fSMIndexStop )
	{
	  //......................................... Get digisEB
	  Handle<EBDigiCollection> pdigisEB;
	  const EBDigiCollection* digisEB = 0;
	  try{
	    iEvent.getByLabel(digiProducer_, EBdigiCollection_, pdigisEB); 
	    digisEB = pdigisEB.product();
	  }catch (std::exception& ex ){
	    std::cerr << "Error! can't get the product " << EBdigiCollection_.c_str() << std::endl;
	  }
	  
	  // Initialize vectors if not already done
	  if ( int(digisEB->size()) > nChannels_ ){nChannels_ = digisEB->size();}
	  
	  //Int_t print_count = 0;
	  if( Int_t(digisEB->end()-digisEB->begin()) >= 0 &&
	      Int_t(digisEB->end()-digisEB->begin()) <=  Int_t(digisEB->size()) )
	    {
	      //..........................................EB
	      //===============================================================================
	      //
	      //                    Loop over Ecal barrel digisEB (Xtals)
	      //
	      //===============================================================================

	      for(EBDigiCollection::const_iterator digiItr = digisEB->begin(); digiItr != digisEB->end(); ++digiItr)
		{
		  EBDetId id_crystal(digiItr->id());
		  // Int_t HashedIndex = id_crystal.hashedIndex();
		  
		  Int_t i0SM = id_crystal.ism() - 1;    //   <============== GET the SM number - 1 here
		  
		  //=============================== Init Ecna EB ===============================
		  if( fMyCnaEBSM[i0SM] == 0 )
		    {
		      fMyCnaEBSM[i0SM] = new TEcnaRun("EB", fNbOfSamples);
		      fMyCnaEBSM[i0SM]->GetReadyToReadData(fAnalysisName,  fRunNumber,
							   fFirstReqEvent, fLastReqEvent, fReqNbOfEvts,
							   i0SM+1,         fRunTypeNumber);
		      
		      std::cout << "*EcnaAnalyzer::analyze(...)> ********* INIT ECNA EB ********* " << endl
				<< "                                   fAnalysisName = " << fAnalysisName << endl
				<< "                                      fRunNumber = " << fRunNumber << endl
				<< "                                  fFirstReqEvent = " << fFirstReqEvent << endl
				<< "                                   fLastReqEvent = " << fLastReqEvent << endl
				<< "                                    fReqNbOfEvts = " << fReqNbOfEvts << endl
				<< "                                              SM = " << i0SM+1 << endl
				<< "                                        run type = " << runtype(fRunTypeNumber) << endl;
		    }
		  //============================================================================

		  if( i0SM >= 0 && i0SM<fMaxTreatedStexCounter )
		    {
		      fStexDigiOK[i0SM]++;
		      if( fStexDigiOK[i0SM] == 1 ){fStexNbOfTreatedEvents[i0SM]++;}

		      if( fStexNbOfTreatedEvents[i0SM] >= 1 && fStexNbOfTreatedEvents[i0SM] <= fReqNbOfEvts )
		      {
			//......................................... date of first event
			if( fStexDigiOK[i0SM] == 1 && fStexNbOfTreatedEvents[i0SM] == 1 &&
			    ( fStexNumber == 0 || i0SM+1 == fStexNumber )
			    )
			  {
			    edm::Timestamp Time = iEvent.time();
			    edm::TimeValue_t t_current_ev_time = (cond::Time_t)Time.value();
			    time_t         i_current_ev_time = (time_t)(t_current_ev_time>>32);
			    const time_t*  p_current_ev_time = &i_current_ev_time;
			    char*          astime            = ctime(p_current_ev_time);
			    fTimeFirst[i0SM] = i_current_ev_time;
			    fDateFirst[i0SM] = astime;
			    std::cout << "*----> beginning of analysis for " << fStexName << i0SM+1
				      << ". First event date : " << astime << endl;
			      //      << " t_current_ev_time = " << t_current_ev_time  << endl
			      //      << " i_current_ev_time = " << i_current_ev_time  << endl
			      //      << " p_current_ev_time = " << p_current_ev_time  << endl
			  }

			  //=============================================> CUT on i0SM value
			  if( (fStexNumber > 0 && i0SM == fStexNumber-1) || (fStexNumber == 0) )
			    { 
			      Int_t iEta = id_crystal.ietaSM(); // ietaSM() : range = [1,85]
			      Int_t iPhi = id_crystal.iphiSM(); // iphiSM() : range = [1,20]
			      
			      Int_t n1SMCrys  =
				(iEta-1)*(fMyEBEcal->MaxTowPhiInSM()*fMyEBEcal->MaxCrysPhiInTow())+iPhi; // range = [1,1700]
			      Int_t n1SMTow   = fMyEBNumbering->Get1SMTowFrom1SMCrys(n1SMCrys);          // range = [1,68]
			      Int_t i0TowEcha = fMyEBNumbering->Get0TowEchaFrom1SMCrys(n1SMCrys);        // range = [0,24]
			      
			      Int_t NbOfSamplesFromDigis = digiItr->size();
			      
			      EBDataFrame df( *digiItr );
			      
			      if( NbOfSamplesFromDigis > 0 && NbOfSamplesFromDigis <= fMyEBEcal->MaxSampADC() )
				{
				  Double_t adcDBLS = (Double_t)0;
				  // Three 1st samples mean value for Dynamic Base Line Substraction (DBLS)
				  if( fDynBaseLineSub == "yes" )
				    {
				      for (Int_t i0Sample=0; i0Sample<3; i0Sample++)  
					{adcDBLS += (Double_t)(df.sample(i0Sample).adc());}
				      adcDBLS /= (Double_t)3;
				    }
				  // Loop over the samples
				  for (Int_t i0Sample=0; i0Sample<fNbOfSamples; i0Sample++)  
				    {
				      Double_t adc = (Double_t)(df.sample(i0Sample).adc()) - adcDBLS;
				      //................................................. Calls to BuildEventDistributions
				      if( fMyCnaEBSM[i0SM]->BuildEventDistributions
					  (fStexNbOfTreatedEvents[i0SM],n1SMTow,i0TowEcha,i0Sample,adc) == kFALSE )
					{fBadBuildCnaRun++;}
				    }
				}
			      else
				{
				  std::cout << "EcnaAnalyzer::analyze(...)> NbOfSamplesFromDigis out of bounds = "
					    << NbOfSamplesFromDigis << endl;
				}
			    } // end of if( (fStexNumber > 0 && i0SM == fStexNumber-1) || (fStexNumber == 0) )
		      } // end of if( fStexNbOfTreatedEvents[i0SM] >= 1 && fStexNbOfTreatedEvents[i0SM] <= fReqNbOfEvts )
		    } // end of if( i0SM >= 0 && i0SM<fMyEBEcal->MaxSMInEB() )
		} // end of for (EBDigiCollection::const_iterator digiItr = digisEB->begin();
	          //             digiItr != digisEB->end(); ++digiItr)
	      
	      for(Int_t i0SM=0; i0SM<fMaxTreatedStexCounter; i0SM++)
		{
		  fStexDigiOK[i0SM] = 0;  // reset fStexDigiOK[i0SM] after loop on digis
		}

	    } // end of if( Int_t(digisEB->end()-digisEB->begin()) >= 0 &&
	      // Int_t(digisEB->end()-digisEB->begin()) <=  Int_t(digisEB->size()) )
	} // end of if( fStexName == "SM" && fSMIndexBegin < fSMIndexStop )
      

      //=============================================================== Record type EE (Dee)
      if( fStexName == "Dee" && fDeeIndexBegin < fDeeIndexStop )
	{
	  //......................................... Get digisEE
	  Handle<EEDigiCollection> pdigisEE;
	  const EEDigiCollection* digisEE = 0;
	  try{
	    iEvent.getByLabel(digiProducer_, EEdigiCollection_, pdigisEE); 
	    digisEE = pdigisEE.product();
	  }catch (std::exception& ex ){
	    std::cerr << "Error! can't get the product " << EEdigiCollection_.c_str() << std::endl;
	  }
      
	  // Initialize vectors if not already done
	  if( int(digisEE->size()) > nChannels_ ){nChannels_ = digisEE->size();}
      
	  //Int_t print_count = 0;
	  if( Int_t(digisEE->end()-digisEE->begin()) >= 0 &&
	      Int_t(digisEE->end()-digisEE->begin()) <= Int_t(digisEE->size()) )
	    {
	      //======================================================================================
	      //
	      //                           Loop over Ecal endcap digisEE (Xtals)
	      //
	      //======================================================================================

	      for (EEDigiCollection::const_iterator digiItr = digisEE->begin(); digiItr != digisEE->end(); ++digiItr) 
		{	  
		  EEDetId id_crystal(digiItr->id());

		  Int_t iX_data = id_crystal.ix();        // iX_data : range = [1,100]
		  Int_t iY_data = id_crystal.iy();        // iY_data : range = [1,100] 
		  Int_t i_quad  = id_crystal.iquadrant(); // iquadrant() : range = [1,4]
		  Int_t i_sgnZ  = id_crystal.zside();     //     zside() : values = -1,+1    
		    
		  Int_t iX = iX_data;
		  Int_t iY = iY_data; // iY : range = [1,100]

		  //.......... See CMSSW/DataFormats/EcalDetId/src/EEDetId.cc::ixQuadrantOne()  [ in which ix() = iX_data ]
		  if( i_quad == 1 || i_quad == 4 ){iX = iX_data - 50;} // iX_data : range = [51,100], iX : range = [1,50]
		  if( i_quad == 3 || i_quad == 2 ){iX = 51 - iX_data;} // iX_data : range = [50,1],   iX : range = [1,50]
		  		    
		  Int_t n1DeeCrys =
		    (iX-1)*(fMyEEEcal->MaxSCIYInDee()*fMyEEEcal->MaxCrysIYInSC()) + iY; // n1DeeCrys: range = [1,5000]

		  Int_t n1DeeNumber = 0;
		  if( i_quad == 1 && i_sgnZ ==  1 ){n1DeeNumber = 2;}
		  if( i_quad == 1 && i_sgnZ == -1 ){n1DeeNumber = 3;}
		  if( i_quad == 2 && i_sgnZ ==  1 ){n1DeeNumber = 1;}
		  if( i_quad == 2 && i_sgnZ == -1 ){n1DeeNumber = 4;}
		  if( i_quad == 3 && i_sgnZ ==  1 ){n1DeeNumber = 1;}
		  if( i_quad == 3 && i_sgnZ == -1 ){n1DeeNumber = 4;}
		  if( i_quad == 4 && i_sgnZ ==  1 ){n1DeeNumber = 2;}
		  if( i_quad == 4 && i_sgnZ == -1 ){n1DeeNumber = 3;}

		  Int_t i0Dee = n1DeeNumber-1;   //   <============== GET the Dee number - 1 here

		  //=============================== Init Ecna EE ===============================
		  if( fMyCnaEEDee[i0Dee] == 0 )
		    {
		      fMyCnaEEDee[i0Dee] = new TEcnaRun("EE", fNbOfSamples);
		      fMyCnaEEDee[i0Dee]->GetReadyToReadData(fAnalysisName,  fRunNumber,
							     fFirstReqEvent, fLastReqEvent, fReqNbOfEvts,
							     i0Dee+1,         fRunTypeNumber);
		      
		      std::cout << "*EcnaAnalyzer::analyze(...)> ********* INIT ECNA EE ********* " << endl
				<< "                                   fAnalysisName = " << fAnalysisName << endl
				<< "                                      fRunNumber = " << fRunNumber << endl
				<< "                                  fFirstReqEvent = " << fFirstReqEvent << endl
				<< "                                   fLastReqEvent = " << fLastReqEvent << endl
				<< "                                    fReqNbOfEvts = " << fReqNbOfEvts << endl
				<< "                                             Dee = " << i0Dee+1 << endl
				<< "                                        run type = " << runtype(fRunTypeNumber) << endl;
		    }
		  //============================================================================

		  if( i0Dee >= 0 && i0Dee<fMaxTreatedStexCounter )
		    {
		      Bool_t cOKForTreatment = kFALSE;

		      if( fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcSPeg12" )
			{
			  if( fFedTcc >= 1 && fFedTcc <= MaxSMAndDS )
			    {
			      fFedDigiOK[fESFromFedTcc[fFedTcc-1]-1]++;

			      if( !(fESFromFedTcc[fFedTcc-1] == 5 || fESFromFedTcc[fFedTcc-1] == 14) )
				{  
				  if( fFedDigiOK[fESFromFedTcc[fFedTcc-1]-1] == 1 )
				    {
				      fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1]++;
				    }
				  if( fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1] >= 1 &&
				      fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1] <= fReqNbOfEvts )
				    {
				      fStexNbOfTreatedEvents[i0Dee] = fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1];
				      cOKForTreatment = kTRUE;
				    }
				}
			      if( fESFromFedTcc[fFedTcc-1] == 5 || fESFromFedTcc[fFedTcc-1] == 14 )
				{
				  if( fFedDigiOK[fESFromFedTcc[fFedTcc-1]-1] == 1 )
				    {
				      fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1]++;
				      fDeeDS5Memo1 = n1DeeNumber;
				      fStexNbOfTreatedEvents[i0Dee] = fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1];
				    }
				  else
				    {
				      if( fDeeDS5Memo2 == 0 )
					{
					  if( n1DeeNumber != fDeeDS5Memo1 )
					    {
					      // change of Dee in Data sector 5
					      fDeeDS5Memo2 = n1DeeNumber;
					      fStexNbOfTreatedEvents[i0Dee] =
						fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1];
					    }
					}
				    }
				  if( fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1] >= 1 &&
				      fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1] <= fReqNbOfEvts )
				    {
				      cOKForTreatment = kTRUE;
				    }
				}
			    } // end of if( fFedTcc >= 1 && fFedTcc <= MaxSMAndDS )
			} // end of if( fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcSPeg12" )
		      else
			{
			  fStexDigiOK[i0Dee]++;
			  if( fStexDigiOK[i0Dee] == 1 ){fStexNbOfTreatedEvents[i0Dee]++;}
			  if( fStexNbOfTreatedEvents[i0Dee] >= 1 &&
			      fStexNbOfTreatedEvents[i0Dee] <= fReqNbOfEvts )
			    {cOKForTreatment = kTRUE;}
			}

		      if( cOKForTreatment== kTRUE )
			{
			  //......................................... date of first event
			  if( ( !(fAnalysisName == "AdcPeg12" || fAnalysisName == "AdcSPeg12") &&
				fStexDigiOK[i0Dee] == 1 &&
				fStexNbOfTreatedEvents[i0Dee] == 1 ) ||
			      (  (fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcSPeg12") && 
				 fFedDigiOK[fESFromFedTcc[fFedTcc-1]-1] == 1 &&
				 fStexNbOfTreatedEvents[i0Dee] == 1 &&
				 fMemoDateFirstEvent[i0Dee] == 0 ) )
			    {
			      edm::Timestamp Time = iEvent.time();
			      edm::TimeValue_t t_current_ev_time = (cond::Time_t)Time.value();
			      time_t         i_current_ev_time = (time_t)(t_current_ev_time>>32);
			      const time_t*  p_current_ev_time = &i_current_ev_time;
			      char*          astime            = ctime(p_current_ev_time);
			      fTimeFirst[i0Dee] = i_current_ev_time;
			      fDateFirst[i0Dee] = astime;
			      std::cout << "----- beginning of analysis for " << fStexName << i0Dee+1 << "-------"  << endl
				//<< " t_current_ev_time = " << t_current_ev_time  << endl
				//<< " i_current_ev_time = " << i_current_ev_time  << endl
				//<< " p_current_ev_time = " << p_current_ev_time  << endl
					<< " First event date  = " << astime << endl
					<< " Nb of selected evts = " << fNbOfSelectedEvents << endl
					<< "---------------------------------------------------------------"  << endl;
			      fMemoDateFirstEvent[i0Dee]++;
			    }

			  //=============================================> cut on i0Dee value
			  if( (fStexNumber > 0 && i0Dee == fStexNumber-1) || (fStexNumber == 0) )
			    {
			      TString sDir = fMyEENumbering->GetDeeDirViewedFromIP(n1DeeNumber);
			      Int_t n1DeeSCEcna = fMyEENumbering->Get1DeeSCEcnaFrom1DeeCrys(n1DeeCrys, sDir);
			      Int_t i0SCEcha = fMyEENumbering->Get1SCEchaFrom1DeeCrys(n1DeeCrys, sDir) - 1;
			      
			      Int_t NbOfSamplesFromDigis = digiItr->size();
			      
			      EEDataFrame df( *digiItr );
			      
			      if( NbOfSamplesFromDigis > 0 && NbOfSamplesFromDigis <= fMyEEEcal->MaxSampADC() )
				{
				  Double_t adcDBLS = (Double_t)0;
				  // Three 1st samples mean value for Dynamic Base Line Substraction (DBLS)
				  if( fDynBaseLineSub == "yes" )
				    {
				      for (Int_t i0Sample=0; i0Sample<3; i0Sample++)  
					{adcDBLS += (Double_t)(df.sample(i0Sample).adc());}
				      adcDBLS /= (Double_t)3;
				    }
				  // Loop over the samples
				  for (Int_t i0Sample = 0; i0Sample < fNbOfSamples; i0Sample++)  
				    {
				      Double_t adc = (Double_t)(df.sample(i0Sample).adc()) - adcDBLS;
				      //................................................. Calls to BuildEventDistributions
				      if( fMyCnaEEDee[i0Dee]->BuildEventDistributions
					  (fStexNbOfTreatedEvents[i0Dee],n1DeeSCEcna,i0SCEcha,i0Sample,adc) == kFALSE )
					{fBadBuildCnaRun++;}
				    }
				}
			      else
				{
				  std::cout << "EcnaAnalyzer::analyze(...)> NbOfSamplesFromDigis out of bounds = "
					    << NbOfSamplesFromDigis << endl;
				}
			    } // end of if( (fStexNumber > 0 && i0Dee == fStexNumber-1) || (fStexNumber == 0) )
			} // end of if( fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1] >= 1 &&
			  // fFedNbOfTreatedEvents[fESFromFedTcc[fFedTcc-1]-1] <= fReqNbOfEvts )
		    } // end of if( i0Dee >= 0 && i0Dee<fMyEEEcal->MaxDeeInEE() )
		} // end of for (EBDigiCollection::const_iterator digiItr = digisEB->begin();
                  //             digiItr != digisEB->end(); ++digiItr)

	      // reset fStexDigiOK[i0Dee] or fFedDigiOK[i0Dee] to zero after loop on digis
	      if( fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcSPeg12" )
		{
		  for(Int_t i0FedES=0; i0FedES<fMaxFedUnitCounter; i0FedES++)
		    {fFedDigiOK[i0FedES] = 0;}
		  
		  // reset fDeeDS5Memo1 and fDeeDS5Memo2 (for Data sector 5 versus Dee number  management)
		  fDeeDS5Memo1 = 0;
		  fDeeDS5Memo2 = 0;
		}
	      else
		{
		  for(Int_t i0Dee=0; i0Dee<fMaxTreatedStexCounter; i0Dee++)
		    {fStexDigiOK[i0Dee] = 0;}
		}
	      
	    } // end of if( Int_t(digisEB->end()-digisEB->begin()) >= 0 &&
	      // Int_t(digisEB->end()-digisEB->begin()) <=  Int_t(digisEB->size()) )

	} // end of if( fStexName == "Dee" && fDeeIndexBegin < fDeeIndexStop )
    }  // end of if( (fStexNumber > 0 && fNbOfTreatedStexs == 0) || (fStexNumber == 0 && fNbOfTreatedStexs < MaxNbOfStex) )

  //=============================================================================================
  //
  //                    Number of treated events. Setting Stex and Fed status.
  //
  //=============================================================================================

  // (take into account the "Accelerating selection with FED number" section - see above -)
  if( fStexName == "SM" ||
      ( fStexName == "Dee" && !( fAnalysisName == "AdcPeg12" || fAnalysisName == "AdcSPeg12" ) )
      )  // one FED = one SM = one Stex
    {
      for(Int_t i0Stex=fStexIndexBegin; i0Stex<fStexIndexStop; i0Stex++)
	{
	  if( fStexStatus[i0Stex] != 2 )  // do not change fStexStatus[i0Stex] if already set to 2
	    // even if fStexNbOfTreatedEvents[i0Stex] == fReqNbOfEvts 
	    {
	      if( fStexNbOfTreatedEvents[i0Stex] == fReqNbOfEvts ){fStexStatus[i0Stex] = 1;}
	      if( fStexNbOfTreatedEvents[i0Stex] >  fReqNbOfEvts ){fStexStatus[i0Stex] = 2;}
	    }
	}
    }

  // one FED = one Data Sector (DS or ES)
  if( fStexName == "Dee" &&
      ( fAnalysisName == "AdcPeg12" || fAnalysisName == "AdcSPeg12" ) )
    {
      for(Int_t i0FedES=0; i0FedES<fMaxFedUnitCounter; i0FedES++)
	{
	  if( fFedStatus[i0FedES] != 2 )  // do not change fFedStatus[i0FedES] if already set to 2
	    // even if fFedNbOfTreatedEvents[i0FedES] == fReqNbOfEvts 
	    {
	      if( fFedNbOfTreatedEvents[i0FedES] == fReqNbOfEvts )
		{
		  fFedStatus[i0FedES] = 1;
		  fTreatedFedOrder++; fFedStatusOrder[i0FedES] = fTreatedFedOrder;
		}
	      if( fFedNbOfTreatedEvents[i0FedES] >  fReqNbOfEvts ){fFedStatus[i0FedES] = 2;}
	    }
	}
    
      Int_t j0Fed = 4;
      //..................................................... Dee 4 (DS 1,2,3,4 ; ES 1,2,3,4)
      for( Int_t i0FedES= 0; i0FedES<= 3; i0FedES++ )
	{if( fFedStatus[i0FedES] == 1 ){fNbOfTreatedFedsInDee[3]++; fFedStatus[i0FedES] = 2;}}
    
      //..................................................... Dee 3, Dee 4 (DS 5 ; ES 5)
      j0Fed = 4;
      if( fFedStatus[j0Fed] == 1 )
	{fNbOfTreatedFedsInDee[3]++; fNbOfTreatedFedsInDee[2]++; fFedStatus[j0Fed] = 2;}
    
       //.................................................... Dee 3 (DS 6,7,8,9 ; ES 6,7,8,9)
      for( Int_t i0FedES= 5; i0FedES<= 8; i0FedES++ )
	{if( fFedStatus[i0FedES] == 1 ){fNbOfTreatedFedsInDee[2]++; fFedStatus[i0FedES] = 2;}}
    
      //..................................................... Dee 1 (DS 1,2,3,4 ; ES 10,11,12,13)
      for( Int_t i0FedES= 9; i0FedES<=12; i0FedES++ )
	{if( fFedStatus[i0FedES] == 1 ){fNbOfTreatedFedsInDee[0]++; fFedStatus[i0FedES] = 2;}}
    
      //..................................................... Dee 1, Dee 2 (DS 5 ; ES 5)
      j0Fed = 13;
      if( fFedStatus[j0Fed] == 1 )
	{fNbOfTreatedFedsInDee[0]++; fNbOfTreatedFedsInDee[1]++; fFedStatus[j0Fed] = 2;}
    
      //..................................................... Dee 2 (DS 6,7,8,9 ; ES 15,16,17,18)
      for( Int_t i0FedES=14; i0FedES<=17; i0FedES++ )
	{if( fFedStatus[i0FedES] == 1 ){fNbOfTreatedFedsInDee[1]++; fFedStatus[i0FedES] = 2;}}
    
      //-----------------------------------------------------
      for( Int_t i0Dee=0; i0Dee<4; i0Dee++ )
	{
	  if( fNbOfTreatedFedsInStex[i0Dee] >= 0 && fNbOfTreatedFedsInStex[i0Dee] < 5 )
	    {fNbOfTreatedFedsInStex[i0Dee] = fNbOfTreatedFedsInDee[i0Dee];}
	  if( fNbOfTreatedFedsInDee[i0Dee] == 5 )
	    {fStexStatus[i0Dee] = 1; fNbOfTreatedFedsInDee[i0Dee] = 0;}
	}
    
    } // end of if( fStexName == "Dee" &&
      // ( fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcSPeg12" ) )

  //----------------------------------------------------------------------------------------------
  for(Int_t i0Stex=fStexIndexBegin; i0Stex<fStexIndexStop; i0Stex++)
    {
      if( fStexStatus[i0Stex] == 1 )
	{
	  fNbOfTreatedStexs++; // increase nb of treated Stex's only if fStexStatus[i0Stex] == 1
	  //....................................................... date of last event
	  edm::Timestamp Time = iEvent.time();
	  edm::TimeValue_t t_current_ev_time = (cond::Time_t)Time.value();
	  time_t         i_current_ev_time = (time_t)(t_current_ev_time>>32);
	  const time_t*  p_current_ev_time = &i_current_ev_time;
	  char*          astime            = ctime(p_current_ev_time);
	  fTimeLast[i0Stex] = i_current_ev_time;
	  fDateLast[i0Stex] = astime;
	  
	  std::cout << "---------- End of analysis for " << fStexName << i0Stex+1 << " -----------" << endl;
	  Int_t n3 = 3; CheckMsg(n3, i0Stex);
	  // std::cout 	   << " t_current_ev_time = " << t_current_ev_time  << endl
	  //<< " i_current_ev_time = " << i_current_ev_time  << endl
	  //<< " p_current_ev_time = " << p_current_ev_time  << endl
	  std::cout 	    << " Last event date           = " << astime << endl
			    << " Number of selected events = " << fNbOfSelectedEvents << endl;
	  std::cout << endl << fNbOfTreatedStexs << " " << fStexName
		    << "'s with " << fReqNbOfEvts << " events analyzed." << endl
		    << "---------------------------------------------------------"  << endl;
	  
	  //================================= WRITE RESULTS FILE
	  if( fStexName == "SM" )
	    {
	      if( fMyCnaEBSM[i0Stex] != 0 )
		{
		  //........................................ register dates 1 and 2
		  fMyCnaEBSM[i0Stex]->StartStopDate(fDateFirst[i0Stex], fDateLast[i0Stex]);
		  fMyCnaEBSM[i0Stex]->StartStopTime(fTimeFirst[i0Stex], fTimeLast[i0Stex]);
		  
		  //........................................ get the sample values in array
		  fMyCnaEBSM[i0Stex]->GetReadyToCompute();
		  fMyCnaEBSM[i0Stex]->SampleValues();
		  
		  //........................................ write the sample values in .root file
		  if( fMyCnaEBSM[i0Stex]->WriteRootFile() == kFALSE )
		    {
		      std::cout << "!EcnaAnalyzer::analyze> PROBLEM with write ROOT file for SM" << i0Stex+1
				<< fTTBELL << endl;
		    }
		}
	      // set pointer to zero in order to avoid recalculation and rewriting at the destructor level
	      delete fMyCnaEBSM[i0Stex];  fMyCnaEBSM[i0Stex] = 0;
	      std::cout << "!EcnaAnalyzer::analyze> Set memory free: delete done for SM " << i0Stex+1 << endl;
	    }

	  if( fStexName == "Dee" )
	    {
	      if( fMyCnaEEDee[i0Stex] != 0 )
		{
		  //........................................ register dates 1 and 2
		  fMyCnaEEDee[i0Stex]->StartStopDate(fDateFirst[i0Stex], fDateLast[i0Stex]);
		  fMyCnaEEDee[i0Stex]->StartStopTime(fTimeFirst[i0Stex], fTimeLast[i0Stex]);
		  
		  //........................................ get the sample values in array
		  fMyCnaEEDee[i0Stex]->GetReadyToCompute();
		  fMyCnaEEDee[i0Stex]->SampleValues();

		  //........................................ write the sample values in .root file
		  if(fMyCnaEEDee[i0Stex]->WriteRootFile() == kFALSE )
		    {
		      std::cout << "!EcnaAnalyzer::analyze> PROBLEM with write ROOT file for Dee" << i0Stex+1
				<< fTTBELL << endl;
		    }
		}
	      // set pointer to zero in order to avoid recalculation and rewriting at the destructor level
	      delete fMyCnaEEDee[i0Stex]; fMyCnaEEDee[i0Stex] = 0;
	      std::cout << "!EcnaAnalyzer::analyze> Set memory free: delete done for Dee " << i0Stex+1 << endl;
	    }

	  fStexStatus[i0Stex] = 2;        // set fStexStatus[i0Stex] to 2 definitively
	  std::cout << "*---------------------------------------------------------------------------- " << endl;

	} // end of if( fStexStatus[i0Stex] == 1 )
    } // end of for(Int_t i0Stex=fStexIndexBegin; i0Stex<fStexIndexStop; i0Stex++)
}
// end of EcnaAnalyzer::analyse(...)

Bool_t EcnaAnalyzer::AnalysisOutcome(const TString s_opt)
{
  //---- STOP if end of analysis

  Bool_t result = kFALSE;

  if( s_opt == "EVT" )
    {
      Int_t MaxNbOfStex = 0;
      if( fStexName == "SM"  ){MaxNbOfStex = fMyEBEcal->MaxSMInEB();}
      if( fStexName == "Dee" ){MaxNbOfStex = fMyEEEcal->MaxDeeInEE();}

      if( ( (fStexNumber > 0 && fNbOfTreatedStexs == 1) ||
	    (fStexNumber == 0 && fNbOfTreatedStexs == MaxNbOfStex) ) &&
	  ( (fLastReqEvent < fFirstReqEvent) ||
	    (fLastReqEvent >= fFirstReqEvent && fCurrentEventNumber <= fLastReqEvent) )
	  )
	{
	  std::cout << endl
		    << "**************************** ANALYSIS REPORT > OK **************************************" << endl
		    << "*EcnaAnalyzer::AnalysisOutcome(...)> The maximum requested number of events and the maximum" << endl
		    << "                                     number of treated " << fStexName
		    << "'s have been reached." << endl
		    << "                                     Analysis successfully ended from EcnaAnalyzer " << endl
		    << "                                     by SIGNAL: kill(getpid(),SIGUSR2)." << endl
		    << "                                     Number of selected events   = " << fNbOfSelectedEvents << endl
		    << "                                     Last requested event number = " << fLastReqEvent << endl
		    << "                                     Current event number        = " << fCurrentEventNumber << endl;

	  Int_t n0 = 0; CheckMsg(n0);

	  std::cout << "****************************************************************************************" << endl
		    << endl;

	  result = kTRUE;
	  kill(getpid(),SIGUSR2);
	}

      if( fLastReqEvent >= fFirstReqEvent && fCurrentEventNumber > fLastReqEvent &&
	  ! ( (fStexNumber > 0 && fNbOfTreatedStexs == 1) ||
	      (fStexNumber == 0 && fNbOfTreatedStexs == MaxNbOfStex) ) )
	{
	  std::cout << endl
		    << "**************************** ANALYSIS REPORT >>> *** WARNING *** WARNING *** WARNING ***" << endl
		    << "*EcnaAnalyzer::AnalysisOutcome(...)> Last event reached before completion of analysis." << endl
		    << "                                     Analysis ended from EcnaAnalyzer " << endl
		    << "                                     by SIGNAL: kill(getpid(),SIGUSR2)." << endl
		    << "                                     Number of selected events   = " << fNbOfSelectedEvents << endl
		    << "                                     Last requested event number = " << fLastReqEvent << endl
		    << "                                     Current event number        = " << fCurrentEventNumber << endl;

	  Int_t n0 = 0; CheckMsg(n0);

	  std::cout << "****************************************************************************************" << endl
		    << endl;
      
	  result = kTRUE;
	  kill(getpid(),SIGUSR2);
	}
    }
  else
    {
      if( s_opt == "ERR_FNEG" )
	{
	  std::cout << endl
		    << "**************************** ANALYSIS REPORT >>> **** ERROR **** ERROR **** ERROR ******" << endl
		    << "*EcnaAnalyzer::AnalysisOutcome(...)> First event number = " << fFirstReqEvent
		    << ". Should be strictly potitive." << endl
		    << "                             Analysis ended from EcnaAnalyzer " << endl
		    << "                             by SIGNAL: kill(getpid(),SIGUSR2)." << endl;

	  std::cout << "****************************************************************************************" << endl
		    << endl;

	  result = kTRUE;
	  kill(getpid(),SIGUSR2);
	}
      if( s_opt == "ERR_LREQ" )
	{
	  std::cout << endl
		    << "**************************** ANALYSIS REPORT >>> **** ERROR **** ERROR **** ERROR ******" << endl
		    << "*EcnaAnalyzer::analyze(...)> Requested number of events = " << fReqNbOfEvts << "." << endl
		    << "                             Too large compared to the event range: "
		    << fFirstReqEvent << " - " << fLastReqEvent << endl
		    << "                             Analysis ended from EcnaAnalyzer " << endl
		    << "                             by SIGNAL: kill(getpid(),SIGUSR2)." << endl;

	  std::cout << "****************************************************************************************" << endl
		    << endl;

	  result = kTRUE;
	  kill(getpid(),SIGUSR2);
	}
    }
  return result;
} // end of EcnaAnalyzer::AnalysisOutcome(const Int_t& n_option)

void EcnaAnalyzer::CheckMsg(const Int_t& MsgNum)
{Int_t nm1 = -1; CheckMsg(MsgNum, nm1);}

void EcnaAnalyzer::CheckMsg(const Int_t& MsgNum, const Int_t& i0Stex)
{
  //------ Cross-check messages

  if( MsgNum == 1 )
    {std::cout << "---------------- CROSS-CHECK A ------------------ " << endl
	       << "**************** CURRENT EVENT ****************** " << endl;}
  if( MsgNum == 2 )
    {std::cout << "---------------- CROSS-CHECK B ------------------ " << endl
	       << "**** FIRST EVENT PASSING USER'S ANALYSIS CUT **** " << endl;}
  if( MsgNum == 3 )
    {std::cout << "---------------- CROSS-CHECK C ------------------ " << endl
	       << "*** CURRENT VALUES BEFORE RESULT FILE WRITING *** " << endl;}
  if( MsgNum == 3 || MsgNum == 4 )
    {std::cout << "          fRecNumber = " << fRecNumber << endl
	       << "          fEvtNumber = " << fEvtNumber << endl;}
  
  std::cout << " fCurrentEventNumber = " << fCurrentEventNumber << endl
	    << " fNbOfSelectedEvents = " << fNbOfSelectedEvents << endl
	    << "          fRunNumber = " << fRunNumber << endl
	    << "     Chozen run type = " << runtype(fChozenRunTypeNumber) << endl
	    << "            Run type = " << runtype(fRunTypeNumber) << endl
	    << "             fFedTcc = " << fFedTcc << endl
	    << "        fFedId(+601) = " << fFedId+601 << endl
	    << "           fStexName = " << fStexName << endl
	    << "         Chozen gain = " << gainvalue(fChozenGainNumber) << endl 
	    << "           Mgpa Gain = " << gainvalue(fMgpaGainNumber) << endl << endl;
	  
  if( fAnalysisName == "AdcPeg12"  || fAnalysisName == "AdcSPeg12" )
    {
      if( fStexName == "SM" )
	{
	  for(Int_t j0Stex=fStexIndexBegin; j0Stex<fStexIndexStop; j0Stex++)
	    {
	      Int_t nStexNbOfTreatedEvents = fStexNbOfTreatedEvents[j0Stex];
	      if( fStexStatus[j0Stex] == 1 ){nStexNbOfTreatedEvents = fStexNbOfTreatedEvents[j0Stex];}
	      if( fStexStatus[j0Stex] == 2 ){nStexNbOfTreatedEvents = fStexNbOfTreatedEvents[j0Stex];}
	      
	      std::cout << fStexName << setw(3) << j0Stex+1 << ": "
			<< setw(5) << nStexNbOfTreatedEvents << " events. "
			<< fStexName << " status: " << fStexStatus[j0Stex];
	      if( j0Stex == i0Stex ){std::cout << " (going to write file for this " << fStexName << ").";}
	      std::cout << endl; 
	    }
	}

      if( fStexName == "Dee" )
	{
	  for(Int_t i0FedES=0; i0FedES<fMaxFedUnitCounter; i0FedES++)
	    {
	      Int_t nFedNbOfTreatedEvents = fFedNbOfTreatedEvents[i0FedES];
	      if( fFedStatus[i0FedES] == 1 ){nFedNbOfTreatedEvents = fFedNbOfTreatedEvents[i0FedES];}
	      if( fFedStatus[i0FedES] == 2 ){nFedNbOfTreatedEvents = fFedNbOfTreatedEvents[i0FedES];}
	      
	      std::cout << "Fed (ES) " << setw(3) << i0FedES+1 << ": "
			<< setw(5) << nFedNbOfTreatedEvents << " events."
			<< " Fed status: " << fFedStatus[i0FedES]
			<< ", order: " << setw(3) << fFedStatusOrder[i0FedES]
			<< " (" << fDeeNumberString[i0FedES] << ")" << endl;
	    }
	  
	  for(Int_t j0Stex=fStexIndexBegin; j0Stex<fStexIndexStop; j0Stex++)
	    {
	      Int_t nStexNbOfTreatedEvents = fStexNbOfTreatedEvents[j0Stex];
	      if( fStexStatus[j0Stex] == 1 ){nStexNbOfTreatedEvents = fStexNbOfTreatedEvents[j0Stex];}
	      if( fStexStatus[j0Stex] == 2 ){nStexNbOfTreatedEvents = fStexNbOfTreatedEvents[j0Stex];}
	      
	      std::cout << fStexName << setw(3) << j0Stex+1 << ": "
			<< setw(5) << fNbOfTreatedFedsInStex[j0Stex] << " analyzed Fed(s). "
			<< fStexName << " status: " << fStexStatus[j0Stex];
	      if( j0Stex == i0Stex ){std::cout << " (going to write file for this " << fStexName << ").";}
	      std::cout << endl; 
	    }
	}

      std::cout << "Number of " << fStexName << "'s with "
		<< fReqNbOfEvts << " events analyzed: " << fNbOfTreatedStexs << endl;
    }
  
  if( MsgNum == 1 || MsgNum == 2 )
    {std::cout << "*---------------------------------------------------------------------------- " << endl;}
  if( MsgNum == 3 )
    {std::cout << "*............................................................................ " << endl;}

} // end of EcnaAnalyzer::CheckMsg(const Int_t& MsgNum, const Int_t& i0Stex)

TString EcnaAnalyzer::runtype(const Int_t& numtype)
{
  TString cType = "?";

  if( numtype ==  0 ){cType = "COSMICS";}
  if( numtype ==  1 ){cType = "BEAMH4";}
  if( numtype ==  2 ){cType = "BEAMH2";}
  if( numtype ==  3 ){cType = "MTCC";}
  if( numtype ==  4 ){cType = "LASER_STD";}
  if( numtype ==  5 ){cType = "LASER_POWER_SCAN";}
  if( numtype ==  6 ){cType = "LASER_DELAY_SCAN";}
  if( numtype ==  7 ){cType = "TESTPULSE_SCAN_MEM";}
  if( numtype ==  8 ){cType = "TESTPULSE_MGPA";}
  if( numtype ==  9 ){cType = "PEDESTAL_STD";}
  if( numtype == 10 ){cType = "PEDESTAL_OFFSET_SCAN";}
  if( numtype == 11 ){cType = "PEDESTAL_25NS_SCAN";}
  if( numtype == 12 ){cType = "LED_STD";}

  if( numtype == 13 ){cType = "PHYSICS_GLOBAL";}
  if( numtype == 14 ){cType = "COSMICS_GLOBAL";}
  if( numtype == 15 ){cType = "HALO_GLOBAL";}

  if( numtype == 16 ){cType = "LASER_GAP";}
  if( numtype == 17 ){cType = "TESTPULSE_GAP";}
  if( numtype == 18 ){cType = "PEDESTAL_GAP";}
  if( numtype == 19 ){cType = "LED_GAP";}

  if( numtype == 20 ){cType = "PHYSICS_LOCAL";}
  if( numtype == 21 ){cType = "COSMICS_LOCAL";}
  if( numtype == 22 ){cType = "HALO_LOCAL";}
  if( numtype == 23 ){cType = "CALIB_LOCAL";}

  if( numtype == 24 ){cType = "PEDSIM";}

  return cType;
}

Int_t EcnaAnalyzer::gainvalue(const Int_t& numgain)
{
  Int_t value = 0;

  if( numgain == 1 ){value = 12;}
  if( numgain == 2 ){value =  6;}
  if( numgain == 3 ){value =  1;}

  return value;
}


