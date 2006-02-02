/** \file
 *
 *  Implementation of CSCQualityTester
 *
 *  $Date: 2006/01/25 16:28:39 $
 *  $Revision: 1.1 $
 *  \author Ilaria Segoni
 */
#include "DQM/CSCMonitorClient/interface/CSCQualityTester.h"
#include "DQM/CSCMonitorClient/interface/CSCQualityTestTypes.h"
#include "DQMServices/QualityTests/interface/QCriterionRoot.h"

#include<iostream>


CSCQualityTester::CSCQualityTester()
{
  
  printout=true;
  logFile.open("CSCDQMClient.log");
  
  if(printout) logFile <<"A CSC Quality Tester Is being Created"<<std::endl;

}


CSCQualityTester::~CSCQualityTester(){ 
  logFile.close();
}

//
void CSCQualityTester::SetupTests(MonitorUserInterface * mui){
    
  if(printout) logFile<<"In CSCQualityTester::SetupTests(...)"<<std::endl;  
  this->SetupTestsFromTextFile(mui);
  this->LinkTeststoME();
  this->AttachRunTests(mui);

}



//
void CSCQualityTester::SetupTestsFromTextFile(MonitorUserInterface * mui){ 


  if(printout) logFile<<"In CSCQualityTester::SetupTestsFromTextFiles(...)"<<std::endl;
  FILE* TestsFile;
  char  fileName[128];
  sprintf(fileName,"QualityTests.db");
  TestsFile= fopen(fileName,"r");

  char  TestType[20];
  char  TestName[20];
  float WarningLevel=0;
  float params[5];
  
  for(int ii=0; ii<5; ii++) params[ii]=0.0;

  int file_status=1;
 
  do{    
     file_status=fscanf(TestsFile,"%s%s%f%f%f%f%f%f\n",
                         &TestType,&TestName,&WarningLevel,
			 &params[0],&params[1],&params[2],&params[3],&params[4]);
			 
    if(file_status != EOF){			 
    
           if(printout)  logFile<<" Reading Quality Tests Configuration: "
                          <<TestType<<" "<<TestName <<" "<<WarningLevel<<" "
	                  <<params[0]<<" "<<params[1]<<" "<<params[2]
	                  <<" "<<params[3]<<" "<<params[4]<< std::endl;
	
	    
           if(!std::strcmp(TestType,dqm::qTestType_csc::XRangeContent)) this->SetContentsXRangeROOTTest(mui,TestName,WarningLevel,params);    
    }

  }while( file_status!=EOF ); 



}


void CSCQualityTester::LinkTeststoME(){ 


  if(printout) logFile<<"in CSCQualityTester::LinkTeststoME"<<std::endl;
  FILE* TestsMEsFile;
  char  fileName[128];
  sprintf(fileName,"TestsMEs.db");
  TestsMEsFile= fopen(fileName,"r");
  
  
  char TestName[20];
  char MEName[20];
  
  int file_status=1;
  
  do{    
      file_status=fscanf(TestsMEsFile,"%s%s\n",&TestName,&MEName);
			 
      if(file_status != EOF){			 
            if(printout) logFile<<"Reading Quality Tests List: " <<TestName <<" "<<MEName << std::endl;
	    
	    std::vector<std::string>  MEList;
	    if( qTestToMEMap.find(TestName) == qTestToMEMap.end()){
 	      MEList.push_back(MEName);
	      qTestToMEMap[TestName]=MEList;
	    }else{
	      MEList=qTestToMEMap[TestName];
	      MEList.push_back(MEName);
	    }	    
      
     }	

  }while( file_status!=EOF ); 



}

void CSCQualityTester::SetContentsXRangeROOTTest(MonitorUserInterface * mui, char  TestName[20], float  WarningLevel, float  params[5]){

  if(printout) logFile<<"In CSCQualityTester::SetContentsXRangeROOTTest, configuring "<<TestName<<" test."<<std::endl;
	
  qTests.push_back(TestName);
  QCriterion * qc1 = mui->createQTest(ContentsXRangeROOT::getAlgoName(),TestName);
  MEContentsXRangeROOT * me_qc1 = (MEContentsXRangeROOT *) qc1;
  me_qc1->setAllowedXRange(params[0],params[1]);
  me_qc1->setWarningProb(WarningLevel);
 
}



void CSCQualityTester::AttachRunTests(MonitorUserInterface * mui){

if(printout) logFile<<"In CSCQualityTester::AttachRunTests(...)"<<std::endl;

for (std::map<std::string, std::vector<std::string> >::iterator testsMap=qTestToMEMap.begin(); 
			   testsMap!=qTestToMEMap.end();++testsMap)	   
    {
     
       std::string testName=testsMap->first;
       std::vector<std::string> MElist=testsMap->second;
    
      	     for(std::vector<std::string>::iterator list = MElist.begin(); list != MElist.end(); ++list){
      	    	  std::string meName = *(list);
      	    	  if(printout) logFile<<"Attaching Test "<< testName <<" to ME "<<meName<<std::endl;
      	    	  mui->useQTest(meName, testName);
      	     }
    }
}



void CSCQualityTester::CheckTests(MonitorUserInterface * mui) 
{
  if(printout) logFile << "In CSCQualityTester::CheckTests(...)" << std::endl;

  
	int status = 1980;
	
	status= mui->getSystemStatus();
        if(printout) logFile << "Possible states: successful "<< dqm::qstatus::STATUS_OK<<", error:  " 
	<< dqm::qstatus::ERROR<<
	", warning:  "<< dqm::qstatus::WARNING<<
	", Other: "<< dqm::qstatus::OTHER<<std::endl;
        if(printout) logFile << "The STATUS IS " << status<<std::endl;
	
	switch(status)
	  {
	  case dqm::qstatus::ERROR:
	    logFile << " Error(s)";
	    break;
	  case dqm::qstatus::WARNING:
	    logFile << " Warning(s)";
	    break;
	  case dqm::qstatus::OTHER:
	    logFile << " Some tests did not run;";
	    break; 
	  default:
	    logFile << " No problems";
	  }

	//MonitorElement * me4 = mui->get("Collector/FU0/C1/C2/histo4");
	//if(me4)
	  //checkTests(me4);
	
  return;
}

  


