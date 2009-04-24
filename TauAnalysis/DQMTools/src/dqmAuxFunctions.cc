#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include <TH1.h>

#include <iostream>

//
//-----------------------------------------------------------------------------------------------------------------------
//

std::string dqmDirectoryName(const std::string& directory)
{
  std::string dirName = directory;
  //if ( dirName == "" || dirName.find_last_of(dqmSeparator) != (dirName.length() - 1) )  dirName.append(dqmSeparator);
//--- add tailing '/'
  if ( dirName != "" && dirName.find_last_of(dqmSeparator) != (dirName.length() - 1) )  dirName.append(dqmSeparator);
  return dirName;
}

std::string dqmSubDirectoryName_merged(const std::string& directory, const std::string& subdirectory)
{ 
  std::string subDirName = subdirectory;
//--- remove characters specifying directory part from name of subdirectory
  if ( subDirName.find(directory) <= 1 ) subDirName.replace(subDirName.find(directory), directory.length(), "");
//--- remove tailing '/'s
  while ( subDirName.find(dqmSeparator) == 0 ) subDirName.replace(subDirName.find(dqmSeparator), dqmSeparator.length(), "");
  return subDirName;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void dqmCheckExistence(DQMStore& dqmStore, const std::string& directoryName, const std::string& meName, int mode, int& errorFlag)
{
  MonitorElement* meOutput = dqmStore.get(dqmDirectoryName(directoryName).append(meName));
  //std::cout << " meOutput = " << meOutput << std::endl;
 
  if ( meOutput ) {
    switch ( mode ) {
    case 1: // print error message 
      edm::LogError ("copyRecursively") << " meName = " << meName << " already exists in directoryName = " << directoryName 
					<< " --> skipping !!";
      errorFlag = 1;
      break;
    case 2:
      break;
    default:
      edm::LogError ("dqmCheckExistence") << " Undefined mode = " << mode << " !!";
      errorFlag = 1;
    }
  }
}

void dqmRegisterHistogram(DQMStore& dqmStore, TH1* histogram, const std::string& name)
{
  //std::cout << "<dqmRegisterHistogram>:" << std::endl;
  //histogram->SetName(std::string(histogram->GetName()).append("_copied").data());
  histogram->SetName(histogram->GetName());  
  if ( TH1F* h = dynamic_cast<TH1F*>(histogram) ) {
    //std::cout << " --> calling DQMStore::book1D" << std::endl;
    dqmStore.book1D(name, h);
  } else if ( TH1S* h = dynamic_cast<TH1S*>(histogram) ) {
    //std::cout << " --> calling DQMStore::book1@" << std::endl;
    dqmStore.book1S(name, h);
  } else if ( TH2F* h = dynamic_cast<TH2F*>(histogram) ) {
    //std::cout << " --> calling DQMStore::book2D" << std::endl;
    dqmStore.book2D(name, h);
  } else if ( TH2S* h = dynamic_cast<TH2S*>(histogram) ) {
    //std::cout << " --> calling DQMStore::book2S" << std::endl;
    dqmStore.book2S(name, h);
  } else if ( TH3F* h = dynamic_cast<TH3F*>(histogram) ) {
    //std::cout << " --> calling DQMStore::book3D" << std::endl;
    dqmStore.book3D(name, h);
  } else if ( TProfile* h = dynamic_cast<TProfile*>(histogram) ) {
    //std::cout << " --> calling DQMStore::bookProfile" << std::endl;
    dqmStore.bookProfile(name, h);
  } else if ( TProfile2D* h = dynamic_cast<TProfile2D*>(histogram) ) {
    //std::cout << " --> calling DQMStore::bookProfile2D" << std::endl;
    dqmStore.bookProfile2D(name, h);
  }
}

void dqmCopyRecursively(DQMStore& dqmStore, const std::string& inputDirectory, const std::string& outputDirectory, 
			double scaleFactor, int mode, bool rmInputDirectory)
{
  //std::cout << "<copyRecursively>:" << std::endl;
  //std::cout << " inputDirectory = " << inputDirectory << std::endl;
  //std::cout << " outputDirectory = " << outputDirectory << std::endl;
  //std::cout << " scaleFactor = " << scaleFactor << std::endl;
  //std::cout << " rmInputDirectory = " << rmInputDirectory << std::endl;

  bool meInput_copied = false;

//--- copy all monitor elements in current inputDirectory to the outputDirectory
  dqmStore.setCurrentFolder(inputDirectory);
  std::vector<std::string> meNames = dqmStore.getMEs();
  for ( std::vector<std::string>::const_iterator meName = meNames.begin();
	meName != meNames.end(); ++meName ) {
    std::string meName_full = dqmDirectoryName(inputDirectory).append(*meName);
    //std::cout << " meName_full = " <<  meName_full << std::endl;

    dqmStore.setCurrentFolder(inputDirectory);
    MonitorElement* meInput = dqmStore.get(meName_full);
    //std::cout << " meInput = " << meInput << std::endl;
    if ( !meInput ) {
      edm::LogError ("copyRecursively") << " Failed to access meName = " << (*meName) << " in DQMStore" 
					<< " --> skipping !!";
      continue;
    }

    //-----------------------------------------------------begin of TauAnalysis specific code 
    //
    // WARNING: the following code is handles MonitorElements of different type differently;
    //          the handling is specific to the TauAnalysis !!
    //
    //           o histograms are scaled then added
    //           o strings are simply copied from input to output
    //             (strings are used to encode run + event numbers in the TauAnalysis)
    //           o integer numbers are added
    //             (integer numbers are used in the TauAnalysis to count the number 
    //              of **unweighted** events passing different stages in the event selection)
    //           o real numbers are scaled then added
    //             (real numbers are used in the TauAnalysis to count the number 
    //              of **weighted** events passing different stages in the event selection)
    //

//--- skip "invalid" MonitorElements
    if ( meInput->kind() == MonitorElement::DQM_KIND_INVALID ) {
      edm::LogWarning ("copyRecursively") << " MonitorElement meName = " << (*meName) << " marked as invalid" 
					  << " --> skipping !!";
      continue;
    }

    if ( meInput->kind() == MonitorElement::DQM_KIND_TH1F      ||
	 meInput->kind() == MonitorElement::DQM_KIND_TH1S      ||
	 meInput->kind() == MonitorElement::DQM_KIND_TH2F      ||
	 meInput->kind() == MonitorElement::DQM_KIND_TH2S      ||
	 meInput->kind() == MonitorElement::DQM_KIND_TH3F      ||
	 meInput->kind() == MonitorElement::DQM_KIND_TPROFILE  ||
	 meInput->kind() == MonitorElement::DQM_KIND_TPROFILE2D ) {
      TH1* histogram = meInput->getTH1();
      //std::cout << " histogram = " << histogram << std::endl;
      if ( !histogram ) {
	edm::LogError ("copyRecursively") << " Failed to access histogram associated to meName = " << (*meName) << " in DQMStore" 
					  << " --> skipping !!";
	continue;
      }

      std::auto_ptr<TH1> clone(dynamic_cast<TH1*>(histogram->Clone()));
      clone->Scale(scaleFactor);

      dqmStore.setCurrentFolder(outputDirectory);

//--- check if output MonitorElement already exists
      int error = 0;
      dqmCheckExistence(dqmStore, outputDirectory, *meName, mode, error);
      if ( error ) continue;

      MonitorElement* meOutput = dqmStore.get(dqmDirectoryName(outputDirectory).append(*meName));
      if ( meOutput ) {
//--- add histogram to outputHistogram
	//std::cout << "--> adding to existing histogram." << std::endl;
	meOutput->getTH1()->Add(clone.get());
      } else {
//--- create new outputHistogram
	//std::cout << "--> registering as new histogram." << std::endl;
	dqmRegisterHistogram(dqmStore, clone.release(), *meName);
      }
    } else if ( meInput->kind() == MonitorElement::DQM_KIND_INT ) {
      int intValue = meInput->getIntValue();

      dqmStore.setCurrentFolder(outputDirectory);

//--- check if output MonitorElement already exists
      int error = 0;
      dqmCheckExistence(dqmStore, outputDirectory, *meName, mode, error);
      if ( error ) continue;

      MonitorElement* meOutput = dqmStore.get(dqmDirectoryName(outputDirectory).append(*meName));
      if ( meOutput ) {
	int intSum = meOutput->getIntValue();
	meOutput->Fill(intSum + intValue);
      } else {
	meOutput = dqmStore.bookInt(*meName);
	meOutput->Fill(intValue);
      }
    } else if ( meInput->kind() == MonitorElement::DQM_KIND_REAL ) {
      double realValue = meInput->getFloatValue();
      realValue *= scaleFactor;

      dqmStore.setCurrentFolder(outputDirectory);

//--- check if output MonitorElement already exists
      int error = 0;
      dqmCheckExistence(dqmStore, outputDirectory, *meName, mode, error);
      if ( error ) continue;

      MonitorElement* meOutput = dqmStore.get(dqmDirectoryName(outputDirectory).append(*meName));
      if ( meOutput ) {
	double realSum = meOutput->getFloatValue();
	meOutput->Fill(realSum + realValue);
      } else {
	meOutput = dqmStore.bookFloat(*meName);
	meOutput->Fill(realValue);
      }
    } else if ( meInput->kind() == MonitorElement::DQM_KIND_STRING ) {
      std::string stringValue = meInput->getStringValue();

      dqmStore.setCurrentFolder(outputDirectory);

//--- check if output MonitorElement already exists
      int error = 0;
      dqmCheckExistence(dqmStore, outputDirectory, *meName, mode, error);
      if ( error ) continue;
      
      //std::cout << " --> calling DQMStore::bookString" << std::endl;
      dqmStore.bookString(*meName, stringValue);
    } else {
      edm::LogError ("copyRecursively") << " MonitorElement meName = " << (*meName) << " of unknown type" 
					<< " --> skipping !!";
      continue;
    }
    //
    //-------------------------------------------------------end of TauAnalysis specific code 

    meInput_copied = true;
  }

//--- call function recursively for all sub-directories
  dqmStore.setCurrentFolder(inputDirectory);
  std::vector<std::string> dirNames = dqmStore.getSubdirs();
  for ( std::vector<std::string>::const_iterator dirName = dirNames.begin();
	dirName != dirNames.end(); ++dirName ) {
    std::string subDirName = dqmSubDirectoryName_merged(inputDirectory, *dirName);
    //std::cout << " subDirName = " << subDirName << std::endl;

    std::string inputDirName_full = dqmDirectoryName(inputDirectory).append(subDirName);
    //std::cout << " inputDirName_full = " << inputDirName_full << std::endl;

    std::string outputDirName_full = dqmDirectoryName(outputDirectory).append(subDirName);
    //std::cout << " outputDirName_full = " << outputDirName_full << std::endl;

    dqmCopyRecursively(dqmStore, inputDirName_full, outputDirName_full, scaleFactor, mode, rmInputDirectory);
  }

//--- delete inputDirectory 
//    (if requested to do so and inputDirectory is **not empty**;
//     otherwise, common parent directories of inputDirectory and outputDirectory might get deleted !!)
  if ( rmInputDirectory && meInput_copied ) dqmStore.rmdir(inputDirectory);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void separateHistogramFromDirectoryName(const std::string& histogramAndDirectoryName, std::string& histogramName, std::string& directoryName)
{
  //std::cout << "<separateHistogramFromDirectoryName>:" << std::endl;

  std::string tempName = histogramAndDirectoryName;

//--- remove DQM root directory from histogram name
  std::string::size_type dqmRootDirectoryPos = tempName.find(dqmRootDirectory);
  if ( dqmRootDirectoryPos != std::string::npos ) {  
    tempName.replace(dqmRootDirectoryPos, dqmRootDirectory.size(), "");  
  }  

  //std::cout << " tempName = " << tempName << std::endl;

//--- extract directory from histogram name
  std::string::size_type lastPos;
  std::string::size_type nextPos = tempName.find(dqmSeparator);  
  do {
    lastPos = nextPos;
    nextPos = tempName.find(dqmSeparator, lastPos + 1);
  } while ( nextPos != std::string::npos );

  histogramName = ( lastPos != std::string::npos ) ? std::string(tempName, lastPos + 1, tempName.length()) : tempName;
  directoryName = ( lastPos != std::string::npos ) ? std::string(tempName, 0, lastPos) : "";

  //std::cout << " histogramName = " << histogramName << std::endl;
  //std::cout << " directoryName = " << directoryName << std::endl;
}

