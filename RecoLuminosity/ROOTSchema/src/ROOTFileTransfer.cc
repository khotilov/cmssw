#include "RecoLuminosity/ROOTSchema/interface/ROOTFileTransfer.h"
#include <sstream>
#include <iostream>

HCAL_HLX::ROOTFileTransfer::ROOTFileTransfer():fileName_(""),
					       dirName_(""),
					       fileType_("RAW")
{}

HCAL_HLX::ROOTFileTransfer::~ROOTFileTransfer(){}


void HCAL_HLX::ROOTFileTransfer::SetFileType( const std::string &fileType ){
  
  fileType_ = fileType;
}

int HCAL_HLX::ROOTFileTransfer::TransferFile(){

  int errorCode;
  std::stringstream commandLine;

  if( fileName_ == "" ){
    // No File set
    errorCode = -1;
  }else{

    //Transfer File to Offline DB
    commandLine.str(std::string());
    commandLine << "lumiTransferScript.sh " << dirName_ << " " << fileName_ << " " << fileType_;
    std::system(commandLine.str().c_str()); 
    
  }
  return 0;
}

