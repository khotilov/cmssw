
/***************************************************************************
                          FIPConfiguration.cc  -  description
                             -------------------
    begin                : Sun Nov 13 2005
    email                : case@ucdhep.ucdavis.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *           FIPConfiguration sub-component of DDD                                *
 *                                                                         *
 ***************************************************************************/

//--------------------------------------------------------------------------
//  Includes
//--------------------------------------------------------------------------
// Parser parts
#include "DetectorDescription/Parser/interface/FIPConfiguration.h"
#include "DetectorDescription/Parser/interface/DDLParser.h"
#include "DetectorDescription/Parser/interface/DDLSAX2ConfigHandler.h"
#include "DetectorDescription/Parser/interface/StrX.h"
#include "DetectorDescription/Base/interface/DDdebug.h"
#include "DetectorDescription/Base/interface/DDException.h"


// Xerces dependencies
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/sax/SAXException.hpp>

// EDM Dependencies
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/EDMException.h"

// STL
#include <string>
#include <iostream>
#include <map>

namespace std{} using namespace std;

//--------------------------------------------------------------------------
//  FIPConfiguration:  Default constructor and destructor.
//--------------------------------------------------------------------------
FIPConfiguration::~FIPConfiguration()
{
  //  parser_->getXMLParser()->setContentHandler(0);  
}

FIPConfiguration::FIPConfiguration() : configHandler_()
{ 
  //  parser_ = DDLParser::instance(); ??
  //  std::cout << "Making a FIPConfiguration with configHandler_ at " << &configHandler_ << std::endl;
}

const std::vector<std::string>&  FIPConfiguration::getFileList(void) const
{
  return files_;
}

const std::vector<std::string>&  FIPConfiguration::getURLList(void) const
{
  return urls_;
}

bool FIPConfiguration::doValidation() const { return configHandler_.doValidation(); }

std::string FIPConfiguration::getSchemaLocation() const { return configHandler_.getSchemaLocation(); }

void FIPConfiguration::dumpFileList(void) const {
  std::cout << "File List:" << std::endl;
  std::cout << "  number of files=" << files_.size() << std::endl;
  for (std::vector<std::string>::const_iterator it = files_.begin(); it != files_.end(); it++)
    std::cout << *it << std::endl;
}

//-----------------------------------------------------------------------
//  Here the Xerces parser is used to process the content of the 
//  configuration file.
//-----------------------------------------------------------------------
int FIPConfiguration::readConfig(const std::string& filename)
{
  DCOUT('P', "FIPConfiguration::ReadConfig(): started");

  // Set the parser to use the handler for the configuration file.
  // This makes sure the Parser is initialized and gets a handle to it.
  // Set these to the flags for the configuration file.
  //parser_->getXMLParser()->setFeature(XMLString::transcode("http://xml.org/sax/features/validation"),true);   // optional
  //parser_->getXMLParser()->setFeature(XMLString::transcode("http://xml.org/sax/features/namespaces"),true);   // optional
  //if (parser_->getXMLParser()->getFeature(XMLString::transcode("http://xml.org/sax/features/validation")) == true)
  // parser_->getXMLParser()->setFeature(XMLString::transcode("http://apache.org/xml/features/validation/dynamic"), true);
//   std::cout << " about to set handler" << std::endl;
//   parser_->getXMLParser()->setContentHandler(&configHandler_);
//   std::cout << " done set handler" << std::endl;
  std::string absoluteFileName (filename);
//   std::cout << " absoluteFileName initialized. " << absoluteFileName << std::endl;
  int toreturn(0);
  try
    {
      XMLPlatformUtils::Initialize();
      //      AlgoInit();
    }

  catch (const XMLException& toCatch)
    {
      std::string e("\nDDLParser(): Error during initialization! Message:");
      e += StrX(toCatch.getMessage()).localForm() + std::string ("\n");
      throw (DDException(e));
    }

  parser_  = XMLReaderFactory::createXMLReader();

  // FIX: Temporarily set validation and namespaces to false always.
  //      Due to ignorance, I did not realize that once set, these can not be
  //      changed for a SAX2XMLReader.  I need to re-think the use of SAX2Parser.
  const XMLCh * tempCh = StrX("http://xml.org/sax/features/validation").xmlChForm();
  parser_->setFeature(tempCh, false);   // optional
  //  delete tempCh[];
  const XMLCh * temp2Ch = StrX("http://xml.org/sax/features/namespaces").xmlChForm();
  parser_->setFeature(temp2Ch, false);   // optional
  //  delete temp2Ch[];
  //  SAX2Parser_->setFeature(XMLString::transcode("http://apache.org/xml/properties/scannerName"), XMLString::transcode("SGXMLScanner"));
  //was not the problem, IGXML was fine!  SAX2Parser_->setProperty(XMLUni::fgXercesScannerName, (void *)XMLUni::fgSGXMLScanner);

  // Specify other parser features, e.g.
  //  SAX2Parser_->setFeature(XMLUni::fgXercesSchemaFullChecking, false);

  parser_->setContentHandler(&configHandler_);

  try {
    edm::FileInPath fp(filename);
    absoluteFileName = fp.fullPath();
    std::cout << "in try..." << fp.fullPath() << std::endl;
  } catch ( const edm::Exception& e ) {
    std::string msg = e.what();
    msg += " caught in readConfig... \nERROR: Could not locate configuration for DetectorDescription " + filename;
    std::cout << msg << std::endl;
    throw DDException(msg);
  }
  std::cout << "Absolute file name is: " << absoluteFileName << std::endl;
  try {
    parser_->parse(absoluteFileName.c_str());
  }
  catch (const XMLException& toCatch) {
    std::cout << "\nXMLException: parsing '" << absoluteFileName << "'\n"
	 << "Exception message is: \n"
	 << std::string(StrX(toCatch.getMessage()).localForm()) << "\n" ;
    toreturn = -1;
  }
  catch (...)
    {
      std::cout << "\nUnexpected exception during parsing: '" << absoluteFileName << "'\n";
      toreturn = 4;
    }
  try {
    const std::vector<std::string>& vURLs = configHandler_.getURLs();
    const std::vector<std::string>& vFiles = configHandler_.getFileNames();
    size_t maxInd = vFiles.size();
    size_t ind = 0;
    for ( ; ind < maxInd ; ind++) {
      edm::FileInPath fp(vURLs[ind] + "/" + vFiles[ind]);
      std::cout << "FileInPath says..." << fp.fullPath() << std::endl;
      files_.push_back(fp.fullPath());
      urls_.push_back("");
    }
  } catch ( const edm::Exception& e ) {
    std::cout << "Caught edm::Exception " << e.what() << std::endl;
  } catch ( ... ) {
    std::cout << "Caught ... exception " << std::endl;
  }


  //   std::vector<std::string> fnames = configHandler_.getFileNames();
  //   std::cout << "there are " << fnames.size() << " files." << std::endl;
  //   for (size_t i = 0; i < fnames.size(); i++)
  //     std::cout << "url=" << configHandler_.getURLs()[i] << " file=" << configHandler_.getFileNames()[i] << std::endl;
  return toreturn;
}
