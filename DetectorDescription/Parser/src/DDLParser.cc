/***************************************************************************
                          DDLParser.cc  -  description
                             -------------------
    begin                : Mon Oct 22 2001
    email                : case@ucdhep.ucdavis.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *           DDDParser sub-component of DDD                                *
 *                                                                         *
 ***************************************************************************/

//--------------------------------------------------------------------------
//  Includes
//--------------------------------------------------------------------------
// Parser parts
#include "DetectorDescription/Parser/interface/DDLParser.h"
#include "DetectorDescription/Parser/interface/DDLDocumentProvider.h"
#include "DetectorDescription/Parser/interface/DDLConfiguration.h"
#include "DetectorDescription/Parser/interface/DDLSAX2Handler.h"
#include "DetectorDescription/Parser/interface/DDLSAX2FileHandler.h"
#include "DetectorDescription/Parser/interface/DDLSAX2ConfigHandler.h"
#include "DetectorDescription/Parser/interface/DDLSAX2ExpressionHandler.h"
#include "DetectorDescription/Parser/interface/DDLElementRegistry.h"

// DDCore Dependencies
#include "DetectorDescription/Base/interface/DDdebug.h"
#include "DetectorDescription/Base/interface/DDException.h"
#include "DetectorDescription/Parser/interface/StrX.h"
#include "DetectorDescription/Algorithm/src/AlgoInit.h"

// Xerces dependencies
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/sax/SAXException.hpp>

// Message logger.
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>
#include <iostream>
#include <map>

#include "SealUtil/SealTimer.h"

namespace std{} using namespace std;

//--------------------------------------------------------------------------
//  DDLParser:  Default constructor and destructor.
//--------------------------------------------------------------------------
/// Destructor terminates the XMLPlatformUtils (as required by Xerces)
DDLParser::~DDLParser()
{ 
  // clean up and leave
  XMLPlatformUtils::Terminate();
  DCOUT_V('P', "DDLParser::~DDLParser(): destruct DDLParser"); 
}

/// Constructor initializes XMLPlatformUtils (as required by Xerces.
DDLParser::DDLParser( )  : nFiles_(0) //, handler_(0) (seal::Context* )
{ 
  // Initialize the XML4C2 system
  static seal::SealTimer tddlparser("DDLParser:Construct", false);

  try
    {
      XMLPlatformUtils::Initialize();
      AlgoInit();
    }

  catch (const XMLException& toCatch)
    {
      std::string e("\nDDLParser(): Error during initialization! Message:");
      e += std::string(StrX(toCatch.getMessage()).localForm()) + std::string ("\n");
      throw (DDException(e));
    }

  SAX2Parser_  = XMLReaderFactory::createXMLReader();

  // FIX: Temporarily set validation and namespaces to false always.
  //      Due to ignorance, I did not realize that once set, these can not be
  //      changed for a SAX2XMLReader.  I need to re-think the use of SAX2Parser.
  SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/validation").xmlChForm(), false);   // optional
  SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/namespaces").xmlChForm(), false);   // optional
  //  SAX2Parser_->setFeature(StrX("http://apache.org/xml/properties/scannerName"), XMLString::transcode("SGXMLScanner"));
  //was not the problem, IGXML was fine!  SAX2Parser_->setProperty(XMLUni::fgXercesScannerName, (void *)XMLUni::fgSGXMLScanner);

  // Specify other parser features, e.g.
  //  SAX2Parser_->setFeature(XMLUni::fgXercesSchemaFullChecking, false);

  expHandler_  = new DDLSAX2ExpressionHandler;
  fileHandler_ = new DDLSAX2FileHandler;
  errHandler_  = new DDLSAX2Handler;
  SAX2Parser_->setErrorHandler(errHandler_); 
  SAX2Parser_->setContentHandler(fileHandler_); 

  //  edm::LogInfo ("DDLParser") << "created SAX2XMLReader at memory " << SAX2Parser_ << std::endl;
  
  DCOUT_V('P', "DDLParser::DDLParser(): new (and only) DDLParser"); 
}

// ---------------------------------------------------------------------------
//  DDLSAX2Parser: Implementation of the DDLParser sitting on top of the Xerces
//  C++ XML Parser
//  
//--------------------------------------------------------------------------
// Initialize singleton pointer.
DDLParser* DDLParser::instance_ = 0;

DDLParser* DDLParser::instance()
{

  if ( instance_ == 0 ) {
    instance_ = new DDLParser();
  }

  return instance_;
}

/**  This method allows external "users" to use the current DDLParser on their own.
 *   by giving them access to the SAX2XMLReader.  This may not be a good idea!  The
 *   reason that I 
 */ 
SAX2XMLReader* DDLParser::getXMLParser() { return SAX2Parser_; }

DDLSAX2FileHandler* DDLParser::getDDLSAX2FileHandler() { 
  return fileHandler_; // dynamic_cast < DDLSAX2FileHandler* > (SAX2Parser_->getContentHandler()); 
}

size_t DDLParser::isFound(const std::string& filename)
{
  FileNameHolder::const_iterator it = fileNames_.begin();
  size_t i = 1;
  bool foundFile = false;
  while (it != fileNames_.end() && !foundFile) //  for (; it != fileNames_.end(); it++) 
    {
      if (it->second.first == filename)
	{
	  foundFile = true;
	}
      else i++;
      it++;
    }
  if (foundFile)
    return i;
  return 0;
}

bool DDLParser::isParsed(const std::string& filename)
{
  size_t found = isFound(filename);
  if (found)
    return parsed_[found];
  return false;
}

bool DDLParser::parseOneFile(const std::string& filename, const std::string& url)
{
  size_t foundFile = isFound(filename);

  if (!foundFile || !parsed_[foundFile])
    {

      int fIndex = foundFile;
      //      edm::LogInfo ("DDLParser") << "fIndex= " << fIndex << std::endl;
      if (!foundFile) //parsed_[foundFile])
	{
	  seal::SealTimer tddlfname("DDLParser:"+filename.substr(filename.rfind('/')+1), false);
	  pair <std::string, std::string> pss;
	  pss.first = filename;
	  if (url.size() && url.substr(url.size() - 1, 1) == "/")
	    pss.second = url+filename;
	  else
	    pss.second = url+ "/" + filename;
	  fIndex = nFiles_++;
	  fileNames_[fIndex] = pss;
	  parsed_[fIndex] = false;
	}
      //      edm::LogInfo ("DDLParser") << "fIndex= " << fIndex << std::endl;
      currFileName_ = fileNames_[fIndex].second;
      try
	{
	  //myExpHandler = new DDLSAX2ExpressionHandler;
	  SAX2Parser_->setContentHandler(expHandler_);
	  //	  edm::LogInfo ("DDLParser") << "Parsing: " << fileNames_[fIndex].second << std::endl;
	  LogDebug ("DDLParser") << "Parsing: " << fileNames_[fIndex].second << std::endl;
	  parseFile ( fIndex );

	  //	  delete myExpHandler;
	}
      catch (const XMLException& toCatch) {
	edm::LogError ("DDLParser") << "\nDDLParser::ParseOneFile, PASS1: XMLException while processing files... \n"
	     << "Exception message is: \n"
	     << StrX(toCatch.getMessage()) << "\n" ;
	//	delete myExpHandler;
	XMLPlatformUtils::Terminate();
	throw (DDException("  See XMLException above. "));
      } catch (...) {
	edm::LogError ("DDLParser") << "Some un-caught exception" << endl;
      }
    
      // PASS 2:

      DCOUT_V('P', "DDLParser::ParseOneFile(): PASS2: Just before setting Xerces content and error handlers... ");

      try
	{ 
	  seal::SealTimer t("DDLParser2:"+filename.substr(filename.rfind('/')+1), false);

	  //	  myFileHandler = new DDLSAX2FileHandler;
	  SAX2Parser_->setContentHandler(fileHandler_);

	  // No need to validate (regardless of user's doValidation
	  // because the files have already been validated on the first pass.
	  // This optimization suggested by Martin Liendl.
// 	  SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/validation"), false);   // optional
// 	  SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/namespaces"), false);   // optional
// 	  SAX2Parser_->setFeature(StrX("http://apache.org/xml/features/validation/dynamic"), false);

	  parseFile ( fIndex );
	  parsed_[fIndex] = true;
	  
	  //	  delete myFileHandler;
	}
      catch (const XMLException& toCatch) {
	edm::LogError ("DDLParser") << "\nDDLParser::ParseOneFile, PASS2: XMLException while processing files... \n"
	     << "Exception message is: \n"
	     << StrX(toCatch.getMessage()) << "\n" ;
	//	 delete myFileHandler;
	XMLPlatformUtils::Terminate();
	throw (DDException("  See XMLException above."));
      }
    }
  else // was found and is parsed...
    {
      DCOUT('P', " WARNING: DDLParser::ParseOneFile() file " + filename 
	   + " was already parsed as " + fileNames_[foundFile].second);
      return true;
    }
  return false;
}

std::vector < std::string >  DDLParser::getFileList(void) 
{
  //  edm::LogInfo ("DDLParser") << "start getFileList" << std::endl;
  std::vector<std::string> flist;
  for (FileNameHolder::const_iterator fit = fileNames_.begin(); fit != fileNames_.end(); fit++)
    {
      //      edm::LogInfo ("DDLParser") << "about to push_back " << std::endl;
      //      edm::LogInfo ("DDLParser") << "fit->second.second = " << fit->second.second << std::endl;
      //      edm::LogInfo ("DDLParser") << "fit->second.first = " << fit->second.first << std::endl;
      flist.push_back(fit->second.first); // was .second (mec: 2003:02:19
    }
  //  edm::LogInfo ("DDLParser") << "about to return the list" << std::endl;
  return flist;
}

void DDLParser::dumpFileList(void) {
  edm::LogInfo ("DDLParser") << "File List:" << std::endl;
  for (FileNameHolder::const_iterator it = fileNames_.begin(); it != fileNames_.end(); it++)
    edm::LogInfo ("DDLParser") << it->second.second << std::endl;
}

void DDLParser::dumpFileList(ostream& co) {
  co << "File List:" << std::endl;
  for (FileNameHolder::const_iterator it = fileNames_.begin(); it != fileNames_.end(); it++)
    co << it->second.second << std::endl;
}

// FIX: CLEAN THIS UP!
int DDLParser::parse(const DDLDocumentProvider& dp)
{
  //  edm::LogInfo ("DDLParser") << "Start Parsing.  Validation is set to " << dp.doValidation() << "." << std::endl;
  edm::LogInfo ("DDLParser") << "Start Parsing.  Validation is set off for the time being." << std::endl;
  // prep for pass 1 through DDD XML
  //  DDLSAX2Handler* errHandler = new DDLSAX2Handler;

  //  edm::LogInfo ("DDLParser") << "after setErrorHandler)" << std::endl;
  if (dp.doValidation())
    { 
      //      seal::SealTimer t("DDLParser:Validation");

      //      std::string tval=dp.getSchemaLocation();
      //      const char* tch = tval.c_str();
      DCOUT_V('P', "WARNING:  PARSER VALIDATION IS TURNED OFF REGARDLESS OF <Schema... ELEMENT");
//       SAX2Parser_->setProperty(StrX("http://apache.org/xml/properties/schema/external-schemaLocation"),XMLString::transcode(tch));
//       SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/validation"), true);   // optional
//       SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/namespaces"), true);   // optional
//       SAX2Parser_->setFeature(StrX("http://apache.org/xml/features/validation/dynamic"), true);
    }
  else
    {
//       seal::SealTimer t("DDLParser:NoValidation");
//       SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/validation"), false);   // optional
//       SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/namespaces"), false);   // optional
//       SAX2Parser_->setFeature(StrX("http://apache.org/xml/features/validation/dynamic"), false);
    }

  //  This need be only done once, so might as well to it here.
  //  DDLSAX2ExpressionHandler* myExpHandler = new DDLSAX2ExpressionHandler;
  size_t fileIndex = 0;
  std::vector<std::string> fullFileName;

  for (; fileIndex < (dp.getFileList()).size(); fileIndex++)
    { 
      std::string ts = dp.getURLList()[fileIndex];
      std::string tf = dp.getFileList()[fileIndex];
      if ( ts.size() > 0 ) {
	if ( ts[ts.size() - 1] == '/') {
	  fullFileName.push_back( ts + tf );
	} else {
	  fullFileName.push_back( ts + "/" + tf );
	}
      } else {
	fullFileName.push_back( tf );
      }
      //      edm::LogInfo ("DDLParser") << "full file name" << fullFileName[fileIndex] << std::endl;
    }

    seal::SealTimer * t = new seal::SealTimer("DDLParser1", false);

    for (std::vector<std::string>::const_iterator fnit = fullFileName.begin(); 
	 fnit != fullFileName.end();
	 fnit++)
      {
	size_t foundFile = isFound(expHandler_->extractFileName( *fnit )); 
	
	if (!foundFile)
	  {
	    pair <std::string, std::string> pss;
	    pss.first = expHandler_->extractFileName( *fnit );
	    pss.second = *fnit;
	    fileNames_[nFiles_++] = pss;
	    parsed_[nFiles_ - 1]=false;
	  }
      }
    
  // Start processing the files found in the config file.
  
  // PASS 1:  This was added later (historically) to implement the DDD
  // requirement for Expressions.
  DCOUT('P', "DDLParser::parse(): PASS1: Just before setting Xerces content and error handlers... ");
  
  
  try
    {
      SAX2Parser_->setContentHandler(expHandler_);
      //edm::LogInfo ("DDLParser") << "1st PASS: Parsing " << fileNames_.size() << " files." << std::endl;
      //std::vector<std::string>::const_iterator currFile = fileNames_.begin(); currFile != fileNames_.end(); currFile++)
      for (size_t i = 0; i < fileNames_.size(); i++)
	{
	  if (!parsed_[i])
	    {
	      //edm::LogInfo ("DDLParser") << "Parsing: " << fileNames_[i].second << std::endl;
	      parseFile(i);
	    }
	}
      expHandler_->dumpElementTypeCounter();
    }
  catch (const XMLException& toCatch) {
    edm::LogInfo ("DDLParser") << "\nPASS1: XMLException while processing files... \n"
	 << "Exception message is: \n"
	 << StrX(toCatch.getMessage()) << "\n" ;
    //    delete myExpHandler;
    XMLPlatformUtils::Terminate();
    // FIX use this after DEPRECATED stuff removed:    throw(DDException("See XML Exception above"));
    return -1;
  }
  catch (DDException& e) {
    edm::LogInfo ("DDLParser") << "unexpected: " << std::endl;
    edm::LogInfo ("DDLParser") << e << std::endl;
    //    delete myExpHandler;
    // FIX use this after DEPRECATED stuff removed    throw(e);
    return 4;
  }
  //  delete myExpHandler;
  delete t;
  // PASS 2:

  DCOUT('P', "DDLParser::parse(): PASS2: Just before setting Xerces content and error handlers... ");
  t = new seal::SealTimer("DDLParser2", false);
  //  DDLSAX2FileHandler* myFileHandler(0);
  try
    {
      //      myFileHandler = new DDLSAX2FileHandler;
      SAX2Parser_->setContentHandler(fileHandler_);
      // *T no need to re-do, all the same!

      // No need to validate (regardless of user's doValidation
      // because the files have already been validated on the first pass.
      // This optimization suggested by Martin Listd::endl.
//       SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/validation"), false);   // optional
//       SAX2Parser_->setFeature(StrX("http://xml.org/sax/features/namespaces"), false);   // optional
//       SAX2Parser_->setFeature(StrX("http://apache.org/xml/features/validation/dynamic"), false);


      // Process files again.
      //edm::LogInfo ("DDLParser") << "Parsing " << fileNames_.size() << " files." << std::endl;
      for (size_t i = 0; i < fileNames_.size(); i++)
	//std::vector<std::string>::const_iterator currFile = fileNames_.begin(); currFile != fileNames_.end(); currFile++)
	{
	  parseFile(i);
	  parsed_[i] = true;
	  pair<std::string, std::string> namePair = fileNames_[i];
	  //	  edm::LogInfo ("DDLParser") << "Completed parsing file " << namePair.second << std::endl;
	  LogDebug ("DDLParser") << "Completed parsing file " << namePair.second << std::endl;
	}
      //myFileHandler->dumpElementTypeCounter();

      //      delete myFileHandler;
    }
  catch (DDException& e) {
    std::string s(e.what());
    s+="\n\t see above:  DDLParser::parse  Exception " ;
    edm::LogError ("DDLParser") << s << std::endl;
    return -1;
    //    throw(e);
  }
  catch (const XMLException& toCatch) {
    edm::LogError ("DDLParser") << "\nPASS2: XMLException while processing files... \n"
	 << "Exception message is: \n"
	 << StrX(toCatch.getMessage()) << "\n" ;
    XMLPlatformUtils::Terminate();
    return -1;
  }
  delete t;
  // FIX remove after DEPRECATED stuff removed
  return 0;
}

void DDLParser::parseFile(const int& numtoproc) 
{
  
  if (!parsed_[numtoproc])
    {
      const std::string & fname = fileNames_[numtoproc].second;
      seal::SealTimer t("DDLParser:"+fname.substr(fname.rfind('/')+1), false);

      try
	{
	  currFileName_ = fname;
	  SAX2Parser_->parse(currFileName_.c_str());
	}
      catch (const XMLException& toCatch)
	{
	  std::string e("\nWARNING: DDLParser::parseFile, File: '");
	  e += currFileName_ + "'\n"
	    + "Exception message is: \n"
	    + std::string(StrX(toCatch.getMessage()).localForm()) + "\n";
	  throw(DDException(e));
	}
      catch (DDException& e)
	{
	  std::string s(e.what());
	  s += "\nERROR: Unexpected exception during parsing: '"
	    + currFileName_ + "'\n"
	    + "Exception message is shown above.\n";
	  throw(DDException(e));
	} catch (...) {
	  edm::LogError ("DDLParser") << "Another un-caught exception!" << endl;
	}
    }
  else
    {
      DCOUT('P', "\nWARNING: File " + fileNames_[numtoproc].first 
	   + " has already been processed as " + fileNames_[numtoproc].second);
    }
}

// Return the name of the Current file being processed by the parser.
std::string DDLParser::getCurrFileName()
{
  return currFileName_;
}

