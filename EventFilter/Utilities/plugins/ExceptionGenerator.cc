#include "ExceptionGenerator.h"

#include <iostream>
#include <typeinfo>
#include <map>
#include <sstream>

#include "xgi/Method.h"
#include "xgi/Utils.h"

#include "cgicc/Cgicc.h"
#include "cgicc/FormEntry.h"
#include "cgicc/HTMLClasses.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "boost/tokenizer.hpp"

#include <stdio.h>

using namespace std;

namespace evf{

    const std::string ExceptionGenerator::menu[menu_items] =  
      {"Sleep x ms", "SleepForever", "Cms Exception", "Exit with error", "Abort", "Unknown Exception", "Endless loop", "Generate Error Message" };

    ExceptionGenerator::ExceptionGenerator( const edm::ParameterSet& pset) : 
      ModuleWeb("ExceptionGenerator"), actionRequired_(false), actionId_(-1)
    {
      

    }
    void ExceptionGenerator::analyze(const edm::Event & e, const edm::EventSetup& c)
    {
      if(actionRequired_) 
	{
	  int ind = 0; 
	  int step = 1; 
	  switch(actionId_)
	    {
	    case 0:
	      ::usleep(intqualifier_*1000);
	      break;
	    case 1:
	      ::sleep(0xFFFFFFF);
	      break;
	    case 2:
	      throw cms::Exception(qualifier_) << "This exception was generated by the ExceptionGenerator";
	      break;
	    case 3:
	      exit(-1);
	      break;
	    case 4:
	      abort();
	      break;
	    case 5:
	      throw qualifier_;
	      break;
	    case 6:
	      while(1){ind+=step; if(ind>1000000) step = -1; if(ind==0) step = 1;}
	      break;
	    case 7:
	      edm::LogError("TestErrorMessage") << qualifier_;
	    }
	}
    }
    
    void ExceptionGenerator::endLuminosityBlock(edm::LuminosityBlock const &lb, edm::EventSetup const &es)
    {

    }
    
    void ExceptionGenerator::defaultWebPage(xgi::Input *in, xgi::Output *out)
    {

      std::string path;
      std::string urn;
      std::string mname;
      std::string query;
      try 
	{
	  cgicc::Cgicc cgi(in);
	  if ( xgi::Utils::hasFormElement(cgi,"exceptiontype") )
	    {
	      actionId_ = xgi::Utils::getFormElement(cgi, "exceptiontype")->getIntegerValue();
	      if(actionId_>0)
		qualifier_ = xgi::Utils::getFormElement(cgi, "qualifier")->getValue();
	      else
		intqualifier_ =  xgi::Utils::getFormElement(cgi, "qualifier")->getIntegerValue();
	      actionRequired_ = true;
	    }
	  if ( xgi::Utils::hasFormElement(cgi,"module") )
	    mname = xgi::Utils::getFormElement(cgi, "module")->getValue();
	  cgicc::CgiEnvironment cgie(in);
	  if(original_referrer_ == "")
	    original_referrer_ = cgie.getReferrer();
	  path = cgie.getPathInfo();
	  query = cgie.getQueryString();
	}
      catch (const std::exception & e) 
	{
	  // don't care if it did not work
	}
      
      
      using std::endl;
      *out << "<html>"                                                   << endl;
      *out << "<head>"                                                   << endl;


      *out << "<STYLE type=\"text/css\"> #T1 {border-width: 2px; border: solid blue; text-align: center} </STYLE> "                                      << endl; 
      *out << "<link type=\"text/css\" rel=\"stylesheet\"";
      *out << " href=\"/" <<  urn
	   << "/styles.css\"/>"                   << endl;

      *out << "<title>" << moduleName_
	   << " MAIN</title>"                                            << endl;

      *out << "</head>"                                                  << endl;
      *out << "<body onload=\"loadXMLDoc()\">"                           << endl;
      *out << "<table border=\"0\" width=\"100%\">"                      << endl;
      *out << "<tr>"                                                     << endl;
      *out << "  <td align=\"left\">"                                    << endl;
      *out << "    <img"                                                 << endl;
      *out << "     align=\"middle\""                                    << endl;
      *out << "     src=\"/evf/images/systemerror.jpg\""	         << endl;
      *out << "     alt=\"main\""                                        << endl;
      *out << "     width=\"90\""                                        << endl;
      *out << "     height=\"64\""                                       << endl;
      *out << "     border=\"\"/>"                                       << endl;
      *out << "    <b>"                                                  << endl;
      *out <<             moduleName_                                    << endl;
      *out << "    </b>"                                                 << endl;
      *out << "  </td>"                                                  << endl;
      *out << "  <td width=\"32\">"                                      << endl;
      *out << "    <a href=\"/urn:xdaq-application:lid=3\">"             << endl;
      *out << "      <img"                                               << endl;
      *out << "       align=\"middle\""                                  << endl;
      *out << "       src=\"/hyperdaq/images/HyperDAQ.jpg\""             << endl;
      *out << "       alt=\"HyperDAQ\""                                  << endl;
      *out << "       width=\"32\""                                      << endl;
      *out << "       height=\"32\""                                     << endl;
      *out << "       border=\"\"/>"                                     << endl;
      *out << "    </a>"                                                 << endl;
      *out << "  </td>"                                                  << endl;
      *out << "  <td width=\"32\">"                                      << endl;
      *out << "  </td>"                                                  << endl;
      *out << "  <td width=\"32\">"                                      << endl;
      *out << "    <a href=\"" << original_referrer_  << "\">"           << endl;
      *out << "      <img"                                               << endl;
      *out << "       align=\"middle\""                                  << endl;
      *out << "       src=\"/evf/images/spoticon.jpg\""			 << endl;
      *out << "       alt=\"main\""                                      << endl;
      *out << "       width=\"32\""                                      << endl;
      *out << "       height=\"32\""                                     << endl;
      *out << "       border=\"\"/>"                                     << endl;
      *out << "    </a>"                                                 << endl;
      *out << "  </td>"                                                  << endl;
      *out << "</tr>"                                                    << endl;
      *out << "</table>"                                                 << endl;

      *out << "<hr/>"                                                    << endl;

  
      *out << cgicc::form().set("method","GET").set("action", path ) 
	   << std::endl;
      boost::char_separator<char> sep("&");
      boost::tokenizer<boost::char_separator<char> > tokens(query, sep);
      for (boost::tokenizer<boost::char_separator<char> >::iterator tok_iter = tokens.begin();
	   tok_iter != tokens.end(); ++tok_iter){
	size_t pos = (*tok_iter).find_first_of("=");
	if(pos != std::string::npos){
	  std::string first  = (*tok_iter).substr(0    ,                        pos);
	  std::string second = (*tok_iter).substr(pos+1, (*tok_iter).length()-pos-1);
	  *out << cgicc::input().set("type","hidden").set("name",first).set("value", second) 
	       << std::endl;
	}
      }

      *out << "Select   "						 << endl;
      *out << cgicc::select().set("name","exceptiontype")     << std::endl;
      char istring[2];

      for(int i = 0; i < menu_items; i++)
	{
	  sprintf(istring,"%d",i);
	  *out << cgicc::option().set("value",istring) << menu[i] << cgicc::option()       << std::endl;
	}
      *out << cgicc::select() 	     << std::endl;
      *out << "<br>"                                                     << endl;
      *out << "Qualifier"      						 << endl;
      *out << cgicc::input().set("type","text").set("name","qualifier")  	     << std::endl;
      *out << cgicc::input().set("type","submit").set("value","Do It !")  	     << std::endl;
      *out << cgicc::form()						   << std::endl;  

      *out << "</body>"                                                  << endl;
      *out << "</html>"                                                  << endl;
    }

} // end namespace evf
