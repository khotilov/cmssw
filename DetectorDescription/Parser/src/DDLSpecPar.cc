/***************************************************************************
                          DDLSpecPar.cc  -  description
                             -------------------
    begin                : Tue Nov 21 2001
    email                : case@ucdhep.ucdavis.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *           DDDParser sub-component of DDD                                *
 *                                                                         *
 ***************************************************************************/

namespace std{} using namespace std;

// Parser parts
#include "DetectorDescription/Parser/interface/DDLSpecPar.h"
#include "DetectorDescription/Parser/interface/DDLElementRegistry.h"
#include "DetectorDescription/Parser/interface/DDXMLElement.h"

// DDCore dependencies
#include "DetectorDescription/Core/interface/DDName.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Core/interface/DDPartSelection.h"
#include "DetectorDescription/Core/interface/DDSpecifics.h"
#include "DetectorDescription/Base/interface/DDdebug.h"
#include "DetectorDescription/Core/interface/DDValue.h"
#include "DetectorDescription/Core/interface/DDValuePair.h"
#include "DetectorDescription/Base/interface/DDException.h"

// CLHEP dependencies
#include "CLHEP/Units/SystemOfUnits.h"
#include "DetectorDescription/ExprAlgo/interface/ExprEvalSingleton.h"

#include <string>
#include <sstream>

// Default constructor
DDLSpecPar::DDLSpecPar()
{
}

// Default destructor
DDLSpecPar::~DDLSpecPar()
{
}

// Process a SpecPar element.  We have to assume that 
// certain things have happened.
void DDLSpecPar::processElement (const string& type, const string& nmspace)
{
  DCOUT_V('P',"DDLSpecPar::processElement started");

  // sends the call to the DDD Core OR does nothing if it is a sub-element

  // What I want to do here is the following:
  // 1.  output PartSelector information.
  // 2.  pass the Path and parameters to DDSpecifics
  // for each of the above, use the name of the SpecPar, since DDL does not
  // provide a name for a PartSelector.

  DDXMLElement* myParameter      = DDLElementRegistry::getElement("Parameter");
  DDXMLElement* myNumeric        = DDLElementRegistry::getElement("Numeric");
  DDXMLElement* myString         = DDLElementRegistry::getElement("String");
  DDXMLElement* myPartSelector   = DDLElementRegistry::getElement("PartSelector");
  DDXMLElement* mySpecParSection = DDLElementRegistry::getElement("SpecParSection");

  // DDPartSelector name comes from DDLSpecPar (this class, there is no analogue to 
  // DDLSpecPar in DDCore)
  vector <string> partsels;
  size_t i;

//    if (getName("name") == "")
//      {
//        cout << "ERROR: no name for SpecPar" << endl;
//        partsels = myPartSelector->getVectorAttribute("path");
//        snames = myParameter->getVectorAttribute("name");
//        cout << "\tParameter Names" << endl;
//        size_t i;
//        for (i = 0; i < snames.size(); i++)
//  	{
//  	  cout << "\t\t" << snames[i] << endl;
//  	}
//        cout << "\tPart Selectors:" << endl;
//        for (i = 0; i < partsels.size(); i++)
//  	{
//  	  cout << "\t\t" << partsels[i] << endl;
//  	}
//      }
//    else 
//      {

  //should i keep this? partsels = myPartSelector->getVectorAttribute("path");
  //otherise I have to do this block...
  for (i = 0; i < myPartSelector->size(); i++)
    partsels.push_back((myPartSelector->getAttributeSet(i).find("path"))->second);

  DDsvalues_type svt;

  // boolean flag to indicate whether the vector<DDValuePair> has been evaluated 
  // using the Evaluator
  typedef map <string, pair<bool,vector<DDValuePair> > > vvvpType;

  vvvpType vvvp;

  /** 08/13/03 doNotEval for Parameter is based on the value of the eval flag.
      For String it is always false and for Numeric it is always true.
      But for "legacy" Parameter, remember, we need to check eval.
      Default is NOT to evaluate.
  **/
  bool doNotEval = true;
  bool doRegex = true;
  {
    // check parent level  
    const DDXMLAttribute & atts = mySpecParSection->getAttributeSet();
    
    if (atts.find("eval") != atts.end() && atts.find("eval")->second == "true")
      doNotEval = false;
    
    if (atts.find("regex") != atts.end() && atts.find("regex")->second == "false")
      doRegex = false;
  }
  {
    // check this level
    const DDXMLAttribute & atts = getAttributeSet();
    
    if (atts.find("eval") != atts.end() && atts.find("eval")->second == "true")
      doNotEval = false;
    else if (atts.find("eval") != atts.end())
      doNotEval = true;
    
    if (atts.find("regex") != atts.end() && atts.find("regex")->second == "false")
      doRegex = false;
    else if (atts.find("regex") != atts.end())
      doRegex = true;
  }

  for (i = 0; i < myParameter->size(); i++)
    {
      const DDXMLAttribute & atts = myParameter->getAttributeSet(i);
      vector <DDValuePair> vvp;
      vvvpType::iterator itv = vvvp.find((atts.find("name")->second));
      if (itv != vvvp.end())
	vvp = itv->second.second;
      
      double tval = 0.0;
      bool isEvaluated = false;

      /** 
	  1.  Check eval flag of each level (SpecParSection, SpecPar
	  and Parameter).
	  2.  Default is the closest specified eval attribute
	  with any value other than "false".
      */
      
      // bool notThis =  doNotEval  myParameter->get(string("eval"), i) != "true";

      // since eval is an optional attribute, we need to check if it exists for
      // the debug statement.  Unfortunately, I see no way to do this internal
      // to the DCOUT_V statement... ts is a "wasted" variable.
      string ts = "** no eval attribute **";
      if (atts.find("eval") != atts.end()) 
	ts = atts.find("eval")->second;

      DCOUT_V('P', string("about to process ") << atts.find("value")->second << string(" eval = ") << ts);

      if ((atts.find("eval") != atts.end() && atts.find("eval")->second !="false")
	  || (atts.find("eval") == atts.end() && !doNotEval))
	{ 
	  tval = ExprEvalSingleton::instance().eval(nmspace, atts.find("value")->second);
	  isEvaluated=true;
	  DCOUT_V('P', string("EVALUATED"));
	}
      else
	{
	  DCOUT_V('P', string("NOT Evaluated"));
	}
      
      DDValuePair vp(atts.find("value")->second, tval);
      vvp.push_back(vp);
      vvvp[atts.find("name")->second] = make_pair(isEvaluated,vvp);
    }

  // Process the String names and values.
  for (i = 0; i < myString->size(); i++)
    {
      const DDXMLAttribute & atts = myString->getAttributeSet(i);
      vector <DDValuePair> vvp;
      vvvpType::iterator itv = vvvp.find(atts.find("name")->second);
      if (itv != vvvp.end())
	vvp = itv->second.second;
      DCOUT_V('P', string("about to process String ") << (atts.find("name")->second) << " = " << (atts.find("value")->second));
      DDValuePair vp(atts.find("value")->second, 0.0);
      vvp.push_back(vp);
      vvvp[atts.find("name")->second] = make_pair(false,vvp);
    }
  
  // Process the Numeric names and values.
  for (i = 0; i < myNumeric->size(); i++)
    {
      const DDXMLAttribute & atts = myNumeric->getAttributeSet(i);
      vector <DDValuePair> vvp;
      vvvpType::iterator itv = vvvp.find(atts.find("name")->second);
      if (itv != vvvp.end())
	vvp = itv->second.second;
      DCOUT_V('P', string("about to process String ") << (atts.find("name")->second) << " = " << (atts.find("value")->second));
      double tval = ExprEvalSingleton::instance().eval(nmspace, atts.find("value")->second);
      DCOUT_V('P', string("EVALUATED"));
      DDValuePair vp(atts.find("value")->second, tval);
      vvp.push_back(vp);
      vvvp[atts.find("name")->second] = make_pair(true,vvp);
    }
  
  for (vvvpType::const_iterator it = vvvp.begin(); it != vvvp.end(); it++)
    {
      DDValue val(it->first, it->second.second);
      bool isEvaluated = it->second.first;
      val.setEvalState(isEvaluated);
      svt[val] = val;
      
    }

  DCOUT_V('p', "DDLSpecPar::processElement\n\tname " << getDDName(nmspace) << "\n\tpartsels.size() = " << myPartSelector->size() << "\n\tsvt " << svt);

  try {
    DDSpecifics ds(getDDName(nmspace), 
		   partsels,
		   svt,
		   doRegex);
  } catch (DDException & e) {
    string msg(e.what());
    msg += "\nDDLSpecPar failed to created a DDSpecifics.";
    throwError(msg);
  }

  myParameter->clear();
  myPartSelector->clear();
  
  // after a SpecPar is done, we can clear
  clear();
  
  DCOUT_V('P',"DDLSpecPar::processElement(...)");
}

