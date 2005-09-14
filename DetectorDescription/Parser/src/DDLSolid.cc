/***************************************************************************
                          DDLSolid.cc  -  description
                             -------------------
    begin                : Wed Oct 3 2002
    email                : case@ucdhep.ucdavis.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *           DDDParser sub-component of DDD                                *
 *                                                                         *
 ***************************************************************************/

namespace std{} using namespace std;

// -------------------------------------------------------------------------
// Includes
// -------------------------------------------------------------------------
#include "DetectorDescription/Parser/interface/DDLSolid.h"
#include "DetectorDescription/Parser/interface/DDLElementRegistry.h"
#include "DetectorDescription/Parser/interface/DDLLogicalPart.h"

// DDCore dependencies
#include "DetectorDescription/Core/interface/DDName.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Base/interface/DDdebug.h"

#include "DetectorDescription/ExprAlgo/interface/ExprEvalSingleton.h"

//#include <strstream>
#include <string>

// Default constructor
DDLSolid::DDLSolid()
{
}

// Default desctructor
DDLSolid::~DDLSolid()
{
}

void DDLSolid::setReference (const std::string& nmspace)
{
  // in case it was a BooleanSolid or a ReflectionSolid, clear rSolid.
  DDXMLElement* myrSolid = DDLElementRegistry::getElement("rSolid");
  myrSolid->clear();

  // Make sure Solid elements are in LogicalPart elements.
  if (parent() == "LogicalPart")
    {
      DDXMLElement* refsol = DDLElementRegistry::getElement("rSolid");
      std::vector<std::string> names;
      std::vector<std::string> values;
      names.push_back("name");
      values.push_back(getAttributeSet().find("name")->second);
      refsol->loadAttributes("rSolid", names, values, nmspace);
    }

  // clear THIS solid's values.
  clear();

  DCOUT_V('P', "DDLSolid::setReference completed");
}
