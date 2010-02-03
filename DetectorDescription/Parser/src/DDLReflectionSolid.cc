/***************************************************************************
                          DDLReflectionSolid.cc  -  description
                             -------------------
    begin                : Mon Mar 4, 2002
    email                : case@ucdhep.ucdavis.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *           DDDParser sub-component of DDD                                *
 *                                                                         *
 ***************************************************************************/



// -------------------------------------------------------------------------
// Includes
// -------------------------------------------------------------------------
#include "DDLReflectionSolid.h"
#include "DDLElementRegistry.h"
#include "DDXMLElement.h"

// DDCore dependencies
#include "DetectorDescription/Core/interface/DDName.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Base/interface/DDdebug.h"



// Default constructor
DDLReflectionSolid::DDLReflectionSolid()
{
}

// Default desctructor
DDLReflectionSolid::~DDLReflectionSolid()
{
}

// Upon starting a ReflectionSolid element, we need to clear all rSolids.
void DDLReflectionSolid::preProcessElement(const std::string& name, const std::string& nmspace)
{
  DDLElementRegistry::getElement("rSolid")->clear();
}

// Upon ending a ReflectionSolid element, call DDCore giving the solid name, and dimensions.
void DDLReflectionSolid::processElement (const std::string& name, const std::string& nmspace)
{
  DCOUT_V('P', "DDLReflectionSolid::processElement started");

  // get solid reference:
  DDXMLElement* myrSolid = DDLElementRegistry::getElement("rSolid");

  if (myrSolid->size() != 1)
    {
      std::cout << "WARNING:  A ReflectionSolid had more than one rSolid.  "
	   << "The first one was used." << std::endl;
      std::cout << "The element to look for is " << getDDName(nmspace) << std::endl;
    }

  DDSolid solid = DDSolid(myrSolid->getDDName(nmspace));
  DDSolid ddreflsol = DDSolidFactory::reflection(getDDName(nmspace), solid);

  DDLSolid::setReference(nmspace);

  DCOUT_V('P', "DDLReflectionSolid::processElement completed");
}
