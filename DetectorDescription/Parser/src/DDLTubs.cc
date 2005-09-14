/***************************************************************************
                          DDLTubs.cpp  -  description
                             -------------------
    begin                : Mon Oct 29 2001
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

// Parser parts
#include "DetectorDescription/Parser/interface/DDLTubs.h"
#include "DetectorDescription/Parser/interface/DDLElementRegistry.h"

// DDCore dependencies
#include "DetectorDescription/Core/interface/DDName.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Base/interface/DDdebug.h"
#include "DetectorDescription/Base/interface/DDException.h"

#include "DetectorDescription/ExprAlgo/interface/ExprEvalSingleton.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <string>
#include <vector>

// Default constructor
DDLTubs::DDLTubs()
{

}

// Default destructor
DDLTubs::~DDLTubs()
{
}

// Upon encountering the end of a Tubs element, call DDCore.
void DDLTubs::processElement (const std::string& name, const std::string& nmspace)
{
  DCOUT_V('P', "DDLTubs::processElement started");

  ExprEvalInterface & ev = ExprEvalSingleton::instance();
  DDXMLAttribute atts = getAttributeSet();

  try {
    if (name == "Tubs")
      {
        DDSolid myTubs = DDSolidFactory::tubs (getDDName(nmspace)
		, ev.eval(nmspace, atts.find("dz")->second)
		, ev.eval(nmspace, atts.find("rMin")->second)
		, ev.eval(nmspace, atts.find("rMax")->second)
		, ev.eval(nmspace, atts.find("startPhi")->second)
		, ev.eval(nmspace, atts.find("deltaPhi")->second)
		);
      }
    else if (name == "Tube")
      {
        DDSolid myTubs = DDSolidFactory::tubs (getDDName(nmspace)
	       , ev.eval(nmspace, atts.find("dz")->second)
	       , ev.eval(nmspace, atts.find("rMin")->second)
	       , ev.eval(nmspace, atts.find("rMax")->second)
	       , 0
	       , 360*deg
	       );
      }
    else if (name == "TruncTubs")
      {      
        DDSolid myTT = DDSolidFactory::truncTubs (getDDName(nmspace)
	      , ev.eval(nmspace, atts.find("zHalf")->second)	
              , ev.eval(nmspace, atts.find("rMin")->second)
	      , ev.eval(nmspace, atts.find("rMax")->second)
	      , 0. // startPhi
	      , ev.eval(nmspace, atts.find("deltaPhi")->second)
	      , ev.eval(nmspace, atts.find("cutAtStart")->second)
	      , ev.eval(nmspace, atts.find("cutAtDelta")->second)
	      , false); // cutInside
      }
    else
      {
	std::string msg = "\nDDLTubs::processElement could not process element.";
	throwError(msg);
      }
  } catch (DDException& e)
    {
      std::string msg(e.what());
      msg += std::string("\nDDLTubs call failed to DDSolidFactory.");
      throwError(msg);
    }
  DDLSolid::setReference(nmspace);

  DCOUT_V('P', "DDLTubs::processElement completed");

}

