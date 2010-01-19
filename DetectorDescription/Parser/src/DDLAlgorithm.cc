/**************************************************************************
      DDLAlgorithm.cc  -  description
                             -------------------
    begin                : Saturday November 29, 2003
    email                : case@ucdhep.ucdavis.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *           DDDParser sub-component of DDD                                *
 *                                                                         *
 ***************************************************************************/



// Parser parts
#include "DDLAlgorithm.h"
#include "DDLVector.h"
#include "DDLMap.h"
#include "DDLElementRegistry.h"

// DDCore dependencies
#include "DetectorDescription/Core/interface/DDName.h"
#include "DetectorDescription/Base/interface/DDdebug.h"
#include "DetectorDescription/Algorithm/interface/DDAlgorithm.h"

#include "DetectorDescription/Algorithm/interface/DDAlgorithmHandler.h"

// CLHEP dependencies
#include "DetectorDescription/ExprAlgo/interface/ExprEvalSingleton.h"

#include <sstream>

DDLAlgorithm::DDLAlgorithm( DDLElementRegistry* myreg ) : DDXMLElement(myreg)
{
}

DDLAlgorithm::~DDLAlgorithm()
{
}

void DDLAlgorithm::preProcessElement (const std::string& name, const std::string& nmspace, DDCompactView& cpv)
{
  myRegistry_->getElement("Vector")->clear();
}

void DDLAlgorithm::processElement (const std::string& name, const std::string& nmspace, DDCompactView& cpv)
{
  DCOUT_V('P',"DDLAlgorithm::processElement started");

  DDXMLElement* myNumeric        = myRegistry_->getElement("Numeric");
  DDXMLElement* myString         = myRegistry_->getElement("String");
  DDXMLElement* myVector         = myRegistry_->getElement("Vector");
  DDXMLElement* myMap            = myRegistry_->getElement("Map");
  DDXMLElement* myrParent        = myRegistry_->getElement("rParent");

  DDName algoName(getDDName(nmspace));  
  DDLogicalPart lp(DDName(myrParent->getDDName(nmspace)));
  DDXMLAttribute atts;

  // handle all Numeric elements in the Algorithm.
  DDNumericArguments nArgs;
  size_t i = 0;
  for (; i < myNumeric->size(); ++i)
    {
      atts = myNumeric->getAttributeSet(i);
      nArgs[atts.find("name")->second] = ExprEvalSingleton::instance().eval(nmspace, atts.find("value")->second);
    }

  DDStringArguments sArgs;
  for (i = 0; i < myString->size(); ++i)
    {
      atts = myString->getAttributeSet(i);
      sArgs[atts.find("name")->second] = atts.find("value")->second;
    }

  DDAlgorithmHandler handler;
  atts = getAttributeSet();
  DDLVector* tv= dynamic_cast<DDLVector*> (myVector);
  DDLMap* tm= dynamic_cast<DDLMap*> (myMap);
  handler.initialize( algoName, lp, nArgs, tv->getMapOfVectors(), tm->getMapOfMaps(), sArgs, tv->getMapOfStrVectors());
  //  handler.execute( cpv );
  DDPositioner pos ( &cpv );
  handler.execute( pos );
  //  handler.execute();

  // clear used/referred to elements.
  myString->clear();
  myNumeric->clear();
  myVector->clear();
  myMap->clear();
  myrParent->clear();
  clear();

  DCOUT_V('P',"DDLAlgorithm::processElement(...)");
}

