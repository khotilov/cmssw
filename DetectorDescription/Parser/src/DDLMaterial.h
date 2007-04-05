#ifndef DDLMaterial_H
#define DDLMaterial_H

// -------------------------------------------------------------------------
// Includes
// -------------------------------------------------------------------------
#include "DetectorDescription/Parser/interface/DDXMLElement.h"

#include <string>

/// DDLMaterial processes Box elements.
/** @class DDLMaterial
 * @author Michael Case
 *                                                                       
 *  DDLMaterial.h  -  description
 *  -------------------
 *  begin: Fri Oct 04 2002
 *  email: case@ucdhep.ucdavis.edu
 *
 *  This class currently serves one purpose only.  That is to create a
 *  reference to the most recently created Material, no matter whether
 *  it is an ElementaryMaterial or CompositeMaterial.
 *                                                                         
 */

class DDLMaterial : public DDXMLElement
{
 public:

  /// Constructor
  DDLMaterial();

  /// Destructor
  virtual ~DDLMaterial();

  virtual void setReference (const std::string& nmspace);

};
#endif
