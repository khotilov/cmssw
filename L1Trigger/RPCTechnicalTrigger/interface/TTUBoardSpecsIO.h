// $Id: $
#ifndef CONFIGCODE_TTUBOARDSPECSIO_H 
#define CONFIGCODE_TTUBOARDSPECSIO_H 1

// Include files
#include "CondFormats/RPCObjects/interface/TTUBoardSpecs.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ios>

/** @class TTUBoardSpecsIO TTUBoardSpecsIO.h ConfigCode/TTUBoardSpecsIO.h
 *  
 *
 *  @author Andres Osorio
 *  @date   2008-12-16
 */
class TTUBoardSpecsIO : public TTUBoardSpecs {
public: 
  /// Standard constructor
  TTUBoardSpecsIO( ); 

  virtual ~TTUBoardSpecsIO( ); ///< Destructor

  friend std::istream& operator>>(std::istream &, TTUBoardConfig &);
  
protected:
  
private:

  
};
#endif // CONFIGCODE_TTUBOARDSPECSIO_H
