// $Id: Numbers.h,v 1.4 2007/05/21 08:29:43 benigno Exp $

/*!
  \file Numbers.h
  \brief Some "id" conversions
  \author B. Gobbo 
  \version $Revision: 1.4 $
  \date $Date: 2007/05/21 08:29:43 $
*/

#ifndef Numbers_H
#define Numbers_H

#include <string>
#include <stdexcept>

#include <DataFormats/EcalDetId/interface/EBDetId.h>
#include <DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h>
#include <DataFormats/EcalDetId/interface/EcalElectronicsId.h>
#include <DataFormats/EcalDetId/interface/EcalPnDiodeDetId.h>
#include <DataFormats/EcalRawData/interface/EcalDCCHeaderBlock.h>

class Numbers {

 public:

  static int         iEB( int ism ) throw( std::runtime_error );

  static std::string sEB( int ism ) throw( std::runtime_error );

  static int         iEE( int ism ) throw( std::runtime_error );

  static std::string sEE( int ism ) throw( std::runtime_error );

  static int         iSM( const EBDetId&               id );

  static int         iSM( const EcalTrigTowerDetId&    id );

  static int         iSM( const EcalElectronicsId&     id );

  static int         iSM( const EcalPnDiodeDetId&      id );

  static int         iSM( const EcalDCCHeaderBlock&    id );

};

#endif // Numbers_H
