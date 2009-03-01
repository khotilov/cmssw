// $Id: Numbers.cc,v 1.64 2008/09/05 13:39:13 emanuele Exp $

/*!
  \file Numbers.cc
  \brief Some "id" conversions
  \author B. Gobbo
  \version $Revision: 1.64 $
  \date $Date: 2008/09/05 13:39:13 $
*/

#include <sstream>
#include <iomanip>

#include <DataFormats/EcalDetId/interface/EBDetId.h>
#include <DataFormats/EcalDetId/interface/EEDetId.h>

#include <DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h>
#include <DataFormats/EcalDetId/interface/EcalElectronicsId.h>
#include <DataFormats/EcalDetId/interface/EcalPnDiodeDetId.h>
#include <DataFormats/EcalRawData/interface/EcalDCCHeaderBlock.h>

#include "FWCore/Framework/interface/NoRecordException.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"

#include "DQM/EcalCommon/interface/Numbers.h"

//-------------------------------------------------------------------------

const EcalElectronicsMapping* Numbers::map = 0;

bool Numbers::init = false;

//-------------------------------------------------------------------------

void Numbers::initGeometry( const edm::EventSetup& setup, bool verbose ) {

  if( Numbers::init ) return;

  if ( verbose ) std::cout << "Initializing EcalElectronicsMapping ..." << std::endl;

  Numbers::init = true;

  edm::ESHandle< EcalElectronicsMapping > handle;
  setup.get< EcalMappingRcd >().get(handle);
  Numbers::map = handle.product();

  if ( verbose ) std::cout << "done." << std::endl;

}

//-------------------------------------------------------------------------

int Numbers::iEB( const int ism ) throw( std::runtime_error ) {

  // EB-
  if( ism >=  1 && ism <= 18 ) return( -ism );

  // EB+
  if( ism >= 19 && ism <= 36 ) return( +ism - 18 );

  std::ostringstream s;
  s << "Wrong SM id determination: iSM = " << ism;
  throw( std::runtime_error( s.str() ) );

}

//-------------------------------------------------------------------------

std::string Numbers::sEB( const int ism  ) {

  int ieb = Numbers::iEB( ism );

  std::ostringstream s;
  s << "EB" << std::setw(3) << std::setfill('0')
    << std::setiosflags( std::ios::showpos )
    << std::setiosflags( std::ios::internal )
    << ieb
    << std::resetiosflags( std::ios::showpos )
    << std::resetiosflags( std::ios::internal );
  return( s.str() );

}

//-------------------------------------------------------------------------

int Numbers::iEE( const int ism ) throw( std::runtime_error ) {

  // EE-
  if( ism ==  1 ) return( -7 );
  if( ism ==  2 ) return( -8 );
  if( ism ==  3 ) return( -9 );
  if( ism ==  4 ) return( -1 );
  if( ism ==  5 ) return( -2 );
  if( ism ==  6 ) return( -3 );
  if( ism ==  7 ) return( -4 );
  if( ism ==  8 ) return( -5 );
  if( ism ==  9 ) return( -6 );

  // EE+
  if( ism == 10 ) return( +7 );
  if( ism == 11 ) return( +8 );
  if( ism == 12 ) return( +9 );
  if( ism == 13 ) return( +1 );
  if( ism == 14 ) return( +2 );
  if( ism == 15 ) return( +3 );
  if( ism == 16 ) return( +4 );
  if( ism == 17 ) return( +5 );
  if( ism == 18 ) return( +6 );

  std::ostringstream s;
  s << "Wrong SM id determination: iSM = " << ism;
  throw( std::runtime_error( s.str() ) );

}

//-------------------------------------------------------------------------

EcalSubdetector Numbers::subDet( const EBDetId& id ) {

  return( id.subdet() );

}

//-------------------------------------------------------------------------

EcalSubdetector Numbers::subDet( const EEDetId& id ) {

  return( id.subdet() );

}

//-------------------------------------------------------------------------

EcalSubdetector Numbers::subDet( const EcalTrigTowerDetId& id ) {

  return( id.subDet() );

}

//-------------------------------------------------------------------------

EcalSubdetector Numbers::subDet( const EcalElectronicsId& id ) {

  return( id.subdet() );

}

//-------------------------------------------------------------------------

EcalSubdetector Numbers::subDet( const EcalPnDiodeDetId& id ) {

  return( (EcalSubdetector) id.iEcalSubDetectorId() );

}

//-------------------------------------------------------------------------

EcalSubdetector Numbers::subDet( const EcalDCCHeaderBlock& id ) throw( std::runtime_error ) {

  int idcc = id.id();

  // EE-
  if ( idcc >=  1 && idcc <=  9 ) return( EcalEndcap );

  // EB-/EB+
  if ( idcc >= 10 && idcc <= 45 ) return( EcalBarrel);

  // EE+
  if ( idcc >= 46 && idcc <= 54 ) return( EcalEndcap );

  std::ostringstream s;
  s << "Wrong DCC id: dcc = " << idcc;
  throw( std::runtime_error( s.str() ) );

}

//-------------------------------------------------------------------------

std::string Numbers::sEE( const int ism  ) {

  int iee = Numbers::iEE( ism );

  std::ostringstream s;
  s << "EE" << std::setw(3) << std::setfill('0')
    << std::setiosflags( std::ios::showpos )
    << std::setiosflags( std::ios::internal )
    << iee
    << std::resetiosflags( std::ios::showpos )
    << std::resetiosflags( std::ios::internal );
  return( s.str() );

}

//-------------------------------------------------------------------------

int Numbers::iSM( const int ism, const EcalSubdetector subdet ) throw( std::runtime_error ) {

  if( subdet == EcalBarrel ) {

    // EB-
    if( ism >=  1 && ism <= 18 ) return( ism+18 );

    // EB+
    if( ism >= 19 && ism <= 36 ) return( ism-18 );

    std::ostringstream s;
    s << "Wrong SM id: iSM = " << ism;
    throw( std::runtime_error( s.str() ) );

  } else if( subdet == EcalEndcap ) {

    // EE-
    if( ism >=  1 && ism <=  9 ) return( ism+9 );

    // EE+
    if (ism >= 10 && ism <= 18 ) return( ism-9 );

    std::ostringstream s;
    s << "Wrong SM id: iSM = " << ism;
    throw( std::runtime_error( s.str() ) );

  } else {

    std::ostringstream s;
    s << "Invalid subdetector: subdet = " << subdet;
    throw( std::runtime_error( s.str() ) );

  }

}

//-------------------------------------------------------------------------

int Numbers::iSM( const EBDetId& id ) throw( std::runtime_error ) {

  if( Numbers::map ) {

    EcalElectronicsId eid = Numbers::map->getElectronicsId(id);
    int idcc = eid.dccId();

    // EB-/EB+
    if( idcc >= 10 && idcc <= 45 ) return( idcc - 9 );

    std::ostringstream s;
    s << "Wrong DCC id: dcc = " << idcc;
    throw( std::runtime_error( s.str() ) );

  } else {

    std::ostringstream s;
    s << "ECAL Geometry not available";
    throw( std::runtime_error( s.str() ) );

  }

}

//-------------------------------------------------------------------------

int Numbers::iSM( const EEDetId& id ) throw( std::runtime_error ) {

  if( Numbers::map ) {

    EcalElectronicsId eid = Numbers::map->getElectronicsId(id);
    int idcc = eid.dccId();

    // EE-
    if( idcc >=  1 && idcc <=  9 ) return( idcc );

    // EE+
    if( idcc >= 46 && idcc <= 54 ) return( idcc - 45 + 9 );

    std::ostringstream s;
    s << "Wrong DCC id: dcc = " << idcc;
    throw( std::runtime_error( s.str() ) );

  } else {

    std::ostringstream s;
    s << "ECAL Geometry not available";
    throw( std::runtime_error( s.str() ) );

  }

}

//-------------------------------------------------------------------------

int Numbers::iSM( const EcalTrigTowerDetId& id ) throw( std::runtime_error ) {

  EcalSubdetector subdet = Numbers::subDet( id );

  if( subdet == EcalBarrel ) {

    if( Numbers::map ) {

      int idcc = Numbers::map->DCCid(id);

      // EB-/EB+
      if( idcc >= 10 && idcc <= 45 ) return( idcc - 9 );

      std::ostringstream s;
      s << "Wrong DCC id: dcc = " << idcc;
      throw( std::runtime_error( s.str() ) );

    } else {

      std::ostringstream s;
      s << "ECAL Geometry not available";
      throw( std::runtime_error( s.str() ) );

    }

  } else if( subdet ==  EcalEndcap) {

    if( Numbers::map ) {

      int idcc = Numbers::map->DCCid(id);

      // EE-
      if( idcc >=  1 && idcc <=  9 ) return( idcc );

      // EE+
      if( idcc >= 46 && idcc <= 54 ) return( idcc - 45 + 9 );

      std::ostringstream s;
      s << "Wrong DCC id: dcc = " << idcc;
      throw( std::runtime_error( s.str() ) );

    } else {

      std::ostringstream s;
      s << "ECAL Geometry not available";
      throw( std::runtime_error( s.str() ) );

    }

  } else {

    std::ostringstream s;
    s << "Invalid subdetector: subdet = " << subdet;
    throw( std::runtime_error( s.str() ) );

  }

}

//-------------------------------------------------------------------------

int Numbers::iSM( const EcalElectronicsId& id ) throw( std::runtime_error ) {

  int idcc = id.dccId();

  // EE-
  if( idcc >=  1 && idcc <=  9 ) return( idcc );

  // EB-/EB+
  if( idcc >= 10 && idcc <= 45 ) return( idcc - 9 );

  // EE+
  if( idcc >= 46 && idcc <= 54 ) return( idcc - 45 + 9 );

  std::ostringstream s;
  s << "Wrong DCC id: dcc = " << idcc;
  throw( std::runtime_error( s.str() ) );

}

//-------------------------------------------------------------------------

int Numbers::iSM( const EcalPnDiodeDetId& id ) throw( std::runtime_error ) {

  int idcc = id.iDCCId();

  // EE-
  if( idcc >=  1 && idcc <=  9 ) return( idcc );

  // EB-/EB+
  if( idcc >= 10 && idcc <= 45 ) return( idcc - 9 );

  // EE+
  if( idcc >= 46 && idcc <= 54 ) return( idcc - 45 + 9 );

  std::ostringstream s;
  s << "Wrong DCC id: dcc = " << idcc;
  throw( std::runtime_error( s.str() ) );

}

//-------------------------------------------------------------------------

int Numbers::iSM( const EcalDCCHeaderBlock& id, const EcalSubdetector subdet ) throw( std::runtime_error ) {

  int idcc = id.id();

  // EE-
  if( idcc >=  1 && idcc <=  9 ) return( idcc );

  // EB-/EB+
  if( idcc >= 10 && idcc <= 45 ) return( idcc - 9 );

  // EE+
  if( idcc >= 46 && idcc <= 54 ) return( idcc - 45 + 9 );

  std::ostringstream s;
  s << "Wrong DCC id: dcc = " << idcc;
  throw( std::runtime_error( s.str() ) );

}

//-------------------------------------------------------------------------

int Numbers::iTT( const int ism, const EcalSubdetector subdet, const int i1, const int i2 ) throw( std::runtime_error ) {

  if( subdet == EcalBarrel ) {

    int iet = 1 + ((i1-1)/5);
    int ipt = 1 + ((i2-1)/5);

    return( (ipt-1) + 4*(iet-1) + 1 );

  } else if( subdet == EcalEndcap ) {

    int iz = 0;

    if( ism >=  1 && ism <=  9 ) iz = -1;
    if( ism >= 10 && ism <= 18 ) iz = +1;

    if( EEDetId::validDetId(i1, i2, iz) ) {

      EEDetId id(i1, i2, iz, EEDetId::XYMODE);

      if( Numbers::iSM( id ) != ism ) return( -1 );

      if( Numbers::map ) {

        EcalElectronicsId eid = Numbers::map->getElectronicsId(id);

        return( eid.towerId() );

      } else {

        std::ostringstream s;
        s << "ECAL Geometry not available";
        throw( std::runtime_error( s.str() ) );

      }

    } else {

      return( -1 );

    }

  } else {

    std::ostringstream s;
    s << "Invalid subdetector: subdet = " << subdet;
    throw( std::runtime_error( s.str() ) );

  }

}

//-------------------------------------------------------------------------

int Numbers::iTT( const EcalTrigTowerDetId& id ) throw( std::runtime_error ) {

  EcalSubdetector subdet = Numbers::subDet( id );

  if( subdet == EcalBarrel ) {

    if( Numbers::map ) {

      return( Numbers::map->iTT(id) );

    } else {

      std::ostringstream s;
      s << "ECAL Geometry not available";
      throw( std::runtime_error( s.str() ) );

    }

  } else if( subdet ==  EcalEndcap) {

    if( Numbers::map ) {

      return( Numbers::map->iTT(id) );

    } else {

      std::ostringstream s;
      s << "ECAL Geometry not available";
      throw( std::runtime_error( s.str() ) );

    }

  } else {

    std::ostringstream s;
    s << "Invalid subdetector: subdet = " << subdet;
    throw( std::runtime_error( s.str() ) );

  }

}

//-------------------------------------------------------------------------

int Numbers::TCCid( const EcalTrigTowerDetId& id ) throw( std::runtime_error ) {

  EcalSubdetector subdet = Numbers::subDet( id );

  if( subdet == EcalBarrel ) {

    if( Numbers::map ) {

      return( Numbers::map->TCCid(id) );

    } else {

      std::ostringstream s;
      s << "ECAL Geometry not available";
      throw( std::runtime_error( s.str() ) );

    }

  } else if( subdet ==  EcalEndcap) {

    if( Numbers::map ) {

      return( Numbers::map->TCCid(id) );

    } else {

      std::ostringstream s;
      s << "ECAL Geometry not available";
      throw( std::runtime_error( s.str() ) );

    }

  } else {

    std::ostringstream s;
    s << "Invalid subdetector: subdet = " << subdet;
    throw( std::runtime_error( s.str() ) );

  }

}

//-------------------------------------------------------------------------

std::vector<DetId> Numbers::crystals( const EcalSubdetector subdet, int itcc, int itt ) throw( std::runtime_error ) {


  if( Numbers::map ) {

    EcalSubdetector sub = Numbers::map->subdet(itcc,1);

    if( subdet == sub ) {
      
      return( Numbers::map->ttConstituents( itcc, itt ) );
      
    } else {

      std::vector<DetId> empty;
      empty.clear();
      return empty;

    }
    
  }  else {
    
    std::ostringstream s;
    s << "ECAL Geometry not available";
    throw( std::runtime_error( s.str() ) );
    
  }

}

//-------------------------------------------------------------------------

std::vector<DetId> Numbers::crystals( const EcalTrigTowerDetId& id ) throw( std::runtime_error ) {

  if( Numbers::map ) {

    int itcc = Numbers::map->TCCid(id);
    int itt = Numbers::map->iTT(id);

    return( Numbers::map->ttConstituents( itcc, itt ) );

  } else {

    std::ostringstream s;
    s << "ECAL Geometry not available";
    throw( std::runtime_error( s.str() ) );

  }

}

//-------------------------------------------------------------------------

int Numbers::RtHalf(const EBDetId& id) {

  int ic = id.ic();
  int ie = (ic-1)/20 + 1;
  int ip = (ic-1)%20 + 1;

  if( ie > 5 && ip < 11 ) return 1;

  return 0;

}

//-------------------------------------------------------------------------

int Numbers::RtHalf(const EEDetId& id) {

  int ix = id.ix();

  int ism = Numbers::iSM( id );

  // EE-05
  if ( ism ==  8 && ix > 50 ) return 1;

  // EE+05
  if ( ism == 17 && ix > 50 ) return 1;

  return 0;

}

//-------------------------------------------------------------------------

std::vector<DetId> Numbers::crystals( const EcalElectronicsId& id ) throw( std::runtime_error ) {

  if( Numbers::map ) {

    int idcc = id.dccId();
    int itt = id.towerId();

    return( Numbers::map->dccTowerConstituents( idcc, itt ) );

  } else {

    std::ostringstream s;
    s << "ECAL Geometry not available";
    throw( std::runtime_error( s.str() ) );

  }

}

//-------------------------------------------------------------------------

int Numbers::indexEB( const int ism, const int ie, const int ip ){

  return( (ip-1) + 20*(ie-1) + 1 );

}

//-------------------------------------------------------------------------

int Numbers::indexEE( const int ism, const int ix, const int iy ){

  int iz = 0;

  if( ism >=  1 && ism <=  9 ) iz = -1;
  if( ism >= 10 && ism <= 18 ) iz = +1;

  if( EEDetId::validDetId(ix, iy, iz) ) {

    return( 1000*ix + iy );

  } else {

    return( -1 );

  }

}

//-------------------------------------------------------------------------

int Numbers::icEB( const int ism, const int ie, const int ip ) {

  return( (ip-1) + 20*(ie-1) + 1 );

}

//-------------------------------------------------------------------------

int Numbers::icEE( const int ism, const int ix, const int iy ) throw( std::runtime_error ) {

  int iz = 0;

  if( ism >=  1 && ism <=  9 ) iz = -1;
  if( ism >= 10 && ism <= 18 ) iz = +1;

  if( EEDetId::validDetId(ix, iy, iz) ) {

    EEDetId id(ix, iy, iz, EEDetId::XYMODE);

    if( Numbers::iSM( id ) != ism ) return( -1 );

    if( Numbers::map ) {

      EcalElectronicsId eid = Numbers::map->getElectronicsId(id);

      int vfe = eid.towerId();
      int strip = eid.stripId();
      int channel = eid.xtalId();

      // EE-05 & EE+05
      if( ism == 8 || ism == 17 ) {
        if( vfe > 17 ) vfe = vfe - 7;
      }

      return ( (vfe-1)*25 + (strip-1)*5 + channel );

    } else {

      std::ostringstream s;
      s << "ECAL Geometry not available";
      throw( std::runtime_error( s.str() ) );

    }

  } else {

    return( -1 );

  }

}

//-------------------------------------------------------------------------

int Numbers::ix0EE( const int ism ) {

  if( ism == 1 || ism == 15 ) return( -  5 );
  if( ism == 2 || ism == 14 ) return( +  0 );
  if( ism == 3 || ism == 13 ) return( + 10 );
  if( ism == 4 || ism == 12 ) return( + 40 );
  if( ism == 5 || ism == 11 ) return( + 50 );
  if( ism == 6 || ism == 10 ) return( + 55 );
  if( ism == 7 || ism == 18 ) return( + 50 );
  if( ism == 8 || ism == 17 ) return( + 25 );
  if( ism == 9 || ism == 16 ) return( +  0 );

  return( + 0 );

}

//-------------------------------------------------------------------------

int Numbers::iy0EE( const int ism ) {

  if( ism == 1 || ism == 10 ) return( + 20 );
  if( ism == 2 || ism == 11 ) return( + 45 );
  if( ism == 3 || ism == 12 ) return( + 55 );
  if( ism == 4 || ism == 13 ) return( + 55 );
  if( ism == 5 || ism == 14 ) return( + 45 );
  if( ism == 6 || ism == 15 ) return( + 20 );
  if( ism == 7 || ism == 16 ) return( +  0 );
  if( ism == 8 || ism == 17 ) return( -  5 );
  if( ism == 9 || ism == 18 ) return( +  0 );

  return( + 0 );

}

//-------------------------------------------------------------------------

bool Numbers::validEE( const int ism, const int ix, const int iy ) {

  int iz = 0;

  if( ism >=  1 && ism <=  9 ) iz = -1;
  if( ism >= 10 && ism <= 18 ) iz = +1;

  if( EEDetId::validDetId(ix, iy, iz) ) {

    EEDetId id(ix, iy, iz, EEDetId::XYMODE);

    if( Numbers::iSM( id ) == ism ) return true;

  }

  return false;

}

//-------------------------------------------------------------------------
