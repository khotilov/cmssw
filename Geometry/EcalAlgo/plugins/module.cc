#include "Geometry/CaloEventSetup/interface/CaloGeometryEP.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"

template class CaloGeometryEP< EcalBarrelGeometry    > ;
template class CaloGeometryEP< EcalEndcapGeometry    > ;
template class CaloGeometryEP< EcalPreshowerGeometry > ;

typedef CaloGeometryEP< EcalBarrelGeometry > EcalBarrelGeometryEP ;
DEFINE_FWK_EVENTSETUP_MODULE(EcalBarrelGeometryEP);

typedef CaloGeometryEP< EcalEndcapGeometry > EcalEndcapGeometryEP ;
DEFINE_FWK_EVENTSETUP_MODULE(EcalEndcapGeometryEP);

typedef CaloGeometryEP< EcalPreshowerGeometry > EcalPreshowerGeometryEP ;
DEFINE_FWK_EVENTSETUP_MODULE(EcalPreshowerGeometryEP);

#include "Geometry/CaloEventSetup/interface/CaloGeometryDBEP.h"


template class CaloGeometryDBEP< EcalBarrelGeometry    , false> ;
template class CaloGeometryDBEP< EcalEndcapGeometry    , false> ;
template class CaloGeometryDBEP< EcalPreshowerGeometry , false> ;

typedef CaloGeometryDBEP< EcalBarrelGeometry , false> EcalBarrelGeometryFromDBEP ;
DEFINE_FWK_EVENTSETUP_MODULE(EcalBarrelGeometryFromDBEP);

typedef CaloGeometryDBEP< EcalEndcapGeometry , false> EcalEndcapGeometryFromDBEP ;
DEFINE_FWK_EVENTSETUP_MODULE(EcalEndcapGeometryFromDBEP);

typedef CaloGeometryDBEP< EcalPreshowerGeometry , false> EcalPreshowerGeometryFromDBEP ;
DEFINE_FWK_EVENTSETUP_MODULE(EcalPreshowerGeometryFromDBEP);
