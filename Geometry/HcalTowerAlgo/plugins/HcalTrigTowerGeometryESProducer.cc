#include "HcalTrigTowerGeometryESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include <memory>

HcalTrigTowerGeometryESProducer::HcalTrigTowerGeometryESProducer( const edm::ParameterSet & config )
{
    setWhatProduced( this );
}

HcalTrigTowerGeometryESProducer::~HcalTrigTowerGeometryESProducer( void ) 
{}

boost::shared_ptr<HcalTrigTowerGeometry>
HcalTrigTowerGeometryESProducer::produce( const CaloGeometryRecord & iRecord )
{
    edm::ESHandle<HcalTopology> hcalTopology;
    iRecord.getRecord<IdealGeometryRecord>().get( hcalTopology );

    m_hcalTrigTowerGeom =
	boost::shared_ptr<HcalTrigTowerGeometry>( new HcalTrigTowerGeometry( &*hcalTopology));

    return m_hcalTrigTowerGeom;
}

DEFINE_FWK_EVENTSETUP_MODULE( HcalTrigTowerGeometryESProducer );
