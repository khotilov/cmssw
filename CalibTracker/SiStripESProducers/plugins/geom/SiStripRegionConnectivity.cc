#include "CalibTracker/SiStripESProducers/plugins/geom/SiStripRegionConnectivity.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"

using namespace sistrip;

SiStripRegionConnectivity::SiStripRegionConnectivity(const edm::ParameterSet& pset) :

  etadivisions_(pset.getUntrackedParameter<unsigned int>("EtaDivisions",10)),
  phidivisions_(pset.getUntrackedParameter<unsigned int>("PhiDivisions",10)),
  etamax_(pset.getUntrackedParameter<double>("EtaMax",2.4))

{
  setWhatProduced(this, &SiStripRegionConnectivity::produceRegionCabling);
}

SiStripRegionConnectivity::~SiStripRegionConnectivity() {}

std::auto_ptr<SiStripRegionCabling> SiStripRegionConnectivity::produceRegionCabling( const SiStripRegionCablingRcd& iRecord ) {

  edm::ESHandle<SiStripDetCabling> detcabling;
  iRecord.getRecord<SiStripDetCablingRcd>().get( detcabling );

  edm::ESHandle<TrackerGeometry> tkgeom;
  iRecord.getRecord<TrackerDigiGeometryRecord>().get( tkgeom );
  
  //here build an object of type SiStripRegionCabling using the information from class SiStripDetCabling **PLUS** the geometry.
  
  //Construct region cabling object
  SiStripRegionCabling* RegionConnections = new SiStripRegionCabling(etadivisions_,phidivisions_,etamax_);
  
  //Construct region cabling map
  SiStripRegionCabling::Cabling regioncabling(etadivisions_*phidivisions_,SiStripRegionCabling::RegionCabling(SiStripRegionCabling::ALLSUBDETS,SiStripRegionCabling::WedgeCabling(SiStripRegionCabling::ALLLAYERS,SiStripRegionCabling::ElementCabling())));
  
  //Loop det cabling
  std::map< uint32_t, std::vector<FedChannelConnection> >::const_iterator idet = detcabling->getDetCabling().begin();
  for (;idet!=detcabling->getDetCabling().end();idet++) {
    if (!idet->first || (idet->first == sistrip::invalid32_)) continue;

    // Check if geom det unit exists
    GeomDetUnit* geom_det = const_cast<GeomDetUnit*>( tkgeom->idToDetUnit(DetId(idet->first)) );
    StripGeomDetUnit* strip_det = dynamic_cast<StripGeomDetUnit*>( geom_det );
    if ( !strip_det ) { continue; }
    
    //Calculate region from geometry
    double eta = tkgeom->idToDet(DetId(idet->first))->position().eta();
    double phi = tkgeom->idToDet(DetId(idet->first))->position().phi().value();
    uint32_t reg = RegionConnections->region(SiStripRegionCabling::Position(eta,phi));
  
    //Find subdet from det-id 
    uint32_t subdet = static_cast<uint32_t>(SiStripRegionCabling::subdetFromDetId(idet->first));
    
    //Find layer from det-id
    uint32_t layer = SiStripRegionCabling::layerFromDetId(idet->first);

    //@@ BELOW IS TEMP FIX TO HANDLE BUG IN DET CABLING
    std::vector<FedChannelConnection> conns = idet->second;
    std::vector<FedChannelConnection>::iterator iconn = conns.begin();
    std::vector<FedChannelConnection>::iterator jconn = conns.end();

    //Update region cabling map
    regioncabling[reg][subdet][layer][idet->first].resize(conns.size());
    for ( ; iconn != jconn; ++iconn ) {
      if ( iconn->apvPairNumber() < conns.size() ) { 
	regioncabling[reg][subdet][layer][idet->first][iconn->apvPairNumber()] = *iconn;
      }
    }
    
  }
  
  //Add map to region cabling object
  RegionConnections->setRegionCabling(regioncabling);
  
  return std::auto_ptr<SiStripRegionCabling>( RegionConnections );
}

