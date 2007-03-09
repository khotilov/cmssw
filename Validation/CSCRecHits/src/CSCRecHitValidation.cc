#include "Validation/CSCRecHits/src/CSCRecHitValidation.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>



CSCRecHitValidation::CSCRecHitValidation(const edm::ParameterSet & ps)
: dbe_( edm::Service<DaqMonitorBEInterface>().operator->() ),
  theOutputFile( ps.getParameter<std::string>("outputFile") ),
  theSimHitMap("MuonCSCHits"),
  theCSCGeometry(0),
  the2DValidation(dbe_, ps.getParameter<edm::InputTag>("recHitLabel") ),
  theSegmentValidation(dbe_, ps.getParameter<edm::InputTag>("segmentLabel") )
{
}


CSCRecHitValidation::~CSCRecHitValidation()
{
  if ( theOutputFile.size() != 0 && dbe_ ) dbe_->save(theOutputFile);
}


void CSCRecHitValidation::endJob() {
  if ( theOutputFile.size() != 0 && dbe_ ) dbe_->save(theOutputFile);
}


void CSCRecHitValidation::analyze(const edm::Event&e, const edm::EventSetup& eventSetup)
{
  theSimHitMap.fill(e);

  // find the geometry & conditions for this event
  edm::ESHandle<CSCGeometry> hGeom;
  eventSetup.get<MuonGeometryRecord>().get( hGeom );
  theCSCGeometry = &*hGeom;

  the2DValidation.setGeometry(theCSCGeometry);
  the2DValidation.setSimHitMap(&theSimHitMap);
  theSegmentValidation.setGeometry(theCSCGeometry);
  theSegmentValidation.setSimHitMap(&theSimHitMap);

  the2DValidation.analyze(e, eventSetup);
  theSegmentValidation.analyze(e, eventSetup);

}
