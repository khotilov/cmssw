#include "CSCGeometryESModule.h"
#include "Geometry/CSCGeometryBuilder/src/CSCGeometryBuilderFromDDD.h"
#include "Geometry/CSCGeometry/interface/CSCChamberSpecs.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/MuonNumberingRecord.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"

// Alignments
#include "CondFormats/Alignment/interface/Alignments.h"
#include "CondFormats/Alignment/interface/AlignmentErrors.h"
#include "CondFormats/DataRecord/interface/CSCAlignmentRcd.h"
#include "CondFormats/DataRecord/interface/CSCAlignmentErrorRcd.h"
#include "Geometry/TrackingGeometryAligner/interface/GeometryAligner.h"
#include "Geometry/MuonNumbering/interface/MuonDDDConstants.h"
#include "DataFormats/TrackingRecHit/interface/AlignmentPositionError.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include <memory>

using namespace edm;

CSCGeometryESModule::CSCGeometryESModule(const edm::ParameterSet & p){

  setWhatProduced(this, dependsOn(&CSCGeometryESModule::geometryCallback_) );

  // Choose wire geometry modelling
  // We now _require_ some wire geometry specification in the CSCOrcaSpec.xml file
  // in the DDD Geometry.
  // Default as of transition to CMSSW is to use real values.
  // Alternative is to use pseudo-values which match reasonably closely
  // the calculated geometry values used up to and including ORCA_8_8_1.
  // (This was the default in ORCA.)

  bool useRealWireGeometry =   p.getParameter<bool>("useRealWireGeometry");

  // Suppress strips altogether in ME1a region of ME11?

  bool useOnlyWiresInME1a =    p.getParameter<bool>("useOnlyWiresInME1a");

  // Allow strips in ME1a region of ME11 but gang them?
  // Default is now to treat ME1a with ganged strips (e.g. in clusterizer)

  bool useGangedStripsInME1a = p.getParameter<bool>("useGangedStripsInME1a");

  if ( useGangedStripsInME1a ) useOnlyWiresInME1a = false; // override possible inconsistentcy

  // Use the backed-out offsets that correct the CTI
  bool useCentreTIOffsets = p.getParameter<bool>("useCentreTIOffsets"); 

  // Switch to apply the alignment corrections
  //  applyAlignment_ = p.getParameter<bool>("applyAlignment");

   applyAlignment_ = p.getUntrackedParameter<bool>("applyAlignment", false);

  // Feed these value to where I need them

  CSCChamberSpecs::setUseRealWireGeometry( useRealWireGeometry );
  CSCChamberSpecs::setOnlyWiresInME1a( useOnlyWiresInME1a );
  CSCChamberSpecs::setGangedStripsInME1a( useGangedStripsInME1a );
  CSCChamberSpecs::setUseCentreTIOffsets( useCentreTIOffsets );

}


CSCGeometryESModule::~CSCGeometryESModule(){}


boost::shared_ptr<CSCGeometry>
CSCGeometryESModule::produce(const MuonGeometryRecord & record) {

  //
  // Called whenever the alignments or alignment errors change
  //                                            
  if ( applyAlignment_ ) {
    edm::ESHandle<Alignments> alignments;
    record.getRecord<CSCAlignmentRcd>().get( alignments );
    edm::ESHandle<AlignmentErrors> alignmentErrors;
    record.getRecord<CSCAlignmentErrorRcd>().get( alignmentErrors );
    GeometryAligner aligner;
    aligner.applyAlignments<CSCGeometry>( &(*_cscGeometry),
                                         &(*alignments), &(*alignmentErrors) );
  }

  return _cscGeometry;

}


//______________________________________________________________________________
void CSCGeometryESModule::geometryCallback_( const MuonNumberingRecord& record )
{

  //
  // Called whenever the muon numbering (or ideal geometry) changes
  //
  edm::ESHandle<DDCompactView> cpv;
  edm::ESHandle<MuonDDDConstants> mdc;
  record.getRecord<IdealGeometryRecord>().get(cpv);
  record.get( mdc );
   CSCGeometryBuilderFromDDD builder;
  _cscGeometry = boost::shared_ptr<CSCGeometry>(builder.build(&(*cpv), *mdc));

}

DEFINE_FWK_EVENTSETUP_MODULE(CSCGeometryESModule);
