
// This is CSCRecHit2DBuilder.cc

#include <RecoLocalMuon/CSCRecHit/src/CSCRecHit2DBuilder.h>
#include <Geometry/CSCGeometry/interface/CSCChamberSpecs.h>
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/CSCDigi/interface/CSCStripDigi.h>
#include <DataFormats/CSCDigi/interface/CSCWireDigi.h>
//#include <DataFormats/CSCDigi/interface/CSCStripDigiCollection.h>
//#include <DataFormats/CSCDigi/interface/CSCWireDigiCollection.h>
//#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>
#include <RecoLocalMuon/CSCRecHit/src/CSCRecHit2DAlgo.h>
#include <RecoLocalMuon/CSCRecHit/src/CSCRecHit2DBuilderPluginFactory.h>

#include <FWCore/Utilities/interface/Exception.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 

// DOMINIQUE:
#include <RecoLocalMuon/CSCRecHit/src/CSCHitFromStripOnly.h>
#include <RecoLocalMuon/CSCRecHit/src/CSCHitFromWireOnly.h>

CSCRecHit2DBuilder::CSCRecHit2DBuilder( const edm::ParameterSet& ps ) : 
   geom_(0), algos_( std::vector<CSCRecHit2DAlgo*>() ) {
	// Receives ParameterSet percolated down from EDProducer
	
  // Find names of algorithms
  std::vector<std::string> algoNames = ps.getParameter<std::vector<std::string> >("algo_types");

  // Find appropriate ParameterSets
  std::vector<edm::ParameterSet> algoPSets = ps.getParameter<std::vector<edm::ParameterSet> >("algo_psets");

  // Find allocation of algorithm to chamber type
  std::vector<int> algoToType = ps.getParameter<std::vector<int> >("algo_per_chamber_type");

  // How many chamber types do we have? This seems hard-wirable, but what the heck
  int ntypes = ps.getParameter<int>("no_of_chamber_types");
  LogDebug("CSC") << "no. of chamber types = " << ntypes;

// DOMINIQUE:  These are the classes for strip and wire only rechit, if we shall produce them...
  Produce1DHits  = ps.getParameter<bool>("CSCproduce1DHits");

  HitsFromStripOnly_ = new CSCHitFromStripOnly( ps );  // --> How can I pass these parameters: Tim ????
  HitsFromWireOnly_  = new CSCHitFromWireOnly();  // --> no parameters needed


  if ( ntypes <= 0 ) {
	 throw cms::Exception("ParameterSetError") << "No. of chamber types=" << ntypes << " is invalid" << std::endl;
  }
  algos_.resize(ntypes);

  // Trap if we don't have enough parameter sets or haven't assigned an algo to every type
   if ( ( algoNames.size() > algoPSets.size() ) || ( algos_.size() != algoToType.size() )) {
	  throw cms::Exception("ParameterSetError") << 
	      "#algos=" << algoNames.size() << ", #psets=" << algoPSets.size() << 
	      ", # chamber types=" << ntypes << ", dim algoToType=" << algoToType.size() << std::endl;
   }

  // Instantiate the requested algorithm(s) and buffer them
  std::vector<CSCRecHit2DAlgo*> algobuf;

  // Ask factory to build this algorithm, giving it appropriate ParameterSet
  for ( size_t i = 0; i < algoNames.size(); ++i ) {
	CSCRecHit2DAlgo* pAlgo = CSCRecHit2DBuilderPluginFactory::get()->
	        create( algoNames[i], algoPSets[i] );
	algobuf.push_back( pAlgo );
	LogDebug("CSC") << "algorithm [" << i << "] named " << algoNames[i] << " has address " << pAlgo;
  }

	// Register appropriate algorithm for each chamber type
  for ( size_t i = 0; i < algoToType.size(); ++i ) {
	algos_[i] = algobuf[ algoToType[i] - 1 ]; // Care! std::vector index is type-1
	LogDebug("CSC") << "address of algorithm for chamber type " << i << " is " << algos_[i];
  }
}

CSCRecHit2DBuilder::~CSCRecHit2DBuilder() {
  //@@ Can I delete the algorithms I had the factory build??
  //	for ( size_t i = 0; i < algos_.size(); ++i ){
  //		delete algos_[i];
  //	}
}


// DOMINIQUE:  added wire/strip only hits

void CSCRecHit2DBuilder::build( const CSCStripDigiCollection* stripdc,
                                const CSCWireDigiCollection* wiredc,
                                      CSCRecHit2DCollection& oc,
                                      CSCWireHitCollection& woc,
                                      CSCStripHitCollection& soc )
{
  // loop over layers with strip digis and run the local reco on each one with wire digis
  for( CSCStripDigiCollection::DigiRangeIterator it = stripdc->begin(); it != stripdc->end(); ++it )
  {
    const CSCDetId& id = (*it).first;  // DetId for strip region

    // Require both strips and wires in a layer before trying to reconstruct...
    // Any wire digis _in this layer_?

    CSCWireDigiCollection::Range rwired = wiredc->get( id );

    LogDebug("CSC") << "found " << rwired.second - rwired.first << " wire digi(s) in layer";

    // Skip if no wire digis in this layer
    // But for ME11, real wire digis are labelled as belonging to ME1b, so that's where ME1a must look
    // (We try ME1a - above - anyway, because simulated wire digis are labelled as ME1a.)
    if ( rwired.second == rwired.first ) {
      if ( id.station()!=1 || id.ring()!=4 ) continue; // not ME1a, skip to next layer
      // It is ME1a but no wire digis there, so try ME1b...
      int endcap  = id.endcap();
      int chamber = id.chamber();
      int layer   = id.layer();
      CSCDetId idw( endcap, 1, 1, chamber, layer ); // Set idw to same layer in ME1b
      CSCWireDigiCollection::Range rwired = wiredc->get( idw );
      if ( rwired.second == rwired.first ) continue; // Nothing there either, skip
    }

    const CSCStripDigiCollection::Range& rstripd = (*it).second;

    LogDebug("CSC") << "found " << rstripd.second - rstripd.first << " strip digi(s) in layer";

    // Build rechits in one layer by running appropriate algorithm
    const CSCLayer* layer = getLayer( id );
    LogDebug("CSC") << "CSCDetId " << id << " adr CSCLayer=" << layer << "\n";

    //    int chamberType = dynamic_cast<const CSCChamberSpecs&>(layer->type()).chamberType(); 
    int chamberType = layer->chamber()->specs()->chamberType();
    LogDebug("CSC") << "chambertype=" << chamberType;

    CSCRecHit2DAlgo* theAlgo = getAlgo(chamberType);
    LogDebug("CSC") << "use algo w. adr=" << theAlgo;

    std::vector<CSCRecHit2D> rhv = theAlgo->run( id, layer, rstripd, rwired );
   
    // Add the rechits to master collection
	oc.put( id, rhv.begin(), rhv.end() );
  }

  // DOMINIQUE: if we allow to produce the 1-D hit production, continue, otherwise exit
  if ( !Produce1DHits ) return;

  // DOMINIQUE: Strip only hits
  for ( CSCStripDigiCollection::DigiRangeIterator it = stripdc->begin(); it != stripdc->end(); ++it )
  {
      const CSCDetId& id = (*it).first;
      const CSCLayer* layer = getLayer( id );
      const CSCStripDigiCollection::Range& rstripd = (*it).second;
      LogDebug("CSC") << "found " << rstripd.second - rstripd.first << " strip digi(s) in layer";
      
      // Skip if no wire digis in this layer
      if ( rstripd.second == rstripd.first ) continue;

      // This is running CSCRecHit1DFromStripOnly.cc
      std::vector<CSCStripHit> rhv = HitsFromStripOnly_->runStrip( id, layer, rstripd );

      // Add the strip hits to master collection
      soc.put( id, rhv.begin(), rhv.end() );
  }

  // DOMINIQUE:  wire only hits !
  for ( CSCWireDigiCollection::DigiRangeIterator it = wiredc->begin(); it != wiredc->end(); ++it )
  {
      const CSCDetId& id = (*it).first;
      const CSCLayer* layer = getLayer( id );
      const CSCWireDigiCollection::Range rwired = wiredc->get( id );

     // Any wire digis _in this layer_?
     LogDebug("CSC") << "found " << rwired.second - rwired.first << " wire digi(s) in layer";
      
    // Skip if no wire digis in this layer  
    if ( rwired.second == rwired.first ) continue;

// This is running CSCRecHit1DFromWireOnly.cc
      std::vector<CSCWireHit> rhv = HitsFromWireOnly_->runWire( id, layer, rwired );   

      // Add the wire hits to master collection 
      woc.put( id, rhv.begin(), rhv.end() );
  }

}

void CSCRecHit2DBuilder::setGeometry( const CSCGeometry* geom ) {
	geom_ = geom;
}

const CSCLayer* CSCRecHit2DBuilder::getLayer( const CSCDetId& detId ) const {
  if ( !geom_ ) throw cms::Exception("MissingGeometry") << 
     "It don't mean a thing if it ain't got that swing, no, geometry" << std::endl;
  return geom_->layer(detId);
}

CSCRecHit2DAlgo* CSCRecHit2DBuilder::getAlgo( int iChamberType ) {
	//@@ No trapping
  return algos_[ iChamberType - 1 ];
}
