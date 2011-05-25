
#include "RecoPixelVertexing/PixelTriplets/interface/QuadrupletSeedMerger.h"


/***

QuadrupletSeedMerger

* merge triplet into quadruplet seeds
  for either SeedingHitSets or TrajectorySeeds

* method: 
    QuadrupletSeedMerger::mergeTriplets( const OrderedSeedingHits&, const edm::EventSetup& ) const
  for use in HLT with: RecoPixelVertexing/PixelTrackFitting/src/PixelTrackReconstruction.cc
  contains the merging functionality

* method:
    QuadrupletSeedMerger::mergeTriplets( const TrajectorySeedCollection&, TrackingRegion&, EventSetup& ) const
  is a wrapper for the former, running on TrajectorySeeds
  for use in iterative tracking with: RecoTracker/TkSeedGenerator/plugins/SeedGeneratorFromRegionHitsEDProducer

***/




/// a binary to unary
/// for remove_if, see QuadrupletSeedMerger::mergeTriplets
struct seedmergerIsEqualLayer : public std::binary_function< SeedMergerPixelLayer, std::pair<SeedMergerPixelLayer,SeedMergerPixelLayer>, bool > {
  bool operator () ( const SeedMergerPixelLayer& l, const std::pair<SeedMergerPixelLayer,SeedMergerPixelLayer>& p ) const {
    return( p.first.getName() == l.getName() || p.second.getName() == l.getName() );
  }
};



///
/// swo comp for sorting quadruplet hits
/// as TransientTrackingRecHit::ConstRecHitPointer
/// according to radius
///
bool QuadrupletSeedMerger::isGreaterHit( const TransientTrackingRecHit::ConstRecHitPointer& a,
					 const TransientTrackingRecHit::ConstRecHitPointer& b  )  {
  const std::pair<const GeomDet*,const GeomDet*> geomDets( theTrackerGeometry_->idToDet( a->hit()->geographicalId() ),
							   theTrackerGeometry_->idToDet( b->hit()->geographicalId() )  );
  const std::pair<double,double> radii( sqrt( pow( geomDets.first->toGlobal(  a->hit()->localPosition() ).x(), 2. ) +
					      pow( geomDets.first->toGlobal(  a->hit()->localPosition() ).y(), 2. )   ),
					sqrt( pow( geomDets.second->toGlobal( b->hit()->localPosition() ).x(), 2. ) +
					      pow( geomDets.second->toGlobal( b->hit()->localPosition() ).y(), 2. )   )  );
  return( radii.first < radii.second );
}



///
///
///
QuadrupletSeedMerger::QuadrupletSeedMerger( const edm::EventSetup& es ) {

  // copy geometry
  es.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry_ );

  // by default, do not..
  // ..merge triplets
  isMergeTriplets_ = false;
  // ..add remaining triplets
  isAddRemainingTriplets_ = false;

  // default is the layer list from plain quadrupletseedmerging_cff
  // unless configured contrarily via setLayerListName()
  layerListName_ = std::string( "PixelSeedMergerQuadruplets" );

}



///
///
///
QuadrupletSeedMerger::~QuadrupletSeedMerger() {

}



///
/// merge triplets into quadruplets
/// INPUT: OrderedSeedingHits
/// OUTPUT: SeedingHitSets
///
/// this method is used in RecoPixelVertexing/PixelTrackFitting/src/PixelTrackReconstruction.cc
/// and contains the basic merger functionality
///
const std::vector<SeedingHitSet> QuadrupletSeedMerger::mergeTriplets( const OrderedSeedingHits& inputTriplets, const edm::EventSetup& es ) {


  // the list of layers on which quadruplets should be formed
  edm::ESHandle<SeedingLayerSetsBuilder> layerBuilder;
  es.get<TrackerDigiGeometryRecord>().get( layerListName_.c_str(), layerBuilder );
  theLayerSets_ = layerBuilder->layers( es ); // this is a vector<vector<SeedingLayer> >

  // make a working copy of the input triplets
  // to be able to remove merged triplets
  std::list<SeedingHitSet> tripletList;
  for( unsigned int it = 0; it < inputTriplets.size(); ++it ) {
    tripletList.push_back( inputTriplets[it] );
  }
  const unsigned int nInputTriplets = tripletList.size();

  // the output
  std::vector<SeedingHitSet> theResult;


  // check if the input is all triplets
  // (code is also used for pairs..)
  // if not, copy the input to the output & return
  bool isAllTriplets = true;
  for( std::list<SeedingHitSet>::iterator triplet = tripletList.begin(); triplet != tripletList.end(); ++triplet ) {
    if( triplet->size() != 3 ) isAllTriplets = false;
  }

  if( !isAllTriplets && isMergeTriplets_ )
    std::cout << "[QuadrupletSeedMerger::mergeTriplets] (in HLT) ** bailing out since non-triplets in input." << std::endl;

  if( !isAllTriplets || !isMergeTriplets_ ) {
    for( std::list<SeedingHitSet>::const_iterator it = tripletList.begin(); it != tripletList.end(); ++it ) {
      theResult.push_back( *it );
    }

    return theResult;
  }


  // loop all possible 4-layer combinations
  // as specified in python/quadrupletseedmerging_cff.py
  for( ctfseeding::SeedingLayerSets::const_iterator lsIt = theLayerSets_.begin(); lsIt < theLayerSets_.end(); ++lsIt ) {
    
    // fill a vector with the layers in this set
    std::vector<SeedMergerPixelLayer> currentLayers;
    for( ctfseeding::SeedingLayers::const_iterator layIt = lsIt->begin(); layIt < lsIt->end(); ++layIt ) {
      currentLayers.push_back( SeedMergerPixelLayer( layIt->name() ) );
    }


    // loop all pair combinations of these 4 layers;
    // strategy is to look for shared hits on such a pair and
    // then merge the remaining hits in both triplets if feasible
    // (SeedingLayers is a vector<SeedingLayer>)
    for( std::vector<SeedMergerPixelLayer>::const_iterator shared1 = currentLayers.begin(); shared1 < currentLayers.end(); ++shared1 ) {
      for( std::vector<SeedMergerPixelLayer>::const_iterator shared2 = shared1 + 1; shared2 < currentLayers.end(); ++shared2 ) {

	// the two layers that the triplets share hits on
	std::pair<SeedMergerPixelLayer, SeedMergerPixelLayer> sharedLayers( *shared1, *shared2 );
	
	// now the other two layers: make a copy of them all
	std::vector<SeedMergerPixelLayer> nonSharedLayers( currentLayers );
	// then erase shared layers
	nonSharedLayers.erase( remove_if( nonSharedLayers.begin(), nonSharedLayers.end(), 
					  std::bind2nd( seedmergerIsEqualLayer(), sharedLayers ) ), 
			       nonSharedLayers.end() );

	// the remaining are those to be merged
	std::pair<SeedMergerPixelLayer, SeedMergerPixelLayer> layersToBeMerged( nonSharedLayers.at( 0 ), nonSharedLayers.at( 1 )  );

	// loop all possible triplet pairs (which come as input)
  	for( std::list<SeedingHitSet>::iterator triplet1 = tripletList.begin(); triplet1 != tripletList.end(); ++triplet1 ) {
  	  for( std::list<SeedingHitSet>::iterator triplet2 = ++std::list<SeedingHitSet>::iterator( triplet1 );
  	       triplet2 != tripletList.end(); ++triplet2 ) {

	    // do both triplets have shared hits on these two layers?
	    std::pair<TransientTrackingRecHit::ConstRecHitPointer,TransientTrackingRecHit::ConstRecHitPointer> sharedHits;
	    if( isTripletsShareHitsOnLayers( *triplet1, *triplet2, sharedLayers, sharedHits ) ) {

	      // are the remaining hits on different layers?
	      std::pair<TransientTrackingRecHit::ConstRecHitPointer,TransientTrackingRecHit::ConstRecHitPointer> nonSharedHits;
	      if( isMergeableHitsInTriplets( *triplet1, *triplet2, layersToBeMerged, nonSharedHits ) ) {

		// create an intermediate vector with all hits
		std::vector<TransientTrackingRecHit::ConstRecHitPointer> unsortedHits;
		unsortedHits.push_back( sharedHits.first );
		unsortedHits.push_back( sharedHits.second );
		unsortedHits.push_back( nonSharedHits.first );
		unsortedHits.push_back( nonSharedHits.second );

		// sort them by increasing radius
		std::sort( unsortedHits.begin(), unsortedHits.end(), boost::bind( &QuadrupletSeedMerger::isGreaterHit, *this, _1, _2 ) );

		// and create the quadruplet
		SeedingHitSet quadruplet;
		for( unsigned int aHit = 0; aHit < unsortedHits.size(); ++aHit ) quadruplet.add( unsortedHits[aHit] );

// 		std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
// 		std::cout << "~~~~ TRIPLET1:" << std::endl;
// 		printNtuplet( *triplet1 );
// 		std::cout << "~~~~ TRIPLET2:" << std::endl;
// 		printNtuplet( *triplet2 );
// 		std::cout << "~~~~ QUADPLET:" << std::endl;
// 		printNtuplet( quadruplet );
// 	        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

		// add to result
		if( isValidQuadruplet( quadruplet, currentLayers ) ) {
		  // insert this quadruplet
		  theResult.push_back( quadruplet );
		  // remove both triplets from the list,
		  // needs this 4-permutation since we're in a double loop
		  tripletList.erase( triplet2 );
		  triplet1 = tripletList.erase( triplet1 );
                  triplet2 = ( triplet1 == tripletList.end() ) ? --triplet1 : triplet1;
		}

	      } // isMergeableHitsInTriplets


 	    } // isTripletsShareHitsOnLayers

	    
	  } // triplet double loop
	}

	
      } // seeding layers double loop
    }
  
  } // seeding layer sets

  
  // add the remaining triplets
  if( isAddRemainingTriplets_ ) {
    for( std::list<SeedingHitSet>::const_iterator it = tripletList.begin(); it != tripletList.end(); ++it ) {
      theResult.push_back( *it );
    }
  }

  // some stats
  std::cout << " [QuadrupletSeedMerger::mergeTriplets] -- Created: " << theResult.size()
	    << " quadruplets from: " << nInputTriplets << " input triplets (" << tripletList.size()
	    << " remaining ";
  std::cout << (isAddRemainingTriplets_?"added":"dropped");
  std::cout << ")." << std::endl;

  return theResult;

}



///
/// merge triplets into quadruplets
/// INPUT: TrajectorySeedCollection
/// OUTPUT: TrajectorySeedCollection
///
/// this is a wrapper for:
/// vector<SeedingHitSet> mergeTriplets( const OrderedSeedingHits& )
/// for use in RecoTracker/TkSeedGenerator/plugins/SeedGeneratorFromRegionHitsEDProducer.cc
/// (iterative tracking)
///
const TrajectorySeedCollection QuadrupletSeedMerger::mergeTriplets( const TrajectorySeedCollection& seedCollection,
								    const TrackingRegion& region,
								    const edm::EventSetup& es ) {

  // ttrh builder for HitSet -> TrajectorySeed conversion;
  // require this to be correctly configured, otherwise -> exception
  es.get<TransientRecHitRecord>().get( theTTRHBuilderLabel_, theTTRHBuilder_ );

  // output collection
  TrajectorySeedCollection theResult;

  // loop to see if we have triplets ONLY
  // if not, copy input -> output and return
  bool isAllTriplets = true;
  for( TrajectorySeedCollection::const_iterator aTrajectorySeed = seedCollection.begin();
       aTrajectorySeed < seedCollection.end(); ++aTrajectorySeed ) {
    if( 3 != aTrajectorySeed->nHits() ) isAllTriplets = false;
  }

  if( !isAllTriplets && isMergeTriplets_ )
    std::cout << " [QuadrupletSeedMerger::mergeTriplets] (in RECO) -- bailing out since non-triplets in input." << std::endl;

  if( !isAllTriplets || !isMergeTriplets_ ) {
    for( TrajectorySeedCollection::const_iterator aTrajectorySeed = seedCollection.begin();
       aTrajectorySeed < seedCollection.end(); ++aTrajectorySeed ) {
      theResult.push_back( *aTrajectorySeed );
    }

    return theResult;
  }


  // all this fiddling here is now about converting
  // TrajectorySeedCollection <-> OrderedSeedingHits

  // create OrderedSeedingHits first;
  OrderedHitTriplets inputTriplets;

  // loop seeds
  for( TrajectorySeedCollection::const_iterator aTrajectorySeed = seedCollection.begin();
       aTrajectorySeed < seedCollection.end(); ++aTrajectorySeed ) {

    std::vector<TransientTrackingRecHit::RecHitPointer> recHitPointers;

    // loop RecHits
    const TrajectorySeed::range theHitsRange = aTrajectorySeed->recHits();
    for( edm::OwnVector<TrackingRecHit>::const_iterator aHit = theHitsRange.first;
	 aHit < theHitsRange.second; ++aHit ) {
      
      // this is a collection of: ReferenceCountingPointer< TransientTrackingRecHit> 
      recHitPointers.push_back( theTTRHBuilder_->build( &(*aHit) ) );

    }
    
    // add to input collection
    inputTriplets.push_back( OrderedHitTriplet( recHitPointers.at( 0 ), recHitPointers.at( 1 ), recHitPointers.at( 2 ) ) );

  }

  // do the real merging..
  const std::vector<SeedingHitSet> quadrupletHitSets = mergeTriplets( inputTriplets, es );
  
  // convert back to TrajectorySeedCollection
  SeedFromConsecutiveHitsCreator seedCreator;
  for( std::vector<SeedingHitSet>::const_iterator aHitSetIt = quadrupletHitSets.begin();
       aHitSetIt < quadrupletHitSets.end(); ++aHitSetIt ) {
    
    // add trajectory seed to result collection
    seedCreator.trajectorySeed( theResult, *aHitSetIt, region, es );
    
  }

  return theResult;
  
}



///
///
///
bool QuadrupletSeedMerger::isEqual( const TrackingRecHit* hit1, const TrackingRecHit* hit2 ) const {

  const double epsilon = 0.00001;

  const GeomDet* geomDet1 = theTrackerGeometry_->idToDet( hit1->geographicalId() );
  GlobalPoint const& gp1 = geomDet1->toGlobal( hit1->localPosition() );
  const GeomDet* geomDet2 = theTrackerGeometry_->idToDet( hit2->geographicalId() );
  GlobalPoint const& gp2 = geomDet2->toGlobal( hit2->localPosition() );

  if( ( fabs( gp1.x() - gp2.x() ) < epsilon ) &&
      ( fabs( gp1.y() - gp2.y() ) < epsilon ) &&
      ( fabs( gp1.z() - gp2.z() ) < epsilon ) ) {
    return true;
  }

  return false;

}



///
///
///
void QuadrupletSeedMerger::printHit( const TransientTrackingRecHit::ConstRecHitPointer& aRecHitPointer ) const {

  printHit( aRecHitPointer->hit() );

}



///
///
///
void QuadrupletSeedMerger::printHit( const TrackingRecHit* aHit ) const {

  const GeomDet* geomDet = theTrackerGeometry_->idToDet( aHit->geographicalId() );
  const double r = geomDet->surface().position().perp();
  const double x = geomDet->toGlobal( aHit->localPosition() ).x();
  const double y = geomDet->toGlobal( aHit->localPosition() ).y();
  const double z = geomDet->toGlobal( aHit->localPosition() ).z();
  std::cout << "<RecHit> x: " << x << " y: " << y << " z: " << z << " r: " << r << std::endl;

}

///
///
///
void QuadrupletSeedMerger::printNtuplet( const SeedingHitSet& aNtuplet ) const {

  std::cout << "DUMPING NTUPLET OF SIZE:";
  std::cout << aNtuplet.size() << std::endl;

  for( unsigned int aHit = 0; aHit < aNtuplet.size(); ++aHit ) {

    const TrackingRecHit* theHit = aNtuplet[aHit]->hit();
    const GeomDet* geomDet = theTrackerGeometry_->idToDet( theHit->geographicalId() );
    const double x = geomDet->toGlobal( theHit->localPosition() ).x();
    const double y = geomDet->toGlobal( theHit->localPosition() ).y();
    const double z = geomDet->toGlobal( theHit->localPosition() ).z();
    const double r = sqrt( x*x + y*y );

    unsigned int layer;
    std::string detName;
    if( PixelSubdetector::PixelBarrel == theHit->geographicalId().subdetId() ) {
      detName = "BPIX ";
      layer = PixelBarrelName::PixelBarrelName( aNtuplet[aHit]->hit()->geographicalId() ).layerName();
    }
    else {
      detName = "FPIX";
      if( z > 0 ) detName += "+";
      else detName += "-";
      layer = PixelEndcapName::PixelEndcapName( theHit->geographicalId() ).diskName();
    }

    std::cout << "<NtupletHit> D: " << detName << " L: " << layer << " x: " << x << " y: " << y << " z: " << z << " r: " << r << std::endl;
    
}

  std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

}



///
///
///
void QuadrupletSeedMerger::setTTRHBuilderLabel( std::string label ) {

  theTTRHBuilderLabel_ = label;
  
}



///
///
///
void QuadrupletSeedMerger::setLayerListName( std::string layerListName ) {
  layerListName_ = layerListName;
}



///
///
///
void QuadrupletSeedMerger::setMergeTriplets( bool isMergeTriplets ) {
  isMergeTriplets_ = isMergeTriplets;
}



///
///
///
void QuadrupletSeedMerger::setAddRemainingTriplets( bool isAddTriplets ) {
  isAddRemainingTriplets_ = isAddTriplets;
}



///
/// check for validity of a (radius-) *sorted* quadruplet:
///  1. after sorting, hits must be on layers according to the 
///     order given in PixelSeedMergerQuadruplets (from cfg)
///
bool QuadrupletSeedMerger::isValidQuadruplet( const SeedingHitSet& quadruplet, const std::vector<SeedMergerPixelLayer>& layers ) const {

  const unsigned int quadrupletSize = quadruplet.size();

  // basic size test..
  if( quadrupletSize != layers.size() ) {
    std::cout << " [QuadrupletSeedMerger::isValidQuadruplet] ** WARNING: size mismatch: "
	      << quadrupletSize << "/" << layers.size() << std::endl;
    return false;
  }

  // go along layers and check if all (parallel) quadruplet hits match
  for( unsigned int index = 0; index < quadruplet.size(); ++index ) {
    if( ! layers[index].isContainsDetector( quadruplet[index]->hit()->geographicalId() ) ) {
      return false;
    }
  }

  return true;

}



///
/// check if both triplets share a hit
/// on either of the sharedLayers
/// and return both hits (unsorted)
///
bool QuadrupletSeedMerger::isTripletsShareHitsOnLayers( const SeedingHitSet& firstTriplet, const SeedingHitSet& secondTriplet, 
  const std::pair<SeedMergerPixelLayer, SeedMergerPixelLayer>& sharedLayers,
  std::pair<TransientTrackingRecHit::ConstRecHitPointer,TransientTrackingRecHit::ConstRecHitPointer>& hits ) const {

  std::pair<bool, bool> isSuccess1( false, false );
  std::pair<bool, bool> isSuccess2( false, false );
  std::pair<TransientTrackingRecHit::ConstRecHitPointer,TransientTrackingRecHit::ConstRecHitPointer> hitsTriplet1, hitsTriplet2;

  // check if firstTriplet and secondTriplet have hits on sharedLayers
  for( unsigned int index = 0; index < 3; ++index ) {

    { // first triplet
      if( ! firstTriplet[index]->isValid() ) return false; // catch invalid TTRH pointers (tbd: erase triplet)
       DetId const& thisDetId = firstTriplet[index]->hit()->geographicalId();
      
      
       if( ! isSuccess1.first ) { // first triplet on shared layer 1
	 if( sharedLayers.first.isContainsDetector( thisDetId ) ) {
	   isSuccess1.first = true;
	   hitsTriplet1.first = firstTriplet[index];
	 }
       }
       
       if ( ! isSuccess1.second ) { // first triplet on shared layer 2
	 if( sharedLayers.second.isContainsDetector( thisDetId ) ) {
	   isSuccess1.second = true;
	   hitsTriplet1.second = firstTriplet[index];
	 }
       }

    }

    { // second triplet

      if( ! secondTriplet[index]->isValid() ) { return false; } // catch invalid TTRH pointers (tbd: erase triplet)
      DetId const& thisDetId = secondTriplet[index]->hit()->geographicalId();

      if( ! isSuccess2.first ) { // second triplet on shared layer 1
	if( sharedLayers.first.isContainsDetector( thisDetId ) ) {
	  isSuccess2.first = true;
	  hitsTriplet2.first = secondTriplet[index];
	}
      }
      
      if( ! isSuccess2.second ) { // second triplet on shared layer 2
	if( sharedLayers.second.isContainsDetector( thisDetId ) ) {
	  isSuccess2.second = true;
	  hitsTriplet2.second = secondTriplet[index];
	}
      }
    }


  }

  // check if these hits are pairwise equal
  if( isSuccess1.first && isSuccess1.second && isSuccess2.first && isSuccess2.second ) {
    
    if( isEqual( hitsTriplet1.first->hit(),  hitsTriplet2.first->hit()  ) &&
	isEqual( hitsTriplet1.second->hit(), hitsTriplet2.second->hit() )    ) {

      // copy to output, take triplet1 since they're equal anyway
      hits.first  = hitsTriplet1.first;
      hits.second = hitsTriplet1.second;

      return true;
    }
    
  }
  
  // empty output, careful
  return false;
  
}



///
/// check if the triplets have hits on the nonSharedLayers
/// triplet1 on layer1 && triplet2 on layer2, or vice versa,
/// and return the hits on those layers (unsorted)
///
bool QuadrupletSeedMerger::isMergeableHitsInTriplets( const SeedingHitSet& firstTriplet, const SeedingHitSet& secondTriplet, 
   const std::pair<SeedMergerPixelLayer, SeedMergerPixelLayer>& nonSharedLayers,
   std::pair<TransientTrackingRecHit::ConstRecHitPointer,TransientTrackingRecHit::ConstRecHitPointer>& hits ) const {

  // check if firstTriplet and secondTriplet have hits on sharedLayers
  for( unsigned int index1 = 0; index1 < 3; ++index1 ) {
    
    { // first triplet on non-shared layer 1
      DetId const& aDetId = firstTriplet[index1]->hit()->geographicalId();
      if( nonSharedLayers.first.isContainsDetector( aDetId ) ) {
	
	// look for hit in other (second) triplet on other layer
	for( unsigned int index2 = 0; index2 < 3; ++index2 ) {
	  
	  DetId const& anotherDetId = secondTriplet[index2]->hit()->geographicalId();
	  if( nonSharedLayers.second.isContainsDetector( anotherDetId ) ) {
	    
	    // ok!
	    hits.first  = firstTriplet[index1];
	    hits.second = secondTriplet[index2];
	    return true;
	    
	  }
	}
      }
    }

    // and vice versa..

    { // second triplet on non-shared layer 1
      DetId const& aDetId = secondTriplet[index1]->hit()->geographicalId();
      if( nonSharedLayers.first.isContainsDetector( aDetId ) ) {
	
	// look for hit in other (second) triplet on other layer
	for( unsigned int index2 = 0; index2 < 3; ++index2 ) {
	  
	  DetId const& anotherDetId = firstTriplet[index2]->hit()->geographicalId();
	  if( nonSharedLayers.second.isContainsDetector( anotherDetId ) ) {
	    
	    // ok!
	    hits.first  = firstTriplet[index1];
	    hits.second = secondTriplet[index2];
	    return true;
	    
	  }
	}
      }
    }

  } // for( index1

  return false;
    
}



///
///
///
SeedMergerPixelLayer::SeedMergerPixelLayer( const std::string& name ) {
  
  if( ! isValidName( name ) ) { 
    std::cerr << " [SeedMergerPixelLayer::SeedMergerPixelLayer] ** ERROR: illegal name: \"" << name << "\"." << std::endl;
    isValid_ = false;
    return;
  }

  // bare name, can be done here
  name_ = name;

  // (output format -> see DataFormats/SiPixelDetId/interface/PixelSubdetector.h)
  if( std::string::npos != name_.find( "BPix" ) ) subdet_ = PixelSubdetector::PixelBarrel;
  else if( std::string::npos != name_.find( "FPix" ) ) subdet_ = PixelSubdetector::PixelEndcap;
  else {
    std::cerr << " [PixelLayerNameParser::subdetector] ** ERROR: something's wrong here.." << std::endl;
  }

  // layer
  layer_ = atoi( name_.substr( 4, 1 ).c_str() );

}



///
///
///
SeedMergerPixelLayer::Side SeedMergerPixelLayer::getSide( void ) const {

  SeedMergerPixelLayer::Side side;

  if( std::string::npos != name_.find( "BPix" ) ) return Undefined; // no side in barrel
  else if( std::string::npos != name_.find( "pos", 6 ) ) return Plus;
  else if( std::string::npos != name_.find( "neg", 6 ) ) return Minus;

  else {
    std::cerr << " [PixelLayerNameParser::side] ** ERROR: something's wrong here.." << std::endl;
    side = SideError;
  }

  return side;

}



///
/// check if we have a name string as expected
///
bool SeedMergerPixelLayer::isValidName( const std::string& name ) {

  const int layer = atoi( name.substr( 4, 1 ).c_str() );

  if( std::string::npos != name.find( "BPix" ) ) {
    if( layer > 0 && layer < 5 ) return true;
  }

  else if( std::string::npos != name.find( "FPix" ) ) {
    if( layer > 0 && layer < 4 ) {
      if( std::string::npos != name.find( "pos", 6 ) || std::string::npos != name.find( "neg", 6 ) ) return true;
    }

  }

  std::cerr << " [SeedMergerPixelLayer::isValidName] ** WARNING: invalid name: \"" << name << "\"." << std::endl;
  return false;

}



///
/// check if the layer or disk described by this object
/// is the one carrying the detector: detId
///
bool SeedMergerPixelLayer::isContainsDetector( const DetId& detId ) const {

  PixelSubdetector::SubDetector subdet = getSubdet();

  // same subdet?
  if( detId.subdetId() == subdet ) {
    
    // same barrel layer?
    if( PixelSubdetector::PixelBarrel == subdet ) {
      
      const int layer = getLayerNumber();
      if( PixelBarrelName::PixelBarrelName( detId ).layerName() == layer ) {
	return true;
      }
      
    }

    // same endcap disk?
    else if( PixelSubdetector::PixelEndcap == subdet ) {
      
      const int layer = getLayerNumber();
      if( PixelEndcapName::PixelEndcapName( detId ).diskName() == layer ) {
	
	// endcap: same side?
	PixelEndcapName::HalfCylinder cylinder = PixelEndcapName::PixelEndcapName( detId ).halfCylinder();
	SeedMergerPixelLayer::Side side = getSide();

	if( ( ( PixelEndcapName::mO == cylinder || PixelEndcapName::mI == cylinder ) && SeedMergerPixelLayer::Minus == side ) ||
	    ( ( PixelEndcapName::pO == cylinder || PixelEndcapName::pI == cylinder ) && SeedMergerPixelLayer::Plus  == side )    ) {
	  
	  return true;

	}
	
      }
      
    }
    
  }
  
  return false;
  
}
