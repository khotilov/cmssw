#include "SimCalorimetry/EcalSimAlgos/interface/EcalHitResponse.h" 
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVSimParameterMap.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloSimParameters.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVShape.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVHitCorrection.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVHitFilter.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVPECorrection.h"
#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h" 
#include<iostream>



EcalHitResponse::EcalHitResponse( const CaloVSimParameterMap* parameterMap ,
				  const CaloVShape*           shape         ) :
   m_parameterMap  ( parameterMap ) ,
   m_shape         ( shape        ) ,
   m_hitCorrection ( 0            ) ,
   m_PECorrection  ( 0            ) ,
   m_hitFilter     ( 0            ) ,
   m_geometry      ( 0            ) ,
   m_RandPoisson   ( 0            ) ,
   m_RandGauss     ( 0            ) ,
   m_minBunch      ( -10          ) ,
   m_maxBunch      (  10          ) ,
   m_phaseShift    ( 1            ) 
{
   edm::Service<edm::RandomNumberGenerator> rng ;
   if ( !rng.isAvailable() ) 
   {
      throw cms::Exception("Configuration")
	 << "EcalHitResponse requires the RandomNumberGeneratorService\n"
	 "which is not present in the configuration file.  You must add the service\n"
	 "in the configuration file or remove the modules that require it.";
   }
   m_RandPoisson = new CLHEP::RandPoissonQ( rng->getEngine() ) ;
   m_RandGauss   = new CLHEP::RandGaussQ(   rng->getEngine() ) ;
}

EcalHitResponse::~EcalHitResponse()
{
   delete m_RandPoisson ;
   delete m_RandGauss   ;
}

CLHEP::RandPoissonQ* 
EcalHitResponse::ranPois() const
{
   return m_RandPoisson ;
}

CLHEP::RandGaussQ* 
EcalHitResponse::ranGauss() const
{
   return m_RandGauss ;
}

const CaloSimParameters*
EcalHitResponse::params( const DetId& detId ) const
{
   assert( 0 != m_parameterMap ) ;
   return &m_parameterMap->simParameters( detId ) ;
}

const CaloVShape*
EcalHitResponse::shape() const
{
   assert( 0 != m_shape ) ;
   return m_shape ;
}

const CaloSubdetectorGeometry*
EcalHitResponse::geometry() const
{
   assert( 0 != m_geometry ) ;
   return m_geometry ;
}

void 
EcalHitResponse::setBunchRange( int minBunch , 
				int maxBunch  ) 
{
   m_minBunch = minBunch ;
   m_maxBunch = maxBunch ;
}

void 
EcalHitResponse::setGeometry( const CaloSubdetectorGeometry* geometry )
{
   m_geometry = geometry ;
}

void 
EcalHitResponse::setPhaseShift( double phaseShift )
{
   m_phaseShift = phaseShift ;
}

double
EcalHitResponse::phaseShift() const
{
   return m_phaseShift ;
}

void 
EcalHitResponse::setHitFilter( const CaloVHitFilter* filter)
{
   m_hitFilter = filter ;
}

void 
EcalHitResponse::setHitCorrection( const CaloVHitCorrection* hitCorrection)
{
   m_hitCorrection = hitCorrection ;
}

void 
EcalHitResponse::setPECorrection( const CaloVPECorrection* peCorrection )
{
   m_PECorrection = peCorrection ;
}

void 
EcalHitResponse::blankOutUsedSamples()  // blank out previously used elements
{
   const unsigned int size ( m_index.size() ) ;

   for( unsigned int i ( 0 ) ; i != size ; ++i )
   {
      vSamAll( m_index[i] )->setZero() ;
   }
   m_index.erase( m_index.begin() ,    // done and make ready to start over
		  m_index.end()    ) ;
}

void 
EcalHitResponse::run( MixCollection<PCaloHit>& hits ) 
{
   if( 0 != m_index.size() ) blankOutUsedSamples() ;

   for( MixCollection<PCaloHit>::MixItr hitItr ( hits.begin() ) ;
	hitItr != hits.end() ; ++hitItr )
   {
      const PCaloHit& hit ( *hitItr ) ;
      const int bunch ( hitItr.bunch() ) ;
      if( m_minBunch <= bunch  &&
	  m_maxBunch >= bunch  &&
	  !isnan( hit.time() ) &&
	  ( 0 == m_hitFilter ||
	    m_hitFilter->accepts( hit ) ) ) putAnalogSignal( hit ) ;
   }
}

void
EcalHitResponse::putAnalogSignal( const PCaloHit& inputHit )
{
   PCaloHit hit ( inputHit ) ;

   if( 0 != m_hitCorrection ) m_hitCorrection->correct( hit ) ;

   const DetId detId ( hit.id() ) ;

   const CaloSimParameters* parameters ( params( detId ) ) ;

   const double signal ( analogSignalAmplitude( hit ) ) ;

   const double jitter ( hit.time() - timeOfFlight( detId ) ) ;

   const double tzero = ( shape()->timeToRise()
			  + parameters->timePhase() 
			  - jitter 
			  - BUNCHSPACE*( parameters->binOfMaximum()
					 - m_phaseShift             ) ) ;
   double binTime ( tzero ) ;

   EcalSamples& result ( *findSignal( detId ) ) ;

   const unsigned int rsize ( result.size() ) ;

   for( unsigned int bin ( 0 ) ; bin != rsize ; ++bin )
   {
      result[ bin ] += (*shape())( binTime )*signal ;
      binTime += BUNCHSPACE ;
   }
}

EcalHitResponse::EcalSamples* 
EcalHitResponse::findSignal( const DetId& detId )
{
   EcalSamples* result ( vSamAll( CaloGenericDetId( detId ).denseIndex() ) ) ;
   if( result->zero() ) m_index.push_back( result - vSamAll(0) ) ;
   return result ;
}

double 
EcalHitResponse::analogSignalAmplitude( const PCaloHit& hit ) const
{
   const DetId& detId ( hit.id() ) ;

   const CaloSimParameters& parameters ( *params( detId ) ) ;

   // OK, the "energy" in the hit could be a real energy, deposited energy,
   // or pe count.  This factor converts to photoelectrons

   double npe ( hit.energy()*parameters.simHitToPhotoelectrons( detId ) ) ;

   // do we need to doPoisson statistics for the photoelectrons?
   if( parameters.doPhotostatistics() ) npe = ranPois()->fire( npe ) ;

   if( 0 != m_PECorrection ) npe = m_PECorrection->correctPE( detId, npe ) ;

   return npe ;
}

double 
EcalHitResponse::timeOfFlight( const DetId& detId ) const 
{
   const CaloCellGeometry* cellGeometry ( geometry()->getGeometry( detId ) ) ;
   assert( 0 != cellGeometry ) ;
   return cellGeometry->getPosition().mag()*cm/c_light ; // Units of c_light: mm/ns
}

void 
EcalHitResponse::add( const EcalSamples* pSam )
{
   EcalSamples& sam ( *findSignal( pSam->id() ) ) ;
   sam += (*pSam) ;
}

int 
EcalHitResponse::minBunch() const 
{
   return m_minBunch ; 
}

int 
EcalHitResponse::maxBunch() const 
{
   return m_maxBunch ; 
}

EcalHitResponse::VecInd& 
EcalHitResponse::index() 
{
   return m_index ; 
}

const CaloVHitFilter* 
EcalHitResponse::hitFilter() const 
{ 
   return m_hitFilter ; 
}
