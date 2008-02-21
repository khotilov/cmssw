// Do not include .h from plugin directory, but locally:
#include "BzeroReferenceTrajectoryFactory.h"
#include "Alignment/ReferenceTrajectories/interface/BzeroReferenceTrajectory.h" 
#include "Alignment/ReferenceTrajectories/interface/TrajectoryFactoryPlugin.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

BzeroReferenceTrajectoryFactory::BzeroReferenceTrajectoryFactory( const edm::ParameterSet & config ) :
  TrajectoryFactoryBase( config )
{
  theMass = config.getParameter< double >( "ParticleMass" );
  theMomentumEstimate = config.getParameter< double >( "MomentumEstimate" );
}

 
BzeroReferenceTrajectoryFactory::~BzeroReferenceTrajectoryFactory( void ) {}


const BzeroReferenceTrajectoryFactory::ReferenceTrajectoryCollection
BzeroReferenceTrajectoryFactory::trajectories( const edm::EventSetup & setup,
					       const ConstTrajTrackPairCollection & tracks ) const
{
  ReferenceTrajectoryCollection trajectories;

  edm::ESHandle< MagneticField > magneticField;
  setup.get< IdealMagneticFieldRecord >().get( magneticField );

  ConstTrajTrackPairCollection::const_iterator itTracks = tracks.begin();

  while ( itTracks != tracks.end() )
  { 
    TrajectoryInput input = this->innermostStateAndRecHits( *itTracks );
    // set the flag for reversing the RecHits to false, since they are already in the correct order.
    trajectories.push_back( ReferenceTrajectoryPtr( new BzeroReferenceTrajectory( input.first, input.second, 
										  false, magneticField.product(),
										  materialEffects(), propagationDirection(),
										  theMass, theMomentumEstimate ) ) );
    ++itTracks;
  }

  return trajectories;
}


const BzeroReferenceTrajectoryFactory::ReferenceTrajectoryCollection
BzeroReferenceTrajectoryFactory::trajectories( const edm::EventSetup & setup,
					       const ConstTrajTrackPairCollection& tracks,
					       const ExternalPredictionCollection& external ) const
{
  ReferenceTrajectoryCollection trajectories;

  if ( tracks.size() != external.size() )
  {
    edm::LogInfo("ReferenceTrajectories") << "@SUB=BzeroReferenceTrajectoryFactory::trajectories"
					  << "Inconsistent input:\n"
					  << "\tnumber of tracks = " << tracks.size()
					  << "\tnumber of external predictions = " << external.size();

    return trajectories;
  }

  edm::ESHandle< MagneticField > magneticField;
  setup.get< IdealMagneticFieldRecord >().get( magneticField );

  ConstTrajTrackPairCollection::const_iterator itTracks = tracks.begin();
  ExternalPredictionCollection::const_iterator itExternal = external.begin();

  while ( itTracks != tracks.end() )
  {
    TrajectoryInput input = innermostStateAndRecHits( *itTracks  );

    if ( (*itExternal).isValid() && sameSurface( (*itExternal).surface(), input.first.surface() ) )
    {
      // set the flag for reversing the RecHits to false, since they are already in the correct order.
      ReferenceTrajectoryPtr refTraj( new BzeroReferenceTrajectory( *itExternal, input.second, false,
								    magneticField.product(), materialEffects(),
								    propagationDirection(), theMass, theMomentumEstimate ) );

      AlgebraicSymMatrix externalParamErrors( asHepMatrix<5>( (*itExternal).localError().matrix() ) );
      refTraj->setParameterErrors( externalParamErrors.sub( 2, 5 ) );
      trajectories.push_back( refTraj );
    }
    else
    {
      trajectories.push_back( ReferenceTrajectoryPtr( new BzeroReferenceTrajectory( input.first, input.second, 
										    false, magneticField.product(),
										    materialEffects(), propagationDirection(),
										    theMass, theMomentumEstimate ) ) );
    }

    ++itTracks;
    ++itExternal;
  }

  return trajectories;
}


DEFINE_EDM_PLUGIN( TrajectoryFactoryPlugin, BzeroReferenceTrajectoryFactory, "BzeroReferenceTrajectoryFactory" );
