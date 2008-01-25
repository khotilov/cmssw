#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "RecoVertex/VertexTools/interface/AnnealingSchedule.h"
#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"
#include "RecoVertex/VertexTools/interface/VertexTrackFactory.h"
#include "RecoVertex/AdaptiveVertexFit/interface/KalmanChiSquare.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <algorithm>

// #define STORE_WEIGHTS
#ifdef STORE_WEIGHTS
#include <dataharvester/Writer.h>
#endif

using namespace std;

namespace {
  static const float initialError = 10000;
  int debug()
  {
    return 0;
  }

  #ifdef STORE_WEIGHTS
  map < RefCountedLinearizedTrackState, int > ids;
  int iter=0;

  int getId ( const RefCountedLinearizedTrackState & r )
  {
    static int ctr=1;
    if ( ids.count(r) == 0 )
    {
      ids[r]=ctr++;
    }
    return ids[r];
  }
  #endif
}

AdaptiveVertexFitter::AdaptiveVertexFitter(
      const AnnealingSchedule & ann,
      const LinearizationPointFinder & linP,
      const VertexUpdator<5> & updator,
      const VertexTrackCompatibilityEstimator<5> & crit,
      const VertexSmoother<5> & smoother,
      const AbstractLTSFactory<5> & ltsf ) :
    theNr(0),
    theLinP(linP.clone()), theUpdator( updator.clone()),
    theSmoother ( smoother.clone() ), theAssProbComputer( ann.clone() ),
    theComp ( crit.clone() ), theLinTrkFactory ( ltsf.clone() )
{
  setParameters();
}

void AdaptiveVertexFitter::setWeightThreshold ( float w )
{
  theWeightThreshold=w;
}

AdaptiveVertexFitter::AdaptiveVertexFitter
                        (const AdaptiveVertexFitter & o ) :
    theMaxShift ( o.theMaxShift ), theMaxLPShift ( o.theMaxLPShift ),
    theMaxStep ( o.theMaxStep ), theWeightThreshold ( o.theWeightThreshold ),
    theNr ( o.theNr ),
    theLinP ( o.theLinP->clone() ), theUpdator ( o.theUpdator->clone() ),
    theSmoother ( o.theSmoother->clone() ),
    theAssProbComputer ( o.theAssProbComputer->clone() ),
    theComp ( o.theComp->clone() ),
    theLinTrkFactory ( o.theLinTrkFactory->clone() )
{}

AdaptiveVertexFitter::~AdaptiveVertexFitter()
{
  delete theLinP;
  delete theUpdator;
  delete theSmoother;
  delete theAssProbComputer;
  delete theComp;
  delete theLinTrkFactory;
}

void AdaptiveVertexFitter::setParameters( double maxshift, double maxlpshift, 
                                          unsigned maxstep, double weightthreshold )
{
  theMaxShift = maxshift;
  theMaxLPShift = maxlpshift;
  theMaxStep = maxstep;
  theWeightThreshold=weightthreshold;
}


void AdaptiveVertexFitter::setParameters ( const edm::ParameterSet & s )
{
    setParameters ( s.getParameter<double>("maxshift"),
                    s.getParameter<double>("maxlpshift"),
                    s.getParameter<int>("maxstep"),
                    s.getParameter<double>("weightthreshold") );
}

CachingVertex<5>
AdaptiveVertexFitter::vertex(const vector<reco::TransientTrack> & tracks) const
{
  if ( tracks.size() < 2 )
  {
    throw VertexException("Supplied fewer than two tracks");
  };
  // Linearization Point
  GlobalPoint linP(0.,0.,0.);
    linP = theLinP->getLinearizationPoint(tracks);
  // Initial vertex seed, with a very large error matrix
  AlgebraicSymMatrix we(3,1);
  GlobalError error( we * initialError );
  VertexState seed (linP, error);
  vector<RefCountedVertexTrack> vtContainer = linearizeTracks(tracks, seed);
  return fit(vtContainer, seed, false);
}

CachingVertex<5>
AdaptiveVertexFitter::vertex(const vector<RefCountedVertexTrack> & tracks) const
{
  if ( tracks.size() < 2 )
  {
    throw VertexException( "Supplied fewer than two tracks" );
  };
  // Initial vertex seed, with a very small weight matrix
  GlobalPoint linP = tracks[0]->linearizedTrack()->linearizationPoint();
  AlgebraicSymMatrix we(3,1);
  GlobalError error( we * initialError );
  VertexState seed (linP, error);
  return fit(tracks, seed, false);
}

CachingVertex<5>
AdaptiveVertexFitter::vertex(const vector<RefCountedVertexTrack> & tracks, const reco::BeamSpot & spot ) const
{
  if ( tracks.size() < 2 )
  {
    throw VertexException( "Supplied fewer than two tracks" );
  };
  VertexState beamSpotState(spot);
  return fit(tracks, beamSpotState, true );
}



/** Fit vertex out of a set of reco::TransientTracks.
 *  Uses the specified linearization point.
 */
CachingVertex<5>
AdaptiveVertexFitter::vertex(const vector<reco::TransientTrack> & tracks,
                  const GlobalPoint& linPoint) const
{
  if ( tracks.size() < 2 )
  {
    throw VertexException( "Supplied fewer than two tracks" );
  };
  // Initial vertex seed, with a very large error matrix
  AlgebraicSymMatrix we(3,1);
  GlobalError error( we * initialError );
  VertexState seed (linPoint, error);
  vector<RefCountedVertexTrack> vtContainer = linearizeTracks(tracks, seed);
  return fit(vtContainer, seed, false);
}


/** Fit vertex out of a set of TransientTracks. 
 *  The specified BeamSpot will be used as priot, but NOT for the linearization.
 *  The specified LinearizationPointFinder will be used to find the linearization point.
 */
CachingVertex<5> 
AdaptiveVertexFitter::vertex(const vector<reco::TransientTrack> & tracks,
			       const reco::BeamSpot& beamSpot) const
{
  if ( tracks.size() < 1 )
  {
    throw VertexException( "Supplied no tracks" );
  };

  VertexState beamSpotState(beamSpot);

  /*
  cout << "[AdaptiveVertexFitter] beamspot vtx state=" << beamSpotState.position()
       << endl
       << beamSpotState.error().matrix()
       << endl;
       */
  vector<RefCountedVertexTrack> vtContainer;

  if (tracks.size() > 1) {
    // Linearization Point search if there are more than 1 track
    GlobalPoint linP(0.,0.,0.);
      linP = theLinP->getLinearizationPoint(tracks);
    AlgebraicSymMatrix we(3,1);
    // AlgebraicSymMatrix33 we;
    // we(0,0)=1; we(1,1)=1; we(2,2);
    GlobalError error(we*10000);
    VertexState lpState(linP, error);
    vtContainer = linearizeTracks(tracks, lpState);
  } else {
    // otherwise take the beamspot position.
    vtContainer = linearizeTracks(tracks, beamSpotState);
  }

  return fit(vtContainer, beamSpotState, true);
}


/** Fit vertex out of a set of reco::TransientTracks.
 *   Uses the position as both the linearization point AND as prior
 *   estimate of the vertex position. The error is used for the
 *   weight of the prior estimate.
 */
CachingVertex<5> AdaptiveVertexFitter::vertex(const vector<reco::TransientTrack> & tracks,
                  const GlobalPoint& priorPos,
                  const GlobalError& priorError) const

{
  if ( tracks.size() < 1 )
  {
    throw VertexException( "Supplied no tracks" );
  };
  VertexState seed (priorPos, priorError);
  vector<RefCountedVertexTrack> vtContainer = linearizeTracks(tracks, seed);
  return fit( vtContainer, seed, true );
}


/** Fit vertex out of a set of VertexTracks
 *   Uses the position and error for the prior estimate of the vertex.
 *   This position is not used to relinearize the tracks.
 */
CachingVertex<5> AdaptiveVertexFitter::vertex(
                const vector<RefCountedVertexTrack> & tracks,
                  const GlobalPoint& priorPos,
                  const GlobalError& priorError) const
{
  if ( tracks.size() < 1 )
  {
    throw VertexException( "Supplied no tracks" );
  };
  VertexState seed (priorPos, priorError);
  return fit(tracks, seed, true);
}


/**
 * Construct a container of VertexTrack from a set of reco::TransientTracks.
 * As this is the first iteration of the adaptive fit, the initial error
 * does not enter in the computation of the weights.
 * This is to avoid that all tracks get the same weight when
 * using a very large initial error matrix.
 */
vector<AdaptiveVertexFitter::RefCountedVertexTrack>
AdaptiveVertexFitter::linearizeTracks(const vector<reco::TransientTrack> & tracks,
                                      const VertexState & seed ) const
{
  GlobalPoint linP = seed.position();
  vector<RefCountedLinearizedTrackState> lTracks;
  for(vector<reco::TransientTrack>::const_iterator i = tracks.begin();
      i != tracks.end(); ++i )
  {
    try {
      RefCountedLinearizedTrackState lTrData
        = theLinTrkFactory->linearizedTrackState(linP, *i);
      lTracks.push_back(lTrData);
    } catch ( exception & e ) {
      cout << "[AdaptiveVertexFitter] Exception " << e.what() << " in ::linearizeTracks."
           << endl;
      cout << "[AdaptiveVertexFitter] Your future vertex has just lost a track."
           << endl;
    };
  }
  AlgebraicSymMatrix we(3,0);
  GlobalError nullError( we );
  VertexState initialSeed (linP, nullError);
  return weightTracks(lTracks, initialSeed);
}

/**
 * Construct new a container of VertexTrack with a new linearization point
 * and vertex seed, from an existing set of VertexTrack, from which only the
 * recTracks will be used.
 */
vector<AdaptiveVertexFitter::RefCountedVertexTrack>
AdaptiveVertexFitter::reLinearizeTracks(
                                const vector<RefCountedVertexTrack> & tracks,
                                const CachingVertex<5> & vertex ) const
{
  VertexState seed = vertex.vertexState();
  GlobalPoint linP = seed.position();
  vector<RefCountedLinearizedTrackState> lTracks;
  for(vector<RefCountedVertexTrack>::const_iterator i = tracks.begin();
    i != tracks.end(); i++)
  {
    try {
      RefCountedLinearizedTrackState lTrData
        = theLinTrkFactory->linearizedTrackState( linP, (**i).linearizedTrack()->track() );
      /*
      RefCountedLinearizedTrackState lTrData =
              (**i).linearizedTrack()->stateWithNewLinearizationPoint(linP);
              */
      lTracks.push_back(lTrData);
    } catch ( exception & e ) {
      cout << "[AdaptiveVertexFitter] Exception " << e.what()
           << " in ::relinearizeTracks." << endl;
      cout << "[AdaptiveVertexFitter] Will not relinearize this track." << endl;
      lTracks.push_back ( (**i).linearizedTrack() );
    };
  };
  return reWeightTracks(lTracks, vertex );
}

AdaptiveVertexFitter * AdaptiveVertexFitter::clone() const
{
  return new AdaptiveVertexFitter( * this );
}

vector<AdaptiveVertexFitter::RefCountedVertexTrack>
AdaptiveVertexFitter::reWeightTracks(
                    const vector<RefCountedLinearizedTrackState> & lTracks,
                    const CachingVertex<5> & vertex ) const
{
  VertexState seed = vertex.vertexState();
  theNr++;
  if (debug() & 4) cout << "Reweighting tracks... " << endl;
  // GlobalPoint pos = seed.position();

  vector<RefCountedVertexTrack> finalTracks;
  VertexTrackFactory<5> vTrackFactory;
  #ifdef STORE_WEIGHTS
  iter++;
  #endif
  for(vector<RefCountedLinearizedTrackState>::const_iterator i
        = lTracks.begin(); i != lTracks.end(); i++)
  {
    double chi2 = theComp->estimate ( vertex, *i );
    double weight = theAssProbComputer->weight(chi2);

    if ( weight > 1.0 )
    {
      cout << "[AdaptiveVertexFitter] Warning: weight " << weight << " > 1.0!"
           << endl;
      weight=1.0;
    };

    if ( weight < 0.0 )
    {
      cout << "[AdaptiveVertexFitter] Warning: weight " << weight << " < 0.0!"
           << endl;
      weight=0.0;
    };

    RefCountedVertexTrack vTrData
       = vTrackFactory.vertexTrack(*i, seed, theAssProbComputer->weight(chi2));

    #ifdef STORE_WEIGHTS
    map < string, dataharvester::MultiType > m;
    m["chi2"]=chi2;
    m["w"]=theAssProbComputer->weight(chi2);
    m["T"]=theAssProbComputer->currentTemp();
    m["n"]=iter;
    m["pos"]="reweight";
    m["id"]=getId ( *i );
    dataharvester::Writer::file("w.txt").save ( m );
    #endif

    finalTracks.push_back(vTrData);
  }
  return finalTracks;
}

vector<AdaptiveVertexFitter::RefCountedVertexTrack>
AdaptiveVertexFitter::weightTracks(
                    const vector<RefCountedLinearizedTrackState> & lTracks,
                    const VertexState & seed ) const
{
  theNr++;
  if (debug() & 4) cout << "Reweighting tracks... " << endl;
  GlobalPoint pos = seed.position();

  vector<RefCountedVertexTrack> finalTracks;
  VertexTrackFactory<5> vTrackFactory;
  KalmanChiSquare computer;
  #ifdef STORE_WEIGHTS
  iter++;
  #endif
  for(vector<RefCountedLinearizedTrackState>::const_iterator i
        = lTracks.begin(); i != lTracks.end(); i++)
  {
    double chi2 = computer.estimate ( pos, *i );
    double weight = theAssProbComputer->weight(chi2);

    if ( weight > 1.0 )
    {
      cout << "[AdaptiveVertexFitter] Warning: weight " << weight << " > 1.0!"
           << endl;
      weight=1.0;
    };

    if ( weight < 0.0 )
    {
      cout << "[AdaptiveVertexFitter] Warning: weight " << weight << " < 0.0!"
           << endl;
      weight=0.0;
    };

    RefCountedVertexTrack vTrData
       = vTrackFactory.vertexTrack(*i, seed, theAssProbComputer->weight(chi2));
    #ifdef STORE_WEIGHTS
    map < string, dataharvester::MultiType > m;
    m["chi2"]=chi2;
    m["w"]=theAssProbComputer->weight(chi2);
    m["T"]=theAssProbComputer->currentTemp();
    m["n"]=iter;
    m["id"]=getId ( *i );
    m["pos"]="weight";
    dataharvester::Writer::file("w.txt").save ( m );
    #endif
    finalTracks.push_back(vTrData);
  }
  return finalTracks;
}

/**
 * Construct new a container of VertexTrack with new weights
 * accounting for vertex error, from an existing set of VertexTracks.
 * From these the LinearizedTracks will be reused.
 */
vector<AdaptiveVertexFitter::RefCountedVertexTrack>
AdaptiveVertexFitter::reWeightTracks(
                            const vector<RefCountedVertexTrack> & tracks,
                            const CachingVertex<5> & seed) const
{
  vector<RefCountedLinearizedTrackState> lTracks;
  for(vector<RefCountedVertexTrack>::const_iterator i = tracks.begin();
    i != tracks.end(); i++)
  {
    lTracks.push_back((**i).linearizedTrack());
  }

  return reWeightTracks(lTracks, seed);
}


/*
 * The method where the vertex fit is actually done!
 */

CachingVertex<5>
AdaptiveVertexFitter::fit(const vector<RefCountedVertexTrack> & tracks,
                          const VertexState & priorSeed,
                          bool withPrior) const
{
  /*
  cout << "[AVF] debug: fitting wp=" << withPrior << " prior czz=" << priorSeed.error().czz() 
       << ", cxx=" << priorSeed.error().cxx() 
       << ", cyy=" << priorSeed.error().cyy() << endl;
       */
  theAssProbComputer->resetAnnealing();
  vector<RefCountedVertexTrack> initialTracks;
  GlobalPoint priorVertexPosition = priorSeed.position();
  GlobalError priorVertexError = priorSeed.error();

  CachingVertex<5> returnVertex( priorVertexPosition,priorVertexError,
                              initialTracks,0);
  if (withPrior)
  {
    returnVertex = CachingVertex<5>(priorVertexPosition,priorVertexError,
                    priorVertexPosition,priorVertexError,initialTracks,0);
  }

  vector<RefCountedVertexTrack> globalVTracks = tracks;

  // main loop through all the VTracks
  int lpStep = 0; int step = 0;

  CachingVertex<5> initialVertex = returnVertex;

  GlobalPoint newPosition = priorVertexPosition;
  GlobalPoint previousPosition = newPosition;

  int ns_trks=0; // number of significant tracks.
  // If we have only two significant tracks, we throw an
  // exception

  do {
    ns_trks=0;
    if (debug() & 8) cout << "lin point convergence step " << lpStep;
    if (debug() & 8) cout << " vtx pos convergence step " << step << endl;
    CachingVertex<5> fVertex = initialVertex;
    if ((previousPosition - newPosition).transverse() > theMaxLPShift)
    {
      if (debug() & 4) cout << "[AdaptiveVertexFitter] Relinearization." << endl;
      // relinearize and reweight.
      // (reLinearizeTracks also reweights tracks)
      globalVTracks = reLinearizeTracks( globalVTracks,
                             returnVertex );
      lpStep++;
    } else if (step) {
      // reweight, if it is not the first step
      globalVTracks = reWeightTracks( globalVTracks,
                                      returnVertex );
    }
    // cout << "-- new iteration" << endl;
    // update sequentially the vertex estimate
    for(vector<RefCountedVertexTrack>::const_iterator i
          = globalVTracks.begin(); i != globalVTracks.end(); i++)
    {
      try {
        // cout << "  -- fVtx pos=" << fVertex.position() << " cxx=" << fVertex.error().cxx() << " ndf=" << fVertex.degreesOfFreedom() << " add " << (**i).weight() << endl;
        fVertex = theUpdator->add( fVertex, *i );
        // cout << "  +- fVtx pos=" << fVertex.position() << " cxx=" << fVertex.error().cxx() << " ndf=" << fVertex.degreesOfFreedom() << endl;
        if ( (**i).weight() >= theWeightThreshold )
        {
          ns_trks++;
        };
      } catch ( exception & e ) {
        cout << "[AdaptiveVertexFitter] warning: updator throws " << e.what() << endl;
        cout << "[AdaptiveVertexFitter] (your vertex might just have lost one good track)"
             << endl;
      };
    }
    previousPosition = newPosition;
    newPosition = fVertex.position();
    returnVertex = fVertex;
    theAssProbComputer->anneal();
    step++;
    if ( step >= theMaxStep ) break;

  } while (
      // repeat as long as
      // - vertex moved too much or
      // - we're not yet annealed
         ( ((previousPosition - newPosition).mag() > theMaxShift) ||
           (!(theAssProbComputer->isAnnealed()) ) ) ) ;

  if (debug() ) cout << "[AdaptiveVertexFitter] debug: steps=" << step
     << " final temp=" << theAssProbComputer->currentTemp() 
     << " lpsteps=" << lpStep << endl;

  if ( theWeightThreshold > 0. &&  ns_trks < 2 && !withPrior ) 
  {
    ostringstream o;
    o << "fewer than two significant tracks (w>" << theWeightThreshold << ")";
    throw VertexException( o.str() );
  };

  #ifdef STORE_WEIGHTS
  map < string, dataharvester::MultiType > m;
  m["chi2"]=chi2;
  m["w"]=theAssProbComputer->weight(chi2);
  m["T"]=theAssProbComputer->currentTemp();
  m["n"]=iter;
  m["id"]=getId ( *i );
  m["pos"]="final";
  dataharvester::Writer::file("w.txt").save ( m );
  #endif
  return theSmoother->smooth( returnVertex );
}
