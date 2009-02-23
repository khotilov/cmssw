#include "RecoVertex/GaussianSumVertexFit/interface/AdaptiveGsfVertexFitter.h"
#include "RecoVertex/VertexTools/interface/DummyVertexSmoother.h"
#include "RecoVertex/GaussianSumVertexFit/interface/GsfVertexSmoother.h"
#include "RecoVertex/GaussianSumVertexFit/interface/MultiPerigeeLTSFactory.h"
#include "RecoVertex/GaussianSumVertexFit/interface/GsfVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"



AdaptiveGsfVertexFitter::AdaptiveGsfVertexFitter(const edm::ParameterSet& pSet,
	const LinearizationPointFinder & linP )
{

  bool limitComponents_ = pSet.getParameter<bool>("limitComponents");
  bool useSmoothing = pSet.getParameter<bool>("smoothTracks");

  DeepCopyPointerByClone<GsfVertexMerger> theMerger;

  if (limitComponents_) {
    edm::ParameterSet mergerPSet = pSet.getParameter<edm::ParameterSet>("GsfMergerParameters");
    theMerger = new GsfVertexMerger(mergerPSet);
  }
  VertexSmoother<5> * sm;

  if (useSmoothing) sm  = new GsfVertexSmoother(limitComponents_, &*theMerger) ;
   else sm = new DummyVertexSmoother<5>();

  theFitter = new AdaptiveVertexFitter(
      GeometricAnnealing(),
      linP,
      GsfVertexUpdator(limitComponents_, &*theMerger),
      GsfVertexTrackCompatibilityEstimator(),
      *sm,
      MultiPerigeeLTSFactory() );
 
  /**
   *   Reads the configurable parameters.
   *   \param maxshift if the vertex moves further than this (in cm),
   *   then we re-iterate.
   *   \param maxlpshift if the vertex moves further than this,
   *   then we re-linearize the tracks.
   *   \param maxstep that's the maximum of iterations that we allow for.
   *   \param weightthreshold that's the minimum track weight
   *   for a track to be considered "significant".
   *   If fewer than two tracks are significant, an exception is thrown.
   */
  theFitter->setParameters ( pSet.getParameter<double>("maxshift"),
                    pSet.getParameter<double>("maxlpshift"),
                    pSet.getParameter<int>("maxstep"),
                    pSet.getParameter<double>("weightthreshold") );

}

AdaptiveGsfVertexFitter::AdaptiveGsfVertexFitter(const AdaptiveGsfVertexFitter & original)
{
  theFitter = original.theFitter->clone();
}

AdaptiveGsfVertexFitter::~AdaptiveGsfVertexFitter()
{
  delete theFitter;
}
