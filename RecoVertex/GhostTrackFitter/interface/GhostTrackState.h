#ifndef RecoBTag_GhostTrackState_h
#define RecoBTag_GhostTrackState_h

#include <utility>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

namespace reco {

class GhostTrackPrediction;

class GhostTrackState {
    public:
	typedef std::pair<GlobalPoint, GlobalError> Vertex;

	GhostTrackState(const TransientTrack &track) :
		track_(track), weight_(1.) {}

	const TransientTrack &track() const { return track_; }
	const TrajectoryStateOnSurface &tsos() const { return tsos_; }

	double lambda() const { return lambda_; }
	double lambdaError(const GhostTrackPrediction &pred,
	                   const GlobalError &pvError = GlobalError()) const;
	bool isValid() const { return tsos_.isValid(); }

	void reset() { tsos_ = TrajectoryStateOnSurface(); }
	bool linearize(const GhostTrackPrediction &pred,
	               bool initial = false, double lambda = 0.);
	bool linearize(const GhostTrackPrediction &pred, double lambda);

	double flightDistance(const GlobalPoint &point,
	                      const GlobalVector &dir) const;
	double axisDistance(const GlobalPoint &point,
	                    const GlobalVector &dir) const;
	double axisDistance(const GhostTrackPrediction &pred) const;

	Vertex vertexStateOnGhostTrack(const GhostTrackPrediction &pred,
	                               bool withRecoTrackError = true) const;
	Vertex vertexStateOnRecoTrack(const GhostTrackPrediction &pred,
	                              bool withGhostTrackError = true) const;

	double weight() const { return weight_; }
	void setWeight(double weight) { weight_ = weight; }

    private:
	TransientTrack			track_;
	TrajectoryStateOnSurface	tsos_;
	double				lambda_;
	double				weight_;
};

}

#endif // RecoBTag_GhostTrackState_h
