#ifndef RecoBTag_SecondaryVertex_TrackSelector_h
#define RecoBTag_SecondaryVertex_TrackSelector_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"

namespace reco {

class TrackSelector {
    public:
	TrackSelector(const edm::ParameterSet &params);
	~TrackSelector() {}

	bool operator() (const reco::Track &track,
	                 const reco::TrackIPTagInfo::TrackIPData &ipData,
	                 const reco::Jet &jet) const;

    private:
	double	jetDeltaR;
	double	sip2dValMin;
	double	sip2dValMax;
	double	sip2dSigMin;
	double	sip2dSigMax;
	double	sip3dValMin;
	double	sip3dValMax;
	double	sip3dSigMin;
	double	sip3dSigMax;
};

} // namespace reco

#endif // RecoBTag_SecondaryVertex_TrackSelector_h
