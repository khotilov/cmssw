#include "RecoTauTag/ImpactParameter/interface/ImpactParameterAlgorithm.h"

#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "RecoBTag/BTagTools/interface/SignedImpactParameter3D.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

ImpactParameterAlgorithm::ImpactParameterAlgorithm(){
        ip_min   = -9999;
        ip_max   = 9999;
        sip_min  = 0;
        use_sign = false;
	use3D    = false; 
}

ImpactParameterAlgorithm::ImpactParameterAlgorithm(const ParameterSet & parameters){
	ip_min       = parameters.getParameter<double>("TauImpactParameterMin");
	ip_max       = parameters.getParameter<double>("TauImpactParameterMax");
	sip_min	     = parameters.getParameter<double>("TauImpactParameterSignificanceMin");
	use_sign     = parameters.getParameter<bool>("UseTauImpactParameterSign");
	use3D        = parameters.getParameter<bool>("UseTau3DImpactParameter");

}

void ImpactParameterAlgorithm::setTransientTrackBuilder(const TransientTrackBuilder * builder) { 
	transientTrackBuilder = builder; 
}


pair<JetTag,TauImpactParameterInfo> ImpactParameterAlgorithm::tag(const IsolatedTauTagInfoRef & tauRef, const Vertex & pv) {

	if(transientTrackBuilder == 0){
	     cout << "Transient track builder is 0. abort!" << endl;
	     abort(); //FIXME: trow an exception here
	}

        TauImpactParameterInfo resultExtended;
	resultExtended.setIsolatedTauTag(tauRef);

	const Jet & jet = tauRef->jet();
	GlobalVector direction(jet.px(),jet.py(),jet.pz());

        const TrackRefVector& tracks = tauRef->selectedTracks();

	RefVector<TrackCollection>::const_iterator iTrack;
	for(iTrack = tracks.begin(); iTrack!= tracks.end(); iTrack++){

          const TransientTrack transientTrack = (transientTrackBuilder->build(&(**iTrack)));

          SignedTransverseImpactParameter stip;
	  Measurement1D ip = stip.apply(transientTrack,direction,pv).second;

	  SignedImpactParameter3D signed_ip3D;
	  Measurement1D ip3D = signed_ip3D.apply(transientTrack,direction,pv).second;
	  //cout << "check pv,ip3d,track z " << pv.z() << " " << ip3D.value() << " " << transientTrack->dz() << endl;
	  if(!use_sign){
	    Measurement1D tmp2D(fabs(ip.value()),ip.error());
	    ip = tmp2D;

	    Measurement1D tmp3D(fabs(ip3D.value()),ip3D.error());
            ip3D = tmp3D;
	  }

          reco::TauImpactParameterTrackData theData;

	  theData.transverseIp = ip;
          theData.ip3D = ip3D;
	  resultExtended.storeTrackData(*iTrack,theData);

	}

	double discriminator = resultExtended.discriminator(ip_min,ip_max,sip_min,use_sign,use3D);

	const JetTracksAssociationRef& jtaRef = tauRef->jtaRef();
	JetTag resultBase(discriminator,jtaRef);

	return pair<JetTag,TauImpactParameterInfo> (resultBase,resultExtended);
}
