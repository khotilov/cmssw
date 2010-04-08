#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GapClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include <algorithm>

using namespace reco;

//
// constructors and destructor
//
PrimaryVertexProducerAlgorithm::PrimaryVertexProducerAlgorithm(const edm::ParameterSet& conf)
  // extract relevant parts of config for components
  : theConfig(conf), 
    theTrackFilter(conf.getParameter<edm::ParameterSet>("TkFilterParameters")), 
    theVertexSelector(VertexDistanceXY(), 
		      conf.getParameter<edm::ParameterSet>("PVSelParameters").getParameter<double>("maxDistanceToBeam"))
{
  edm::LogInfo("PVDebugInfo") 
    << "PVSelParameters::maxDistanceToBeam = " 
    << conf.getParameter<edm::ParameterSet>("PVSelParameters").getParameter<double>("maxDistanceToBeam") << "\n";


  fUseBeamConstraint = conf.getParameter<bool>("useBeamConstraint");
  fVerbose           = conf.getUntrackedParameter<bool>("verbose", false);
  fMinNdof           = conf.getParameter<double>("minNdof");
  fFailsafe          = true; //conf.getUntrackedParameter<bool>("failsafe",true);


  // select and configure the track clusterizer
  std::string clusteringAlgorithm=conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<std::string>("algorithm");
  if (clusteringAlgorithm=="gap"){
    theTrackClusterizer = new GapClusterizerInZ(conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkGapClusParameters"));
  }else if(clusteringAlgorithm=="DA"){
    theTrackClusterizer = new DAClusterizerInZ(conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  }else{
    throw VertexException("PrimaryVertexProducerAlgorithm: unknown clustering algorithm: " + clusteringAlgorithm);  
  }

  // select and configure the vertex fitter
  std::string algorithm = conf.getParameter<std::string>("algorithm");
  fapply_finder = false;
  if (algorithm == "TrimmedKalmanFinder") {
    fapply_finder = true;
    theFinder.setParameters(conf.getParameter<edm::ParameterSet>("VtxFinderParameters"));
  } else if (algorithm=="KalmanVertexFitter") {
    theFitter=new KalmanVertexFitter();
  } else if( algorithm=="AdaptiveVertexFitter") {
    theFitter=new AdaptiveVertexFitter();
  } else {
    throw VertexException("PrimaryVertexProducerAlgorithm: unknown algorithm: " + algorithm);  
  }

  edm::LogInfo("PVDebugInfo") 
    << "Using " << algorithm << "\n";
  edm::LogInfo("PVDebugInfo") 
    << "beam-constraint  " << fUseBeamConstraint << "\n"; 

  edm::LogInfo("PVDebugInfo") 
    << "PV producer algorithm initialization: done" << "\n";

}


PrimaryVertexProducerAlgorithm::~PrimaryVertexProducerAlgorithm() 
{
  if (theFitter) delete theFitter;
  if (theTrackClusterizer) delete theTrackClusterizer;
}


//
// member functions
//

// obsolete method, unfortunately required throgh inheritance from  VertexReconstructor
vector<TransientVertex> 
PrimaryVertexProducerAlgorithm::vertices(const vector<reco::TransientTrack> & tracks) const
{

   throw VertexException("PrimaryVertexProducerAlgorithm: cannot make a Primary Vertex without a beam spot constraint " );

  /*  std::cout<< "PrimaryVertexProducer::vertices> Obsolete function, using dummy beamspot " << std::endl;
    reco::BeamSpot dummyBeamSpot;
    dummyBeamSpot.dummy();
    return vertices(tracks,dummyBeamSpot); */
   return vector<TransientVertex>();
}


vector<TransientVertex> 
PrimaryVertexProducerAlgorithm::vertices(const vector<reco::TransientTrack> & tracks,
					 const reco::BeamSpot & beamSpot) const
{
  bool validBS = true;
  VertexState beamVertexState(beamSpot);
  if ( (beamVertexState.error().cxx() <= 0.) || 
       (beamVertexState.error().cyy() <= 0.) ||
       (beamVertexState.error().czz() <= 0.) ) {
    validBS = false;
    edm::LogError("UnusableBeamSpot") << "Beamspot with invalid errors "<<beamVertexState.error().matrix();
  }

  if ( fapply_finder) {
        return theFinder.vertices( tracks );
  }
  vector<TransientVertex> pvs;


  // select tracks
  vector<TransientTrack> seltks;

  if (validBS){
    for (vector<reco::TransientTrack>::const_iterator itk = tracks.begin();
	 itk != tracks.end(); itk++) {
      if (theTrackFilter(*itk)) seltks.push_back(*itk);
    }
  } else {
    seltks = tracks;
  }



  // clusterize tracks in Z
  vector< vector<reco::TransientTrack> > clusters =  theTrackClusterizer->clusterize(seltks);
  if (fVerbose){cout <<  " clustering returned  "<< clusters.size() << " clusters  from " << seltks.size() << " selected tracks" <<endl;}


  // look for primary vertices in each cluster
  vector<TransientVertex> pvCand;
  int nclu=0;
  for (vector< vector<reco::TransientTrack> >::const_iterator iclus
	 = clusters.begin(); iclus != clusters.end(); iclus++) {


    TransientVertex v;
    if( fUseBeamConstraint && validBS &&((*iclus).size()>1) ){
      if (fVerbose){cout <<  " constrained fit with "<< (*iclus).size() << " tracks"  <<endl;}
      v = theFitter->vertex(*iclus, beamSpot);
      if (v.isValid() && (v.degreesOfFreedom()>=fMinNdof)) pvCand.push_back(v);

      if (fVerbose){
	if (v.isValid()) cout << "x,y,z=" << v.position().x() <<" " << v.position().y() << " " <<  v.position().z() << endl;
	else cout <<"Invalid fitted vertex\n";
      }

    }else if((*iclus).size()>1){
      if (fVerbose){cout <<  " unconstrained fit with "<< (*iclus).size() << " tracks"  << endl;}

      v = theFitter->vertex(*iclus); 
      if (v.isValid() && (v.degreesOfFreedom()>=fMinNdof)) pvCand.push_back(v);

      if (fVerbose){
	if (v.isValid()) cout << "x,y,z=" << v.position().x() <<" " << v.position().y() << " " <<  v.position().z() << endl;
	else cout <<"Invalid fitted vertex\n";
      }

    }

    nclu++;

  }// end of cluster loop

  if(fVerbose){
    cout << "PrimaryVertexProducerAlgorithm::vertices  candidates =" << pvCand.size() << endl;
  }



  // select vertices compatible with beam
  int npv=0;
  for (vector<TransientVertex>::const_iterator ipv = pvCand.begin();
       ipv != pvCand.end(); ipv++) {
    if(fVerbose){
      cout << "PrimaryVertexProducerAlgorithm::vertices cand " << npv++ << " sel=" <<
	(validBS && theVertexSelector(*ipv,beamVertexState)) << "   z="  << ipv->position().z() << endl;
    }
    if (!validBS || theVertexSelector(*ipv,beamVertexState)) pvs.push_back(*ipv);
  }


  if(pvs.size()>0){

    // sort vertices by pt**2  vertex (aka signal vertex tagging)
    sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());

  }else{

    if (    fFailsafe 
	    && (seltks.size()>1) 
	    && ( (clusters.size()!=1)  ||  ( (clusters.size()==1) && (clusters.begin()->size()<seltks.size())) )
        )
      { 
	// if no vertex was found, try fitting all selected tracks, unless this has already been tried
	// in low/no pile-up situations with low multiplicity vertices, this can recover vertices lost in clustering
	// with very small zSep. Only makes sense when used with a robust fitter, like the AdaptiveVertexFitter
	TransientVertex v;
	if( fUseBeamConstraint && validBS ){
	  v = theFitter->vertex(seltks, beamSpot);
	}else{
	  v = theFitter->vertex(seltks); 
	}
	if (fVerbose){ cout << "PrimaryVertexProducerAlgorithm: failsafe vertex "
			    <<"  tracks="  << seltks.size()
			    <<"  valid()=" << v.isValid() << " ndof()=" << v.degreesOfFreedom() 
			    <<"  selected="<< theVertexSelector(v,beamVertexState) << endl; }
	if ( v.isValid() && (v.degreesOfFreedom()>=fMinNdof) && (theVertexSelector(v,beamVertexState)) )  pvs.push_back(v);
      }
    
  }


  return pvs;
  
}
