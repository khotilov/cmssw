#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/InputSort.h"
#include "RecoVertex/LinearizationPointFinders/interface/DefaultLinearizationPointFinder.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertexFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"



KinematicConstrainedVertexFitter::KinematicConstrainedVertexFitter()
{
 finder = new DefaultLinearizationPointFinder();
 vCons = new VertexKinematicConstraint();
 updator = new KinematicConstrainedVertexUpdator();
 tBuilder = new ConstrainedTreeBuilder;
 defaultParameters();
}

KinematicConstrainedVertexFitter:: KinematicConstrainedVertexFitter(const LinearizationPointFinder& fnd)
{

 finder = fnd.clone();
 vCons = new VertexKinematicConstraint();
 updator = new KinematicConstrainedVertexUpdator();
 tBuilder = new ConstrainedTreeBuilder;
 defaultParameters();
}

KinematicConstrainedVertexFitter::~KinematicConstrainedVertexFitter()
{
 delete finder;
 delete vCons;
 delete updator;
 delete tBuilder;
}

void KinematicConstrainedVertexFitter::setParameters(const edm::ParameterSet& pSet)
{
  theMaxDiff = pSet.getParameter<double>("maxDistance");
  theMaxStep = pSet.getParameter<int>("maxNbrOfIterations");;
  theMaxInitial = pSet.getParameter<double>("maxOfInitialValue");
}

void KinematicConstrainedVertexFitter::defaultParameters()
{
  theMaxDiff = 0.0001;
  theMaxStep = 1000;
  theMaxInitial = 9999.; //dummy value
}

RefCountedKinematicTree KinematicConstrainedVertexFitter::fit(vector<RefCountedKinematicParticle> part,
                                                             MultiTrackKinematicConstraint * cs)const
{
 if(part.size()<2) throw VertexException("KinematicConstrainedVertexFitter::input states are less than 2");

//sorting out the input particles
 InputSort iSort;
 pair<vector<RefCountedKinematicParticle>, vector<FreeTrajectoryState> > input = iSort.sort(part);
 const vector<RefCountedKinematicParticle> & particles  = input.first;
 const vector<FreeTrajectoryState> & fStates = input.second;

// linearization point:
 GlobalPoint linPoint  = finder->getLinearizationPoint(fStates);

//initial parameters:
 int vSize = particles.size();
 AlgebraicVector inPar(3 + 7*vSize,0);

//final parameters
 AlgebraicVector finPar(3 + 7*vSize,0);

//initial covariance
 AlgebraicMatrix inCov(3 + 7*vSize,3 + 7*vSize,0);

//making initial vector of parameters and initial particle-related covariance
 int nSt = 0;
 vector<KinematicState> inStates;
 for(vector<RefCountedKinematicParticle>::const_iterator i = particles.begin(); i!=particles.end(); i++)
 {
  KinematicState state = (*i)->stateAtPoint(linPoint);
  if (!state.isValid()) {
      LogDebug("KinematicConstrainedVertexFitter")
       << "State is invalid at point: "<<linPoint<<endl;
      return ReferenceCountingPointer<KinematicTree>(new KinematicTree());
  }
  AlgebraicVector prPar = asHepVector<7>(state.kinematicParameters().vector());
  for(int j = 1; j<8; j++){inPar(3 + 7*nSt + j) = prPar(j);}
  AlgebraicSymMatrix l_cov  = asHepMatrix<7>(state.kinematicParametersError().matrix());
  inCov.sub(4 + 7*nSt,4 + 7*nSt ,l_cov);
  inStates.push_back(state);
  ++nSt;
 }

//initial vertex error matrix components (huge error method)
//and vertex related initial vector components
 double in_er = 100.;
 inCov(1,1) = in_er;
 inCov(2,2) = in_er;
 inCov(3,3) = in_er;

 inPar(1) = linPoint.x();
 inPar(2) = linPoint.y();
 inPar(3) = linPoint.z();

//constraint equations value and number of iterations
 double eq;
 int nit = 0;

 vector<KinematicState> lStates = inStates;
 GlobalPoint lPoint  = linPoint;
 RefCountedKinematicVertex rVtx;
 AlgebraicMatrix refCCov;

//iterarions over the updator: each time updated parameters
//are taken as new linearization point
 do{
  eq = 0.;
  pair< pair< vector<KinematicState>, AlgebraicMatrix >,RefCountedKinematicVertex> lRes =
                                      updator->update(inPar,inCov,lStates,lPoint,cs);
  lStates = lRes.first.first;
  if (particles.size() != lStates.size()) {
    LogDebug("KinematicConstrainedVertexFitter")
	<< "updator failure\n";
    return ReferenceCountingPointer<KinematicTree>(new KinematicTree());
  }
  rVtx = lRes.second;
  lPoint = rVtx->position();

  AlgebraicVector vValue = vCons->value(lStates, lPoint);
  for(int i = 1; i<vValue.num_row();++i)
  {eq += abs(vValue(i));}
  if(cs !=0)
  {
   AlgebraicVector cVal = cs->value(lStates, lPoint);
   for(int i = 1; i<cVal.num_row();++i)
   {eq += abs(cVal(i));}
  }
  if (nit == 0) {
    if (eq>theMaxInitial) break;
  }
  refCCov = lRes.first.second;
  nit++;
 }while(nit<theMaxStep && eq>theMaxDiff);

// cout<<"number of relinearizations "<<nit<<endl;
// cout<<"value obtained: "<<eq<<endl;

  return  tBuilder->buildTree(particles, lStates, rVtx, refCCov);

}
