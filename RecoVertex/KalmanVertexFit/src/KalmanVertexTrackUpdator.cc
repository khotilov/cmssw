#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackUpdator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
//#include "Utilities/GenUtil/interface/ReferenceCountingPointer.h"
#include "RecoVertex/VertexPrimitives/interface/VertexTrack.h"
#include "RecoVertex/VertexPrimitives/interface/CachingVertex.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"


RefCountedVertexTrack KalmanVertexTrackUpdator::update
	(const CachingVertex & vertex , RefCountedVertexTrack track) const

{
  pair<RefCountedRefittedTrackState, AlgebraicMatrix> thePair = 
  	trackRefit(vertex.vertexState(), track->linearizedTrack());
  return theVTFactory.vertexTrack(track->linearizedTrack(),
  	vertex.vertexState(), thePair.first, thePair.second, 
	track->weight());
}

pair<RefCountedRefittedTrackState, AlgebraicMatrix> 
KalmanVertexTrackUpdator::trackRefit(const VertexState & vertex,
	 RefCountedLinearizedTrackState linTrackState) const

{
  //Vertex position 
  GlobalPoint vertexPosition = vertex.position();

  AlgebraicVector vertexCoord(3);
  vertexCoord[0] = vertexPosition.x();
  vertexCoord[1] = vertexPosition.y();
  vertexCoord[2] = vertexPosition.z();
  AlgebraicSymMatrix vertexErrorMatrix = vertex.error().matrix();

//track information
  AlgebraicMatrix a = linTrackState->positionJacobian();
  AlgebraicMatrix b = linTrackState->momentumJacobian();

  AlgebraicVector trackParameters = 
  	linTrackState->predictedStateParameters();

  AlgebraicSymMatrix trackParametersWeight = 
  	linTrackState->predictedStateWeight();

  AlgebraicSymMatrix s = trackParametersWeight.similarityT(b);
  
  int ifail;
  s.invert(ifail);
  if(ifail !=0) throw VertexException
  	("KalmanVertexTrackUpdator::S matrix inversion failed");
   
  AlgebraicVector newTrackMomentumP =  s * b.T() * trackParametersWeight * 
    (trackParameters - linTrackState->constantTerm() - a*vertexCoord);

  AlgebraicMatrix refittedPositionMomentumConvariance = 
    -vertexErrorMatrix * a.T() * trackParametersWeight * b * s;

  AlgebraicSymMatrix refittedMomentumConvariance = s +  
     vertex.weight().matrix().similarityT(refittedPositionMomentumConvariance);

  
  int matrixSize = 3+refittedMomentumConvariance.num_col();
  AlgebraicMatrix  covMatrix(matrixSize, matrixSize);
  covMatrix.sub(1, 4, refittedPositionMomentumConvariance);
  covMatrix.sub(4, 1, refittedPositionMomentumConvariance.T());
  covMatrix.sub(1, 1, vertexErrorMatrix);
  covMatrix.sub(4, 4, refittedMomentumConvariance);
  AlgebraicSymMatrix covSymMatrix;
  covSymMatrix.assign(covMatrix);

  RefCountedRefittedTrackState refittedTrackState = linTrackState->
	createRefittedTrackState(vertexPosition, newTrackMomentumP, covSymMatrix);

  return pair<RefCountedRefittedTrackState, AlgebraicMatrix>
  		(refittedTrackState, refittedPositionMomentumConvariance);
} 
