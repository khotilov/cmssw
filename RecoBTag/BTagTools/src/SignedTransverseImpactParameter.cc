#include "Geometry/Surface/interface/Line.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "Geometry/CommonDetAlgo/interface/Measurement1D.h"
#include "RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertError.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Matrix/Vector.h"
#include <string>


pair<bool,Measurement1D> SignedTransverseImpactParameter::apply(const Track & aRecTrack,
 const GlobalVector & direction, const  Vertex & vertex)
 const {

  TrajectoryStateOnSurface stateOnSurface = TransientTrack(aRecTrack).impactPointState();

  FreeTrajectoryState * FTS =  stateOnSurface.freeState(); //aRecTrack.stateAtFirstPoint().freeTrajectoryState());

  GlobalPoint vertexPosition(vertex.x(),vertex.y(),vertex.z());
 
  double theValue;
  double theError;
  TransverseImpactPointExtrapolator TIPE;
  TrajectoryStateOnSurface TSOS = TIPE.extrapolate(*FTS, vertexPosition);

  if(!TSOS.isValid()){
    theValue=0.;
    theError=0.;
  }
  else{
    
    GlobalPoint D0(TSOS.globalPosition());
    
    GlobalVector DD0(D0.x()-vertex.x(),D0.y()-vertex.y(),0.);
    GlobalVector JetDir(direction);
    double ps = DD0.dot(JetDir);
    theValue = DD0.mag()*(ps/abs(ps));

    //error calculation
    
    HepVector deriv(6);
    HepVector deriv_v(3);
    GlobalVector dd0 = DD0.unit();//check
    
    deriv_v[0] = - dd0.x();
    deriv_v[1] = - dd0.y();
    deriv_v[2] = - dd0.z();
    
    deriv[0] = dd0.x();
    deriv[1] = dd0.y();
    deriv[2] = dd0.z();
    deriv[3] =  0.;
    deriv[4] =  0.;
    deriv[5] =  0.;
    
    double E1 = (TSOS.cartesianError().matrix()).similarity(deriv);
    double E2 = RecoVertex::convertError(vertex.error()).matrix().similarity(deriv_v);
             // (aJet.vertex().positionError().matrix()).similarity(deriv_v);
    theError = sqrt(E1+E2);
  }//end if
  
  bool x = true;
  
  Measurement1D A(theValue, theError);
  return pair<bool,Measurement1D>(x,A);
}// end constructor declaration



pair<bool,Measurement1D> SignedTransverseImpactParameter::zImpactParameter ( const Track & aRecTrack, 
       const GlobalVector & direction, const  Vertex & vertex) const {
  
  TransverseImpactPointExtrapolator TIPE ;
  TrajectoryStateOnSurface TSOS  = TransientTrack(aRecTrack).impactPointState();

  if ( !TSOS.isValid() ) {
    cout << "====>>>> SignedTransverseImpactParameter::zImpactParameter : TSOS not valid" << endl ;
    return pair<bool,Measurement1D> (false,Measurement1D(0.0,0.0)) ;
  }

  GlobalPoint PV(vertex.x(),vertex.y(),vertex.z());
  TrajectoryStateOnSurface statePCA = TIPE.extrapolate( TSOS , PV ) ;

  if ( !statePCA.isValid() ) {
    cout << "====>>>> SignedTransverseImpactParameter::zImpactParameter : statePCA not valid" << endl ;
    return pair<bool,Measurement1D> (false,Measurement1D(0.0,0.0)) ;
  }

  GlobalPoint PCA = statePCA.globalPosition() ;

  // sign as in rphi
  GlobalVector PVPCA ( PCA.x()-PV.x() , PCA.y()-PV.y() , 0. );
  GlobalVector JetDir(direction);
  double sign = PVPCA.dot(JetDir);
  sign /= fabs(sign) ;

  // z IP
  double deltaZ = fabs(PCA.z()-PV.z()) * sign ;
  
  // error
  double errPvZ2 = RecoVertex::convertError(vertex.error()).czz() ;
  //CW  cout << "CW number or rows and columns : " << statePCA.cartesianError().matrix().num_row() << " , "
  //CW                                             << statePCA.cartesianError().matrix().num_col() << endl ;
  double errTrackZ2 =  statePCA.cartesianError().matrix()[2][2] ;
  double errZ = sqrt ( errPvZ2 + errTrackZ2 ) ;

  // CW alt. -> gives the same!!
  //  double errTZAlt = statePCA.localError().matrix()[4][4] ;
  //CW
  
  return pair<bool,Measurement1D> ( true , Measurement1D( deltaZ , errZ ) ) ;
}









