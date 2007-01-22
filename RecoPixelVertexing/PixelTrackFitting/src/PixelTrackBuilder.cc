#include "PixelTrackBuilder.h"


#include "Geometry/Surface/interface/LocalError.h"
#include "Geometry/Surface/interface/BoundPlane.h"
#include "Geometry/Vector/interface/GlobalPoint.h"

#include "TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/LocalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include <sstream>
using namespace std;
using namespace reco;


template <class T> T sqr( T t) {return t*t;}

reco::Track * PixelTrackBuilder::build(
      const Measurement1D & pt,
      const Measurement1D & phi, 
      const Measurement1D & cotTheta,
      const Measurement1D & tip,  
      const Measurement1D & zip,
      float chi2,
      int   charge,
      const std::vector<const TrackingRecHit* >& hits,
      const MagneticField * mf) const 
{

  LogDebug("PixelTrackBuilder::build");
  LogTrace("")<<"reconstructed TRIPLET kinematics:\n"<<print(pt,phi,cotTheta,tip,zip,chi2,charge);

  float sinTheta = 1/sqrt(1+sqr(cotTheta.value()));
  float cosTheta = cotTheta.value()*sinTheta;
  int tipSign = tip.value() > 0 ? 1 : -1;

  AlgebraicSymMatrix m(5,0);
  float invPtErr = 1./sqr(pt.value()) * pt.error();
  m[0][0] = sqr(sinTheta) * (
              sqr(invPtErr)
            + sqr(cotTheta.error()/pt.value()*cosTheta * sinTheta)
            );
  m[0][2] = sqr( cotTheta.error()) * cosTheta * sqr(sinTheta) / pt.value();
  m[1][1] = sqr( phi.error() );
  m[2][2] = sqr( cotTheta.error());
  m[3][3] = sqr( tip.error() );
  m[4][4] = sqr( zip.error() );
  LocalTrajectoryError error(m);

  LocalTrajectoryParameters lpar(
    LocalPoint(tipSign*tip.value(), -tipSign*zip.value(), 0),
    LocalVector(0., -tipSign*pt.value()*cotTheta.value(), pt.value()),
    charge);

  Surface::RotationType rot(
      sin(phi.value())*tipSign, -cos(phi.value())*tipSign,             0,
                     0,                 0,     -1*tipSign,
      cos(phi.value()),          sin(phi.value()),             0);
  BoundPlane * impPointPlane = new BoundPlane(GlobalPoint(0.,0.,0.), rot);

  TrajectoryStateOnSurface impactPointState( lpar , error, *impPointPlane, mf, 1.0);
  
//  checkState(impactPointState,mf);
  LogTrace("")<<"constructed TSOS :\n"<<print(impactPointState);

  int ndof = 2*hits.size()-5;
  GlobalPoint vv = impactPointState.globalPosition();
  math::XYZPoint  pos( vv.x(), vv.y(), vv.z() );
  GlobalVector pp = impactPointState.globalMomentum();
  math::XYZVector mom( pp.x(), pp.y(), pp.z() );

  reco::Track * track = new reco::Track( chi2, ndof, pos, mom, 
        impactPointState.charge(), impactPointState.curvilinearError());

  LogTrace("") <<"RECONSTRUCTED TRACK:\n"<< print(*track)<<std::endl;;

  return track;
}

std::string PixelTrackBuilder::print(
    const Measurement1D & pt,
    const Measurement1D & phi,
    const Measurement1D & cotTheta,
    const Measurement1D & tip,
    const Measurement1D & zip,
    float chi2,
    int   charge) const
{
    ostringstream str;
    str <<"\t pt: "  << pt.value() <<"+/-"<<pt.error()  
        <<"\t phi: " << phi.value() <<"+/-"<<phi.error()
        <<"\t cot: " << cotTheta.value() <<"+/-"<<cotTheta.error()
        <<"\t tip: " << tip.value() <<"+/-"<<tip.error()
        <<"\t zip: " << zip.value() <<"+/-"<<zip.error()
        <<"\t charge: " << charge;
    return str.str();
}

std::string PixelTrackBuilder::print(const reco::Track & track) const
{

  Measurement1D phi( track.phi(), track.phiError());

  float theta = track.theta();
  float cosTheta = cos(theta);
  float sinTheta = sin(theta);
  float errLambda2 = sqr( track.lambdaError() );
  Measurement1D cotTheta(cosTheta/sinTheta, sqrt(errLambda2)/sqr(sinTheta));

  float pt_v = track.pt();
  float errInvP2 = sqr(track.qoverpError());
  float covIPtTheta = track.covariance(TrackBase::i_qoverp, TrackBase::i_lambda);
  float errInvPt2 = (   errInvP2
                      + sqr(cosTheta/pt_v)*errLambda2
                      + 2*(cosTheta/pt_v)*covIPtTheta
                     ) / sqr(sinTheta);
  Measurement1D pt(pt_v, sqr(pt_v)*sqrt(errInvPt2));

  Measurement1D tip(track.d0(), track.d0Error());

  Measurement1D zip(track.dz(), track.dzError());

  return print(pt, phi, cotTheta, tip, zip, track.chi2(),  track.charge());
}

std::string PixelTrackBuilder::print(const TrajectoryStateOnSurface & state) const
{

  float pt_v = state.globalMomentum().perp();
  float phi_v = state.globalMomentum().phi();
  float theta_v = state.globalMomentum().theta();

  CurvilinearTrajectoryError curv = state.curvilinearError();
  float errPhi2 = curv.matrix()(3,3);
  float errLambda2 = curv.matrix()(2,2);
  float errInvP2 = curv.matrix()(1,1);
  float covIPtTheta = curv.matrix()(1,2);
  float cosTheta = cos(theta_v);
  float sinTheta = sin(theta_v);
  float errInvPt2 = (   errInvP2
                      + sqr(cosTheta/pt_v)*errLambda2
                      + 2*(cosTheta/pt_v)*covIPtTheta) / sqr(sinTheta);
  float errCotTheta = sqrt(errLambda2)/sqr(sinTheta) ;
  Measurement1D pt(pt_v, sqr(pt_v) * sqrt(errInvPt2));
  Measurement1D phi(phi_v, sqrt(errPhi2) );
  Measurement1D cotTheta(cosTheta/sinTheta, errCotTheta);

  float zip_v = state.globalPosition().z();
  float zip_e = sqrt( state.localError().matrix()(5,5));
  Measurement1D zip(zip_v, zip_e);

  float tip_v  = state.localPosition().x(); 
  int tip_sign = (state.localMomentum().y()*cotTheta.value() > 0) ? -1 : 1;
  float tip_e  = sqrt( state.localError().matrix()(4,4) );
  Measurement1D tip( tip_sign*tip_v, tip_e);

  return print(pt, phi, cotTheta, tip, zip, 0., state.charge());
}

void PixelTrackBuilder::checkState(const TrajectoryStateOnSurface & state, const MagneticField* mf) const
{
  LogTrace("")<<" *** PixelTrackBuilder::checkState: ";
  LogTrace("")<<"INPUT,  ROTATION" << endl<<state.surface().rotation();
  LogTrace("")<<"INPUT,  TSOS:"<<endl<<state;
  
  TransverseImpactPointExtrapolator tipe(mf);
  TrajectoryStateOnSurface test= tipe.extrapolate(state, GlobalPoint(0,0,0) );
  LogTrace("")<<"CHECK-1 ROTATION" << endl<<"\n"<<test.surface().rotation();
  LogTrace("")<<"CHECK-1 TSOS" << endl<<test;

  TSCPBuilderNoMaterial tscpBuilder;
  TrajectoryStateClosestToPoint tscp =
      tscpBuilder(*(state.freeState()), GlobalPoint(0,0,0) );
  FreeTrajectoryState fs = tscp.theState();
  LogTrace("")<<"CHECK-2 FTS: " << fs;
}
