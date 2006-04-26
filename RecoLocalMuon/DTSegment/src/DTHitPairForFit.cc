/** \file
 *
 * $Date: 2006/04/21 16:01:05 $
 * $Revision: 1.6 $
 * \author Stefano Lacaprara - INFN Legnaro <stefano.lacaprara@pd.infn.it>
 * \author Riccardo Bellan - INFN TO <riccardo.bellan@cern.ch>
 */

/* This Class Header */
#include "RecoLocalMuon/DTSegment/src/DTHitPairForFit.h"

/* Collaborating Class Header */
#include "FWCore/Utilities/interface/Exception.h"

/* C++ Headers */
#include <iostream>
using namespace std;

/* ====================================================================== */

/// Constructor
DTHitPairForFit::DTHitPairForFit(const DTRecHit1DPair& pair,
                                 const DTSuperLayer& sl,
                                 const edm::ESHandle<DTGeometry>& dtGeom) {

  theWireId = pair.wireId();
  theDigiTime = pair.digiTime();
  
  const DTLayer* layer = dtGeom->layer(theWireId.layerId());

  // transform the Local position in Layer-rf in a SL local position
  theLeftPos =
    sl.toLocal(layer->toGlobal(pair.componentRecHit(DTEnums::Left)->localPosition()));
  theRightPos =
    sl.toLocal(layer->toGlobal(pair.componentRecHit(DTEnums::Right)->localPosition()));

  // TODO how do I transform an error from local to global?
  theError = pair.componentRecHit(DTEnums::Left)->localPositionError();
  // theError =
  //   layer->surface().toLocal(sl.surface().toGlobal(pair.componentRecHit(DTEnums::Left)->localPositionError()));
  
}

/// Destructor
DTHitPairForFit::~DTHitPairForFit() {
}

/* Operations */ 
LocalPoint DTHitPairForFit::localPosition(DTEnums::DTCellSide s) const {
  if (s==DTEnums::Left) return theLeftPos;
  else if (s==DTEnums::Right) return theRightPos;
  else{ 
    throw cms::Exception("DTHitPairForFit")<<" localPosition called with undef LR code"<<endl;
    return LocalPoint();
  }
}

pair<bool,bool> 
DTHitPairForFit::isCompatible(const LocalPoint& posIni,
                              const LocalVector& dirIni) const {
  return std::pair<bool,bool>(isCompatible(posIni, dirIni, DTEnums::Left),
                              isCompatible(posIni, dirIni, DTEnums::Right));
}

bool DTHitPairForFit::isCompatible(const LocalPoint& posIni,
                                   const LocalVector& dirIni,
                                   DTEnums::DTCellSide code) const {
  // all is in SL frame
  LocalPoint pos= localPosition(code);
  LocalError err= localPositionError();

  const float errorScale=10.; // to be tuned!!

  LocalPoint segPosAtZ=
    posIni+dirIni*(pos.z()-posIni.z())/dirIni.z();

  // cout << "segPosAtZ     " << segPosAtZ << endl;
  // cout << "segPosInLayer " << pos<< endl;
  // cout << "errInLayer (" << err.xx() << "," << 
  //    err.xy() << "," << err.yy() << ")" << endl;

  float dx=pos.x()-segPosAtZ.x();
  // cout << "Dx " << dx << " vs " << sqrt(err.xx())*errorScale << endl;

  return fabs(dx)<sqrt(err.xx())*errorScale;

}

bool DTHitPairForFit::operator<(const DTHitPairForFit& hit) const {
  //SL if same layer use x() for strict ordering
  if (id()==hit.id()) 
    return (theLeftPos.x() < hit.localPosition(DTEnums::Left).x());
  return (id() < hit.id());
}

bool DTHitPairForFit::operator==(const DTHitPairForFit& hit) const {
  return  (id() == hit.id() && fabs(digiTime() - hit.digiTime()) < 0.1 );
}

ostream& operator<<(ostream& out, const DTHitPairForFit& hit) {
  out << hit.leftPos() << " " << hit.rightPos() ;
  return out;
}
