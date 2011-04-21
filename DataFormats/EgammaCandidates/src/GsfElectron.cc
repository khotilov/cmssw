#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

using namespace reco ;

GsfElectronCoreRef GsfElectron::core() const { return core_ ; }

void GsfElectron::init()
 {
  passCutBasedPreselection_ = false ;
  passMvaPreslection_ = false ;
  ambiguous_ = true ;
  fbrem_ = 0 ;
  class_ = UNKNOWN ;
 }

GsfElectron::GsfElectron()
 { init() ; }

GsfElectron::GsfElectron( const GsfElectronCoreRef & core )
 : core_(core)
 { init() ; }

GsfElectron::GsfElectron
 ( int charge, const ChargeInfo & chargeInfo,
   const GsfElectronCoreRef & core,
   const TrackClusterMatching & tcm, const TrackExtrapolations & te,
   const ClosestCtfTrack & ctfInfo,
   const FiducialFlags & ff, const ShowerShape & ss,
   const ConversionRejection & crv,
   float fbrem
 )
 : chargeInfo_(chargeInfo),
   core_(core),
   trackClusterMatching_(tcm), trackExtrapolations_(te),
   //closestCtfTrack_(ctfInfo),
   fiducialFlags_(ff), showerShape_(ss), conversionRejection_(crv)
 {
  init() ;
  fbrem_ = fbrem ;
  setCharge(charge) ;
  setVertex(math::XYZPoint(te.positionAtVtx.x(),te.positionAtVtx.y(),te.positionAtVtx.z())) ;
  setPdgId(-11*charge) ;
  /*if (ecalDrivenSeed())*/ corrections_.ecalEnergy = superCluster()->energy() ;
  assert(ctfInfo.ctfTrack==(GsfElectron::core()->ctfTrack())) ;
  assert(ctfInfo.shFracInnerHits==(GsfElectron::core()->ctfGsfOverlap())) ;
 }

GsfElectron::GsfElectron
 ( const GsfElectron & electron,
   const GsfElectronCoreRef & core )
 : RecoCandidate(electron),
   chargeInfo_(electron.chargeInfo_),
   core_(core),
   trackClusterMatching_(electron.trackClusterMatching_),
   trackExtrapolations_(electron.trackExtrapolations_),
   //closestCtfTrack_(electron.closestCtfTrack_),
   fiducialFlags_(electron.fiducialFlags_),
   showerShape_(electron.showerShape_),
   dr03_(electron.dr03_), dr04_(electron.dr04_),
   conversionRejection_(electron.conversionRejection_),
   passCutBasedPreselection_(electron.passCutBasedPreselection_),
   passMvaPreslection_(electron.passMvaPreslection_),
   ambiguous_(electron.ambiguous_),
   ambiguousGsfTracks_(electron.ambiguousGsfTracks_),
   //mva_(electron.mva_),
   fbrem_(electron.fbrem_),
   class_(electron.class_),
   corrections_(electron.corrections_)
 {
  assert(electron.core()->ctfTrack()==core->ctfTrack()) ;
  assert(electron.core()->ctfGsfOverlap()==core->ctfGsfOverlap()) ;
 }

GsfElectron::GsfElectron
 ( const GsfElectron & electron,
   const GsfElectronCoreRef & core,
   const CaloClusterPtr & electronCluster,
   const TrackRef & closestCtfTrack,
   const TrackBaseRef & conversionPartner,
   const GsfTrackRefVector & ambiguousTracks )
 : RecoCandidate(electron),
   chargeInfo_(electron.chargeInfo_),
   core_(core),
   trackClusterMatching_(electron.trackClusterMatching_),
   trackExtrapolations_(electron.trackExtrapolations_),
   //closestCtfTrack_(electron.closestCtfTrack_),
   fiducialFlags_(electron.fiducialFlags_),
   showerShape_(electron.showerShape_),
   dr03_(electron.dr03_), dr04_(electron.dr04_),
   conversionRejection_(electron.conversionRejection_),
   passCutBasedPreselection_(electron.passCutBasedPreselection_),
   passMvaPreslection_(electron.passMvaPreslection_),
   ambiguous_(electron.ambiguous_),
   ambiguousGsfTracks_(ambiguousTracks),
   //mva_(electron.mva_),
   fbrem_(electron.fbrem_),
   class_(electron.class_),
   corrections_(electron.corrections_)
 {
  trackClusterMatching_.electronCluster = electronCluster ;
  //closestCtfTrack_.ctfTrack = closestCtfTrack ;
  conversionRejection_.partner = conversionPartner ;
  assert(closestCtfTrack==core->ctfTrack()) ;
  assert(electron.core()->ctfGsfOverlap()==core->ctfGsfOverlap()) ;
  // TO BE DONE
  // Check that the new edm references are really
  // the clones of the former references, and therefore other attributes
  // stay valid :
  // * electron.core_ ~ core ?
  // * electron.trackClusterMatching_.electronCluster ~ electronCluster ?
  // * electron.closestCtfTrack_.ctfTrack ~ closestCtfTrack ?
  // * electron.ambiguousGsfTracks_ ~ ambiguousTracks ?
 }

bool GsfElectron::overlap( const Candidate & c ) const {
  const RecoCandidate * o = dynamic_cast<const RecoCandidate *>( & c );
  return ( o != 0 &&
	   ( checkOverlap( gsfTrack(), o->gsfTrack() ) ||
	     checkOverlap( superCluster(), o->superCluster() ) )
	   );
  //?? return false;
}

GsfElectron * GsfElectron::clone() const
 { return new GsfElectron(*this) ; }

GsfElectron * GsfElectron::clone
 (
  const GsfElectronCoreRef & core,
  const CaloClusterPtr & electronCluster,
  const TrackRef & closestCtfTrack,
  const TrackBaseRef & conversionPartner,
  const GsfTrackRefVector & ambiguousTracks
 ) const
 { return new GsfElectron(*this,core,electronCluster,closestCtfTrack,conversionPartner,ambiguousTracks) ; }

bool GsfElectron::ecalDriven() const
 {
  if (!passingCutBasedPreselection()&&!passingMvaPreselection())
   {
    edm::LogWarning("GsfElectronAlgo|UndefinedPreselectionInfo")
      <<"All preselection flags are false,"
      <<" either the data is too old or electrons were not preselected." ;
   }
  return (ecalDrivenSeed()&&passingCutBasedPreselection()) ;
 }

void GsfElectron::setCorrectedEcalEnergy( float newEnergy, float newEnergyError )
 {
  math::XYZTLorentzVectorD momentum = p4() ;
  momentum *= newEnergy/momentum.e() ;
  setP4(momentum) ;
  showerShape_.hcalDepth1OverEcal *= corrections_.ecalEnergy/newEnergy ;
  showerShape_.hcalDepth2OverEcal *= corrections_.ecalEnergy/newEnergy ;
  trackClusterMatching_.eSuperClusterOverP *= newEnergy/corrections_.ecalEnergy ;
  trackClusterMatching_.eSeedClusterOverP *= newEnergy/corrections_.ecalEnergy ;
  trackClusterMatching_.eEleClusterOverPout *= newEnergy/corrections_.ecalEnergy ;
  corrections_.ecalEnergy = newEnergy ;
  corrections_.ecalEnergyError = newEnergyError ;
  corrections_.isEcalEnergyCorrected = true ;
 }

void GsfElectron::setTrackMomentumError( float trackErr )
 { corrections_.trackMomentumError = trackErr ; }

void GsfElectron::setP4
 ( P4Kind kind, const reco::Candidate::LorentzVector & p4, float error, bool setCandidate )
 {
  switch(kind)
   {
    case P4_FROM_SUPER_CLUSTER:
      corrections_.fromSuperClusterP4 = p4 ;
      corrections_.fromSuperClusterP4Error = error ;
      break ;
    case P4_COMBINATION:
      corrections_.combinedP4 = p4 ;
      corrections_.combinedP4Error = error ;
      break ;
    case P4_PFLOW_COMBINATION:
      corrections_.pflowP4 = p4 ;
      corrections_.pflowP4Error = error ;
      break ;
    default:
      throw cms::Exception("GsfElectron")<<"unexpected p4 kind: "<<kind ;
   }
  if (setCandidate)
   {
    setP4(p4) ;
    corrections_.candidateP4Kind = kind ;
   }
 }

const Candidate::LorentzVector &  GsfElectron::p4( P4Kind kind ) const
 {
  switch(kind)
   {
    case P4_FROM_SUPER_CLUSTER: return corrections_.fromSuperClusterP4 ;
    case P4_COMBINATION: return corrections_.combinedP4 ;
    case P4_PFLOW_COMBINATION: return corrections_.pflowP4 ;
    default: throw cms::Exception("GsfElectron")<<"unexpected p4 kind: "<<kind ;
   }
 }


float GsfElectron::p4Error( P4Kind kind ) const
 {
  switch(kind)
   {
    case P4_FROM_SUPER_CLUSTER: return corrections_.fromSuperClusterP4Error ;
    case P4_COMBINATION: return corrections_.combinedP4Error ;
    case P4_PFLOW_COMBINATION: return corrections_.pflowP4Error ;
    default: throw cms::Exception("GsfElectron")<<"unexpected p4 kind: "<<kind ;
   }
 }

