#include "CaloExtractorByAssociator.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Utilities/Timing/interface/TimingReport.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace muonisolation;

CaloExtractorByAssociator::CaloExtractorByAssociator(const ParameterSet& par) :
  theUseRecHitsFlag(par.getParameter<bool>("UseRecHitsFlag")),
  theDepositLabel(par.getUntrackedParameter<string>("DepositLabel")),
  theDepositInstanceLabels(par.getParameter<std::vector<std::string> >("DepositInstanceLabels")),
  thePropagatorName(par.getParameter<std::string>("PropagatorName")),
  theThreshold_E(par.getParameter<double>("Threshold_E")),
  theThreshold_H(par.getParameter<double>("Threshold_H")),
  theThreshold_HO(par.getParameter<double>("Threshold_HO")),
  theDR_Veto_E(par.getParameter<double>("DR_Veto_E")),
  theDR_Veto_H(par.getParameter<double>("DR_Veto_H")),
  theDR_Veto_HO(par.getParameter<double>("DR_Veto_HO")),
  theDR_Max(par.getParameter<double>("DR_Max")),
  theNoise_EB(par.getParameter<double>("Noise_EB")),
  theNoise_EE(par.getParameter<double>("Noise_EE")),
  theNoise_HB(par.getParameter<double>("Noise_HB")),
  theNoise_HE(par.getParameter<double>("Noise_HE")),
  theNoise_HO(par.getParameter<double>("Noise_HO")),	
  theNoiseTow_EB(par.getParameter<double>("NoiseTow_EB")),
  theNoiseTow_EE(par.getParameter<double>("NoiseTow_EE")),
  theAssociator(0),
  thePropagator(0),
  thePrintTimeReport(par.getUntrackedParameter<bool>("PrintTimeReport"))
{
  theAssociatorParameters = new TrackAssociatorParameters(par.getParameter<edm::ParameterSet>("TrackAssociatorParameters"));
  theAssociator = new TrackDetectorAssociator();
}

CaloExtractorByAssociator::~CaloExtractorByAssociator(){
  if (thePrintTimeReport) TimingReport::current()->dump(std::cout);
  if (theAssociatorParameters) delete theAssociatorParameters;
  if (theAssociator) delete theAssociator;
  if (thePropagator) delete thePropagator;
}

void CaloExtractorByAssociator::fillVetos(const edm::Event& event, const edm::EventSetup& eventSetup, const TrackCollection& muons)
{
//   LogWarning("CaloExtractorByAssociator")
//     <<"fillVetos does nothing now: MuIsoDeposit provides enough functionality\n"
//     <<"to remove a deposit at/around given (eta, phi)";

}

MuIsoDeposit CaloExtractorByAssociator::deposit( const Event & event, const EventSetup& eventSetup, const Track & muon) const
{
  MuIsoDeposit::Direction muonDir(muon.eta(), muon.phi());
  MuIsoDeposit dep(theDepositLabel, muonDir );

//   LogWarning("CaloExtractorByAssociator")
//     <<"single deposit is not an option here\n"
//     <<"use ::deposits --> extract all and reweight as necessary";

  return dep;

}


//Make separate deposits: for ECAL, HCAL, HO
std::vector<MuIsoDeposit> CaloExtractorByAssociator::deposits( const Event & event, const EventSetup& eventSetup, const Track & muon) const
{
  if (thePropagator == 0){
    ESHandle<Propagator> prop;
    eventSetup.get<TrackingComponentsRecord>().get(thePropagatorName, prop);
    thePropagator = prop->clone();
    theAssociator->setPropagator(thePropagator);
  }
  if (theDepositInstanceLabels.size() != 3){
    LogError("MuonIsolation")<<"Configuration is inconsistent: Need 3 deposit instance labels";
  }
  if (! theDepositInstanceLabels[0].compare(0,1, std::string("e")) == 0
      || ! theDepositInstanceLabels[1].compare(0,1, std::string("h")) == 0
      || ! theDepositInstanceLabels[2].compare(0,2, std::string("ho")) == 0){
    LogWarning("MuonIsolation")<<"Deposit instance labels do not look like  (e*, h*, ho*):"
			       <<"proceed at your own risk. The extractor interprets lab0=from ecal; lab1=from hcal; lab2=from ho";
  }

  typedef MuIsoDeposit::Veto Veto;
  MuIsoDeposit::Direction muonDir(muon.eta(), muon.phi());
  
  MuIsoDeposit depEcal(theDepositInstanceLabels[0], muonDir);
  MuIsoDeposit depHcal(theDepositInstanceLabels[1], muonDir);
  MuIsoDeposit depHOcal(theDepositInstanceLabels[2], muonDir);

  edm::ESHandle<MagneticField> bField;
  eventSetup.get<IdealMagneticFieldRecord>().get(bField);


  reco::TransientTrack tMuon(muon, &*bField);
  FreeTrajectoryState iFTS = tMuon.initialFreeState();
  TrackDetMatchInfo mInfo = theAssociator->associate(event, eventSetup, iFTS, *theAssociatorParameters);

  depEcal.setVeto(Veto(Direction(mInfo.trkGlobPosAtEcal.eta(), mInfo.trkGlobPosAtEcal.phi()),
		       theDR_Veto_E));
  depHcal.setVeto(Veto(Direction(mInfo.trkGlobPosAtHcal.eta(), mInfo.trkGlobPosAtHcal.phi()),
		       theDR_Veto_H));
  depHOcal.setVeto(Veto(Direction(mInfo.trkGlobPosAtHO.eta(), mInfo.trkGlobPosAtHO.phi()),
			theDR_Veto_HO));

  if (theUseRecHitsFlag){
    edm::ESHandle<CaloGeometry> caloGeom;
    eventSetup.get<IdealGeometryRecord>().get(caloGeom);

    //Ecal
    std::vector<EcalRecHit>::const_iterator eHitCI = mInfo.ecalRecHits.begin();
    for (; eHitCI != mInfo.ecalRecHits.end(); ++eHitCI){
      GlobalPoint eHitPos = caloGeom->getPosition(eHitCI->detid());
      double deltar0 = deltaR(muon, eHitPos);
      double cosTheta = 1./cosh(eHitPos.eta());
      double energy = eHitCI->energy();
      double et = energy*cosTheta; 
      if (deltar0 > theDR_Max 
	  || ! (et > theThreshold_E && energy > 3*noiseRecHit(eHitCI->detid()))) continue;

      bool vetoHit = false;
      double deltar = deltaR(mInfo.trkGlobPosAtEcal, eHitPos);
      if (deltar < theDR_Veto_E ){
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Veto ECAL hit: Calo deltaR= " << deltar;
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Calo eta phi ethcal: " << eHitPos.eta() << " " << eHitPos.phi() << " " << et;
	vetoHit = true;
      }

      if (vetoHit ){
	depEcal.addMuonEnergy(et);
      } else {
	depEcal.addDeposit(Direction(eHitPos.eta(), eHitPos.phi()), et);      
      }
    }

    //Hcal
    std::vector<HBHERecHit>::const_iterator hHitCI = mInfo.hcalRecHits.begin();
    for (; hHitCI != mInfo.hcalRecHits.end(); ++hHitCI){
      GlobalPoint hHitPos = caloGeom->getPosition(hHitCI->detid());
      double deltar0 = deltaR(muon, hHitPos);
      double cosTheta = 1./cosh(hHitPos.eta());
      double energy = hHitCI->energy();
      double et = energy*cosTheta;
      if (deltar0 > theDR_Max 
	  || ! (et > theThreshold_H && energy > 3*noiseRecHit(hHitCI->detid()))) continue;

      bool vetoHit = false;
      double deltar = deltaR(mInfo.trkGlobPosAtHcal, hHitPos);
      if (deltar < theDR_Veto_H ){
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Veto HBHE hit: Calo deltaR= " << deltar;
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Calo eta phi ethcal: " << hHitPos.eta() << " " << hHitPos.phi() << " " << et;
	vetoHit = true;
      }

      if (vetoHit ){
	depHcal.addMuonEnergy(et);
      } else {
	depHcal.addDeposit(Direction(hHitPos.eta(), hHitPos.phi()), et);      
      }
    }

    //HOcal
    std::vector<HORecHit>::const_iterator hoHitCI = mInfo.hoRecHits.begin();
    for (; hoHitCI != mInfo.hoRecHits.end(); ++hoHitCI){
      GlobalPoint hoHitPos = caloGeom->getPosition(hoHitCI->detid());
      double deltar0 = deltaR(muon, hoHitPos);
      double cosTheta = 1./cosh(hoHitPos.eta());
      double energy = hoHitCI->energy();
      double et = energy*cosTheta;
      if (deltar0 > theDR_Max 
	  || ! (et > theThreshold_HO && energy > 3*noiseRecHit(hoHitCI->detid()))) continue;

      bool vetoHit = false;
      double deltar = deltaR(mInfo.trkGlobPosAtHO, hoHitPos);
      if (deltar < theDR_Veto_HO ){
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Veto HO hit: Calo deltaR= " << deltar;
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Calo eta phi ethcal: " << hoHitPos.eta() << " " << hoHitPos.phi() << " " << et;
	vetoHit = true;
      }

      if (vetoHit ){
	depHOcal.addMuonEnergy(et);
      } else {
	depHOcal.addDeposit(Direction(hoHitPos.eta(), hoHitPos.phi()), et);      	
      }
    }


  } else {
    //use calo towers    
    CaloTowerCollection::const_iterator calCI = mInfo.towers.begin();
    for (; calCI != mInfo.towers.end(); ++calCI){
      double deltar0 = deltaR(muon,*calCI);
      if (deltar0>theDR_Max) continue;
    
      //even more copy-pasting .. need to refactor
      double etecal = calCI->emEt();
      double eecal = calCI->emEnergy();
      bool doEcal = etecal>theThreshold_E && eecal>3*noiseEcal(*calCI);
      double ethcal = calCI->hadEt();
      double ehcal = calCI->hadEnergy();
      bool doHcal = ethcal>theThreshold_H && ehcal>3*noiseHcal(*calCI);
      double ethocal = calCI->outerEt();
      double ehocal = calCI->outerEnergy();
      bool doHOcal = ethocal>theThreshold_HO && ehocal>3*noiseHOcal(*calCI);
      if ((!doEcal) && (!doHcal) && (!doHcal)) continue;
    
      bool vetoTowerEcal = false;
      double deltarEcal = deltaR(mInfo.trkGlobPosAtEcal, *calCI);
      if (deltarEcal < theDR_Veto_E ){
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Veto ecal tower: Calo deltaR= " << deltarEcal;
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Calo eta phi ethcal: " << calCI->eta() << " " << calCI->phi() << " " << ethcal;
	vetoTowerEcal = true;
      }
      bool vetoTowerHcal = false;
      double deltarHcal = deltaR(mInfo.trkGlobPosAtHcal, *calCI);
      if (deltarHcal < theDR_Veto_H ){
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Veto hcal tower: Calo deltaR= " << deltarHcal;
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Calo eta phi ethcal: " << calCI->eta() << " " << calCI->phi() << " " << ethcal;
	vetoTowerHcal = true;
      }
      bool vetoTowerHOCal = false;
      double deltarHOcal = deltaR(mInfo.trkGlobPosAtHO, *calCI);
      if (deltarHOcal < theDR_Veto_HO ){
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Veto HO tower: Calo deltaR= " << deltarHOcal;
	LogDebug("RecoMuon|CaloExtractorByAssociator")
	  << " >>> Calo eta phi ethcal: " << calCI->eta() << " " << calCI->phi() << " " << ethcal;
	vetoTowerHOCal = true;
      }

      Direction towerDir(calCI->eta(), calCI->phi());
      if (doEcal){
	if (vetoTowerEcal) depEcal.addMuonEnergy(etecal);
	else depEcal.addDeposit(towerDir, etecal);
      }
      if (doHcal){
	if (vetoTowerHcal) depHcal.addMuonEnergy(ethcal);
	else depHcal.addDeposit(towerDir, ethcal);
      }
      if (doHOcal){
	if (vetoTowerHOCal) depHOcal.addMuonEnergy(ethocal);
	else depHOcal.addDeposit(towerDir, ethocal);
      }
    }
  }

  std::vector<MuIsoDeposit> resultDeps;    
  resultDeps.push_back(depEcal);
  resultDeps.push_back(depHcal);
  resultDeps.push_back(depHOcal);

  return resultDeps;

}

double CaloExtractorByAssociator::PhiInRange(const double& phi) {
      double phiout = phi;

      if( phiout > 2*M_PI || phiout < -2*M_PI) {
            phiout = fmod( phiout, 2*M_PI);
      }
      if (phiout <= -M_PI) phiout += 2*M_PI;
      else if (phiout >  M_PI) phiout -= 2*M_PI;

      return phiout;
}

template <class T, class U>
double CaloExtractorByAssociator::deltaR(const T& t, const U& u) {
      return sqrt(pow(t.eta()-u.eta(),2) +pow(PhiInRange(t.phi()-u.phi()),2));
}

double CaloExtractorByAssociator::noiseEcal(const CaloTower& tower) const {
      double noise = theNoiseTow_EB;
      double eta = tower.eta();
      if (fabs(eta)>1.479) noise = theNoiseTow_EE;
      return noise;
}

double CaloExtractorByAssociator::noiseHcal(const CaloTower& tower) const {
  double noise = fabs(tower.eta())> 1.479 ? theNoise_HE : theNoise_HB;      
  return noise;
}

double CaloExtractorByAssociator::noiseHOcal(const CaloTower& tower) const {
      double noise = theNoise_HO;
      return noise;
}


double CaloExtractorByAssociator::noiseRecHit(const DetId& detId) const {
  double  noise = 100;
  DetId::Detector det = detId.det();
  if (det == DetId::Ecal){
    EcalSubdetector subDet = (EcalSubdetector)(detId.subdetId());
    if (subDet == EcalBarrel){
      noise = theNoise_EB;
    } else if (subDet == EcalEndcap){
      noise = theNoise_EE;
    }
  } else if (det == DetId::Hcal){
    HcalSubdetector subDet = (HcalSubdetector)(detId.subdetId());
    if (subDet == HcalBarrel){
      noise = theNoise_HB;
    } else if (subDet == HcalEndcap){
      noise = theNoise_HE;      
    } else if (subDet == HcalOuter){
      noise = theNoise_HO;
    }
  }
  return noise;
}

