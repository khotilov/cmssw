#include "RecoLocalCalo/CaloTowersCreator/interface/CaloTowersCreationAlgo.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoLocalCalo/CaloTowersCreator/interface/HcalMaterials.h"
#include "Math/Interpolator.h"
#include<iostream>

CaloTowersCreationAlgo::CaloTowersCreationAlgo()
 : theEBthreshold(-1000.),
   theEEthreshold(-1000.),
   theHcalThreshold(-1000.),
   theHBthreshold(-1000.),
   theHESthreshold(-1000.),
   theHEDthreshold(-1000.),
   theHOthreshold(-1000.),
   theHF1threshold(-1000.),
   theHF2threshold(-1000.),
   theEBGrid(std::vector<double>(5,10.)),
   theEBWeights(std::vector<double>(5,1.)),
   theEEGrid(std::vector<double>(5,10.)),
   theEEWeights(std::vector<double>(5,1.)),
   theHBGrid(std::vector<double>(5,10.)),
   theHBWeights(std::vector<double>(5,1.)),
   theHESGrid(std::vector<double>(5,10.)),
   theHESWeights(std::vector<double>(5,1.)),
   theHEDGrid(std::vector<double>(5,10.)),
   theHEDWeights(std::vector<double>(5,1.)),
   theHOGrid(std::vector<double>(5,10.)),
   theHOWeights(std::vector<double>(5,1.)),
   theHF1Grid(std::vector<double>(5,10.)),
   theHF1Weights(std::vector<double>(5,1.)),
   theHF2Grid(std::vector<double>(5,10.)),
   theHF2Weights(std::vector<double>(5,1.)),
   theEBweight(1.),
   theEEweight(1.),
   theHBweight(1.),
   theHESweight(1.),
   theHEDweight(1.),
   theHOweight(1.),
   theHF1weight(1.),
   theHF2weight(1.),
   theEcutTower(-1000.),
   theEBSumThreshold(-1000.),
   theEESumThreshold(-1000.),
   theHcalTopology(0),
   theGeometry(0),
   theTowerConstituentsMap(0),
   theHOIsUsed(true),
   // (for momentum reconstruction algorithm)
   theMomConstrMethod(0),
   theMomEmDepth(0.),
   theMomHadDepth(0.),
   theMomTotDepth(0.)
{
}

CaloTowersCreationAlgo::CaloTowersCreationAlgo(double EBthreshold, double EEthreshold, double HcalThreshold,
    double HBthreshold, double HESthreshold, double  HEDthreshold,
    double HOthreshold, double HF1threshold, double HF2threshold,
    double EBweight, double EEweight,
    double HBweight, double HESweight, double HEDweight,
    double HOweight, double HF1weight, double HF2weight,
    double EcutTower, double EBSumThreshold, double EESumThreshold,
    bool useHO,
    // (momentum reconstruction algorithm)
    int momConstrMethod,
    double momEmDepth,
    double momHadDepth,
    double momTotDepth)

 : theEBthreshold(EBthreshold),
   theEEthreshold(EEthreshold),
   theHcalThreshold(HcalThreshold),
   theHBthreshold(HBthreshold),
   theHESthreshold(HESthreshold),
   theHEDthreshold(HEDthreshold),
   theHOthreshold(HOthreshold),
   theHF1threshold(HF1threshold),
   theHF2threshold(HF2threshold),
   theEBGrid(std::vector<double>(5,10.)),
   theEBWeights(std::vector<double>(5,1.)),
   theEEGrid(std::vector<double>(5,10.)),
   theEEWeights(std::vector<double>(5,1.)),
   theHBGrid(std::vector<double>(5,10.)),
   theHBWeights(std::vector<double>(5,1.)),
   theHESGrid(std::vector<double>(5,10.)),
   theHESWeights(std::vector<double>(5,1.)),
   theHEDGrid(std::vector<double>(5,10.)),
   theHEDWeights(std::vector<double>(5,1.)),
   theHOGrid(std::vector<double>(5,10.)),
   theHOWeights(std::vector<double>(5,1.)),
   theHF1Grid(std::vector<double>(5,10.)),
   theHF1Weights(std::vector<double>(5,1.)),
   theHF2Grid(std::vector<double>(5,10.)),
   theHF2Weights(std::vector<double>(5,1.)),
   theEBweight(EBweight),
   theEEweight(EEweight),
   theHBweight(HBweight),
   theHESweight(HESweight),
   theHEDweight(HEDweight),
   theHOweight(HOweight),
   theHF1weight(HF1weight),
   theHF2weight(HF2weight),
   theEcutTower(EcutTower),
   theEBSumThreshold(EBSumThreshold),
   theEESumThreshold(EESumThreshold),
   theHOIsUsed(useHO),
   // (momentum reconstruction algorithm)
   theMomConstrMethod(momConstrMethod),
   theMomEmDepth(momEmDepth),
   theMomHadDepth(momHadDepth),
   theMomTotDepth(momTotDepth)
{
}

CaloTowersCreationAlgo::CaloTowersCreationAlgo(double EBthreshold, double EEthreshold, double HcalThreshold,
    double HBthreshold, double HESthreshold, double  HEDthreshold,
    double HOthreshold, double HF1threshold, double HF2threshold,
    std::vector<double> EBGrid, std::vector<double> EBWeights,
    std::vector<double> EEGrid, std::vector<double> EEWeights,
    std::vector<double> HBGrid, std::vector<double> HBWeights,
    std::vector<double> HESGrid, std::vector<double> HESWeights,
    std::vector<double> HEDGrid, std::vector<double> HEDWeights,
    std::vector<double> HOGrid, std::vector<double> HOWeights,
    std::vector<double> HF1Grid, std::vector<double> HF1Weights,
    std::vector<double> HF2Grid, std::vector<double> HF2Weights,
    double EBweight, double EEweight,
    double HBweight, double HESweight, double HEDweight,
    double HOweight, double HF1weight, double HF2weight,
    double EcutTower, double EBSumThreshold, double EESumThreshold,
    bool useHO,
    // (for the momentum construction algorithm)
    int momConstrMethod,
    double momEmDepth,
    double momHadDepth,
    double momTotDepth
    )

 : theEBthreshold(EBthreshold),
   theEEthreshold(EEthreshold),
   theHcalThreshold(HcalThreshold),
   theHBthreshold(HBthreshold),
   theHESthreshold(HESthreshold),
   theHEDthreshold(HEDthreshold),
   theHOthreshold(HOthreshold),
   theHF1threshold(HF1threshold),
   theHF2threshold(HF2threshold),
   theEBGrid(EBGrid),
   theEBWeights(EBWeights),
   theEEGrid(EEGrid),
   theEEWeights(EEWeights),
   theHBGrid(HBGrid),
   theHBWeights(HBWeights),
   theHESGrid(HESGrid),
   theHESWeights(HESWeights),
   theHEDGrid(HEDGrid),
   theHEDWeights(HEDWeights),
   theHOGrid(HOGrid),
   theHOWeights(HOWeights),
   theHF1Grid(HF1Grid),
   theHF1Weights(HF1Weights),
   theHF2Grid(HF2Grid),
   theHF2Weights(HF2Weights),
   theEBweight(EBweight),
   theEEweight(EEweight),
   theHBweight(HBweight),
   theHESweight(HESweight),
   theHEDweight(HEDweight),
   theHOweight(HOweight),
   theHF1weight(HF1weight),
   theHF2weight(HF2weight),
   theEcutTower(EcutTower),
   theEBSumThreshold(EBSumThreshold),
   theEESumThreshold(EESumThreshold),
   theHOIsUsed(useHO),
   // (momentum reconstruction algorithm)
   theMomConstrMethod(momConstrMethod),
   theMomEmDepth(momEmDepth),
   theMomHadDepth(momHadDepth),
   theMomTotDepth(momTotDepth)
{
}


void CaloTowersCreationAlgo::setGeometry(const CaloTowerConstituentsMap* ctt, const HcalTopology* topo, const CaloGeometry* geo) {
  theTowerConstituentsMap=ctt;
  theHcalTopology = topo;
  theGeometry = geo;
  theTowerGeometry=geo->getSubdetectorGeometry(DetId::Calo,CaloTowerDetId::SubdetId);
}

void CaloTowersCreationAlgo::begin() {
  theTowerMap.clear();
}

void CaloTowersCreationAlgo::process(const HBHERecHitCollection& hbhe) { 
  for(HBHERecHitCollection::const_iterator hbheItr = hbhe.begin();
      hbheItr != hbhe.end(); ++hbheItr)
    assignHit(&(*hbheItr));
}

void CaloTowersCreationAlgo::process(const HORecHitCollection& ho) { 
  for(HORecHitCollection::const_iterator hoItr = ho.begin();
      hoItr != ho.end(); ++hoItr)
    assignHit(&(*hoItr));
}  

void CaloTowersCreationAlgo::process(const HFRecHitCollection& hf) { 
  for(HFRecHitCollection::const_iterator hfItr = hf.begin();
      hfItr != hf.end(); ++hfItr)  
    assignHit(&(*hfItr));
}

void CaloTowersCreationAlgo::process(const EcalRecHitCollection& ec) { 
  for(EcalRecHitCollection::const_iterator ecItr = ec.begin();
      ecItr != ec.end(); ++ecItr)  
    assignHit(&(*ecItr));
}
void CaloTowersCreationAlgo::process(const CaloTowerCollection& ctc) {
  for(CaloTowerCollection::const_iterator ctcItr = ctc.begin();
      ctcItr != ctc.end(); ++ctcItr) { 
    rescale(&(*ctcItr));
    }
}



void CaloTowersCreationAlgo::finish(CaloTowerCollection& result) {
  // now copy this map into the final collection
  for(MetaTowerMap::const_iterator mapItr = theTowerMap.begin();
      mapItr != theTowerMap.end(); ++ mapItr) {
    CaloTower ct=convert(mapItr->first,mapItr->second);
    if (ct.constituentsSize()>0 && ct.energy()>theEcutTower) {
      result.push_back(ct);
    }
  }
  theTowerMap.clear(); // save the memory
}


void CaloTowersCreationAlgo::rescaleTowers(const CaloTowerCollection& ctc, CaloTowerCollection& ctcResult) {

    for (CaloTowerCollection::const_iterator ctcItr = ctc.begin();
      ctcItr != ctc.end(); ++ctcItr) { 
      
      CaloTowerDetId  twrId = ctcItr->id(); 
      double newE_em    = ctcItr->emEnergy();
      double newE_had   = ctcItr->hadEnergy();
      double newE_outer = ctcItr->outerEnergy(); 

      double threshold = 0.0; // not used: we do not change thresholds
      double weight    = 1.0;

      // HF
      if (ctcItr->ietaAbs()>=30) {
        double E_short = 0.5 * newE_had;             // from the definitions for HF
        double E_long  = newE_em + 0.5 * newE_had;   //
        // scale
        E_long  *= theHF1weight;
        E_short *= theHF2weight;
        // convert
        newE_em  = E_long - E_short;
        newE_had = 2.0 * E_short;
      }

      else {   // barrel/endcap

        // find if its in EB, or EE; determine from first ecal constituent found
        for (uint iConst = 0; iConst < ctcItr->constituentsSize(); ++iConst) {
          DetId constId = ctcItr->constituent(iConst);
          if (constId.det()!=DetId::Ecal) continue;
          getThresholdAndWeight(constId, threshold, weight);
          newE_em *= weight;
          break;
        }
        // HO
        for (uint iConst = 0; iConst < ctcItr->constituentsSize(); ++iConst) {
          DetId constId = ctcItr->constituent(iConst);
          if (constId.det()!=DetId::Hcal) continue;
          if (HcalDetId(constId).subdet()!=HcalOuter) continue;
          getThresholdAndWeight(constId, threshold, weight);
          newE_outer *= weight;
          break;
        }
        // HB/HE
        for (uint iConst = 0; iConst < ctcItr->constituentsSize(); ++iConst) {
          DetId constId = ctcItr->constituent(iConst);
          if (constId.det()!=DetId::Hcal) continue;
          if (HcalDetId(constId).subdet()==HcalOuter) continue;
          getThresholdAndWeight(constId, threshold, weight);
          newE_had *= weight;
          if (ctcItr->ietaAbs()>16) newE_outer *= weight;
          break;
        }
        
    }   // barrel/endcap region

    // now make the new tower

    double newE_hadTot = (theHOIsUsed &&  twrId.ietaAbs()<16)? newE_had+newE_outer : newE_had;

    GlobalPoint  emPoint = ctcItr->emPosition(); 
    GlobalPoint hadPoint = ctcItr->emPosition(); 

    double f_em  = 1.0/cosh(emPoint.eta());
    double f_had = 1.0/cosh(hadPoint.eta());

    CaloTower::LorentzVector towerP4;

    if (ctcItr->ietaAbs()<30) {
      if (newE_em>0)     towerP4 += CaloTower::PolarLorentzVector(newE_em*f_em,   emPoint.eta(),  emPoint.phi(),  0); 
      if (newE_hadTot>0) towerP4 += CaloTower::PolarLorentzVector(newE_hadTot*f_had, hadPoint.eta(), hadPoint.phi(), 0); 
    }
    else {
      double newE_tot = newE_em + newE_had;
      // for HF we use common point for ecal, hcal shower positions regardless of the method
      if (newE_tot>0) towerP4 += CaloTower::PolarLorentzVector(newE_tot*f_had, hadPoint.eta(), hadPoint.phi(), 0);
    }



    CaloTower rescaledTower(twrId, newE_em, newE_had, newE_outer, -1, -1, towerP4, emPoint, hadPoint);
    // copy the timings, have to convert back to int, 1 unit = 0.01 ns
    rescaledTower.setEcalTime( int(ctcItr->ecalTime()*100.0 + 0.5) );
    rescaledTower.setHcalTime( int(ctcItr->hcalTime()*100.0 + 0.5) );

    std::vector<DetId> contains;
    for (uint iConst = 0; iConst < ctcItr->constituentsSize(); ++iConst) {
      contains.push_back(ctcItr->constituent(iConst));
    }
    rescaledTower.addConstituents(contains);

    ctcResult.push_back(rescaledTower);

    } // end of loop over towers


}







void CaloTowersCreationAlgo::assignHit(const CaloRecHit * recHit) {
  DetId detId = recHit->detid();
  double threshold, weight;
  getThresholdAndWeight(detId, threshold, weight);

  // SPECIAL handling of tower 28/depth 3 --> half into tower 28 and half into tower 29
  if (detId.det()==DetId::Hcal && 
      HcalDetId(detId).subdet()==HcalEndcap &&
      HcalDetId(detId).depth()==3 &&
      HcalDetId(detId).ietaAbs()==28) {

    CaloTowerDetId towerDetId = theTowerConstituentsMap->towerOf(detId);
    if (towerDetId.null()) return;    
    MetaTower & tower28 = find(towerDetId);    
    CaloTowerDetId towerDetId29 = CaloTowerDetId(towerDetId.ieta()+
						 towerDetId.zside(),
						 towerDetId.iphi());
    MetaTower & tower29 = find(towerDetId29);

    double energy = recHit->energy()/2; // NOTE DIVIDE BY 2!!!
    if(energy >= threshold) {
      double e28 = energy * weight;
      double e29 = energy * weight;
      
      tower28.E_had += e28;
      tower28.E += e28;
      std::pair<DetId,double> mc(detId,e28);
      tower28.metaConstituents.push_back(mc);

      tower29.E_had += e29;
      tower29.E += e29;
      tower29.metaConstituents.push_back(mc);

      // store the energy in layer 3 also in E_outer
      tower28.E_outer += e28;
      tower29.E_outer += e29;

    }
  } else {
    CaloTowerDetId towerDetId = theTowerConstituentsMap->towerOf(detId);
    if (towerDetId.null()) return;    
    MetaTower & tower = find(towerDetId);

    double energy = recHit->energy();
    if(energy >= threshold) {
      // TODO calculate crystal by crystal
      double e = energy * weight;
      
      DetId::Detector det = detId.det();
      if(det == DetId::Ecal) {
        tower.E_em += e;
        tower.E += e;
        // do not use "recovered" hits in time calculation
        if (recHit->time() != EcalRecHit::kRECOVERED) {
          tower.emSumTimeTimesE += ( e * recHit->time() );
          tower.emSumEForTime   += e;  // see above
        }
      }
      // HCAL
      else {
        HcalDetId hcalDetId(detId);
        if(hcalDetId.subdet() == HcalOuter) {
          tower.E_outer += e;
          if(theHOIsUsed) tower.E += e;
        } 
        // HF calculates EM fraction differently
        else if(hcalDetId.subdet() == HcalForward) {
          if(hcalDetId.depth() == 1) {
            // long fiber, so E_EM = E(Long) - E(Short)
            tower.E_em += e;
          } 
          else {
            // short fiber, EHAD = 2 * E(Short)
            tower.E_em -= e;
            tower.E_had += 2. * e;
          }
          tower.E += e;

	  // put the timing only in HCAL
	  tower.hadSumTimeTimesE += ( e * recHit->time() );
	  tower.hadSumEForTime   += e;
        }
        else {
          // HCAL situation normal
          tower.E_had += e;
          tower.E += e;

      	  // time info
          tower.hadSumTimeTimesE += ( e * recHit->time() );
          tower.hadSumEForTime   += e;

          // store energy for depth 2 for towers 18-27
          if (HcalDetId(detId).subdet()==HcalEndcap & HcalDetId(detId).depth()==2 &&
            HcalDetId(detId).ietaAbs()>=18 && HcalDetId(detId).ietaAbs()<=27) {
              tower.E_outer += e;
          }
        }
      }
      std::pair<DetId,double> mc(detId,e);
      tower.metaConstituents.push_back(mc);
    } 
  }
}


// This method is not flexible enough for the new CaloTower format. 
// For now make a quick compatibility "fix" : WILL NOT WORK CORRECTLY with anything 
// except the default simple p4 assignment!!!
// Must be rewritten for full functionality.
void CaloTowersCreationAlgo::rescale(const CaloTower * ct) {
  double threshold, weight;

  CaloTowerDetId towerDetId = ct->id();
//  if (towerDetId.null()) return;    
  MetaTower & tower = find(towerDetId);

  tower.E_em = 0.;
  tower.E_had = 0.;
  tower.E_outer = 0.;
  for (unsigned int i=0; i<ct->constituentsSize(); i++) {
    DetId detId = ct->constituent(i);
    getThresholdAndWeight(detId, threshold, weight);
    DetId::Detector det = detId.det();
    if(det == DetId::Ecal) {
      tower.E_em = ct->emEnergy()*weight;
    }
    else {
      HcalDetId hcalDetId(detId);
      if(hcalDetId.subdet() == HcalForward) {
        if (hcalDetId.depth()==1) tower.E_em = ct->emEnergy()*weight;
        if (hcalDetId.depth()==2) tower.E_had = ct->hadEnergy()*weight;
      }
      else if(hcalDetId.subdet() == HcalOuter) {
        tower.E_outer = ct->outerEnergy()*weight;
      }
      else {
        tower.E_had = ct->hadEnergy()*weight;
      }
    }
    tower.E = tower.E_had+tower.E_em+tower.E_outer;

    // this is to be compliant with the new MetaTower setup
    // used only for the default simple vector assignment
    std::pair<DetId, double> mc(detId, 0);
    tower.metaConstituents.push_back(mc);
  }

  // preserve time inforamtion
  tower.emSumTimeTimesE  = ct->ecalTime();
  tower.hadSumTimeTimesE = ct->hcalTime();
  tower.emSumEForTime = 1.0;
  tower.hadSumEForTime = 1.0;
}

CaloTowersCreationAlgo::MetaTower::MetaTower() : 
E(0),E_em(0),E_had(0),E_outer(0), emSumTimeTimesE(0), hadSumTimeTimesE(0), emSumEForTime(0), hadSumEForTime(0) { }


CaloTowersCreationAlgo::MetaTower & CaloTowersCreationAlgo::find(const CaloTowerDetId & detId) {
  MetaTowerMap::iterator itr = theTowerMap.find(detId);
  if(itr == theTowerMap.end()) {
    // need to build a new tower
    MetaTower t;

    // store it in the map
    theTowerMap.insert(std::pair<CaloTowerDetId, CaloTowersCreationAlgo::MetaTower>(detId, t));
    itr = theTowerMap.find(detId);
  }
  return itr->second;
}

CaloTower CaloTowersCreationAlgo::convert(const CaloTowerDetId& id, const MetaTower& mt) {

    double ecalThres=(id.ietaAbs()<=17)?(theEBSumThreshold):(theEESumThreshold);
    double E=mt.E;
    double E_em=mt.E_em;
    double E_had=mt.E_had;
    double E_outer=mt.E_outer;

    // Note: E_outer is used to save HO energy OR energy in the outermost depths in endcap region
    // In the methods with separate treatment of EM and HAd components:
    //  - HO is not used to determine direction, however HO energy is added to get "total had energy"
    //  => Check if the tower is within HO coverage before adding E_outer to the "total had" energy
    //     else the energy will be double counted
    // When summing up the energy of the tower these checks are performed in the loops over RecHits


    float  ecalTime = (mt.emSumEForTime>0)?   mt.emSumTimeTimesE/mt.emSumEForTime  : -9999;
    float  hcalTime = (mt.hadSumEForTime>0)?  mt.hadSumTimeTimesE/mt.hadSumEForTime : -9999;

    std::vector<std::pair<DetId,double> > metaContains=mt.metaConstituents;

    if (id.ietaAbs()<=29 && E_em<ecalThres) { // ignore EM threshold in HF
      E-=E_em;
      E_em=0;
      std::vector<std::pair<DetId,double> > metaContains_noecal;

    for (std::vector<std::pair<DetId,double> >::iterator i=metaContains.begin(); i!=metaContains.end(); ++i) 
	        if (i->first.det()!=DetId::Ecal) metaContains_noecal.push_back(*i);
      metaContains.swap(metaContains_noecal);
    }

    if (id.ietaAbs()<=29 && E_had<theHcalThreshold) {
      E-=E_had;

      if (theHOIsUsed && id.ietaAbs()<16)  E-=E_outer; // not subtracted before, think it should be done
     
      E_had=0;
      E_outer=0;
      std::vector<std::pair<DetId,double> > metaContains_nohcal;

      for (std::vector<std::pair<DetId,double> >::iterator i=metaContains.begin(); i!=metaContains.end(); ++i) 
        if (i->first.det()!=DetId::Hcal) metaContains_nohcal.push_back(*i);
      metaContains.swap(metaContains_nohcal);
    }

    // create CaloTower using the selected algorithm

    GlobalPoint emPoint, hadPoint;
    CaloTower::LorentzVector towerP4;

  switch (theMomConstrMethod) {

  case 0 :
    {  // Simple 4-momentum assignment
      GlobalPoint p=theTowerGeometry->getGeometry(id)->getPosition();

      double pf=1.0/cosh(p.eta());
      if (E>0) towerP4 = CaloTower::PolarLorentzVector(E*pf, p.eta(), p.phi(), 0);
      
      emPoint  = p;   
      hadPoint = p;
    }  // end case 0
    break;

  case 1 :
    {   // separate 4-vectors for ECAL, HCAL, add to get the 4-vector of the tower (=>tower has mass!)
      if (id.ietaAbs()<=29) {
        if (E_em>0) {
          emPoint   = emShwrPos(metaContains, theMomEmDepth, E_em);
          double emPf = 1.0/cosh(emPoint.eta());
          towerP4 += CaloTower::PolarLorentzVector(E_em*emPf, emPoint.eta(), emPoint.phi(), 0); 
        }
        if (E_had>0) {
          double E_had_tot = (theHOIsUsed && id.ietaAbs()<16)? E_had+E_outer : E_had;
//          hadPoint  = hadShwrPos(metaContains, theMomHadDepth, E_had_tot);
          hadPoint  = hadShwrPos(id, theMomHadDepth);
          double hadPf = 1.0/cosh(hadPoint.eta());
          towerP4 += CaloTower::PolarLorentzVector(E_had_tot*hadPf, hadPoint.eta(), hadPoint.phi(), 0); 
        }
      }
      else {  // forward detector: use the CaloTower position 
        GlobalPoint p=theTowerGeometry->getGeometry(id)->getPosition();
        double pf=1.0/cosh(p.eta());
        if (E>0) towerP4 = CaloTower::PolarLorentzVector(E*pf, p.eta(), p.phi(), 0);  // simple momentum assignment, same position
        emPoint  = p;   
        hadPoint = p;
      }
    }  // end case 1
    break;

  case 2:
    {   // use ECAL, HCAL shower position to get weighted overall position, assign full energy to get the 4-vector of the tower (massless tower)
      if (id.ietaAbs()<=29) {
        if (E_em>0)  emPoint = emShwrPos(metaContains, theMomEmDepth, E_em);
        double E_had_tot = (theHOIsUsed && id.ietaAbs()<16)? E_had+E_outer : E_had; 
        if (E_had>0) hadPoint  = hadShwrPos(metaContains, theMomHadDepth, E_had_tot);

        // common point for EM/HAD components based on predefined depths
        GlobalPoint p = GlobalPoint( 
          (E_em*emPoint.x()+E_had_tot*hadPoint.x())/E,
          (E_em*emPoint.y()+E_had_tot*hadPoint.y())/E,
          (E_em*emPoint.z()+E_had_tot*hadPoint.z())/E
          );

        double sumPf = 1.0/cosh(p.eta());
        if (E>0) towerP4 = CaloTower::PolarLorentzVector(E*sumPf, p.eta(), p.phi(), 0); 

        emPoint  = p;   
        hadPoint = p;

      }
      else {  // forward detector: use the CaloTower position 
        GlobalPoint p=theTowerGeometry->getGeometry(id)->getPosition();
        double pf=1.0/cosh(p.eta());
        if (E>0) towerP4 = CaloTower::PolarLorentzVector(E*pf, p.eta(), p.phi(), 0);  // simple momentum assignment, same position
        emPoint  = p;   
        hadPoint = p;
      }
    }  // end case 2
    break;

  case 3:
    {   // separate 4-vectors for ECAL crystals, and HCAL; add to get the 4-vector of the tower (tower has mass!)
      if (id.ietaAbs()<=29) {
        if (E_em>0) {
          
          emPoint   = emShwrPos(metaContains, theMomEmDepth, E_em);
          for (std::vector<std::pair<DetId,double> >::iterator mc_it = metaContains.begin(); 
            mc_it!=metaContains.end(); ++mc_it) {
              if (mc_it->first.det() == DetId::Ecal) {
                GlobalPoint p = emCrystalShwrPos(mc_it->first.det(), theMomEmDepth);
                double pf=1.0/cosh(p.eta());
                double e = mc_it->second;
                towerP4 += CaloTower::PolarLorentzVector(e*pf, p.eta(), p.phi(), 0); 
              }
          }                    
          emPoint   = emShwrPos(metaContains, theMomEmDepth, E_em);
        }
        if (E_had>0) {
          double E_had_tot = (theHOIsUsed && id.ietaAbs()<16)? E_had+E_outer : E_had; 
//          hadPoint = hadShwrPos(metaContains, theMomHadDepth, E_had_tot);
          hadPoint  = hadShwrPos(id, theMomHadDepth);
          double hadPf = 1.0/cosh(hadPoint.eta());
          towerP4 += CaloTower::PolarLorentzVector(E_had_tot*hadPf, hadPoint.eta(), hadPoint.phi(), 0); 
        }
      }
      else {  // forward detector: use the CaloTower position 
        GlobalPoint p=theTowerGeometry->getGeometry(id)->getPosition();
        double pf=1.0/cosh(p.eta());
        towerP4 = CaloTower::PolarLorentzVector(E*pf, p.eta(), p.phi(), 0);  // simple momentum assignment, same position
        emPoint  = p;   
        hadPoint = p;
      }
    }  // end case 3
    break;

  case 4:
    {   // use ECAL position for the tower (when E_cal>0), else default CaloTower position (massless tower)
      if (id.ietaAbs()<=29) {
        if (E_em>0)  emPoint = emShwrLogWeightPos(metaContains, theMomEmDepth, E_em);
        else emPoint = theTowerGeometry->getGeometry(id)->getPosition();

        double sumPf = 1.0/cosh(emPoint.eta());
        if (E>0) towerP4 = CaloTower::PolarLorentzVector(E*sumPf, emPoint.eta(), emPoint.phi(), 0); 
        
        hadPoint = emPoint;
      }
      else {  // forward detector: use the CaloTower position 
        GlobalPoint p=theTowerGeometry->getGeometry(id)->getPosition();
        double pf=1.0/cosh(p.eta());
        if (E>0) towerP4 = CaloTower::PolarLorentzVector(E*pf, p.eta(), p.phi(), 0);  // simple momentum assignment, same position
        emPoint  = p;   
        hadPoint = p;
      }
    }   // end case 4
    break;

  }  // end of decision on p4 reconstruction method


//    CaloTower::LorentzVector lv = caloTowerMomentum(id, metaContains, E, E_em, E_had, E_outer);

    CaloTower retval(id, E_em, E_had, E_outer, -1, -1, towerP4, emPoint, hadPoint);

    // set the timings
    retval.setEcalTime(compactTime(ecalTime));
    retval.setHcalTime(compactTime(hcalTime));

    std::vector<DetId> contains;
    for (std::vector<std::pair<DetId,double> >::iterator i=metaContains.begin(); i!=metaContains.end(); ++i) 
        contains.push_back(i->first);

    retval.addConstituents(contains);
    return retval;
} 



void CaloTowersCreationAlgo::getThresholdAndWeight(const DetId & detId, double & threshold, double & weight) const {
  DetId::Detector det = detId.det();
  weight=0; // in case the hit is not identified

  if(det == DetId::Ecal) {
    // may or may not be EB.  We'll find out.

    EcalSubdetector subdet = (EcalSubdetector)(detId.subdetId());
    if(subdet == EcalBarrel) {
      threshold = theEBthreshold;
      weight = theEBweight;
      if (weight <= 0.) {
        ROOT::Math::Interpolator my(theEBGrid,theEBWeights,ROOT::Math::Interpolation::AKIMA);
        weight = my.Eval(theEBEScale);
      }
    }
    else if(subdet == EcalEndcap) {
      threshold = theEEthreshold;
      weight = theEEweight;
      if (weight <= 0.) {
        ROOT::Math::Interpolator my(theEEGrid,theEEWeights,ROOT::Math::Interpolation::AKIMA);
        weight = my.Eval(theEEEScale);
      }
    }
  }
  else if(det == DetId::Hcal) {
    HcalDetId hcalDetId(detId);
    HcalSubdetector subdet = hcalDetId.subdet();
    
    if(subdet == HcalBarrel) {
      threshold = theHBthreshold;
      weight = theHBweight;
      if (weight <= 0.) {
        ROOT::Math::Interpolator my(theHBGrid,theHBWeights,ROOT::Math::Interpolation::AKIMA);
        weight = my.Eval(theHBEScale);
      }
    }
    
    else if(subdet == HcalEndcap) {
      // check if it's single or double tower
      if(hcalDetId.ietaAbs() < theHcalTopology->firstHEDoublePhiRing()) {
        threshold = theHESthreshold;
        weight = theHESweight;
        if (weight <= 0.) {
          ROOT::Math::Interpolator my(theHESGrid,theHESWeights,ROOT::Math::Interpolation::AKIMA);
          weight = my.Eval(theHESEScale);
        }
      }
      else {
        threshold = theHEDthreshold;
        weight = theHEDweight;
        if (weight <= 0.) {
          ROOT::Math::Interpolator my(theHEDGrid,theHEDWeights,ROOT::Math::Interpolation::AKIMA);
          weight = my.Eval(theHEDEScale);
        }
      }
    } else if(subdet == HcalOuter) {
      threshold = theHOthreshold;
      weight = theHOweight;
      if (weight <= 0.) {
        ROOT::Math::Interpolator my(theHOGrid,theHOWeights,ROOT::Math::Interpolation::AKIMA);
        weight = my.Eval(theHOEScale);
      }
    } else if(subdet == HcalForward) {
      if(hcalDetId.depth() == 1) {
        threshold = theHF1threshold;
        weight = theHF1weight;
        if (weight <= 0.) {
          ROOT::Math::Interpolator my(theHF1Grid,theHF1Weights,ROOT::Math::Interpolation::AKIMA);
          weight = my.Eval(theHF1EScale);
        }
      } else {
        threshold = theHF2threshold;
        weight = theHF2weight;
        if (weight <= 0.) {
          ROOT::Math::Interpolator my(theHF2Grid,theHF2Weights,ROOT::Math::Interpolation::AKIMA);
          weight = my.Eval(theHF2EScale);
        }
      }
    }
  }
  else {
    std::cout << "BAD CELL det " << det << std::endl;
  }
}

void CaloTowersCreationAlgo::setEBEScale(double scale){
  if (scale>0.00001) *&theEBEScale = scale;
  else *&theEBEScale = 50.;
}

void CaloTowersCreationAlgo::setEEEScale(double scale){
  if (scale>0.00001) *&theEEEScale = scale;
  else *&theEEEScale = 50.;
}

void CaloTowersCreationAlgo::setHBEScale(double scale){
  if (scale>0.00001) *&theHBEScale = scale;
  else *&theHBEScale = 50.;
}

void CaloTowersCreationAlgo::setHESEScale(double scale){
  if (scale>0.00001) *&theHESEScale = scale;
  else *&theHESEScale = 50.;
}

void CaloTowersCreationAlgo::setHEDEScale(double scale){
  if (scale>0.00001) *&theHEDEScale = scale;
  else *&theHEDEScale = 50.;
}

void CaloTowersCreationAlgo::setHOEScale(double scale){
  if (scale>0.00001) *&theHOEScale = scale;
  else *&theHOEScale = 50.;
}

void CaloTowersCreationAlgo::setHF1EScale(double scale){
  if (scale>0.00001) *&theHF1EScale = scale;
  else *&theHF1EScale = 50.;
}

void CaloTowersCreationAlgo::setHF2EScale(double scale){
  if (scale>0.00001) *&theHF2EScale = scale;
  else *&theHF2EScale = 50.;
}


GlobalPoint CaloTowersCreationAlgo::emCrystalShwrPos(DetId detId, float fracDepth) {
   const CaloCellGeometry* cellGeometry = theGeometry->getGeometry(detId);
   GlobalPoint point = cellGeometry->getPosition();  // face of the cell

   if      (fracDepth<0) fracDepth=0;
   else if (fracDepth>1) fracDepth=1;

   if (fracDepth>0.0) {
     CaloCellGeometry::CornersVec cv = cellGeometry->getCorners();
     GlobalPoint backPoint = GlobalPoint( 0.25*( cv[4].x() + cv[5].x() + cv[6].x() + cv[7].x() ),
                                          0.25*( cv[4].y() + cv[5].y() + cv[6].y() + cv[7].y() ),
                                          0.25*( cv[4].z() + cv[5].z() + cv[6].z() + cv[7].z() ) );
     point += fracDepth * (backPoint-point);
   }

   return point;
}

GlobalPoint CaloTowersCreationAlgo::hadSegmentShwrPos(DetId detId, float fracDepth) {
   const CaloCellGeometry* cellGeometry = theGeometry->getGeometry(detId);
   GlobalPoint point = cellGeometry->getPosition();  // face of the cell

   if      (fracDepth<0) fracDepth=0;
   else if (fracDepth>1) fracDepth=1;

   if (fracDepth>0.0) {
     CaloCellGeometry::CornersVec cv = cellGeometry->getCorners();
     GlobalPoint backPoint = GlobalPoint( 0.25*( cv[4].x() + cv[5].x() + cv[6].x() + cv[7].x() ),
                                          0.25*( cv[4].y() + cv[5].y() + cv[6].y() + cv[7].y() ),
                                          0.25*( cv[4].z() + cv[5].z() + cv[6].z() + cv[7].z() ) );
     point += fracDepth * (backPoint-point);
   }

   return point;
}


GlobalPoint CaloTowersCreationAlgo::hadShwrPos(std::vector<std::pair<DetId,double> >& metaContains,
                                               float fracDepth, double hadE) {
                                                  
  // this is based on available RecHits, can lead to different actual depths if
  // hits in multi-depth towers are not all there
  if (hadE<=0) return GlobalPoint(0,0,0);

  double hadX = 0.0;
  double hadY = 0.0;
  double hadZ = 0.0;

  int nConst = 0;

  std::vector<std::pair<DetId,double> >::iterator mc_it = metaContains.begin();
  for (; mc_it!=metaContains.end(); ++mc_it) {
    if (mc_it->first.det() != DetId::Hcal) continue;
    // do not use HO for deirection calculations for now
    if (HcalDetId(mc_it->first).subdet() == HcalOuter) continue;
    ++nConst;

    GlobalPoint p = hadSegmentShwrPos(mc_it->first, fracDepth);

    // longitudinal segmentation: do not weight by energy,
    // get the geometrical position
    hadX += p.x();
    hadY += p.y();
    hadZ += p.z();
  }

   return GlobalPoint(hadX/nConst, hadY/nConst, hadZ/nConst);
}


GlobalPoint CaloTowersCreationAlgo::hadShwrPos(CaloTowerDetId towerId, float fracDepth) {

  // set depth using geometry of cells that are associated with the
  // tower (regrdeleess if they have non-zero energies)

//  if (hadE <= 0) return GlobalPoint(0, 0, 0);

  if (fracDepth < 0) fracDepth = 0;
  else if (fracDepth > 1) fracDepth = 1;

  GlobalPoint point(0,0,0);

  int iEta = towerId.ieta();
  int iPhi = towerId.iphi();

  HcalDetId frontCellId, backCellId;

  if (towerId.ietaAbs() <= 14) {
    // barrel, one depth only
    frontCellId = HcalDetId(HcalBarrel, iEta, iPhi, 1);
    backCellId  = HcalDetId(HcalBarrel, iEta, iPhi, 1);
  }
  else if (towerId.ietaAbs() == 15) {
    // barrel, two depths
    frontCellId = HcalDetId(HcalBarrel, iEta, iPhi, 1);
    backCellId  = HcalDetId(HcalBarrel, iEta, iPhi, 2);
  }
  else if (towerId.ietaAbs() == 16) {
    // barrel and endcap: two depths HB, one depth HE 
    frontCellId = HcalDetId(HcalBarrel, iEta, iPhi, 1);
    backCellId  = HcalDetId(HcalEndcap, iEta, iPhi, 3);  // this cell is in endcap!
  }
  else if (towerId.ietaAbs() == 17) {
    // endcap, one depth only
   frontCellId = HcalDetId(HcalEndcap, iEta, iPhi, 1);
   backCellId  = HcalDetId(HcalEndcap, iEta, iPhi, 1);
  }
  else if (towerId.ietaAbs() >= 18 && towerId.ietaAbs() <= 26) {
  // endcap: two depths
  frontCellId = HcalDetId(HcalEndcap, iEta, iPhi, 1);
  backCellId  = HcalDetId(HcalEndcap, iEta, iPhi, 2);
  }
  else if (towerId.ietaAbs() <= 29) {
  // endcap: three depths
  frontCellId = HcalDetId(HcalEndcap, iEta, iPhi, 1);
  // there is no iEta=29 for depth 3
  if (iEta ==  29) iEta =  28;
  if (iEta == -29) iEta = -28;
  backCellId  = HcalDetId(HcalEndcap, iEta, iPhi, 3);
  }
  else if (towerId.ietaAbs() >= 30) {
  // forward, take the goemetry for long fibers
  frontCellId = HcalDetId(HcalForward, iEta, iPhi, 1);
  backCellId  = HcalDetId(HcalForward, iEta, iPhi, 1);
  }
  else {
    // should not get here
    return point;
  }

  point = hadShwPosFromCells(DetId(frontCellId), DetId(backCellId), fracDepth);

  return point;
}

GlobalPoint CaloTowersCreationAlgo::hadShwPosFromCells(DetId frontCellId, DetId backCellId, float fracDepth) {

   // uses the "front" and "back" cells
   // to determine the axis. point set by the predefined depth.

    const CaloCellGeometry* frontCellGeometry = theGeometry->getGeometry(DetId(frontCellId));
    const CaloCellGeometry* backCellGeometry  = theGeometry->getGeometry(DetId(backCellId));

    GlobalPoint point = frontCellGeometry->getPosition();

    CaloCellGeometry::CornersVec cv = backCellGeometry->getCorners();

    GlobalPoint backPoint = GlobalPoint(0.25 * (cv[4].x() + cv[5].x() + cv[6].x() + cv[7].x()),
      0.25 * (cv[4].y() + cv[5].y() + cv[6].y() + cv[7].y()),
      0.25 * (cv[4].z() + cv[5].z() + cv[6].z() + cv[7].z()));

    point += fracDepth * (backPoint - point);

    return point;
}


GlobalPoint CaloTowersCreationAlgo::emShwrPos(std::vector<std::pair<DetId,double> >& metaContains, 
                                              float fracDepth, double emE) {

  if (emE<=0) return GlobalPoint(0,0,0);

  double emX = 0.0;
  double emY = 0.0;
  double emZ = 0.0;

  double eSum = 0;

  std::vector<std::pair<DetId,double> >::iterator mc_it = metaContains.begin();
  for (; mc_it!=metaContains.end(); ++mc_it) {
    if (mc_it->first.det() != DetId::Ecal) continue;
    GlobalPoint p = emCrystalShwrPos(mc_it->first, fracDepth);
    double e = mc_it->second;

    if (e>0) {
      emX += p.x() * e;
      emY += p.y() * e;
      emZ += p.z() * e;
      eSum += e;
    }

  }

   return GlobalPoint(emX/eSum, emY/eSum, emZ/eSum);
}


GlobalPoint CaloTowersCreationAlgo::emShwrLogWeightPos(std::vector<std::pair<DetId,double> >& metaContains, 
                               float fracDepth, double emE) {

  double emX = 0.0;
  double emY = 0.0;
  double emZ = 0.0;

  double weight = 0;
  double sumWeights = 0;
  double sumEmE = 0;  // add crystals with E/E_EM > 1.5%
  double crystalThresh = 0.015 * emE;

  std::vector<std::pair<DetId,double> >::iterator mc_it = metaContains.begin();
  for (; mc_it!=metaContains.end(); ++mc_it) {
    if (mc_it->first.det() == DetId::Ecal && mc_it->second > crystalThresh) sumEmE += mc_it->second;
  }

  for (mc_it = metaContains.begin(); mc_it!=metaContains.end(); ++mc_it) {
    
    if (mc_it->first.det() != DetId::Ecal || mc_it->second < crystalThresh) continue;
    
    GlobalPoint p = emCrystalShwrPos(mc_it->first, fracDepth);

    weight = 4.2 + log(mc_it->second/sumEmE);
    sumWeights += weight;
      
    emX += p.x() * weight;
    emY += p.y() * weight;
    emZ += p.z() * weight;
  }

   return GlobalPoint(emX/sumWeights, emY/sumWeights, emZ/sumWeights);
}





int CaloTowersCreationAlgo::compactTime(float time) {

  const float timeUnit = 0.01; // discretization (ns)

  if (time>  300.0) return  30000;
  if (time< -300.0) return -30000;

  return int(time/timeUnit + 0.5);

}
