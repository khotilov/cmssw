#ifndef RECOLOCALCALO_CALOTOWERSCREATOR_CALOTOWERSCREATIONALGO_H
#define RECOLOCALCALO_CALOTOWERSCREATOR_CALOTOWERSCREATIONALGO_H 1

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// channel status
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"


#include <map>
class HcalTopology;
class CaloGeometry;
class CaloSubdetectorGeometry;
class CaloTowerConstituentsMap;
class CaloRecHit;
class DetId;

/** \class CaloTowersCreationAlgo
  *  
  * $Date: 2008/09/03 20:49:30 $
  * $Revision: 1.13 $
  * \author R. Wilkinson - Caltech
  */

//
// Modify MetaTower to save energy of rechits for use in tower 4-momentum assignment,
// added containers for timing assignment and for holding status information.
// Anton Anastassov (Northwestern)
//

class CaloTowersCreationAlgo {
public:
  CaloTowersCreationAlgo();

  CaloTowersCreationAlgo(double EBthreshold, double EEthreshold, double HcalThreshold,
    double HBthreshold, double HESthreshold, double HEDthreshold,
    double HOthreshold, double HF1threshold, double HF2threshold,
    double EBweight, double EEweight,
    double HBweight, double HESweight, double HEDweight, 
    double HOweight, double HF1weight, double HF2weight,
    double EcutTower, double EBSumThreshold, double EESumThreshold, bool useHO,
    // (for momentum reconstruction algorithm)
    int momConstrMethod,
    double momHBDepth,
    double momHEDepth,
    double momEBDepth,
    double momEEDepth
    );
  
  CaloTowersCreationAlgo(double EBthreshold, double EEthreshold, double HcalThreshold,
    double HBthreshold, double HESthreshold, double HEDthreshold,
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
    double EcutTower, double EBSumThreshold, double EESumThreshold, bool useHO,
    // (for momentum reconstruction algorithm)
    int momConstrMethod,
    double momHBDepth,
    double momHEDepth,
    double momEBDepth,
    double momEEDepth
);
  
  void setGeometry(const CaloTowerConstituentsMap* cttopo, const HcalTopology* htopo, const CaloGeometry* geo);

  // pass the containers of channels status from the event record (stored in DB)
  // these are called in  CaloTowersCreator
  void setHcalChStatusFromDB(const HcalChannelQuality* s) { theHcalChStatus = s; }
  void setEcalChStatusFromDB(const EcalChannelStatus* s) { theEcalChStatus = s; }

  // make maps based on information found in the DB
  void makeHcalDeadChMap();
  void makeEcalDeadChMap();

  void begin();
  void process(const HBHERecHitCollection& hbhe);
  void process(const HORecHitCollection& ho);
  void process(const HFRecHitCollection& hf); 
  void process(const EcalRecHitCollection& ecal); 
  
  
  void process(const CaloTowerCollection& ctc);

  void finish(CaloTowerCollection& destCollection);

  // modified rescale method
  void rescaleTowers(const CaloTowerCollection& ctInput, CaloTowerCollection& ctResult);

  void setEBEScale(double scale);
  void setEEEScale(double scale);
  void setHBEScale(double scale);
  void setHESEScale(double scale);
  void setHEDEScale(double scale);
  void setHOEScale(double scale);
  void setHF1EScale(double scale);
  void setHF2EScale(double scale);


  // Assign to categories based on info from DB and RecHit status
  // Called in assignHit to check if the energy should be added to
  // calotower, and how to flag the channel
  uint hbheChanStatusForCaloTower(const CaloRecHit* hit);
  uint hfChanStatusForCaloTower(const CaloRecHit* hit);
  uint hoChanStatusForCaloTower(const CaloRecHit* hit);
  uint ecalChanStatusForCaloTower(const CaloRecHit* hit);

  // Channel flagging is based on acceptable severity levels specified in the
  // configuration file. These methods are used to pass the values read in
  // CaloTowersCreator
  // 
  // from DB
  void setHbheAcceptSevLevelDb(uint level) {theHbheAcceptSevLevelDb = level;} 
  void setHfAcceptSevLevelDb(uint level) {theHfAcceptSevLevelDb = level;} 
  void setHoAcceptSevLevelDb(uint level)  {theHoAcceptSevLevelDb = level;} 
  void setEcalAcceptSevLevelDb(uint level)  {theEcalAcceptSevLevelDb = level;} 
  // from the RecHit
  void setHbheAcceptSevLevelRecHit(uint level) {theHbheAcceptSevLevelRecHit = level; } 
  void setHfAcceptSevLevelRecHit(uint level) {theHfAcceptSevLevelRecHit = level; } 
  void setHoAcceptSevLevelRecHit(uint level) {theHoAcceptSevLevelRecHit = level; } 
  void setEcalAcceptSevLevelRecHit(uint level){theEcalAcceptSevLevelRecHit = level; } 
  // flag to use recovered hits
  void setRecovHbheIsUsed(bool flag) {theRecovHbheIsUsed = flag; };
  void setRecovHoIsUsed(bool flag) {theRecovHoIsUsed = flag; };
  void setRecovHfIsUsed(bool flag) {theRecovHfIsUsed = flag; };
  void setRecovEcalIsUsed(bool flag) {theRecovEcalIsUsed = flag; };


  // Add methods to get the seperate positions for ECAL/HCAL 
  // used in constructing the 4-vectors using new methods
  GlobalPoint emCrystalShwrPos (DetId detId, float fracDepth); 
  GlobalPoint hadSegmentShwrPos(DetId detId, float fracDepth);
  // "effective" point for the EM/HAD shower in CaloTower
  //  position based on non-zero energy cells
  GlobalPoint hadShwrPos(std::vector<std::pair<DetId,double> >& metaContains,
    float fracDepth, double hadE);
  GlobalPoint emShwrPos(std::vector<std::pair<DetId,double> >& metaContains, 
    float fracDepth, double totEmE);

  // overloaded function to get had position based on all had cells in the tower
  GlobalPoint hadShwrPos(CaloTowerDetId id, float fracDepth);
  GlobalPoint hadShwPosFromCells(DetId frontCell, DetId backCell, float fracDepth);

  // for Chris
  GlobalPoint emShwrLogWeightPos(std::vector<std::pair<DetId,double> >& metaContains, 
    float fracDepth, double totEmE);


private:

  struct MetaTower {
    MetaTower();
    double E, E_em, E_had, E_outer;
    // contains also energy of RecHit
    std::vector< std::pair<DetId, double> > metaConstituents;
    double emSumTimeTimesE, hadSumTimeTimesE, emSumEForTime, hadSumEForTime; // Sum(Energy x Timing) : intermediate container

    // needed to set CaloTower status word
    int numBadEcalCells, numRecEcalCells, numProbEcalCells, numBadHcalCells, numRecHcalCells, numProbHcalCells; 

 };

  /// adds a single hit to the tower
  void assignHit(const CaloRecHit * recHit);

  void rescale(const CaloTower * ct);

  /// looks for a given tower in the internal cache.  If it can't find it, it makes it.
  MetaTower & find(const CaloTowerDetId & id);
  
  /// helper method to look up the appropriate threshold & weight
  void getThresholdAndWeight(const DetId & detId, double & threshold, double & weight) const;
  
  double theEBthreshold, theEEthreshold, theHcalThreshold;
  double theHBthreshold, theHESthreshold,  theHEDthreshold; 
  double theHOthreshold, theHF1threshold, theHF2threshold;
  std::vector<double> theEBGrid, theEBWeights;
  std::vector<double> theEEGrid, theEEWeights;
  std::vector<double> theHBGrid, theHBWeights;
  std::vector<double> theHESGrid, theHESWeights;
  std::vector<double> theHEDGrid, theHEDWeights;
  std::vector<double> theHOGrid, theHOWeights;
  std::vector<double> theHF1Grid, theHF1Weights;
  std::vector<double> theHF2Grid, theHF2Weights;
  double theEBweight, theEEweight;
  double theHBweight, theHESweight, theHEDweight, theHOweight, theHF1weight, theHF2weight;
  double theEcutTower, theEBSumThreshold, theEESumThreshold;

  double theEBEScale;
  double theEEEScale;
  double theHBEScale;
  double theHESEScale;
  double theHEDEScale;
  double theHOEScale;
  double theHF1EScale;
  double theHF2EScale;
  const HcalTopology* theHcalTopology;
  const CaloGeometry* theGeometry;
  const CaloTowerConstituentsMap* theTowerConstituentsMap;
  const CaloSubdetectorGeometry* theTowerGeometry;

  // for checking the status of ECAL and HCAL channels stored in the DB 
  const EcalChannelStatus* theEcalChStatus;
  const HcalChannelQuality* theHcalChStatus;


  // fields that hold the information passed from the CaloTowersCreator configuration file:
  // controll what is considered bad/recovered/problematic channel for CaloTower purposes 
  //
  // from DB
  uint theHbheAcceptSevLevelDb;
  uint theHfAcceptSevLevelDb;
  uint theHoAcceptSevLevelDb;
  uint theEcalAcceptSevLevelDb;
  // from the RecHit
  uint theHbheAcceptSevLevelRecHit;
  uint theHfAcceptSevLevelRecHit;
  uint theHoAcceptSevLevelRecHit;
  uint theEcalAcceptSevLevelRecHit;
  // flag to use recovered hits
  bool theRecovHbheIsUsed;
  bool theRecovHoIsUsed;
  bool theRecovHfIsUsed;
  bool theRecovEcalIsUsed;


  /// only affects energy and ET calculation.  HO is still recorded in the tower
  bool theHOIsUsed;

  // Switches and paramters for CaloTower 4-momentum assignment
  // "depth" variables do not affect all algorithms 
  int theMomConstrMethod;
  double theMomEmDepth;
  double theMomHadDepth;

  double theMomHBDepth;
  double theMomHEDepth;
  double theMomEBDepth;
  double theMomEEDepth;

  // compactify timing info
  int compactTime(float time);

  CaloTower convert(const CaloTowerDetId& id, const MetaTower& mt);

  // internal map
  typedef std::map<CaloTowerDetId, MetaTower> MetaTowerMap;
  MetaTowerMap theTowerMap;

  // maps of number of dead channels in towers based in information in the DB 
  // The number of bad channels can be reduced/incremented in the loop over RecHits if the channel
  // was "recovered" (using a HCAL/ECAL recovery algorithm), or had a problem with
  // high severety an event-by-event basis
  //  
  std::map<CaloTowerDetId, int> hcalDeadChMap;
  std::map<CaloTowerDetId, int> ecalDeadChMap;

  // vector of DetId's of channel that were used to make hcalDeadChMap and ecalDeadChMap
  // These vectors are used to prevent double counting: if for some reason such channels
  // make it to the RecHit collections and are reasigned (for example are "recovered")
  // the count of "bad" channels in a tower is reduced.
  //
  std::vector<DetId> hcalBadChanIdInDB;
  std::vector<DetId> ecalBadChanIdInDB;
  
  // clasification of channels in tower construction: the category definition is
  // affected by the setting in the configuration file
  // 
  enum ctHitCategory {GoodChan = 0, BadChan = 1, RecoveredChan = 2, ProblematicChan = 3 };



};

#endif
