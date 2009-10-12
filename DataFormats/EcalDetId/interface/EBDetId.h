#ifndef ECALDETID_EBDETID_H
#define ECALDETID_EBDETID_H

#include <ostream>
#include <cmath>
#include <cstdlib>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"


/** \class EBDetId
 *  Crystal identifier class for the ECAL barrel
 *
 *
 *  $Id: EBDetId.h,v 1.24 2009/01/28 16:23:16 elmer Exp $
 */


class EBDetId : public DetId {
 public:
  enum { Subdet=EcalBarrel};
  /** Constructor of a null id */
  EBDetId() {}
  /** Constructor from a raw value */
  EBDetId(uint32_t rawid) : DetId(rawid) {}
  /** Constructor from crystal ieta and iphi 
      or from SM# and crystal# */
  EBDetId(int index1, int index2, int mode = ETAPHIMODE);
  /** Constructor from a generic cell id */
  EBDetId(const DetId& id);
  /** Assignment operator from cell id */
  EBDetId& operator=(const DetId& id);

  /// get the subdetector .i.e EcalBarrel (what else?)
  // EcalSubdetector subdet() const { return EcalSubdetector(subdetId()); }
  static EcalSubdetector subdet() { return EcalBarrel;}

  /// get the z-side of the crystal (1/-1)
  int zside() const { return (id_&0x10000)?(1):(-1); }
  /// get the absolute value of the crystal ieta
  int ietaAbs() const { return (id_>>9)&0x7F; }
  /// get the crystal ieta
  int ieta() const { return zside()*ietaAbs(); }
  /// get the crystal iphi
  int iphi() const { return id_&0x1FF; }
  /// get the HCAL/trigger ieta of this crystal
  int tower_ieta() const { return ((ietaAbs()-1)/5+1)*zside(); }
  /// get the HCAL/trigger iphi of this crystal
  int tower_iphi() const;
  /// get the HCAL/trigger iphi of this crystal
  EcalTrigTowerDetId tower() const { return EcalTrigTowerDetId(zside(),EcalBarrel,abs(tower_ieta()),tower_iphi()); }
  /// get the ECAL/SM id
  int ism() const;
  /// get the number of module inside the SM (1-4)
  int im() const;
  /// get ECAL/crystal number inside SM
  int ic() const;
  /// get the crystal ieta in the SM convention (1-85)
  int ietaSM() const { return ietaAbs(); }
  /// get the crystal iphi (1-20)
  int iphiSM() const { return (( ic() -1 ) % kCrystalsInPhi ) + 1; }
  
  // is z positive?
  bool positiveZ() const { return id_&0x10000;}
  // crystal number in eta-phi grid
  int numberByEtaPhi() const { 
    return (MAX_IETA + (positiveZ() ? ietaAbs()-1 : -ietaAbs()) )*MAX_IPHI+ iphi()-1;
  }
  // index numbering crystal by SM
  int numberBySM() const; 
  /// get a compact index for arrays
  int hashedIndex() const { return numberByEtaPhi(); }

  uint32_t denseIndex() const { return hashedIndex() ; }

  static bool validDenseIndex( uint32_t din ) { return ( din < kSizeForDenseIndexing ) ; }

  static EBDetId detIdFromDenseIndex( uint32_t di ) { return unhashIndex( di ) ; }

  /// get a DetId from a compact index for arrays
  static EBDetId unhashIndex( int hi ) ;

  static bool validHashIndex(int i) { return !(i<MIN_HASH || i>MAX_HASH); }

  /// check if a valid index combination
  static bool validDetId(int i, int j) ;

  static bool isNextToBoundary(EBDetId id);

  static bool isNextToEtaBoundary(EBDetId id);

  static bool isNextToPhiBoundary(EBDetId id);

  //return the distance in eta units between two EBDetId
  static int distanceEta(const EBDetId& a,const EBDetId& b); 
  //return the distance in phi units between two EBDetId
  static int distancePhi(const EBDetId& a,const EBDetId& b); 

  /// range constants
  static const int MIN_IETA = 1;
  static const int MIN_IPHI = 1;
  static const int MAX_IETA = 85;
  static const int MAX_IPHI = 360;
  static const int kChannelsPerCard = 5;
  static const int kTowersInPhi = 4;  // per SM
  static const int kModulesPerSM = 4;
  static const int kModuleBoundaries[4] ;
  static const int kCrystalsInPhi = 20; // per SM
  static const int kCrystalsInEta = 85; // per SM
  static const int kCrystalsPerSM = 1700;
  static const int MIN_SM = 1;
  static const int MAX_SM = 36;
  static const int MIN_C = 1;
  static const int MAX_C = kCrystalsPerSM;
  static const int MIN_HASH =  0; // always 0 ...
  static const int MAX_HASH =  2*MAX_IPHI*MAX_IETA-1;

  enum { kSizeForDenseIndexing = MAX_HASH + 1 } ;
  

  // function modes for (int, int) constructor
  static const int ETAPHIMODE = 0;
  static const int SMCRYSTALMODE = 1;
};

std::ostream& operator<<(std::ostream& s,const EBDetId& id);


#endif
