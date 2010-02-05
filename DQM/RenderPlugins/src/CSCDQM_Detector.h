/*
 * =====================================================================================
 *
 *       Filename:  CSCDQM_Detector.h
 *
 *    Description:  CSC detector functions.
 *
 *        Version:  1.0
 *        Created:  05/19/2008 10:52:21 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Valdas Rapsevicius (VR), Valdas.Rapsevicius@cern.ch
 *        Company:  CERN, CH
 *
 * =====================================================================================
 */

#ifndef CSCDQM_Detector_H
#define CSCDQM_Detector_H

#include <math.h>
#include <float.h>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>

#ifdef CSC_RENDER_PLUGIN
#include "CSCDQM_Utility.h"
#else
#include "DQM/CSCMonitorModule/interface/CSCDQM_Utility.h"
#endif

namespace cscdqm {

/**
 * Number of Detector Components.
 */

#define N_SIDES    2
#define N_STATIONS 4
#define N_RINGS    3
#define N_CHAMBERS 36
#define N_LAYERS   6
#define N_CFEBS    5
#define N_HVS      5

/** Size of the address (number of components) */
#define ADDR_SIZE  7

/** Number of addressing elements in detector */
#define N_ELEMENTS 9540
//(7740 + 1800)

/**
 * Partition function shortcuts
 */

#define PARTITION_INDEX(x,y)  (x * partitions_y + y)
#define PARTITION_STEP_X      (5.0 / partitions_x)
#define PARTITION_STEP_Y      ((2.0 * 3.14159) / partitions_y)

/**
 * @brief  Mask of the address which is used to switch on and off appropriate Address fields.
 */
struct AddressMask {
  bool side;
  bool station;
  bool ring;
  bool chamber;
  bool layer;
  bool cfeb;
  bool hv;
};

/**
 * @brief  Structure to store detector addresses of any granularity: from
 * whole detector to the single HV element.
 */
struct Address {

  unsigned int side;
  unsigned int station;
  unsigned int ring;
  unsigned int chamber;
  unsigned int layer;
  unsigned int cfeb;
  unsigned int hv;

  AddressMask mask;

  bool operator== (const Address& a) const {
    if (mask.side    == a.mask.side    && mask.side    == true && side    != a.side)    return false;
    if (mask.station == a.mask.station && mask.station == true && station != a.station) return false;
    if (mask.ring    == a.mask.ring    && mask.ring    == true && ring    != a.ring)    return false;
    if (mask.chamber == a.mask.chamber && mask.chamber == true && chamber != a.chamber) return false;
    if (mask.layer   == a.mask.layer   && mask.layer   == true && layer   != a.layer)   return false;
    if (mask.cfeb    == a.mask.cfeb    && mask.cfeb    == true && cfeb    != a.cfeb)    return false;
    if (mask.hv      == a.mask.hv      && mask.hv      == true && hv      != a.hv)      return false;
    return true;
  };

  const Address* operator= (const Address& a) {
    mask.side    = a.mask.side;
    side         = a.side;
    mask.station = a.mask.station;
    station      = a.station;
    mask.ring    = a.mask.ring;
    ring         = a.ring;
    mask.chamber = a.mask.chamber;
    chamber      = a.chamber;
    mask.layer   = a.mask.layer;
    layer        = a.layer;
    mask.cfeb    = a.mask.cfeb;
    cfeb         = a.cfeb;
    mask.hv      = a.mask.hv;
    hv           = a.hv;
    return this;
  };

};

/**
 * @brief  Area covered by Address in eta/phy space
 */
struct AddressBox {
  Address adr;
  float xmin;
  float xmax;
  float ymin;
  float ymax;
};


/** Map of partitions and partition covering adresses indexes type */
typedef std::map<unsigned int, std::vector<unsigned int> > PartitionMap;

/** Iterator type of PartitionMap */
typedef PartitionMap::iterator PartitionMapIterator;

/**
 * @class Detector
 * @brief Detector geometry and addressing related imformation and routines
 */
class Detector {

  public:

    Detector(unsigned int p_partitions_x = 0, unsigned int p_partitions_y = 0);

    bool NextAddress(unsigned int& i, const Address*& adr, const Address& mask) const;
    bool NextAddressBox(unsigned int& i, const AddressBox*& box, const Address& mask) const;
    //bool NextAddressBoxByPartition(unsigned int& i, unsigned int& px, unsigned int& py, const AddressBox*& box, const Address& mask, float xmin, float xmax, float ymin, float ymax);
    bool NextAddressBoxByPartition (unsigned int& i, unsigned int px, unsigned int py, AddressBox*& box);

    float Area(unsigned int station) const;
    float Area(const Address& adr) const;

    void PrintAddress(const Address& adr) const;
    const std::string AddressName(const Address& adr) const;
    bool AddressFromString(const std::string str_address, Address& adr) const;

    unsigned int NumberOfRings(unsigned int station) const;
    unsigned int NumberOfChambers(unsigned int station, unsigned int ring) const;
    unsigned int NumberOfChamberCFEBs(unsigned int station, unsigned int ring) const;
    unsigned int NumberOfChamberHVs(unsigned int station, unsigned int ring) const;

  private:

    float Eta(float r, float z) const;
    float EtaToX(float eta) const;
    float PhiToY(float phi) const;
    float Z(const int station, const int ring) const;
    float RMinHV(const int station, const int ring, const int n_hv) const;
    float RMaxHV(const int station, const int ring, const int n_hv) const;
    float PhiMinCFEB(const int station, const int ring, const int chamber, const int cfeb) const;
    float PhiMaxCFEB(const int station, const int ring, const int chamber, const int cfeb) const;

    /** Address boxes in epa/phi space */
    AddressBox boxes[N_ELEMENTS];

    /** Station areas precalculated */
    float station_area[N_STATIONS];

    /** Number of partitions in X axis */
    unsigned int partitions_x;

    /** Number of partitions in Y axis */
    unsigned int partitions_y;

    /** Map of partitions and list of it covering addresses indexes */
    PartitionMap partitions;

};

}

#endif
