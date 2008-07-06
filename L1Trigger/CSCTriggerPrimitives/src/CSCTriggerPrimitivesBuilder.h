#ifndef CSCTriggerPrimitives_CSCTriggerPrimitivesBuilder_h
#define CSCTriggerPrimitives_CSCTriggerPrimitivesBuilder_h

/** \class CSCTriggerPrimitivesBuilder
 *
 * Algorithm to build anode, cathode, and correlated LCTs from wire and
 * comparator digis in endcap muon CSCs by implementing a 'build' function
 * required by CSCTriggerPrimitivesProducer.
 *
 * Configured via the Producer's ParameterSet.
 *
 * \author Slava Valuev, UCLA.
 *
 * $Date: 2007/04/18 16:08:55 $
 * $Revision: 1.4 $
 *
 */

#include <DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCWireDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

class CSCL1TPParameters;
class CSCMotherboard;
class CSCMuonPortCard;

class CSCTriggerPrimitivesBuilder
{
 public:

  /** Configure the algorithm via constructor.
   *  Receives ParameterSet percolated down from EDProducer which owns this
   *  Builder.
   */
  explicit CSCTriggerPrimitivesBuilder(const edm::ParameterSet&);

  ~CSCTriggerPrimitivesBuilder();

  /** Sets configuration parameters obtained via EventSetup mechanism. */
  void setConfigParameters(const CSCL1TPParameters* conf);

  /** Build anode, cathode, and correlated LCTs in each chamber and fill
   *  them into output collections.  Select up to three best correlated LCTs
   *  in each (sub)sector and put them into an output collection as well. */
  void build(const CSCWireDigiCollection* wiredc,
	     const CSCComparatorDigiCollection* compdc,
	     CSCALCTDigiCollection& oc_alct, CSCCLCTDigiCollection& oc_clct,
	     CSCCorrelatedLCTDigiCollection& oc_lct,
	     CSCCorrelatedLCTDigiCollection& oc_sorted_lct);

  /** Max values of trigger labels for all CSCs; used to construct TMB
   *  processors. */
  enum trig_cscs {MAX_ENDCAPS = 2, MAX_STATIONS = 4, MAX_SECTORS = 6,
		  MAX_SUBSECTORS = 2, MAX_CHAMBERS = 9};
 private:

  /** Min and max allowed values for various CSC elements, defined in
   *  CSCDetId and CSCTriggerNumbering classes. */
  static const int min_endcap;    // endcaps
  static const int max_endcap;
  static const int min_station;   // stations per endcap
  static const int max_station;
  static const int min_sector;    // trigger sectors per station
  static const int max_sector;
  static const int min_subsector; // trigger subsectors per sector
  static const int max_subsector;
  static const int min_chamber;   // chambers per trigger subsector
  static const int max_chamber;

  int m_minBX, m_maxBX; // min and max BX to sort.

  /** Pointers to TMB processors for all possible chambers. */
  CSCMotherboard*
    tmb_[MAX_ENDCAPS][MAX_STATIONS][MAX_SECTORS][MAX_SUBSECTORS][MAX_CHAMBERS];

  /** Pointer to MPC processor. */
  CSCMuonPortCard* m_muonportcard;
};

#endif
