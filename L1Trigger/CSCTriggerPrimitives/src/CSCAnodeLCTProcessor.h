#ifndef CSCTriggerPrimitives_CSCAnodeLCTProcessor_h
#define CSCTriggerPrimitives_CSCAnodeLCTProcessor_h

/** \class CSCAnodeLCTProcessor
 *
 * This class simulates the functionality of the anode LCT card. It is run by
 * the MotherBoard and returns up to two AnodeLCTs.  It can be run either in a
 * test mode, where it is passed an array of wire times, or in normal mode
 * where it determines the wire times from the wire digis.
 *
 * \author Benn Tannenbaum  benn@physics.ucla.edu 13 July 1999
 * Numerous later improvements by Jason Mumford and Slava Valuev (see cvs
 * in ORCA).
 * Porting from ORCA by S. Valuev (Slava.Valuev@cern.ch), May 2006.
 *
 * $Date: 2007/02/19 14:59:46 $
 * $Revision: 1.8 $
 *
 */

#include <vector>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/CSCDigi/interface/CSCWireDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCALCTDigi.h>
#include <CondFormats/L1TObjects/interface/L1CSCTPParameters.h>
#include <L1Trigger/CSCCommonTrigger/interface/CSCConstants.h>

class CSCAnodeLCTProcessor
{
 public:
  /** Normal constructor. */
  CSCAnodeLCTProcessor(unsigned endcap, unsigned station, unsigned sector,
		       unsigned subsector, unsigned chamber,
		       const edm::ParameterSet& conf);

  /** Default constructor. Used for testing. */
  CSCAnodeLCTProcessor();

  /** Sets configuration parameters obtained via EventSetup mechanism. */
  void setConfigParameters(const L1CSCTPParameters* conf);

  /** Clears the LCT containers. */
  void clear();

  /** Runs the LCT processor code. Called in normal running -- gets info from
      a collection of wire digis. */
  std::vector<CSCALCTDigi> run(const CSCWireDigiCollection* wiredc);

  /** Runs the LCT processor code. Called in normal running or in testing
      mode. */
  void run(const int wire[CSCConstants::NUM_LAYERS][CSCConstants::MAX_NUM_WIRES]);

  /** Access routine to wire digis. */
  bool getDigis(const CSCWireDigiCollection* wiredc);

  /** Best LCT in this chamber, as found by the processor. */
  CSCALCTDigi bestALCT;

  /** Second best LCT in this chamber, as found by the processor. */
  CSCALCTDigi secondALCT;

  /** Returns vector of found ALCTs, if any. */
  std::vector<CSCALCTDigi> getALCTs();

  /** Access to times on wires on any layer. */
  std::vector<int> wireHits(const int layer) const;

  /** Access to time on single wire on any layer. */
  int wireHit(const int layer, const int wire) const;

  /** Pre-defined patterns. */
  enum {NUM_PATTERN_WIRES = 14};
  static const int pattern_envelope[CSCConstants::NUM_ALCT_PATTERNS][NUM_PATTERN_WIRES];
  static const int pattern_mask[CSCConstants::NUM_ALCT_PATTERNS][NUM_PATTERN_WIRES];

 private:
  /** Verbosity level: 0: no print (default).
   *                   1: print only ALCTs found.
   *                   2: info at every step of the algorithm.
   *                   3: add special-purpose prints. */
  int infoV;

  /** Chamber id (trigger-type labels). */
  const unsigned theEndcap;
  const unsigned theStation;
  const unsigned theSector;
  const unsigned theSubsector;
  const unsigned theTrigChamber;

  int numWireGroups;
  int MESelection;

  int first_bx[CSCConstants::MAX_NUM_WIRES];
  int quality[CSCConstants::MAX_NUM_WIRES][3];
  std::vector<CSCWireDigi> digiV[CSCConstants::NUM_LAYERS];
  unsigned int pulse[CSCConstants::NUM_LAYERS][CSCConstants::MAX_NUM_WIRES];

  std::vector<int> theWireHits[CSCConstants::NUM_LAYERS];

  /** Configuration parameters. */
  unsigned int fifo_tbins, fifo_pretrig, bx_width, drift_delay;
  unsigned int nph_thresh, nph_pattern;
  unsigned int trig_mode, alct_amode, l1a_window;

  /** Set default values for configuration parameters. */
  void setDefaultConfigParameters();

  /** Make sure that the parameter values are within the allowed range. */
  void checkConfigParameters() const;

  /** Clears the quality for a given wire and pattern if it is a ghost. */
  void clear(const int wire, const int pattern);

  /** ALCT algorithm methods. */
  void readWireDigis(int wire[CSCConstants::NUM_LAYERS][CSCConstants::MAX_NUM_WIRES]);
  bool pulseExtension(const int wire[CSCConstants::NUM_LAYERS][CSCConstants::MAX_NUM_WIRES]);
  bool preTrigger(const int key_wire);
  bool patternDetection(const int key_wire);
  void ghostCancellationLogic();
  void lctSearch();
  void trigMode(const int key_wire);
  void alctAmode(const int key_wire);

  std::vector<CSCALCTDigi>
    bestTrackSelector(const std::vector<CSCALCTDigi>& all_alcts);
  bool isBetterALCT(const CSCALCTDigi& lhsALCT, const CSCALCTDigi& rhsALCT);

  /** Dump ALCT configuration parameters. */
  void dumpConfigParams() const;

  /** Dump digis on wire groups. */
  void dumpDigis(const int wire[CSCConstants::NUM_LAYERS][CSCConstants::MAX_NUM_WIRES]) const;

  /** Set times on all layers for all wires. */
  void saveAllHits(const int wires[CSCConstants::NUM_LAYERS][CSCConstants::MAX_NUM_WIRES]);

  /** Set times on wires on any layer. */
  void setWireHits(const int layer, const std::vector<int>& wireHits);

  /** Set time on single wire on any layer. */
  void setWireHit(const int layer, const int wire, const int hit);

  void showPatterns(const int key_wire);
};

#endif
