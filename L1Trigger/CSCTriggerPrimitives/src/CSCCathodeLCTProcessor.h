#ifndef CSCTriggerPrimitives_CSCCathodeLCTProcessor_h
#define CSCTriggerPrimitives_CSCCathodeLCTProcessor_h

/** \class CSCCathodeLCTProcessor
 *
 * This class simulates the functionality of the cathode LCT card. It is run by
 * the MotherBoard and returns up to two CathodeLCTs.  It can be run either in
 * a test mode, where it is passed an array of comparator times and comparator
 * values, or in normal mode where it determines the time and comparator
 * information from the comparator digis.
 *
 * The CathodeLCTs come in distrip and halfstrip flavors; they are sorted
 * (from best to worst) as follows: 6/6H, 5/6H, 6/6D, 4/6H, 5/6D, 4/6D.
 *
 * \date May 2001  Removed the card boundaries.  Changed the Pretrigger to
 * emulate the hardware electronic logic.  Also changed the keylayer to be
 * the 4th layer in a chamber instead of the 3rd layer from the interaction
 * region. The code is a more realistic simulation of hardware LCT logic now.
 * -Jason Mumford.
 *
 * \author Benn Tannenbaum  UCLA 13 July 1999 benn@physics.ucla.edu
 * Numerous later improvements by Jason Mumford and Slava Valuev (see cvs
 * in ORCA).
 * Porting from ORCA by S. Valuev (Slava.Valuev@cern.ch), May 2006.
 *
 */

#include <vector>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCLCTDigi.h>
#include <L1Trigger/CSCCommonTrigger/interface/CSCConstants.h>

class CSCCathodeLCTProcessor
{
 public:
  /** Normal constructor. */
#ifndef L1CSC_STANDALONE
  CSCCathodeLCTProcessor(unsigned endcap, unsigned station, unsigned sector,
			 unsigned subsector, unsigned chamber,
			 const edm::ParameterSet& conf);
#endif

  /** Default constructor. Used for testing. */
  CSCCathodeLCTProcessor();

  /** Runs the LCT processor code. Called in normal running -- gets info from
      a collection of comparator digis. */
  std::vector<CSCCLCTDigi> run(const CSCComparatorDigiCollection* compdc);

  /** Called in test mode and by the run() function; does the actual LCT
      finding. */
  void run(int triad[CSCConstants::NUM_LAYERS][CSCConstants::MAX_NUM_STRIPS],
	   int time[CSCConstants::NUM_LAYERS][CSCConstants::MAX_NUM_STRIPS],
	   int digiNum[CSCConstants::NUM_LAYERS][CSCConstants::MAX_NUM_STRIPS]);

  /** Clears the LCT containers. */
  void clear();

  /** Turns on the debug flag for this class. */
  static void setDebug() {debug = true;}

  /** Turns off the debug flag for this class (default). */
  static void setNoDebug() {debug = false;}

  /** Best LCT in this chamber, as found by the processor. */
  CSCCLCTDigi bestCLCT;

  /** Second best LCT in this chamber, as found by the processor. */
  CSCCLCTDigi secondCLCT;

#ifndef L1CSC_STANDALONE
  /** Access routine to comparator digis. */
  void getDigis(const CSCComparatorDigiCollection* compdc);
#endif

  /** This version used for L1CSC_STANDALONE mode. */
  void getDigis(const std::vector<std::vector<CSCComparatorDigi> > digis);

  /** Returns vector of found CLCTs, if any. */
  std::vector<CSCCLCTDigi> getCLCTs();

  /** Access to times on halfstrips on any layer. */
  std::vector<int> halfStripHits(const int layer) const;

  /** Access to time on single halfstrip on any layer. */
  int halfStripHit(const int layer, const int strip) const;

  /** Access to times on distrips on any layer. */
  std::vector<int> diStripHits(const int layer) const;

  /** Access to time on single distrip on any layer. */
  int diStripHit(const int layer, const int strip) const;

  static void distripStagger(int stag_triad[CSCConstants::MAX_NUM_STRIPS],
			     int stag_time[CSCConstants::MAX_NUM_STRIPS],
			     int stag_digi[CSCConstants::MAX_NUM_STRIPS],
			     int i_distrip);

  /** Pre-defined patterns. */
  enum {NUM_PATTERN_STRIPS = 36};
  static const int pre_hit_pattern[2][NUM_PATTERN_STRIPS];
  static const int pattern[CSCConstants::NUM_CLCT_PATTERNS][NUM_PATTERN_STRIPS+1];

  /** Maximum number of cathode front-end boards (move to CSCConstants?). */
  enum {NUM_CFEBS = 5};

  // we use these next ones to address the various bits inside the array that's
  // used to make the cathode LCTs.
  enum CLCT_INDICES {CLCT_PATTERN, CLCT_BEND, CLCT_STRIP, CLCT_BX,
		     CLCT_STRIP_TYPE, CLCT_QUALITY};

 private:
  const unsigned theEndcap;
  const unsigned theStation;
  const unsigned theSector;
  const unsigned theSubsector;
  const unsigned theTrigChamber;

  int numStrips;
  int stagger[CSCConstants::NUM_LAYERS];

  std::vector<CSCComparatorDigi> digiV[CSCConstants::NUM_LAYERS];

  std::vector<int> theHalfStripHits[CSCConstants::NUM_LAYERS];
  std::vector<int> theDiStripHits[CSCConstants::NUM_LAYERS];

  static int test_iteration;
  static bool debug;

  /** Set times on all layers for distrips and halfstrips. */
  void saveAllHits(const int distrip[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS], 
		   const int halfstrip[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS]);

  /** Set times on halfstrips on any layer. */
  void setHalfStripHits(const int layer, const std::vector<int>& hStripHits);

  /** Set time on single halfstrip on any layer. */
  void setHalfStripHit(const int layer, const int strip, const int hit);

  /** Set times on distrips on any layer. */
  void setDiStripHits(const int layer, const std::vector<int>& dStripHits);

  /** Set time on single distrip on any layer. */
  void setDiStripHit(const int layer, const int strip, const int hit);

  //----------------------- Default ORCA Fcns ---------------------------------
  bool hitIsGood(int hitTime, int BX);
  std::vector<CSCCLCTDigi> findLCTs(const int strip[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS],
				    int width, int numStrips);
  bool preTrigger(const int strip[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS],
		  int nStrips, int& first_bx);
  void getKeyStripData(const int strip[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS],
		       int keystrip_data[CSCConstants::NUM_HALF_STRIPS][6],
		       int nStrips, int first_bx, int& best_strip,
		       int stripType);
  void getPattern(int pattern_num, int strip_value[NUM_PATTERN_STRIPS],
		  int bx_time, int &quality, int &bend);

  //----------------------- Test Beam Fcns below ----------------------------
  std::vector<CSCCLCTDigi> findLCTs(const int halfstrip[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS],
				    const int distrip[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS]);
  bool preTrigger(const int strip[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS],
		  unsigned long int pulse[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS], 
		  const int nStrips, int& first_bx);
  bool preTrigLookUp(const unsigned long int pulse[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS],
		     const int nStrips, const int bx_time);
  void latchLCTs(const unsigned long int pulse[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS],
		 int keyStrip[NUM_CFEBS], int nhits[NUM_CFEBS],
		 const int nStrips, const int bx_time);
  void priorityEncode(const int h_keyStrip[NUM_CFEBS],
		      const int h_nhits[NUM_CFEBS],
		      const int d_keyStrip[NUM_CFEBS],
		      const int d_nhits[NUM_CFEBS], int keystrip_data[2][6]);
  void getKeyStripData(const unsigned long int h_pulse[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS],
		       const unsigned long int d_pulse[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS],
		       int keystrip_data[2][6],
		       const int latch_bx, const int drift_delay);
  void getPattern(int pattern_num, const int strip_value[NUM_PATTERN_STRIPS],
		  int& quality, int& bend);
  //-------------------------------------------------------------------------

  void testDistripStagger();
  void testLCTs();
  void printPatterns();
  void testPatterns();
  int findNumLayersHit(int stripsHit[CSCConstants::NUM_LAYERS][CSCConstants::NUM_HALF_STRIPS]);
};

#endif
