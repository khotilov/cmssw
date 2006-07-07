#ifndef MuonDetId_CSCDetId_h
#define MuonDetId_CSCDetId_h

/** \class CSCDetId
 * Identifier class for hierarchy of Endcap Muon detector components.
 *
 * Ported from MuEndDetectorId but now derived from DetId and updated accordingly.
 *
 * Allows access to hardware integer labels of the subcomponents
 * of the Muon Endcap CSC detector system.
 *
 * The STATIC member functions can be used to translate back and
 * forth between a MuEndLayer 'rawId' and the %set of subdetector labels.
 *
 * \warning EVERY LABEL COUNTS FROM ONE NOT ZERO.
 *
 */

#include <iosfwd>
#include <DataFormats/DetId/interface/DetId.h>
#include <DataFormats/MuonDetId/interface/MuonSubdetId.h>

class CSCDetId;

std::ostream& operator<<( std::ostream& os, const CSCDetId& id );

class CSCDetId : public DetId {

public:

  /// Default constructor; fills the common part in the base
  /// and leaves 0 in all other fields
  CSCDetId();

  /// Construct from a packed id. It is required that the Detector part of
  /// id is Muon and the SubDet part is CSC, otherwise an exception is thrown.
  explicit CSCDetId(uint32_t id);


  /// Construct from fully qualified identifier.
  /// Input values are required to be within legal ranges, otherwise an
  /// exception is thrown.
  CSCDetId( int iendcap, int istation, 
	    int iring, int ichamber, 
	    int ilayer );

  /** Copy ctor.
   */
  CSCDetId( const CSCDetId& id )
     : DetId(id.id_) { }  

  /** Chamber CSCDetId from a Layer CSCDetId
   */
  CSCDetId chamberId() const {
    // build chamber id by removing layer bits
    return CSCDetId( id_ - layer() ) ; }

  /**
   * Return Layer label.
   *
   */
  int layer() const {
    return (id_ & MASK_LAYER); } 

  /**
   * Return Chamber label.
   *
   */
   int chamber() const {
     return (  (id_>>START_CHAMBER) & MASK_CHAMBER ); }

  /**
   * Return %Ring label.
   *
   */
   int ring() const {
     if (((id_>>START_STATION) & MASK_STATION) == 1)
       return (  detIdToInt((id_>>START_RING) & MASK_RING )); 
     else
       return (((id_>>START_RING) & MASK_RING )); 
   }

  /**
   * Return Station label.
   *
   */
   int station() const {
     return (  (id_>>START_STATION) & MASK_STATION ); }

  /**
   * Return Endcap label.
   *
   */
   int endcap() const {
     return (  (id_>>START_ENDCAP) & MASK_ENDCAP ); }


  // static methods
  // Used when we need information about subdetector labels.

  /** 
   * Returns the unique integer 'rawId' which labels each CSC layer.
   *
   * The arguments are the integer labels for, respectively,  <br>
   * endcap, station, ring, chamber, layer.
   *
   * \warning The input int args are expected to be .ge. 1 and there
   * is no sanity-checking for their upper limits.
   *
   * \warning The returned integers are not necessarily consecutive, i.e.
   * there are gaps. This is to permit computational efficiency
   * starting from the component ids.
   *
   */

  // Tim dislikes the necessity of this ugly code - magic numbers included
  // Thanks a lot, CMSSW

   static int rawIdMaker( int iendcap, int istation, int iring, 
               int ichamber, int ilayer ) {
     return ((DetId::Muon&0xF)<<28)|((MuonSubdetId::CSC&0x7)<<25)|
               init(iendcap, istation, iring, ichamber, ilayer) ; }

   /**
    * Return Layer label for supplied CSCDetId index.
    *
    */
   static int layer( int index ) {
     return (index & MASK_LAYER); }

   /**
    * Return Chamber label for supplied CSCDetId index.
    *
    */
   static int chamber( int index ) {
     return (  (index>>START_CHAMBER) & MASK_CHAMBER ); }

   /**
    * Return Ring label for supplied CSCDetId index.
    *
    */
   static int ring( int index ) {
     if (((index>>START_STATION) & MASK_STATION) == 1)
       return (  detIdToInt((index>>START_RING) & MASK_RING )); 
     else
       return (( index>>START_RING) & MASK_RING ); 
   }

   /**
    * Return Station label for supplied CSCDetId index.
    *
    */
   static int station( int index ) {
     return (  (index>>START_STATION) & MASK_STATION ); }

   /**
    * Return Endcap label for supplied CSCDetId index.
    *
    */
   static int endcap( int index ) {
     return (  (index>>START_ENDCAP) & MASK_ENDCAP ); }

   /**
    * Return trigger-level sector id for an Endcap Muon chamber.
    *
    * This method encapsulates the information about which chambers
    * are in which sectors, and may need updating according to
    * hardware changes, or software chamber indexing.
    *
    * Station 1 has 3 rings of 10-degree chambers. <br>
    * Stations 2, 3, 4 have an inner ring of 20-degree chambers
    * and an outer ring of 10-degree chambers. <br>
    *
    * Sectors are 60 degree slices of a station, covering both rings. <br>
    * For Station 1, there are subsectors of 30 degrees: 9 10-degree
    * chambers (3 each from ME1/1, ME1/2, ME1/3.) <br>
    * 
    * The first sector starts at phi = 15 degrees so it matches Barrel Muon sectors.
    * We count from one not zero.
    *
    */
   int triggerSector() const;

   /**
    * Return trigger-level CSC id  within a sector for an Endcap Muon chamber.
    *
    * This id is an index within a sector such that the 3 inner ring chambers 
    * (20 degrees each) are 1, 2, 3 (increasing counterclockwise) and the 6 outer ring 
    * chambers (10 degrees each) are 4, 5, 6, 7, 8, 9 (again increasing counter-clockwise.) 
    *
    * This method knows which chambers are part of which sector and returns
    * the chamber label/index/identifier accordingly.
    * Beware that this information is liable to change according to hardware
    * and software changes.
    *
    */
   int triggerCscId() const;

   /**
    * Lower and upper counts for the subdetector hierarchy
    */
   static int minEndcapId()  { return MIN_ENDCAP; }
   static int maxEndcapId()  { return MAX_ENDCAP; }
   static int minStationId() { return MIN_STATION; }
   static int maxStationId() { return MAX_STATION; }
   static int minRingId()    { return MIN_RING; }
   static int maxRingId()    { return MAX_RING; }
   static int minChamberId() { return MIN_CHAMBER; }
   static int maxChamberId() { return MAX_CHAMBER; }
   static int minLayerId()   { return MIN_LAYER; }
   static int maxLayerId()   { return MAX_LAYER; }

private:
 
  /**
   * Method for initialization within ctors.
   *
   */
  static uint32_t init( int iendcap, int istation, 
			int iring, int ichamber, int ilayer ) {
    
    if (istation == 1)
      iring = intToDetId(iring);

     return
         (ilayer   & MASK_LAYER)                      |
       ( (ichamber & MASK_CHAMBER) << START_CHAMBER ) |
       ( (iring    & MASK_RING)    << START_RING )    |
       ( (istation & MASK_STATION) << START_STATION ) | 
       ( (iendcap  & MASK_ENDCAP)  << START_ENDCAP ) ; }

  static int intToDetId(int iring) {
    int i = (iring+1)%4;
    if (i == 0)
      i = 4;
    return i;
  }
  static int detIdToInt(int iring) {
    int i = (iring-1)%4;
    if (i == 0)
      i = 4;
    return i;
  }
 
  // The following define the bit-packing implementation...

  // The maximum numbers of various parts
  enum eMaxNum{ MAX_ENDCAP=2, MAX_STATION=4, MAX_RING=4, MAX_CHAMBER=36, MAX_LAYER=6 };
  // We count from 1 
  enum eMinNum{ MIN_ENDCAP=1, MIN_STATION=1, MIN_RING=1, MIN_CHAMBER=1, MIN_LAYER=1 };

  // BITS_det is no. of binary bits required to label 'det' but allow 0 as a wild-card character
  enum eNumBitDet{ BITS_ENDCAP=2, BITS_STATION=3,  BITS_RING=3, BITS_CHAMBER=7, BITS_LAYER=3 };

  // MASK_det is binary bits set to pick off the bits for 'det' (defined as octal)
  enum eMaskBitDet{ MASK_ENDCAP=03, MASK_STATION=07, MASK_RING=07, MASK_CHAMBER=0177, MASK_LAYER=07 };

  // START_det is bit position (counting from zero) at which bits for 'det' start in 'rawId' word
  enum eStartBitDet{ START_CHAMBER=BITS_LAYER, START_RING=START_CHAMBER+BITS_CHAMBER,
          START_STATION=START_RING+BITS_RING, START_ENDCAP=START_STATION+BITS_STATION };
};

#endif


