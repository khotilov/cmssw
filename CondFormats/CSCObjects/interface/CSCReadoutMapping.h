#ifndef CondFormats_CSCReadoutMapping_h
#define CondFormats_CSCReadoutMapping_h

/** 
 * \class CSCReadoutMapping
 * \author Tim Cox
 * Abstract class to define mapping between CSC readout hardware ids and other labels.
 *
 * Defines the ids and labels in the mapping and supplies tramslation interface.
 * A derived class must define how hardware labels map to a unique integer.
 * A derived, concrete, class must define from where the mapping information comes.
 */

//@@ FIXME This whole design would better suit a Factory/Builder pattern

#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <vector>
#include <map>

class CSCReadoutMapping {
 public:

  /// Default constructor
  CSCReadoutMapping();

  /// Destructor
  virtual ~CSCReadoutMapping();

  /**
   * Instead of a set of vectors of int use one vector of a set of ints
   */  
  struct CSCLabel{
    CSCLabel(){}
    CSCLabel( int endcap, int station, int ring, int chamber,  
	      int vmecrate, int dmb, int tmb, int tsector, int cscid )
      : endcap_( endcap ), station_( station ), ring_( ring ), chamber_( chamber ),
         vmecrate_( vmecrate ), dmb_( dmb ), tmb_( tmb ), 
	 tsector_( tsector ), cscid_( cscid ) {}
    ~CSCLabel(){}

    int endcap_;
    int station_;
    int ring_;
    int chamber_;
    int vmecrate_;
    int dmb_;
    int tmb_;
    int tsector_;
    int cscid_;
  };

   /**
    * Return CSCDetId for layer corresponding to readout ids vme, tmb, and dmb for given endcap
    * and layer no. 1-6, or for chamber if no layer no. supplied.
    * Args: endcap = 1 (+z), 2 (-z), station, vme crate number, dmb slot number, tmb slot number, layer#
    */
    // layer at end so it can have default arg
    CSCDetId detId( int endcap, int station, int vmecrate, int dmb, int tmb, int layer = 0 ) const;

   /** 
    * Return chamber label corresponding to readout ids vme, tmb and dmb for given endcap
    *  endcap = 1 (+z), 2 (-z), vme crate number, tmb slot number, dmb slot number
    */
    int chamber( int endcap, int station, int vmecrate, int dmb, int tmb ) const;

   /** 
    * Fill mapping store
    */
    virtual void fill( void ) = 0;

   /**
    * Add one record of info to mapping
    */
    void addRecord( int endcap, int station, int ring, int chamber, 
		    int vmecrate, int dmb, int tmb, int tsector, int cscid ); 

    /**
     * Set debug printout flag
     */
    void setDebugV( bool dbg ) { debugV_ = dbg; }

    /**
     * Status of debug printout flag
     */
    bool debugV( void ) const { return debugV_; }

    /**
     * Return class name
     */
    const std::string& myName( void ) const { return myName_; }

 private: 

    /**
     * Build a unique integer out of the readout electronics labels.
     *
     * In general this must depend on endcap and station, as well as
     * vme crate number and dmb slot number. And possibly tmb slot?
     */
    virtual int hwId( int endcap, int station, int vme, int dmb, int tmb ) const = 0;

    /**
     * Build a unique integer out of chamber labels.
     *
     * We'll probably use rawId of CSCDetId... You know it makes sense!
     */
    int swId( int endcap, int station, int ring, int chamber) const;

    std::string myName_;
    bool debugV_;
    std::vector< CSCLabel > mapping_;
    std::map< int, int > hw2sw_;
};

#endif
