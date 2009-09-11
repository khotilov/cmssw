#ifndef MuonReco_Muon_h
#define MuonReco_Muon_h
/** \class reco::Muon Muon.h DataFormats/MuonReco/interface/Muon.h
 *  
 * A reconstructed Muon.
 * contains reference to three fits:
 *  - tracker alone
 *  - muon detector alone
 *  - combined muon plus tracker
 *
 * \author Luca Lista, Claudio Campagnari, Dmytro Kovalskyi, Jake Ribnik
 *
 * \version $Id: Muon.h,v 1.53 2009/09/09 20:55:46 aeverett Exp $
 *
 */
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"

namespace reco {
 
  class Muon : public RecoCandidate {
  public:
    Muon();
    /// constructor from values
    Muon(  Charge, const LorentzVector &, const Point & = Point( 0, 0, 0 ) );
    /// create a clone
    Muon * clone() const;
    
    ///
    /// ====================== TRACK BLOCK ===========================
    ///
    /// reference to Track reconstructed in the tracker only
    using reco::RecoCandidate::track;
    virtual TrackRef innerTrack() const { return innerTrack_; }
    virtual TrackRef track() const { return innerTrack(); }
    /// reference to Track reconstructed in the muon detector only
    virtual TrackRef outerTrack() const { return outerTrack_; }
    virtual TrackRef standAloneMuon() const { return outerTrack(); }
    /// reference to Track reconstructed in both tracked and muon detector
    virtual TrackRef globalTrack() const { return globalTrack_; }
    virtual TrackRef combinedMuon() const { return globalTrack(); }
    /// set reference to Track
    virtual void setInnerTrack( const TrackRef & t );
    virtual void setTrack( const TrackRef & t );
    /// set reference to Track
    virtual void setOuterTrack( const TrackRef & t );
    virtual void setStandAlone( const TrackRef & t );
    /// set reference to Track
    virtual void setGlobalTrack( const TrackRef & t );
    virtual void setCombined( const TrackRef & t );
    
    ///
    /// ====================== ENERGY BLOCK ===========================
    ///
    /// energy deposition
    bool isEnergyValid() const { return energyValid_; }
    /// get energy deposition information
    MuonEnergy calEnergy() const { return calEnergy_; }
    /// set energy deposition information
    void setCalEnergy( const MuonEnergy& calEnergy ) { calEnergy_ = calEnergy; energyValid_ = true; }
    
    ///
    /// ====================== Quality BLOCK ===========================
    ///
    /// energy deposition
    bool isQualityValid() const { return qualityValid_; }
    /// get energy deposition information
    MuonQuality combinedQuality() const { return combinedQuality_; }
    /// set energy deposition information
    void setCombinedQuality( const MuonQuality& combinedQuality ) { combinedQuality_ = combinedQuality; qualityValid_ = true; }

    ///
    /// ====================== TIMING BLOCK ===========================
    ///
    /// timing information
    bool isTimeValid() const { return (time_.nDof>0); }
    /// get timing information
    MuonTime time() const { return time_; }
    /// set timing information
    void setTime( const MuonTime& time ) { time_ = time; }
     
    ///
    /// ====================== MUON MATCH BLOCK ===========================
    ///
    bool isMatchesValid() const { return matchesValid_; }
    /// get muon matching information
    std::vector<MuonChamberMatch>& matches() { return muMatches_;}
    const std::vector<MuonChamberMatch>& matches() const { return muMatches_;	}
    /// set muon matching information
    void setMatches( const std::vector<MuonChamberMatch>& matches ) { muMatches_ = matches; matchesValid_ = true; }
     
    ///
    /// ====================== MUON COMPATIBILITY BLOCK ===========================
    ///
    /// Relative likelihood based on ECAL, HCAL, HO energy defined as
    /// L_muon/(L_muon+L_not_muon)
    float caloCompatibility() const { return caloCompatibility_; }
    void  setCaloCompatibility(float input){ caloCompatibility_ = input; }
    bool  isCaloCompatibilityValid() const { return caloCompatibility_>=0; } 
    
    ///
    /// ====================== ISOLATION BLOCK ===========================
    ///
    /// Summary of muon isolation information 
    const MuonIsolation& isolationR03() const { return isolationR03_; }
    const MuonIsolation& isolationR05() const { return isolationR05_; }
    void setIsolation( const MuonIsolation& isoR03, const MuonIsolation& isoR05 );
    bool isIsolationValid() const { return isolationValid_; }

    /// define arbitration schemes
    enum ArbitrationType { NoArbitration, SegmentArbitration, SegmentAndTrackArbitration };
    
    ///
    /// ====================== USEFUL METHODs ===========================
    ///
    /// number of chambers
    int numberOfChambers() const { return muMatches_.size(); }
    /// get number of chambers with matched segments
    int numberOfMatches( ArbitrationType type = SegmentArbitration ) const;
    /// get bit map of stations with matched segments
    /// bits 0-1-2-3 = DT stations 1-2-3-4
    /// bits 4-5-6-7 = CSC stations 1-2-3-4
    unsigned int stationMask( ArbitrationType type = SegmentArbitration ) const;
    /// get bit map of stations with tracks within
    /// given distance (in cm) of chamber edges 
    /// bit assignments are same as above
    unsigned int stationGapMaskDistance( float distanceCut = 10. ) const;
    /// same as above for given number of sigmas
    unsigned int stationGapMaskPull( float sigmaCut = 3. ) const;
     
    /// muon type - type of the algorithm that reconstructed this muon
    /// multiple algorithms can reconstruct the same muon
    static const unsigned int GlobalMuon     =  1<<1;
    static const unsigned int TrackerMuon    =  1<<2;
    static const unsigned int StandAloneMuon =  1<<3;
    static const unsigned int CaloMuon =  1<<4;
    void setType( unsigned int type ) { type_ = type; }
    unsigned int type() const { return type_; }
    // override of method in base class reco::Candidate
    bool isMuon() const { return true; }
    bool isGlobalMuon()     const { return type_ & GlobalMuon; }
    bool isTrackerMuon()    const { return type_ & TrackerMuon; }
    bool isStandAloneMuon() const { return type_ & StandAloneMuon; }
    bool isCaloMuon() const { return type_ & CaloMuon; }
    
  private:
    /// check overlap with another candidate
    virtual bool overlap( const Candidate & ) const;
    /// reference to Track reconstructed in the tracker only
    TrackRef innerTrack_;
    /// reference to Track reconstructed in the muon detector only
    TrackRef outerTrack_;
    /// reference to Track reconstructed in both tracked and muon detector
    TrackRef globalTrack_;
    /// energy deposition 
    MuonEnergy calEnergy_;
    /// quality block
    MuonQuality combinedQuality_;
    /// Information on matching between tracks and segments
    std::vector<MuonChamberMatch> muMatches_;
    /// timing
    MuonTime time_;
     
    bool energyValid_;
    bool matchesValid_;
    bool isolationValid_;
    bool qualityValid_;
    /// muon hypothesis compatibility with observer calorimeter energy
    float caloCompatibility_;
    /// Isolation information for two cones with dR=0.3 and dR=0.5
    MuonIsolation isolationR03_;
    MuonIsolation isolationR05_;
    /// muon type mask
    unsigned int type_;

    // FixMe: Still missing trigger information

    /// get vector of muon chambers for given station and detector
    const std::vector<const MuonChamberMatch*> chambers( int station, int muonSubdetId ) const;
    /// get pointers to best segment and corresponding chamber in vector of chambers
    std::pair<const MuonChamberMatch*,const MuonSegmentMatch*> pair( const std::vector<const MuonChamberMatch*> &,
									ArbitrationType type = SegmentArbitration ) const;
     
   public:
     /// get number of segments
     int numberOfSegments( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     /// get deltas between (best) segment and track
     /// If no chamber or no segment returns 999999
     float dX       ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float dY       ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float dDxDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float dDyDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float pullX    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration, bool includeSegmentError = false ) const;
     float pullY    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration, bool includeSegmentError = false ) const;
     float pullDxDz ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration, bool includeSegmentError = false ) const;
     float pullDyDz ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration, bool includeSegmentError = false ) const;
     /// get (best) segment information
     /// If no segment returns 999999
     float segmentX       ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float segmentY       ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float segmentDxDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float segmentDyDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float segmentXErr    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float segmentYErr    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float segmentDxDzErr ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float segmentDyDzErr ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     /// get track information in chamber that contains (best) segment
     /// If no segment, get track information in chamber that has the most negative distance between the track
     /// and the nearest chamber edge (the chamber with the deepest track)
     /// If no chamber returns 999999
     float trackEdgeX   ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackEdgeY   ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackX       ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackY       ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackDxDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackDyDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackXErr    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackYErr    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackDxDzErr ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackDyDzErr ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackDist    ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     float trackDistErr ( int station, int muonSubdetId, ArbitrationType type = SegmentArbitration ) const;
     
     float t0(int n=0) {
	int i = 0;
	for( std::vector<MuonChamberMatch>::const_iterator chamber = muMatches_.begin();
	     chamber != muMatches_.end(); ++chamber )
	  for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->segmentMatches.begin();
		segment != chamber->segmentMatches.end(); ++segment )
	    {
	       if (i==n) return segment->t0;
	       ++i;
	    }
	return 0;
     }
  };

}


#endif
