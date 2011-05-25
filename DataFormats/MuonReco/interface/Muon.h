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
 * \version $Id: Muon.h,v 1.59 2011/01/29 00:32:44 jribnik Exp $
 *
 */
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
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



    ///set reference to PFCandidate
    ///
    /// ====================== PF BLOCK ===========================
    ///

    PFCandidateRef pfMuon() {return pfMuon_;}
    virtual void setPFMuon( const PFCandidateRef & muon ); 

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

    const MuonPFIsolation& pfIsolationR03() const { return pfIsolationR03_; }
    const MuonPFIsolation& pfIsolationR04() const { return pfIsolationR04_; }

    void setIsolation( const MuonIsolation& isoR03, const MuonIsolation& isoR05 );
    bool isIsolationValid() const { return isolationValid_; }

    void setPFIsolation( const MuonPFIsolation& isoR03, const MuonPFIsolation& isoR04 );
    bool isPFIsolationValid() const { return pfIsolationValid_; }


    /// define arbitration schemes
    enum ArbitrationType { NoArbitration, SegmentArbitration, SegmentAndTrackArbitration, SegmentAndTrackArbitrationCleaned };
    
    ///
    /// ====================== USEFUL METHODs ===========================
    ///
    /// number of chambers
    int numberOfChambers() const { return muMatches_.size(); }
    /// get number of chambers with matched segments
    int numberOfMatches( ArbitrationType type = SegmentAndTrackArbitration ) const;
    /// get number of stations with matched segments
    /// just adds the bits returned by stationMask
    int numberOfMatchedStations( ArbitrationType type = SegmentAndTrackArbitration ) const;
    /// get bit map of stations with matched segments
    /// bits 0-1-2-3 = DT stations 1-2-3-4
    /// bits 4-5-6-7 = CSC stations 1-2-3-4
    unsigned int stationMask( ArbitrationType type = SegmentAndTrackArbitration ) const;
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
    static const unsigned int PFMuon =  1<<5;

    void setType( unsigned int type ) { type_ = type; }
    unsigned int type() const { return type_; }
    // override of method in base class reco::Candidate
    bool isMuon() const { return true; }
    bool isGlobalMuon()     const { return type_ & GlobalMuon; }
    bool isTrackerMuon()    const { return type_ & TrackerMuon; }
    bool isStandAloneMuon() const { return type_ & StandAloneMuon; }
    bool isCaloMuon() const { return type_ & CaloMuon; }
    bool isPFMuon() {return type_ & PFMuon;} //fix me ! Has to go to type
    
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
    bool pfIsolationValid_;
    bool qualityValid_;
    /// muon hypothesis compatibility with observer calorimeter energy
    float caloCompatibility_;
    /// Isolation information for two cones with dR=0.3 and dR=0.5
    MuonIsolation isolationR03_;
    MuonIsolation isolationR05_;

    /// PF Isolation information for two cones with dR=0.3 and dR=0.4
    MuonPFIsolation pfIsolationR03_;
    MuonPFIsolation pfIsolationR04_;

    /// muon type mask
    unsigned int type_;

    //PF muon ref
    PFCandidateRef pfMuon_;

    // FixMe: Still missing trigger information

    /// get vector of muon chambers for given station and detector
    const std::vector<const MuonChamberMatch*> chambers( int station, int muonSubdetId ) const;
    /// get pointers to best segment and corresponding chamber in vector of chambers
    std::pair<const MuonChamberMatch*,const MuonSegmentMatch*> pair( const std::vector<const MuonChamberMatch*> &,
									ArbitrationType type = SegmentAndTrackArbitration ) const;
     
   public:
     /// get number of segments
     int numberOfSegments( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     /// get deltas between (best) segment and track
     /// If no chamber or no segment returns 999999
     float dX       ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float dY       ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float dDxDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float dDyDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float pullX    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration, bool includeSegmentError = true ) const;
     float pullY    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration, bool includeSegmentError = true ) const;
     float pullDxDz ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration, bool includeSegmentError = true ) const;
     float pullDyDz ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration, bool includeSegmentError = true ) const;
     /// get (best) segment information
     /// If no segment returns 999999
     float segmentX       ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float segmentY       ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float segmentDxDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float segmentDyDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float segmentXErr    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float segmentYErr    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float segmentDxDzErr ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float segmentDyDzErr ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     /// get track information in chamber that contains (best) segment
     /// If no segment, get track information in chamber that has the most negative distance between the track
     /// and the nearest chamber edge (the chamber with the deepest track)
     /// If no chamber returns 999999
     float trackEdgeX   ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackEdgeY   ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackX       ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackY       ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackDxDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackDyDz    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackXErr    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackYErr    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackDxDzErr ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackDyDzErr ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackDist    ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     float trackDistErr ( int station, int muonSubdetId, ArbitrationType type = SegmentAndTrackArbitration ) const;
     
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
