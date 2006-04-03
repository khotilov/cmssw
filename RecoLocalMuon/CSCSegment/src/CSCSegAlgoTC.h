#ifndef CSCSegment_CSCSegAlgoTC_h
#define CSCSegment_CSCSegAlgoTC_h

/**
 * \class CSCSegAlgoTC
 *
 * This is an alternative algorithm for building endcap muon track segments
 * out of the rechit's in a CSCChamber. cf. CSCSegmentizerSK.<BR>
 * 'TC' = 'Tim Cox' = Try (all) Combinations <BR>
 *
 * A CSCSegment isa BasicRecHit4D, and is built from
 * CSCRhit objects, each of which isa BasicRecHit2D. <BR>
 *
 * This class is used by the CSCSegmentRecDet. <BR>
 * Alternative algorithms can be used for the segment building
 * by writing classes like this, and then selecting which one is actually
 * used via the CSCSegmentizerBuilder in CSCDetector. <BR>
 *
 * \author Matteo Sani
 *  Porting of the original algorithm from ORCA to CMSSW 2006-04-03
 */

/* #include "Muon/MAbsAlgo/interface/AbsAlgorithm.h" */
/* #include "CommonDet/BasicDet/interface/Det.h" */
/* #include "Utilities/UI/interface/SimpleConfigurable.h" */
/* #include "Utilities/UI/interface/Verbosity.h" */

#include <RecoLocalMuon/CSCSegment/src/CSCSegmentAlgorithm.h>

#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>

#include <deque>
#include <vector>

//class CSCChamber;
//class CSCChamberSpecs;
//class RecHit;
//class CSCSegment;


class CSCSegAlgoTC : public CSCSegmentAlgorithm {
public:

    // Tim tried using map as basic container of all (space-point) RecHit's in a chamber:
    // The 'key' is a pseudo-layer number (1-6 but with 1 always closest to IP).
    // The 'value' is a vector of the RecHit's on that layer.
    // Using the layer number like this removes the need to sort in global z.
    // Instead we just have to ensure the layer index is correctly adjusted 
    // to enforce the requirement that 'layer 1' is closest in the chamber
    // to the IP.
    // However a map is very painful to use when handling the original 'SK' algorithm,
    // particularly when we need to flag a hit as 'used'. Because of this I am now
    // favouring a pair of vectors, vector<RecHit> and vector<int> for the layer id.

    typedef std::vector<CSCRecHit2D> ChamberHitContainer;
    typedef ChamberHitContainer::const_iterator ChamberHitContainerCIt;

    // We need to be able to flag a hit as 'used' and so need a container
    // of bool's. Naively, this would be vector<bool>... but AVOID that since it's
    // non-standard i.e. packed-bit implementation which is not a standard STL container. 
    // We don't need what it offers and it could lead to unexpected trouble in the future

    typedef std::deque<bool> BoolContainer;

    explicit CSCSegAlgoTC(const edm::ParameterSet& ps);
    virtual ~CSCSegAlgoTC() {};

    /**
     * Build track segments in this chamber (this is where the actual
     * segment-building algorithm hides.)
     */
    CSCSegmentCollection buildSegments(ChamberHitContainer rechits);

    /**
     * Here we must implement the algorithm
     */
    CSCSegmentCollection run(const CSCChamber* aChamber, ChamberHitContainer rechits);

private:

    // Utility functions 
    
    bool addHit(const CSCRecHit2D& aHit, int layer);
    bool replaceHit(const CSCRecHit2D& h, int layer);
    void compareProtoSegment(const CSCRecHit2D& h, int layer);
    void increaseProtoSegment(const CSCRecHit2D& h, int layer);
    AlgebraicSymMatrix calculateError() const;
    HepMatrix derivativeMatrix() const;
    AlgebraicSymMatrix weightMatrix() const;

    /**
     * Return true if the difference in (local) x of two hits is < dRPhiMax
     */
    bool areHitsCloseInLocalX(const CSCRecHit2D& h1, const CSCRecHit2D& h2) const;

    /**
     * Return true if the difference in (global) phi of two hits is < dPhiMax
     */
    bool areHitsCloseInGlobalPhi(const CSCRecHit2D& h1, const CSCRecHit2D& h2) const;
    
    bool hasHitOnLayer(int layer) const;

    /**
     * Return true if hit is near segment.
     * 'Near' means deltaphi and rxy*deltaphi are within ranges
     * specified by orcarc parameters dPhiFineMax and dRPhiFineMax,
     * where rxy = sqrt(x**2+y**2) of the hit in global coordinates.
     */
    bool isHitNearSegment(const CSCRecHit2D& h) const;

    /**
      * Dump global and local coordinate of each rechit in chamber after sort in z
      */
    void dumpHits(const ChamberHitContainer& rechits) const;

    /**
     * Try adding non-used hits to segment
     */
     void tryAddingHitsToSegment(const ChamberHitContainer& rechits, 
        const ChamberHitContainerCIt i1, 
        const ChamberHitContainerCIt i2);

    /**
     * Return true if segment is good.
     * In this algorithm, this means it shares no hits with any other segment.
     */
    bool isSegmentGood(const CSCSegment& seg, 
        const ChamberHitContainer& rechitsInChamber,  
        BoolContainer& used) const;

    /**
     * Flag hits on segment as used
     */
    void flagHitsAsUsed(const CSCSegment& seg, 
        const ChamberHitContainer& rechitsInChamber, BoolContainer& used) const;

    /**
     * Order segments by quality (chi2/#hits) and select the best,
     * requiring that they have unique hits.
     */
    void pruneTheSegments(CSCSegmentCollection& segments, 
        const ChamberHitContainer& rechitsInChamber);
    /**
     * Sort criterion for segment quality, for use in pruneTheSegments.
     */   
    void segmentSort(CSCSegmentCollection& segs);  
    
    float phiAtZ(float z) const;  
    void fillLocalDirection();
    void fillChiSquared();
    void fitSlopes();
    void updateParameters();

    /**
     * Pass final segments to RecDet (the CSCChamber)
     */
//    void passToRecDet(const SegmentContainer& segments, 
//         const CSCChamber* theChamber ) const;


    // Member variables
    // ================


    const CSCChamber* theChamber;
    ChamberHitContainer proto_segment;
    double theChi2;
    LocalPoint theOrigin;
    LocalVector theDirection;
    float uz, vz;
    
    /** max segment chi squared
     */
    float chi2Max;

    /** max hit deviation in r-phi from the segment axis.
     *  Function hitNearSegment requires rxy*abs(deltaphi) < dRPhiFineMax.
     */
    float dRPhiFineMax;

    /** max hit deviation in global phi from the segment axis.
     *  Function hitNearSegment requires abs(deltaphi) < dPhiFineMax.
     */
    float dPhiFineMax;

    /** max distance in local x between hits in one segment
     * @@ The name is historical!
     */
    float dRPhiMax;

    /** max distance in global phi between hits in one segment
     */
    float dPhiMax;

    /** Require end-points of segment are at least minLayersApart 
     */
    int minLayersApart;

    /** Name of this class
     */
    const std::string myName;
    bool debugInfo;
};

#endif
