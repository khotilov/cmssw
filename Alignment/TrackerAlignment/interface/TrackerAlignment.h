#ifndef Alignment_TrackerAlignment_TrackerAlignment_H
#define Alignment_TrackerAlignment_TrackerAlignment_H

/** \class TrackerAlignment
 *  The TrackerAlignment helper class for alignment jobs
 *  - Rotates and Translates any module for the tracker based on rawId
 *
 *  \author Nhan Tran
 */

#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"

class TrackerAlignment{
 
public:
  TrackerAlignment( const edm::EventSetup& setup );
	
  ~TrackerAlignment();
	
  AlignableTracker* getAlignableTracker() { return theAlignableTracker; }
	
  void moveAlignablePixelEndCaps( int rawId, std::vector<float> localDisplacements, std::vector<float> localRotations  );
  void moveAlignableEndCaps( int rawId, std::vector<float> localDisplacements, std::vector<float> localRotations  );
  void moveAlignablePixelHalfBarrels( int rawId, std::vector<float> localDisplacements, std::vector<float> localRotations  );
  void moveAlignableInnerHalfBarrels( int rawId, std::vector<float> localDisplacements, std::vector<float> localRotations  );
  void moveAlignableOuterHalfBarrels( int rawId, std::vector<float> localDisplacements, std::vector<float> localRotations  );	
  void moveAlignableTIDs( int rawId, std::vector<float> localDisplacements, std::vector<float> localRotations  );
  
  void saveToDB();
  int rawid;
  
  std::vector<float> local_displacements;
  
  std::vector<float> local_rotations;
	
  
private:  
  AlignableTracker* theAlignableTracker;

  std::string theAlignRecordName, theErrorRecordName;
  
};
#endif //TrackerAlignment_H
