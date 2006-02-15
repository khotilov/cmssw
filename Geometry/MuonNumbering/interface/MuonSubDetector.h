#ifndef MuonNumbering_MuonSubDetector_h
#define MuonNumbering_MuonSubDetector_h

/** \class MuonSubDetector
 *
 * class to handle muon sensitive detectors,
 * possible arguments for constructor:
 * "MuonDTHits", "MuonCSCHits", "MuonRPCHits"
 *
 * the function suIdName() returns the detector SuId
 * for the ROU factory
 *  
 *  $Date: 2005/10/26 17:45:00 $
 *  $Revision: 1.2 $
 * \author Arno Straessner, CERN <arno.straessner@cern.ch>
 *
 */

#include<string>

class MuonSubDetector {
 public:

  /*  
   * possible arguments for constructor:
   * "MuonDTHits", "MuonCSCHits", "MuonRPCHits"
   */

  MuonSubDetector(std::string name);
  ~MuonSubDetector(){};

  bool isBarrel();
  bool isEndcap();
  bool isRpc();
  std::string name();
  std::string suIdName();
      
 private:
  enum subDetector {barrel,endcap,rpc,nodef};
  subDetector detector;
  std::string detectorName;
};

#endif
