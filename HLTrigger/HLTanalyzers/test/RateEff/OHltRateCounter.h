//////////////////////////////////////////////////////////
//
// Class to store  and process rate counts
//
//////////////////////////////////////////////////////////

#ifndef OHltRateCounter_h
#define OHltRateCounter_h

#include <vector>
#include <libconfig.h++>
#include <TMath.h>

using namespace std;
using namespace libconfig;

class OHltRateCounter {
 public:

  OHltRateCounter(unsigned int size);
  virtual ~OHltRateCounter(){};

  // Helper functions
  static inline float eff(int a, int b){ 
    if (b==0.){return -1.;}
    float af = float(a),bf = float(b),effi = af/bf;
    return effi;
  }
  static inline float effErr(int a, int b){
    if (b==0.){return -1.;}
    float af = float(a),bf = float(b),r = af/bf;
    float unc = sqrt(af + (r*r*bf) )/bf;
    return unc;
  }
  static inline float eff(float a, float b){ 
    if (b==0.){return -1.;}
    float af = float(a),bf = float(b),effi = af/bf;
    return effi;
  }
  static inline float effErr(float a, float b){
    if (b==0.){return -1.;}
    float af = float(a),bf = float(b),r = af/bf;
    float unc = sqrt(af + (r*r*bf) )/bf;
    return unc;
  }
  
  
  // Data
  vector<int> iCount;
  vector<int> sPureCount;
  vector<int> pureCount;
  vector< vector<int> > overlapCount;

};
#endif
