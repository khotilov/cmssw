#include "RecoLocalTracker/SiStripZeroSuppression/interface/PercentileCMNSubtractor.h"

void PercentileCMNSubtractor::subtract(const uint32_t& detId,std::vector<int16_t>& digis) {subtract_(detId,digis);}
void PercentileCMNSubtractor::subtract(const uint32_t& detId,std::vector<float>& digis) {subtract_(detId,digis);}

template<typename T> 
inline
void PercentileCMNSubtractor::
subtract_(const uint32_t& detId,std::vector<T>& digis){
  
  std::vector<T> tmp;  tmp.reserve(128);  
  typename std::vector<T>::iterator  
    strip( digis.begin() ), 
    end(   digis.end()   ),
    endAPV;
  
  while( strip < end ) {
    endAPV = strip+128; tmp.clear();
    tmp.insert(tmp.end(),strip,endAPV);
    const float offset = percentile(tmp,percentile_);

    while (strip < endAPV) {
      *strip = static_cast<T>(*strip-offset);
      strip++;
    }

  }
}


template<typename T>
inline
float PercentileCMNSubtractor::
percentile( std::vector<T>& sample, double pct) {
  typename std::vector<T>::iterator mid = sample.begin() + int(sample.size()*pct/100.0);
  std::nth_element(sample.begin(), mid, sample.end());
  return *mid;
} 

