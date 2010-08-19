#include "RecoLocalTracker/SiStripZeroSuppression/interface/IteratedMedianCMNSubtractor.h"

#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CondFormats/DataRecord/interface/SiStripNoisesRcd.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
#include <cmath>

void IteratedMedianCMNSubtractor::init(const edm::EventSetup& es){
  uint32_t n_cache_id = es.get<SiStripNoisesRcd>().cacheIdentifier();
  uint32_t q_cache_id = es.get<SiStripQualityRcd>().cacheIdentifier();

  if(n_cache_id != noise_cache_id) {
    es.get<SiStripNoisesRcd>().get( noiseHandle );
    noise_cache_id = n_cache_id;
  }
  if(q_cache_id != quality_cache_id) {
    es.get<SiStripQualityRcd>().get( qualityHandle );
    quality_cache_id = q_cache_id;
  }
}

void IteratedMedianCMNSubtractor::subtract(const uint32_t& detId,std::vector<int16_t>& digis){ subtract_(detId,digis);}
void IteratedMedianCMNSubtractor::subtract(const uint32_t& detId,std::vector<float>& digis){ subtract_(detId,digis);}

template<typename T>
inline
void IteratedMedianCMNSubtractor::
subtract_(const uint32_t& detId,std::vector<T>& digis){

  SiStripNoises::Range detNoiseRange = noiseHandle->getRange(detId);
  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId);

  typename std::vector<T>::iterator fs,ls;
  float offset = 0;  
  std::vector< std::pair<float,float> > subset;
  subset.reserve(128); 

  for( uint16_t APV=0; APV< digis.size()/128; ++APV)
  {
    subset.clear();
    // fill subset vector with all good strips and their noises
    for (uint16_t istrip=APV*128; istrip<(APV+1)*128; ++istrip)
    {
      if ( !qualityHandle->IsStripBad(detQualityRange,istrip) )
      {
        std::pair<float,float> pin((float)digis[istrip], (float)noiseHandle->getNoise(istrip,detNoiseRange));
        subset.push_back( pin );
      }
    }

    // caluate offset for all good strips (first iteration)
    if (subset.size() != 0)
      offset = pairMedian(subset);

    // for second, third... iterations, remove strips over threshold
    // and recalculate offset on remaining strips
    for ( int ii = 0; ii<iterations_-1; ++ii )
    {
      std::vector< std::pair<float,float> >::iterator si = subset.begin();
      while(  si != subset.end() )
      {
        if( si->first-offset < cut_to_avoid_signal_*si->second )  
          si = subset.erase(si);
        else
          ++si;
      }
      if ( subset.size() == 0 ) break;
      offset = pairMedian(subset);
    }        

    // remove offset
    fs = digis.begin()+APV*128;
    ls = digis.begin()+(APV+1)*128;
    while (fs < ls)
      *fs++ = static_cast<T>(*fs-offset);

  }
}



inline float IteratedMedianCMNSubtractor::pairMedian( std::vector<std::pair<float,float> >& sample) {
  std::vector<std::pair<float,float> >::iterator mid = sample.begin() + sample.size()/2;
  std::nth_element(sample.begin(), mid, sample.end());
  if( sample.size() & 1 ) //odd size
    return (*mid).first;
  return ( (*std::max_element(sample.begin(), mid)).first + (*mid).first ) / 2.;
}

