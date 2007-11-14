#ifndef ANALYSISDATAFORMATS_SISTRIPCLUSTERINFO_H
#define ANALYSISDATAFORMATS_SISTRIPCLUSTERINFO_H

#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include <vector>
#include <sstream>
#include "boost/cstdint.hpp"

class SiStripClusterInfo {
public:

  SiStripClusterInfo() : detId_(0) {};

  SiStripClusterInfo(const SiStripCluster& cluster):
    detId_(cluster.geographicalId()),
    FirstStrip(cluster.firstStrip()),
    StripAmplitudes(cluster.amplitudes()),    
    Position(cluster.barycenter()),
    Width(cluster.amplitudes().size()),
    MaxCharge(0),
    MaxPosition(0)
    {
      RawDigiAmplitudesL=std::vector<float>(0);
      RawDigiAmplitudesR=std::vector<float>(0);
      StripNoises=std::vector<float>(0);
      ApvGains=std::vector<float>(1);
    };
  

  uint16_t firstStrip() const {return FirstStrip;}
  unsigned int geographicalId() const {return detId_;}
  const std::vector<uint16_t>& stripAmplitudes()    const {return StripAmplitudes;}
  const std::vector<float>&    rawdigiAmplitudesL() const {return RawDigiAmplitudesL;}
  const std::vector<float>&    rawdigiAmplitudesR() const {return RawDigiAmplitudesR;}
  const std::vector<float>&    stripNoises()        const {return StripNoises;}
  const std::vector<float>&    apvGains()           const {return ApvGains;}
  
  float charge()    const {return Charge;}
  float noise()     const {return Noise;}
  float position()  const {return Position;}
  float width()     const {return Width;}
  float maxCharge() const {return MaxCharge;}
  uint16_t maxPos() const {return MaxPosition;}
  float chargeL()   const {return ChargeL;}
  float chargeR()   const {return ChargeR;}
  float getGainForStripNb(uint16_t istrip) const;

  void setCharge(const float& value) {Charge=value;}  
  void setNoise(const float& value) {Noise=value;}
  void setStripNoises(std::vector<float>& value) {StripNoises=value;}
  void setApvGains(std::vector<float>& value) {ApvGains=value;}
  void setMaxCharge(const float& value) {MaxCharge=value;}
  void setMaxPos(const uint16_t& value) {MaxPosition=value;}
  void setChargeL(const float& value) {ChargeL=value;}
  void setChargeR(const float& value) {ChargeR=value;}
  void setRawDigiAmplitudesL(const std::vector<float>& value){RawDigiAmplitudesL=value;}
  void setRawDigiAmplitudesR(const std::vector<float>& value){RawDigiAmplitudesR=value;}

  void print(std::stringstream &ss);

private:

  uint32_t                detId_;
  uint16_t                FirstStrip;
  std::vector<uint16_t>   StripAmplitudes;
  std::vector<float>      RawDigiAmplitudesL;
  std::vector<float>      RawDigiAmplitudesR;
  std::vector<float>      StripNoises;
  std::vector<float>      ApvGains;

  float     Charge;
  float     Noise;
  float     Position;
  uint16_t  Width;
  float     MaxCharge;
  uint16_t  MaxPosition;
  float     ChargeL;
  float     ChargeR;
};

// Comparison operators
inline bool operator<( const SiStripClusterInfo& one, const SiStripClusterInfo& other) {
  if(one.geographicalId() == other.geographicalId()) {
    float StoN_one = (one.noise()!= 0) ? one.charge()/one.noise() : .0;
    float StoN_other = (other.noise()!= 0) ? other.charge()/other.noise() : .0;
 
    return StoN_one < StoN_other;
  }
  return one.geographicalId() < other.geographicalId();
}
#endif // ANALYSISDATAFORMATS_SISTRIPCLUSTER_H
