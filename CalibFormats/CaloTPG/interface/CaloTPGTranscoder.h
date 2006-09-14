#ifndef CALIBFORMATS_CALOTPG_CALOTPGTRANSCODER_H
#define CALIBFORMATS_CALOTPG_CALOTPGTRANSCODER_H 1

#include <boost/shared_ptr.hpp>
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveSample.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveSample.h"

class HcalTPGCompressor;
class EcalTPGCompressor;

namespace edm {
  class EventSetup; 
}

/** \class CaloTPGTranscoder
  *  
  * Abstract interface for the mutual transcoder required for compressing
  * and uncompressing the ET stored in HCAL and ECAL Trigger Primitives
  * 
  * $Date: 2006/09/14 16:24:10 $
  * $Revision: 1.1 $
  * \author J. Mans - Minnesota
  */
class CaloTPGTranscoder {
public:
  CaloTPGTranscoder();
  virtual ~CaloTPGTranscoder(); 

  enum Mode { All=0, RCT=1, HcalTPG=2, EcalTPG=3 };
  /// Obtain any needed objects from the EventSetup
  virtual void setup(const edm::EventSetup& es, Mode mode=All);
  /// Release any objects obtained from the EventSetup
  virtual void releaseSetup();
  /** \brief Compression from linear samples+fine grain in the HTR */
  virtual HcalTriggerPrimitiveSample hcalCompress(const HcalTrigTowerDetId& id, unsigned int sample, bool fineGrain) const = 0;
  /** \brief Compression from linear samples+fine grain in the ECAL */
  virtual EcalTriggerPrimitiveSample ecalCompress(const EcalTrigTowerDetId& id, unsigned int sample, bool fineGrain) const = 0;
  /** \brief Uncompression for the Electron/Photon path in the RCT */
  virtual void rctEGammaUncompress(const HcalTrigTowerDetId& hid, const HcalTriggerPrimitiveSample& hc,
				   const EcalTrigTowerDetId& eid, const EcalTriggerPrimitiveSample& ec, 
				   unsigned int& et, bool& egVecto, bool& activity) const = 0;
  /** \brief Uncompression for the JET path in the RCT */
  virtual void rctJetUncompress(const HcalTrigTowerDetId& hid, const HcalTriggerPrimitiveSample& hc,
				   const EcalTrigTowerDetId& eid, const EcalTriggerPrimitiveSample& ec, 
				   unsigned int& et) const = 0;

  boost::shared_ptr<HcalTPGCompressor> getHcalCompressor() const { return hccompress_; }
  boost::shared_ptr<EcalTPGCompressor> getEcalCompressor() const { return eccompress_; }
private:
  boost::shared_ptr<HcalTPGCompressor> hccompress_;
  boost::shared_ptr<EcalTPGCompressor> eccompress_;
};

#endif
