#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "DataFormats/HcalDetId/interface/HcalOtherDetId.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include <iostream>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace HcalUnpacker_impl {
  template <class DigiClass>
  const HcalQIESample* unpack(const HcalQIESample* startPoint, const HcalQIESample* limit, DigiClass& digi, int presamples, const HcalElectronicsId& eid, int startSample, int endSample) {
    // set parameters
    digi.setPresamples(presamples);
    digi.setReadoutIds(eid);

    // what is my sample number?
    int myFiberChan=startPoint->fiberAndChan();
    int ncurr=0;
    const HcalQIESample* qie_work=startPoint;
    while (qie_work!=limit && qie_work->fiberAndChan()==myFiberChan) {
      if (ncurr>=startSample && ncurr<=endSample) {
	digi.setSample(digi.size(),*qie_work);
	digi.setSize(digi.size()+1);
      }
      ncurr++;
      qie_work++;
    }
    return qie_work;
  }
}


void HcalUnpacker::unpack(const FEDRawData& raw, const HcalElectronicsMap& emap,
			  Collections& colls, HcalUnpackerReport& report) {

  if (raw.size()<16) {
    edm::LogWarning("Invalid Data") << "Empty/invalid DCC data, size = " << raw.size();
    return;
  }

  // get the DCC header
  const HcalDCCHeader* dccHeader=(const HcalDCCHeader*)(raw.data());
  int dccid=dccHeader->getSourceId()-sourceIdOffset_;

  // check the summary status
  
  // walk through the HTR data...
  HcalHTRData htr;
  const unsigned short* daq_first, *daq_last, *tp_first, *tp_last;
  const HcalQIESample* qie_begin, *qie_end, *qie_work;
  const HcalTriggerPrimitiveSample *tp_begin, *tp_end, *tp_work; 
  for (int spigot=0; spigot<HcalDCCHeader::SPIGOT_COUNT; spigot++) {
    if (!dccHeader->getSpigotPresent(spigot)) continue;

    dccHeader->getSpigotData(spigot,htr);
    // check
    if (!htr.check()) {
      edm::LogWarning("Invalid Data") << "Invalid HTR data observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
      report.countSpigotFormatError();
      continue;
    }
    if (htr.isHistogramEvent()) {
      edm::LogWarning("Invalid Data") << "Histogram data passed to non-histogram unpacker on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
      continue;

    }
    // calculate "real" number of presamples
    int nps=htr.getNPS()-startSample_;
    
    // get pointers
    htr.dataPointers(&daq_first,&daq_last,&tp_first,&tp_last);
    unsigned int smid=htr.getSubmodule();
    int htr_tb=smid&0x1;
    int htr_slot=(smid>>1)&0x1F;
    int htr_cr=(smid>>6)&0x1F;
    
    tp_begin=(HcalTriggerPrimitiveSample*)tp_first;
    tp_end=(HcalTriggerPrimitiveSample*)(tp_last+1); // one beyond last..
    
    /// work through the samples
    int currFiberChan=0x3F; // invalid fiber+channel...
    int ncurr=0;
    bool valid=false;
    
    /*
      Unpack the trigger primitives
    */
    for (tp_work=tp_begin; tp_work!=tp_end; tp_work++) {
      if (tp_work->raw()==0xFFFF) continue; // filler word
      if (tp_work->fiberAndChan()!=currFiberChan) { // start new set
	currFiberChan=tp_work->fiberAndChan();
	// lookup the right channel
	HcalElectronicsId eid(tp_work->fiberChan(),tp_work->fiber(),spigot,dccid);
	eid.setHTR(htr_cr,htr_slot,htr_tb);
	DetId did=emap.lookupTrigger(eid);
	if (did.null()) {
	  report.countUnmappedTPDigi();
	  if (unknownIdsTrig_.find(eid)==unknownIdsTrig_.end()) {
	    edm::LogWarning("HCAL") << "HcalUnpacker: No trigger primitive match found for electronics id :" << eid;
	    unknownIdsTrig_.insert(eid);
	  }
	  valid=false;
	  continue;
	} else if (did==HcalTrigTowerDetId::Undefined || 
		   (did.det()==DetId::Hcal && did.subdetId()==0)) {
	  // known to be unmapped
	  valid=false;
	  continue;
	}
	HcalTrigTowerDetId id(did);
	colls.tpCont->push_back(HcalTriggerPrimitiveDigi(id));
	// set the various bits
	colls.tpCont->back().setPresamples(nps);
	// no hits recorded for current
	ncurr=0;
	valid=true;
      }
      // add the word (if within settings) [ TODO: correct behavior when just one TP... ]
      if (valid && ncurr>=startSample_ && ncurr<=endSample_) {
	colls.tpCont->back().setSample(colls.tpCont->back().size(),*tp_work);
	colls.tpCont->back().setSize(colls.tpCont->back().size()+1);
      }
      ncurr++;
    }


    qie_begin=(HcalQIESample*)daq_first;
    qie_end=(HcalQIESample*)(daq_last+1); // one beyond last..

    /// work through the samples
    currFiberChan=0x3F; // invalid fiber+channel...
    ncurr=0;
    valid=false;
    
    for (qie_work=qie_begin; qie_work!=qie_end; ) {
      if (qie_work->raw()==0xFFFF) {
	qie_work++;
	continue; // filler word
      }
      // always at the beginning ...
      currFiberChan=qie_work->fiberAndChan();

      // lookup the right channel
      HcalElectronicsId eid(qie_work->fiberChan(),qie_work->fiber(),spigot,dccid);
      eid.setHTR(htr_cr,htr_slot,htr_tb);
      DetId did=emap.lookup(eid);

      if (!did.null()) {
	switch (((HcalSubdetector)did.subdetId())) {
	case (HcalBarrel):
	case (HcalEndcap): {
	  colls.hbheCont->push_back(HBHEDataFrame(HcalDetId(did)));
	  qie_work=HcalUnpacker_impl::unpack<HBHEDataFrame>(qie_work, qie_end, colls.hbheCont->back(), nps, eid, startSample_, endSample_);
	} break;
	case (HcalOuter): {
	  colls.hoCont->push_back(HODataFrame(HcalDetId(did)));
	  qie_work=HcalUnpacker_impl::unpack<HODataFrame>(qie_work, qie_end, colls.hoCont->back(), nps, eid, startSample_, endSample_);
	} break;
	case (HcalForward): {
	  colls.hfCont->push_back(HFDataFrame(HcalDetId(did)));
	  qie_work=HcalUnpacker_impl::unpack<HFDataFrame>(qie_work, qie_end, colls.hfCont->back(), nps, eid, startSample_, endSample_);
	} break;
	case (HcalOther) : {
	  HcalOtherDetId odid(did);
	  if (odid.subdet()==HcalCalibration) {
	    colls.calibCont->push_back(HcalCalibDataFrame(HcalCalibDetId(did)));
	    qie_work=HcalUnpacker_impl::unpack<HcalCalibDataFrame>(qie_work, qie_end, colls.calibCont->back(), nps, eid, startSample_, endSample_); 
	  } else if (odid.subdet()==HcalZDC) {
	    colls.zdcCont->push_back(ZDCDataFrame(HcalZDCDetId(did)));
	    qie_work=HcalUnpacker_impl::unpack<ZDCDataFrame>(qie_work, qie_end, colls.zdcCont->back(), nps, eid, startSample_, endSample_); 
	  }
	} break;
	case (HcalEmpty): 
	default: {
	  for (int fiberC=qie_work->fiberAndChan();
	       qie_work!=qie_end && qie_work->fiberAndChan()==fiberC;
	       qie_work++);
	}
	  break;
	}
      } else {
	report.countUnmappedDigi();
	if (unknownIds_.find(eid)==unknownIds_.end()) {
	  edm::LogWarning("HCAL") << "HcalUnpacker: No match found for electronics id :" << eid;
	  unknownIds_.insert(eid);
	}
	for (int fiberC=qie_work->fiberAndChan();
	     qie_work!=qie_end && qie_work->fiberAndChan()==fiberC;
	     qie_work++);
      }
    }
  }
}

HcalUnpacker::Collections::Collections() {
  hbheCont=0;
  hoCont=0;
  hfCont=0;
  tpCont=0;
  zdcCont=0;
  calibCont=0;
}

void HcalUnpacker::unpack(const FEDRawData& raw, const HcalElectronicsMap& emap, std::vector<HBHEDataFrame>& container, std::vector<HcalTriggerPrimitiveDigi>& tp) {
  Collections c;
  c.hbheCont=&container;
  c.tpCont=&tp;
  HcalUnpackerReport r;
  unpack(raw,emap,c,r);
}

void HcalUnpacker::unpack(const FEDRawData& raw, const HcalElectronicsMap& emap, std::vector<HODataFrame>& container, std::vector<HcalTriggerPrimitiveDigi>& tp) {
  Collections c;
  c.hoCont=&container;
  c.tpCont=&tp;
  HcalUnpackerReport r;
  unpack(raw,emap,c,r);
}

void HcalUnpacker::unpack(const FEDRawData& raw, const HcalElectronicsMap& emap, std::vector<HFDataFrame>& container, std::vector<HcalTriggerPrimitiveDigi>& tp) {
  Collections c;
  c.hfCont=&container;
  c.tpCont=&tp;
  HcalUnpackerReport r;
  unpack(raw,emap,c,r);
}

void HcalUnpacker::unpack(const FEDRawData& raw, const HcalElectronicsMap& emap, std::vector<HcalHistogramDigi>& histoDigis) {

  // get the DCC header
  const HcalDCCHeader* dccHeader=(const HcalDCCHeader*)(raw.data());
  int dccid=dccHeader->getSourceId()-sourceIdOffset_;
  
  // check the summary status
  
  // walk through the HTR data...
  HcalHTRData htr;
  for (int spigot=0; spigot<HcalDCCHeader::SPIGOT_COUNT; spigot++) {
    if (!dccHeader->getSpigotPresent(spigot)) continue;
    
    dccHeader->getSpigotData(spigot,htr);
    // check
    if (!htr.check()) {
      edm::LogWarning("Invalid Data") << "Invalid HTR data observed on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
      continue;
    }
    if (!htr.isHistogramEvent()) {
      edm::LogWarning("Invalid Data") << "Non-histogram data passed to histogram unpacker on spigot " << spigot << " of DCC with source id " << dccHeader->getSourceId();
      continue;
    }

    unsigned int smid=htr.getSubmodule();
    int htr_tb=smid&0x1;
    int htr_slot=(smid>>1)&0x1F;
    int htr_cr=(smid>>6)&0x1F;
    
    // find out the fibers
    int f[2],fc;
    htr.getHistogramFibers(f[0],f[1]);

    for (int nf=0; nf<2; nf++) {
      if (f[nf]<0 || nf==1 && f[0]==f[1]) continue; // skip if invalid or the same
      for (fc=0; fc<=2; fc++) {
	HcalElectronicsId eid(fc,f[nf],spigot,dccid);	  
	eid.setHTR(htr_cr,htr_slot,htr_tb);
	DetId did=emap.lookup(eid);

	if (did.null() || did.det()!=DetId::Hcal || did.subdetId()==0) {
	  if (unknownIds_.find(eid)==unknownIds_.end()) {
	    edm::LogWarning("HCAL") << "HcalHistogramUnpacker: No match found for electronics id :" << eid;
	    unknownIds_.insert(eid);
	  }	  
	  continue;
	}
	histoDigis.push_back(HcalHistogramDigi(HcalDetId(did))); // add it!
	HcalHistogramDigi& digi=histoDigis.back();
	
	// unpack the four capids
	for (int capid=0; capid<4; capid++) 
	  htr.unpackHistogram(f[nf],fc,capid,digi.getArray(capid));
	
      }
    }
  }
}      

