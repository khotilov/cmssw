#include "EventFilter/CastorRawToDigi/interface/CastorPacker.h"
#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "EventFilter/Utilities/interface/Crc.h"

CastorPacker::Collections::Collections() {
  castorCont=0;
  tpCont=0;
  calibCont=0;
}

template <class Coll, class DetIdClass> 
int process(const Coll* pt, const DetId& did, unsigned short* buffer, int& presamples) {
  if (pt==0) return 0;
  int size=0;
  typename Coll::const_iterator i=pt->find(DetIdClass(did));
  if (i!=pt->end()) {
    presamples=i->presamples();
    size=i->size();
    for (int j=0; j<size; j++) 
      buffer[j]=(*i)[j].raw();
  }
  return size;
}

static unsigned char processTrig(const HcalTrigPrimDigiCollection* pt, const HcalTrigTowerDetId& tid, unsigned short* buffer) {
  if (pt==0) return 0;
  int size=0;
  HcalTrigPrimDigiCollection::const_iterator i=pt->find(tid);
  if (i!=pt->end()) {
    int presamples=i->presamples();
    int samples=i->size();

    for (int j=0; j<samples; j++) {
       buffer[j]=(*i)[j].raw();
       if (j==presamples) buffer[j]|=0x0200;
    }
  }
  return size;
}

int CastorPacker::findSamples(const DetId& did, const Collections& inputs,
			    unsigned short* buffer, int &presamples) {
  if (did.det()!=DetId::Hcal) return 0;
  int size=0;
  HcalGenericDetId genId(did);
  size=process<CastorDigiCollection,HcalDetId>(inputs.castorCont,did,buffer,presamples);
  
/*  switch (genId.genericSubdet()) {
  case(HcalGenericDetId::HcalGenBarrel):
  case(HcalGenericDetId::HcalGenEndcap):
    size=process<HBHEDigiCollection,HcalDetId>(inputs.hbhe,did,buffer,presamples);
    break;
  case(HcalGenericDetId::HcalGenOuter):
    size=process<HODigiCollection,HcalDetId>(inputs.hoCont,did,buffer,presamples);
    break;
  case(HcalGenericDetId::HcalGenForward):
    size=process<HFDigiCollection,HcalDetId>(inputs.hfCont,did,buffer,presamples);
    break;
  case(HcalGenericDetId::HcalGenZDC):
    size=process<ZDCDigiCollection,HcalZDCDetId>(inputs.zdcCont,did,buffer,presamples);
    break;
  case(HcalGenericDetId::HcalGenCalibration):
    size=process<HcalCalibDigiCollection,HcalCalibDetId>(inputs.calibCont,did,buffer,presamples);
    break;
  default: size=0;
  } */
  return size;
}

void CastorPacker::pack(int fedid, int dccnumber,
		      int nl1a, int orbitn, int bcn,
		      const Collections& inputs, 
		      const CastorElectronicsMap& emap,
		      FEDRawData& output) {
  std::vector<unsigned short> precdata(HcalHTRData::CHANNELS_PER_SPIGOT*HcalHTRData::MAXIMUM_SAMPLES_PER_CHANNEL);
  std::vector<unsigned short> trigdata(HcalHTRData::CHANNELS_PER_SPIGOT*HcalHTRData::MAXIMUM_SAMPLES_PER_CHANNEL);
  std::vector<unsigned char> preclen(HcalHTRData::CHANNELS_PER_SPIGOT);
  std::vector<unsigned char> triglen(HcalHTRData::CHANNELS_PER_SPIGOT);
  static const int HTRFormatVersion=3;

  HcalHTRData spigots[15];
  // loop over all valid channels in the given dcc, spigot by spigot.
  for (int spigot=0; spigot<15; spigot++) {
    spigots[spigot].allocate(HTRFormatVersion);
    CastorElectronicsId exampleEId;
    int npresent=0;
    int presamples=-1, samples=-1;
    for (int fiber=1; fiber<=8; fiber++) 
      for (int fiberchan=0; fiberchan<3; fiberchan++) {
	int linear=(fiber-1)*3+fiberchan;
	HcalQIESample chanSample(0,0,fiber,fiberchan,false,false);
	unsigned short chanid=chanSample.raw()&0xF800;
	preclen[linear]=0;

	CastorElectronicsId partialEid(fiberchan,fiber,spigot,dccnumber);
	// does this partial id exist?
	CastorElectronicsId fullEid;
	HcalGenericDetId genId;
	if (!emap.lookup(partialEid,fullEid,genId)) continue;


	// next, see if there is a digi with this id
	unsigned short* database=&(precdata[linear*HcalHTRData::MAXIMUM_SAMPLES_PER_CHANNEL]);
	int mypresamples;
	int mysamples=findSamples(genId,inputs,database,mypresamples);

	if (mysamples>0) {
	  if (samples<0) samples=mysamples;
	  else if (samples!=mysamples) {
	    edm::LogError("HCAL") << "Mismatch of samples in a single HTR (unsupported) " << mysamples << " != " << samples;
	    continue;
	  }
	  if (presamples<0) {
	    presamples=mypresamples;
	    exampleEId=fullEid;
	  } else if (mypresamples!=presamples) {
	    edm::LogError("HCAL") << "Mismatch of presamples in a single HTR (unsupported) " << mypresamples << " != " << presamples;
	    continue;	    
	  }
	  for (int ii=0; ii<samples; ii++)
	    database[ii]=(database[ii]&0x7FF)|chanid;
	  preclen[linear]=(unsigned char)(samples);
	  npresent++;
	}	
      }
    for (int slb=1; slb<=6; slb++) 
      for (int slbchan=0; slbchan<=3; slbchan++) {
	int linear=(slb-1)*4+slbchan;
	HcalTriggerPrimitiveSample idCvt(0,0,slb,slbchan);
	unsigned short chanid=idCvt.raw()&0xF800;
	triglen[linear]=0;

	CastorElectronicsId partialEid(slbchan,slb,spigot,dccnumber,0,0,0);
	// does this partial id exist?
	CastorElectronicsId fullEid;
	HcalTrigTowerDetId tid;
	if (!emap.lookup(partialEid,fullEid,tid)) {
//	  std::cout << "TPGPACK : no match for " << partialEid << std::endl;
	  continue;
	}  //else std::cout << "TPGPACK : converted " << partialEid << " to " << fullEid << "/" << tid << std::endl;
          
	// finally, what about a trigger channel?
	if (!tid.null()) {
	  unsigned short* trigbase=&(trigdata[linear*HcalHTRData::MAXIMUM_SAMPLES_PER_CHANNEL]);
	  triglen[linear]=processTrig(inputs.tpCont,tid,trigbase);
	  
	  for (unsigned char q=0; q<triglen[linear]; q++)
	    trigbase[q]=(trigbase[q]&0x7FF)|chanid;
	}
      }
    /// pack into HcalHTRData
    if (npresent>0) {
      spigots[spigot].pack(&(preclen[0]),&(precdata[0]),
			   &(triglen[0]),&(trigdata[0]),
			   false);
      static const int pipeline=0x22;
      static const int firmwareRev=0;
      int submodule=exampleEId.htrTopBottom()&0x1;
      submodule|=(exampleEId.htrSlot()&0x1F)<<1;
      submodule|=(exampleEId.readoutVMECrateId()&0x1f)<<6;
      spigots[spigot].packHeaderTrailer(nl1a,
					bcn,
					submodule,
					orbitn,
					pipeline,
					samples,
					presamples,
					firmwareRev);
      
    }
  }
  // calculate the total length, and resize the FEDRawData
  int theSize=0;
  for (int spigot=0; spigot<15; spigot++) {
    theSize+=spigots[spigot].getRawLength()*sizeof(unsigned short);
  }
  theSize+=sizeof(HcalDCCHeader)+8; // 8 for trailer
  theSize+=(8-(theSize%8))%8; // even number of 64-bit words.
  output.resize(theSize);
  
  // construct the bare DCC Header
  HcalDCCHeader* dcc=(HcalDCCHeader*)(output.data());
  dcc->clear();
  dcc->setHeader(fedid,bcn,nl1a,orbitn);

  // pack the HTR data into the FEDRawData block using HcalDCCHeader
  for (int spigot=0; spigot<15; spigot++) {
    if (spigots[spigot].getRawLength()>0)
      dcc->copySpigotData(spigot,spigots[spigot],true,0);
  }
  // trailer
  FEDTrailer fedTrailer(output.data()+(output.size()-8));
  fedTrailer.set(output.data()+(output.size()-8),
    output.size()/8,
    evf::compute_crc(output.data(),output.size()), 0, 0);
}
