#include "EventFilter/ESRawToDigi/interface/ESUnpacker.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

const int ESUnpacker::bDHEAD    = 2;
const int ESUnpacker::bDH       = 6;
const int ESUnpacker::bDEL      = 24;
const int ESUnpacker::bDERR     = 8;
const int ESUnpacker::bDRUN     = 24;
const int ESUnpacker::bDRUNTYPE = 32;
const int ESUnpacker::bDTRGTYPE = 16;
const int ESUnpacker::bDCOMFLAG = 8;
const int ESUnpacker::bDORBIT   = 32;
const int ESUnpacker::bDVMAJOR  = 8;
const int ESUnpacker::bDVMINOR  = 8;
const int ESUnpacker::bDCH      = 4; 
const int ESUnpacker::bDOPTO    = 8;

const int ESUnpacker::sDHEAD    = 26;
const int ESUnpacker::sDH       = 24;
const int ESUnpacker::sDEL      = 0;
const int ESUnpacker::sDERR     = bDEL + sDEL;
const int ESUnpacker::sDRUN     = 0;
const int ESUnpacker::sDRUNTYPE = 0;
const int ESUnpacker::sDTRGTYPE = 0;
const int ESUnpacker::sDCOMFLAG = bDTRGTYPE + sDTRGTYPE;
const int ESUnpacker::sDORBIT   = 0;
const int ESUnpacker::sDVMINOR  = 8;
const int ESUnpacker::sDVMAJOR  = bDVMINOR + sDVMINOR;
const int ESUnpacker::sDCH      = 0;
const int ESUnpacker::sDOPTO    = 16;

const int ESUnpacker::bKEC    = 8;   // KCHIP packet event counter
const int ESUnpacker::bKFLAG2 = 8;
const int ESUnpacker::bKBC    = 12;  // KCHIP packet bunch counter
const int ESUnpacker::bKFLAG1 = 4;
const int ESUnpacker::bKET    = 1;
const int ESUnpacker::bKCRC   = 1;
const int ESUnpacker::bKCE    = 1;
const int ESUnpacker::bKID    = 11;
const int ESUnpacker::bFIBER  = 6;   // Fiber number
const int ESUnpacker::bKHEAD1 = 2;
const int ESUnpacker::bKHEAD2 = 2;

const int ESUnpacker::sKEC    = 0;  
const int ESUnpacker::sKFLAG2 = bKEC + sKEC;
const int ESUnpacker::sKBC    = bKFLAG2 + sKFLAG2; 
const int ESUnpacker::sKFLAG1 = bKBC + sKBC;
const int ESUnpacker::sKET    = bKFLAG1 + sKFLAG1;
const int ESUnpacker::sKCRC   = bKET + sKET;
const int ESUnpacker::sKCE    = bKCRC + sKCRC;
const int ESUnpacker::sKID    = bKCE + sKCE + 5;
const int ESUnpacker::sFIBER  = bKID + sKID + 1;  
const int ESUnpacker::sKHEAD1 = bFIBER + sFIBER + 2;
const int ESUnpacker::sKHEAD2 = bKHEAD1 + sKHEAD1;

const int ESUnpacker::bADC0  = 16;
const int ESUnpacker::bADC1  = 16;
const int ESUnpacker::bADC2  = 16;
const int ESUnpacker::bPACE  = 2;
const int ESUnpacker::bSTRIP = 5;
const int ESUnpacker::bE0    = 1;
const int ESUnpacker::bE1    = 1;
const int ESUnpacker::bHEAD  = 2;

const int ESUnpacker::sADC0  = 0;
const int ESUnpacker::sADC1  = bADC0 + sADC0;
const int ESUnpacker::sADC2  = bADC1 + sADC1;
const int ESUnpacker::sPACE  = bADC2 + sADC2;
const int ESUnpacker::sSTRIP = bPACE + sPACE; 
const int ESUnpacker::sE0    = bSTRIP + sSTRIP + 1;
const int ESUnpacker::sE1    = bE0 + sE0;
const int ESUnpacker::sHEAD  = bE1 + sE1 + 4;

ESUnpacker::ESUnpacker(const ParameterSet& ps) 
  : pset_(ps), fedId_(0), run_number_(0), orbit_number_(0), bx_(0), lv1_(0), trgtype_(0)
{

  debug_ = pset_.getUntrackedParameter<bool>("debugMode", false);

}

ESUnpacker::~ESUnpacker() {
}

void ESUnpacker::interpretRawData(int fedId, const FEDRawData & rawData, ESDigiCollection & digis) {
  
  int nWords = rawData.size()/sizeof(Word64);
  if (nWords==0) return;
  int dccWords = 6;
  
  // Event header
  const Word64* header = reinterpret_cast<const Word64* >(rawData.data()); --header;
  bool moreHeaders = true;
  while (moreHeaders) {
    ++header;
    FEDHeader ESHeader( reinterpret_cast<const unsigned char*>(header) );
    if ( !ESHeader.check() ) {
      if (debug_) edm::LogWarning("Invalid Data")<<"ES : Failed header check !";
      return;
    } 

    fedId_ = ESHeader.sourceID();
    lv1_   = ESHeader.lvl1ID();
    bx_    = ESHeader.bxID();

    if (debug_) {
      cout<<"[ESUnpacker]: FED Header candidate. Is header? "<< ESHeader.check();
      if (ESHeader.check())
        cout <<". BXID: "<<bx_<<" SourceID : "<<fedId_<<" L1ID: "<<lv1_<<endl;
      else cout<<" WARNING!, this is not a ES Header"<<endl;
    }

    moreHeaders = ESHeader.moreHeaders();
  }
  if (fedId_ != fedId) {
    if (debug_) edm::LogWarning("Invalid Data")<<"Invalid ES data with source id " <<fedId_;
    return;
  }

  // Event trailer
  const Word64* trailer = reinterpret_cast<const Word64* >(rawData.data())+(nWords-1); ++trailer;
  bool moreTrailers = true;
  while (moreTrailers) {
    --trailer;
    FEDTrailer ESTrailer(reinterpret_cast<const unsigned char*>(trailer));
    if ( !ESTrailer.check()) { ++trailer; break; } 

    if ( ESTrailer.lenght() != nWords) {
      if (debug_) edm::LogWarning("Invalid Data")<<"Invalid ES data : the length is not correct !";
      return;
    }
    if ( ESTrailer.lenght() < 8) {
      if (debug_) edm::LogWarning("Invalid Data")<<"Invalid ES data : the length is not correct !";
      return;
    }

    if (debug_)  {
      cout<<"[ESUnpacker]: FED Trailer candidate. Is trailer? "<<ESTrailer.check();
      if (ESTrailer.check())
        cout<<". Length of the ES event: "<<ESTrailer.lenght()<<endl;
      else cout<<" WARNING!, this is not a ES Trailer"<<endl;
    }

    moreTrailers = ESTrailer.moreTrailers();
  }

  // Check ES data format version
  int vmajor = 0;
  int vminor = 0;
  int dccLineCount = 0;
  for (const Word64* word=(header+1); word!=(header+dccWords+1); ++word) {
    dccLineCount++;
    if (dccLineCount==3) {
      vminor = (*word >> 40) & 0x00ff;
      vmajor = (*word >> 48) & 0x00ff;
    }
  }
  //cout<<"ES version : "<<vmajor<<" "<<vminor<<endl;
  // DCC data
  int dccHeaderCount = 0;
  dccLineCount = 0;
  int dccHead, dccLine;
  Word64 m4  = ~(~Word64(0) << 4);
  for (const Word64* word=(header+1); word!=(header+dccWords+1); ++word) {
    if (debug_) cout<<"DCC   : "<<print(*word)<<endl;
    dccHead = (*word >> 60) & m4;
    if (dccHead == 3) dccHeaderCount++;
    dccLine = (*word >> 56) & m4;
    dccLineCount++;
    if (dccLine != dccLineCount && vmajor != 1 && vminor !=1) {
      if (debug_) edm::LogWarning("Invalid Data")<<"Invalid ES data : DCC header order is not correct !";
      return; 
    }
  }
  if (dccHeaderCount != 6 && vmajor != 1 && vminor !=1) {
    if (debug_) edm::LogWarning("Invalid Data")<<"Invalid ES data : DCC header lines are "<<dccHeaderCount;
    return;
  }
  
  // Event data
  static const Word64 mHEAD = ~(~Word64(0) << bHEAD);
  static const Word64 mKID  = ~(~Word64(0) << bKID);
  int head;
  int kchip = 0;
  for (const Word64* word=(header+dccWords+1); word!=trailer; ++word) {
    if (debug_) cout<<"Event : "<<print(*word)<<endl;

    head = (*word >> sHEAD) & mHEAD;
    if (head == 1) 
      kchip = (*word >> sKID) & mKID;
    else if (head == 2) 
      word2digi(kchip, *word, digis);
    
  }

}

void ESUnpacker::word2digi(int kchip, const Word64 & word, ESDigiCollection & digis) 
{

  //static const Word64 mHEAD  = ~(~Word64(0) << bHEAD);
  //static const Word64 mE1    = ~(~Word64(0) << bE1);
  //static const Word64 mE0    = ~(~Word64(0) << bE0);
  static const Word64 mSTRIP = ~(~Word64(0) << bSTRIP);
  static const Word64 mPACE  = ~(~Word64(0) << bPACE);
  static const Word64 mADC2  = ~(~Word64(0) << bADC2);
  static const Word64 mADC1  = ~(~Word64(0) << bADC1);
  static const Word64 mADC0  = ~(~Word64(0) << bADC0);
  
  //static const Word64 mKHEAD2  = ~(~Word64(0) << bKHEAD2);
  //static const Word64 mKHEAD1  = ~(~Word64(0) << bKHEAD1);
  //static const Word64 mFIBER   = ~(~Word64(0) << bFIBER);
  //static const Word64 mKID     = ~(~Word64(0) << bKID);
  //static const Word64 mKCE     = ~(~Word64(0) << bKCE);
  //static const Word64 mKCRC    = ~(~Word64(0) << bKCRC);
  //static const Word64 mKET     = ~(~Word64(0) << bKET);
  //static const Word64 mKFLAG1  = ~(~Word64(0) << bKFLAG1);
  //static const Word64 mKBC     = ~(~Word64(0) << bKBC);
  //static const Word64 mKFLAG2  = ~(~Word64(0) << bKFLAG2);
  //static const Word64 mKEC     = ~(~Word64(0) << bKEC);

  int adc[3];
  int strip = (word >> sSTRIP) & mSTRIP;
  int pace  = (word >> sPACE) & mPACE;
  adc[2]    = (word >> sADC2) & mADC2;
  adc[1]    = (word >> sADC1) & mADC1;
  adc[0]    = (word >> sADC0) & mADC0;

  if (debug_) cout<<kchip<<" "<<strip<<" "<<pace<<" "<<adc[2]<<" "<<adc[1]<<" "<<adc[0]<<endl;

  int zside = 0;
  int plane = 0;
  int ix = 0;
  int iy = 0;

  if (kchip>=0 && kchip<400) {
    zside = 1; plane = 1;
  } else if (kchip>=400 && kchip<800) {
    zside = 1; plane = 2; kchip -= 400;
  } else if (kchip>=800 && kchip<1200) {
    zside = -1; plane = 1; kchip -= 800;
  } else if (kchip>=1200 && kchip<1600) {
    zside = -1; plane = 2; kchip -= 1200;
  }

  ix = kchip % 20;
  iy = kchip / 20;

  if (pace == 0) {
    ix = (ix+1)*2-1;
    iy = (iy+1)*2-1;
  } else if (pace == 1) {
    ix = (ix+1)*2;
    iy = (iy+1)*2-1;
  } else if (pace == 2) {
    ix = (ix+1)*2-1;
    iy = (iy+1)*2; 
  } else if (pace == 3) {
    ix = (ix+1)*2;
    iy = (iy+1)*2;
  }

  if (debug_) cout<<"DetId : "<<zside<<" "<<plane<<" "<<ix<<" "<<iy<<" "<<strip+1<<endl;

  ESDetId detId(strip+1, ix, iy, plane, zside);
  ESDataFrame df(detId);
  df.setSize(3);

  for (int i=0; i<3; i++) df.setSample(i, adc[i]);  

  digis.push_back(df);

  if (debug_) 
    cout<<"Si : "<<detId.zside()<<" "<<detId.plane()<<" "<<detId.six()<<" "<<detId.siy()<<" "<<detId.strip()<<" ("<<kchip<<","<<pace<<") "<<df.sample(0).adc()<<" "<<df.sample(1).adc()<<" "<<df.sample(2).adc()<<endl;

}

string ESUnpacker::print(const  Word64 & word) const
{
  ostringstream str;
  str << "Word64:  " << reinterpret_cast<const bitset<64>&> (word);
  return str.str();
}

