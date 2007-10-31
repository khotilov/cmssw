#ifndef DEUTILS_H
#define DEUTILS_H

/*\class template DEutils
 *\description data|emulation auxiliary template
               collection operations struct
 *\author Nuno Leonardo (CERN)
 *\date 07.04
 */

#include "L1Trigger/HardwareValidation/interface/DEtrait.h"

template <typename T> 
struct DEutils {

  typedef typename T::size_type col_sz;
  typedef typename T::const_iterator col_cit;
  typedef typename T::iterator col_it;
  typedef DEtrait<T> de_trait;
  typedef typename de_trait::cand_type cand_type;
  typedef typename de_trait::coll_type coll_type;

  public:
  
  DEutils() {
    if(de_type()>16)
      throw cms::Exception("ERROR") 
	<< "DEutils::DEutils() :: "
	<< "specialization is still missing for collection of type:" 
	<< de_type() << std::endl;
  }
  ~DEutils(){}
  
  inline int de_type() const {return de_trait::de_type();}
  bool   de_equal     (const cand_type&, const cand_type&);
  bool   de_equal_loc (const cand_type&, const cand_type&);
  bool   de_nequal    (const cand_type&, const cand_type&);
  bool   de_nequal_loc(const cand_type&, const cand_type&);
  col_it de_find      ( col_it, col_it,  const cand_type&);
  //col_it de_find_loc  ( col_it, col_it,  const cand_type&);

  std::string print(col_cit) const;
  bool is_empty(col_cit) const;
  std::string GetName(int) const;

  L1DataEmulDigi DEDigi(col_cit itd, col_cit itm, int ctype);
  
};


/// --- form de-digi ---

template <typename T> 
L1DataEmulDigi DEutils<T>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  ///return empty digi by default
  return L1DataEmulDigi();
}

template<> inline L1DataEmulDigi 
DEutils<EcalTrigPrimDigiCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  //fill data if flagged, otherwise emulator
  double x1 = (aflag!=4) ? itd->id().iphi() : itm->id().iphi();
  double x2 = (aflag!=4) ? itd->id().ieta() : itm->id().ieta();
  L1DataEmulDigi digi(dedefs::ETP,cid, x1,x2,0, errt);
  unsigned int dw = (aflag==4)?0:itd->sample(itd->sampleOfInterest()).raw();
  unsigned int ew = (aflag==3)?0:itm->sample(itm->sampleOfInterest()).raw();
  //also available: uint32_t id().rawId() ... merge words ?
  dw &= 0x0fff; ew &= 0x0fff; //12-bit
  digi.setData(dw,ew);
  int de = (aflag==4)?0:itd->compressedEt() ;
  int ee = (aflag==3)?0:itm->compressedEt() ;
  digi.setRank((float)de,(float)ee);
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<HcalTrigPrimDigiCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1 = (aflag!=4) ? itd->id().iphi() : itm->id().iphi();
  double x2 = (aflag!=4) ? itd->id().ieta() : itm->id().ieta();
  L1DataEmulDigi digi(dedefs::HTP,cid, x1,x2,0, errt);
  unsigned int dw = (aflag==4)?0:itd->t0().raw();
  unsigned int ew = (aflag==3)?0:itm->t0().raw();
  dw &= 0xf9ff; ew &= 0xf9ff; //16-bit (bits 10, 11 not set !?)
  digi.setData(dw,ew); 
  int de = (aflag==4)?0:itd->SOI_compressedEt();
  int ee = (aflag==3)?0:itm->SOI_compressedEt();
  digi.setRank((float)de,(float)ee);
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<L1CaloEmCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1, x2, x3(0.);
  // global index ieta (0-21), iphi (0-17)
  x1 = (aflag!=4) ? itd->regionId().iphi() : itm->regionId().iphi();
  x2 = (aflag!=4) ? itd->regionId().ieta() : itm->regionId().ieta();
  //alternative coordinates: rctCrate(), rctCard(), index()
  L1DataEmulDigi digi(dedefs::RCT,cid, x1,x2,x3, errt);
  unsigned int dw = (aflag==4)?0:itd->raw();
  unsigned int ew = (aflag==3)?0:itm->raw();
  dw &= 0x03ff; ew &= 0x03ff; //10-bit
  digi.setData(dw,ew); 
  int de = (aflag==4)?0:itd->rank();
  int ee = (aflag==3)?0:itm->rank();
  digi.setRank((float)de,(float)ee);
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<L1CaloRegionCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1 = (aflag!=4) ? itd->rctCrate() : itm->rctCrate(); //rctPhi()
  double x2 = (aflag!=4) ? itd->rctCard () : itm->rctCard (); //rctEta()
  double x3 = (aflag!=4) ? itd->rctRegionIndex() : itm->rctRegionIndex();
  L1DataEmulDigi digi(dedefs::RCT,cid, x1,x2,x3, errt);
  //Missing raw data accessor! (pack method private;)
  //dw &= 0x2fff; de &= 0x2fff; //14-bit 
  int de = (aflag==4)?0:itd->et();
  int ee = (aflag==3)?0:itm->et();
  digi.setRank((float)de,(float)ee);
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<L1GctEmCandCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  //phi: 0..17; eta: -6..-0,+0..+6; eta sign:1(z-),0(z+)
  // bring it to global coordinates 0..21 below
  double x1 = (aflag!=4) ? itd->phiIndex() : itm->phiIndex();
  unsigned deta = (itd->etaSign()==1 ? 10-(itd->etaIndex()&0x7) : 11+(itd->etaIndex()&0x7) );
  unsigned eeta = (itm->etaSign()==1 ? 10-(itm->etaIndex()&0x7) : 11+(itm->etaIndex()&0x7) );
  double x2 = (aflag!=4) ? deta : eeta;
  L1DataEmulDigi digi(dedefs::GCT,cid, x1,x2,0., errt);
  unsigned int dw = (aflag==4)?0:itd->raw();
  unsigned int ew = (aflag==3)?0:itm->raw();
  dw &= 0x7fff; ew &= 0x7fff; //15-bit
  digi.setData(dw,ew); 
  int de = (aflag==4)?0:itd->rank();
  int ee = (aflag==3)?0:itm->rank();
  digi.setRank((float)de,(float)ee);
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<L1GctJetCandCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  //phi: 0..17; eta: -6..-0,+0..+6; eta sign:1(z-),0(z+)
  // bring it to global coordinates 0..21 below
  double x1 = (aflag!=4) ? itd->phiIndex() : itm->phiIndex();
  unsigned deta(0), eeta(0);
  if (!itd->isForward()) deta=(itd->etaSign()==1?10-(itd->etaIndex()&0x7):(itd->etaIndex()&0x7)+11);
  else                   deta=(itd->etaSign()==1? 3-(itd->etaIndex()&0x7):(itd->etaIndex()&0x7)+18 );
  if (!itm->isForward()) eeta=(itm->etaSign()==1?10-(itm->etaIndex()&0x7):(itm->etaIndex()&0x7)+11);
  else                   eeta=(itm->etaSign()==1? 3-(itm->etaIndex()&0x7):(itm->etaIndex()&0x7)+18 );
  double x2 = (aflag!=4) ? deta : eeta;
  L1DataEmulDigi digi(dedefs::GCT,cid, x1,x2,0., errt);
  unsigned int dw = (aflag==4)?0:itd->raw();
  unsigned int ew = (aflag==3)?0:itm->raw();
  dw &= 0x7fff; ew &= 0x7fff; //15-bit
  digi.setData(dw,ew); 
  int de = (aflag==4)?0:itd->rank();
  int ee = (aflag==3)?0:itm->rank();
  digi.setRank((float)de,(float)ee);
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<L1MuRegionalCandCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int sid;
  switch(itd->type_idx()) { // 0 DT, 1 bRPC, 2 CSC, 3 fRPC
  case 0:  sid=dedefs::DTF; break;
  case 1:  sid=dedefs::RPC; break;
  case 2:  sid=dedefs::CTF; break;
  case 3:  sid=dedefs::RPC; break;
  default: sid=-1;
  }
  int cid = de_type();
  int errt = aflag;
  double x1 = (aflag!=4) ? itd->phiValue() : itm->phiValue();
  double x2 = (aflag!=4) ? itd->etaValue() : itm->etaValue();
  L1DataEmulDigi digi(sid,cid, x1,x2,0, errt);
  unsigned int dw = (aflag==4)?0 : itd->getDataWord();
  unsigned int ew = (aflag==3)?0 : itm->getDataWord();
  dw &= 0xffffffff; ew &= 0xffffffff; //32-bit
  digi.setData(dw,ew);
  int de = (aflag==4)?0:itd->pt_packed();//ptValue();
  int ee = (aflag==3)?0:itm->pt_packed();//ptValue();
  digi.setRank((float)de,(float)ee);
  //note: phi,eta,pt 'values' not always set for all muon tf systems
  //(under discussion) need universal mechanism for setting up physical units
  if(0) //check print
    std::cout << "L1DataEmulDigi DEutils<L1MuRegionalCandCollection>] dedigi info"
	      << " phivalue:" << itd->phiValue()   << "," << itm->phiValue()
	      << " etavalue:" << itd->etaValue()   << "," << itm->etaValue()
	      << " phipackd:" << itd->phi_packed() << "," << itm->phi_packed()
	      << " etapackd:" << itd->eta_packed() << "," << itm->eta_packed()
	      << std::endl;
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<L1MuGMTCandCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1 = (aflag!=4) ? itd->phiValue() : itm->phiValue();
  double x2 = (aflag!=4) ? itd->etaValue() : itm->etaValue();
  L1DataEmulDigi digi(dedefs::GMT,cid, x1,x2,0, errt);
  unsigned int dw = (aflag==4)?0 : itd->getDataWord();
  unsigned int ew = (aflag==3)?0 : itm->getDataWord();
  dw &= 0x3ffffff; ew &= 0x3ffffff; //26-bit
  digi.setData(dw,ew);
  int de = (aflag==4)?0:itd->ptIndex();//ptValue();
  int ee = (aflag==3)?0:itm->ptIndex();//ptValue();
  digi.setRank((float)de,(float)ee);
  if(0) //check print
  std::cout << "l1dataemuldigi l1mugmtcandcoll type:" << cid 
	    << " eta:" << itd->etaValue() << ", " << itm->etaValue()
	    << " phi:" << itd->phiValue() << ", " << itm->phiValue()
	    << std::hex << " word d:" << dw << "e:" << ew << std::dec 
	    << std::endl;
  return digi;
}

template<> 
inline L1DataEmulDigi DEutils<L1MuDTChambPhDigiCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1 = (aflag!=4) ? itd->scNum() : itm->scNum();
  double x2 = (aflag!=4) ? itd->whNum() : itm->whNum();
  double x3 = (aflag!=4) ? itd->stNum() : itm->stNum();
  L1DataEmulDigi digi(dedefs::DTP,cid, x1,x2,x3, errt);
  //other coordinate methods phi(), phiB()
  //note: no data word defined for candidate
  int dr = (aflag==4)?0:itd->code();
  int er = (aflag==3)?0:itm->code();
  digi.setRank((float)dr,(float)er);
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<L1MuDTChambThDigiCollection>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1 = (aflag!=4) ? itd->scNum() : itm->scNum();
  double x2 = (aflag!=4) ? itd->whNum() : itm->whNum();
  double x3 = (aflag!=4) ? itd->stNum() : itm->stNum();
  L1DataEmulDigi digi(dedefs::DTP,cid, x1,x2,x3, errt);
  //note: no data word defined for candidate
  int dr(0), er(0);
  for(int i=0; i<7;i++){
    if(itd->code(i)>=dr) dr=itd->quality(i);
    if(itm->code(i)>=er) er=itm->quality(i);
  }
  //alternatives: code() = quality() + positions()
  dr = (aflag==4)?0:dr;
  er = (aflag==3)?0:er;
  digi.setRank((float)dr,(float)er);
  return digi;
}


template<> inline L1DataEmulDigi 
DEutils<CSCCorrelatedLCTDigiCollection_>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1 = (aflag!=4) ? itd->getTrknmb() : itm->getTrknmb();
  double x2 = (aflag!=4) ? itd->getKeyWG () : itm->getKeyWG ();
  //multiple subsystem ctp,ctf
  L1DataEmulDigi digi(-1,cid, x1,x2,0, errt);
  //note: no data word and rank defined for candidate
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<CSCALCTDigiCollection_>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1 = (aflag!=4) ? itd->getTrknmb() : itm->getTrknmb();
  double x2 = (aflag!=4) ? itd->getKeyWG () : itm->getKeyWG ();
  L1DataEmulDigi digi(dedefs::CTP,cid, x1,x2,0, errt);
  //note: no data word and rank defined for candidate
  return digi;
}
template<> inline L1DataEmulDigi 
DEutils<CSCCLCTDigiCollection_>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1 = (aflag!=4) ? itd->getTrknmb() : itm->getTrknmb();
  L1DataEmulDigi digi(dedefs::CTP,cid, x1,0,0, errt);
  //note: no data word and rank defined for candidate
  return digi;
}

template<> inline L1DataEmulDigi 
DEutils<L1CSCSPStatusDigiCollection_>::DEDigi(col_cit itd,  col_cit itm, int aflag) {
  int cid = de_type();
  int errt = aflag;
  double x1; //sector/slot
  x1 = (aflag!=4) ? itd->slot() : itm->slot();
  //sector-slot map to be in principle to be provided from event setup
  //int de_cscstatus_slot2sector[22] = 
  // {0,0,0,0,0, 0,1,2,3,4, 5,6,0,0,0, 0,7,8,9,10,  11,12};
  //x1 = (aflag!=4) ? slot2sector[itd->slot()] : slot2sector[itm->slot()];
  L1DataEmulDigi digi(dedefs::CTF,cid, x1,0,0, errt);
  //note: no data word and rank defined for candidate
  return digi;
}

/// --- find candidate ---

template <typename T> typename 
DEutils<T>::col_it DEutils<T>::de_find( col_it first, col_it last, const cand_type& value ) {
  for ( ;first!=last; first++) 
    if ( de_equal(*first,value) ) break;
  return first;
}

/*
template <typename T> typename 
DEutils<T>::col_it DEutils<T>::de_find_loc( col_it first, col_it last, const cand_type& value ) {
  for ( ;first!=last; first++) 
    if ( de_equal_loc(*first,value) ) break;
  return first;
}
*/

/// --- candidate match definition ---

template <typename T>
bool DEutils<T>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  //declare candidate matching by default
  return true;
}
template <typename T>
bool DEutils<T>::de_nequal(const cand_type& lhs, const cand_type& rhs) {
  return !de_equal(lhs,rhs);
}

template <> inline bool 
DEutils<EcalTrigPrimDigiCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs[lhs.sampleOfInterest()].raw() == rhs[rhs.sampleOfInterest()].raw());
  val &= (lhs.id().rawId()                  == rhs.id().rawId());
  return val;
}

template <> inline bool 
DEutils<HcalTrigPrimDigiCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.t0().raw()     == rhs.t0().raw());
  val &= (lhs.id().rawId()   == rhs.id().rawId());
  return val;
}

template <> inline bool 
DEutils<L1CaloEmCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.raw()      == rhs.raw()     );
  val &= (lhs.rctCrate() == rhs.rctCrate());
  val &= (lhs.isolated() == rhs.isolated());
  val &= (lhs.index()    == rhs.index()   );
  //val &= (lhs.bx()       == rhs.bx()      );
  return val;
}

template <> inline bool 
DEutils<L1CaloRegionCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.et()        == rhs.et()       );
  val &= (lhs.rctCrate()  == rhs.rctCrate() );	
  val &= (lhs.rctRegionIndex() == rhs.rctRegionIndex());
  val &= (lhs.id().isHf() == rhs.id().isHf());  
  if (!lhs.id().isHf()){
    val &= (lhs.overFlow()  == rhs.overFlow() );
    val &= (lhs.tauVeto()   == rhs.tauVeto()  );
    val &= (lhs.mip()       == rhs.mip()      );
    val &= (lhs.quiet()     == rhs.quiet()    );
    val &= (lhs.rctCard()   == rhs.rctCard()  );
  } else {
    val &= (lhs.fineGrain() == rhs.fineGrain());
  }
  return val;
}

template <> inline bool 
DEutils<L1GctEmCandCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  return lhs==rhs;
}

template <> inline bool 
DEutils<L1GctJetCandCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  return lhs==rhs;
}

template <> inline bool 
DEutils<L1MuDTChambPhDigiCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.whNum() ==rhs.whNum() );
  val &= (lhs.scNum() ==rhs.scNum() );
  val &= (lhs.stNum() ==rhs.stNum() );
  val &= (lhs.phi()   ==rhs.phi()   );
  val &= (lhs.phiB()  ==rhs.phiB()  );
  val &= (lhs.code()  ==rhs.code()  );
  val &= (lhs.Ts2Tag()==rhs.Ts2Tag());
  //val &= (lhs.BxCnt() ==rhs.BxCnt() ); 
  //val &= (lhs.bxNum() ==rhs.bxNum() );
  return val;
}

template <> inline bool 
DEutils<L1MuDTChambThDigiCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.whNum() ==rhs.whNum() );
  val &= (lhs.scNum() ==rhs.scNum() );
  val &= (lhs.stNum() ==rhs.stNum() );
  //for(int i=0; i<7; i++) {
  //  val &= (lhs.code(i)    ==rhs.code(i)    );
  //  val &= (lhs.position(i)==rhs.position(i));
  //  val &= (lhs.quality(i) ==rhs.quality(i) );
  //}
  //val &= (lhs.bxNum() ==rhs.bxNum() );
  return val;
}

template <> inline bool 
DEutils<L1MuRegionalCandCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.getDataWord() == rhs.getDataWord() );
  //check whether collections being compared refer to same system and bx!
  val &= (lhs.type_idx() == rhs.type_idx());
  val &= (lhs.bx()       == rhs.bx());
  return val;
}

template <> inline bool 
DEutils<L1MuGMTCandCollection>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.getDataWord() == rhs.getDataWord() );
  return val;
  //return lhs==rhs; //(dataword,bx..)
  }

template <> inline bool 
DEutils<CSCCorrelatedLCTDigiCollection_>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.getTrknmb()  == rhs.getTrknmb() );
  val &= (lhs.isValid()    == rhs.isValid()   );
  val &= (lhs.getQuality() == rhs.getQuality());
  val &= (lhs.getKeyWG()   == rhs.getKeyWG()  );
  val &= (lhs.getStrip()   == rhs.getStrip()  );
  val &= (lhs.getPattern() == rhs.getPattern());
  val &= (lhs.getBend()    == rhs.getBend()   );
  val &= (lhs.getMPCLink() == rhs.getMPCLink()); 
  val &= (lhs.getBX()      == rhs.getBX()     );    
  return val;
  //return lhs==rhs;
}
template <> inline bool 
DEutils<CSCALCTDigiCollection_>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  return lhs==rhs;
}
template <> inline bool 
DEutils<CSCCLCTDigiCollection_>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  return lhs==rhs;
}
template <> inline bool 
DEutils<L1CSCSPStatusDigiCollection_>::de_equal(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.slot() == rhs.slot());
  val &= (lhs.BXN () == rhs.BXN ());
  val &= (lhs.FMM () == rhs.FMM ());
  val &= (lhs.SEs () == rhs.SEs ());
  val &= (lhs.SMs () == rhs.SMs ());
  val &= (lhs.BXs () == rhs.BXs ());
  val &= (lhs.AFs () == rhs.AFs ());
  val &= (lhs.VPs () == rhs.VPs ());
  return val;
}

/// --- candidate location-match definition ---

template <typename T>
bool DEutils<T>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  //declare candidate matching by default
  return true;
}
template <typename T>
bool DEutils<T>::de_nequal_loc(const cand_type& lhs, const cand_type& rhs) {
  return !de_equal_loc(lhs,rhs);
}


template <> inline bool 
DEutils<EcalTrigPrimDigiCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.id().zside()   == rhs.id().zside()  );
  val &= (lhs.id().ietaAbs() == rhs.id().ietaAbs());
  val &= (lhs.id().iphi()    == rhs.id().iphi()   );
  return val;
}

template <> inline bool 
DEutils<HcalTrigPrimDigiCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.id().zside()   == rhs.id().zside()  );
  val &= (lhs.id().ietaAbs() == rhs.id().ietaAbs());
  val &= (lhs.id().iphi()    == rhs.id().iphi()   );
  return val;
}

template <> inline bool 
DEutils<L1CaloEmCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.rctCrate()  == rhs.rctCrate());
  val &= (lhs.rctCard()   == rhs.rctCard());
  val &= (lhs.rctRegion() == rhs.rctRegion());
  return val;
}

template <> inline bool 
DEutils<L1CaloRegionCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.rctCrate()  == rhs.rctCrate() );	
  val &= (lhs.id().isHf() == rhs.id().isHf());  
  if (!lhs.id().isHf())
    val &= (lhs.rctCard() == rhs.rctCard()  );
  val &= (lhs.rctRegionIndex() == rhs.rctRegionIndex());
  return val;
}

template <> inline bool 
DEutils<L1GctEmCandCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.etaIndex() == rhs.etaIndex());
  val &= (lhs.phiIndex() == rhs.phiIndex());
  return val;
}
template <> inline bool 
DEutils<L1GctJetCandCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.etaIndex() == rhs.etaIndex());
  val &= (lhs.phiIndex() == rhs.phiIndex());
  return val;
}

template <> inline bool 
DEutils<L1MuRegionalCandCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.phi_packed() ==rhs.phi_packed() );
  val &= (lhs.eta_packed() ==rhs.eta_packed() );
  //val &= (lhs.type_idx() == rhs.type_idx());
  //val &= (lhs.bx()       == rhs.bx());
  return val;
}

template <> inline bool 
DEutils<L1MuGMTCandCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.phiIndex() ==rhs.phiIndex() );
  val &= (lhs.etaIndex() ==rhs.etaIndex() );
  return val;
}

template <> inline bool 
DEutils<L1MuDTChambPhDigiCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.whNum() ==rhs.whNum() );
  val &= (lhs.scNum() ==rhs.scNum() );
  val &= (lhs.stNum() ==rhs.stNum() );
  //val &= (lhs.phi()   ==rhs.phi()   );
  //val &= (lhs.phiB()  ==rhs.phiB()  );
  //val &= (lhs.bxNum() ==rhs.bxNum() );
  return val;
}

template <> inline bool 
DEutils<L1MuDTChambThDigiCollection>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.whNum() ==rhs.whNum() );
  val &= (lhs.scNum() ==rhs.scNum() );
  val &= (lhs.stNum() ==rhs.stNum() );
  //val &= (lhs.bxNum() ==rhs.bxNum() );
  return val;
}

template <> inline bool 
DEutils<CSCCorrelatedLCTDigiCollection_>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.getTrknmb()  == rhs.getTrknmb() );
  val &= (lhs.getKeyWG()   == rhs.getKeyWG()  );
  return val;
}

template <> inline bool 
DEutils<CSCALCTDigiCollection_>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.getTrknmb()  == rhs.getTrknmb() );
  val &= (lhs.getKeyWG()   == rhs.getKeyWG()  );
  return val;
}
template <> inline bool 
DEutils<CSCCLCTDigiCollection_>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.getTrknmb()  == rhs.getTrknmb() );
  return val;
}
template <> inline bool 
DEutils<L1CSCSPStatusDigiCollection_>::de_equal_loc(const cand_type& lhs, const cand_type& rhs) {
  bool val = true;
  val &= (lhs.slot() == rhs.slot());
  return val;
}
 
/// --- candidate emptiness definition ---

template <typename T> 
bool DEutils<T>::is_empty(col_cit it) const { 
  //declare candidate non-empty by default
  return false; 
}

template<>
inline bool DEutils<EcalTrigPrimDigiCollection>::is_empty(col_cit it) const { 
  return ( it->size()==0 || it->sample(it->sampleOfInterest()).raw()==0);
}

template<>
inline bool DEutils<HcalTrigPrimDigiCollection>::is_empty(col_cit it) const { 
  return (  it->size()==0 || it->t0().raw()==0 || it->SOI_compressedEt()==0 );
}

template<>
inline bool DEutils<L1CaloEmCollection>::is_empty(col_cit it) const { 
    return  ((it->rank())==0);
    //return  ((it->raw())==0);
}

template<>
inline bool DEutils<L1CaloRegionCollection>::is_empty(col_cit it) const { 
    return  ((it->et())==0);
    //note: missing accessors in dataformats
}

template<>
inline bool DEutils<L1GctEmCandCollection>::is_empty(col_cit it) const { 
  return (it->empty());
}

template<>
inline bool DEutils<L1GctJetCandCollection>::is_empty(col_cit it) const { 
    return  (it->empty());
}

template<>
inline bool DEutils<L1MuDTChambPhDigiCollection>::is_empty(col_cit it) const { 
  return (it->code() == 7); 
  //return (it->qualityCode() == 7); 
  //return  false;
}
template<>
inline bool DEutils<L1MuDTChambThDigiCollection>::is_empty(col_cit it) const { 
  return (it->whNum()==0 && it->scNum()==0 && it->stNum()==0);//tmp!
  //return  false;
}

template<>
inline bool DEutils<L1MuRegionalCandCollection>::is_empty(col_cit it) const { 
  //note: following call used to give trouble sometimes
  return (it->empty()); 
  //virtual bool empty() const { return readDataField( PT_START, PT_LENGTH) == 0; }
  //return  (it->getDataWord()==0);
  //return  (it->pt_packed()==0);
}

template<>
inline bool DEutils<L1MuGMTCandCollection>::is_empty(col_cit it) const { 
  return (it->empty());
  //return (it->ptIndex()==0);
  //return  (it->getDataWord()==0);
}

template<>
inline bool DEutils<CSCCorrelatedLCTDigiCollection_>::is_empty(col_cit it) const { 
  return !(it->isValid());
}
template<>
inline bool DEutils<CSCALCTDigiCollection_>::is_empty(col_cit it) const { 
  return !(it->isValid());
}
template<>
inline bool DEutils<CSCCLCTDigiCollection_>::is_empty(col_cit it) const { 
  return !(it->isValid());
}

template<>
inline bool DEutils<L1CSCSPStatusDigiCollection_>::is_empty(col_cit it) const { 
  unsigned data = 
    it->slot() | it->BXN () | it->FMM () | it->SEs () |
    it->SMs () | it->BXs () | it->AFs () | it->VPs ();
  return data==0;
}

/// --- print candidate ---

template <typename T> 
std::string DEutils<T>::print(col_cit it) const {
  std::stringstream ss;
  ss << "[DEutils<T>::print()] specialization still missing for collection!";
  //ss << *it; // default
  ss << std::endl;
  return ss.str();
}

template <> 
inline std::string DEutils<EcalTrigPrimDigiCollection>::print(col_cit it) const {
  std::stringstream ss;
  ss << "0x" << std::setw(4) << std::setfill('0') << std::hex 
     << it->sample(it->sampleOfInterest()).raw()
     << std::setfill(' ') << std::dec 
     << ", et:"   << std::setw(3) << it->compressedEt() 
     << ", fg:"   << std::setw(1) << it->fineGrain()
     << ", ttf:"  << std::setw(2) << it->ttFlag()
     << ", sdet:" << ((it->id().subDet()==EcalBarrel)?("Barrel"):("Endcap")) 
     << ", iz:"   << ((it->id().zside()>0)?("+"):("-")) 
     << ", ieta:" << std::setw(2) << it->id().ietaAbs()
     << ", iphi:" << std::setw(2) << it->id().iphi()
    //<< "\n\t: " << *it 
     << std::endl;
  return ss.str();
}

template <> 
inline std::string DEutils<HcalTrigPrimDigiCollection>::print(col_cit it) const {
  std::stringstream ss;
  ss << "0x" << std::setw(4) << std::setfill('0') << std::hex 
     << it->t0().raw()
     << std::setfill(' ') << std::dec 
     << ", et:"   << std::setw(3) << it->SOI_compressedEt()
     << ", fg:"   << std::setw(1) << it->SOI_fineGrain()
     << ", sdet:" << it->id().subdet()
     << ", iz:"   << ((it->id().zside()>0)?("+"):("-")) 
     << ", ieta:" << std::setw(2) << it->id().ietaAbs()
     << ", iphi:" << std::setw(2) << it->id().iphi()
     << std::endl;
  return ss.str();
}

template <> 
inline std::string DEutils<L1CaloEmCollection>::print(col_cit it) const {
  std::stringstream ss;
  ss << "0x" << std::setw(4) << std::setfill('0') << std::hex << it->raw() 
     << ", rank=0x"<< std::setw(2) << std::hex << it->rank()
     << std::setfill(' ') << std::dec 
     << ", region="<< std::setw(1) << it->rctRegion() 
     << ", card="  << std::setw(1) << it->rctCard() 
     << ", crate=" << std::setw(2) << it->rctCrate()
     << ", iso="   << std::setw(1) << it->isolated()
     << ", index=" << std::setw(1) << it->index() 
     << ", bx="    << it->bx()
     << std::endl;
  //ss << *it;
  return ss.str();
}

template <> 
inline std::string DEutils<L1CaloRegionCollection>::print(col_cit it) const {
  std::stringstream ss;
  ss << *it;
  //note: raw() data accessor missing in dataformats
  return ss.str();
}

template <> 
inline std::string DEutils<L1GctEmCandCollection>::print(col_cit it) const {
  std::stringstream ss;
  ss << "0x" << std::setw(4) << std::setfill('0') << std::hex << it->raw() 
     << *it << std::dec << std::endl;
  return ss.str();
}
template <> 
inline std::string DEutils<L1GctJetCandCollection>::print(col_cit it) const {
  std::stringstream ss;
  ss << "0x" << std::setw(4) << std::setfill('0') << std::hex << it->raw() 
     << *it << std::dec; 
  return ss.str();
}

template <> 
inline std::string DEutils<L1MuDTChambPhDigiCollection>::print(col_cit it) const {
  std::stringstream ss;
  ss << ""
     << " bxNum:"  << it->bxNum()  
     << " whNum:"  << it->whNum()  
     << " scNum:"  << it->scNum()  
     << " stNum:"  << it->stNum()  
     << " phi:"    << it->phi()    
     << " phiB:"   << it->phiB()   
     << " code:"   << it->code()   
     << " Ts2Tag:" << it->Ts2Tag() 
     << " BxCnt:"  << it->BxCnt()  
     << std::endl;
  //nb: operator << not implemented in base class L1MuDTChambPhDigi
  return ss.str();
}

template <> 
inline std::string DEutils<L1MuDTChambThDigiCollection>::print(col_cit it) const {
  std::stringstream ss;
  ss << ""
     << " bxNum:"  << it->bxNum()  
     << " whNum:"  << it->whNum()  
     << " scNum:"  << it->scNum()  
     << " stNum:"  << it->stNum()  
     << std::endl;
  //nb: operator << not implemented in base class L1MuDTChambThDigi
  return ss.str();
}

template <> 
inline std::string DEutils<L1MuRegionalCandCollection>::print(col_cit it) const {
  std::stringstream ss;
  const float noval = -10; //L1MuRegionalCand::m_invalidValue;
  ss << std::setiosflags(std::ios::showpoint | std::ios::fixed | std::ios::right | std::ios::adjustfield);
  ss   << std::hex << std::setfill('0')    
       << " 0x"    << std::setw(8) << it->getDataWord();
  if(it->phiValue()==noval || it->etaValue()==noval || it->ptValue()==noval )
    ss << std::hex << std::setfill('0')    
       << " pt:0x" << std::setw(2) << it->pt_packed() 
       << " phi:0x"<< std::setw(2) << it->phi_packed()
       << " eta:0x"<< std::setw(2) << it->eta_packed();
  else
    ss << std::dec << std::setfill(' ')
       << " pt:"   << std::setw(5) << std::setprecision(1) << it->ptValue() <<"[GeV]"
       << " phi:"  << std::setw(5) << std::setprecision(3) << it->phiValue()<<"[rad]"
       << " eta:"  << std::setw(6) << std::setprecision(3) << it->etaValue(); 
  ss   << std::dec << std::setfill(' ')
       << " qua:"  << std::setw(1) << it->quality() 
       << " cha:"  << std::setw(2) << it->chargeValue() 
       << " chav:" << std::setw(1) << it->chargeValid() 
       << " fh:"   << std::setw(1) << it->isFineHalo() 
       << " bx:"   << std::setw(4) << it->bx() 
       << " [id:"  << std::setw(1) << it->type_idx() << "]" // 0 DT, 1 bRPC, 2 CSC, 3 fRPC
       << std::endl;  
  //ss << it->print() 
  return ss.str();
}

template <> 
inline std::string DEutils<L1MuGMTCandCollection>::print(col_cit it) const {
  std::stringstream ss;
  ss << std::setiosflags(std::ios::showpoint | std::ios::fixed | std::ios::right | std::ios::adjustfield);
  const float noval = -10; //L1MuGMTCand::m_invalidValue;
  ss   << std::hex << std::setfill('0')    
       << " 0x"    << std::setw(7) << it->getDataWord();
  if(it->phiValue()==noval || it->etaValue()==noval || it->ptValue()==noval)
    ss << std::hex << std::setfill('0')    
       << " pt:0x" << std::setw(2) << it->ptIndex()
       << " eta:0x"<< std::setw(2) << it->etaIndex()
       << " phi:0x"<< std::setw(3) << it->phiIndex();
  else
    ss << std::dec << std::setfill(' ')    
       << " pt:"   << std::setw(5) << std::setprecision(1) << it->ptValue() <<"[GeV]"
       << " phi:"  << std::setw(5) << std::setprecision(3) << it->phiValue()<<"[rad]"
       << " eta:"  << std::setw(6) << std::setprecision(2) << it->etaValue();
  ss   << std::dec << std::setfill(' ')
       << " cha:"  << std::setw(2) << it->charge()  
       << " qua:"  << std::setw(3) << it->quality() 
       << " iso:"  << std::setw(1) << it->isol()    
       << " mip:"  << std::setw(1) << it->mip() 
       << " bx:"                   << it->bx() 
       << std::endl;
  //ss << it->print() 
  return ss.str();
}

template <> 
inline std::string DEutils<CSCCorrelatedLCTDigiCollection_>::print(col_cit it) const {
  std::stringstream ss;
  ss 
    << " ltc#:"     << it->getTrknmb()
    << " val:"      << it->isValid()
    << " qua:"      << it->getQuality() 
    << " mpc-link:" << it->getMPCLink()
    << " strip:"    << it->getStrip()
    << "("          << ((it->getStripType() == 0) ? 'D' : 'H') << ")"
    << " bend:"     << ((it->getBend() == 0) ? 'L' : 'R')
    << " patt:"     << it->getCLCTPattern()
    <<"  key-wire:" << it->getKeyWG()
    << " bx:"       << it->getBX()
    << std::endl;
  //ss << *it;
  return ss.str();
}

template <> 
inline std::string DEutils<CSCALCTDigiCollection_>::print(col_cit it) const {
  std::stringstream ss;
  ss
    << *it
    << std::endl;
  return ss.str();
}

template <> 
inline std::string DEutils<CSCCLCTDigiCollection_>::print(col_cit it) const {
  std::stringstream ss;
  ss 
    << *it
    << std::endl;
  return ss.str();
}

template <> 
inline std::string DEutils<L1CSCSPStatusDigiCollection_>::print(col_cit it) const {
  std::stringstream ss;
  ss 
    << " slot:"<< it->slot()
    << " bxn:" << it->BXN ()
    << " fmm:" << it->FMM ()
    << " ses:" << it->SEs ()
    << " sms:" << it->SMs ()
    << " bxs:" << it->BXs ()
    << " afs:" << it->AFs ()
    << " vps:" << it->VPs ()
    << std::endl;
  return ss.str();
}

/// --- name candidate ---

template <typename T> 
std::string DEutils<T>::GetName(int i=0) const {

  const int nlabel = 16;
  if(!(i<nlabel)) 
    return                  "un-defined" ;
  std::string str[nlabel]= {"un-registered"};

  switch(de_type()) {
  case dedefs::ECALtp:
    str[0] = "ECAL tp";
    str[1] = "EcalTrigPrimDigiCollection";
    str[2] = "EcalTriggerPrimitiveDigi";
  break;
  case dedefs::HCALtp:
    str[0] = "HCAL tp";
    str[1] = "HcalTrigPrimDigiCollection";
    str[2] = "HcalTriggerPrimitiveDigi";
  break;
  case dedefs::RCTem:
    str[0] = "RCT em";
    str[1] = "L1CaloEmCollection";
    str[2] = "L1CaloEmCand";
  break;
  case dedefs::RCTrgn:
    str[0] = "RCT region";
    str[1] = "L1CaloRegionCollection";
    str[2] = "L1CaloRegion";
    break;
  case dedefs::GCTem:
    str[0] = "GCT em";
    str[1] = "L1GctEmCandCollection";
    str[2] = "L1GctEmCand";
   break;
  case dedefs::GCTjet:
    str[0] = "GCT jet";
    str[1] = "L1GctJetCandCollection";
    str[2] = "L1GctJetCand";
   break;
  case dedefs::DTtpPh:
    str[0] = "DT tp phi";
    str[1] = "L1MuDTChambPhDigiCollection";
    str[2] = "L1MuDTChambPhDigi";
   break;
  case dedefs::DTtpTh:
    str[0] = "DT tp theta";
    str[1] = "L1MuDTChambThDigiCollection";
    str[2] = "L1MuDTChambThDigi";
   break;
  case dedefs::CSCtpa:
    str[0] = "CSC tpa";
    str[1] = "CSCALCTDigiCollection";
    str[2] = "CSCALCTDigi";
   break;
  case dedefs::CSCtpc:
    str[0] = "CSC tpc";
    str[1] = "CSCCLCTDigiCollection";
    str[2] = "CSCCLCTDigi";
   break;
  case dedefs::CSCtpl:
    str[0] = "CSC tp";
    str[1] = "CSCCorrelatedLCTDigiCollection";
    str[2] = "CSCCorrelatedLCTDigi";
   break;
  case dedefs::CSCsta:
    str[0] = "CSC tf status";
    str[1] = "L1CSCSPStatusDigiCollection_";
    str[2] = "L1CSCSPStatusDigi";
   break;
  case dedefs::MUrtf:
    str[0] = "Mu reg tf";
    str[1] = "L1MuRegionalCandCollection";
    str[2] = "L1MuRegionalCand";
   break;
  case dedefs::LTCi:
    str[0] = "LTC";
    str[1] = "LTCDigiCollection";
    str[2] = "LTCDigi";
    break;
  case dedefs::GMTcnd:
    str[0] = "GMT cand";
    str[1] = "L1MuGMTCandCollection";
    str[2] = "L1MuGMTCand";
    break;
  case dedefs::GMTrdt:
    str[0] = "GMT record";
    str[1] = "L1MuGMTReadoutRecordCollection";
    str[2] = "L1MuGMTReadoutRecord";
    break;
  case dedefs::GTdword:
    str[0] = "";
    str[1] = "";
    str[2] = "";
    break;
    //default:
  }
  return str[i];
}

/// --- order candidates ---

template <typename T>
struct de_rank : public DEutils<T> , public std::binary_function<typename DEutils<T>::cand_type, typename DEutils<T>::cand_type, bool> {
  typedef DEtrait<T> de_trait;
  typedef typename de_trait::cand_type cand_type;
  bool operator()(const cand_type& x, const cand_type& y) const {
    return true; //default
  }
};

template <> inline bool de_rank<EcalTrigPrimDigiCollection>::operator()(const cand_type& x, const cand_type& y) const { return x.compressedEt() > y.compressedEt(); }
template <> inline bool de_rank<HcalTrigPrimDigiCollection>::operator()(const cand_type& x, const cand_type& y) const { return x.SOI_compressedEt() > y.SOI_compressedEt(); }

template <> 
inline bool de_rank<L1CaloEmCollection>::operator() 
     (const cand_type& x, const cand_type& y) const {
  if       (x.rank()      != y.rank())     {
    return (x.rank()      <  y.rank())     ;
  } else if(x.isolated()  != y.isolated()) {
    return (x.isolated())?1:0;
  } else if(x.rctRegion() != y.rctRegion()){
    return (x.rctRegion() <  y.rctRegion());
  } else if(x.rctCrate()  != y.rctCrate()) {
    return (x.rctCrate()  <  y.rctCrate()) ;
  } else if(x.rctCard()   != y.rctCard())  {
    return (x.rctCard()   <  y.rctCard())  ;
  } else {
    return x.raw() < y.raw();
  }
}

template <> inline bool de_rank<L1CaloRegionCollection>::operator()(const cand_type& x, const cand_type& y) const { return x.et() < y.et(); }

template <> inline bool de_rank<L1GctEmCandCollection>::operator()(const cand_type& x, const cand_type& y)const { if(x.rank()!=y.rank()){return x.rank() < y.rank();} else{if(x.etaIndex()!=y.etaIndex()){return y.etaIndex() < x.etaIndex();}else{ return x.phiIndex() < y.phiIndex();}}}
template <> inline bool de_rank<L1GctJetCandCollection>::operator()(const cand_type& x, const cand_type& y)const { if(x.rank()!=y.rank()){return x.rank() < y.rank();} else{if(x.etaIndex()!=y.etaIndex()){return y.etaIndex() < x.etaIndex();}else{ return x.phiIndex() < y.phiIndex();}}}

template <> inline bool de_rank<L1MuDTChambPhDigiCollection>::operator()(const cand_type& x, const cand_type& y)const { if(x.whNum()!=y.whNum()){return x.whNum() < y.whNum();} else{if(x.scNum()!=y.scNum()){return y.scNum() < x.scNum();}else{ return x.stNum() < y.stNum();}}}
template <> inline bool de_rank<L1MuDTChambThDigiCollection>::operator()(const cand_type& x, const cand_type& y)const { if(x.whNum()!=y.whNum()){return x.whNum() < y.whNum();} else{if(x.scNum()!=y.scNum()){return y.scNum() < x.scNum();}else{ return x.stNum() < y.stNum();}}}

template <> inline bool de_rank<L1MuRegionalCandCollection>::operator()(const cand_type& x, const cand_type& y)const {if(x.phi_packed()!=y.phi_packed()){return x.phi_packed() < y.phi_packed();} else{if(x.eta_packed()!=y.eta_packed()){return y.eta_packed() < x.eta_packed();}else{ return x.quality_packed() < y.quality_packed();}}}

template <> inline bool de_rank<L1MuGMTCandCollection>::operator()(const cand_type& x, const cand_type& y)const {
  if(x.bx()!=y.bx()){return x.bx() < y.bx();} 
  else if(x.ptIndex()!=y.ptIndex()){return x.ptIndex() < y.ptIndex();}
  else{ return x.quality() < y.quality();}
}

template <> inline bool de_rank<CSCCorrelatedLCTDigiCollection_>::operator()(const cand_type& x, const cand_type& y)const {if(x.getTrknmb()!=y.getTrknmb()){return x.getTrknmb() < y.getTrknmb();} else{if(x.getKeyWG()!=y.getKeyWG()){return y.getKeyWG() < x.getKeyWG();} else{ return x.getQuality() < y.getQuality();}}}

#endif
