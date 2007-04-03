#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEmCand.h"

#include <iostream>

#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

using std::ostream;
using std::string;
using std::hex;
using std::dec;

// default constructor
L1GctEmCand::L1GctEmCand() :
  m_data(0),
  m_iso(false),
  m_source(0)
{ 

}

// construct from raw data, no source - used in GT
L1GctEmCand::L1GctEmCand(uint16_t data, bool iso) :
  m_data(data),
  m_iso(iso),
  m_source(0),
  m_bx(0)
 {

 }

// construct from raw data with source - used in GCT unpacker
 L1GctEmCand::L1GctEmCand(uint16_t data, bool iso, uint16_t block, uint16_t index, int16_t bx) :
   m_data(data),
   m_iso(iso),
   m_source( ((block&0x7f)<<9) + (index&0x1ff) ),
   m_bx(bx)
 {

 }

// construct from content - used in GCT emulator
// eta = -6 to -0, +0 to +6. Sign is bit 3, 1 means -ve Z, 0 means +ve Z
L1GctEmCand::L1GctEmCand(unsigned rank, unsigned eta, unsigned phi, bool iso) : 
  m_iso(iso),
  m_bx(0)
 
{
  construct(rank, eta, phi);
  m_source=0;
}

// construct from content, with source - will be used in GCT emulator one day?
// eta = -6 to -0, +0 to +6. Sign is bit 3, 1 means -ve Z, 0 means +ve Z
L1GctEmCand::L1GctEmCand(unsigned rank, unsigned eta, unsigned phi, bool iso, uint16_t block, uint16_t index, int16_t bx) : 
  m_iso(iso),
  m_source( ((block&0x7f)<<9) + (index&0x1ff) ),
  m_bx(bx)
{
  construct(rank, eta, phi);
}

// construct from RCT output candidate
L1GctEmCand::L1GctEmCand(L1CaloEmCand& c) :
  m_iso(c.isolated()),
  m_source(0),
  m_bx(c.bx())
{
  construct(c.rank(), c.regionId().gctEta(), c.regionId().gctPhi());
}

// destructor
L1GctEmCand::~L1GctEmCand() { } 

// name of candidate type
string L1GctEmCand::name() const {
  return (isolated() ? "iso EM" : "non iso EM" ); 
}

// was a candidate found
bool L1GctEmCand::empty() const { 
  return (rank() == 0); 
}

// pretty print
ostream& operator<<(ostream& s, const L1GctEmCand& cand) {
  s << "L1GctEmCand : ";
  s << "rank=" << hex << cand.rank();
  s << ", etaSign=" << cand.etaSign() << ", ieta=" << (cand.etaIndex()&0x7) << ", iphi=" << cand.phiIndex();
  s << ", iso=" << cand.isolated() << dec;
  s << hex << " cap block=" << cand.capBlock() << ", index=" << cand.capIndex() << ", BX=" << cand.bx() << dec;
  return s;
}

// return region object
L1CaloRegionDetId L1GctEmCand::regionId() const {
  // get global eta
  unsigned eta = ( etaSign()==1 ? 10-(etaIndex()&0x7) : 11+(etaIndex()&0x7) );
  return L1CaloRegionDetId(eta, phiIndex());
}

// construct from rank, eta, phi
void L1GctEmCand::construct(unsigned rank, unsigned eta, unsigned phi) {
  m_data = (rank & 0x3f) + ((eta & 0xf)<<6) + ((phi & 0x1f)<<10);
}
