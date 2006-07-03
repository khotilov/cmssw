#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJet.h"

#include "FWCore/Utilities/interface/Exception.h"  

//DEFINE STATICS
const unsigned L1GctJet::RAWSUM_BITWIDTH = 10;
const unsigned L1GctJet::N_RGN_ETA = 22;
const unsigned L1GctJet::N_RGN_PHI = 18;


L1GctJet::L1GctJet(uint16_t rawsum, unsigned eta, unsigned phi, bool tauVeto,
		   L1GctJetEtCalibrationLut* lut) :
  m_rawsum(rawsum),
  m_id(eta, phi),
  m_tauVeto(tauVeto),
  m_jetEtCalibrationLut(lut)
{

}

L1GctJet::~L1GctJet()
{
}

std::ostream& operator << (std::ostream& os, const L1GctJet& cand)
{
  os << "L1 Gct jet";
  os << " energy sum " << cand.m_rawsum;
  os << " Eta " << cand.globalEta();
  os << " Phi " << cand.globalPhi();
  os << " Tau " << cand.m_tauVeto;
  if (cand.m_jetEtCalibrationLut == 0) {
    os << " using default lut!" << std::endl;
  } else {
    os << " rank " << cand.rank();
    os << " lut address " << cand.m_jetEtCalibrationLut << std::endl;
  }

  return os;
}	

void L1GctJet::setupJet(uint16_t rawsum, unsigned eta, unsigned phi, bool tauVeto)
{
  L1CaloRegionDetId temp(eta, phi);
  m_rawsum = rawsum;
  m_id = temp;
  m_tauVeto = tauVeto;
}

/// Methods to return the jet rank
uint16_t L1GctJet::rank()      const
{
  uint16_t result;
  // If no lut setup, just return the MSB of the rawsum as the rank
  if (m_jetEtCalibrationLut==0) {
    result = std::min(63, m_rawsum >> (RAWSUM_BITWIDTH - 6));
  } else {
    result = m_jetEtCalibrationLut->convertToSixBitRank(m_rawsum, m_id.ieta());
  }
  return result;
}

uint16_t L1GctJet::rankForHt() const
{
  uint16_t result;
  // If no lut setup, just return the MSB of the rawsum as the rank
  if (m_jetEtCalibrationLut==0) {
    result = std::min(1023, m_rawsum >> (RAWSUM_BITWIDTH - 10));
  } else {
    result = m_jetEtCalibrationLut->convertToTenBitRank(m_rawsum, m_id.ieta());
  }
  return result;
}

/// convert to central jet digi
L1GctJetCand L1GctJet::makeJetCand() {
  return L1GctJetCand(this->rank(), this->hwEta(), this->hwPhi(), this->isTauJet(), this->isForwardJet());
}

/// eta value as encoded in hardware at the GCT output
unsigned L1GctJet::hwEta() const
{
  // Force into 4 bits.
  // Count eta bins separately for central and forward jets. Set MSB to indicate the Wheel
  return (((m_id.rctEta() % 7) & 0x7) | (m_id.ieta()<11 ? 0x10 : 0));
}

/// phi value as encoded in hardware at the GCT output
unsigned L1GctJet::hwPhi() const
{
  // Force into 5 bits.
  return m_id.iphi() & 0x1f;
}
