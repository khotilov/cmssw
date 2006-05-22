#include "L1Trigger/GlobalCaloTrigger/interface/L1GctGlobalEnergyAlgos.h"

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelEnergyFpga.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelJetFpga.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetFinalStage.h"

L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() :
  m_jcValPlusWheel(12),
  m_jcVlMinusWheel(12),
  m_jcBoundaryJets(12),
  m_outputJetCounts(12)
{
}

L1GctGlobalEnergyAlgos::~L1GctGlobalEnergyAlgos()
{
}

std::ostream& operator << (std::ostream& os, const L1GctGlobalEnergyAlgos& fpga)
{

  os << "Output Emiss " << fpga.m_outputEtMiss << std::endl;
  os << "Output Emiss Phi " << fpga.m_outputEtMissPhi << std::endl;
  os << "Output EtSum " << fpga.m_outputEtSum << std::endl;
  os << "Output EtHad " << fpga.m_outputEtHad << std::endl;
  os << "Output Jet counts " << std::endl;
  for(unsigned i=0; i < fpga.m_outputJetCounts.size(); i++)
    {
      os << fpga.m_outputJetCounts[i];
    } 

  return os;
}

// clear internal data
void L1GctGlobalEnergyAlgos::reset()
{
  m_exValPlusWheel.reset();
  m_exVlMinusWheel.reset();
  m_eyValPlusWheel.reset();
  m_eyVlMinusWheel.reset();
  m_etValPlusWheel.reset();
  m_etVlMinusWheel.reset();
  m_htValPlusWheel.reset();
  m_htVlMinusWheel.reset();
  m_htBoundaryJets.reset();
  for (int i=0; i<12; i++) {
    m_jcValPlusWheel[i].reset();
    m_jcVlMinusWheel[i].reset();
    m_jcBoundaryJets[i].reset();
  }
  //
  m_outputEtMiss.reset();
  m_outputEtMissPhi.reset();
  m_outputEtSum.reset();
  m_outputEtHad.reset();
  for (int i=0; i<12; i++) {
    m_outputJetCounts[i].reset();
  }
}

void L1GctGlobalEnergyAlgos::fetchInput() {
  // input from WheelEnergyFpgas
  m_exValPlusWheel = m_plusWheelFpga->outputEx();
  m_eyValPlusWheel = m_plusWheelFpga->outputEy();
  m_etValPlusWheel = m_plusWheelFpga->outputEt();
  m_htValPlusWheel = m_plusWheelJetFpga->outputHt();
  
  m_exVlMinusWheel = m_minusWheelFpga->outputEx();
  m_eyVlMinusWheel = m_minusWheelFpga->outputEy();
  m_etVlMinusWheel = m_minusWheelFpga->outputEt();
  m_htVlMinusWheel = m_minusWheelJetFpga->outputHt();

  m_htBoundaryJets = m_jetFinalStage->outputHt();
  //
  for (unsigned i=0; i<12; i++) {
    m_jcValPlusWheel[i] = m_plusWheelJetFpga->outputJc(i);
    m_jcVlMinusWheel[i] = m_minusWheelJetFpga->outputJc(i);
    m_jcBoundaryJets[i] = m_jetFinalStage->outputJc(i);
  }
}


// process the event
void L1GctGlobalEnergyAlgos::process()
{
  L1GctEtComponent ExSum, EySum;
  L1GctGlobalEnergyAlgos::etmiss_vec EtMissing;

  //
  //-----------------------------------------------------------------------------
  // Form the Ex and Ey sums
  ExSum = m_exValPlusWheel + m_exVlMinusWheel;
  EySum = m_eyValPlusWheel + m_eyVlMinusWheel;
  // Execute the missing Et algorithm
  EtMissing = calculate_etmiss_vec(ExSum, EySum);
  //
  m_outputEtMiss    = EtMissing.mag;
  m_outputEtMissPhi = EtMissing.phi;

  //
  //-----------------------------------------------------------------------------
  // Form the Et and Ht sums
  m_outputEtSum = m_etValPlusWheel + m_etVlMinusWheel;
  m_outputEtHad = m_htValPlusWheel + m_htVlMinusWheel + m_htBoundaryJets;

  //
  //-----------------------------------------------------------------------------
  // Add the jet counts.
  for (int i=0; i<12; i++) {
    m_outputJetCounts[i] =
      L1GctJcFinalType(m_jcValPlusWheel[i]) +
      L1GctJcFinalType(m_jcVlMinusWheel[i]) +
      L1GctJcFinalType(m_jcBoundaryJets[i]);
  }
}

//----------------------------------------------------------------------------------------------
// set input data sources
//
void L1GctGlobalEnergyAlgos::setPlusWheelEnergyFpga (L1GctWheelEnergyFpga* fpga)
{
  m_plusWheelFpga = fpga;
}

void L1GctGlobalEnergyAlgos::setMinusWheelEnergyFpga(L1GctWheelEnergyFpga* fpga)
{
  m_minusWheelFpga = fpga;
}

void L1GctGlobalEnergyAlgos::setPlusWheelJetFpga (L1GctWheelJetFpga* fpga)
{
  m_plusWheelJetFpga = fpga;
}

void L1GctGlobalEnergyAlgos::setMinusWheelJetFpga(L1GctWheelJetFpga* fpga)
{
  m_minusWheelJetFpga = fpga;
}

void L1GctGlobalEnergyAlgos::setJetFinalStage(L1GctJetFinalStage* jfs)
{
    m_jetFinalStage = jfs;
}

//----------------------------------------------------------------------------------------------
// set input data per wheel: x component of missing Et
//
void L1GctGlobalEnergyAlgos::setInputWheelEx(unsigned wheel, int energy, bool overflow)
{
  if (wheel==0) {
    m_exValPlusWheel.setValue(energy);
    m_exValPlusWheel.setOverFlow(overflow);
  } else if (wheel==1) {
    m_exVlMinusWheel.setValue(energy);
    m_exVlMinusWheel.setOverFlow(overflow);
  }
}

//----------------------------------------------------------------------------------------------
// set input data per wheel: y component of missing Et
//
void L1GctGlobalEnergyAlgos::setInputWheelEy(unsigned wheel, int energy, bool overflow)
{
  if (wheel==0) {
    m_eyValPlusWheel.setValue(energy);
    m_eyValPlusWheel.setOverFlow(overflow);
  } else if (wheel==1) {
    m_eyVlMinusWheel.setValue(energy);
    m_eyVlMinusWheel.setOverFlow(overflow);
  }
}

//----------------------------------------------------------------------------------------------
// set input data per wheel: scalar sum of Et
//
void L1GctGlobalEnergyAlgos::setInputWheelEt(unsigned wheel, unsigned energy, bool overflow)
{
  if (wheel==0) {
    m_etValPlusWheel.setValue(energy);
    m_etValPlusWheel.setOverFlow(overflow);
  } else if (wheel==1) {
    m_etVlMinusWheel.setValue(energy);
    m_etVlMinusWheel.setOverFlow(overflow);
  }
}

//----------------------------------------------------------------------------------------------
// set input data per wheel: sum of transverse energy in jets (Ht)
//
void L1GctGlobalEnergyAlgos::setInputWheelHt(unsigned wheel, unsigned energy, bool overflow)
{
  if (wheel==0) {
    m_htValPlusWheel.setValue(energy);
    m_htValPlusWheel.setOverFlow(overflow);
  } else if (wheel==1) {
    m_htVlMinusWheel.setValue(energy);
    m_htVlMinusWheel.setOverFlow(overflow);
  }
}


//----------------------------------------------------------------------------------------------
// An extra contribution to Ht from jets at
// the boundary between wheels
//
void L1GctGlobalEnergyAlgos::setInputBoundaryHt(unsigned energy, bool overflow)
{
  m_htBoundaryJets.setValue(energy);
  m_htBoundaryJets.setOverFlow(overflow);
}


//----------------------------------------------------------------------------------------------
// Set the jet count input values
//
void L1GctGlobalEnergyAlgos::setInputWheelJc(unsigned wheel, unsigned jcnum, unsigned count)
{
  if (wheel==0) {
    m_jcValPlusWheel[jcnum].setValue(count);
  } else if (wheel==1) {
    m_jcVlMinusWheel[jcnum].setValue(count);
  }

}


//----------------------------------------------------------------------------------------------
// Extra contributions to jet counts from jets at
// the boundary between wheels
//
void L1GctGlobalEnergyAlgos::setInputBoundaryJc(unsigned jcnum, unsigned count)
{
  m_jcBoundaryJets[jcnum].setValue(count);

}

//----------------------------------------------------------------------------------------------
//
// PRIVATE MEMBER FUNCTION
//
// Here's the Etmiss calculation
//
//-----------------------------------------------------------------------------------
L1GctGlobalEnergyAlgos::etmiss_vec
L1GctGlobalEnergyAlgos::calculate_etmiss_vec (L1GctEtComponent ex, L1GctEtComponent ey)
{
  //---------------------------------------------------------------------------------
  //
  // Calculates magnitude and direction of missing Et, given measured Ex and Ey.
  //
  // The algorithm used is suitable for implementation in hardware, using integer
  // multiplication, addition and comparison and bit shifting operations.
  //
  // Proceed in two stages. The first stage gives a result that lies between
  // 92% and 100% of the true Et, with the direction measured in 45 degree bins.
  // The final precision depends on the number of factors used in corrFact.
  // The present version with eleven factors gives a precision of 1% on Et, and
  // finds the direction to the nearest 5 degrees.
  //
  //---------------------------------------------------------------------------------
  etmiss_vec result;

  unsigned eneCoarse, phiCoarse;
  unsigned eneCorect, phiCorect;

  const unsigned root2fact = 181;
  const unsigned corrFact[11] = {24, 39, 51, 60, 69, 77, 83, 89, 95, 101, 106};
  const unsigned corrDphi[11] = { 0,  1,  2,  2,  3,  3,  3,  3,  4,   4,   4};

  std::vector<bool> s(3);
  unsigned Mx, My, Mw;

  unsigned Dx, Dy;
  unsigned eFact;

  unsigned b,phibin;
  bool midphi;

  // Here's the coarse calculation, with just one multiply operation
  //
  My = static_cast<unsigned>(abs(ey.value()));
  Mx = static_cast<unsigned>(abs(ex.value()));
  Mw = ((Mx+My)*root2fact)>>8;

  s[0] = (ey.value()<0);
  s[1] = (ex.value()<0);
  s[2] = (My>Mx);

  phibin = 0; b = 0;
  for (int i=0; i<3; i++) {
    if (s[i]) { b=1-b;} phibin = 2*phibin + b;
  }

  eneCoarse = std::max(std::max(Mx, My), Mw);
  phiCoarse = phibin*9;

  // For the fine calculation we multiply both input components
  // by all the factors in the corrFact list in order to find
  // the required corrections to the energy and angle
  //
  for (eFact=0; eFact<10; eFact++) {
    Dx = (Mx*corrFact[eFact])>>8;
    Dy = (My*corrFact[eFact])>>8;
    if         ((Dx>My) || (Dy>Mx))         {midphi=false; break;}
    if ((Mx+Dx)>(My-Dy) && (My+Dy)>(Mx-Dx)) {midphi=true;  break;}
  }
  eneCorect = (eneCoarse*(128+eFact))>>7;
  if (midphi ^ (b==1)) {
    phiCorect = phiCoarse + 8 - corrDphi[eFact];
  } else {
    phiCorect = phiCoarse + corrDphi[eFact];
  }

  // Store the result of the calculation
  //
  result.mag.setValue(eneCorect);
  result.phi.setValue(phiCorect);

  result.mag.setOverFlow( result.mag.overFlow() || ex.overFlow() || ey.overFlow() );

  return result;
}

