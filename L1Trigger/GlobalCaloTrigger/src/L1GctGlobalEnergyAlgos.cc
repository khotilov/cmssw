#include "L1Trigger/GlobalCaloTrigger/interface/L1GctGlobalEnergyAlgos.h"

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelEnergyFpga.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelJetFpga.h"

#include "FWCore/Utilities/interface/Exception.h"

using std::ostream;
using std::endl;
using std::vector;
using std::max;

L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos(L1GctWheelEnergyFpga* plusWheelFpga,
					       L1GctWheelEnergyFpga* minusWheelFpga,
					       L1GctWheelJetFpga* plusWheelJetFpga,
					       L1GctWheelJetFpga* minusWheelJetFpga) :
  m_plusWheelFpga(plusWheelFpga),
  m_minusWheelFpga(minusWheelFpga),
  m_plusWheelJetFpga(plusWheelJetFpga),
  m_minusWheelJetFpga(minusWheelJetFpga),
  m_jcValPlusWheel(12),
  m_jcVlMinusWheel(12),
  m_outputJetCounts(12)
{
    if(m_plusWheelFpga == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() has been incorrectly constructed!\n"
      << "Plus Wheel Fpga pointer has not been set!\n";
    }
    if(m_minusWheelFpga == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() has been incorrectly constructed!\n"
      << "Minus Wheel Fpga pointer has not been set!\n";
    }
    if(m_plusWheelJetFpga == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() has been incorrectly constructed!\n"
      << "Plus Wheel Jet Fpga pointer has not been set!\n";
    }
    if(m_minusWheelJetFpga == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctGlobalEnergyAlgos::L1GctGlobalEnergyAlgos() has been incorrectly constructed!\n"
      << "Minus Wheel Jet Fpga pointer has not been set!\n";
    }
}

L1GctGlobalEnergyAlgos::~L1GctGlobalEnergyAlgos()
{
}

ostream& operator << (ostream& os, const L1GctGlobalEnergyAlgos& fpga)
{
  os << "=== Global Energy Algos ===" << endl;
  os << "Inputs from Plus wheel:" << endl;
  os << "  Ex " << fpga.m_exValPlusWheel << "  Ey " << fpga.m_eyValPlusWheel;
  os << "  Et " << fpga.m_etValPlusWheel << "  Ht " << fpga.m_htValPlusWheel << endl; 
  os << "Inputs from Minus wheel:" << endl;
  os << "  Ex " << fpga.m_exVlMinusWheel << "  Ey " << fpga.m_eyVlMinusWheel;
  os << "  Et " << fpga.m_etVlMinusWheel << "  Ht " << fpga.m_htVlMinusWheel << endl; 
  os << "Input Jet counts " << endl;
  for(unsigned i=0; i < fpga.m_jcValPlusWheel.size(); i++)
    {
      os << "  Plus wheel  " << i << ": " << fpga.m_jcValPlusWheel[i];
      os << "  Minus wheel " << i << ": " << fpga.m_jcVlMinusWheel[i];
    } 
  os << endl;
  // No endlines here since the overloaded << for
  // the UnsignedInt data type includes it.
  os << "Output Etmiss " << fpga.m_outputEtMiss;
  os << "Output Etmiss Phi " << fpga.m_outputEtMissPhi;
  os << "Output EtSum " << fpga.m_outputEtSum;
  os << "Output EtHad " << fpga.m_outputEtHad;
  os << "Output Jet counts " << endl;
  for(unsigned i=0; i < fpga.m_outputJetCounts.size(); i++)
    {
      os << i << ": " << fpga.m_outputJetCounts[i];
    } 
  os << endl;

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
  for (int i=0; i<12; i++) {
    m_jcValPlusWheel[i].reset();
    m_jcVlMinusWheel[i].reset();
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
  m_exValPlusWheel = m_plusWheelFpga->getOutputEx();
  m_eyValPlusWheel = m_plusWheelFpga->getOutputEy();
  m_etValPlusWheel = m_plusWheelFpga->getOutputEt();
  m_htValPlusWheel = m_plusWheelJetFpga->getOutputHt();
  
  m_exVlMinusWheel = m_minusWheelFpga->getOutputEx();
  m_eyVlMinusWheel = m_minusWheelFpga->getOutputEy();
  m_etVlMinusWheel = m_minusWheelFpga->getOutputEt();
  m_htVlMinusWheel = m_minusWheelJetFpga->getOutputHt();

  //
  for (unsigned i=0; i<12; i++) {
    m_jcValPlusWheel[i] = m_plusWheelJetFpga->getOutputJc(i);
    m_jcVlMinusWheel[i] = m_minusWheelJetFpga->getOutputJc(i);
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
  m_outputEtHad = m_htValPlusWheel + m_htVlMinusWheel;

  //
  //-----------------------------------------------------------------------------
  // Add the jet counts.
  for (int i=0; i<12; i++) {
    m_outputJetCounts[i] =
      L1GctJcFinalType(m_jcValPlusWheel[i]) +
      L1GctJcFinalType(m_jcVlMinusWheel[i]);
  }
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

  vector<bool> s(3);
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

  eneCoarse = max(max(Mx, My), Mw);
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

