#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelEnergyFpga.h"

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetLeafCard.h"

#include "FWCore/Utilities/interface/Exception.h"

using std::vector;
using std::ostream;
using std::endl;

L1GctWheelEnergyFpga::L1GctWheelEnergyFpga(int id, vector<L1GctJetLeafCard*> leafCards) :
	m_id(id),
        m_inputLeafCards(leafCards),
	m_inputEx(3),
	m_inputEy(3),
	m_inputEt(3)
{
  //Check wheelEnergyFpga setup
  if(m_id != 0 && m_id != 1)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctWheelEnergyFpga::L1GctWheelEnergyFpga() : Wheel Energy Fpga ID " << m_id << " has been incorrectly constructed!\n"
    << "ID number should be 0 or 1.\n";
  } 
  
  if(m_inputLeafCards.size() != 3)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctWheelEnergyFpga::L1GctWheelEnergyFpga() : Wheel Energy Fpga ID " << m_id << " has been incorrectly constructed!\n"
    << "This class needs 3 leaf card pointers, yet only " << m_inputLeafCards.size()
    << " leaf card pointers are present.\n";
  }
  
  for(unsigned int i = 0; i < m_inputLeafCards.size(); ++i)
  {
    if(m_inputLeafCards[i] == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctWheelEnergyFpga::L1GctWheelEnergyFpga() : Wheel Energy Fpga ID " << m_id << " has been incorrectly constructed!\n"
      << "Input Leaf card pointer " << i << " has not been set!\n";
    }
  }
}

L1GctWheelEnergyFpga::~L1GctWheelEnergyFpga()
{
}

ostream& operator << (ostream& os, const L1GctWheelEnergyFpga& fpga)
{
  os << "===L1GctWheelEnergyFPGA===" << endl;
  os << "ID : " << fpga.m_id << endl;
  os << "No. of Input Leaf Cards " << fpga.m_inputLeafCards.size() << endl;
  for(unsigned i=0; i < fpga.m_inputLeafCards.size(); i++)
    {
      os << "LeafCard* " << i << " = " << fpga.m_inputLeafCards[i] << endl;
    } 
  os << "Input Ex " << endl;
  for(unsigned i=0; i < fpga.m_inputEx.size(); i++)
    {
      os << fpga.m_inputEx[i] << endl;
    } 
  os << "Input Ey " << endl;
  for(unsigned i=0; i < fpga.m_inputEy.size(); i++)
    {
      os << fpga.m_inputEy[i] << endl;
    } 
  os << "Input Et " << endl;
  for(unsigned i=0; i < fpga.m_inputEt.size(); i++)
    {
      os << fpga.m_inputEt[i] << endl;
    } 
  os << "Output Ex " << fpga.m_outputEx << endl;
  os << "Output Ey " << fpga.m_outputEy << endl;
  os << "Output Et " << fpga.m_outputEt << endl;
  os << endl;
  return os;
}

void L1GctWheelEnergyFpga::reset()
{
  for (int i=0; i<3; i++) {
    m_inputEx[i].reset();
    m_inputEy[i].reset();
    m_inputEt[i].reset();
  }
  m_outputEx.reset();
  m_outputEy.reset();
  m_outputEt.reset();
}

void L1GctWheelEnergyFpga::fetchInput()
{
  // Fetch the output values from each of our input leaf cards.
  for (int i=0; i<3; i++) {
    m_inputEx[i] = m_inputLeafCards[i]->getOutputEx();
    m_inputEy[i] = m_inputLeafCards[i]->getOutputEy();
    m_inputEt[i] = m_inputLeafCards[i]->getOutputEt();
  }
}

void L1GctWheelEnergyFpga::process()
{

  m_outputEx = m_inputEx[0] + m_inputEx[1] + m_inputEx[2];
  m_outputEy = m_inputEy[0] + m_inputEy[1] + m_inputEy[2];
  m_outputEt = m_inputEt[0] + m_inputEt[1] + m_inputEt[2];

}


///
/// set input data
void L1GctWheelEnergyFpga::setInputEnergy(int i, int ex, int ey, unsigned et)
{
  // Set the three input values from this Leaf card
  if (i>=0 && i<3) {
    m_inputEx[i].setValue(ex);
    m_inputEy[i].setValue(ey);
    m_inputEt[i].setValue(et);
  }

}

