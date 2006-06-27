#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetFinderBase.h"
 
#include "FWCore/Utilities/interface/Exception.h"  

#include <iostream>
using namespace std;

//DEFINE STATICS
const unsigned int L1GctJetFinderBase::MAX_JETS_OUT = 6;
const unsigned int L1GctJetFinderBase::MAX_SOURCE_CARDS = 9;
const unsigned int L1GctJetFinderBase::COL_OFFSET = ((L1GctMap::N_RGN_ETA)/2)+1;
const unsigned int L1GctJetFinderBase::N_JF_PER_WHEEL = ((L1GctMap::N_RGN_PHI)/2);

const unsigned int L1GctJetFinderBase::MAX_REGIONS_IN = L1GctJetFinderBase::COL_OFFSET*L1GctJetFinderBase::N_COLS;
const int L1GctJetFinderBase::N_COLS = 2;
const unsigned int L1GctJetFinderBase::CENTRAL_COL0 = 0;


L1GctJetFinderBase::L1GctJetFinderBase(int id, vector<L1GctSourceCard*> sourceCards,
				       L1GctJetEtCalibrationLut* jetEtCalLut):
  m_id(id),
  m_sourceCards(sourceCards),
  m_neighbourJetFinders(2),
  m_gotNeighbourPointers(false),
  m_jetEtCalLut(jetEtCalLut),
  m_inputRegions(MAX_REGIONS_IN),
  m_outputJets(MAX_JETS_OUT)
{
  map = L1GctMap::getMap();
  // Call reset to initialise vectors for input and output
  this->reset();
  //Check jetfinder setup
  if(m_id < 0 || m_id > 17)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctJetFinderBase::L1GctJetFinderBase() : Jet Finder ID " << m_id << " has been incorrectly constructed!\n"
    << "ID number should be between the range of 0 to 17\n";
  } 
  
  if(m_sourceCards.size() != MAX_SOURCE_CARDS)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctJetFinderBase::L1GctJetFinderBase() : Jet Finder ID " << m_id << " has been incorrectly constructed!\n"
    << "This class needs " << MAX_SOURCE_CARDS << " source card pointers, yet only " << m_sourceCards.size()
    << " source card pointers are present.\n";
  }
  
  for(unsigned int i = 0; i < m_sourceCards.size(); ++i)
  {
    if(m_sourceCards.at(i) == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctJetFinderBase::L1GctJetFinderBase() : Jet Finder ID " << m_id << " has been incorrectly constructed!\n"
      << "Source card pointer " << i << " has not been set!\n";
    }
  }
  
  if(m_jetEtCalLut == 0)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctJetFinderBase::L1GctJetFinderBase() : Jet Finder ID " << m_id << " has been incorrectly constructed!\n"
    << "The jet Et calibration LUT pointer has not been set!\n";  
  }
}

L1GctJetFinderBase::~L1GctJetFinderBase()
{
}

/// Set pointers to neighbours
void L1GctJetFinderBase::setNeighbourJetFinders(std::vector<L1GctJetFinderBase*> neighbours)
{
  if (neighbours.size()==2) {
    m_neighbourJetFinders = neighbours;
  } else {
    throw cms::Exception("L1GctSetupError")
      << "L1GctJetFinderBase::setNeighbourJetFinders() : In Jet Finder ID " << m_id 
      << " size of input vector should be 2, but is in fact " << neighbours.size() << "\n";
  }
  if (m_neighbourJetFinders.at(0) == 0) {
    throw cms::Exception("L1GctSetupError")
      << "L1GctJetFinderBase::setNeighbourJetFinders() : In Jet Finder ID " << m_id 
      << " first neighbour pointer is set to zero\n";
  }
  if (m_neighbourJetFinders.at(1) == 0) {
    throw cms::Exception("L1GctSetupError")
      << "L1GctJetFinderBase::setNeighbourJetFinders() : In Jet Finder ID " << m_id 
      << " second neighbour pointer is set to zero\n";
  }
  m_gotNeighbourPointers = true;
}

ostream& operator << (ostream& os, const L1GctJetFinderBase& algo)
{
  os << "ID = " << algo.m_id << endl;
  os << "No of Source cards " << algo.m_sourceCards.size() << endl;
  for (unsigned i=0; i<algo.m_sourceCards.size(); i++) {
    os << "SourceCard* " << i << " = " << algo.m_sourceCards.at(i)<< endl;
    os << "No of regions from this sourceCard " << algo.m_sourceCards.at(i)->getRegions().size() << endl;
  }
  os << "JetEtCalibrationLut* = " <<  algo.m_jetEtCalLut << endl;
  os << "No of input regions " << algo.m_inputRegions.size() << endl;
//   for(unsigned i=0; i < algo.m_inputRegions.size(); ++i)
//     {
//       os << algo.m_inputRegions.at(i); 
//     }
  os << "No of output jets " << algo.m_outputJets.size() << endl;
//   for(unsigned i=0; i < algo.m_outputJets.size(); ++i)
//     {
//       os << algo.m_outputJets.at(i); 
//     }
  os << "Output Et strip 0 " << algo.m_outputEtStrip0 << endl;
  os << "Output Et strip 1 " << algo.m_outputEtStrip1 << endl;
  os << "Output Ht " << algo.m_outputHt << endl;
  os << endl;

  return os;
}


void L1GctJetFinderBase::reset()
{
  m_inputRegions.clear();
  m_inputRegions.resize(this->maxRegionsIn());
  m_outputJets.clear();
  m_outputJets.resize(MAX_JETS_OUT);
  
  // All jets need a pointer to the calibration LUT,
  // so set one up here.
  JetVector::iterator currentJet;  
  for(currentJet = m_outputJets.begin(); currentJet != m_outputJets.end(); ++currentJet)
  {
    currentJet->setLut(m_jetEtCalLut);
  }

  m_outputEtStrip0 = 0;
  m_outputEtStrip1 = 0;
  m_outputHt = 0;
}

void L1GctJetFinderBase::setInputRegion(unsigned i, L1GctRegion region)
{
  if(i >= 0 && i < this->maxRegionsIn())
  {
    m_inputRegions.at(i) = region;
  }
  else
  {
    throw cms::Exception("L1GctInputError")
    << "L1GctJetFinderBase::setInputRegion() : In Jet Finder ID " << m_id << ", inputted region " 
    << i << " is outside input index range of 0 to " << (this->maxRegionsIn()-1) << "\n";
  }
}

// PROTECTED METHODS BELOW
/// Get the input regions for the 2x11 search window plus eta=0 neighbours
void L1GctJetFinderBase::fetchCentreStripsInput() {
  fetchScInput(m_sourceCards.at(0), 2, this->centralCol0());
  fetchScInput(m_sourceCards.at(1), 3, this->centralCol0());
  fetchNeighbourScInput(m_sourceCards.at(2), this->centralCol0());
}

/// Get the input regions for adjacent 2x11 search windows plus eta=0 neighbours
void L1GctJetFinderBase::fetchEdgeStripsInput() {
  fetchScInput(m_sourceCards.at(3), 2, (this->centralCol0()+2));
  fetchScInput(m_sourceCards.at(4), 3, (this->centralCol0()+2));
  fetchNeighbourScInput(m_sourceCards.at(5), (this->centralCol0()+2));
  fetchScInput(m_sourceCards.at(6), 2, (this->centralCol0()-2));
  fetchScInput(m_sourceCards.at(7), 3, (this->centralCol0()-2));
  fetchNeighbourScInput(m_sourceCards.at(8), (this->centralCol0()-2));
}

/// Copy the input regions from one source card into the m_inputRegions vector
void L1GctJetFinderBase::fetchScInput(L1GctSourceCard* sourceCard, unsigned scType, int col0) {
  for (unsigned i=0; i<sourceCard->getRegions().size(); ++i) {
    unsigned localEta = map->rctEta(scType, i);
    unsigned localPhi = map->rctPhi(scType, i);

    int col = col0+localPhi;
    if (col>=0 && col<this->nCols()) {
      unsigned offset = col*COL_OFFSET + localEta + 1;
      m_inputRegions.at(offset) = sourceCard->getRegions().at(i);
    }
  }
}

/// Copy the input regions from one eta=0 neighbour source card
void L1GctJetFinderBase::fetchNeighbourScInput(L1GctSourceCard* sourceCard, int col0) {
  for (unsigned iphi=0; iphi<2; ++iphi) {
    int col = col0+iphi;
    if (col>=0 && col<this->nCols()) {
      unsigned scOutput = map->sourceCardOutput(0, iphi);

      unsigned offset = col*COL_OFFSET;
      m_inputRegions.at(offset) = sourceCard->getRegions().at(scOutput);
    }
  }
}

/// Sort the found jets. All jetFinders should call this in process().
void L1GctJetFinderBase::sortJets()
{
  //presort the jets into decending order of energy
  sort(m_outputJets.begin(), m_outputJets.end(), L1GctJet::rankGreaterThan());
}
   
/// Fill the Et strip sums and Ht sum. All jetFinders should call this in process().
void L1GctJetFinderBase::doEnergySums()
{
  //calculate the raw Et strip sums
  m_outputEtStrip0 = calcEtStrip(0);
  m_outputEtStrip1 = calcEtStrip(1);

  //calculate the Ht
  m_outputHt = calcHt();
    
  return;
}


// Calculates total (raw) energy in a phi strip
L1GctScalarEtVal L1GctJetFinderBase::calcEtStrip(const UShort strip) const
{
  if (strip !=0 && strip != 1) {
    throw cms::Exception("L1GctProcessingError")
      << "L1GctJetFinderBase::calcEtStrip() has been called with strip number "
      << strip << "; should be 0 or 1 \n";
  } 
  // Add the Et values from regions 13 to 23 for strip 0,
  //     the Et values from regions 25 to 35 for strip 1.
  unsigned et = 0;
  bool of = false;
  unsigned offset = COL_OFFSET * (strip+centralCol0());
  for (UShort i=1; i < COL_OFFSET; ++i) {
    offset++;
    et += m_inputRegions.at(offset).et();
    of |= m_inputRegions.at(offset).overFlow();
  }
  L1GctScalarEtVal temp(et);
  temp.setOverFlow(temp.overFlow() || of);
  return temp;
}

// Calculates total calibrated energy in jets (Ht) sum
L1GctScalarEtVal L1GctJetFinderBase::calcHt() const
{    
  unsigned ht = 0;
  for(UShort i=0; i < MAX_JETS_OUT; ++i)
  {
    // Only sum Ht for valid jets
    if (!m_outputJets.at(i).isNullJet()) {
      ht += static_cast<unsigned>(m_outputJets.at(i).rankForHt());
    }
  }
  L1GctScalarEtVal temp(ht);
  return temp;
}
