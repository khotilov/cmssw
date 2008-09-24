#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetFinderBase.h"

#include "CondFormats/L1TObjects/interface/L1GctJetFinderParams.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetEtCalibrationLut.h"

#include "FWCore/Utilities/interface/Exception.h"  

//DEFINE STATICS
const unsigned int L1GctJetFinderBase::MAX_JETS_OUT = 6;
const unsigned int L1GctJetFinderBase::COL_OFFSET = ((L1CaloRegionDetId::N_ETA)/2)+1;
const unsigned int L1GctJetFinderBase::N_JF_PER_WHEEL = ((L1CaloRegionDetId::N_PHI)/2);

const unsigned int L1GctJetFinderBase::MAX_REGIONS_IN = L1GctJetFinderBase::COL_OFFSET*L1GctJetFinderBase::N_COLS;
const unsigned int L1GctJetFinderBase::N_COLS = 2;
const unsigned int L1GctJetFinderBase::CENTRAL_COL0 = 0;


L1GctJetFinderBase::L1GctJetFinderBase(int id):
  L1GctProcessor(),
  m_id(id),
  m_neighbourJetFinders(2),
  m_gotNeighbourPointers(false),
  m_gotJetFinderParams(false),
  m_CenJetSeed(0), m_FwdJetSeed(0), m_TauJetSeed(0), m_EtaBoundry(0),
  m_jetEtCalLut(0),
  m_inputRegions(MAX_REGIONS_IN),
  m_sentProtoJets(MAX_JETS_OUT), m_rcvdProtoJets(MAX_JETS_OUT), m_keptProtoJets(MAX_JETS_OUT),
  m_outputJets(MAX_JETS_OUT), m_sortedJets(MAX_JETS_OUT),
  m_outputEtStrip0(0), m_outputEtStrip1(0),
  m_outputHtStrip0(0), m_outputHtStrip1(0),
  m_outputHfSums(),
  m_outputJetsPipe(MAX_JETS_OUT)
{
  // Call reset to initialise vectors for input and output
  this->reset();
  //Check jetfinder setup
  if(m_id < 0 || m_id >= static_cast<int>(L1CaloRegionDetId::N_PHI))
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctJetFinderBase::L1GctJetFinderBase() : Jet Finder ID " << m_id << " has been incorrectly constructed!\n"
    << "ID number should be between the range of 0 to " << L1CaloRegionDetId::N_PHI-1 << "\n";
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

/// Set pointer to parameters - needed to complete the setup
void L1GctJetFinderBase::setJetFinderParams(const L1GctJetFinderParams* jfpars)
{
  m_CenJetSeed = jfpars->CENTRAL_JET_SEED;
  m_FwdJetSeed = jfpars->FORWARD_JET_SEED;
  m_TauJetSeed = jfpars->TAU_JET_SEED;
  m_EtaBoundry = jfpars->CENTRAL_FORWARD_ETA_BOUNDARY;
  m_gotJetFinderParams = true;
}

/// Set pointer to calibration Lut - needed to complete the setup
void L1GctJetFinderBase::setJetEtCalibrationLut(const L1GctJetEtCalibrationLut* lut)
{
  m_jetEtCalLut = lut;
}

std::ostream& operator << (std::ostream& os, const L1GctJetFinderBase& algo)
{
  using std::endl;
  os << "ID = " << algo.m_id << endl;
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
  os << "Output Ht strip 0 " << algo.m_outputHtStrip0 << endl;
  os << "Output Ht strip 1 " << algo.m_outputHtStrip1 << endl;
  os << endl;

  return os;
}


void L1GctJetFinderBase::resetProcessor()
{
  m_inputRegions.clear();
  m_inputRegions.resize(this->maxRegionsIn());
  m_outputJets.clear();
  m_outputJets.resize(MAX_JETS_OUT);
  m_sortedJets.clear();
  m_sortedJets.resize(MAX_JETS_OUT);
  
  m_sentProtoJets.clear();
  m_sentProtoJets.resize(MAX_JETS_OUT);
  m_rcvdProtoJets.clear();
  m_rcvdProtoJets.resize(MAX_JETS_OUT);
  m_keptProtoJets.clear();
  m_keptProtoJets.resize(MAX_JETS_OUT);

  m_outputEtStrip0 = 0;
  m_outputEtStrip1 = 0;
  m_outputHtStrip0 = 0;
  m_outputHtStrip1 = 0;

  m_outputHfSums.reset();
}

void L1GctJetFinderBase::resetPipelines()
{
  m_outputJetsPipe.reset(numOfBx());
}

/// Initialise inputs with null objects for the correct bunch crossing
/// If no other input candidates "arrive", we have the correct
/// bunch crossing to propagate through the processing.
void L1GctJetFinderBase::setupObjects()
{
  /// Create a null input region with the right bunch crossing, 
  /// and fill the input candidates with copies of this.
  L1GctRegion tempRgn;
  tempRgn.setBx(bxAbs());
  m_inputRegions.assign(this->maxRegionsIn(), tempRgn);

  /// The same for the lists of pre-clustered jets
  /// passed between neighbour jetFinders
  m_sentProtoJets.assign(MAX_JETS_OUT, tempRgn);
  m_rcvdProtoJets.assign(MAX_JETS_OUT, tempRgn);
  m_keptProtoJets.assign(MAX_JETS_OUT, tempRgn);

  /// The same for the lists of output jets
  L1GctJet tempJet;
  tempJet.setBx(bxAbs());
  m_outputJets.assign(MAX_JETS_OUT, tempJet);
}

// This is how the regions from the RCT get into the GCT for processing 
void L1GctJetFinderBase::setInputRegion(const L1CaloRegion& region)
{
  static const unsigned NPHI = L1CaloRegionDetId::N_PHI;
  unsigned crate = region.rctCrate();
  // Find the column for this region in a global (eta,phi) array
  // Note the column numbers here are not the same as region->gctPhi()
  // because the RCT crates are not numbered from phi=0.
  unsigned colAbsolute = crate*2 + region.rctPhi();
  unsigned colRelative = ((colAbsolute+NPHI) - m_minColThisJf) % NPHI;
  if (colRelative < this->nCols()) {
    // We are in the right range in phi
    // Now check we are in the right wheel (positive or negative eta)
    if ( (crate/N_JF_PER_WHEEL) == (m_id/N_JF_PER_WHEEL) ) {
      unsigned i = colRelative*COL_OFFSET + region.rctEta() + 1;
      L1GctRegion temp(region);
      m_inputRegions.at(i) = temp;
    } else {
      // Accept neighbouring regions from the other wheel
      if (region.rctEta() == 0) {
	unsigned i = colRelative*COL_OFFSET;
        L1GctRegion temp(region);
	m_inputRegions.at(i) = temp;
      }
    }
  }
}

// PROTECTED METHODS BELOW
/// fetch the protoJets from neighbour jetFinder
void L1GctJetFinderBase::fetchProtoJetsFromNeighbour(const fetchType ft)
{
  switch (ft) {
  case TOP : 
    m_rcvdProtoJets = m_neighbourJetFinders.at(0)->getSentProtoJets(); break;
  case BOT :
    m_rcvdProtoJets = m_neighbourJetFinders.at(1)->getSentProtoJets(); break;
  case TOPBOT :
    // Copy half the jets from each neighbour
    static const unsigned int MAX_TOPBOT_JETS = MAX_JETS_OUT/2;
    unsigned j=0;
    RegionsVector temp;
    temp = m_neighbourJetFinders.at(0)->getSentProtoJets();
    for ( ; j<MAX_TOPBOT_JETS; ++j) {
      m_rcvdProtoJets.at(j) = temp.at(j);
    } 
    temp = m_neighbourJetFinders.at(1)->getSentProtoJets();
    for ( ; j<MAX_JETS_OUT; ++j) {
      m_rcvdProtoJets.at(j) = temp.at(j);
    }     
    break;
  }
}


/// Sort the found jets. All jetFinders should call this in process().
void L1GctJetFinderBase::sortJets()
{
  //transform the jets to the final GCT output format
  for (unsigned j=0; j<MAX_JETS_OUT; ++j) {
    m_sortedJets.at(j) = m_outputJets.at(j).jetCand(m_jetEtCalLut);
  }
  //presort the jets into descending order of energy
  sort(m_sortedJets.begin(), m_sortedJets.end(), rankGreaterThan());
  //store jets in "pipeline memory" for checking
  m_outputJetsPipe.store(m_outputJets, bxRel());
}
   
/// Fill the Et strip sums and Ht sum. All jetFinders should call this in process().
void L1GctJetFinderBase::doEnergySums()
{
  //calculate the raw Et strip sums
  m_outputEtStrip0 = calcEtStrip(0);
  m_outputEtStrip1 = calcEtStrip(1);

  //calculate the Ht
  m_outputHtStrip0 = calcHtStrip(0);
  m_outputHtStrip1 = calcHtStrip(1);

  //calculate the Hf tower Et sums and tower-over-threshold counts
  m_outputHfSums = calcHfSums();
    
  return;
}


// Calculates total (raw) energy in a phi strip
L1GctJetFinderBase::etTotalType L1GctJetFinderBase::calcEtStrip(const UShort strip) const
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
  etTotalType temp(et);
  temp.setOverFlow(temp.overFlow() || of);
  return temp;
}

// Calculates total calibrated energy in jets (Ht) sum
L1GctJetFinderBase::etTotalType L1GctJetFinderBase::calcHtStrip(const UShort strip) const
{    
  unsigned ht = 0;
  bool of = false;
  for(UShort i=0; i < MAX_JETS_OUT; ++i)
  {
    // Only sum Ht for valid jets
    if (!m_outputJets.at(i).isNullJet()) {
      if (m_outputJets.at(i).rctPhi() == strip) {
	ht += m_outputJets.at(i).calibratedEt(m_jetEtCalLut);
	of |= m_outputJets.at(i).overFlow();
      }
    }
  }
  etHadType temp(ht);
  temp.setOverFlow(temp.overFlow() || of);
  return temp;
}

// Calculates Hf inner rings Et sum, and counts number of "fineGrain" bits set
L1GctJetFinderBase::hfTowerSumsType L1GctJetFinderBase::calcHfSums() const
{
  static const UShort NUMBER_OF_FRWRD_RINGS = 4;
  static const UShort NUMBER_OF_INNER_RINGS = 2;
  std::vector<unsigned> et(NUMBER_OF_INNER_RINGS, 0);
  std::vector<bool>     of(NUMBER_OF_INNER_RINGS, false);
  std::vector<unsigned> nt(NUMBER_OF_INNER_RINGS, 0);

  UShort offset = COL_OFFSET*(centralCol0() + 1);
  for (UShort i=0; i < NUMBER_OF_FRWRD_RINGS; ++i) {
    offset--;

    // Sum HF Et and count jets above threshold over "inner rings"
    if (i<NUMBER_OF_INNER_RINGS) {
      et.at(i) += m_inputRegions.at(offset).et();
      of.at(i) = of.at(i) || m_inputRegions.at(offset).overFlow();

      et.at(i) += m_inputRegions.at(offset+COL_OFFSET).et();
      of.at(i) = of.at(i) || m_inputRegions.at(offset+COL_OFFSET).overFlow();

      if (m_inputRegions.at(offset).fineGrain()) nt.at(i)++;
      if (m_inputRegions.at(offset+COL_OFFSET).fineGrain()) nt.at(i)++;
    }
  }
  hfTowerSumsType temp(et.at(0), et.at(1), nt.at(0), nt.at(1));
  temp.etSum0.setOverFlow(temp.etSum0.overFlow() || of.at(0));
  temp.etSum1.setOverFlow(temp.etSum1.overFlow() || of.at(1));
  return temp;
}
