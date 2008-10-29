/*
  Test distributor class
  Validates test patterns at the distributor
  level (i.e. for lumi sections)
  Copyright John Jones 2007
*/

// Prototype class definition
#include "TestDistributor.hh"

// Various generic and type-specific exception implementations
#include "ICException.hh"
#include "HardwareAccessException.hh"
#include "MemoryAllocationException.hh"
//#include "ArgumentOutOfRangeException.hh"

// HCAL HLX namespace
namespace HCAL_HLX
{

  // Define namespaces in use
  using namespace ICCoreUtils;
  using namespace std;
  
  // Default constructor
  TestDistributor::TestDistributor() {
    mNumGoodSections = 0;
    mNumBadSections = 0;
  }

  // Destructor deletes the hardware interface
  TestDistributor::~TestDistributor() {}

  bool TestDistributor::ProcessSection(const LUMI_SECTION & lumiSection) {
    try {

      for ( u16 k = 0 ; k != lumiSection.hdr.numHLXs ; k++ ) {
	u32 i, value;
	bool bExit = false;
	for ( i = 0 ; i != lumiSection.hdr.numBunches ; i++ ) {
	  for ( u32 j = 0 ; j != 6 ; j++ ) {
	    if ( lumiSection.occupancy[k].hdr.numNibbles == 0 ) {
	      cout << "Occupancy " << dec << k << " lumi section contains no nibbles!!!" << endl;
	      break;
	    }
	    value = i * lumiSection.hdr.numOrbits / lumiSection.occupancy[k].hdr.numNibbles;
	    if ( value > 0xFFFF ) value = 0xFFFF;
	    value *= lumiSection.occupancy[k].hdr.numNibbles;
	    if ( lumiSection.occupancy[k].data[j][i] != value ) {
	      cout << "occ" << static_cast<u16>(k) << "\t" << j << "\t" << i << "\t" << lumiSection.occupancy[k].data[j][i] << "\t" << value << endl;
	      bExit = true;
	      break;
	    }
	  }
	  if ( bExit ) break;
	}
	if ( i != lumiSection.hdr.numBunches ) {
	  mNumBadSections++;
	  return true;
	}
	
	for ( i = 0 ; i != lumiSection.hdr.numBunches ; i++ ) {
	  if ( lumiSection.lhc[k].hdr.numNibbles == 0 ) {
	    cout << "LHC " << dec << k << " lumi section contains no nibbles!" << endl;
	    break;
	  }
	  value = i * lumiSection.hdr.numOrbits / lumiSection.lhc[k].hdr.numNibbles;
	  if ( value > 0xFFFF ) value = 0xFFFF;
	  value *= lumiSection.lhc[k].hdr.numNibbles;
	  if ( lumiSection.lhc[k].data[i] != value ) {
	    cout << "lhc" << static_cast<u16>(k) << "\t" << i << "\t" << lumiSection.lhc[k].data[i] << "\t" << value << endl;
	    //break;
	  }
	}	
	if ( i != lumiSection.hdr.numBunches ) {
	  mNumBadSections++;
	  return true;
	}

	for ( i = 0 ; i != lumiSection.hdr.numBunches ; i++ ) {
	  if ( lumiSection.etSum[k].hdr.numNibbles == 0 ) {
	    cout << "etSum " << dec << k << " lumi section contains no nibbles!!!" << endl;
	    break;
	  }
	  value = i * lumiSection.hdr.numOrbits;
	  if ( lumiSection.etSum[k].data[i] != value ) {
	    cout << "et" << static_cast<u16>(k) << "\t" << i << "\t" << lumiSection.etSum[k].data[i] << "\t" << value << endl;
	    break;
	  }
	}
	if ( i != lumiSection.hdr.numBunches ) {
	  mNumBadSections++;
	} else {
	  mNumGoodSections++;
	}
      }
      return true;
    } catch (ICException & aExc) {
      cerr << aExc.what() << endl;
      return false;
      //RETHROW(aExc);
    }
  }


} // namespace HCAL_HLX
