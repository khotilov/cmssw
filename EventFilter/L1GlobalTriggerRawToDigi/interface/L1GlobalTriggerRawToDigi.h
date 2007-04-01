#ifndef EventFilter_L1GlobalTriggerRawToDigi_L1GlobalTriggerRawToDigi_h
#define EventFilter_L1GlobalTriggerRawToDigi_L1GlobalTriggerRawToDigi_h

/**
 * \class L1GlobalTriggerRawToDigi
 * 
 * 
 * Description: unpack raw data into digitized data.  
 *
 * Implementation:
 *    <TODO: enter implementation details>
 *   
 * \author: Vasile Mihai Ghete - HEPHY Vienna -  GT 
 * \author: Ivan Mikulec       - HEPHY Vienna - GMT
 * 
 * $Date:$
 * $Revision:$
 *
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// forward declarations
class L1GtfeWord;
class L1GtFdlWord;
class L1GtPsbWord;

class L1MuGMTReadoutCollection;

// class declaration
class L1GlobalTriggerRawToDigi : public edm::EDProducer
{

public:

    /// constructor(s)
    explicit L1GlobalTriggerRawToDigi(const edm::ParameterSet&);

    /// destructor
    virtual ~L1GlobalTriggerRawToDigi();

private:

    virtual void beginJob(const edm::EventSetup&);

    virtual void produce(edm::Event&, const edm::EventSetup&);

    /// block unpackers

    /// unpack header
    void unpackHeader(const unsigned char*);

    /// unpack the GTFE block
    /// gives the number of bunch crosses in the event, as well as the active boards
    /// records for inactive boards are not written in the GT DAQ record
    void unpackGTFE(const edm::EventSetup&, const unsigned char*, L1GtfeWord*);

    /// unpack FDL blocks for various bunch crosses
    void unpackFDL(const edm::EventSetup&, const unsigned char*, L1GtFdlWord&);

    /// unpack PSB blocks
    /// unpacking is done in PSB class format
    /// methods are given later to translate from the PSB format
    /// to the physical input of the PSB
    void unpackPSB(const edm::EventSetup&, const unsigned char*, L1GtPsbWord&);

    /// unpack the GMT record
    void unpackGMT(const unsigned char*, std::auto_ptr<L1MuGMTReadoutCollection>&);

    /// unpack trailer word
    void unpackTrailer(const unsigned char*);

    ///
    virtual void endJob();

private:

    L1GtfeWord* m_gtfeWord;
    L1GtPsbWord* m_gtPsbWord;
    L1GtFdlWord* m_gtFdlWord;

    /// total Bx's in the event, obtained from GTFE block    
    int m_totalBxInEvent;

};

#endif // EventFilter_L1GlobalTriggerRawToDigi_L1GlobalTriggerRawToDigi_h
