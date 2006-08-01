#ifndef STREAMER_EVENTSTREAMHTTPREADER_H
#define STREAMER_EVENTSTREAMHTTPREADER_H

#include "IOPool/Streamer/interface/EventBuffer.h"
#include "IOPool/Streamer/interface/Utilities.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ProductRegistry.h"
#include "FWCore/Framework/interface/InputSource.h"

#include <vector>
#include <memory>
#include <string>
#include <fstream>

namespace edmtestp
{
  struct ReadData;

  class EventStreamHttpReader : public edm::InputSource
  {
  public:
    typedef std::vector<char> Buf;

    EventStreamHttpReader(edm::ParameterSet const& pset,
		 edm::InputSourceDescription const& desc);
    virtual ~EventStreamHttpReader();

    virtual std::auto_ptr<edm::EventPrincipal> read();
    virtual std::auto_ptr<edm::SendJobHeader> readHeader();

  private:  
    std::string sourceurl_;
    char eventurl_[256];
    char headerurl_[256];
    Buf buf_;
    int events_read_;      // use this until we inherent from the right input source
    //edm::EventDecoder decoder_;
    int hltBitCount;
    int l1BitCount;
  };

}
#endif

