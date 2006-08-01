/*
  Input source for event consumers that will get events from the
  Storage Manager Event Server. This does uses a HTTP get using the
  cURL library. The Storage Manager Event Server responses with
  a binary octet-stream.  The product registry is also obtained
  through a HTTP get.
      There is currently no test of the product registry against
  the consumer client product registry within the code. It should
  already be done if this was inherenting from the standard
  framework input source. Currently we inherit from InputSource.

  17 Mar 2006 - HWKC - initial version for testing
  30 Mar 2006 - HWKC - first proof of principle version that can
                wait for the the product registry and events from
                the Storage Manager. Only way to stop the cmsRun
                using this input source is to kill the Storage
                Manager or specify a maximum number of events for
                the client to read through a maxEvents parameter.
*/

#include "EventFilter/StorageManager/src/EventStreamHttpReader.h"
#include "IOPool/Streamer/interface/BufferArea.h"
#include "FWCore/Utilities/interface/DebugMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "IOPool/Streamer/interface/ClassFiller.h"
#include "IOPool/Streamer/interface/Utilities.h"
#include "IOPool/Streamer/interface/StreamTranslator.h"
#include "IOPool/Streamer/interface/OtherMessage.h"

#include <algorithm>
#include <iterator>
#include "curl/curl.h"
#include <string>

#include <wait.h>

using namespace std;
using namespace edm;

namespace edmtestp
{  
  struct ReadData
  {
    std::string d_;
  };  

  size_t func(void* buf,size_t size, size_t nmemb, void* userp)
  {
    ReadData* rdata = (ReadData*)userp;
    size_t sz = size * nmemb;
    char* cbuf = (char*)buf;
    rdata->d_.insert(rdata->d_.end(),cbuf,cbuf+sz);
    return sz;
  }

  template <class Han, class Opt, class Par>
  int setopt(Han han,Opt opt,Par par)
  {
    if(curl_easy_setopt(han,opt,par)!=0)
      {
        cerr << "could not setopt " << opt << endl;
        abort();
      }
    return 0;
  }

  // ----------------------------------

  EventStreamHttpReader::EventStreamHttpReader(edm::ParameterSet const& ps,
					       edm::InputSourceDescription const& desc):
    edm::InputSource(ps, desc),
    sourceurl_(ps.getParameter<string>("sourceURL")),
    buf_(1000*1000*7), 
    events_read_(0)
  {
    std::string evturl = sourceurl_ + "/geteventdata";
    int stlen = evturl.length();
    for (int i=0; i<stlen; i++) eventurl_[i]=evturl[i];
    eventurl_[stlen] = '\0';

    std::string header = sourceurl_ + "/getregdata";
    stlen = header.length();
    for (int i=0; i<stlen; i++) headerurl_[i]=header[i];
    headerurl_[stlen] = '\0';

    std::auto_ptr<SendJobHeader> p = readHeader();
    SendDescs & descs = p->descs_;
    edm::mergeWithRegistry(descs, productRegistry());

    // next taken from IOPool/Streamer/EventStreamFileReader
    // jbk - the next line should not be needed
    edm::declareStreamers(descs);
    edm::buildClassCache(descs);
    loadExtraClasses();
  }

  EventStreamHttpReader::~EventStreamHttpReader()
  {
  }

  std::auto_ptr<edm::EventPrincipal> EventStreamHttpReader::read()
  {
    // repeat a http get every 5 seconds until we get an event
    // wait for Storage Manager event server buffer to not be empty
    // only way to stop is specify a maxEvents parameter
    // or kill the STorage Manager so the http get fails.
    // do it like test for the proof of principle test

    // see if already read maxEvents
    if(maxEvents() > 0 && events_read_ >= maxEvents()) 
      return std::auto_ptr<edm::EventPrincipal>();

    ReadData data;
    do {
      CURL* han = curl_easy_init();

      if(han==0)
      {
        cerr << "could not create handle" << endl;
        // this will end cmsRun 
        return std::auto_ptr<edm::EventPrincipal>();
      }

      setopt(han,CURLOPT_URL,eventurl_);
      setopt(han,CURLOPT_WRITEFUNCTION,func);
      setopt(han,CURLOPT_WRITEDATA,&data);

      if(curl_easy_perform(han)!=0)
      {
        cerr << "curl perform failed for event" << endl;
        // this will end cmsRun 
        return std::auto_ptr<edm::EventPrincipal>();
      }
      curl_easy_cleanup(han);
      if(data.d_.length() == 0)
      {
        std::cout << "...waiting for event from Storage Manager..." << std::endl;
        // sleep for 5 seconds
        sleep(5);
      }
    } while (data.d_.length() == 0);

    int len = data.d_.length();
    FDEBUG(9) << "EventStreamHttpReader received len = " << len << std::endl;
    buf_.resize(len);
    for (int i=0; i<len ; i++) buf_[i] = data.d_[i];

    // first check if done message
    // need to use this BUT does EventMsgView EVER look like a MsgCode DONE message!!!
    //edm::MsgCode msgtest(&buf_[0],len);
    //if(msgtest.getCode() == MsgCode::DONE) {
    // OtherMessageView class not working
    OtherMessageView msgView(&buf_[0]);
    //std::cout << "received other message code = " << msgView.code()
    //          << " and size = " << msgView.size()
    //          << " and check against " << Header::DONE << endl;

    if (msgView.code() == Header::DONE) {
      // this will end cmsRun 
      std::cout << "Storage Manager has halted - ending run" << std::endl;
      return std::auto_ptr<edm::EventPrincipal>();
    } else {
      events_read_++;
      //edm::EventMsg msg(&buf_[0],len);
      //return decoder_.decodeEvent(msg,productRegistry());
      EventMsgView eventView(&buf_[0],hltBitCount,l1BitCount);
      return StreamTranslator::deserializeEvent(eventView,productRegistry());
    }
  }

  std::auto_ptr<SendJobHeader> EventStreamHttpReader::readHeader()
  {
    // repeat a http get every 5 seconds until we get the registry
    // do it like this for the proof of principle test
    ReadData data;
    do {
      CURL* han = curl_easy_init();

      if(han==0)
        {
          cerr << "could not create handle" << endl;
          //return 0; //or use this?
          throw cms::Exception("readHeader","EventStreamHttpReader")
            << "Could not get header: probably XDAQ not running on Storage Manager "
            << "\n";
        }

      setopt(han,CURLOPT_URL,headerurl_);
      setopt(han,CURLOPT_WRITEFUNCTION,func);
      setopt(han,CURLOPT_WRITEDATA,&data);

      if(curl_easy_perform(han)!=0)
      {
        cerr << "curl perform failed for header" << endl;
        //return 0; //or use this?
        throw cms::Exception("readHeader","EventStreamHttpReader")
          << "Could not get header: probably XDAQ not running on Storage Manager "
          << "\n";
      }
      curl_easy_cleanup(han);
      if(data.d_.length() == 0)
      {
        std::cout << "...waiting for Storage Manager..." << std::endl;
        // sleep for 5 seconds
        sleep(5);
      }
    } while (data.d_.length() == 0);

    //JobHeaderDecoder hdecoder;
    std::vector<char> regdata(1000*1000);

    // rely on http transfer string of correct length!
    int len = data.d_.length();
    FDEBUG(9) << "EventStreamHttpReader received registry len = " << len << std::endl;
    regdata.resize(len);
    for (int i=0; i<len ; i++) regdata[i] = data.d_[i];
    //edm::InitMsg msg(&regdata[0],len);
    InitMsgView initView(&regdata[0]);
    //hltBitCount = initView.hltTriggerCount();
    //l1BitCount = initView.l1TriggerCount();
    hltBitCount = initView.get_hlt_bit_cnt();
    l1BitCount = initView.get_l1_bit_cnt();
    // 21-Jun-2006, KAB:  catch (and re-throw) any exceptions decoding
    // the job header so that we can display the returned HTML and
    // (hopefully) give the user a hint as to the cause of the problem.
    std::auto_ptr<SendJobHeader> p;
    try {
      //p = hdecoder.decodeJobHeader(msg);
      p = StreamTranslator::deserializeRegistry(initView);
    }
    catch (cms::Exception excpt) {
      const unsigned int MAX_DUMP_LENGTH = 1000;
      std::cout << "========================================" << std::endl;
      std::cout << "* Exception decoding the getregdata response from the storage manager!" << std::endl;
      if (data.d_.length() <= MAX_DUMP_LENGTH) {
        std::cout << "* Here is the raw text that was returned:" << std::endl;
        std::cout << data.d_ << std::endl;
      }
      else {
        std::cout << "* Here are the first " << MAX_DUMP_LENGTH <<
          " characters of the raw text that was returned:" << std::endl;
        std::cout << (data.d_.substr(0, MAX_DUMP_LENGTH)) << std::endl;
      }
      std::cout << "========================================" << std::endl;
      throw excpt;
    }
    return p;
  }
}
