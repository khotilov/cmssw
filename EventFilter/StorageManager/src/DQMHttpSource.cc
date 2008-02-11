/** 
 *  An input source for DQM consumers run in cmsRun that connect to
 *  the StorageManager or SMProxyServer to get DQM data.
 *
 *  $Id: DQMHttpSource.cc,v 1.9 2008/01/22 19:28:36 muzaffar Exp $
 */

#include "EventFilter/StorageManager/src/DQMHttpSource.h"
#include "EventFilter/StorageManager/interface/SMCurlInterface.h"
#include "FWCore/Utilities/interface/DebugMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "IOPool/Streamer/interface/OtherMessage.h"
#include "IOPool/Streamer/interface/ConsRegMessage.h"
#include "IOPool/Streamer/interface/DQMEventMessage.h"
#include "IOPool/Streamer/interface/StreamDQMDeserializer.h"

#include <iostream>
#include <sys/time.h>
#include "curl/curl.h"
#include <wait.h>

using namespace edm;
using namespace std;

namespace edm
{  
  DQMHttpSource::DQMHttpSource(const ParameterSet& pset, 
                                         const InputSourceDescription& desc) :
    edm::RawInputSource(pset,desc), 
    updatesCounter_(0),
    sourceurl_(pset.getUntrackedParameter<string>("sourceURL")),
    buf_(1000*1000*7), 
    events_read_(0),
    consumerTopFolderName_(pset.getUntrackedParameter<string>("topLevelFolderName")),
    alreadySaidHalted_(false)
  {
    std::string evturl = sourceurl_ + "/getDQMeventdata";
    int stlen = evturl.length();
    for (int i=0; i<stlen; i++) DQMeventurl_[i]=evturl[i];
    DQMeventurl_[stlen] = '\0';

    std::string regurl = sourceurl_ + "/registerDQMConsumer";
    stlen = regurl.length();
    for (int i=0; i<stlen; i++) DQMsubscriptionurl_[i]=regurl[i];
    DQMsubscriptionurl_[stlen] = '\0';

    const double MAX_REQUEST_INTERVAL = 300.0;  // seconds
    DQMconsumerName_ = pset.getUntrackedParameter<string>("DQMconsumerName","Unknown");
    DQMconsumerPriority_ = pset.getUntrackedParameter<string>("DQMconsumerPriority","normal");
    headerRetryInterval_ = pset.getUntrackedParameter<int>("headerRetryInterval",5);
    double maxEventRequestRate = pset.getUntrackedParameter<double>("maxDQMEventRequestRate",1.0);
    if (maxEventRequestRate < (1.0 / MAX_REQUEST_INTERVAL)) {
      minDQMEventRequestInterval_ = MAX_REQUEST_INTERVAL;
    }
    else {
      minDQMEventRequestInterval_ = 1.0 / maxEventRequestRate;  // seconds
    }
    lastDQMRequestTime_.tv_sec = 0;
    lastDQMRequestTime_.tv_usec = 0;

    // register this DQM consumer with the DQMevent server of the Storage Manager
    DQMconsumerId_ = (time(0) & 0xffffff);  // temporary - will get from ES later
    registerWithDQMEventServer();
    // when running Async it seems bei_ is not NULL at the start after default ctor
    bei_ = NULL;
  }


  std::auto_ptr<Event> DQMHttpSource::readOneEvent()
  {
    // repeat a http get every X seconds until we get a DQMevent
    // only way to stop is specify a maxEvents parameter
    // or kill the Storage Manager XDAQ application so the http get fails.

    // try to get an event repeat until we get one, this allows
    // re-registration is the SM is halted or stopped

    bool gotEvent = false;
    std::auto_ptr<Event> result(0);
    while (!gotEvent)
    { 
       result = getOneDQMEvent();
       if(result.get() != NULL) gotEvent = true;
    } 
    return result;
  }

  std::auto_ptr<Event> DQMHttpSource::getOneDQMEvent()
  {
    // repeat a http get every X seconds until we get a DQMevent
    // only way to stop is specify a maxEvents parameter
    // or kill the Storage Manager XDAQ application so the http get fails.

    // check if we need to sleep (to enforce the allowed request rate)
    struct timeval now;
    struct timezone dummyTZ;
    gettimeofday(&now, &dummyTZ);
    double timeDiff = (double) now.tv_sec;
    timeDiff -= (double) lastDQMRequestTime_.tv_sec;
    timeDiff += ((double) now.tv_usec / 1000000.0);
    timeDiff -= ((double) lastDQMRequestTime_.tv_usec / 1000000.0);
    if (timeDiff < minDQMEventRequestInterval_)
    {
      double sleepTime = minDQMEventRequestInterval_ - timeDiff;
      // trim off a little sleep time to account for the time taken by
      // calling gettimeofday again
      sleepTime -= 0.01;
      if (sleepTime < 0.0) {sleepTime = 0.0;}
      usleep(static_cast<int>(1000000 * sleepTime));
      gettimeofday(&lastDQMRequestTime_, &dummyTZ);
    }
    else
    {
      lastDQMRequestTime_ = now;
    }

    stor::ReadData data;
    bool alreadySaidWaiting = false;
    do {
      CURL* han = curl_easy_init();

      if(han==0)
      {
        cerr << "DQMHttpSOurce: could not create handle" << endl;
        throw cms::Exception("getOneEvent","DQMHttpSource")
            << "Unable to create curl handle\n";
        // this will end cmsRun
      }

      stor::setopt(han,CURLOPT_URL,DQMeventurl_);
      stor::setopt(han,CURLOPT_WRITEFUNCTION,stor::func);
      stor::setopt(han,CURLOPT_WRITEDATA,&data);

      // send our consumer ID as part of the event request
      char msgBuff[100];
      OtherMessageBuilder requestMessage(&msgBuff[0], Header::DQMEVENT_REQUEST,
                                         sizeof(char_uint32));
      uint8 *bodyPtr = requestMessage.msgBody();
      char_uint32 convertedId;
      convert(DQMconsumerId_, convertedId);
      for (unsigned int idx = 0; idx < sizeof(char_uint32); idx++) {
        bodyPtr[idx] = convertedId[idx];
      }
      stor::setopt(han, CURLOPT_POSTFIELDS, requestMessage.startAddress());
      stor::setopt(han, CURLOPT_POSTFIELDSIZE, requestMessage.size());
      struct curl_slist *headers=NULL;
      headers = curl_slist_append(headers, "Content-Type: application/octet-stream");
      headers = curl_slist_append(headers, "Content-Transfer-Encoding: binary");
      stor::setopt(han, CURLOPT_HTTPHEADER, headers);

      // send the HTTP POST, read the reply, and cleanup before going on
      CURLcode messageStatus = curl_easy_perform(han);
      curl_slist_free_all(headers);
      curl_easy_cleanup(han);

      if(messageStatus!=0)
      {
        cerr << "curl perform failed for DQMevent, messageStatus = "
             << messageStatus << endl;
        throw cms::Exception("getOneDQMEvent","DQMHttpSource")
            << "Could not get event: probably XDAQ not running on Storage Manager "
            << "\n";
        // this will end cmsRun
      }
      if(data.d_.length() == 0)
      {
        if(!alreadySaidWaiting) {
          std::cout << "...waiting for DQMevent from Storage Manager..." << std::endl;
          alreadySaidWaiting = true;
        }
        // sleep for the standard request interval
        usleep(static_cast<int>(1000000 * minDQMEventRequestInterval_));
      }
    } while (data.d_.length() == 0);

    int len = data.d_.length();
    FDEBUG(9) << "DQMHttpSource received len = " << len << std::endl;
    buf_.resize(len);
    for (int i=0; i<len ; i++) buf_[i] = data.d_[i];

    OtherMessageView msgView(&buf_[0]);

    RunNumber_t iRun = 0;
    EventNumber_t iEvent = 0;
    TimeValue_t tStamp = 1;
    Timestamp timeStamp (tStamp);

    if (msgView.code() == Header::DONE) {
      // Continue past run boundaries (SM halt)
      // no need to register again as the SM/EventServer is kept alive on a stopAction
     if(!alreadySaidHalted_) {
       alreadySaidHalted_ = true;
       std::cout << "Storage Manager has halted - waiting for restart" << std::endl;
     }
     return std::auto_ptr<edm::Event>();
    } else {
      // counting the updates
      ++updatesCounter_;
      ++events_read_;
      DQMEventMsgView dqmEventView(&buf_[0]);
      iRun = dqmEventView.runNumber();
      iEvent = dqmEventView.eventNumberAtUpdate();
      timeStamp = dqmEventView.timeStamp();

      FDEBUG(8) << "  DQM Message data:" << std::endl;
      FDEBUG(8) << "    protocol version = "
                << dqmEventView.protocolVersion() << std::endl;
      FDEBUG(8) << "    header size = "
                << dqmEventView.headerSize() << std::endl;
      FDEBUG(8) << "    run number = "
                << dqmEventView.runNumber() << std::endl;
      FDEBUG(8) << "    event number = "
                << dqmEventView.eventNumberAtUpdate() << std::endl;
      FDEBUG(8) << "    lumi section = "
                << dqmEventView.lumiSection() << std::endl;
      FDEBUG(8) << "    update number = "
                << dqmEventView.updateNumber() << std::endl;
      FDEBUG(8) << "    compression flag = "
                << dqmEventView.compressionFlag() << std::endl;
      FDEBUG(8) << "    reserved word = "
                << dqmEventView.reserved() << std::endl;
      FDEBUG(8) << "    release tag = "
                << dqmEventView.releaseTag() << std::endl;
      FDEBUG(8) << "    top folder name = "
                << dqmEventView.topFolderName() << std::endl;
      FDEBUG(8) << "    sub folder count = "
                << dqmEventView.subFolderCount() << std::endl;

      // deserialize and stick into DQM backend
      // need both types of interfaces as the extractObject I use is
      // only in DaqMonitorBEInterface
      if (bei_ == NULL) {
        bei_ = dynamic_cast<DaqMonitorROOTBackEnd*>(edm::Service<DaqMonitorBEInterface>().operator->());
      }
      if (bei_ == NULL) {
        throw cms::Exception("readOneEvent", "DQMHttpSource")
          << "Unable to lookup the DaqMonitorBEInterface service!\n";
      }

      edm::StreamDQMDeserializer deserializeWorker;
      std::auto_ptr<DQMEvent::TObjectTable> toTablePtr =
          deserializeWorker.deserializeDQMEvent(dqmEventView);

      unsigned int count = 0;
      DQMEvent::TObjectTable::const_iterator toIter;
      for (toIter = toTablePtr->begin();
           toIter != toTablePtr->end(); toIter++) {
        std::string subFolderName = toIter->first;
        std::vector<TObject *> toList = toIter->second;
        MonitorElementRootFolder * dqmdir = 
          dynamic_cast<DaqMonitorBEInterface*>(bei_)->makeDirectory(subFolderName);  // fetch or create
        const bool fromRemoteNode = true; // Pretend this is like connecting to DQM Collector
        for (int tdx = 0; tdx < (int) toList.size(); tdx++) {
          TObject *toPtr = toList[tdx];
          std::string cls = toPtr->IsA()->GetName();
          std::string nm = toPtr->GetName();
          FDEBUG(8) << "    TObject class = " << cls << ", name = " << nm << std::endl;
          // special handling for TObjString - try to use DQM core code for this:
          if(bei_->isInt(toPtr) || bei_->isFloat(toPtr) ||
             bei_->isString(toPtr) || bei_->isQReport(toPtr)) {
            std::string tos_name, tos_value;
            bool extractOK = bei_->extractObject(toPtr, fromRemoteNode, tos_name, tos_value);
            if(extractOK) {
              addMonitorable(tos_name, subFolderName);
              setIsDesired(dqmdir, tos_name, true);
            }
          } else {
            addMonitorable(nm, subFolderName);
            setIsDesired(dqmdir, nm, true);
          }
          // for TObjString of type String this will complain
          bool success = bei_->extractObject(toPtr, dqmdir,fromRemoteNode);
          if(success) ++count; // success looks to be hardwired to be true on return!
        }
      }

      // clean up memory by spinning through the DQMEvent::TObjectTable map and
      // deleting each TObject in the std::vector<TObject *> later we will
      // change map to use std::vector< boost::shared_ptr<TObject> >
      DQMEvent::TObjectTable::iterator ti(toTablePtr->begin()), te(toTablePtr->end());
      for ( ; ti != te; ++ti) {
        std::vector<TObject *>::iterator vi(ti->second.begin()), ve(ti->second.end());
        for ( ; vi != ve; ++vi) delete *vi;
      }
    }

    setRunNumber(iRun);
    EventID eventId(iRun,iEvent);

    // make a fake event containing no data but the evId and runId from DQMEvent
    // and the time stamp from the event at update
    std::auto_ptr<Event> e = makeEvent(eventId,timeStamp);
  
    return e;
  }

  void DQMHttpSource::registerWithDQMEventServer()
  {
    stor::ReadData data;
    uint32 registrationStatus;
    bool alreadySaidWaiting = false;
    do {
      data.d_.clear();
      CURL* han = curl_easy_init();
      if(han==0)
        {
          cerr << "could not create handle" << endl;
          throw cms::Exception("registerWithDQMEventServer","DQMHttpSource")
            << "Unable to create curl handle\n";
        }

      // set the standard http request options
      stor::setopt(han,CURLOPT_URL,DQMsubscriptionurl_);
      stor::setopt(han,CURLOPT_WRITEFUNCTION,stor::func);
      stor::setopt(han,CURLOPT_WRITEDATA,&data);

      // build the registration request message to send to the storage manager
      const int BUFFER_SIZE = 2000;
      char msgBuff[BUFFER_SIZE];
      ConsRegRequestBuilder requestMessage(msgBuff, BUFFER_SIZE, DQMconsumerName_,
                                       DQMconsumerPriority_, consumerTopFolderName_);

      // add the request message as a http post
      stor::setopt(han, CURLOPT_POSTFIELDS, requestMessage.startAddress());
      stor::setopt(han, CURLOPT_POSTFIELDSIZE, requestMessage.size());
      struct curl_slist *headers=NULL;
      headers = curl_slist_append(headers, "Content-Type: application/octet-stream");
      headers = curl_slist_append(headers, "Content-Transfer-Encoding: binary");
      stor::setopt(han, CURLOPT_HTTPHEADER, headers);

      // send the HTTP POST, read the reply, and cleanup before going on
      CURLcode messageStatus = curl_easy_perform(han);
      curl_slist_free_all(headers);
      curl_easy_cleanup(han);

      if(messageStatus!=0)
      {
        cerr << "curl perform failed for DQM registration" << endl;
        throw cms::Exception("registerWithDQMEventServer","DQMHttpSource")
          << "Could not register: probably XDAQ not running or no Storage Manager/SMProxyServer loaded"
          << "\n";
      }
      registrationStatus = ConsRegResponseBuilder::ES_NOT_READY;
      if(data.d_.length() > 0)
      {
        int len = data.d_.length();
        FDEBUG(9) << "DQMHttpSource received len = " << len << std::endl;
        buf_.resize(len);
        for (int i=0; i<len ; i++) buf_[i] = data.d_[i];

        try {
          ConsRegResponseView respView(&buf_[0]);
          registrationStatus = respView.getStatus();
          DQMconsumerId_ = respView.getConsumerId();
        }
        catch (cms::Exception excpt) {
          const unsigned int MAX_DUMP_LENGTH = 1000;
          std::cout << "========================================" << std::endl;
          std::cout << "* Exception decoding the registerWithEventServer response!" << std::endl;
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
      }

      if (registrationStatus == ConsRegResponseBuilder::ES_NOT_READY)
      {
        if(!alreadySaidWaiting) {
          std::cout << "...waiting for DQM registration response from StorageManager or SMProxyServer..." 
                    << std::endl;
          alreadySaidWaiting = true;
        }
        // sleep for desired amount of time
        sleep(headerRetryInterval_);
      }
    } while (registrationStatus == ConsRegResponseBuilder::ES_NOT_READY);

    FDEBUG(9) << "Consumer ID = " << DQMconsumerId_ << endl;
  }

  void DQMHttpSource::addMonitorable(std::string name, std::string dir_path)
  {
    if (bei_ == NULL) {
       std::cerr << " *** DQMHttpSource::addMonitorable: No backend interface defined! " << std::endl;
       std::cerr << " Ignoring addMonitorable operation... " << std::endl;
       return;
    }
    if(!bei_->findObject(name, dir_path)) bei_->addElement(name, dir_path);

    bei_->lock();
    // make the monitorable string, <dir_path>:<obj1>,<obj2>,...:<tag> where tag is optional
    // see StringUtil:::unpackDirFormat
    std::string new_name = dir_path + ":" + name;
    bei_->addedMonitorable.push_back(new_name);
    bei_->unlock();
  }

  void DQMHttpSource::setIsDesired(MonitorElementRootFolder *folder, std::string ME_name, bool flag)
  {
    if(!folder)
    {
      std::cerr << " *** Null MonitorElementRootFolder in DQMHttpSource::setIsDesired"
           << std::endl;
      return;
    }

    if(!folder->hasMonitorable(ME_name))
    {
      std::cerr << " *** DQMHttpSource::setIsDesired: Object " << ME_name << " does not exist in "
           << folder->getPathname() << std::endl;
      std::cerr << " Ignoring setIsDesired operation... " << std::endl;
      return;
    }
  
    folder->isDesired[ME_name] = flag;
  }
}
